/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "davidson"

   Step: initialize subspace V
*/

#include "davidson.h"

typedef struct {
  PetscInt k;                 /* desired initial subspace size */
  PetscInt user;              /* number of user initial vectors */
  void     *old_initV_data;   /* old initV data */
} dvdInitV;

static PetscErrorCode dvd_initV_classic_0(dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdInitV       *data = (dvdInitV*)d->initV_data;
  PetscInt       i,user = PetscMin(data->user,d->eps->mpd), l,k;

  PetscFunctionBegin;
  ierr = BVGetActiveColumns(d->eps->V,&l,&k);CHKERRQ(ierr);
  /* User vectors are added at the beginning, so no active column should be in V */
  if (data->user>0&&l>0) SETERRQ(PETSC_COMM_SELF,1,"Consistency broken");
  /* Generate a set of random initial vectors and orthonormalize them */
  for (i=l+user;i<l+data->k && i<d->eps->ncv && i-l<d->eps->mpd;i++) {
    ierr = BVSetRandomColumn(d->eps->V,i);CHKERRQ(ierr);
  }
  d->V_tra_s = 0; d->V_tra_e = 0;
  d->V_new_s = 0; d->V_new_e = i-l;

  /* After that the user vectors will be destroyed */
  data->user = 0;
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_initV_krylov_0(dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdInitV       *data = (dvdInitV*)d->initV_data;
  PetscInt       i,user = PetscMin(data->user,d->eps->mpd),l,k;
  Vec            av,v1,v2;

  PetscFunctionBegin;
  ierr = BVGetActiveColumns(d->eps->V,&l,&k);CHKERRQ(ierr);
  /* User vectors are added at the beginning, so no active column should be in V */
  if (data->user>0&&l>0) SETERRQ(PETSC_COMM_SELF,1,"Consistency broken");

  /* If needed, generate a random vector for starting the arnoldi method */
  if (user == 0) {
    ierr = BVSetRandomColumn(d->eps->V,l);CHKERRQ(ierr);
    user = 1;
  }

  /* Perform k steps of Arnoldi with the operator K^{-1}*(t[1]*A-t[2]*B) */
  ierr = dvd_orthV(d->eps->V,l,l+user);CHKERRQ(ierr);
  for (i=l+user;i<l+data->k && i<d->eps->ncv && i-l<d->eps->mpd;i++) {
    /* aux <- theta[1]A*in - theta[0]*B*in */
    ierr = BVGetColumn(d->eps->V,i,&v1);CHKERRQ(ierr);
    ierr = BVGetColumn(d->eps->V,i-user,&v2);CHKERRQ(ierr);
    ierr = BVGetColumn(d->auxBV,0,&av);CHKERRQ(ierr);
    if (d->B) {
      ierr = MatMult(d->A,v2,v1);CHKERRQ(ierr);
      ierr = MatMult(d->B,v2,av);CHKERRQ(ierr);
      ierr = VecAXPBY(av,d->target[1],-d->target[0],v1);CHKERRQ(ierr);
    } else {
      ierr = MatMult(d->A,v2,av);CHKERRQ(ierr);
      ierr = VecAXPBY(av,-d->target[0],d->target[1],v2);CHKERRQ(ierr);
    }
    ierr = d->improvex_precond(d,0,av,v1);CHKERRQ(ierr);
    ierr = BVRestoreColumn(d->eps->V,i,&v1);CHKERRQ(ierr);
    ierr = BVRestoreColumn(d->eps->V,i-user,&v2);CHKERRQ(ierr);
    ierr = BVRestoreColumn(d->auxBV,0,&av);CHKERRQ(ierr);
    ierr = dvd_orthV(d->eps->V,i,i+1);CHKERRQ(ierr);
  }

  d->V_tra_s = 0; d->V_tra_e = 0;
  d->V_new_s = 0; d->V_new_e = i-l;

  /* After that the user vectors will be destroyed */
  data->user = 0;
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_initV_d(dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdInitV       *data = (dvdInitV*)d->initV_data;

  PetscFunctionBegin;
  /* Restore changes in dvdDashboard */
  d->initV_data = data->old_initV_data;

  /* Free local data */
  ierr = PetscFree(data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_initV(dvdDashboard *d, dvdBlackboard *b, PetscInt k,PetscInt user, PetscBool krylov)
{
  PetscErrorCode ierr;
  dvdInitV       *data;

  PetscFunctionBegin;
  /* Setting configuration constrains */
  b->max_size_V = PetscMax(b->max_size_V, k);

  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    ierr = PetscNewLog(d->eps,&data);CHKERRQ(ierr);
    data->k = k;
    data->user = user;
    data->old_initV_data = d->initV_data;
    d->initV_data = data;
    if (krylov) d->initV = dvd_initV_krylov_0;
    else d->initV = dvd_initV_classic_0;
    ierr = EPSDavidsonFLAdd(&d->destroyList,dvd_initV_d);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_orthV(BV V,PetscInt V_new_s,PetscInt V_new_e)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  for (i=V_new_s;i<V_new_e;i++) {
    ierr = BVOrthonormalizeColumn(V,i,PETSC_TRUE,NULL,NULL);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

