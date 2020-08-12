/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/vecimpl.h>

#if defined(__WITH_MPI__)
#define __SUF__(A) A##_MPI
#else
#define __SUF__(A) A##_Seq
#endif
#define __QUOTEME_(x) #x
#define __QUOTEME(x) __QUOTEME_(x)
#define __SUF_C__(A) __QUOTEME(__SUF__(A))

PetscErrorCode __SUF__(VecDot_Comp)(Vec a,Vec b,PetscScalar *z)
{
  PetscScalar    sum = 0.0,work;
  PetscInt       i;
  PetscErrorCode ierr;
  Vec_Comp       *as = (Vec_Comp*)a->data,*bs = (Vec_Comp*)b->data;

  PetscFunctionBegin;
  SlepcValidVecComp(a,1);
  SlepcValidVecComp(b,2);
  if (as->x[0]->ops->dot_local) {
    for (i=0,sum=0.0;i<as->n->n;i++) {
      ierr = as->x[i]->ops->dot_local(as->x[i],bs->x[i],&work);CHKERRQ(ierr);
      sum += work;
    }
#if defined(__WITH_MPI__)
    work = sum;
    ierr = MPI_Allreduce(&work,&sum,1,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)a));CHKERRQ(ierr);
#endif
  } else {
    for (i=0,sum=0.0;i<as->n->n;i++) {
      ierr = VecDot(as->x[i],bs->x[i],&work);CHKERRQ(ierr);
      sum += work;
    }
  }
  *z = sum;
  PetscFunctionReturn(0);
}

PetscErrorCode __SUF__(VecMDot_Comp)(Vec a,PetscInt n,const Vec b[],PetscScalar *z)
{
  PetscScalar    *work,*work0,*r;
  PetscErrorCode ierr;
  Vec_Comp       *as = (Vec_Comp*)a->data;
  Vec            *bx;
  PetscInt       i,j;

  PetscFunctionBegin;
  SlepcValidVecComp(a,1);
  for (i=0;i<n;i++) SlepcValidVecComp(b[i],3);

  if (as->n->n == 0) {
    *z = 0;
    PetscFunctionReturn(0);
  }

  ierr = PetscMalloc2(n,&work0,n,&bx);CHKERRQ(ierr);

#if defined(__WITH_MPI__)
  if (as->x[0]->ops->mdot_local) {
    r = work0; work = z;
  } else
#endif
  {
    r = z; work = work0;
  }

  /* z[i] <- a.x' * b[i].x */
  for (i=0;i<n;i++) r[i] = 0.0;
  for (j=0;j<as->n->n;j++) {
    for (i=0;i<n;i++) bx[i] = ((Vec_Comp*)b[i]->data)->x[j];
    if (as->x[0]->ops->mdot_local) {
      ierr = as->x[j]->ops->mdot_local(as->x[j],n,bx,work);CHKERRQ(ierr);
    } else {
      ierr = VecMDot(as->x[j],n,bx,work);CHKERRQ(ierr);
    }
    for (i=0;i<n;i++) r[i] += work[i];
  }

  /* If def(__WITH_MPI__) and exists mdot_local */
#if defined(__WITH_MPI__)
  if (as->x[0]->ops->mdot_local) {
    /* z[i] <- Allreduce(work[i]) */
    ierr = MPI_Allreduce(r,z,n,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)a));CHKERRQ(ierr);
  }
#endif

  ierr = PetscFree2(work0,bx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode __SUF__(VecTDot_Comp)(Vec a,Vec b,PetscScalar *z)
{
  PetscScalar    sum = 0.0,work;
  PetscInt       i;
  PetscErrorCode ierr;
  Vec_Comp       *as = (Vec_Comp*)a->data,*bs = (Vec_Comp*)b->data;

  PetscFunctionBegin;
  SlepcValidVecComp(a,1);
  SlepcValidVecComp(b,2);
  if (as->x[0]->ops->tdot_local) {
    for (i=0,sum=0.0;i<as->n->n;i++) {
      ierr = as->x[i]->ops->tdot_local(as->x[i],bs->x[i],&work);CHKERRQ(ierr);
      sum += work;
    }
#if defined(__WITH_MPI__)
    work = sum;
    ierr = MPI_Allreduce(&work,&sum,1,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)a));CHKERRQ(ierr);
#endif
  } else {
    for (i=0,sum=0.0;i<as->n->n;i++) {
      ierr = VecTDot(as->x[i],bs->x[i],&work);CHKERRQ(ierr);
      sum += work;
    }
  }
  *z = sum;
  PetscFunctionReturn(0);
}

PetscErrorCode __SUF__(VecMTDot_Comp)(Vec a,PetscInt n,const Vec b[],PetscScalar *z)
{
  PetscScalar    *work,*work0,*r;
  PetscErrorCode ierr;
  Vec_Comp       *as = (Vec_Comp*)a->data;
  Vec            *bx;
  PetscInt       i,j;

  PetscFunctionBegin;
  SlepcValidVecComp(a,1);
  for (i=0;i<n;i++) SlepcValidVecComp(b[i],3);

  if (as->n->n == 0) {
    *z = 0;
    PetscFunctionReturn(0);
  }

  ierr = PetscMalloc2(n,&work0,n,&bx);CHKERRQ(ierr);

#if defined(__WITH_MPI__)
  if (as->x[0]->ops->mtdot_local) {
    r = work0; work = z;
  } else
#endif
  {
    r = z; work = work0;
  }

  /* z[i] <- a.x' * b[i].x */
  for (i=0;i<n;i++) r[i] = 0.0;
  for (j=0;j<as->n->n;j++) {
    for (i=0;i<n;i++) bx[i] = ((Vec_Comp*)b[i]->data)->x[j];
    if (as->x[0]->ops->mtdot_local) {
      ierr = as->x[j]->ops->mtdot_local(as->x[j],n,bx,work);CHKERRQ(ierr);
    } else {
      ierr = VecMTDot(as->x[j],n,bx,work);CHKERRQ(ierr);
    }
    for (i=0;i<n;i++) r[i] += work[i];
  }

  /* If def(__WITH_MPI__) and exists mtdot_local */
#if defined(__WITH_MPI__)
  if (as->x[0]->ops->mtdot_local) {
    /* z[i] <- Allreduce(work[i]) */
    ierr = MPI_Allreduce(r,z,n,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)a));CHKERRQ(ierr);
  }
#endif

  ierr = PetscFree2(work0,bx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode __SUF__(VecNorm_Comp)(Vec a,NormType t,PetscReal *norm)
{
  PetscReal      work[3],s=0.0;
  PetscErrorCode ierr;
  Vec_Comp       *as = (Vec_Comp*)a->data;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(a,1);
  /* Initialize norm */
  switch (t) {
    case NORM_1: case NORM_INFINITY: *norm = 0.0; break;
    case NORM_2: case NORM_FROBENIUS: *norm = 1.0; s = 0.0; break;
    case NORM_1_AND_2: norm[0] = 0.0; norm[1] = 1.0; s = 0.0; break;
  }
  for (i=0;i<as->n->n;i++) {
    if (as->x[0]->ops->norm_local) {
      ierr = as->x[0]->ops->norm_local(as->x[i],t,work);CHKERRQ(ierr);
    } else {
      ierr = VecNorm(as->x[i],t,work);CHKERRQ(ierr);
    }
    /* norm+= work */
    switch (t) {
      case NORM_1: *norm+= *work; break;
      case NORM_2: case NORM_FROBENIUS: AddNorm2(norm,&s,*work); break;
      case NORM_1_AND_2: norm[0]+= work[0]; AddNorm2(&norm[1],&s,work[1]); break;
      case NORM_INFINITY: *norm = PetscMax(*norm,*work); break;
    }
  }

  /* If def(__WITH_MPI__) and exists norm_local */
#if defined(__WITH_MPI__)
  if (as->x[0]->ops->norm_local) {
    PetscReal work0[3];
    /* norm <- Allreduce(work) */
    switch (t) {
    case NORM_1:
      work[0] = *norm;
      ierr = MPI_Allreduce(work,norm,1,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)a));CHKERRQ(ierr);
      break;
    case NORM_2: case NORM_FROBENIUS:
      work[0] = *norm; work[1] = s;
      ierr = MPI_Allreduce(work,work0,1,MPIU_NORM2,MPIU_NORM2_SUM,PetscObjectComm((PetscObject)a));CHKERRQ(ierr);
      *norm = GetNorm2(work0[0],work0[1]);
      break;
    case NORM_1_AND_2:
      work[0] = norm[0]; work[1] = norm[1]; work[2] = s;
      ierr = MPI_Allreduce(work,work0,1,MPIU_NORM1_AND_2,MPIU_NORM2_SUM,PetscObjectComm((PetscObject)a));CHKERRQ(ierr);
      norm[0] = work0[0];
      norm[1] = GetNorm2(work0[1],work0[2]);
      break;
    case NORM_INFINITY:
      work[0] = *norm;
      ierr = MPI_Allreduce(work,norm,1,MPIU_REAL,MPIU_MAX,PetscObjectComm((PetscObject)a));CHKERRQ(ierr);
      break;
    }
  }
#else
  /* Norm correction */
  switch (t) {
    case NORM_2: case NORM_FROBENIUS: *norm = GetNorm2(*norm,s); break;
    case NORM_1_AND_2: norm[1] = GetNorm2(norm[1],s); break;
    default: ;
  }
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode __SUF__(VecDotNorm2_Comp)(Vec v,Vec w,PetscScalar *dp,PetscScalar *nm)
{
  PetscErrorCode    ierr;
  PetscScalar       dp0,nm0,dp1,nm1;
  const PetscScalar *vx,*wx;
  Vec_Comp          *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data;
  PetscInt          i,n;
  PetscBool         t0,t1;
#if defined(__WITH_MPI__)
  PetscScalar       work[4];
#endif

  PetscFunctionBegin;
  /* Compute recursively the local part */
  dp0 = nm0 = 0.0;
  ierr = PetscObjectTypeCompare((PetscObject)v,VECCOMP,&t0);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)w,VECCOMP,&t1);CHKERRQ(ierr);
  if (t0 && t1) {
    SlepcValidVecComp(v,1);
    SlepcValidVecComp(w,2);
    for (i=0;i<vs->n->n;i++) {
      ierr = VecDotNorm2_Comp_Seq(vs->x[i],ws->x[i],&dp1,&nm1);CHKERRQ(ierr);
      dp0 += dp1;
      nm0 += nm1;
    }
  } else if (!t0 && !t1) {
    ierr = VecGetLocalSize(v,&n);CHKERRQ(ierr);
    ierr = VecGetArrayRead(v,&vx);CHKERRQ(ierr);
    ierr = VecGetArrayRead(w,&wx);CHKERRQ(ierr);
    for (i=0;i<n;i++) {
      dp0 += vx[i]*PetscConj(wx[i]);
      nm0 += wx[i]*PetscConj(wx[i]);
    }
    ierr = VecRestoreArrayRead(v,&vx);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(w,&wx);CHKERRQ(ierr);
  } else SETERRQ(PetscObjectComm((PetscObject)v),PETSC_ERR_ARG_INCOMP,"Incompatible vector types");

#if defined(__WITH_MPI__)
    /* [dp, nm] <- Allreduce([dp0, nm0]) */
    work[0] = dp0; work[1] = nm0;
    ierr = MPI_Allreduce(work,&work[2],2,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)v));CHKERRQ(ierr);
    *dp = work[2]; *nm = work[3];
#else
    *dp = dp0; *nm = nm0;
#endif
  PetscFunctionReturn(0);
}

#undef __SUF__
#undef __QUOTEME
#undef __SUF_C__

