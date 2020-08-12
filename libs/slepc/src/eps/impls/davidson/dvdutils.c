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

   Some utils
*/

#include "davidson.h"

typedef struct {
  PC pc;
} dvdPCWrapper;

/*
  Configure the harmonics.
  switch (mode) {
  DVD_HARM_RR:    harmonic RR
  DVD_HARM_RRR:   relative harmonic RR
  DVD_HARM_REIGS: rightmost eigenvalues
  DVD_HARM_LEIGS: largest eigenvalues
  }
  fixedTarged, if true use the target instead of the best eigenvalue
  target, the fixed target to be used
*/
typedef struct {
  PetscScalar Wa,Wb;       /* span{W} = span{Wa*AV - Wb*BV} */
  PetscScalar Pa,Pb;       /* H=W'*(Pa*AV - Pb*BV), G=W'*(Wa*AV - Wb*BV) */
  PetscBool   withTarget;
  HarmType_t  mode;
} dvdHarmonic;

typedef struct {
  Vec diagA, diagB;
} dvdJacobiPrecond;

static PetscErrorCode dvd_improvex_precond_d(dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdPCWrapper   *dvdpc = (dvdPCWrapper*)d->improvex_precond_data;

  PetscFunctionBegin;
  /* Free local data */
  ierr = PCDestroy(&dvdpc->pc);CHKERRQ(ierr);
  ierr = PetscFree(d->improvex_precond_data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_static_precond_PC_0(dvdDashboard *d,PetscInt i,Vec x,Vec Px)
{
  PetscErrorCode ierr;
  dvdPCWrapper   *dvdpc = (dvdPCWrapper*)d->improvex_precond_data;

  PetscFunctionBegin;
  ierr = PCApply(dvdpc->pc,x,Px);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Create a trivial preconditioner
*/
static PetscErrorCode dvd_precond_none(dvdDashboard *d,PetscInt i,Vec x,Vec Px)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecCopy(x,Px);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Create a static preconditioner from a PC
*/
PetscErrorCode dvd_static_precond_PC(dvdDashboard *d,dvdBlackboard *b,PC pc)
{
  PetscErrorCode ierr;
  dvdPCWrapper   *dvdpc;
  Mat            P;
  PetscBool      t0,t1,t2;

  PetscFunctionBegin;
  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    /* If the preconditioner is valid */
    if (pc) {
      ierr = PetscNewLog(d->eps,&dvdpc);CHKERRQ(ierr);
      dvdpc->pc = pc;
      ierr = PetscObjectReference((PetscObject)pc);CHKERRQ(ierr);
      d->improvex_precond_data = dvdpc;
      d->improvex_precond = dvd_static_precond_PC_0;

      /* PC saves the matrix associated with the linear system, and it has to
         be initialize to a valid matrix */
      ierr = PCGetOperatorsSet(pc,NULL,&t0);CHKERRQ(ierr);
      ierr = PetscObjectTypeCompare((PetscObject)pc,PCNONE,&t1);CHKERRQ(ierr);
      ierr = PetscObjectTypeCompare((PetscObject)pc,PCSHELL,&t2);CHKERRQ(ierr);
      if (t0 && !t1) {
        ierr = PCGetOperators(pc,NULL,&P);CHKERRQ(ierr);
        ierr = PetscObjectReference((PetscObject)P);CHKERRQ(ierr);
        ierr = PCSetOperators(pc,P,P);CHKERRQ(ierr);
        ierr = PCSetReusePreconditioner(pc,PETSC_TRUE);CHKERRQ(ierr);
        ierr = MatDestroy(&P);CHKERRQ(ierr);
      } else if (t2) {
        ierr = PCSetOperators(pc,d->A,d->A);CHKERRQ(ierr);
        ierr = PCSetReusePreconditioner(pc,PETSC_TRUE);CHKERRQ(ierr);
      } else {
        d->improvex_precond = dvd_precond_none;
      }

      ierr = EPSDavidsonFLAdd(&d->destroyList,dvd_improvex_precond_d);CHKERRQ(ierr);

    /* Else, use no preconditioner */
    } else d->improvex_precond = dvd_precond_none;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_harm_d(dvdDashboard *d)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Free local data */
  ierr = PetscFree(d->calcpairs_W_data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_harm_transf(dvdHarmonic *dvdh,PetscScalar t)
{
  PetscFunctionBegin;
  switch (dvdh->mode) {
  case DVD_HARM_RR:    /* harmonic RR */
    dvdh->Wa = 1.0; dvdh->Wb = t;   dvdh->Pa = 0.0; dvdh->Pb = -1.0;
    break;
  case DVD_HARM_RRR:   /* relative harmonic RR */
    dvdh->Wa = 1.0; dvdh->Wb = t;   dvdh->Pa = 1.0; dvdh->Pb = 0.0;
    break;
  case DVD_HARM_REIGS: /* rightmost eigenvalues */
    dvdh->Wa = 1.0; dvdh->Wb = t;   dvdh->Pa = 1.0; dvdh->Pb = -PetscConj(t);
    break;
  case DVD_HARM_LEIGS: /* largest eigenvalues */
    dvdh->Wa = 0.0; dvdh->Wb = 1.0; dvdh->Pa = 1.0; dvdh->Pb = 0.0;
    break;
  case DVD_HARM_NONE:
  default:
    SETERRQ(PETSC_COMM_SELF,1,"Harmonic type not supported");
  }

  /* Check the transformation does not change the sign of the imaginary part */
#if !defined(PETSC_USE_COMPLEX)
  if (dvdh->Pb*dvdh->Wa - dvdh->Wb*dvdh->Pa < 0.0) {
    dvdh->Pa *= -1.0;
    dvdh->Pb *= -1.0;
  }
#endif
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_harm_updateW(dvdDashboard *d)
{
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;
  PetscErrorCode ierr;
  PetscInt       l,k;
  BV             BX = d->BX?d->BX:d->eps->V;

  PetscFunctionBegin;
  /* Update the target if it is necessary */
  if (!data->withTarget) {
    ierr = dvd_harm_transf(data,d->eigr[0]);CHKERRQ(ierr);
  }

  /* W(i) <- Wa*AV(i) - Wb*BV(i) */
  ierr = BVGetActiveColumns(d->eps->V,&l,&k);CHKERRQ(ierr);
  if (k != l+d->V_new_s) SETERRQ(PETSC_COMM_SELF,1,"Consistency broken");
  ierr = BVSetActiveColumns(d->W,l+d->V_new_s,l+d->V_new_e);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(d->AX,l+d->V_new_s,l+d->V_new_e);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(BX,l+d->V_new_s,l+d->V_new_e);CHKERRQ(ierr);
  ierr = BVCopy(d->AX,d->W);CHKERRQ(ierr);
  ierr = BVScale(d->W,data->Wa);CHKERRQ(ierr);
  ierr = BVMult(d->W,-data->Wb,1.0,BX,NULL);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(d->W,l,k);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(d->AX,l,k);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(BX,l,k);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_harm_proj(dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;
  PetscInt       i,j,l0,l,k,ld;
  PetscScalar    h,g,*H,*G;

  PetscFunctionBegin;
  ierr = BVGetActiveColumns(d->eps->V,&l0,&k);CHKERRQ(ierr);
  l = l0 + d->V_new_s;
  k = l0 + d->V_new_e;
  ierr = MatGetSize(d->H,&ld,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(d->H,&H);CHKERRQ(ierr);
  ierr = MatDenseGetArray(d->G,&G);CHKERRQ(ierr);
  /* [H G] <- [Pa*H - Pb*G, Wa*H - Wb*G] */
  /* Right part */
  for (i=l;i<k;i++) {
    for (j=l0;j<k;j++) {
      h = H[ld*i+j];
      g = G[ld*i+j];
      H[ld*i+j] = data->Pa*h - data->Pb*g;
      G[ld*i+j] = data->Wa*h - data->Wb*g;
    }
  }
  /* Left part */
  for (i=l0;i<l;i++) {
    for (j=l;j<k;j++) {
      h = H[ld*i+j];
      g = G[ld*i+j];
      H[ld*i+j] = data->Pa*h - data->Pb*g;
      G[ld*i+j] = data->Wa*h - data->Wb*g;
    }
  }
  ierr = MatDenseRestoreArray(d->H,&H);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(d->G,&G);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_harm_updateproj(dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;
  PetscInt       i,j,l,k,ld;
  PetscScalar    h,g,*H,*G;

  PetscFunctionBegin;
  ierr = BVGetActiveColumns(d->eps->V,&l,&k);CHKERRQ(ierr);
  k = l + d->V_tra_s;
  ierr = MatGetSize(d->H,&ld,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(d->H,&H);CHKERRQ(ierr);
  ierr = MatDenseGetArray(d->G,&G);CHKERRQ(ierr);
  /* [H G] <- [Pa*H - Pb*G, Wa*H - Wb*G] */
  /* Right part */
  for (i=l;i<k;i++) {
    for (j=0;j<l;j++) {
      h = H[ld*i+j];
      g = G[ld*i+j];
      H[ld*i+j] = data->Pa*h - data->Pb*g;
      G[ld*i+j] = data->Wa*h - data->Wb*g;
    }
  }
  /* Lower triangular part */
  for (i=0;i<l;i++) {
    for (j=l;j<k;j++) {
      h = H[ld*i+j];
      g = G[ld*i+j];
      H[ld*i+j] = data->Pa*h - data->Pb*g;
      G[ld*i+j] = data->Wa*h - data->Wb*g;
    }
  }
  ierr = MatDenseRestoreArray(d->H,&H);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(d->G,&G);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_harm_backtrans(dvdHarmonic *data,PetscScalar *ar,PetscScalar *ai)
{
  PetscScalar xr;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar xi, k;
#endif

  PetscFunctionBegin;
  xr = *ar;
#if !defined(PETSC_USE_COMPLEX)
  xi = *ai;
  if (xi != 0.0) {
    k = (data->Pa - data->Wa*xr)*(data->Pa - data->Wa*xr) + data->Wa*data->Wa*xi*xi;
    *ar = (data->Pb*data->Pa - (data->Pb*data->Wa + data->Wb*data->Pa)*xr + data->Wb*data->Wa*(xr*xr + xi*xi))/k;
    *ai = (data->Pb*data->Wa - data->Wb*data->Pa)*xi/k;
  } else
#endif
    *ar = (data->Pb - data->Wb*xr) / (data->Pa - data->Wa*xr);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_harm_eig_backtrans(dvdDashboard *d,PetscScalar ar,PetscScalar ai,PetscScalar *br,PetscScalar *bi)
{
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = dvd_harm_backtrans(data,&ar,&ai);CHKERRQ(ierr);
  *br = ar;
  *bi = ai;
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_harm_eigs_trans(dvdDashboard *d)
{
  dvdHarmonic    *data = (dvdHarmonic*)d->calcpairs_W_data;
  PetscInt       i,l,k;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = BVGetActiveColumns(d->eps->V,&l,&k);CHKERRQ(ierr);
  for (i=0;i<k-l;i++) {
    ierr = dvd_harm_backtrans(data,&d->eigr[i],&d->eigi[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_harm_conf(dvdDashboard *d,dvdBlackboard *b,HarmType_t mode,PetscBool fixedTarget,PetscScalar t)
{
  PetscErrorCode ierr;
  dvdHarmonic    *dvdh;

  PetscFunctionBegin;
  /* Set the problem to GNHEP:
     d->G maybe is upper triangular due to biorthogonality of V and W */
  d->sEP = d->sA = d->sB = 0;

  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    ierr = PetscNewLog(d->eps,&dvdh);CHKERRQ(ierr);
    dvdh->withTarget = fixedTarget;
    dvdh->mode = mode;
    if (fixedTarget) dvd_harm_transf(dvdh, t);
    d->calcpairs_W_data = dvdh;
    d->calcpairs_W = dvd_harm_updateW;
    d->calcpairs_proj_trans = dvd_harm_proj;
    d->calcpairs_eigs_trans = dvd_harm_eigs_trans;
    d->calcpairs_eig_backtrans = dvd_harm_eig_backtrans;

    ierr = EPSDavidsonFLAdd(&d->destroyList,dvd_harm_d);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

