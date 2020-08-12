/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Shift spectral transformation, applies (A + sigma I) as operator, or
   inv(B)(A + sigma B) for generalized problems
*/

#include <slepc/private/stimpl.h>

PetscErrorCode STBackTransform_Shift(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscInt j;

  PetscFunctionBegin;
  for (j=0;j<n;j++) {
    eigr[j] += st->sigma;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode STPostSolve_Shift(ST st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (st->matmode == ST_MATMODE_INPLACE) {
    if (st->nmat>1) {
      ierr = MatAXPY(st->A[0],st->sigma,st->A[1],st->str);CHKERRQ(ierr);
    } else {
      ierr = MatShift(st->A[0],st->sigma);CHKERRQ(ierr);
    }
    st->Astate[0] = ((PetscObject)st->A[0])->state;
    st->state   = ST_STATE_INITIAL;
    st->opready = PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}

/*
   Operator (shift):
               Op               P         M
   if nmat=1:  A-sI             NULL      A-sI
   if nmat=2:  B^-1 (A-sB)      B         A-sB
*/
PetscErrorCode STComputeOperator_Shift(ST st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  st->usesksp = (st->nmat>1)? PETSC_TRUE: PETSC_FALSE;
  ierr = PetscObjectReference((PetscObject)st->A[1]);CHKERRQ(ierr);
  ierr = MatDestroy(&st->T[1]);CHKERRQ(ierr);
  st->T[1] = st->A[1];
  ierr = STMatMAXPY_Private(st,-st->sigma,0.0,0,NULL,PetscNot(st->state==ST_STATE_UPDATED),&st->T[0]);CHKERRQ(ierr);
  if (st->nmat>1) { ierr = PetscObjectReference((PetscObject)st->T[1]);CHKERRQ(ierr); }
  ierr = MatDestroy(&st->P);CHKERRQ(ierr);
  st->P = (st->nmat>1)? st->T[1]: NULL;
  st->M = st->T[0];
  PetscFunctionReturn(0);
}

PetscErrorCode STSetUp_Shift(ST st)
{
  PetscErrorCode ierr;
  PetscInt       k,nc,nmat=st->nmat;
  PetscScalar    *coeffs=NULL;

  PetscFunctionBegin;
  if (nmat>1) {
    ierr = STSetWorkVecs(st,1);CHKERRQ(ierr);
  }
  if (nmat>2) {  /* set-up matrices for polynomial eigenproblems */
    if (st->transform) {
      nc = (nmat*(nmat+1))/2;
      ierr = PetscMalloc1(nc,&coeffs);CHKERRQ(ierr);
      /* Compute coeffs */
      ierr = STCoeffs_Monomial(st,coeffs);CHKERRQ(ierr);
      /* T[n] = A_n */
      k = nmat-1;
      ierr = PetscObjectReference((PetscObject)st->A[k]);CHKERRQ(ierr);
      ierr = MatDestroy(&st->T[k]);CHKERRQ(ierr);
      st->T[k] = st->A[k];
      for (k=0;k<nmat-1;k++) {
        ierr = STMatMAXPY_Private(st,nmat>2?st->sigma:-st->sigma,0.0,k,coeffs?coeffs+((nmat-k)*(nmat-k-1))/2:NULL,PetscNot(st->state==ST_STATE_UPDATED),&st->T[k]);CHKERRQ(ierr);
      }
      ierr = PetscFree(coeffs);CHKERRQ(ierr);
      ierr = PetscObjectReference((PetscObject)st->T[nmat-1]);CHKERRQ(ierr);
      ierr = MatDestroy(&st->P);CHKERRQ(ierr);
      st->P = st->T[nmat-1];
      if (!st->ksp) { ierr = STGetKSP(st,&st->ksp);CHKERRQ(ierr); }
      ierr = STCheckFactorPackage(st);CHKERRQ(ierr);
      ierr = KSPSetOperators(st->ksp,st->P,st->P);CHKERRQ(ierr);
    } else {
      for (k=0;k<nmat;k++) {
        ierr = PetscObjectReference((PetscObject)st->A[k]);CHKERRQ(ierr);
        ierr = MatDestroy(&st->T[k]);CHKERRQ(ierr);
        st->T[k] = st->A[k];
      }
    }
  }
  if (st->P) {
    ierr = KSPSetUp(st->ksp);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode STSetShift_Shift(ST st,PetscScalar newshift)
{
  PetscErrorCode ierr;
  PetscInt       k,nc,nmat=PetscMax(st->nmat,2);
  PetscScalar    *coeffs=NULL;

  PetscFunctionBegin;
  if (st->transform) {
    if (st->matmode == ST_MATMODE_COPY && nmat>2) {
      nc = (nmat*(nmat+1))/2;
      ierr = PetscMalloc1(nc,&coeffs);CHKERRQ(ierr);
      /* Compute coeffs */
      ierr = STCoeffs_Monomial(st,coeffs);CHKERRQ(ierr);
    }
    for (k=0;k<nmat-1;k++) {
      ierr = STMatMAXPY_Private(st,nmat>2?newshift:-newshift,nmat>2?st->sigma:-st->sigma,k,coeffs?coeffs+((nmat-k)*(nmat-k-1))/2:NULL,PETSC_FALSE,&st->T[k]);CHKERRQ(ierr);
    }
    if (st->matmode == ST_MATMODE_COPY && nmat>2) {
      ierr = PetscFree(coeffs);CHKERRQ(ierr);
    }
    if (st->nmat<=2) st->M = st->T[0];
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode STCreate_Shift(ST st)
{
  PetscFunctionBegin;
  st->usesksp = PETSC_TRUE;

  st->ops->apply           = STApply_Generic;
  st->ops->applytrans      = STApplyTranspose_Generic;
  st->ops->backtransform   = STBackTransform_Shift;
  st->ops->setshift        = STSetShift_Shift;
  st->ops->getbilinearform = STGetBilinearForm_Default;
  st->ops->setup           = STSetUp_Shift;
  st->ops->computeoperator = STComputeOperator_Shift;
  st->ops->postsolve       = STPostSolve_Shift;
  st->ops->setdefaultksp   = STSetDefaultKSP_Default;
  PetscFunctionReturn(0);
}
