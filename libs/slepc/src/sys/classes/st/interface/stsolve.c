/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   ST interface routines, callable by users
*/

#include <slepc/private/stimpl.h>            /*I "slepcst.h" I*/

PetscErrorCode STApply_Generic(ST st,Vec x,Vec y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (st->M && st->P) {
    ierr = MatMult(st->M,x,st->work[0]);CHKERRQ(ierr);
    ierr = STMatSolve(st,st->work[0],y);CHKERRQ(ierr);
  } else if (st->M) {
    ierr = MatMult(st->M,x,y);CHKERRQ(ierr);
  } else {
    ierr = STMatSolve(st,x,y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   STApply - Applies the spectral transformation operator to a vector, for
   instance (A - sB)^-1 B in the case of the shift-and-invert transformation
   and generalized eigenproblem.

   Collective on st

   Input Parameters:
+  st - the spectral transformation context
-  x  - input vector

   Output Parameter:
.  y - output vector

   Level: developer

.seealso: STApplyTranspose(), STApplyHermitianTranspose()
@*/
PetscErrorCode STApply(ST st,Vec x,Vec y)
{
  PetscErrorCode ierr;
  Mat            Op;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  if (x == y) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  ierr = VecSetErrorIfLocked(y,3);CHKERRQ(ierr);
  if (!st->ops->apply) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"ST does not have apply");
  ierr = STGetOperator_Private(st,&Op);CHKERRQ(ierr);
  ierr = MatMult(Op,x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STApplyTranspose_Generic(ST st,Vec x,Vec y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (st->M && st->P) {
    ierr = STMatSolveTranspose(st,x,st->work[0]);CHKERRQ(ierr);
    ierr = MatMultTranspose(st->M,st->work[0],y);CHKERRQ(ierr);
  } else if (st->M) {
    ierr = MatMultTranspose(st->M,x,y);CHKERRQ(ierr);
  } else {
    ierr = STMatSolveTranspose(st,x,y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   STApplyTranspose - Applies the transpose of the operator to a vector, for
   instance B^T(A - sB)^-T in the case of the shift-and-invert transformation
   and generalized eigenproblem.

   Collective on st

   Input Parameters:
+  st - the spectral transformation context
-  x  - input vector

   Output Parameter:
.  y - output vector

   Level: developer

.seealso: STApply(), STApplyHermitianTranspose()
@*/
PetscErrorCode STApplyTranspose(ST st,Vec x,Vec y)
{
  PetscErrorCode ierr;
  Mat            Op;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  if (x == y) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  ierr = VecSetErrorIfLocked(y,3);CHKERRQ(ierr);
  if (!st->ops->applytrans) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"ST does not have applytrans");
  ierr = STGetOperator_Private(st,&Op);CHKERRQ(ierr);
  ierr = MatMultTranspose(Op,x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   STApplyHermitianTranspose - Applies the hermitian-transpose of the operator
   to a vector, for instance B^H(A - sB)^-H in the case of the shift-and-invert
   transformation and generalized eigenproblem.

   Collective on st

   Input Parameters:
+  st - the spectral transformation context
-  x  - input vector

   Output Parameter:
.  y - output vector

   Note:
   Currently implemented via STApplyTranspose() with appropriate conjugation.

   Level: developer

.seealso: STApply(), STApplyTranspose()
@*/
PetscErrorCode STApplyHermitianTranspose(ST st,Vec x,Vec y)
{
  PetscErrorCode ierr;
  Mat            Op;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  if (x == y) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  ierr = VecSetErrorIfLocked(y,3);CHKERRQ(ierr);
  if (!st->ops->applytrans) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"ST does not have applytrans");
  ierr = STGetOperator_Private(st,&Op);CHKERRQ(ierr);
  ierr = MatMultHermitianTranspose(Op,x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   STGetBilinearForm - Returns the matrix used in the bilinear form with a
   generalized problem with semi-definite B.

   Not collective, though a parallel Mat may be returned

   Input Parameters:
.  st - the spectral transformation context

   Output Parameter:
.  B - output matrix

   Notes:
   The output matrix B must be destroyed after use. It will be NULL in
   case of standard eigenproblems.

   Level: developer
@*/
PetscErrorCode STGetBilinearForm(ST st,Mat *B)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  PetscValidPointer(B,2);
  STCheckMatrices(st,1);
  ierr = (*st->ops->getbilinearform)(st,B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STGetBilinearForm_Default(ST st,Mat *B)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (st->nmat==1) *B = NULL;
  else {
    *B = st->A[1];
    ierr = PetscObjectReference((PetscObject)*B);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_STOperator(Mat Op,Vec x,Vec y)
{
  PetscErrorCode ierr;
  ST             st;

  PetscFunctionBegin;
  ierr = MatShellGetContext(Op,(void**)&st);CHKERRQ(ierr);
  ierr = STSetUp(st);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(ST_Apply,st,x,y,0);CHKERRQ(ierr);
  if (st->D) { /* with balancing */
    ierr = VecPointwiseDivide(st->wb,x,st->D);CHKERRQ(ierr);
    ierr = (*st->ops->apply)(st,st->wb,y);CHKERRQ(ierr);
    ierr = VecPointwiseMult(y,y,st->D);CHKERRQ(ierr);
  } else {
    ierr = (*st->ops->apply)(st,x,y);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(ST_Apply,st,x,y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMultTranspose_STOperator(Mat Op,Vec x,Vec y)
{
  PetscErrorCode ierr;
  ST             st;

  PetscFunctionBegin;
  ierr = MatShellGetContext(Op,(void**)&st);CHKERRQ(ierr);
  ierr = STSetUp(st);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(ST_ApplyTranspose,st,x,y,0);CHKERRQ(ierr);
  if (st->D) { /* with balancing */
    ierr = VecPointwiseMult(st->wb,x,st->D);CHKERRQ(ierr);
    ierr = (*st->ops->applytrans)(st,st->wb,y);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(y,y,st->D);CHKERRQ(ierr);
  } else {
    ierr = (*st->ops->applytrans)(st,x,y);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(ST_ApplyTranspose,st,x,y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_USE_COMPLEX)
static PetscErrorCode MatMultHermitianTranspose_STOperator(Mat Op,Vec x,Vec y)
{
  PetscErrorCode ierr;
  ST             st;

  PetscFunctionBegin;
  ierr = MatShellGetContext(Op,(void**)&st);CHKERRQ(ierr);
  ierr = STSetUp(st);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(ST_ApplyTranspose,st,x,y,0);CHKERRQ(ierr);
  if (!st->wht) {
    ierr = MatCreateVecs(st->A[0],&st->wht,NULL);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)st,(PetscObject)st->wht);CHKERRQ(ierr);
  }
  ierr = VecCopy(x,st->wht);CHKERRQ(ierr);
  ierr = VecConjugate(st->wht);CHKERRQ(ierr);
  if (st->D) { /* with balancing */
    ierr = VecPointwiseMult(st->wb,st->wht,st->D);CHKERRQ(ierr);
    ierr = (*st->ops->applytrans)(st,st->wb,y);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(y,y,st->D);CHKERRQ(ierr);
  } else {
    ierr = (*st->ops->applytrans)(st,st->wht,y);CHKERRQ(ierr);
  }
  ierr = VecConjugate(y);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(ST_ApplyTranspose,st,x,y,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

PetscErrorCode STGetOperator_Private(ST st,Mat *Op)
{
  PetscErrorCode ierr;
  PetscInt       m,n,M,N;

  PetscFunctionBegin;
  if (!st->Op) {
    if (Op) *Op = NULL;
    ierr = MatGetLocalSize(st->A[0],&m,&n);CHKERRQ(ierr);
    ierr = MatGetSize(st->A[0],&M,&N);CHKERRQ(ierr);
    ierr = MatCreateShell(PetscObjectComm((PetscObject)st),m,n,M,N,st,&st->Op);CHKERRQ(ierr);
    ierr = MatShellSetOperation(st->Op,MATOP_MULT,(void(*)(void))MatMult_STOperator);CHKERRQ(ierr);
    ierr = MatShellSetOperation(st->Op,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMultTranspose_STOperator);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = MatShellSetOperation(st->Op,MATOP_MULT_HERMITIAN_TRANSPOSE,(void(*)(void))MatMultHermitianTranspose_STOperator);CHKERRQ(ierr);
#else
    ierr = MatShellSetOperation(st->Op,MATOP_MULT_HERMITIAN_TRANSPOSE,(void(*)(void))MatMultTranspose_STOperator);CHKERRQ(ierr);
#endif
    ierr = STComputeOperator(st);CHKERRQ(ierr);
  }
  if (Op) *Op = st->Op;
  PetscFunctionReturn(0);
}

/*@
   STGetOperator - Returns a shell matrix that represents the operator of the
   spectral transformation.

   Collective on st

   Input Parameter:
.  st - the spectral transformation context

   Output Parameter:
.  Op - operator matrix

   Notes:
   The operator is defined in linear eigenproblems only, not in polynomial ones,
   so the call will fail if more than 2 matrices were passed in STSetMatrices().

   The returned shell matrix is essentially a wrapper to the STApply() and
   STApplyTranspose() operations. The operator can often be expressed as

.vb
      Op = D*inv(K)*M*inv(D)
.ve

   where D is the balancing matrix, and M and K are two matrices corresponding
   to the numerator and denominator for spectral transformations that represent
   a rational matrix function. In the case of STSHELL, the inner part inv(K)*M
   is replaced by the user-provided operation from STShellSetApply().

   The preconditioner matrix K typically depends on the value of the shift, and
   its inverse is handled via an internal KSP object. Normal usage does not
   require explicitly calling STGetOperator(), but it can be used to force the
   creation of K and M, and then K is passed to the KSP. This is useful for
   setting options associated with the PCFactor (to set MUMPS options, for instance).

   The returned matrix must NOT be destroyed by the user. Instead, when no
   longer needed it must be returned with STRestoreOperator(). In particular,
   this is required before modifying the ST matrices or the shift.

   A NULL pointer can be passed in Op in case the matrix is not required but we
   want to force its creation. In this case, STRestoreOperator() should not be
   called.

   Level: advanced

.seealso: STApply(), STApplyTranspose(), STSetBalanceMatrix(), STShellSetApply(),
          STGetKSP(), STSetShift(), STRestoreOperator(), STSetMatrices()
@*/
PetscErrorCode STGetOperator(ST st,Mat *Op)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  STCheckNotSeized(st,1);
  if (st->nmat>2) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONGSTATE,"The operator is not defined in polynomial eigenproblems");
  ierr = STGetOperator_Private(st,Op);CHKERRQ(ierr);
  if (Op) st->opseized = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*@
   STRestoreOperator - Restore the previously seized operator matrix.

   Collective on st

   Input Parameters:
+  st - the spectral transformation context
-  Op - operator matrix

   Notes:
   The arguments must match the corresponding call to STGetOperator().

   Level: advanced

.seealso: STGetOperator()
@*/
PetscErrorCode STRestoreOperator(ST st,Mat *Op)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidPointer(Op,2);
  PetscValidHeaderSpecific(*Op,MAT_CLASSID,2);
  if (!st->opseized) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONGSTATE,"Must be called after STGetOperator()");
  *Op = NULL;
  st->opseized = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*
   STComputeOperator - Computes the matrices that constitute the operator

      Op = D*inv(K)*M*inv(D).

   K and M are computed here (D is user-provided) from the system matrices
   and the shift sigma (whenever these are changed, this function recomputes
   K and M). This is used only in linear eigenproblems (nmat<3).

   K is the "preconditioner matrix": it is the denominator in rational operators,
   e.g. (A-sigma*B) in shift-and-invert. In non-rational transformations such
   as STFILTER, K=NULL which means identity. After computing K, it is passed to
   the internal KSP object via KSPSetOperators.

   M is the numerator in rational operators. If unused it is set to NULL (e.g.
   in STPRECOND).

   STSHELL does not compute anything here, but sets the flag as if it was ready.
*/
PetscErrorCode STComputeOperator(ST st)
{
  PetscErrorCode ierr;
  PC             pc;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  if (!st->opready && st->ops->computeoperator) {
    STCheckMatrices(st,1);
    if (!st->T) {
      ierr = PetscCalloc1(PetscMax(2,st->nmat),&st->T);CHKERRQ(ierr);
      ierr = PetscLogObjectMemory((PetscObject)st,PetscMax(2,st->nmat)*sizeof(Mat));CHKERRQ(ierr);
    }
    ierr = PetscLogEventBegin(ST_ComputeOperator,st,0,0,0);CHKERRQ(ierr);
    ierr = (*st->ops->computeoperator)(st);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(ST_ComputeOperator,st,0,0,0);CHKERRQ(ierr);
    if (st->usesksp) {
      if (!st->ksp) { ierr = STGetKSP(st,&st->ksp);CHKERRQ(ierr); }
      if (st->P) {
        ierr = STSetDefaultKSP(st);CHKERRQ(ierr);
        ierr = STCheckFactorPackage(st);CHKERRQ(ierr);
        ierr = KSPSetOperators(st->ksp,st->P,st->P);CHKERRQ(ierr);
      } else {
        /* STPRECOND defaults to PCNONE if st->P is empty */
        ierr = KSPGetPC(st->ksp,&pc);CHKERRQ(ierr);
        ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
      }
    }
  }
  st->opready = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*@
   STSetUp - Prepares for the use of a spectral transformation.

   Collective on st

   Input Parameter:
.  st - the spectral transformation context

   Level: advanced

.seealso: STCreate(), STApply(), STDestroy()
@*/
PetscErrorCode STSetUp(ST st)
{
  PetscInt       i,n,k;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  STCheckMatrices(st,1);
  switch (st->state) {
    case ST_STATE_INITIAL:
      ierr = PetscInfo(st,"Setting up new ST\n");CHKERRQ(ierr);
      if (!((PetscObject)st)->type_name) {
        ierr = STSetType(st,STSHIFT);CHKERRQ(ierr);
      }
      break;
    case ST_STATE_SETUP:
      PetscFunctionReturn(0);
    case ST_STATE_UPDATED:
      ierr = PetscInfo(st,"Setting up updated ST\n");CHKERRQ(ierr);
      break;
  }
  ierr = PetscLogEventBegin(ST_SetUp,st,0,0,0);CHKERRQ(ierr);
  if (st->state!=ST_STATE_UPDATED) {
    if (!(st->nmat<3 && st->opready)) {
      if (st->T) {
        for (i=0;i<PetscMax(2,st->nmat);i++) {
          ierr = MatDestroy(&st->T[i]);CHKERRQ(ierr);
        }
      }
      ierr = MatDestroy(&st->P);CHKERRQ(ierr);
    }
  }
  if (st->D) {
    ierr = MatGetLocalSize(st->A[0],NULL,&n);CHKERRQ(ierr);
    ierr = VecGetLocalSize(st->D,&k);CHKERRQ(ierr);
    if (n != k) SETERRQ2(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_SIZ,"Balance matrix has wrong dimension %D (should be %D)",k,n);
    if (!st->wb) {
      ierr = VecDuplicate(st->D,&st->wb);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)st,(PetscObject)st->wb);CHKERRQ(ierr);
    }
  }
  if (st->nmat<3 && st->transform) {
    ierr = STComputeOperator(st);CHKERRQ(ierr);
  } else {
    if (!st->T) {
      ierr = PetscCalloc1(PetscMax(2,st->nmat),&st->T);CHKERRQ(ierr);
      ierr = PetscLogObjectMemory((PetscObject)st,PetscMax(2,st->nmat)*sizeof(Mat));CHKERRQ(ierr);
    }
  }
  if (st->ops->setup) { ierr = (*st->ops->setup)(st);CHKERRQ(ierr); }
  st->state = ST_STATE_SETUP;
  ierr = PetscLogEventEnd(ST_SetUp,st,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Computes coefficients for the transformed polynomial,
   and stores the result in argument S.

   alpha - value of the parameter of the transformed polynomial
   beta - value of the previous shift (only used in inplace mode)
   k - index of first matrix included in the computation
   coeffs - coefficients of the expansion
   initial - true if this is the first time (only relevant for shell mode)
*/
PetscErrorCode STMatMAXPY_Private(ST st,PetscScalar alpha,PetscScalar beta,PetscInt k,PetscScalar *coeffs,PetscBool initial,Mat *S)
{
  PetscErrorCode ierr;
  PetscInt       *matIdx=NULL,nmat,i,ini=-1;
  PetscScalar    t=1.0,ta,gamma;
  PetscBool      nz=PETSC_FALSE;

  PetscFunctionBegin;
  nmat = st->nmat-k;
  switch (st->matmode) {
  case ST_MATMODE_INPLACE:
    if (st->nmat>2) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"ST_MATMODE_INPLACE not supported for polynomial eigenproblems");
    if (initial) {
      ierr = PetscObjectReference((PetscObject)st->A[0]);CHKERRQ(ierr);
      *S = st->A[0];
      gamma = alpha;
    } else gamma = alpha-beta;
    if (gamma != 0.0) {
      if (st->nmat>1) {
        ierr = MatAXPY(*S,gamma,st->A[1],st->str);CHKERRQ(ierr);
      } else {
        ierr = MatShift(*S,gamma);CHKERRQ(ierr);
      }
    }
    break;
  case ST_MATMODE_SHELL:
    if (initial) {
      if (st->nmat>2) {
        ierr = PetscMalloc1(nmat,&matIdx);CHKERRQ(ierr);
        for (i=0;i<nmat;i++) matIdx[i] = k+i;
      }
      ierr = STMatShellCreate(st,alpha,nmat,matIdx,coeffs,S);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)st,(PetscObject)*S);CHKERRQ(ierr);
      if (st->nmat>2) { ierr = PetscFree(matIdx);CHKERRQ(ierr); }
    } else {
      ierr = STMatShellShift(*S,alpha);CHKERRQ(ierr);
    }
    break;
  case ST_MATMODE_COPY:
    if (coeffs) {
      for (i=0;i<nmat && ini==-1;i++) {
        if (coeffs[i]!=0.0) ini = i;
        else t *= alpha;
      }
      if (coeffs[ini] != 1.0) nz = PETSC_TRUE;
      for (i=ini+1;i<nmat&&!nz;i++) if (coeffs[i]!=0.0) nz = PETSC_TRUE;
    } else { nz = PETSC_TRUE; ini = 0; }
    if ((alpha == 0.0 || !nz) && t==1.0) {
      ierr = PetscObjectReference((PetscObject)st->A[k+ini]);CHKERRQ(ierr);
      ierr = MatDestroy(S);CHKERRQ(ierr);
      *S = st->A[k+ini];
    } else {
      if (*S && *S!=st->A[k+ini]) {
        ierr = MatSetOption(*S,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
        ierr = MatCopy(st->A[k+ini],*S,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
      } else {
        ierr = MatDestroy(S);CHKERRQ(ierr);
        ierr = MatDuplicate(st->A[k+ini],MAT_COPY_VALUES,S);CHKERRQ(ierr);
        ierr = MatSetOption(*S,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
        ierr = PetscLogObjectParent((PetscObject)st,(PetscObject)*S);CHKERRQ(ierr);
      }
      if (coeffs && coeffs[ini]!=1.0) {
        ierr = MatScale(*S,coeffs[ini]);CHKERRQ(ierr);
      }
      for (i=ini+k+1;i<PetscMax(2,st->nmat);i++) {
        t *= alpha;
        ta = t;
        if (coeffs) ta *= coeffs[i-k];
        if (ta!=0.0) {
          if (st->nmat>1) {
            ierr = MatAXPY(*S,ta,st->A[i],st->str);CHKERRQ(ierr);
          } else {
            ierr = MatShift(*S,ta);CHKERRQ(ierr);
          }
        }
      }
    }
  }
  ierr = MatSetOption(*S,MAT_SYMMETRIC,st->asymm);CHKERRQ(ierr);
  ierr = MatSetOption(*S,MAT_HERMITIAN,(PetscImaginaryPart(st->sigma)==0.0)?st->aherm:PETSC_FALSE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Computes the values of the coefficients required by STMatMAXPY_Private
   for the case of monomial basis.
*/
PetscErrorCode STCoeffs_Monomial(ST st, PetscScalar *coeffs)
{
  PetscInt  k,i,ini,inip;

  PetscFunctionBegin;
  /* Compute binomial coefficients */
  ini = (st->nmat*(st->nmat-1))/2;
  for (i=0;i<st->nmat;i++) coeffs[ini+i]=1.0;
  for (k=st->nmat-1;k>=1;k--) {
    inip = ini+1;
    ini = (k*(k-1))/2;
    coeffs[ini] = 1.0;
    for (i=1;i<k;i++) coeffs[ini+i] = coeffs[ini+i-1]+coeffs[inip+i-1];
  }
  PetscFunctionReturn(0);
}

/*@
   STPostSolve - Optional post-solve phase, intended for any actions that must
   be performed on the ST object after the eigensolver has finished.

   Collective on st

   Input Parameters:
.  st  - the spectral transformation context

   Level: developer

.seealso: EPSSolve()
@*/
PetscErrorCode STPostSolve(ST st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  if (st->ops->postsolve) {
    ierr = (*st->ops->postsolve)(st);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   STBackTransform - Back-transformation phase, intended for
   spectral transformations which require to transform the computed
   eigenvalues back to the original eigenvalue problem.

   Not Collective

   Input Parameters:
+  st   - the spectral transformation context
   eigr - real part of a computed eigenvalue
-  eigi - imaginary part of a computed eigenvalue

   Level: developer

.seealso: STIsInjective()
@*/
PetscErrorCode STBackTransform(ST st,PetscInt n,PetscScalar* eigr,PetscScalar* eigi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  if (st->ops->backtransform) {
    ierr = (*st->ops->backtransform)(st,n,eigr,eigi);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   STIsInjective - Ask if this spectral transformation is injective or not
   (that is, if it corresponds to a one-to-one mapping). If not, then it
   does not make sense to call STBackTransform().

   Not collective

   Input Parameter:
.  st   - the spectral transformation context

   Output Parameter:
.  is - the answer

   Level: developer

.seealso: STBackTransform()
@*/
PetscErrorCode STIsInjective(ST st,PetscBool* is)
{
  PetscErrorCode ierr;
  PetscBool      shell;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidType(st,1);
  PetscValidBoolPointer(is,2);

  ierr = PetscObjectTypeCompare((PetscObject)st,STSHELL,&shell);CHKERRQ(ierr);
  if (shell) {
    ierr = STIsInjective_Shell(st,is);CHKERRQ(ierr);
  } else *is = st->ops->backtransform? PETSC_TRUE: PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@
   STMatSetUp - Build the preconditioner matrix used in STMatSolve().

   Collective on st

   Input Parameters:
+  st     - the spectral transformation context
.  sigma  - the shift
-  coeffs - the coefficients (may be NULL)

   Note:
   This function is not intended to be called by end users, but by SLEPc
   solvers that use ST. It builds matrix st->P as follows, then calls KSPSetUp().
.vb
    If (coeffs):  st->P = Sum_{i=0:nmat-1} coeffs[i]*sigma^i*A_i.
    else          st->P = Sum_{i=0:nmat-1} sigma^i*A_i
.ve

   Level: developer

.seealso: STMatSolve()
@*/
PetscErrorCode STMatSetUp(ST st,PetscScalar sigma,PetscScalar *coeffs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveScalar(st,sigma,2);
  STCheckMatrices(st,1);

  ierr = PetscLogEventBegin(ST_MatSetUp,st,0,0,0);CHKERRQ(ierr);
  ierr = STMatMAXPY_Private(st,sigma,0.0,0,coeffs,PETSC_TRUE,&st->P);CHKERRQ(ierr);
  if (!st->ksp) { ierr = STGetKSP(st,&st->ksp);CHKERRQ(ierr); }
  ierr = STCheckFactorPackage(st);CHKERRQ(ierr);
  ierr = KSPSetOperators(st->ksp,st->P,st->P);CHKERRQ(ierr);
  ierr = KSPSetUp(st->ksp);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(ST_MatSetUp,st,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   STSetWorkVecs - Sets a number of work vectors into the ST object.

   Collective on st

   Input Parameters:
+  st - the spectral transformation context
-  nw - number of work vectors to allocate

   Developers Note:
   This is SLEPC_EXTERN because it may be required by shell STs.

   Level: developer
@*/
PetscErrorCode STSetWorkVecs(ST st,PetscInt nw)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,nw,2);
  if (nw <= 0) SETERRQ1(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"nw must be > 0: nw = %D",nw);
  if (st->nwork < nw) {
    ierr = VecDestroyVecs(st->nwork,&st->work);CHKERRQ(ierr);
    st->nwork = nw;
    ierr = PetscMalloc1(nw,&st->work);CHKERRQ(ierr);
    for (i=0;i<nw;i++) { ierr = STMatCreateVecs(st,&st->work[i],NULL);CHKERRQ(ierr); }
    ierr = PetscLogObjectParents(st,nw,st->work);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

