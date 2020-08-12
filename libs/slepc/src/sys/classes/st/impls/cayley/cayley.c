/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Implements the Cayley spectral transform
*/

#include <slepc/private/stimpl.h>          /*I "slepcst.h" I*/

typedef struct {
  PetscScalar nu;
  PetscBool   nu_set;
} ST_CAYLEY;

static PetscErrorCode MatMult_Cayley(Mat B,Vec x,Vec y)
{
  PetscErrorCode ierr;
  ST             st;
  ST_CAYLEY      *ctx;
  PetscScalar    nu;

  PetscFunctionBegin;
  ierr = MatShellGetContext(B,(void**)&st);CHKERRQ(ierr);
  ctx = (ST_CAYLEY*)st->data;
  nu = ctx->nu;

  if (st->matmode == ST_MATMODE_INPLACE) { nu = nu + st->sigma; };

  if (st->nmat>1) {
    /* generalized eigenproblem: y = (A + tB)x */
    ierr = MatMult(st->A[0],x,y);CHKERRQ(ierr);
    ierr = MatMult(st->A[1],x,st->work[1]);CHKERRQ(ierr);
    ierr = VecAXPY(y,nu,st->work[1]);CHKERRQ(ierr);
  } else {
    /* standard eigenproblem: y = (A + tI)x */
    ierr = MatMult(st->A[0],x,y);CHKERRQ(ierr);
    ierr = VecAXPY(y,nu,x);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMultTranspose_Cayley(Mat B,Vec x,Vec y)
{
  PetscErrorCode ierr;
  ST             st;
  ST_CAYLEY      *ctx;
  PetscScalar    nu;

  PetscFunctionBegin;
  ierr = MatShellGetContext(B,(void**)&st);CHKERRQ(ierr);
  ctx = (ST_CAYLEY*)st->data;
  nu = ctx->nu;

  if (st->matmode == ST_MATMODE_INPLACE) { nu = nu + st->sigma; };
  nu = PetscConj(nu);

  if (st->nmat>1) {
    /* generalized eigenproblem: y = (A + tB)x */
    ierr = MatMultTranspose(st->A[0],x,y);CHKERRQ(ierr);
    ierr = MatMultTranspose(st->A[1],x,st->work[1]);CHKERRQ(ierr);
    ierr = VecAXPY(y,nu,st->work[1]);CHKERRQ(ierr);
  } else {
    /* standard eigenproblem: y = (A + tI)x */
    ierr = MatMultTranspose(st->A[0],x,y);CHKERRQ(ierr);
    ierr = VecAXPY(y,nu,x);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode STGetBilinearForm_Cayley(ST st,Mat *B)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = STSetUp(st);CHKERRQ(ierr);
  *B = st->T[0];
  ierr = PetscObjectReference((PetscObject)*B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STBackTransform_Cayley(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  ST_CAYLEY   *ctx = (ST_CAYLEY*)st->data;
  PetscInt    j;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar t,i,r;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  for (j=0;j<n;j++) {
    if (eigi[j] == 0.0) eigr[j] = (ctx->nu + eigr[j] * st->sigma) / (eigr[j] - 1.0);
    else {
      r = eigr[j];
      i = eigi[j];
      r = st->sigma * (r * r + i * i - r) + ctx->nu * (r - 1);
      i = - st->sigma * i - ctx->nu * i;
      t = i * i + r * (r - 2.0) + 1.0;
      eigr[j] = r / t;
      eigi[j] = i / t;
    }
  }
#else
  for (j=0;j<n;j++) {
    eigr[j] = (ctx->nu + eigr[j] * st->sigma) / (eigr[j] - 1.0);
  }
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode STPostSolve_Cayley(ST st)
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
   Operator (cayley):
               Op                  P         M
   if nmat=1:  (A-sI)^-1 (A+tI)    A-sI      A+tI
   if nmat=2:  (A-sB)^-1 (A+tB)    A-sB      A+tI
*/
PetscErrorCode STComputeOperator_Cayley(ST st)
{
  PetscErrorCode ierr;
  PetscInt       n,m;
  ST_CAYLEY      *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  /* if the user did not set the shift, use the target value */
  if (!st->sigma_set) st->sigma = st->defsigma;

  if (!ctx->nu_set) ctx->nu = st->sigma;
  if (ctx->nu == 0.0 && st->sigma == 0.0) SETERRQ(PetscObjectComm((PetscObject)st),1,"Values of shift and antishift cannot be zero simultaneously");
  if (ctx->nu == -st->sigma) SETERRQ(PetscObjectComm((PetscObject)st),1,"It is not allowed to set the antishift equal to minus the shift (the target)");

  /* T[0] = A+nu*B */
  if (st->matmode==ST_MATMODE_INPLACE) {
    ierr = MatGetLocalSize(st->A[0],&n,&m);CHKERRQ(ierr);
    ierr = MatCreateShell(PetscObjectComm((PetscObject)st),n,m,PETSC_DETERMINE,PETSC_DETERMINE,st,&st->T[0]);CHKERRQ(ierr);
    ierr = MatShellSetOperation(st->T[0],MATOP_MULT,(void(*)(void))MatMult_Cayley);CHKERRQ(ierr);
    ierr = MatShellSetOperation(st->T[0],MATOP_MULT_TRANSPOSE,(void(*)(void))MatMultTranspose_Cayley);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)st,(PetscObject)st->T[0]);CHKERRQ(ierr);
  } else {
    ierr = STMatMAXPY_Private(st,ctx->nu,0.0,0,NULL,PetscNot(st->state==ST_STATE_UPDATED),&st->T[0]);CHKERRQ(ierr);
  }
  st->M = st->T[0];

  /* T[1] = A-sigma*B */
  ierr = STMatMAXPY_Private(st,-st->sigma,0.0,0,NULL,PetscNot(st->state==ST_STATE_UPDATED),&st->T[1]);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)st->T[1]);CHKERRQ(ierr);
  ierr = MatDestroy(&st->P);CHKERRQ(ierr);
  st->P = st->T[1];
  PetscFunctionReturn(0);
}

PetscErrorCode STSetUp_Cayley(ST st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (st->nmat>2) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_SUP,"Cayley transform cannot be used in polynomial eigenproblems");
  ierr = STSetWorkVecs(st,2);CHKERRQ(ierr);
  ierr = KSPSetUp(st->ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STSetShift_Cayley(ST st,PetscScalar newshift)
{
  PetscErrorCode ierr;
  ST_CAYLEY      *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  if (newshift==0.0 && (!ctx->nu_set || (ctx->nu_set && ctx->nu==0.0))) SETERRQ(PetscObjectComm((PetscObject)st),1,"Values of shift and antishift cannot be zero simultaneously");
  if (ctx->nu == -newshift) SETERRQ(PetscObjectComm((PetscObject)st),1,"It is not allowed to set the shift equal to minus the antishift");

  if (!ctx->nu_set) {
    if (st->matmode!=ST_MATMODE_INPLACE) {
      ierr = STMatMAXPY_Private(st,newshift,ctx->nu,0,NULL,PETSC_FALSE,&st->T[0]);CHKERRQ(ierr);
    }
    ctx->nu = newshift;
  }
  ierr = STMatMAXPY_Private(st,-newshift,-st->sigma,0,NULL,PETSC_FALSE,&st->T[1]);CHKERRQ(ierr);
  if (st->P!=st->T[1]) {
    ierr = PetscObjectReference((PetscObject)st->T[1]);CHKERRQ(ierr);
    ierr = MatDestroy(&st->P);CHKERRQ(ierr);
    st->P = st->T[1];
  }
  ierr = KSPSetOperators(st->ksp,st->P,st->P);CHKERRQ(ierr);
  ierr = KSPSetUp(st->ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STSetFromOptions_Cayley(PetscOptionItems *PetscOptionsObject,ST st)
{
  PetscErrorCode ierr;
  PetscScalar    nu;
  PetscBool      flg;
  ST_CAYLEY      *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"ST Cayley Options");CHKERRQ(ierr);

    ierr = PetscOptionsScalar("-st_cayley_antishift","Value of the antishift","STCayleySetAntishift",ctx->nu,&nu,&flg);CHKERRQ(ierr);
    if (flg) { ierr = STCayleySetAntishift(st,nu);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STCayleySetAntishift_Cayley(ST st,PetscScalar newshift)
{
  PetscErrorCode ierr;
  ST_CAYLEY *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  if (ctx->nu != newshift) {
    STCheckNotSeized(st,1);
    if (st->state && st->matmode!=ST_MATMODE_INPLACE) {
      ierr = STMatMAXPY_Private(st,newshift,ctx->nu,0,NULL,PETSC_FALSE,&st->T[0]);CHKERRQ(ierr);
    }
    ctx->nu = newshift;
  }
  ctx->nu_set = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*@
   STCayleySetAntishift - Sets the value of the anti-shift for the Cayley
   spectral transformation.

   Logically Collective on st

   Input Parameters:
+  st  - the spectral transformation context
-  nu  - the anti-shift

   Options Database Key:
.  -st_cayley_antishift - Sets the value of the anti-shift

   Level: intermediate

   Note:
   In the generalized Cayley transform, the operator can be expressed as
   OP = inv(A - sigma B)*(A + nu B). This function sets the value of nu.
   Use STSetShift() for setting sigma. The value nu=-sigma is not allowed.

.seealso: STSetShift(), STCayleyGetAntishift()
@*/
PetscErrorCode STCayleySetAntishift(ST st,PetscScalar nu)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveScalar(st,nu,2);
  ierr = PetscTryMethod(st,"STCayleySetAntishift_C",(ST,PetscScalar),(st,nu));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STCayleyGetAntishift_Cayley(ST st,PetscScalar *nu)
{
  ST_CAYLEY *ctx = (ST_CAYLEY*)st->data;

  PetscFunctionBegin;
  *nu = ctx->nu;
  PetscFunctionReturn(0);
}

/*@
   STCayleyGetAntishift - Gets the value of the anti-shift used in the Cayley
   spectral transformation.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
.  nu  - the anti-shift

   Level: intermediate

.seealso: STGetShift(), STCayleySetAntishift()
@*/
PetscErrorCode STCayleyGetAntishift(ST st,PetscScalar *nu)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidScalarPointer(nu,2);
  ierr = PetscUseMethod(st,"STCayleyGetAntishift_C",(ST,PetscScalar*),(st,nu));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STView_Cayley(ST st,PetscViewer viewer)
{
  PetscErrorCode ierr;
  char           str[50];
  ST_CAYLEY      *ctx = (ST_CAYLEY*)st->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = SlepcSNPrintfScalar(str,50,ctx->nu,PETSC_FALSE);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  antishift: %s\n",str);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode STDestroy_Cayley(ST st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(st->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STCayleySetAntishift_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STCayleyGetAntishift_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode STCreate_Cayley(ST st)
{
  PetscErrorCode ierr;
  ST_CAYLEY      *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(st,&ctx);CHKERRQ(ierr);
  st->data = (void*)ctx;

  st->usesksp = PETSC_TRUE;

  st->ops->apply           = STApply_Generic;
  st->ops->applytrans      = STApplyTranspose_Generic;
  st->ops->backtransform   = STBackTransform_Cayley;
  st->ops->setshift        = STSetShift_Cayley;
  st->ops->getbilinearform = STGetBilinearForm_Cayley;
  st->ops->setup           = STSetUp_Cayley;
  st->ops->computeoperator = STComputeOperator_Cayley;
  st->ops->setfromoptions  = STSetFromOptions_Cayley;
  st->ops->postsolve       = STPostSolve_Cayley;
  st->ops->destroy         = STDestroy_Cayley;
  st->ops->view            = STView_Cayley;
  st->ops->checknullspace  = STCheckNullSpace_Default;
  st->ops->setdefaultksp   = STSetDefaultKSP_Default;

  ierr = PetscObjectComposeFunction((PetscObject)st,"STCayleySetAntishift_C",STCayleySetAntishift_Cayley);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STCayleyGetAntishift_C",STCayleyGetAntishift_Cayley);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

