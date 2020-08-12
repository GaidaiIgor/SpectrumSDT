/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Rational function  r(x) = p(x)/q(x), where p(x) and q(x) are polynomials
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  PetscScalar *pcoeff;    /* numerator coefficients */
  PetscInt    np;         /* length of array pcoeff, p(x) has degree np-1 */
  PetscScalar *qcoeff;    /* denominator coefficients */
  PetscInt    nq;         /* length of array qcoeff, q(x) has degree nq-1 */
} FN_RATIONAL;

PetscErrorCode FNEvaluateFunction_Rational(FN fn,PetscScalar x,PetscScalar *y)
{
  FN_RATIONAL *ctx = (FN_RATIONAL*)fn->data;
  PetscInt    i;
  PetscScalar p,q;

  PetscFunctionBegin;
  if (!ctx->np) p = 1.0;
  else {
    p = ctx->pcoeff[0];
    for (i=1;i<ctx->np;i++)
      p = ctx->pcoeff[i]+x*p;
  }
  if (!ctx->nq) *y = p;
  else {
    q = ctx->qcoeff[0];
    for (i=1;i<ctx->nq;i++)
      q = ctx->qcoeff[i]+x*q;
    if (q==0.0) SETERRQ(PETSC_COMM_SELF,1,"Function not defined in the requested value");
    *y = p/q;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode FNEvaluateFunctionMat_Rational_Private(FN fn,PetscScalar *Aa,PetscScalar *Ba,PetscInt m,PetscBool firstonly)
{
  PetscErrorCode ierr;
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscBLASInt   n,k,ld,*ipiv,info;
  PetscInt       i,j;
  PetscScalar    *W,*P,*Q,one=1.0,zero=0.0;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ld = n;
  k  = firstonly? 1: n;
  if (Aa==Ba) {
    ierr = PetscMalloc4(m*m,&P,m*m,&Q,m*m,&W,ld,&ipiv);CHKERRQ(ierr);
  } else {
    P = Ba;
    ierr = PetscMalloc3(m*m,&Q,m*m,&W,ld,&ipiv);CHKERRQ(ierr);
  }
  ierr = PetscArrayzero(P,m*m);CHKERRQ(ierr);
  if (!ctx->np) {
    for (i=0;i<m;i++) P[i+i*ld] = 1.0;
  } else {
    for (i=0;i<m;i++) P[i+i*ld] = ctx->pcoeff[0];
    for (j=1;j<ctx->np;j++) {
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&one,P,&ld,Aa,&ld,&zero,W,&ld));
      ierr = PetscArraycpy(P,W,m*m);CHKERRQ(ierr);
      for (i=0;i<m;i++) P[i+i*ld] += ctx->pcoeff[j];
    }
    ierr = PetscLogFlops(2.0*n*n*n*(ctx->np-1));CHKERRQ(ierr);
  }
  if (ctx->nq) {
    ierr = PetscArrayzero(Q,m*m);CHKERRQ(ierr);
    for (i=0;i<m;i++) Q[i+i*ld] = ctx->qcoeff[0];
    for (j=1;j<ctx->nq;j++) {
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&one,Q,&ld,Aa,&ld,&zero,W,&ld));
      ierr = PetscArraycpy(Q,W,m*m);CHKERRQ(ierr);
      for (i=0;i<m;i++) Q[i+i*ld] += ctx->qcoeff[j];
    }
    PetscStackCallBLAS("LAPACKgesv",LAPACKgesv_(&n,&k,Q,&ld,ipiv,P,&ld,&info));
    SlepcCheckLapackInfo("gesv",info);
    ierr = PetscLogFlops(2.0*n*n*n*(ctx->nq-1)+2.0*n*n*n/3.0+2.0*n*n*k);CHKERRQ(ierr);
  }
  if (Aa==Ba) {
    ierr = PetscArraycpy(Aa,P,m*k);CHKERRQ(ierr);
    ierr = PetscFree4(P,Q,W,ipiv);CHKERRQ(ierr);
  } else {
    ierr = PetscFree3(Q,W,ipiv);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FNEvaluateFunctionMat_Rational(FN fn,Mat A,Mat B)
{
  PetscErrorCode ierr;
  PetscInt       m;
  PetscScalar    *Aa,*Ba;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseGetArray(B,&Ba);CHKERRQ(ierr);
  ierr = MatGetSize(A,&m,NULL);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMat_Rational_Private(fn,Aa,Ba,m,PETSC_FALSE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(B,&Ba);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNEvaluateFunctionMatVec_Rational(FN fn,Mat A,Vec v)
{
  PetscErrorCode ierr;
  PetscInt       m;
  PetscScalar    *Aa,*Ba;
  Mat            B;

  PetscFunctionBegin;
  ierr = FN_AllocateWorkMat(fn,A,&B);CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseGetArray(B,&Ba);CHKERRQ(ierr);
  ierr = MatGetSize(A,&m,NULL);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMat_Rational_Private(fn,Aa,Ba,m,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(B,&Ba);CHKERRQ(ierr);
  ierr = MatGetColumnVector(B,v,0);CHKERRQ(ierr);
  ierr = FN_FreeWorkMat(fn,&B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNEvaluateDerivative_Rational(FN fn,PetscScalar x,PetscScalar *yp)
{
  FN_RATIONAL *ctx = (FN_RATIONAL*)fn->data;
  PetscInt    i;
  PetscScalar p,q,pp,qp;

  PetscFunctionBegin;
  if (!ctx->np) {
    p = 1.0;
    pp = 0.0;
  } else {
    p = ctx->pcoeff[0];
    pp = 0.0;
    for (i=1;i<ctx->np;i++) {
      pp = p+x*pp;
      p = ctx->pcoeff[i]+x*p;
    }
  }
  if (!ctx->nq) *yp = pp;
  else {
    q = ctx->qcoeff[0];
    qp = 0.0;
    for (i=1;i<ctx->nq;i++) {
      qp = q+x*qp;
      q = ctx->qcoeff[i]+x*q;
    }
    if (q==0.0) SETERRQ(PETSC_COMM_SELF,1,"Derivative not defined in the requested value");
    *yp = (pp*q-p*qp)/(q*q);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FNView_Rational(FN fn,PetscViewer viewer)
{
  PetscErrorCode ierr;
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscBool      isascii;
  PetscInt       i;
  char           str[50];

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    if (fn->alpha!=(PetscScalar)1.0 || fn->beta!=(PetscScalar)1.0) {
      ierr = SlepcSNPrintfScalar(str,50,fn->alpha,PETSC_FALSE);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"  Scale factors: alpha=%s,",str);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
      ierr = SlepcSNPrintfScalar(str,50,fn->beta,PETSC_FALSE);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer," beta=%s\n",str);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    }
    if (!ctx->nq) {
      if (!ctx->np) {
        ierr = PetscViewerASCIIPrintf(viewer,"  Constant: 1.0\n");CHKERRQ(ierr);
      } else if (ctx->np==1) {
        ierr = SlepcSNPrintfScalar(str,50,ctx->pcoeff[0],PETSC_FALSE);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"  Constant: %s\n",str);CHKERRQ(ierr);
      } else {
        ierr = PetscViewerASCIIPrintf(viewer,"  Polynomial: ");CHKERRQ(ierr);
        ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
        for (i=0;i<ctx->np-1;i++) {
          ierr = SlepcSNPrintfScalar(str,50,ctx->pcoeff[i],PETSC_TRUE);CHKERRQ(ierr);
          ierr = PetscViewerASCIIPrintf(viewer,"%s*x^%1D",str,ctx->np-i-1);CHKERRQ(ierr);
        }
        ierr = SlepcSNPrintfScalar(str,50,ctx->pcoeff[ctx->np-1],PETSC_TRUE);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"%s\n",str);CHKERRQ(ierr);
        ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
      }
    } else if (!ctx->np) {
      ierr = PetscViewerASCIIPrintf(viewer,"  Inverse polinomial: 1 / (");CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
      for (i=0;i<ctx->nq-1;i++) {
        ierr = SlepcSNPrintfScalar(str,50,ctx->qcoeff[i],PETSC_TRUE);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"%s*x^%1D",str,ctx->nq-i-1);CHKERRQ(ierr);
      }
      ierr = SlepcSNPrintfScalar(str,50,ctx->qcoeff[ctx->nq-1],PETSC_TRUE);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%s)\n",str);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"  Rational function: (");CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
      for (i=0;i<ctx->np-1;i++) {
        ierr = SlepcSNPrintfScalar(str,50,ctx->pcoeff[i],PETSC_TRUE);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"%s*x^%1D",str,ctx->np-i-1);CHKERRQ(ierr);
      }
      ierr = SlepcSNPrintfScalar(str,50,ctx->pcoeff[ctx->np-1],PETSC_TRUE);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%s) / (",str);CHKERRQ(ierr);
      for (i=0;i<ctx->nq-1;i++) {
        ierr = SlepcSNPrintfScalar(str,50,ctx->qcoeff[i],PETSC_TRUE);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"%s*x^%1D",str,ctx->nq-i-1);CHKERRQ(ierr);
      }
      ierr = SlepcSNPrintfScalar(str,50,ctx->qcoeff[ctx->nq-1],PETSC_TRUE);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%s)\n",str);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode FNRationalSetNumerator_Rational(FN fn,PetscInt np,PetscScalar *pcoeff)
{
  PetscErrorCode ierr;
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (np<0) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"Argument np cannot be negative");
  ctx->np = np;
  ierr = PetscFree(ctx->pcoeff);CHKERRQ(ierr);
  if (np) {
    ierr = PetscMalloc1(np,&ctx->pcoeff);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)fn,np*sizeof(PetscScalar));CHKERRQ(ierr);
    for (i=0;i<np;i++) ctx->pcoeff[i] = pcoeff[i];
  }
  PetscFunctionReturn(0);
}

/*@C
   FNRationalSetNumerator - Sets the parameters defining the numerator of the
   rational function.

   Logically Collective on fn

   Input Parameters:
+  fn     - the math function context
.  np     - number of coefficients
-  pcoeff - coefficients (array of scalar values)

   Notes:
   Let the rational function r(x) = p(x)/q(x), where p(x) and q(x) are polynomials.
   This function provides the coefficients of the numerator p(x).
   Hence, p(x) is of degree np-1.
   If np is zero, then the numerator is assumed to be p(x)=1.

   In polynomials, high order coefficients are stored in the first positions
   of the array, e.g. to represent x^2-3 use {1,0,-3}.

   Level: intermediate

.seealso: FNRationalSetDenominator(), FNRationalGetNumerator()
@*/
PetscErrorCode FNRationalSetNumerator(FN fn,PetscInt np,PetscScalar *pcoeff)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveInt(fn,np,2);
  if (np) PetscValidScalarPointer(pcoeff,3);
  ierr = PetscTryMethod(fn,"FNRationalSetNumerator_C",(FN,PetscInt,PetscScalar*),(fn,np,pcoeff));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode FNRationalGetNumerator_Rational(FN fn,PetscInt *np,PetscScalar *pcoeff[])
{
  PetscErrorCode ierr;
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (np) *np = ctx->np;
  if (pcoeff) {
    if (!ctx->np) *pcoeff = NULL;
    else {
      ierr = PetscMalloc1(ctx->np,pcoeff);CHKERRQ(ierr);
      for (i=0;i<ctx->np;i++) (*pcoeff)[i] = ctx->pcoeff[i];
    }
  }
  PetscFunctionReturn(0);
}

/*@C
   FNRationalGetNumerator - Gets the parameters that define the numerator of the
   rational function.

   Not Collective

   Input Parameter:
.  fn     - the math function context

   Output Parameters:
+  np     - number of coefficients
-  pcoeff - coefficients (array of scalar values, length nq)

   Notes:
   The values passed by user with FNRationalSetNumerator() are returned (or null
   pointers otherwise).
   The pcoeff array should be freed by the user when no longer needed.

   Level: intermediate

.seealso: FNRationalSetNumerator()
@*/
PetscErrorCode FNRationalGetNumerator(FN fn,PetscInt *np,PetscScalar *pcoeff[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  ierr = PetscUseMethod(fn,"FNRationalGetNumerator_C",(FN,PetscInt*,PetscScalar**),(fn,np,pcoeff));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode FNRationalSetDenominator_Rational(FN fn,PetscInt nq,PetscScalar *qcoeff)
{
  PetscErrorCode ierr;
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (nq<0) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"Argument nq cannot be negative");
  ctx->nq = nq;
  ierr = PetscFree(ctx->qcoeff);CHKERRQ(ierr);
  if (nq) {
    ierr = PetscMalloc1(nq,&ctx->qcoeff);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)fn,nq*sizeof(PetscScalar));CHKERRQ(ierr);
    for (i=0;i<nq;i++) ctx->qcoeff[i] = qcoeff[i];
  }
  PetscFunctionReturn(0);
}

/*@C
   FNRationalSetDenominator - Sets the parameters defining the denominator of the
   rational function.

   Logically Collective on fn

   Input Parameters:
+  fn     - the math function context
.  nq     - number of coefficients
-  qcoeff - coefficients (array of scalar values)

   Notes:
   Let the rational function r(x) = p(x)/q(x), where p(x) and q(x) are polynomials.
   This function provides the coefficients of the denominator q(x).
   Hence, q(x) is of degree nq-1.
   If nq is zero, then the function is assumed to be polynomial, r(x) = p(x).

   In polynomials, high order coefficients are stored in the first positions
   of the array, e.g. to represent x^2-3 use {1,0,-3}.

   Level: intermediate

.seealso: FNRationalSetNumerator(), FNRationalGetDenominator()
@*/
PetscErrorCode FNRationalSetDenominator(FN fn,PetscInt nq,PetscScalar *qcoeff)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveInt(fn,nq,2);
  if (nq) PetscValidScalarPointer(qcoeff,3);
  ierr = PetscTryMethod(fn,"FNRationalSetDenominator_C",(FN,PetscInt,PetscScalar*),(fn,nq,qcoeff));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode FNRationalGetDenominator_Rational(FN fn,PetscInt *nq,PetscScalar *qcoeff[])
{
  PetscErrorCode ierr;
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (nq) *nq = ctx->nq;
  if (qcoeff) {
    if (!ctx->nq) *qcoeff = NULL;
    else {
      ierr = PetscMalloc1(ctx->nq,qcoeff);CHKERRQ(ierr);
      for (i=0;i<ctx->nq;i++) (*qcoeff)[i] = ctx->qcoeff[i];
    }
  }
  PetscFunctionReturn(0);
}

/*@C
   FNRationalGetDenominator - Gets the parameters that define the denominator of the
   rational function.

   Not Collective

   Input Parameter:
.  fn     - the math function context

   Output Parameters:
+  nq     - number of coefficients
-  qcoeff - coefficients (array of scalar values, length nq)

   Notes:
   The values passed by user with FNRationalSetDenominator() are returned (or a null
   pointer otherwise).
   The qcoeff array should be freed by the user when no longer needed.

   Level: intermediate

.seealso: FNRationalSetDenominator()
@*/
PetscErrorCode FNRationalGetDenominator(FN fn,PetscInt *nq,PetscScalar *qcoeff[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  ierr = PetscUseMethod(fn,"FNRationalGetDenominator_C",(FN,PetscInt*,PetscScalar**),(fn,nq,qcoeff));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNSetFromOptions_Rational(PetscOptionItems *PetscOptionsObject,FN fn)
{
  PetscErrorCode ierr;
#define PARMAX 10
  PetscScalar    array[PARMAX];
  PetscInt       i,k;
  PetscBool      flg;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"FN Rational Options");CHKERRQ(ierr);

    k = PARMAX;
    for (i=0;i<k;i++) array[i] = 0;
    ierr = PetscOptionsScalarArray("-fn_rational_numerator","Numerator coefficients (one or more scalar values separated with a comma without spaces)","FNRationalSetNumerator",array,&k,&flg);CHKERRQ(ierr);
    if (flg) { ierr = FNRationalSetNumerator(fn,k,array);CHKERRQ(ierr); }

    k = PARMAX;
    for (i=0;i<k;i++) array[i] = 0;
    ierr = PetscOptionsScalarArray("-fn_rational_denominator","Denominator coefficients (one or more scalar values separated with a comma without spaces)","FNRationalSetDenominator",array,&k,&flg);CHKERRQ(ierr);
    if (flg) { ierr = FNRationalSetDenominator(fn,k,array);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNDuplicate_Rational(FN fn,MPI_Comm comm,FN *newfn)
{
  PetscErrorCode ierr;
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data,*ctx2 = (FN_RATIONAL*)(*newfn)->data;
  PetscInt       i;

  PetscFunctionBegin;
  ctx2->np = ctx->np;
  if (ctx->np) {
    ierr = PetscMalloc1(ctx->np,&ctx2->pcoeff);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)(*newfn),ctx->np*sizeof(PetscScalar));CHKERRQ(ierr);
    for (i=0;i<ctx->np;i++) ctx2->pcoeff[i] = ctx->pcoeff[i];
  }
  ctx2->nq = ctx->nq;
  if (ctx->nq) {
    ierr = PetscMalloc1(ctx->nq,&ctx2->qcoeff);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)(*newfn),ctx->nq*sizeof(PetscScalar));CHKERRQ(ierr);
    for (i=0;i<ctx->nq;i++) ctx2->qcoeff[i] = ctx->qcoeff[i];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FNDestroy_Rational(FN fn)
{
  PetscErrorCode ierr;
  FN_RATIONAL    *ctx = (FN_RATIONAL*)fn->data;

  PetscFunctionBegin;
  ierr = PetscFree(ctx->pcoeff);CHKERRQ(ierr);
  ierr = PetscFree(ctx->qcoeff);CHKERRQ(ierr);
  ierr = PetscFree(fn->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNRationalSetNumerator_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNRationalGetNumerator_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNRationalSetDenominator_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNRationalGetDenominator_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode FNCreate_Rational(FN fn)
{
  PetscErrorCode ierr;
  FN_RATIONAL    *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(fn,&ctx);CHKERRQ(ierr);
  fn->data = (void*)ctx;

  fn->ops->evaluatefunction          = FNEvaluateFunction_Rational;
  fn->ops->evaluatederivative        = FNEvaluateDerivative_Rational;
  fn->ops->evaluatefunctionmat[0]    = FNEvaluateFunctionMat_Rational;
  fn->ops->evaluatefunctionmatvec[0] = FNEvaluateFunctionMatVec_Rational;
  fn->ops->setfromoptions            = FNSetFromOptions_Rational;
  fn->ops->view                      = FNView_Rational;
  fn->ops->duplicate                 = FNDuplicate_Rational;
  fn->ops->destroy                   = FNDestroy_Rational;
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNRationalSetNumerator_C",FNRationalSetNumerator_Rational);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNRationalGetNumerator_C",FNRationalGetNumerator_Rational);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNRationalSetDenominator_C",FNRationalSetDenominator_Rational);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNRationalGetDenominator_C",FNRationalGetDenominator_Rational);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

