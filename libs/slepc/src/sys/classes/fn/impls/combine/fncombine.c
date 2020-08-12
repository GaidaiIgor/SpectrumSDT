/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   A function that is obtained by combining two other functions (either by
   addition, multiplication, division or composition)

      addition:          f(x) = f1(x)+f2(x)
      multiplication:    f(x) = f1(x)*f2(x)
      division:          f(x) = f1(x)/f2(x)      f(A) = f2(A)\f1(A)
      composition:       f(x) = f2(f1(x))
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  FN            f1,f2;    /* functions */
  FNCombineType comb;     /* how the functions are combined */
} FN_COMBINE;

PetscErrorCode FNEvaluateFunction_Combine(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;
  PetscScalar    a,b;

  PetscFunctionBegin;
  ierr = FNEvaluateFunction(ctx->f1,x,&a);CHKERRQ(ierr);
  switch (ctx->comb) {
    case FN_COMBINE_ADD:
      ierr = FNEvaluateFunction(ctx->f2,x,&b);CHKERRQ(ierr);
      *y = a+b;
      break;
    case FN_COMBINE_MULTIPLY:
      ierr = FNEvaluateFunction(ctx->f2,x,&b);CHKERRQ(ierr);
      *y = a*b;
      break;
    case FN_COMBINE_DIVIDE:
      ierr = FNEvaluateFunction(ctx->f2,x,&b);CHKERRQ(ierr);
      if (b==0.0) SETERRQ(PETSC_COMM_SELF,1,"Function not defined in the requested value");
      *y = a/b;
      break;
    case FN_COMBINE_COMPOSE:
      ierr = FNEvaluateFunction(ctx->f2,a,y);CHKERRQ(ierr);
      break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FNEvaluateDerivative_Combine(FN fn,PetscScalar x,PetscScalar *yp)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;
  PetscScalar    a,b,ap,bp;

  PetscFunctionBegin;
  switch (ctx->comb) {
    case FN_COMBINE_ADD:
      ierr = FNEvaluateDerivative(ctx->f1,x,&ap);CHKERRQ(ierr);
      ierr = FNEvaluateDerivative(ctx->f2,x,&bp);CHKERRQ(ierr);
      *yp = ap+bp;
      break;
    case FN_COMBINE_MULTIPLY:
      ierr = FNEvaluateDerivative(ctx->f1,x,&ap);CHKERRQ(ierr);
      ierr = FNEvaluateDerivative(ctx->f2,x,&bp);CHKERRQ(ierr);
      ierr = FNEvaluateFunction(ctx->f1,x,&a);CHKERRQ(ierr);
      ierr = FNEvaluateFunction(ctx->f2,x,&b);CHKERRQ(ierr);
      *yp = ap*b+a*bp;
      break;
    case FN_COMBINE_DIVIDE:
      ierr = FNEvaluateDerivative(ctx->f1,x,&ap);CHKERRQ(ierr);
      ierr = FNEvaluateDerivative(ctx->f2,x,&bp);CHKERRQ(ierr);
      ierr = FNEvaluateFunction(ctx->f1,x,&a);CHKERRQ(ierr);
      ierr = FNEvaluateFunction(ctx->f2,x,&b);CHKERRQ(ierr);
      if (b==0.0) SETERRQ(PETSC_COMM_SELF,1,"Derivative not defined in the requested value");
      *yp = (ap*b-a*bp)/(b*b);
      break;
    case FN_COMBINE_COMPOSE:
      ierr = FNEvaluateFunction(ctx->f1,x,&a);CHKERRQ(ierr);
      ierr = FNEvaluateDerivative(ctx->f1,x,&ap);CHKERRQ(ierr);
      ierr = FNEvaluateDerivative(ctx->f2,a,yp);CHKERRQ(ierr);
      *yp *= ap;
      break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FNEvaluateFunctionMat_Combine(FN fn,Mat A,Mat B)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;
  PetscScalar    *Aa,*Ba,*Wa,*Za,one=1.0,zero=0.0;
  PetscBLASInt   n,ld,ld2,inc=1,*ipiv,info;
  PetscInt       m;
  Mat            W,Z;

  PetscFunctionBegin;
  ierr = FN_AllocateWorkMat(fn,A,&W);CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseGetArray(B,&Ba);CHKERRQ(ierr);
  ierr = MatDenseGetArray(W,&Wa);CHKERRQ(ierr);
  ierr = MatGetSize(A,&m,NULL);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ld  = n;
  ld2 = ld*ld;

  switch (ctx->comb) {
    case FN_COMBINE_ADD:
      ierr = FNEvaluateFunctionMat_Private(ctx->f1,A,W,PETSC_FALSE);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat_Private(ctx->f2,A,B,PETSC_FALSE);CHKERRQ(ierr);
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&ld2,&one,Wa,&inc,Ba,&inc));
      ierr = PetscLogFlops(1.0*n*n);CHKERRQ(ierr);
      break;
    case FN_COMBINE_MULTIPLY:
      ierr = FN_AllocateWorkMat(fn,A,&Z);CHKERRQ(ierr);
      ierr = MatDenseGetArray(Z,&Za);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat_Private(ctx->f1,A,W,PETSC_FALSE);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat_Private(ctx->f2,A,Z,PETSC_FALSE);CHKERRQ(ierr);
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&one,Wa,&ld,Za,&ld,&zero,Ba,&ld));
      ierr = PetscLogFlops(2.0*n*n*n);CHKERRQ(ierr);
      ierr = MatDenseRestoreArray(Z,&Za);CHKERRQ(ierr);
      ierr = FN_FreeWorkMat(fn,&Z);CHKERRQ(ierr);
      break;
    case FN_COMBINE_DIVIDE:
      ierr = FNEvaluateFunctionMat_Private(ctx->f2,A,W,PETSC_FALSE);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat_Private(ctx->f1,A,B,PETSC_FALSE);CHKERRQ(ierr);
      ierr = PetscMalloc1(ld,&ipiv);CHKERRQ(ierr);
      PetscStackCallBLAS("LAPACKgesv",LAPACKgesv_(&n,&n,Wa,&ld,ipiv,Ba,&ld,&info));
      SlepcCheckLapackInfo("gesv",info);
      ierr = PetscLogFlops(2.0*n*n*n/3.0+2.0*n*n*n);CHKERRQ(ierr);
      ierr = PetscFree(ipiv);CHKERRQ(ierr);
      break;
    case FN_COMBINE_COMPOSE:
      ierr = FNEvaluateFunctionMat_Private(ctx->f1,A,W,PETSC_FALSE);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat_Private(ctx->f2,W,B,PETSC_FALSE);CHKERRQ(ierr);
      break;
  }

  ierr = MatDenseRestoreArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(B,&Ba);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(W,&Wa);CHKERRQ(ierr);
  ierr = FN_FreeWorkMat(fn,&W);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNEvaluateFunctionMatVec_Combine(FN fn,Mat A,Vec v)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;
  PetscScalar    *va,*Za;
  PetscBLASInt   n,ld,*ipiv,info,one=1;
  PetscInt       m;
  Mat            Z;
  Vec            w;

  PetscFunctionBegin;
  ierr = MatGetSize(A,&m,NULL);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ld = n;

  switch (ctx->comb) {
    case FN_COMBINE_ADD:
      ierr = VecDuplicate(v,&w);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMatVec(ctx->f1,A,w);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMatVec(ctx->f2,A,v);CHKERRQ(ierr);
      ierr = VecAXPY(v,1.0,w);CHKERRQ(ierr);
      ierr = VecDestroy(&w);CHKERRQ(ierr);
      break;
    case FN_COMBINE_MULTIPLY:
      ierr = VecDuplicate(v,&w);CHKERRQ(ierr);
      ierr = FN_AllocateWorkMat(fn,A,&Z);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat_Private(ctx->f1,A,Z,PETSC_FALSE);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMatVec_Private(ctx->f2,A,w,PETSC_FALSE);CHKERRQ(ierr);
      ierr = MatMult(Z,w,v);CHKERRQ(ierr);
      ierr = FN_FreeWorkMat(fn,&Z);CHKERRQ(ierr);
      ierr = VecDestroy(&w);CHKERRQ(ierr);
      break;
    case FN_COMBINE_DIVIDE:
      ierr = VecDuplicate(v,&w);CHKERRQ(ierr);
      ierr = FN_AllocateWorkMat(fn,A,&Z);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat_Private(ctx->f2,A,Z,PETSC_FALSE);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMatVec_Private(ctx->f1,A,v,PETSC_FALSE);CHKERRQ(ierr);
      ierr = PetscMalloc1(ld,&ipiv);CHKERRQ(ierr);
      ierr = MatDenseGetArray(Z,&Za);CHKERRQ(ierr);
      ierr = VecGetArray(v,&va);CHKERRQ(ierr);
      PetscStackCallBLAS("LAPACKgesv",LAPACKgesv_(&n,&one,Za,&ld,ipiv,va,&ld,&info));
      SlepcCheckLapackInfo("gesv",info);
      ierr = PetscLogFlops(2.0*n*n*n/3.0+2.0*n*n);CHKERRQ(ierr);
      ierr = VecRestoreArray(v,&va);CHKERRQ(ierr);
      ierr = MatDenseRestoreArray(Z,&Za);CHKERRQ(ierr);
      ierr = PetscFree(ipiv);CHKERRQ(ierr);
      ierr = FN_FreeWorkMat(fn,&Z);CHKERRQ(ierr);
      ierr = VecDestroy(&w);CHKERRQ(ierr);
      break;
    case FN_COMBINE_COMPOSE:
      ierr = FN_AllocateWorkMat(fn,A,&Z);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat_Private(ctx->f1,A,Z,PETSC_FALSE);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMatVec_Private(ctx->f2,Z,v,PETSC_FALSE);CHKERRQ(ierr);
      ierr = FN_FreeWorkMat(fn,&Z);CHKERRQ(ierr);
      break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FNView_Combine(FN fn,PetscViewer viewer)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    switch (ctx->comb) {
      case FN_COMBINE_ADD:
        ierr = PetscViewerASCIIPrintf(viewer,"  Two added functions f1+f2\n");CHKERRQ(ierr);
        break;
      case FN_COMBINE_MULTIPLY:
        ierr = PetscViewerASCIIPrintf(viewer,"  Two multiplied functions f1*f2\n");CHKERRQ(ierr);
        break;
      case FN_COMBINE_DIVIDE:
        ierr = PetscViewerASCIIPrintf(viewer,"  A quotient of two functions f1/f2\n");CHKERRQ(ierr);
        break;
      case FN_COMBINE_COMPOSE:
        ierr = PetscViewerASCIIPrintf(viewer,"  Two composed functions f2(f1(.))\n");CHKERRQ(ierr);
        break;
    }
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = FNView(ctx->f1,viewer);CHKERRQ(ierr);
    ierr = FNView(ctx->f2,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode FNCombineSetChildren_Combine(FN fn,FNCombineType comb,FN f1,FN f2)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;

  PetscFunctionBegin;
  ctx->comb = comb;
  ierr = PetscObjectReference((PetscObject)f1);CHKERRQ(ierr);
  ierr = FNDestroy(&ctx->f1);CHKERRQ(ierr);
  ctx->f1 = f1;
  ierr = PetscLogObjectParent((PetscObject)fn,(PetscObject)ctx->f1);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)f2);CHKERRQ(ierr);
  ierr = FNDestroy(&ctx->f2);CHKERRQ(ierr);
  ctx->f2 = f2;
  ierr = PetscLogObjectParent((PetscObject)fn,(PetscObject)ctx->f2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   FNCombineSetChildren - Sets the two child functions that constitute this
   combined function, and the way they must be combined.

   Logically Collective on fn

   Input Parameters:
+  fn   - the math function context
.  comb - how to combine the functions (addition, multiplication, division or composition)
.  f1   - first function
-  f2   - second function

   Level: intermediate

.seealso: FNCombineGetChildren()
@*/
PetscErrorCode FNCombineSetChildren(FN fn,FNCombineType comb,FN f1,FN f2)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveEnum(fn,comb,2);
  PetscValidHeaderSpecific(f1,FN_CLASSID,3);
  PetscValidHeaderSpecific(f2,FN_CLASSID,4);
  ierr = PetscTryMethod(fn,"FNCombineSetChildren_C",(FN,FNCombineType,FN,FN),(fn,comb,f1,f2));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode FNCombineGetChildren_Combine(FN fn,FNCombineType *comb,FN *f1,FN *f2)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;

  PetscFunctionBegin;
  if (comb) *comb = ctx->comb;
  if (f1) {
    if (!ctx->f1) {
      ierr = FNCreate(PetscObjectComm((PetscObject)fn),&ctx->f1);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)fn,(PetscObject)ctx->f1);CHKERRQ(ierr);
    }
    *f1 = ctx->f1;
  }
  if (f2) {
    if (!ctx->f2) {
      ierr = FNCreate(PetscObjectComm((PetscObject)fn),&ctx->f2);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)fn,(PetscObject)ctx->f2);CHKERRQ(ierr);
    }
    *f2 = ctx->f2;
  }
  PetscFunctionReturn(0);
}

/*@
   FNCombineGetChildren - Gets the two child functions that constitute this
   combined function, and the way they are combined.

   Not Collective

   Input Parameter:
.  fn   - the math function context

   Output Parameters:
+  comb - how to combine the functions (addition, multiplication, division or composition)
.  f1   - first function
-  f2   - second function

   Level: intermediate

.seealso: FNCombineSetChildren()
@*/
PetscErrorCode FNCombineGetChildren(FN fn,FNCombineType *comb,FN *f1,FN *f2)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  ierr = PetscUseMethod(fn,"FNCombineGetChildren_C",(FN,FNCombineType*,FN*,FN*),(fn,comb,f1,f2));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNDuplicate_Combine(FN fn,MPI_Comm comm,FN *newfn)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data,*ctx2 = (FN_COMBINE*)(*newfn)->data;

  PetscFunctionBegin;
  ctx2->comb = ctx->comb;
  ierr = FNDuplicate(ctx->f1,comm,&ctx2->f1);CHKERRQ(ierr);
  ierr = FNDuplicate(ctx->f2,comm,&ctx2->f2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNDestroy_Combine(FN fn)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx = (FN_COMBINE*)fn->data;

  PetscFunctionBegin;
  ierr = FNDestroy(&ctx->f1);CHKERRQ(ierr);
  ierr = FNDestroy(&ctx->f2);CHKERRQ(ierr);
  ierr = PetscFree(fn->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNCombineSetChildren_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNCombineGetChildren_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode FNCreate_Combine(FN fn)
{
  PetscErrorCode ierr;
  FN_COMBINE     *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(fn,&ctx);CHKERRQ(ierr);
  fn->data = (void*)ctx;

  fn->ops->evaluatefunction          = FNEvaluateFunction_Combine;
  fn->ops->evaluatederivative        = FNEvaluateDerivative_Combine;
  fn->ops->evaluatefunctionmat[0]    = FNEvaluateFunctionMat_Combine;
  fn->ops->evaluatefunctionmatvec[0] = FNEvaluateFunctionMatVec_Combine;
  fn->ops->view                      = FNView_Combine;
  fn->ops->duplicate                 = FNDuplicate_Combine;
  fn->ops->destroy                   = FNDestroy_Combine;
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNCombineSetChildren_C",FNCombineSetChildren_Combine);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)fn,"FNCombineGetChildren_C",FNCombineGetChildren_Combine);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

