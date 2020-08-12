/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic FN routines
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/
#include <slepcblaslapack.h>

PetscFunctionList FNList = 0;
PetscBool         FNRegisterAllCalled = PETSC_FALSE;
PetscClassId      FN_CLASSID = 0;
PetscLogEvent     FN_Evaluate = 0;
static PetscBool  FNPackageInitialized = PETSC_FALSE;

const char *FNParallelTypes[] = {"REDUNDANT","SYNCHRONIZED","FNParallelType","FN_PARALLEL_",0};

/*@C
   FNFinalizePackage - This function destroys everything in the Slepc interface
   to the FN package. It is called from SlepcFinalize().

   Level: developer

.seealso: SlepcFinalize()
@*/
PetscErrorCode FNFinalizePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListDestroy(&FNList);CHKERRQ(ierr);
  FNPackageInitialized = PETSC_FALSE;
  FNRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
  FNInitializePackage - This function initializes everything in the FN package.
  It is called from PetscDLLibraryRegister() when using dynamic libraries, and
  on the first call to FNCreate() when using static libraries.

  Level: developer

.seealso: SlepcInitialize()
@*/
PetscErrorCode FNInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (FNPackageInitialized) PetscFunctionReturn(0);
  FNPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  ierr = PetscClassIdRegister("Math Function",&FN_CLASSID);CHKERRQ(ierr);
  /* Register Constructors */
  ierr = FNRegisterAll();CHKERRQ(ierr);
  /* Register Events */
  ierr = PetscLogEventRegister("FNEvaluate",FN_CLASSID,&FN_Evaluate);CHKERRQ(ierr);
  /* Process Info */
  classids[0] = FN_CLASSID;
  ierr = PetscInfoProcessClass("fn",1,&classids[0]);CHKERRQ(ierr);
  /* Process summary exclusions */
  ierr = PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt);CHKERRQ(ierr);
  if (opt) {
    ierr = PetscStrInList("fn",logList,',',&pkg);CHKERRQ(ierr);
    if (pkg) { ierr = PetscLogEventDeactivateClass(FN_CLASSID);CHKERRQ(ierr); }
  }
  /* Register package finalizer */
  ierr = PetscRegisterFinalize(FNFinalizePackage);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   FNCreate - Creates an FN context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  newfn - location to put the FN context

   Level: beginner

.seealso: FNDestroy(), FN
@*/
PetscErrorCode FNCreate(MPI_Comm comm,FN *newfn)
{
  FN             fn;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidPointer(newfn,2);
  *newfn = 0;
  ierr = FNInitializePackage();CHKERRQ(ierr);
  ierr = SlepcHeaderCreate(fn,FN_CLASSID,"FN","Math Function","FN",comm,FNDestroy,FNView);CHKERRQ(ierr);

  fn->alpha    = 1.0;
  fn->beta     = 1.0;
  fn->method   = 0;

  fn->nw       = 0;
  fn->cw       = 0;
  fn->data     = NULL;

  *newfn = fn;
  PetscFunctionReturn(0);
}

/*@C
   FNSetOptionsPrefix - Sets the prefix used for searching for all
   FN options in the database.

   Logically Collective on fn

   Input Parameters:
+  fn - the math function context
-  prefix - the prefix string to prepend to all FN option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: FNAppendOptionsPrefix()
@*/
PetscErrorCode FNSetOptionsPrefix(FN fn,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)fn,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   FNAppendOptionsPrefix - Appends to the prefix used for searching for all
   FN options in the database.

   Logically Collective on fn

   Input Parameters:
+  fn - the math function context
-  prefix - the prefix string to prepend to all FN option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: FNSetOptionsPrefix()
@*/
PetscErrorCode FNAppendOptionsPrefix(FN fn,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  ierr = PetscObjectAppendOptionsPrefix((PetscObject)fn,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   FNGetOptionsPrefix - Gets the prefix used for searching for all
   FN options in the database.

   Not Collective

   Input Parameters:
.  fn - the math function context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: FNSetOptionsPrefix(), FNAppendOptionsPrefix()
@*/
PetscErrorCode FNGetOptionsPrefix(FN fn,const char *prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidPointer(prefix,2);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)fn,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   FNSetType - Selects the type for the FN object.

   Logically Collective on fn

   Input Parameter:
+  fn   - the math function context
-  type - a known type

   Notes:
   The default is FNRATIONAL, which includes polynomials as a particular
   case as well as simple functions such as f(x)=x and f(x)=constant.

   Level: intermediate

.seealso: FNGetType()
@*/
PetscErrorCode FNSetType(FN fn,FNType type)
{
  PetscErrorCode ierr,(*r)(FN);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidCharPointer(type,2);

  ierr = PetscObjectTypeCompare((PetscObject)fn,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);

  ierr =  PetscFunctionListFind(FNList,type,&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested FN type %s",type);

  if (fn->ops->destroy) { ierr = (*fn->ops->destroy)(fn);CHKERRQ(ierr); }
  ierr = PetscMemzero(fn->ops,sizeof(struct _FNOps));CHKERRQ(ierr);

  ierr = PetscObjectChangeTypeName((PetscObject)fn,type);CHKERRQ(ierr);
  ierr = (*r)(fn);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   FNGetType - Gets the FN type name (as a string) from the FN context.

   Not Collective

   Input Parameter:
.  fn - the math function context

   Output Parameter:
.  name - name of the math function

   Level: intermediate

.seealso: FNSetType()
@*/
PetscErrorCode FNGetType(FN fn,FNType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidPointer(type,2);
  *type = ((PetscObject)fn)->type_name;
  PetscFunctionReturn(0);
}

/*@
   FNSetScale - Sets the scaling parameters that define the matematical function.

   Logically Collective on fn

   Input Parameters:
+  fn    - the math function context
.  alpha - inner scaling (argument)
-  beta  - outer scaling (result)

   Notes:
   Given a function f(x) specified by the FN type, the scaling parameters can
   be used to realize the function beta*f(alpha*x). So when these values are given,
   the procedure for function evaluation will first multiply the argument by alpha,
   then evaluate the function itself, and finally scale the result by beta.
   Likewise, these values are also considered when evaluating the derivative.

   If you want to provide only one of the two scaling factors, set the other
   one to 1.0.

   Level: intermediate

.seealso: FNGetScale(), FNEvaluateFunction()
@*/
PetscErrorCode FNSetScale(FN fn,PetscScalar alpha,PetscScalar beta)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveScalar(fn,alpha,2);
  PetscValidLogicalCollectiveScalar(fn,beta,2);
  if (PetscAbsScalar(alpha)==0.0 || PetscAbsScalar(beta)==0.0) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_WRONG,"Scaling factors must be nonzero");
  fn->alpha = alpha;
  fn->beta  = beta;
  PetscFunctionReturn(0);
}

/*@
   FNGetScale - Gets the scaling parameters that define the matematical function.

   Not Collective

   Input Parameter:
.  fn    - the math function context

   Output Parameters:
+  alpha - inner scaling (argument)
-  beta  - outer scaling (result)

   Level: intermediate

.seealso: FNSetScale()
@*/
PetscErrorCode FNGetScale(FN fn,PetscScalar *alpha,PetscScalar *beta)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  if (alpha) *alpha = fn->alpha;
  if (beta)  *beta  = fn->beta;
  PetscFunctionReturn(0);
}

/*@
   FNSetMethod - Selects the method to be used to evaluate functions of matrices.

   Logically Collective on fn

   Input Parameter:
+  fn   - the math function context
-  meth - an index indentifying the method

   Options Database Key:
.  -fn_method <meth> - Sets the method

   Notes:
   In some FN types there are more than one algorithms available for computing
   matrix functions. In that case, this function allows choosing the wanted method.

   If meth is currently set to 0 (the default) and the input argument A of
   FNEvaluateFunctionMat() is a symmetric/Hermitian matrix, then the computation
   is done via the eigendecomposition of A, rather than with the general algorithm.

   Level: intermediate

.seealso: FNGetMethod(), FNEvaluateFunctionMat()
@*/
PetscErrorCode FNSetMethod(FN fn,PetscInt meth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveInt(fn,meth,2);
  if (meth<0) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"The method must be a non-negative integer");
  if (meth>FN_MAX_SOLVE) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"Too large value for the method");
  fn->method = meth;
  PetscFunctionReturn(0);
}

/*@
   FNGetMethod - Gets the method currently used in the FN.

   Not Collective

   Input Parameter:
.  fn - the math function context

   Output Parameter:
.  meth - identifier of the method

   Level: intermediate

.seealso: FNSetMethod()
@*/
PetscErrorCode FNGetMethod(FN fn,PetscInt *meth)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidIntPointer(meth,2);
  *meth = fn->method;
  PetscFunctionReturn(0);
}

/*@
   FNSetParallel - Selects the mode of operation in parallel runs.

   Logically Collective on fn

   Input Parameter:
+  fn    - the math function context
-  pmode - the parallel mode

   Options Database Key:
.  -fn_parallel <mode> - Sets the parallel mode, either 'redundant' or 'synchronized'

   Notes:
   This is relevant only when the function is evaluated on a matrix, with
   either FNEvaluateFunctionMat() or FNEvaluateFunctionMatVec().

   In the 'redundant' parallel mode, all processes will make the computation
   redundantly, starting from the same data, and producing the same result.
   This result may be slightly different in the different processes if using a
   multithreaded BLAS library, which may cause issues in ill-conditioned problems.

   In the 'synchronized' parallel mode, only the first MPI process performs the
   computation and then the computed matrix is broadcast to the other
   processes in the communicator. This communication is done automatically at
   the end of FNEvaluateFunctionMat() or FNEvaluateFunctionMatVec().

   Level: advanced

.seealso: FNEvaluateFunctionMat() or FNEvaluateFunctionMatVec(), FNGetParallel()
@*/
PetscErrorCode FNSetParallel(FN fn,FNParallelType pmode)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidLogicalCollectiveEnum(fn,pmode,2);
  fn->pmode = pmode;
  PetscFunctionReturn(0);
}

/*@
   FNGetParallel - Gets the mode of operation in parallel runs.

   Not Collective

   Input Parameter:
.  fn - the math function context

   Output Parameter:
.  pmode - the parallel mode

   Level: advanced

.seealso: FNSetParallel()
@*/
PetscErrorCode FNGetParallel(FN fn,FNParallelType *pmode)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidPointer(pmode,2);
  *pmode = fn->pmode;
  PetscFunctionReturn(0);
}

/*@
   FNEvaluateFunction - Computes the value of the function f(x) for a given x.

   Not collective

   Input Parameters:
+  fn - the math function context
-  x  - the value where the function must be evaluated

   Output Parameter:
.  y  - the result of f(x)

   Note:
   Scaling factors are taken into account, so the actual function evaluation
   will return beta*f(alpha*x).

   Level: intermediate

.seealso: FNEvaluateDerivative(), FNEvaluateFunctionMat(), FNSetScale()
@*/
PetscErrorCode FNEvaluateFunction(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscErrorCode ierr;
  PetscScalar    xf,yf;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidType(fn,1);
  PetscValidScalarPointer(y,3);
  ierr = PetscLogEventBegin(FN_Evaluate,fn,0,0,0);CHKERRQ(ierr);
  xf = fn->alpha*x;
  ierr = (*fn->ops->evaluatefunction)(fn,xf,&yf);CHKERRQ(ierr);
  *y = fn->beta*yf;
  ierr = PetscLogEventEnd(FN_Evaluate,fn,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   FNEvaluateDerivative - Computes the value of the derivative f'(x) for a given x.

   Not Collective

   Input Parameters:
+  fn - the math function context
-  x  - the value where the derivative must be evaluated

   Output Parameter:
.  y  - the result of f'(x)

   Note:
   Scaling factors are taken into account, so the actual derivative evaluation will
   return alpha*beta*f'(alpha*x).

   Level: intermediate

.seealso: FNEvaluateFunction()
@*/
PetscErrorCode FNEvaluateDerivative(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscErrorCode ierr;
  PetscScalar    xf,yf;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidType(fn,1);
  PetscValidScalarPointer(y,3);
  ierr = PetscLogEventBegin(FN_Evaluate,fn,0,0,0);CHKERRQ(ierr);
  xf = fn->alpha*x;
  ierr = (*fn->ops->evaluatederivative)(fn,xf,&yf);CHKERRQ(ierr);
  *y = fn->alpha*fn->beta*yf;
  ierr = PetscLogEventEnd(FN_Evaluate,fn,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode FNEvaluateFunctionMat_Sym_Private(FN fn,PetscScalar *As,PetscScalar *Bs,PetscInt m,PetscBool firstonly)
{
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscBLASInt   n,k,ld,lwork,info;
  PetscScalar    *Q,*W,*work,a,x,y,one=1.0,zero=0.0;
  PetscReal      *eig,dummy;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork,rdummy;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ld = n;
  k = firstonly? 1: n;

  /* workspace query and memory allocation */
  lwork = -1;
#if defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,As,&ld,&dummy,&a,&lwork,&rdummy,&info));
  ierr = PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork);CHKERRQ(ierr);
  ierr = PetscMalloc5(m,&eig,m*m,&Q,m*k,&W,lwork,&work,PetscMax(1,3*m-2),&rwork);CHKERRQ(ierr);
#else
  PetscStackCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,As,&ld,&dummy,&a,&lwork,&info));
  ierr = PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork);CHKERRQ(ierr);
  ierr = PetscMalloc4(m,&eig,m*m,&Q,m*k,&W,lwork,&work);CHKERRQ(ierr);
#endif

  /* compute eigendecomposition */
  for (j=0;j<n;j++) for (i=j;i<n;i++) Q[i+j*ld] = As[i+j*ld];
#if defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,Q,&ld,eig,work,&lwork,rwork,&info));
#else
  PetscStackCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n,Q,&ld,eig,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("syev",info);

  /* W = f(Lambda)*Q' */
  for (i=0;i<n;i++) {
    x = fn->alpha*eig[i];
    ierr = (*fn->ops->evaluatefunction)(fn,x,&y);CHKERRQ(ierr);  /* y = f(x) */
    for (j=0;j<k;j++) W[i+j*ld] = PetscConj(Q[j+i*ld])*fn->beta*y;
  }
  /* Bs = Q*W */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&k,&n,&one,Q,&ld,W,&ld,&zero,Bs,&ld));
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscFree5(eig,Q,W,work,rwork);CHKERRQ(ierr);
#else
  ierr = PetscFree4(eig,Q,W,work);CHKERRQ(ierr);
#endif
  ierr = PetscLogFlops(9.0*n*n*n+2.0*n*n*n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   FNEvaluateFunctionMat_Sym_Default - given a symmetric matrix A,
   compute the matrix function as f(A)=Q*f(D)*Q' where the spectral
   decomposition of A is A=Q*D*Q'
*/
static PetscErrorCode FNEvaluateFunctionMat_Sym_Default(FN fn,Mat A,Mat B)
{
  PetscErrorCode ierr;
  PetscInt       m;
  PetscScalar    *As,*Bs;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  ierr = MatDenseGetArray(B,&Bs);CHKERRQ(ierr);
  ierr = MatGetSize(A,&m,NULL);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMat_Sym_Private(fn,As,Bs,m,PETSC_FALSE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(B,&Bs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNEvaluateFunctionMat_Private(FN fn,Mat A,Mat B,PetscBool sync)
{
  PetscErrorCode ierr;
  PetscBool      set,flg,symm=PETSC_FALSE;
  PetscInt       m,n;
  PetscMPIInt    size,rank;
  PetscScalar    *pF;
  Mat            M,F;

  PetscFunctionBegin;
  /* destination matrix */
  F = B?B:A;

  /* check symmetry of A */
  ierr = MatIsHermitianKnown(A,&set,&flg);CHKERRQ(ierr);
  symm = set? flg: PETSC_FALSE;

  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)fn),&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)fn),&rank);CHKERRQ(ierr);
  if (size==1 || fn->pmode==FN_PARALLEL_REDUNDANT || (fn->pmode==FN_PARALLEL_SYNCHRONIZED && !rank)) {

    ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
    if (symm && !fn->method) {  /* prefer diagonalization */
      ierr = PetscInfo(fn,"Computing matrix function via diagonalization\n");CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat_Sym_Default(fn,A,F);CHKERRQ(ierr);
    } else {
      /* scale argument */
      if (fn->alpha!=(PetscScalar)1.0) {
        ierr = FN_AllocateWorkMat(fn,A,&M);CHKERRQ(ierr);
        ierr = MatScale(M,fn->alpha);CHKERRQ(ierr);
      } else M = A;
      if (fn->ops->evaluatefunctionmat[fn->method]) {
        ierr = (*fn->ops->evaluatefunctionmat[fn->method])(fn,M,F);CHKERRQ(ierr);
      } else if (!fn->method) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_SUP,"Matrix functions not implemented in this FN type");
      else SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"The specified method number does not exist for this FN type");
      if (fn->alpha!=(PetscScalar)1.0) {
        ierr = FN_FreeWorkMat(fn,&M);CHKERRQ(ierr);
      }
      /* scale result */
      ierr = MatScale(F,fn->beta);CHKERRQ(ierr);
    }
    ierr = PetscFPTrapPop();CHKERRQ(ierr);
  }
  if (size>1 && fn->pmode==FN_PARALLEL_SYNCHRONIZED && sync) {  /* synchronize */
    ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
    ierr = MatDenseGetArray(F,&pF);CHKERRQ(ierr);
    ierr = MPI_Bcast(pF,n*n,MPIU_SCALAR,0,PetscObjectComm((PetscObject)fn));CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(F,&pF);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   FNEvaluateFunctionMat - Computes the value of the function f(A) for a given
   matrix A, where the result is also a matrix.

   Logically Collective on fn

   Input Parameters:
+  fn - the math function context
-  A  - matrix on which the function must be evaluated

   Output Parameter:
.  B  - (optional) matrix resulting from evaluating f(A)

   Notes:
   Matrix A must be a square sequential dense Mat, with all entries equal on
   all processes (otherwise each process will compute different results).
   If matrix B is provided, it must also be a square sequential dense Mat, and
   both matrices must have the same dimensions. If B is NULL (or B=A) then the
   function will perform an in-place computation, overwriting A with f(A).

   If A is known to be real symmetric or complex Hermitian then it is
   recommended to set the appropriate flag with MatSetOption(), because
   symmetry can sometimes be exploited by the algorithm.

   Scaling factors are taken into account, so the actual function evaluation
   will return beta*f(alpha*A).

   Level: advanced

.seealso: FNEvaluateFunction(), FNEvaluateFunctionMatVec(), FNSetMethod()
@*/
PetscErrorCode FNEvaluateFunctionMat(FN fn,Mat A,Mat B)
{
  PetscErrorCode ierr;
  PetscBool      match,inplace=PETSC_FALSE;
  PetscInt       m,n,n1;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidType(fn,1);
  PetscValidType(A,2);
  if (B) {
    PetscValidHeaderSpecific(B,MAT_CLASSID,3);
    PetscValidType(B,3);
  } else inplace = PETSC_TRUE;
  ierr = PetscObjectTypeCompare((PetscObject)A,MATSEQDENSE,&match);CHKERRQ(ierr);
  if (!match) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_SUP,"Mat A must be of type seqdense");
  ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
  if (m!=n) SETERRQ2(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Mat A is not square (has %D rows, %D cols)",m,n);
  if (!inplace) {
    ierr = PetscObjectTypeCompare((PetscObject)B,MATSEQDENSE,&match);CHKERRQ(ierr);
    if (!match) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_SUP,"Mat B must be of type seqdense");
    n1 = n;
    ierr = MatGetSize(B,&m,&n);CHKERRQ(ierr);
    if (m!=n) SETERRQ2(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Mat B is not square (has %D rows, %D cols)",m,n);
    if (n1!=n) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Matrices A and B must have the same dimension");
  }

  /* evaluate matrix function */
  ierr = PetscLogEventBegin(FN_Evaluate,fn,0,0,0);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMat_Private(fn,A,B,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(FN_Evaluate,fn,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   FNEvaluateFunctionMatVec_Default - computes the full matrix f(A)
   and then copies the first column.
*/
static PetscErrorCode FNEvaluateFunctionMatVec_Default(FN fn,Mat A,Vec v)
{
  PetscErrorCode ierr;
  Mat            F;

  PetscFunctionBegin;
  ierr = FN_AllocateWorkMat(fn,A,&F);CHKERRQ(ierr);
  if (fn->ops->evaluatefunctionmat[fn->method]) {
    ierr = (*fn->ops->evaluatefunctionmat[fn->method])(fn,A,F);CHKERRQ(ierr);
  } else if (!fn->method) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_SUP,"Matrix functions not implemented in this FN type");
  else SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_OUTOFRANGE,"The specified method number does not exist for this FN type");
  ierr = MatGetColumnVector(F,v,0);CHKERRQ(ierr);
  ierr = FN_FreeWorkMat(fn,&F);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   FNEvaluateFunctionMatVec_Sym_Default - given a symmetric matrix A,
   compute the matrix function as f(A)=Q*f(D)*Q' where the spectral
   decomposition of A is A=Q*D*Q'. Only the first column is computed.
*/
static PetscErrorCode FNEvaluateFunctionMatVec_Sym_Default(FN fn,Mat A,Vec v)
{
  PetscErrorCode ierr;
  PetscInt       m;
  PetscScalar    *As,*vs;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  ierr = VecGetArray(v,&vs);CHKERRQ(ierr);
  ierr = MatGetSize(A,&m,NULL);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMat_Sym_Private(fn,As,vs,m,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  ierr = VecRestoreArray(v,&vs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNEvaluateFunctionMatVec_Private(FN fn,Mat A,Vec v,PetscBool sync)
{
  PetscErrorCode ierr;
  PetscBool      set,flg,symm=PETSC_FALSE;
  PetscInt       m,n;
  Mat            M;
  PetscMPIInt    size,rank;
  PetscScalar    *pv;

  PetscFunctionBegin;
  /* check symmetry of A */
  ierr = MatIsHermitianKnown(A,&set,&flg);CHKERRQ(ierr);
  symm = set? flg: PETSC_FALSE;

  /* evaluate matrix function */
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)fn),&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)fn),&rank);CHKERRQ(ierr);
  if (size==1 || fn->pmode==FN_PARALLEL_REDUNDANT || (fn->pmode==FN_PARALLEL_SYNCHRONIZED && !rank)) {
    ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
    if (symm && !fn->method) {  /* prefer diagonalization */
      ierr = PetscInfo(fn,"Computing matrix function via diagonalization\n");CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMatVec_Sym_Default(fn,A,v);CHKERRQ(ierr);
    } else {
      /* scale argument */
      if (fn->alpha!=(PetscScalar)1.0) {
        ierr = FN_AllocateWorkMat(fn,A,&M);CHKERRQ(ierr);
        ierr = MatScale(M,fn->alpha);CHKERRQ(ierr);
      } else M = A;
      if (fn->ops->evaluatefunctionmatvec[fn->method]) {
        ierr = (*fn->ops->evaluatefunctionmatvec[fn->method])(fn,M,v);CHKERRQ(ierr);
      } else {
        ierr = FNEvaluateFunctionMatVec_Default(fn,M,v);CHKERRQ(ierr);
      }
      if (fn->alpha!=(PetscScalar)1.0) {
        ierr = FN_FreeWorkMat(fn,&M);CHKERRQ(ierr);
      }
      /* scale result */
      ierr = VecScale(v,fn->beta);CHKERRQ(ierr);
    }
    ierr = PetscFPTrapPop();CHKERRQ(ierr);
  }

  /* synchronize */
  if (size>1 && fn->pmode==FN_PARALLEL_SYNCHRONIZED && sync) {
    ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
    ierr = VecGetArray(v,&pv);CHKERRQ(ierr);
    ierr = MPI_Bcast(pv,n,MPIU_SCALAR,0,PetscObjectComm((PetscObject)fn));CHKERRQ(ierr);
    ierr = VecRestoreArray(v,&pv);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   FNEvaluateFunctionMatVec - Computes the first column of the matrix f(A)
   for a given matrix A.

   Logically Collective on fn

   Input Parameters:
+  fn - the math function context
-  A  - matrix on which the function must be evaluated

   Output Parameter:
.  v  - vector to hold the first column of f(A)

   Notes:
   This operation is similar to FNEvaluateFunctionMat() but returns only
   the first column of f(A), hence saving computations in most cases.

   Level: advanced

.seealso: FNEvaluateFunction(), FNEvaluateFunctionMat(), FNSetMethod()
@*/
PetscErrorCode FNEvaluateFunctionMatVec(FN fn,Mat A,Vec v)
{
  PetscErrorCode ierr;
  PetscBool      match;
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidHeaderSpecific(v,VEC_CLASSID,3);
  PetscValidType(fn,1);
  PetscValidType(A,2);
  PetscValidType(v,3);
  ierr = PetscObjectTypeCompare((PetscObject)A,MATSEQDENSE,&match);CHKERRQ(ierr);
  if (!match) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_SUP,"Mat A must be of type seqdense");
  ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
  if (m!=n) SETERRQ2(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Mat A is not square (has %D rows, %D cols)",m,n);
  ierr = VecGetSize(v,&m);CHKERRQ(ierr);
  if (m!=n) SETERRQ(PetscObjectComm((PetscObject)fn),PETSC_ERR_ARG_SIZ,"Matrix A and vector v must have the same size");
  ierr = PetscLogEventBegin(FN_Evaluate,fn,0,0,0);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMatVec_Private(fn,A,v,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(FN_Evaluate,fn,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   FNSetFromOptions - Sets FN options from the options database.

   Collective on fn

   Input Parameters:
.  fn - the math function context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner
@*/
PetscErrorCode FNSetFromOptions(FN fn)
{
  PetscErrorCode ierr;
  char           type[256];
  PetscScalar    array[2];
  PetscInt       k,meth;
  PetscBool      flg;
  FNParallelType pmode;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  ierr = FNRegisterAll();CHKERRQ(ierr);
  ierr = PetscObjectOptionsBegin((PetscObject)fn);CHKERRQ(ierr);
    ierr = PetscOptionsFList("-fn_type","Math function type","FNSetType",FNList,(char*)(((PetscObject)fn)->type_name?((PetscObject)fn)->type_name:FNRATIONAL),type,256,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = FNSetType(fn,type);CHKERRQ(ierr);
    } else if (!((PetscObject)fn)->type_name) {
      ierr = FNSetType(fn,FNRATIONAL);CHKERRQ(ierr);
    }

    k = 2;
    array[0] = 0.0; array[1] = 0.0;
    ierr = PetscOptionsScalarArray("-fn_scale","Scale factors (one or two scalar values separated with a comma without spaces)","FNSetScale",array,&k,&flg);CHKERRQ(ierr);
    if (flg) {
      if (k<2) array[1] = 1.0;
      ierr = FNSetScale(fn,array[0],array[1]);CHKERRQ(ierr);
    }

    ierr = PetscOptionsInt("-fn_method","Method to be used for computing matrix functions","FNSetMethod",fn->method,&meth,&flg);CHKERRQ(ierr);
    if (flg) { ierr = FNSetMethod(fn,meth);CHKERRQ(ierr); }

    ierr = PetscOptionsEnum("-fn_parallel","Operation mode in parallel runs","FNSetParallel",FNParallelTypes,(PetscEnum)fn->pmode,(PetscEnum*)&pmode,&flg);CHKERRQ(ierr);
    if (flg) { ierr = FNSetParallel(fn,pmode);CHKERRQ(ierr); }

    if (fn->ops->setfromoptions) {
      ierr = (*fn->ops->setfromoptions)(PetscOptionsObject,fn);CHKERRQ(ierr);
    }
    ierr = PetscObjectProcessOptionsHandlers(PetscOptionsObject,(PetscObject)fn);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   FNView - Prints the FN data structure.

   Collective on fn

   Input Parameters:
+  fn - the math function context
-  viewer - optional visualization context

   Note:
   The available visualization contexts include
+     PETSC_VIEWER_STDOUT_SELF - standard output (default)
-     PETSC_VIEWER_STDOUT_WORLD - synchronized standard
         output where only the first processor opens
         the file.  All other processors send their
         data to the first processor to print.

   The user can open an alternative visualization context with
   PetscViewerASCIIOpen() - output to a specified file.

   Level: beginner
@*/
PetscErrorCode FNView(FN fn,PetscViewer viewer)
{
  PetscBool      isascii;
  PetscErrorCode ierr;
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)fn),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(fn,1,viewer,2);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscObjectPrintClassNamePrefixType((PetscObject)fn,viewer);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)fn),&size);CHKERRQ(ierr);
    if (size>1) {
      ierr = PetscViewerASCIIPrintf(viewer,"  parallel operation mode: %s\n",FNParallelTypes[fn->pmode]);CHKERRQ(ierr);
    }
    if (fn->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*fn->ops->view)(fn,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@C
   FNViewFromOptions - View from options

   Collective on FN

   Input Parameters:
+  fn   - the math function context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: FNView(), FNCreate()
@*/
PetscErrorCode FNViewFromOptions(FN fn,PetscObject obj,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  ierr = PetscObjectViewFromOptions((PetscObject)fn,obj,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   FNDuplicate - Duplicates a math function, copying all parameters, possibly with a
   different communicator.

   Collective on fn

   Input Parameters:
+  fn   - the math function context
-  comm - MPI communicator

   Output Parameter:
.  newfn - location to put the new FN context

   Note:
   In order to use the same MPI communicator as in the original object,
   use PetscObjectComm((PetscObject)fn).

   Level: developer

.seealso: FNCreate()
@*/
PetscErrorCode FNDuplicate(FN fn,MPI_Comm comm,FN *newfn)
{
  PetscErrorCode ierr;
  FNType         type;
  PetscScalar    alpha,beta;
  PetscInt       meth;
  FNParallelType ptype;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(fn,FN_CLASSID,1);
  PetscValidType(fn,1);
  PetscValidPointer(newfn,3);
  ierr = FNCreate(comm,newfn);CHKERRQ(ierr);
  ierr = FNGetType(fn,&type);CHKERRQ(ierr);
  ierr = FNSetType(*newfn,type);CHKERRQ(ierr);
  ierr = FNGetScale(fn,&alpha,&beta);CHKERRQ(ierr);
  ierr = FNSetScale(*newfn,alpha,beta);CHKERRQ(ierr);
  ierr = FNGetMethod(fn,&meth);CHKERRQ(ierr);
  ierr = FNSetMethod(*newfn,meth);CHKERRQ(ierr);
  ierr = FNGetParallel(fn,&ptype);CHKERRQ(ierr);
  ierr = FNSetParallel(*newfn,ptype);CHKERRQ(ierr);
  if (fn->ops->duplicate) {
    ierr = (*fn->ops->duplicate)(fn,comm,newfn);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   FNDestroy - Destroys FN context that was created with FNCreate().

   Collective on fn

   Input Parameter:
.  fn - the math function context

   Level: beginner

.seealso: FNCreate()
@*/
PetscErrorCode FNDestroy(FN *fn)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  if (!*fn) PetscFunctionReturn(0);
  PetscValidHeaderSpecific(*fn,FN_CLASSID,1);
  if (--((PetscObject)(*fn))->refct > 0) { *fn = 0; PetscFunctionReturn(0); }
  if ((*fn)->ops->destroy) { ierr = (*(*fn)->ops->destroy)(*fn);CHKERRQ(ierr); }
  for (i=0;i<(*fn)->nw;i++) {
    ierr = MatDestroy(&(*fn)->W[i]);CHKERRQ(ierr);
  }
  ierr = PetscHeaderDestroy(fn);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   FNRegister - Adds a mathematical function to the FN package.

   Not collective

   Input Parameters:
+  name - name of a new user-defined FN
-  function - routine to create context

   Notes:
   FNRegister() may be called multiple times to add several user-defined functions.

   Level: advanced

.seealso: FNRegisterAll()
@*/
PetscErrorCode FNRegister(const char *name,PetscErrorCode (*function)(FN))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = FNInitializePackage();CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&FNList,name,function);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

