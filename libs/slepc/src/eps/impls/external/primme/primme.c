/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the PRIMME package
*/

#include <slepc/private/epsimpl.h>    /*I "slepceps.h" I*/

#include <primme.h>

#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#define PRIMME_DRIVER cprimme
#else
#define PRIMME_DRIVER zprimme
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
#define PRIMME_DRIVER sprimme
#else
#define PRIMME_DRIVER dprimme
#endif
#endif

#if defined(PRIMME_VERSION_MAJOR) && PRIMME_VERSION_MAJOR*100+PRIMME_VERSION_MINOR >= 202
#define SLEPC_HAVE_PRIMME2p2
#endif

typedef struct {
  primme_params        primme;    /* param struc */
  PetscInt             bs;        /* block size */
  primme_preset_method method;    /* primme method */
  Mat                  A,B;       /* problem matrices */
  KSP                  ksp;       /* linear solver and preconditioner */
  Vec                  x,y;       /* auxiliary vectors */
  double               target;    /* a copy of eps's target */
} EPS_PRIMME;

static void par_GlobalSumReal(void *sendBuf,void *recvBuf,int *count,primme_params *primme,int *ierr)
{
  if (sendBuf == recvBuf) {
    *ierr = MPI_Allreduce(MPI_IN_PLACE,recvBuf,*count,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)primme->commInfo));CHKERRABORT(PetscObjectComm((PetscObject)primme->commInfo),*ierr);
  } else {
    *ierr = MPI_Allreduce(sendBuf,recvBuf,*count,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)primme->commInfo));CHKERRABORT(PetscObjectComm((PetscObject)primme->commInfo),*ierr);
  }
}

#if defined(SLEPC_HAVE_PRIMME3)
static void par_broadcastReal(void *buf,int *count,primme_params *primme,int *ierr)
{
  *ierr = MPI_Bcast(buf,*count,MPIU_REAL,0/*root*/,PetscObjectComm((PetscObject)primme->commInfo));CHKERRABORT(PetscObjectComm((PetscObject)primme->commInfo),*ierr);
}
#endif

#if defined(SLEPC_HAVE_PRIMME2p2)
static void convTestFun(double *eval,void *evec,double *resNorm,int *isconv,primme_params *primme,int *err)
{
  PetscErrorCode ierr;
  EPS            eps = (EPS)primme->commInfo;
  PetscScalar    eigvr = eval?*eval:0.0;
  PetscReal      r = resNorm?*resNorm:HUGE_VAL,errest;

  *err = 1;
  ierr = (*eps->converged)(eps,eigvr,0.0,r,&errest,eps->convergedctx);CHKERRABORT(PetscObjectComm((PetscObject)eps),ierr);
  *isconv = (errest<=eps->tol?1:0);
  *err = 0;
}

static void monitorFun(void *basisEvals,int *basisSize,int *basisFlags,int *iblock,int *blockSize,void *basisNorms,int *numConverged,void *lockedEvals,int *numLocked,int *lockedFlags,void *lockedNorms,int *inner_its,void *LSRes,
#if defined(SLEPC_HAVE_PRIMME3)
                       const char *msg,double *time,
#endif
                       primme_event *event,struct primme_params *primme,int *err)
{
  PetscErrorCode ierr;
  EPS            eps = (EPS)primme->commInfo;
  PetscInt       i,k,nerrest;

  *err = 1;
  switch (*event) {
    case primme_event_outer_iteration:
      /* Update EPS */
      eps->its = primme->stats.numOuterIterations;
      eps->nconv = primme->initSize;
      k=0;
      if (lockedEvals && numLocked) for (i=0; i<*numLocked && k<eps->ncv; i++) eps->eigr[k++] = ((PetscReal*)lockedEvals)[i];
      nerrest = k;
      if (iblock && blockSize) {
        for (i=0; i<*blockSize && k+iblock[i]<eps->ncv; i++) eps->errest[k+iblock[i]] = ((PetscReal*)basisNorms)[i];
        nerrest = k+(*blockSize>0?1+iblock[*blockSize-1]:0);
      }
      if (basisEvals && basisSize) for (i=0; i<*basisSize && k<eps->ncv; i++) eps->eigr[k++] = ((PetscReal*)basisEvals)[i];
      /* Show progress */
      ierr = EPSMonitor(eps,eps->its,numConverged?*numConverged:0,eps->eigr,eps->eigi,eps->errest,nerrest);CHKERRABORT(PetscObjectComm((PetscObject)eps),ierr);
      break;
#if defined(SLEPC_HAVE_PRIMME3)
    case primme_event_message:
      /* Print PRIMME information messages */
      ierr = PetscInfo(eps,msg);CHKERRABORT(PetscObjectComm((PetscObject)eps),ierr);
      break;
#endif
    default:
      break;
  }
  *err = 0;
}
#endif /* SLEPC_HAVE_PRIMME2p2 */

static void matrixMatvec_PRIMME(void *xa,PRIMME_INT *ldx,void *ya,PRIMME_INT *ldy,int *blockSize,struct primme_params *primme,int *ierr)
{
  PetscInt   i;
  EPS_PRIMME *ops = (EPS_PRIMME*)primme->matrix;
  Vec        x = ops->x,y = ops->y;
  Mat        A = ops->A;

  PetscFunctionBegin;
  for (i=0;i<*blockSize;i++) {
    *ierr = VecPlaceArray(x,(PetscScalar*)xa+(*ldx)*i);CHKERRABORT(PetscObjectComm((PetscObject)A),*ierr);
    *ierr = VecPlaceArray(y,(PetscScalar*)ya+(*ldy)*i);CHKERRABORT(PetscObjectComm((PetscObject)A),*ierr);
    *ierr = MatMult(A,x,y);CHKERRABORT(PetscObjectComm((PetscObject)A),*ierr);
    *ierr = VecResetArray(x);CHKERRABORT(PetscObjectComm((PetscObject)A),*ierr);
    *ierr = VecResetArray(y);CHKERRABORT(PetscObjectComm((PetscObject)A),*ierr);
  }
  PetscFunctionReturnVoid();
}

#if defined(SLEPC_HAVE_PRIMME3)
static void massMatrixMatvec_PRIMME(void *xa,PRIMME_INT *ldx,void *ya,PRIMME_INT *ldy,int *blockSize,struct primme_params *primme,int *ierr)
{
  PetscInt   i;
  EPS_PRIMME *ops = (EPS_PRIMME*)primme->massMatrix;
  Vec        x = ops->x,y = ops->y;
  Mat        B = ops->B;

  PetscFunctionBegin;
  for (i=0;i<*blockSize;i++) {
    *ierr = VecPlaceArray(x,(PetscScalar*)xa+(*ldx)*i);CHKERRABORT(PetscObjectComm((PetscObject)B),*ierr);
    *ierr = VecPlaceArray(y,(PetscScalar*)ya+(*ldy)*i);CHKERRABORT(PetscObjectComm((PetscObject)B),*ierr);
    *ierr = MatMult(B,x,y);CHKERRABORT(PetscObjectComm((PetscObject)B),*ierr);
    *ierr = VecResetArray(x);CHKERRABORT(PetscObjectComm((PetscObject)B),*ierr);
    *ierr = VecResetArray(y);CHKERRABORT(PetscObjectComm((PetscObject)B),*ierr);
  }
  PetscFunctionReturnVoid();
}
#endif

static void applyPreconditioner_PRIMME(void *xa,PRIMME_INT *ldx,void *ya,PRIMME_INT *ldy,int *blockSize,struct primme_params *primme,int *ierr)
{
  PetscInt   i;
  EPS_PRIMME *ops = (EPS_PRIMME*)primme->matrix;
  Vec        x = ops->x,y = ops->y;

  PetscFunctionBegin;
  for (i=0;i<*blockSize;i++) {
    *ierr = VecPlaceArray(x,(PetscScalar*)xa+(*ldx)*i);CHKERRABORT(PetscObjectComm((PetscObject)ops->ksp),*ierr);
    *ierr = VecPlaceArray(y,(PetscScalar*)ya+(*ldy)*i);CHKERRABORT(PetscObjectComm((PetscObject)ops->ksp),*ierr);
    *ierr = KSPSolve(ops->ksp,x,y);CHKERRABORT(PetscObjectComm((PetscObject)ops->ksp),*ierr);
    *ierr = VecResetArray(x);CHKERRABORT(PetscObjectComm((PetscObject)ops->ksp),*ierr);
    *ierr = VecResetArray(y);CHKERRABORT(PetscObjectComm((PetscObject)ops->ksp),*ierr);
  }
  PetscFunctionReturnVoid();
}

PetscErrorCode EPSSetUp_PRIMME(EPS eps)
{
  PetscErrorCode ierr;
  PetscMPIInt    numProcs,procID;
  EPS_PRIMME     *ops = (EPS_PRIMME*)eps->data;
  primme_params  *primme = &ops->primme;
  PetscBool      istrivial,flg;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)eps),&numProcs);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)eps),&procID);CHKERRQ(ierr);

  /* Check some constraints and set some default values */
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PETSC_MAX_INT;
  ierr = STGetMatrix(eps->st,0,&ops->A);CHKERRQ(ierr);
  if (!eps->ishermitian) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"PRIMME is only available for Hermitian problems");
  if (eps->isgeneralized) {
#if defined(SLEPC_HAVE_PRIMME3)
    ierr = STGetMatrix(eps->st,1,&ops->B);CHKERRQ(ierr);
#else
    SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This version of PRIMME is not available for generalized problems");
#endif
  }
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");
  if (eps->stopping!=EPSStoppingBasic) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"External packages do not support user-defined stopping test");
  if (!eps->which) eps->which = EPS_LARGEST_REAL;
#if !defined(SLEPC_HAVE_PRIMME2p2)
  if (eps->converged != EPSConvergedAbsolute) { ierr = PetscInfo(eps,"Warning: using absolute convergence test\n");CHKERRQ(ierr); }
#endif
  ierr = RGIsTrivial(eps->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support region filtering");

  /* Transfer SLEPc options to PRIMME options */
  primme_free(primme);
  primme_initialize(primme);
  primme->n                             = eps->n;
  primme->nLocal                        = eps->nloc;
  primme->numEvals                      = eps->nev;
  primme->matrix                        = ops;
  primme->matrixMatvec                  = matrixMatvec_PRIMME;
#if defined(SLEPC_HAVE_PRIMME3)
  if (eps->isgeneralized) {
    primme->massMatrix                  = ops;
    primme->massMatrixMatvec            = massMatrixMatvec_PRIMME;
  }
#endif
  primme->commInfo                      = eps;
  primme->maxOuterIterations            = eps->max_it;
#if !defined(SLEPC_HAVE_PRIMME2p2)
  primme->eps                           = eps->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL:eps->tol;
#endif
  primme->numProcs                      = numProcs;
  primme->procID                        = procID;
  primme->printLevel                    = 1;
  primme->correctionParams.precondition = 1;
  primme->globalSumReal                 = par_GlobalSumReal;
#if defined(SLEPC_HAVE_PRIMME3)
  primme->broadcastReal                 = par_broadcastReal;
#endif
#if defined(SLEPC_HAVE_PRIMME2p2)
  primme->convTestFun                   = convTestFun;
  primme->monitorFun                    = monitorFun;
#endif
  if (ops->bs > 0) primme->maxBlockSize = ops->bs;

  switch (eps->which) {
    case EPS_LARGEST_REAL:
      primme->target = primme_largest;
      break;
    case EPS_SMALLEST_REAL:
      primme->target = primme_smallest;
      break;
    case EPS_LARGEST_MAGNITUDE:
      primme->target = primme_largest_abs;
      ops->target = 0.0;
      primme->numTargetShifts = 1;
      primme->targetShifts = &ops->target;
      break;
    case EPS_SMALLEST_MAGNITUDE:
      primme->target = primme_closest_abs;
      ops->target = 0.0;
      primme->numTargetShifts = 1;
      primme->targetShifts = &ops->target;
      break;
    case EPS_TARGET_MAGNITUDE:
    case EPS_TARGET_REAL:
      primme->target = primme_closest_abs;
      primme->numTargetShifts = 1;
      ops->target = PetscRealPart(eps->target);
      primme->targetShifts = &ops->target;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"'which' value not supported by PRIMME");
      break;
  }

  switch (eps->extraction) {
    case EPS_RITZ:
      primme->projectionParams.projection = primme_proj_RR;
      break;
    case EPS_HARMONIC:
      primme->projectionParams.projection = primme_proj_harmonic;
      break;
    case EPS_REFINED:
      primme->projectionParams.projection = primme_proj_refined;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"'extraction' value not supported by PRIMME");
      break;
  }

  /* If user sets mpd or ncv, maxBasisSize is modified */
  if (eps->mpd!=PETSC_DEFAULT) {
    primme->maxBasisSize = eps->mpd;
    if (eps->ncv!=PETSC_DEFAULT) { ierr = PetscInfo(eps,"Warning: 'ncv' is ignored by PRIMME\n");CHKERRQ(ierr); }
  } else if (eps->ncv!=PETSC_DEFAULT) primme->maxBasisSize = eps->ncv;

  if (primme_set_method(ops->method,primme) < 0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"PRIMME method not valid");

  eps->mpd = primme->maxBasisSize;
  eps->ncv = (primme->locking?eps->nev:0)+primme->maxBasisSize;
  ops->bs  = primme->maxBlockSize;

  /* Set workspace */
  ierr = EPSAllocateSolution(eps,0);CHKERRQ(ierr);

  /* Setup the preconditioner */
  if (primme->correctionParams.precondition) {
    ierr = STGetKSP(eps->st,&ops->ksp);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)ops->ksp,KSPPREONLY,&flg);CHKERRQ(ierr);
    if (!flg) { ierr = PetscInfo(eps,"Warning: ignoring KSP, should use KSPPREONLY\n");CHKERRQ(ierr); }
    primme->preconditioner = NULL;
    primme->applyPreconditioner = applyPreconditioner_PRIMME;
  }

  /* Prepare auxiliary vectors */
  if (!ops->x) {
    ierr = MatCreateVecsEmpty(ops->A,&ops->x,&ops->y);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)ops->x);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)ops->y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_PRIMME(EPS eps)
{
  PetscErrorCode ierr;
  EPS_PRIMME     *ops = (EPS_PRIMME*)eps->data;
  PetscScalar    *a;
  PetscInt       i,ierrprimme;
  PetscReal      *evals,*rnorms;

  PetscFunctionBegin;
  /* Reset some parameters left from previous runs */
#if defined(SLEPC_HAVE_PRIMME2p2)
  ops->primme.aNorm    = 0.0;
#else
  /* Force PRIMME to stop by absolute error */
  ops->primme.aNorm    = 1.0;
#endif
  ops->primme.initSize = eps->nini;
  ops->primme.iseed[0] = -1;
  ops->primme.iseed[1] = -1;
  ops->primme.iseed[2] = -1;
  ops->primme.iseed[3] = -1;

  /* Call PRIMME solver */
  ierr = BVGetArray(eps->V,&a);CHKERRQ(ierr);
  ierr = PetscMalloc2(eps->ncv,&evals,eps->ncv,&rnorms);CHKERRQ(ierr);
  ierrprimme = PRIMME_DRIVER(evals,a,rnorms,&ops->primme);
  for (i=0;i<eps->ncv;i++) eps->eigr[i] = evals[i];
  for (i=0;i<eps->ncv;i++) eps->errest[i] = rnorms[i];
  ierr = PetscFree2(evals,rnorms);CHKERRQ(ierr);
  ierr = BVRestoreArray(eps->V,&a);CHKERRQ(ierr);

  eps->nconv  = ops->primme.initSize >= 0 ? ops->primme.initSize : 0;
  eps->reason = eps->nconv >= eps->nev ? EPS_CONVERGED_TOL: EPS_DIVERGED_ITS;
  eps->its    = ops->primme.stats.numOuterIterations;

  /* Process PRIMME error code */
  if (ierrprimme == 0) {
    /* no error */
  } else if (ierrprimme == -1) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"PRIMME library failed with error code=%d: unexpected error",ierrprimme);
  else if (ierrprimme == -2) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"PRIMME library failed with error code=%d: allocation error",ierrprimme);
  else if (ierrprimme == -3) {
    /* stop by maximum number of iteration or matvecs */
  } else if (ierrprimme >= -39) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"PRIMME library failed with error code=%d: configuration error; check PRIMME's manual",ierrprimme);
  else SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"PRIMME library failed with error code=%d: runtime error; check PRIMME's manual",ierrprimme);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_PRIMME(EPS eps)
{
  PetscErrorCode ierr;
  EPS_PRIMME     *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  primme_free(&ops->primme);
  ierr = VecDestroy(&ops->x);CHKERRQ(ierr);
  ierr = VecDestroy(&ops->y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_PRIMME(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMESetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMESetMethod_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMEGetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMEGetMethod_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_PRIMME(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isascii;
  EPS_PRIMME     *ctx = (EPS_PRIMME*)eps->data;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  block size=%D\n",ctx->bs);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  solver method: %s\n",EPSPRIMMEMethods[(EPSPRIMMEMethod)ctx->method]);CHKERRQ(ierr);

    /* Display PRIMME params */
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)eps),&rank);CHKERRQ(ierr);
    if (!rank) primme_display_params(ctx->primme);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_PRIMME(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode  ierr;
  EPS_PRIMME      *ctx = (EPS_PRIMME*)eps->data;
  PetscInt        bs;
  EPSPRIMMEMethod meth;
  PetscBool       flg;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS PRIMME Options");CHKERRQ(ierr);

    ierr = PetscOptionsInt("-eps_primme_blocksize","Maximum block size","EPSPRIMMESetBlockSize",ctx->bs,&bs,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSPRIMMESetBlockSize(eps,bs);CHKERRQ(ierr); }

    ierr = PetscOptionsEnum("-eps_primme_method","Method for solving the eigenproblem","EPSPRIMMESetMethod",EPSPRIMMEMethods,(PetscEnum)ctx->method,(PetscEnum*)&meth,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSPRIMMESetMethod(eps,meth);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPRIMMESetBlockSize_PRIMME(EPS eps,PetscInt bs)
{
  EPS_PRIMME *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  if (bs == PETSC_DEFAULT) ops->bs = 0;
  else if (bs <= 0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"PRIMME: block size must be positive");
  else ops->bs = bs;
  PetscFunctionReturn(0);
}

/*@
   EPSPRIMMESetBlockSize - The maximum block size that PRIMME will try to use.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  bs - block size

   Options Database Key:
.  -eps_primme_blocksize - Sets the max allowed block size value

   Notes:
   If the block size is not set, the value established by primme_initialize
   is used.

   The user should set the block size based on the architecture specifics
   of the target computer, as well as any a priori knowledge of multiplicities.
   The code does NOT require bs > 1 to find multiple eigenvalues. For some
   methods, keeping bs = 1 yields the best overall performance.

   Level: advanced

.seealso: EPSPRIMMEGetBlockSize()
@*/
PetscErrorCode EPSPRIMMESetBlockSize(EPS eps,PetscInt bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,bs,2);
  ierr = PetscTryMethod(eps,"EPSPRIMMESetBlockSize_C",(EPS,PetscInt),(eps,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPRIMMEGetBlockSize_PRIMME(EPS eps,PetscInt *bs)
{
  EPS_PRIMME *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  *bs = ops->bs;
  PetscFunctionReturn(0);
}

/*@
   EPSPRIMMEGetBlockSize - Get the maximum block size the code will try to use.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  bs - returned block size

   Level: advanced

.seealso: EPSPRIMMESetBlockSize()
@*/
PetscErrorCode EPSPRIMMEGetBlockSize(EPS eps,PetscInt *bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(bs,2);
  ierr = PetscUseMethod(eps,"EPSPRIMMEGetBlockSize_C",(EPS,PetscInt*),(eps,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPRIMMESetMethod_PRIMME(EPS eps,EPSPRIMMEMethod method)
{
  EPS_PRIMME *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  ops->method = (primme_preset_method)method;
  PetscFunctionReturn(0);
}

/*@
   EPSPRIMMESetMethod - Sets the method for the PRIMME library.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  method - method that will be used by PRIMME

   Options Database Key:
.  -eps_primme_method - Sets the method for the PRIMME library

   Note:
   If not set, the method defaults to EPS_PRIMME_DEFAULT_MIN_TIME.

   Level: advanced

.seealso: EPSPRIMMEGetMethod(), EPSPRIMMEMethod
@*/
PetscErrorCode EPSPRIMMESetMethod(EPS eps,EPSPRIMMEMethod method)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,method,2);
  ierr = PetscTryMethod(eps,"EPSPRIMMESetMethod_C",(EPS,EPSPRIMMEMethod),(eps,method));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPRIMMEGetMethod_PRIMME(EPS eps,EPSPRIMMEMethod *method)
{
  EPS_PRIMME *ops = (EPS_PRIMME*)eps->data;

  PetscFunctionBegin;
  *method = (EPSPRIMMEMethod)ops->method;
  PetscFunctionReturn(0);
}

/*@
   EPSPRIMMEGetMethod - Gets the method for the PRIMME library.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  method - method that will be used by PRIMME

   Level: advanced

.seealso: EPSPRIMMESetMethod(), EPSPRIMMEMethod
@*/
PetscErrorCode EPSPRIMMEGetMethod(EPS eps,EPSPRIMMEMethod *method)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(method,2);
  ierr = PetscUseMethod(eps,"EPSPRIMMEGetMethod_C",(EPS,EPSPRIMMEMethod*),(eps,method));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_PRIMME(EPS eps)
{
  PetscErrorCode ierr;
  EPS_PRIMME     *primme;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&primme);CHKERRQ(ierr);
  eps->data = (void*)primme;

  primme_initialize(&primme->primme);
  primme->primme.globalSumReal = par_GlobalSumReal;
#if defined(SLEPC_HAVE_PRIMME3)
  primme->primme.broadcastReal = par_broadcastReal;
#endif
#if defined(SLEPC_HAVE_PRIMME2p2)
  primme->primme.convTestFun = convTestFun;
  primme->primme.monitorFun = monitorFun;
#endif
  primme->method = (primme_preset_method)EPS_PRIMME_DEFAULT_MIN_TIME;

  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_PRIMME;
  eps->ops->setup          = EPSSetUp_PRIMME;
  eps->ops->setfromoptions = EPSSetFromOptions_PRIMME;
  eps->ops->destroy        = EPSDestroy_PRIMME;
  eps->ops->reset          = EPSReset_PRIMME;
  eps->ops->view           = EPSView_PRIMME;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_GMRES;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMESetBlockSize_C",EPSPRIMMESetBlockSize_PRIMME);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMESetMethod_C",EPSPRIMMESetMethod_PRIMME);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMEGetBlockSize_C",EPSPRIMMEGetBlockSize_PRIMME);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPRIMMEGetMethod_C",EPSPRIMMEGetMethod_PRIMME);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

