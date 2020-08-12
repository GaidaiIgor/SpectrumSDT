/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the PRIMME SVD solver
*/

#include <slepc/private/svdimpl.h>    /*I "slepcsvd.h" I*/

#include <primme.h>

#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#define PRIMME_DRIVER cprimme_svds
#else
#define PRIMME_DRIVER zprimme_svds
#endif
#else
#if defined(PETSC_USE_REAL_SINGLE)
#define PRIMME_DRIVER sprimme_svds
#else
#define PRIMME_DRIVER dprimme_svds
#endif
#endif

#if defined(PRIMME_VERSION_MAJOR) && PRIMME_VERSION_MAJOR*100+PRIMME_VERSION_MINOR >= 202
#define SLEPC_HAVE_PRIMME2p2
#endif

typedef struct {
  primme_svds_params        primme;   /* param struct */
  PetscInt                  bs;       /* block size */
  primme_svds_preset_method method;   /* primme method */
  SVD                       svd;      /* reference to the solver */
  Vec                       x,y;      /* auxiliary vectors */
} SVD_PRIMME;

static void multMatvec_PRIMME(void*,PRIMME_INT*,void*,PRIMME_INT*,int*,int*,struct primme_svds_params*,int*);

static void par_GlobalSumReal(void *sendBuf,void *recvBuf,int *count,primme_svds_params *primme,int *ierr)
{
  if (sendBuf == recvBuf) {
    *ierr = MPI_Allreduce(MPI_IN_PLACE,recvBuf,*count,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)primme->commInfo));CHKERRABORT(PetscObjectComm((PetscObject)primme->commInfo),*ierr);
  } else {
    *ierr = MPI_Allreduce(sendBuf,recvBuf,*count,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)primme->commInfo));CHKERRABORT(PetscObjectComm((PetscObject)primme->commInfo),*ierr);
  }
}

#if defined(SLEPC_HAVE_PRIMME3)
static void par_broadcastReal(void *buf,int *count,primme_svds_params *primme,int *ierr)
{
  *ierr = MPI_Bcast(buf,*count,MPIU_REAL,0/*root*/,PetscObjectComm((PetscObject)primme->commInfo));CHKERRABORT(PetscObjectComm((PetscObject)primme->commInfo),*ierr);
}
#endif

#if defined(SLEPC_HAVE_PRIMME2p2)
static void convTestFun(double *sval,void *leftsvec,void *rightsvec,double *resNorm,
#if defined(SLEPC_HAVE_PRIMME3)
                        int *method,
#endif
                        int *isconv,struct primme_svds_params *primme,int *err)
{
  PetscErrorCode ierr;
  SVD            svd = (SVD)primme->commInfo;
  PetscReal      sigma = sval?*sval:0.0;
  PetscReal      r = resNorm?*resNorm:HUGE_VAL,errest;

  *err = 1;
  ierr = (*svd->converged)(svd,sigma,r,&errest,svd->convergedctx);CHKERRABORT(PetscObjectComm((PetscObject)svd),ierr);
  *isconv = (errest<=svd->tol?1:0);
  *err = 0;
}

static void monitorFun(void *basisSvals,int *basisSize,int *basisFlags,int *iblock,int *blockSize,void *basisNorms,int *numConverged,void *lockedSvals,int *numLocked,int *lockedFlags,void *lockedNorms,int *inner_its,void *LSRes,
#if defined(SLEPC_HAVE_PRIMME3)
                       const char *msg,double *time,
#endif
                       primme_event *event,int *stage,struct primme_svds_params *primme,int *err)
{

  PetscErrorCode ierr;
  SVD            svd = (SVD)primme->commInfo;
  PetscInt       i,k,nerrest;

  *err = 1;
  switch (*event) {
    case primme_event_outer_iteration:
      /* Update SVD */
      svd->its = primme->stats.numOuterIterations;
      if (numConverged) svd->nconv = *numConverged;
      k = 0;
      if (lockedSvals && numLocked) for (i=0; i<*numLocked && k<svd->ncv; i++) svd->sigma[k++] = ((PetscReal*)lockedSvals)[i];
      nerrest = k;
      if (iblock && blockSize) {
        for (i=0; i<*blockSize && k+iblock[i]<svd->ncv; i++) svd->errest[k+iblock[i]] = ((PetscReal*)basisNorms)[i];
        nerrest = k+(*blockSize>0?1+iblock[*blockSize-1]:0);
      }
      if (basisSvals && basisSize) for (i=0; i<*basisSize && k<svd->ncv; i++) svd->sigma[k++] = ((PetscReal*)basisSvals)[i];
      /* Show progress */
      ierr = SVDMonitor(svd,svd->its,numConverged?*numConverged:0,svd->sigma,svd->errest,nerrest);CHKERRABORT(PetscObjectComm((PetscObject)svd),ierr);
      break;
#if defined(SLEPC_HAVE_PRIMME3)
    case primme_event_message:
      /* Print PRIMME information messages */
      ierr = PetscInfo(svd,msg);CHKERRABORT(PetscObjectComm((PetscObject)svd),ierr);
      break;
#endif
    default:
      break;
  }
  *err = 0;
}
#endif /* SLEPC_HAVE_PRIMME2p2 */

static void multMatvec_PRIMME(void *xa,PRIMME_INT *ldx,void *ya,PRIMME_INT *ldy,int *blockSize,int *transpose,struct primme_svds_params *primme,int *ierr)
{
  PetscInt   i;
  SVD_PRIMME *ops = (SVD_PRIMME*)primme->matrix;
  Vec        x = ops->x,y = ops->y;
  SVD        svd = ops->svd;

  PetscFunctionBegin;
  for (i=0;i<*blockSize;i++) {
    if (*transpose) {
      *ierr = VecPlaceArray(y,(PetscScalar*)xa+(*ldx)*i);CHKERRABORT(PetscObjectComm((PetscObject)svd),*ierr);
      *ierr = VecPlaceArray(x,(PetscScalar*)ya+(*ldy)*i);CHKERRABORT(PetscObjectComm((PetscObject)svd),*ierr);
      *ierr = SVDMatMult(svd,PETSC_TRUE,y,x);CHKERRABORT(PetscObjectComm((PetscObject)svd),*ierr);
    } else {
      *ierr = VecPlaceArray(x,(PetscScalar*)xa+(*ldx)*i);CHKERRABORT(PetscObjectComm((PetscObject)svd),*ierr);
      *ierr = VecPlaceArray(y,(PetscScalar*)ya+(*ldy)*i);CHKERRABORT(PetscObjectComm((PetscObject)svd),*ierr);
      *ierr = SVDMatMult(svd,PETSC_FALSE,x,y);CHKERRABORT(PetscObjectComm((PetscObject)svd),*ierr);
    }
    *ierr = VecResetArray(x);CHKERRABORT(PetscObjectComm((PetscObject)svd),*ierr);
    *ierr = VecResetArray(y);CHKERRABORT(PetscObjectComm((PetscObject)svd),*ierr);
  }
  PetscFunctionReturnVoid();
}

PetscErrorCode SVDSetUp_PRIMME(SVD svd)
{
  PetscErrorCode     ierr;
  PetscMPIInt        numProcs,procID;
  PetscInt           n,m,nloc,mloc;
  SVD_PRIMME         *ops = (SVD_PRIMME*)svd->data;
  primme_svds_params *primme = &ops->primme;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)svd),&numProcs);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)svd),&procID);CHKERRQ(ierr);

  /* Check some constraints and set some default values */
  ierr = SVDMatGetSize(svd,&m,&n);CHKERRQ(ierr);
  ierr = SVDMatGetLocalSize(svd,&mloc,&nloc);CHKERRQ(ierr);
  ierr = SVDSetDimensions_Default(svd);CHKERRQ(ierr);
  if (svd->max_it==PETSC_DEFAULT) svd->max_it = PETSC_MAX_INT;
  svd->leftbasis = PETSC_TRUE;
  if (svd->stopping!=SVDStoppingBasic) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"External packages do not support user-defined stopping test");
#if !defined(SLEPC_HAVE_PRIMME2p2)
  if (svd->converged != SVDConvergedAbsolute) { ierr = PetscInfo(svd,"Warning: using absolute convergence test\n");CHKERRQ(ierr); }
#endif

  /* Transfer SLEPc options to PRIMME options */
  primme_svds_free(primme);
  primme_svds_initialize(primme);
  primme->m             = m;
  primme->n             = n;
  primme->mLocal        = mloc;
  primme->nLocal        = nloc;
  primme->numSvals      = svd->nsv;
  primme->matrix        = ops;
  primme->commInfo      = svd;
  primme->maxMatvecs    = svd->max_it;
#if !defined(SLEPC_HAVE_PRIMME2p2)
  primme->eps           = svd->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL:svd->tol;
#endif
  primme->numProcs      = numProcs;
  primme->procID        = procID;
  primme->printLevel    = 1;
  primme->matrixMatvec  = multMatvec_PRIMME;
  primme->globalSumReal = par_GlobalSumReal;
#if defined(SLEPC_HAVE_PRIMME3)
  primme->broadcastReal = par_broadcastReal;
#endif
#if defined(SLEPC_HAVE_PRIMME2p2)
  primme->convTestFun   = convTestFun;
  primme->monitorFun    = monitorFun;
#endif
  if (ops->bs > 0) primme->maxBlockSize = ops->bs;

  switch (svd->which) {
    case SVD_LARGEST:
      primme->target = primme_svds_largest;
      break;
    case SVD_SMALLEST:
      primme->target = primme_svds_smallest;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"'which' value not supported by PRIMME");
      break;
  }

  /* If user sets mpd or ncv, maxBasisSize is modified */
  if (svd->mpd!=PETSC_DEFAULT) {
    primme->maxBasisSize = svd->mpd;
    if (svd->ncv!=PETSC_DEFAULT) { ierr = PetscInfo(svd,"Warning: 'ncv' is ignored by PRIMME\n");CHKERRQ(ierr); }
  } else if (svd->ncv!=PETSC_DEFAULT) primme->maxBasisSize = svd->ncv;

  if (primme_svds_set_method(ops->method,(primme_preset_method)EPS_PRIMME_DEFAULT_MIN_TIME,PRIMME_DEFAULT_METHOD,primme) < 0) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"PRIMME method not valid");

  svd->mpd = primme->maxBasisSize;
  svd->ncv = (primme->locking?svd->nsv:0)+primme->maxBasisSize;
  ops->bs  = primme->maxBlockSize;

  /* Set workspace */
  ierr = SVDAllocateSolution(svd,0);CHKERRQ(ierr);

  /* Prepare auxiliary vectors */
  if (!ops->x) {
    ierr = MatCreateVecsEmpty(svd->A,&ops->x,&ops->y);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)ops->x);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)ops->y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSolve_PRIMME(SVD svd)
{
  PetscErrorCode ierr;
  SVD_PRIMME     *ops = (SVD_PRIMME*)svd->data;
  PetscScalar    *svecs, *a;
  PetscInt       i,ierrprimme;
  PetscReal      *svals,*rnorms;

  PetscFunctionBegin;
  /* Reset some parameters left from previous runs */
  ops->primme.aNorm    = 0.0;
  ops->primme.initSize = svd->nini;
  ops->primme.iseed[0] = -1;
  ops->primme.iseed[1] = -1;
  ops->primme.iseed[2] = -1;
  ops->primme.iseed[3] = -1;

  /* Allocating left and right singular vectors contiguously */
  ierr = PetscCalloc1(ops->primme.numSvals*(ops->primme.mLocal+ops->primme.nLocal),&svecs);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)svd,sizeof(PetscReal)*ops->primme.numSvals*(ops->primme.mLocal+ops->primme.nLocal));CHKERRQ(ierr);

  /* Call PRIMME solver */
  ierr = PetscMalloc2(svd->ncv,&svals,svd->ncv,&rnorms);CHKERRQ(ierr);
  ierrprimme = PRIMME_DRIVER(svals,svecs,rnorms,&ops->primme);
  for (i=0;i<svd->ncv;i++) svd->sigma[i] = svals[i];
  for (i=0;i<svd->ncv;i++) svd->errest[i] = rnorms[i];
  ierr = PetscFree2(svals,rnorms);CHKERRQ(ierr);

  /* Copy left and right singular vectors into svd */
  ierr = BVGetArray(svd->U,&a);CHKERRQ(ierr);
  ierr = PetscArraycpy(a,svecs,ops->primme.mLocal*ops->primme.initSize);CHKERRQ(ierr);
  ierr = BVRestoreArray(svd->U,&a);CHKERRQ(ierr);

  ierr = BVGetArray(svd->V,&a);CHKERRQ(ierr);
  ierr = PetscArraycpy(a,svecs+ops->primme.mLocal*ops->primme.initSize,ops->primme.nLocal*ops->primme.initSize);CHKERRQ(ierr);
  ierr = BVRestoreArray(svd->V,&a);CHKERRQ(ierr);

  ierr = PetscFree(svecs);CHKERRQ(ierr);

  svd->nconv  = ops->primme.initSize >= 0 ? ops->primme.initSize : 0;
  svd->reason = svd->nconv >= svd->nsv ? SVD_CONVERGED_TOL: SVD_DIVERGED_ITS;
  svd->its    = ops->primme.stats.numOuterIterations;

  /* Process PRIMME error code */
  if (ierrprimme == 0) {
    /* no error */
  } else if (ierrprimme%100 == -1) SETERRQ1(PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"PRIMME library failed with error code=%d: unexpected error",ierrprimme);
  else if (ierrprimme%100 == -2) SETERRQ1(PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"PRIMME library failed with error code=%d: allocation error",ierrprimme);
  else if (ierrprimme%100 == -3) {
    /* stop by maximum number of iteration or matvecs */
  } else if (ierrprimme%100 >= -39) SETERRQ1(PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"PRIMME library failed with error code=%d: configuration error; check PRIMME's manual",ierrprimme);
  else SETERRQ1(PetscObjectComm((PetscObject)svd),PETSC_ERR_LIB,"PRIMME library failed with error code=%d: runtime error; check PRIMME's manual",ierrprimme);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDReset_PRIMME(SVD svd)
{
  PetscErrorCode ierr;
  SVD_PRIMME     *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  primme_svds_free(&ops->primme);
  ierr = VecDestroy(&ops->x);CHKERRQ(ierr);
  ierr = VecDestroy(&ops->y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDDestroy_PRIMME(SVD svd)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(svd->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMESetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMEGetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMESetMethod_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMEGetMethod_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDView_PRIMME(SVD svd,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isascii;
  SVD_PRIMME     *ctx = (SVD_PRIMME*)svd->data;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  block size=%D\n",ctx->bs);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  solver method: %s\n",SVDPRIMMEMethods[(SVDPRIMMEMethod)ctx->method]);CHKERRQ(ierr);

    /* Display PRIMME params */
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)svd),&rank);CHKERRQ(ierr);
    if (!rank) primme_svds_display_params(ctx->primme);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSetFromOptions_PRIMME(PetscOptionItems *PetscOptionsObject,SVD svd)
{
  PetscErrorCode  ierr;
  SVD_PRIMME      *ctx = (SVD_PRIMME*)svd->data;
  PetscInt        bs;
  SVDPRIMMEMethod meth;
  PetscBool       flg;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"SVD PRIMME Options");CHKERRQ(ierr);

    ierr = PetscOptionsInt("-svd_primme_blocksize","Maximum block size","SVDPRIMMESetBlockSize",ctx->bs,&bs,&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDPRIMMESetBlockSize(svd,bs);CHKERRQ(ierr); }

    ierr = PetscOptionsEnum("-svd_primme_method","Method for solving the singular value problem","SVDPRIMMESetMethod",SVDPRIMMEMethods,(PetscEnum)ctx->method,(PetscEnum*)&meth,&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDPRIMMESetMethod(svd,meth);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDPRIMMESetBlockSize_PRIMME(SVD svd,PetscInt bs)
{
  SVD_PRIMME *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  if (bs == PETSC_DEFAULT) ops->bs = 0;
  else if (bs <= 0) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"PRIMME: block size must be positive");
  else ops->bs = bs;
  PetscFunctionReturn(0);
}

/*@
   SVDPRIMMESetBlockSize - The maximum block size that PRIMME will try to use.

   Logically Collective on svd

   Input Parameters:
+  svd - the singular value solver context
-  bs - block size

   Options Database Key:
.  -svd_primme_blocksize - Sets the max allowed block size value

   Notes:
   If the block size is not set, the value established by primme_svds_initialize
   is used.

   The user should set the block size based on the architecture specifics
   of the target computer, as well as any a priori knowledge of multiplicities.
   The code does NOT require bs > 1 to find multiple eigenvalues. For some
   methods, keeping bs = 1 yields the best overall performance.

   Level: advanced

.seealso: SVDPRIMMEGetBlockSize()
@*/
PetscErrorCode SVDPRIMMESetBlockSize(SVD svd,PetscInt bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveInt(svd,bs,2);
  ierr = PetscTryMethod(svd,"SVDPRIMMESetBlockSize_C",(SVD,PetscInt),(svd,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDPRIMMEGetBlockSize_PRIMME(SVD svd,PetscInt *bs)
{
  SVD_PRIMME *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  *bs = ops->bs;
  PetscFunctionReturn(0);
}

/*@
   SVDPRIMMEGetBlockSize - Get the maximum block size the code will try to use.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  bs - returned block size

   Level: advanced

.seealso: SVDPRIMMESetBlockSize()
@*/
PetscErrorCode SVDPRIMMEGetBlockSize(SVD svd,PetscInt *bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidIntPointer(bs,2);
  ierr = PetscUseMethod(svd,"SVDPRIMMEGetBlockSize_C",(SVD,PetscInt*),(svd,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDPRIMMESetMethod_PRIMME(SVD svd,SVDPRIMMEMethod method)
{
  SVD_PRIMME *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  ops->method = (primme_svds_preset_method)method;
  PetscFunctionReturn(0);
}

/*@
   SVDPRIMMESetMethod - Sets the method for the PRIMME SVD solver.

   Logically Collective on svd

   Input Parameters:
+  svd - the singular value solver context
-  method - method that will be used by PRIMME

   Options Database Key:
.  -svd_primme_method - Sets the method for the PRIMME SVD solver

   Note:
   If not set, the method defaults to SVD_PRIMME_HYBRID.

   Level: advanced

.seealso: SVDPRIMMEGetMethod(), SVDPRIMMEMethod
@*/
PetscErrorCode SVDPRIMMESetMethod(SVD svd,SVDPRIMMEMethod method)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,method,2);
  ierr = PetscTryMethod(svd,"SVDPRIMMESetMethod_C",(SVD,SVDPRIMMEMethod),(svd,method));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDPRIMMEGetMethod_PRIMME(SVD svd,SVDPRIMMEMethod *method)
{
  SVD_PRIMME *ops = (SVD_PRIMME*)svd->data;

  PetscFunctionBegin;
  *method = (SVDPRIMMEMethod)ops->method;
  PetscFunctionReturn(0);
}

/*@
   SVDPRIMMEGetMethod - Gets the method for the PRIMME SVD solver.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  method - method that will be used by PRIMME

   Level: advanced

.seealso: SVDPRIMMESetMethod(), SVDPRIMMEMethod
@*/
PetscErrorCode SVDPRIMMEGetMethod(SVD svd,SVDPRIMMEMethod *method)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(method,2);
  ierr = PetscUseMethod(svd,"SVDPRIMMEGetMethod_C",(SVD,SVDPRIMMEMethod*),(svd,method));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_PRIMME(SVD svd)
{
  PetscErrorCode ierr;
  SVD_PRIMME     *primme;

  PetscFunctionBegin;
  ierr = PetscNewLog(svd,&primme);CHKERRQ(ierr);
  svd->data = (void*)primme;

  primme_svds_initialize(&primme->primme);
  primme->bs = 0;
  primme->method = (primme_svds_preset_method)SVD_PRIMME_HYBRID;
  primme->svd = svd;

  svd->ops->solve          = SVDSolve_PRIMME;
  svd->ops->setup          = SVDSetUp_PRIMME;
  svd->ops->setfromoptions = SVDSetFromOptions_PRIMME;
  svd->ops->destroy        = SVDDestroy_PRIMME;
  svd->ops->reset          = SVDReset_PRIMME;
  svd->ops->view           = SVDView_PRIMME;

  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMESetBlockSize_C",SVDPRIMMESetBlockSize_PRIMME);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMEGetBlockSize_C",SVDPRIMMEGetBlockSize_PRIMME);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMESetMethod_C",SVDPRIMMESetMethod_PRIMME);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDPRIMMEGetMethod_C",SVDPRIMMEGetMethod_PRIMME);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

