/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the BLOPEX package
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include "blopex.h"
#include <lobpcg.h>
#include <interpreter.h>
#include <multivector.h>
#include <temp_multivector.h>

PetscInt slepc_blopex_useconstr = -1;

typedef struct {
  lobpcg_Tolerance           tol;
  lobpcg_BLASLAPACKFunctions blap_fn;
  mv_InterfaceInterpreter    ii;
  ST                         st;
  Vec                        w;
  PetscInt                   bs;     /* block size */
} EPS_BLOPEX;

static void Precond_FnSingleVector(void *data,void *x,void *y)
{
  PetscErrorCode ierr;
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)data;
  MPI_Comm       comm = PetscObjectComm((PetscObject)blopex->st);
  KSP            ksp;

  PetscFunctionBegin;
  ierr = STGetKSP(blopex->st,&ksp);CHKERRABORT(comm,ierr);
  ierr = KSPSolve(ksp,(Vec)x,(Vec)y);CHKERRABORT(comm,ierr);
  PetscFunctionReturnVoid();
}

static void Precond_FnMultiVector(void *data,void *x,void *y)
{
  EPS_BLOPEX *blopex = (EPS_BLOPEX*)data;

  PetscFunctionBegin;
  blopex->ii.Eval(Precond_FnSingleVector,data,x,y);
  PetscFunctionReturnVoid();
}

static void OperatorASingleVector(void *data,void *x,void *y)
{
  PetscErrorCode ierr;
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)data;
  MPI_Comm       comm = PetscObjectComm((PetscObject)blopex->st);
  Mat            A,B;
  PetscScalar    sigma;
  PetscInt       nmat;

  PetscFunctionBegin;
  ierr = STGetNumMatrices(blopex->st,&nmat);CHKERRABORT(comm,ierr);
  ierr = STGetMatrix(blopex->st,0,&A);CHKERRABORT(comm,ierr);
  if (nmat>1) { ierr = STGetMatrix(blopex->st,1,&B);CHKERRABORT(comm,ierr); }
  ierr = MatMult(A,(Vec)x,(Vec)y);CHKERRABORT(comm,ierr);
  ierr = STGetShift(blopex->st,&sigma);CHKERRABORT(comm,ierr);
  if (sigma != 0.0) {
    if (nmat>1) {
      ierr = MatMult(B,(Vec)x,blopex->w);CHKERRABORT(comm,ierr);
    } else {
      ierr = VecCopy((Vec)x,blopex->w);CHKERRABORT(comm,ierr);
    }
    ierr = VecAXPY((Vec)y,-sigma,blopex->w);CHKERRABORT(comm,ierr);
  }
  PetscFunctionReturnVoid();
}

static void OperatorAMultiVector(void *data,void *x,void *y)
{
  EPS_BLOPEX *blopex = (EPS_BLOPEX*)data;

  PetscFunctionBegin;
  blopex->ii.Eval(OperatorASingleVector,data,x,y);
  PetscFunctionReturnVoid();
}

static void OperatorBSingleVector(void *data,void *x,void *y)
{
  PetscErrorCode ierr;
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)data;
  MPI_Comm       comm = PetscObjectComm((PetscObject)blopex->st);
  Mat            B;

  PetscFunctionBegin;
  ierr = STGetMatrix(blopex->st,1,&B);CHKERRABORT(comm,ierr);
  ierr = MatMult(B,(Vec)x,(Vec)y);CHKERRABORT(comm,ierr);
  PetscFunctionReturnVoid();
}

static void OperatorBMultiVector(void *data,void *x,void *y)
{
  EPS_BLOPEX *blopex = (EPS_BLOPEX*)data;

  PetscFunctionBegin;
  blopex->ii.Eval(OperatorBSingleVector,data,x,y);
  PetscFunctionReturnVoid();
}

PetscErrorCode EPSSetDimensions_BLOPEX(EPS eps,PetscInt nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscErrorCode ierr;
  EPS_BLOPEX     *ctx = (EPS_BLOPEX*)eps->data;
  PetscInt       k;

  PetscFunctionBegin;
  k = ((eps->nev-1)/ctx->bs+1)*ctx->bs;
  if (*ncv!=PETSC_DEFAULT) { /* ncv set */
    if (*ncv<k) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv is not sufficiently large");
  } else *ncv = k;
  if (*mpd==PETSC_DEFAULT) *mpd = *ncv;
  else { ierr = PetscInfo(eps,"Warning: given value of mpd ignored\n");CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetUp_BLOPEX(EPS eps)
{
  PetscErrorCode ierr;
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)eps->data;
  PetscBool      istrivial,flg;
  KSP            ksp;

  PetscFunctionBegin;
  if (!eps->ishermitian || (eps->isgeneralized && !eps->ispositive)) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"blopex only works for Hermitian problems");
  if (!blopex->bs) blopex->bs = PetscMin(16,eps->nev);
  ierr = EPSSetDimensions_BLOPEX(eps,eps->nev,&eps->ncv,&eps->mpd);CHKERRQ(ierr);
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) eps->which = EPS_SMALLEST_REAL;
  if (eps->which!=EPS_SMALLEST_REAL) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Wrong value of eps->which");
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");
  if (eps->stopping!=EPSStoppingBasic) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"External packages do not support user-defined stopping test");
  if (eps->extraction) { ierr = PetscInfo(eps,"Warning: extraction type ignored\n");CHKERRQ(ierr); }
  ierr = RGIsTrivial(eps->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support region filtering");

  blopex->st = eps->st;

  if (eps->converged == EPSConvergedRelative) {
    blopex->tol.absolute = 0.0;
    blopex->tol.relative = eps->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL:eps->tol;
  } else if (eps->converged == EPSConvergedAbsolute) {
    blopex->tol.absolute = eps->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL:eps->tol;
    blopex->tol.relative = 0.0;
  } else SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Convergence test not supported in this solver");

  SLEPCSetupInterpreter(&blopex->ii);

  ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&flg);CHKERRQ(ierr);
  if (!flg) { ierr = PetscInfo(eps,"Warning: ignoring KSP, should use KSPPREONLY\n");CHKERRQ(ierr); }

  /* allocate memory */
  if (!eps->V) { ierr = EPSGetBV(eps,&eps->V);CHKERRQ(ierr); }
  ierr = PetscObjectTypeCompareAny((PetscObject)eps->V,&flg,BVVECS,BVCONTIGUOUS,"");CHKERRQ(ierr);
  if (!flg) {  /* blopex only works with BVVECS or BVCONTIGUOUS */
    ierr = BVSetType(eps->V,BVCONTIGUOUS);CHKERRQ(ierr);
  }
  ierr = EPSAllocateSolution(eps,0);CHKERRQ(ierr);
  if (!blopex->w) {
    ierr = BVCreateVec(eps->V,&blopex->w);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)blopex->w);CHKERRQ(ierr);
  }

#if defined(PETSC_USE_COMPLEX)
  blopex->blap_fn.zpotrf = PETSC_zpotrf_interface;
  blopex->blap_fn.zhegv = PETSC_zsygv_interface;
#else
  blopex->blap_fn.dpotrf = PETSC_dpotrf_interface;
  blopex->blap_fn.dsygv = PETSC_dsygv_interface;
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_BLOPEX(EPS eps)
{
  EPS_BLOPEX        *blopex = (EPS_BLOPEX*)eps->data;
  PetscScalar       sigma,*eigr=NULL;
  PetscReal         *errest=NULL;
  int               i,j,info,its,nconv;
  double            *residhist=NULL;
  PetscErrorCode    ierr;
  mv_MultiVectorPtr eigenvectors,constraints;
#if defined(PETSC_USE_COMPLEX)
  komplex           *lambda=NULL,*lambdahist=NULL;
#else
  double            *lambda=NULL,*lambdahist=NULL;
#endif

  PetscFunctionBegin;
  ierr = STGetShift(eps->st,&sigma);CHKERRQ(ierr);
  ierr = PetscMalloc1(blopex->bs,&lambda);CHKERRQ(ierr);
  if (eps->numbermonitors>0) {
    ierr = PetscMalloc4(blopex->bs*(eps->max_it+1),&lambdahist,eps->ncv,&eigr,blopex->bs*(eps->max_it+1),&residhist,eps->ncv,&errest);CHKERRQ(ierr);
  }

  /* Complete the initial basis with random vectors */
  for (i=0;i<eps->nini;i++) {  /* in case the initial vectors were also set with VecSetRandom */
    ierr = BVSetRandomColumn(eps->V,eps->nini);CHKERRQ(ierr);
  }
  for (i=eps->nini;i<eps->ncv;i++) {
    ierr = BVSetRandomColumn(eps->V,i);CHKERRQ(ierr);
  }

  while (eps->reason == EPS_CONVERGED_ITERATING) {

    /* Create multivector of constraints from leading columns of V */
    ierr = PetscObjectComposedDataSetInt((PetscObject)eps->V,slepc_blopex_useconstr,1);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(eps->V,0,eps->nconv);CHKERRQ(ierr);
    constraints = mv_MultiVectorCreateFromSampleVector(&blopex->ii,eps->nds+eps->nconv,eps->V);

    /* Create multivector where eigenvectors of this run will be stored */
    ierr = PetscObjectComposedDataSetInt((PetscObject)eps->V,slepc_blopex_useconstr,0);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(eps->V,eps->nconv,eps->nconv+blopex->bs);CHKERRQ(ierr);
    eigenvectors = mv_MultiVectorCreateFromSampleVector(&blopex->ii,blopex->bs,eps->V);

#if defined(PETSC_USE_COMPLEX)
    info = lobpcg_solve_complex(eigenvectors,blopex,OperatorAMultiVector,
          eps->isgeneralized?blopex:NULL,eps->isgeneralized?OperatorBMultiVector:NULL,
          blopex,Precond_FnMultiVector,constraints,
          blopex->blap_fn,blopex->tol,eps->max_it,0,&its,
          lambda,lambdahist,blopex->bs,eps->errest+eps->nconv,residhist,blopex->bs);
#else
    info = lobpcg_solve_double(eigenvectors,blopex,OperatorAMultiVector,
          eps->isgeneralized?blopex:NULL,eps->isgeneralized?OperatorBMultiVector:NULL,
          blopex,Precond_FnMultiVector,constraints,
          blopex->blap_fn,blopex->tol,eps->max_it,0,&its,
          lambda,lambdahist,blopex->bs,eps->errest+eps->nconv,residhist,blopex->bs);
#endif
    if (info>0) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"BLOPEX failed with exit code=%d",info);
    mv_MultiVectorDestroy(constraints);
    mv_MultiVectorDestroy(eigenvectors);

    for (j=0;j<blopex->bs;j++) {
#if defined(PETSC_USE_COMPLEX)
      eps->eigr[eps->nconv+j] = PetscCMPLX(lambda[j].real,lambda[j].imag);
#else
      eps->eigr[eps->nconv+j] = lambda[j];
#endif
    }

    if (eps->numbermonitors>0) {
      for (i=0;i<its;i++) {
        nconv = 0;
        for (j=0;j<blopex->bs;j++) {
#if defined(PETSC_USE_COMPLEX)
          eigr[eps->nconv+j] = PetscCMPLX(lambdahist[j+i*blopex->bs].real,lambdahist[j+i*blopex->bs].imag);
#else
          eigr[eps->nconv+j] = lambdahist[j+i*blopex->bs];
#endif
          errest[eps->nconv+j] = residhist[j+i*blopex->bs];
          if (residhist[j+i*blopex->bs]<=eps->tol) nconv++;
        }
        ierr = EPSMonitor(eps,eps->its+i,eps->nconv+nconv,eigr,eps->eigi,errest,eps->nconv+blopex->bs);CHKERRQ(ierr);
      }
    }

    eps->its += its;
    if (info==-1) {
      eps->reason = EPS_DIVERGED_ITS;
      break;
    } else {
      for (i=0;i<blopex->bs;i++) {
        if (sigma != 0.0) eps->eigr[eps->nconv+i] += sigma;
      }
      eps->nconv += blopex->bs;
      if (eps->nconv>=eps->nev) eps->reason = EPS_CONVERGED_TOL;
    }
  }

  ierr = PetscFree(lambda);CHKERRQ(ierr);
  if (eps->numbermonitors>0) {
    ierr = PetscFree4(lambdahist,eigr,residhist,errest);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSBLOPEXSetBlockSize_BLOPEX(EPS eps,PetscInt bs)
{
  EPS_BLOPEX *ctx = (EPS_BLOPEX*)eps->data;

  PetscFunctionBegin;
  if (bs==PETSC_DEFAULT) {
    ctx->bs    = 0;
    eps->state = EPS_STATE_INITIAL;
  } else {
    if (bs<=0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Block size must be >0");
    ctx->bs = bs;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSBLOPEXSetBlockSize - Sets the block size of the BLOPEX solver.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  bs  - the block size

   Options Database Key:
.  -eps_blopex_blocksize - Sets the block size

   Level: advanced

.seealso: EPSBLOPEXGetBlockSize()
@*/
PetscErrorCode EPSBLOPEXSetBlockSize(EPS eps,PetscInt bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,bs,2);
  ierr = PetscTryMethod(eps,"EPSBLOPEXSetBlockSize_C",(EPS,PetscInt),(eps,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSBLOPEXGetBlockSize_BLOPEX(EPS eps,PetscInt *bs)
{
  EPS_BLOPEX *ctx = (EPS_BLOPEX*)eps->data;

  PetscFunctionBegin;
  *bs = ctx->bs;
  PetscFunctionReturn(0);
}

/*@
   EPSBLOPEXGetBlockSize - Gets the block size used in the BLOPEX solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  bs - the block size

   Level: advanced

.seealso: EPSBLOPEXSetBlockSize()
@*/
PetscErrorCode EPSBLOPEXGetBlockSize(EPS eps,PetscInt *bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(bs,2);
  ierr = PetscUseMethod(eps,"EPSBLOPEXGetBlockSize_C",(EPS,PetscInt*),(eps,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_BLOPEX(EPS eps)
{
  PetscErrorCode ierr;
  EPS_BLOPEX     *blopex = (EPS_BLOPEX*)eps->data;

  PetscFunctionBegin;
  ierr = VecDestroy(&blopex->w);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_BLOPEX(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  LOBPCG_DestroyRandomContext();
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBLOPEXSetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBLOPEXGetBlockSize_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_BLOPEX(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  EPS_BLOPEX     *ctx = (EPS_BLOPEX*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  block size %D\n",ctx->bs);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_BLOPEX(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      flg;
  PetscInt       bs;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS BLOPEX Options");CHKERRQ(ierr);

    ierr = PetscOptionsInt("-eps_blopex_blocksize","Block size","EPSBLOPEXSetBlockSize",20,&bs,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSBLOPEXSetBlockSize(eps,bs);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);

  LOBPCG_SetFromOptionsRandomContext();
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_BLOPEX(EPS eps)
{
  EPS_BLOPEX     *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&ctx);CHKERRQ(ierr);
  eps->data = (void*)ctx;

  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_BLOPEX;
  eps->ops->setup          = EPSSetUp_BLOPEX;
  eps->ops->setfromoptions = EPSSetFromOptions_BLOPEX;
  eps->ops->destroy        = EPSDestroy_BLOPEX;
  eps->ops->reset          = EPSReset_BLOPEX;
  eps->ops->view           = EPSView_BLOPEX;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_GMRES;

  LOBPCG_InitRandomContext(PetscObjectComm((PetscObject)eps),NULL);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBLOPEXSetBlockSize_C",EPSBLOPEXSetBlockSize_BLOPEX);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBLOPEXGetBlockSize_C",EPSBLOPEXGetBlockSize_BLOPEX);CHKERRQ(ierr);
  if (slepc_blopex_useconstr < 0) { ierr = PetscObjectComposedDataRegister(&slepc_blopex_useconstr);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

