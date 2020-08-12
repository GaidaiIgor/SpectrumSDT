/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the BLZPACK package
*/

#include <slepc/private/epsimpl.h>    /*I "slepceps.h" I*/
#include "blzpack.h"

const char* blzpack_error[33] = {
  "",
  "illegal data, LFLAG ",
  "illegal data, dimension of (U), (V), (X) ",
  "illegal data, leading dimension of (U), (V), (X) ",
  "illegal data, leading dimension of (EIG) ",
  "illegal data, number of required eigenpairs ",
  "illegal data, Lanczos algorithm block size ",
  "illegal data, maximum number of steps ",
  "illegal data, number of starting vectors ",
  "illegal data, number of eigenpairs provided ",
  "illegal data, problem type flag ",
  "illegal data, spectrum slicing flag ",
  "illegal data, eigenvectors purification flag ",
  "illegal data, level of output ",
  "illegal data, output file unit ",
  "illegal data, LCOMM (MPI or PVM) ",
  "illegal data, dimension of ISTOR ",
  "illegal data, convergence threshold ",
  "illegal data, dimension of RSTOR ",
  "illegal data on at least one PE ",
  "ISTOR(3:14) must be equal on all PEs ",
  "RSTOR(1:3) must be equal on all PEs ",
  "not enough space in ISTOR to start eigensolution ",
  "not enough space in RSTOR to start eigensolution ",
  "illegal data, number of negative eigenvalues ",
  "illegal data, entries of V ",
  "illegal data, entries of X ",
  "failure in computational subinterval ",
  "file I/O error, blzpack.__.BQ ",
  "file I/O error, blzpack.__.BX ",
  "file I/O error, blzpack.__.Q ",
  "file I/O error, blzpack.__.X ",
  "parallel interface error "
};

PetscErrorCode EPSSetUp_BLZPACK(EPS eps)
{
  PetscErrorCode ierr;
  PetscInt       listor,lrstor,ncuv,k1,k2,k3,k4;
  EPS_BLZPACK    *blz = (EPS_BLZPACK*)eps->data;
  PetscBool      issinv,istrivial,flg;

  PetscFunctionBegin;
  if (eps->ncv!=PETSC_DEFAULT) {
    if (eps->ncv < PetscMin(eps->nev+10,eps->nev*2)) SETERRQ(PetscObjectComm((PetscObject)eps),0,"Warning: BLZpack recommends that ncv be larger than min(nev+10,nev*2)");
  } else eps->ncv = PetscMin(eps->nev+10,eps->nev*2);
  if (eps->mpd!=PETSC_DEFAULT) { ierr = PetscInfo(eps,"Warning: parameter mpd ignored\n");CHKERRQ(ierr); }
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(1000,eps->n);

  if (!blz->block_size) blz->block_size = 3;
  if (!eps->ishermitian) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Requested method is only available for Hermitian problems");
  if (eps->which==EPS_ALL) {
    if (eps->inta==0.0 && eps->intb==0.0) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Must define a computational interval when using EPS_ALL");
    blz->slice = 1;
  }
  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&issinv);CHKERRQ(ierr);
  if (blz->slice || eps->isgeneralized) {
    if (!issinv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Shift-and-invert ST is needed for generalized problems or spectrum slicing");
  }
  if (blz->slice) {
    if (eps->intb >= PETSC_MAX_REAL) { /* right-open interval */
      if (eps->inta <= PETSC_MIN_REAL) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The defined computational interval should have at least one of their sides bounded");
      ierr = STSetDefaultShift(eps->st,eps->inta);CHKERRQ(ierr);
    } else {
      ierr = STSetDefaultShift(eps->st,eps->intb);CHKERRQ(ierr);
    }
  }
  if (!eps->which) {
    if (issinv) eps->which = EPS_TARGET_REAL;
    else eps->which = EPS_SMALLEST_REAL;
  }
  if ((issinv && eps->which!=EPS_TARGET_REAL && eps->which!=EPS_TARGET_MAGNITUDE && eps->which!=EPS_ALL) || (!issinv && eps->which!=EPS_SMALLEST_REAL)) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Wrong value of eps->which");
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");

  k1 = PetscMin(eps->n,180);
  k2 = blz->block_size;
  k4 = PetscMin(eps->ncv,eps->n);
  k3 = 484+k1*(13+k1*2+k2+PetscMax(18,k2+2))+k2*k2*3+k4*2;

  listor = 123+k1*12;
  ierr = PetscFree(blz->istor);CHKERRQ(ierr);
  ierr = PetscMalloc1((17+listor),&blz->istor);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,(17+listor)*sizeof(PetscBLASInt));CHKERRQ(ierr);
  ierr = PetscBLASIntCast(listor,&blz->istor[14]);CHKERRQ(ierr);

  if (blz->slice) lrstor = eps->nloc*(k2*4+k1*2+k4)+k3;
  else lrstor = eps->nloc*(k2*4+k1)+k3;
lrstor*=10;
  ierr = PetscFree(blz->rstor);CHKERRQ(ierr);
  ierr = PetscMalloc1((4+lrstor),&blz->rstor);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,(4+lrstor)*sizeof(PetscReal));CHKERRQ(ierr);
  blz->rstor[3] = lrstor;

  ncuv = PetscMax(3,blz->block_size);
  ierr = PetscFree(blz->u);CHKERRQ(ierr);
  ierr = PetscMalloc1(ncuv*eps->nloc,&blz->u);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,ncuv*eps->nloc*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscFree(blz->v);CHKERRQ(ierr);
  ierr = PetscMalloc1(ncuv*eps->nloc,&blz->v);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,ncuv*eps->nloc*sizeof(PetscScalar));CHKERRQ(ierr);

  ierr = PetscFree(blz->eig);CHKERRQ(ierr);
  ierr = PetscMalloc1(2*eps->ncv,&blz->eig);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,2*eps->ncv*sizeof(PetscReal));CHKERRQ(ierr);

  if (eps->extraction) { ierr = PetscInfo(eps,"Warning: extraction type ignored\n");CHKERRQ(ierr); }

  ierr = EPSAllocateSolution(eps,0);CHKERRQ(ierr);
  ierr = EPS_SetInnerProduct(eps);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)eps->V,BVVECS,&flg);CHKERRQ(ierr);
  if (flg) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver requires a BV with contiguous storage");
  ierr = RGIsTrivial(eps->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support region filtering");
  if (eps->stopping!=EPSStoppingBasic) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"External packages do not support user-defined stopping test");
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_BLZPACK(EPS eps)
{
  PetscErrorCode ierr;
  EPS_BLZPACK    *blz = (EPS_BLZPACK*)eps->data;
  PetscInt       nn;
  PetscBLASInt   i,nneig,lflag,nvopu;
  Vec            x,y;
  PetscScalar    sigma,*pV;
  Mat            A,M;
  KSP            ksp;
  PC             pc;

  PetscFunctionBegin;
  ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
  ierr = MatCreateVecsEmpty(A,&x,&y);CHKERRQ(ierr);
  ierr = EPSGetStartVector(eps,0,NULL);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,0,0);CHKERRQ(ierr);  /* just for deflation space */
  ierr = BVGetArray(eps->V,&pV);CHKERRQ(ierr);

  if (eps->isgeneralized && !blz->slice) {
    ierr = STGetShift(eps->st,&sigma);CHKERRQ(ierr); /* shift of origin */
    blz->rstor[0]  = sigma;        /* lower limit of eigenvalue interval */
    blz->rstor[1]  = sigma;        /* upper limit of eigenvalue interval */
  } else {
    sigma = 0.0;
    blz->rstor[0]  = eps->inta;    /* lower limit of eigenvalue interval */
    blz->rstor[1]  = eps->intb;    /* upper limit of eigenvalue interval */
  }
  nneig = 0;                       /* no. of eigs less than sigma */

  ierr = PetscBLASIntCast(eps->nloc,&blz->istor[0]);CHKERRQ(ierr); /* no. of rows of U, V, X */
  ierr = PetscBLASIntCast(eps->nloc,&blz->istor[1]);CHKERRQ(ierr); /* leading dim of U, V, X */
  ierr = PetscBLASIntCast(eps->nev,&blz->istor[2]);CHKERRQ(ierr);  /* required eigenpairs */
  ierr = PetscBLASIntCast(eps->ncv,&blz->istor[3]);CHKERRQ(ierr);  /* working eigenpairs */
  blz->istor[4]  = blz->block_size;    /* number of vectors in a block */
  blz->istor[5]  = blz->nsteps;        /* maximun number of steps per run */
  blz->istor[6]  = 1;                  /* number of starting vectors as input */
  blz->istor[7]  = 0;                  /* number of eigenpairs given as input */
  blz->istor[8]  = (blz->slice || eps->isgeneralized) ? 1 : 0;   /* problem type */
  blz->istor[9]  = blz->slice;         /* spectrum slicing */
  blz->istor[10] = eps->isgeneralized ? 1 : 0;   /* solutions refinement (purify) */
  blz->istor[11] = 0;                  /* level of printing */
  blz->istor[12] = 6;                  /* file unit for output */
  ierr = PetscBLASIntCast(MPI_Comm_c2f(PetscObjectComm((PetscObject)eps)),&blz->istor[13]);CHKERRQ(ierr);

  blz->rstor[2]  = eps->tol;           /* threshold for convergence */

  lflag = 0;           /* reverse communication interface flag */

  do {
    BLZpack_(blz->istor,blz->rstor,&sigma,&nneig,blz->u,blz->v,&lflag,&nvopu,blz->eig,pV);

    switch (lflag) {
    case 1:
      /* compute v = OP u */
      for (i=0;i<nvopu;i++) {
        ierr = VecPlaceArray(x,blz->u+i*eps->nloc);CHKERRQ(ierr);
        ierr = VecPlaceArray(y,blz->v+i*eps->nloc);CHKERRQ(ierr);
        if (blz->slice || eps->isgeneralized) {
          ierr = STMatSolve(eps->st,x,y);CHKERRQ(ierr);
        } else {
          ierr = STApply(eps->st,x,y);CHKERRQ(ierr);
        }
        ierr = BVOrthogonalizeVec(eps->V,y,NULL,NULL,NULL);CHKERRQ(ierr);
        ierr = VecResetArray(x);CHKERRQ(ierr);
        ierr = VecResetArray(y);CHKERRQ(ierr);
      }
      /* monitor */
      eps->nconv  = BLZistorr_(blz->istor,"NTEIG",5);
      ierr = EPSMonitor(eps,eps->its,eps->nconv,blz->rstor+BLZistorr_(blz->istor,"IRITZ",5),eps->eigi,blz->rstor+BLZistorr_(blz->istor,"IRITZ",5)+BLZistorr_(blz->istor,"JT",2),BLZistorr_(blz->istor,"NRITZ",5));CHKERRQ(ierr);
      eps->its = eps->its + 1;
      if (eps->its >= eps->max_it || eps->nconv >= eps->nev) lflag = 5;
      break;
    case 2:
      /* compute v = B u */
      for (i=0;i<nvopu;i++) {
        ierr = VecPlaceArray(x,blz->u+i*eps->nloc);CHKERRQ(ierr);
        ierr = VecPlaceArray(y,blz->v+i*eps->nloc);CHKERRQ(ierr);
        ierr = BVApplyMatrix(eps->V,x,y);CHKERRQ(ierr);
        ierr = VecResetArray(x);CHKERRQ(ierr);
        ierr = VecResetArray(y);CHKERRQ(ierr);
      }
      break;
    case 3:
      /* update shift */
      ierr = PetscInfo1(eps,"Factorization update (sigma=%g)\n",(double)sigma);CHKERRQ(ierr);
      ierr = STSetShift(eps->st,sigma);CHKERRQ(ierr);
      ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr);
      ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
      ierr = PCFactorGetMatrix(pc,&M);CHKERRQ(ierr);
      ierr = MatGetInertia(M,&nn,NULL,NULL);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(nn,&nneig);CHKERRQ(ierr);
      break;
    case 4:
      /* copy the initial vector */
      ierr = VecPlaceArray(x,blz->v);CHKERRQ(ierr);
      ierr = BVCopyVec(eps->V,0,x);CHKERRQ(ierr);
      ierr = VecResetArray(x);CHKERRQ(ierr);
      break;
    }

  } while (lflag > 0);

  ierr = BVRestoreArray(eps->V,&pV);CHKERRQ(ierr);

  eps->nconv  = BLZistorr_(blz->istor,"NTEIG",5);
  eps->reason = EPS_CONVERGED_TOL;
  if (blz->slice) eps->nev = eps->nconv;
  for (i=0;i<eps->nconv;i++) eps->eigr[i]=blz->eig[i];

  if (lflag!=0) {
    char msg[2048] = "";
    for (i = 0; i < 33; i++) {
      if (blz->istor[15] & (1 << i)) { ierr = PetscStrcat(msg,blzpack_error[i]);CHKERRQ(ierr); }
    }
    SETERRQ2(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Error in BLZPACK (code=%d): '%s'",blz->istor[15],msg);
  }
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSBackTransform_BLZPACK(EPS eps)
{
  PetscErrorCode ierr;
  EPS_BLZPACK    *blz = (EPS_BLZPACK*)eps->data;

  PetscFunctionBegin;
  if (!blz->slice && !eps->isgeneralized) {
    ierr = EPSBackTransform_Default(eps);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_BLZPACK(EPS eps)
{
  PetscErrorCode ierr;
  EPS_BLZPACK    *blz = (EPS_BLZPACK*)eps->data;

  PetscFunctionBegin;
  ierr = PetscFree(blz->istor);CHKERRQ(ierr);
  ierr = PetscFree(blz->rstor);CHKERRQ(ierr);
  ierr = PetscFree(blz->u);CHKERRQ(ierr);
  ierr = PetscFree(blz->v);CHKERRQ(ierr);
  ierr = PetscFree(blz->eig);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_BLZPACK(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBlzpackSetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBlzpackGetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBlzpackSetNSteps_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBlzpackGetNSteps_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_BLZPACK(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  EPS_BLZPACK    *blz = (EPS_BLZPACK*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  block size=%d\n",blz->block_size);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum number of steps per run=%d\n",blz->nsteps);CHKERRQ(ierr);
    if (blz->slice) {
      ierr = PetscViewerASCIIPrintf(viewer,"  computational interval [%f,%f]\n",eps->inta,eps->intb);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_BLZPACK(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode ierr;
  EPS_BLZPACK    *blz = (EPS_BLZPACK*)eps->data;
  PetscInt       bs,n;
  PetscBool      flg;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS BLZPACK Options");CHKERRQ(ierr);

    bs = blz->block_size;
    ierr = PetscOptionsInt("-eps_blzpack_blocksize","Block size","EPSBlzpackSetBlockSize",bs,&bs,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSBlzpackSetBlockSize(eps,bs);CHKERRQ(ierr); }

    n = blz->nsteps;
    ierr = PetscOptionsInt("-eps_blzpack_nsteps","Number of steps","EPSBlzpackSetNSteps",n,&n,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSBlzpackSetNSteps(eps,n);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSBlzpackSetBlockSize_BLZPACK(EPS eps,PetscInt bs)
{
  PetscErrorCode ierr;
  EPS_BLZPACK    *blz = (EPS_BLZPACK*)eps->data;

  PetscFunctionBegin;
  if (bs == PETSC_DEFAULT) blz->block_size = 3;
  else if (bs <= 0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Block size must be positive");
  else {
    ierr = PetscBLASIntCast(bs,&blz->block_size);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   EPSBlzpackSetBlockSize - Sets the block size for the BLZPACK package.

   Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  bs - block size

   Options Database Key:
.  -eps_blzpack_blocksize - Sets the value of the block size

   Level: advanced
@*/
PetscErrorCode EPSBlzpackSetBlockSize(EPS eps,PetscInt bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,bs,2);
  ierr = PetscTryMethod(eps,"EPSBlzpackSetBlockSize_C",(EPS,PetscInt),(eps,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSBlzpackGetBlockSize_BLZPACK(EPS eps,PetscInt *bs)
{
  EPS_BLZPACK *blz = (EPS_BLZPACK*)eps->data;

  PetscFunctionBegin;
  *bs = blz->block_size;
  PetscFunctionReturn(0);
}

/*@
   EPSBlzpackGetBlockSize - Gets the block size used in the BLZPACK solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  bs - the block size

   Level: advanced

.seealso: EPSBlzpackSetBlockSize()
@*/
PetscErrorCode EPSBlzpackGetBlockSize(EPS eps,PetscInt *bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(bs,2);
  ierr = PetscUseMethod(eps,"EPSBlzpackGetBlockSize_C",(EPS,PetscInt*),(eps,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSBlzpackSetNSteps_BLZPACK(EPS eps,PetscInt nsteps)
{
  PetscErrorCode ierr;
  EPS_BLZPACK    *blz = (EPS_BLZPACK*)eps->data;

  PetscFunctionBegin;
  if (nsteps == PETSC_DEFAULT) blz->nsteps = 0;
  else {
    ierr = PetscBLASIntCast(nsteps,&blz->nsteps);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   EPSBlzpackSetNSteps - Sets the maximum number of steps per run for the BLZPACK
   package.

   Collective on eps

   Input Parameters:
+  eps     - the eigenproblem solver context
-  nsteps  - maximum number of steps

   Options Database Key:
.  -eps_blzpack_nsteps - Sets the maximum number of steps per run

   Level: advanced

@*/
PetscErrorCode EPSBlzpackSetNSteps(EPS eps,PetscInt nsteps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nsteps,2);
  ierr = PetscTryMethod(eps,"EPSBlzpackSetNSteps_C",(EPS,PetscInt),(eps,nsteps));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSBlzpackGetNSteps_BLZPACK(EPS eps,PetscInt *nsteps)
{
  EPS_BLZPACK *blz = (EPS_BLZPACK*)eps->data;

  PetscFunctionBegin;
  *nsteps = blz->nsteps;
  PetscFunctionReturn(0);
}

/*@
   EPSBlzpackGetNSteps - Gets the maximum number of steps per run in the BLZPACK solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  nsteps - the maximum number of steps

   Level: advanced

.seealso: EPSBlzpackSetNSteps()
@*/
PetscErrorCode EPSBlzpackGetNSteps(EPS eps,PetscInt *nsteps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(nsteps,2);
  ierr = PetscUseMethod(eps,"EPSBlzpackGetNSteps_C",(EPS,PetscInt*),(eps,nsteps));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_BLZPACK(EPS eps)
{
  PetscErrorCode ierr;
  EPS_BLZPACK    *blzpack;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&blzpack);CHKERRQ(ierr);
  eps->data = (void*)blzpack;

  eps->ops->solve          = EPSSolve_BLZPACK;
  eps->ops->setup          = EPSSetUp_BLZPACK;
  eps->ops->setfromoptions = EPSSetFromOptions_BLZPACK;
  eps->ops->destroy        = EPSDestroy_BLZPACK;
  eps->ops->reset          = EPSReset_BLZPACK;
  eps->ops->view           = EPSView_BLZPACK;
  eps->ops->backtransform  = EPSBackTransform_BLZPACK;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBlzpackSetBlockSize_C",EPSBlzpackSetBlockSize_BLZPACK);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBlzpackGetBlockSize_C",EPSBlzpackGetBlockSize_BLZPACK);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBlzpackSetNSteps_C",EPSBlzpackSetNSteps_BLZPACK);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSBlzpackGetNSteps_C",EPSBlzpackGetNSteps_BLZPACK);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

