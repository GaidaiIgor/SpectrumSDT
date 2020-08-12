/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "arnoldi"

   Method: Explicitly Restarted Arnoldi

   Algorithm:

       Arnoldi method with explicit restart and deflation.

   References:

       [1] "Arnoldi Methods in SLEPc", SLEPc Technical Report STR-4,
           available at https://slepc.upv.es.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/

typedef struct {
  PetscBool delayed;
} EPS_ARNOLDI;

PetscErrorCode EPSSetUp_Arnoldi(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = EPSSetDimensions_Default(eps,eps->nev,&eps->ncv,&eps->mpd);CHKERRQ(ierr);
  if (eps->ncv>eps->nev+eps->mpd) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv must not be larger than nev+mpd");
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) { ierr = EPSSetWhichEigenpairs_Default(eps);CHKERRQ(ierr); }
  if (eps->which==EPS_ALL || (eps->ishermitian && eps->ispositive && (eps->which==EPS_LARGEST_IMAGINARY || eps->which==EPS_SMALLEST_IMAGINARY))) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Wrong value of eps->which");

  if (!eps->extraction) {
    ierr = EPSSetExtraction(eps,EPS_RITZ);CHKERRQ(ierr);
  }
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");

  ierr = EPSAllocateSolution(eps,1);CHKERRQ(ierr);
  ierr = EPS_SetInnerProduct(eps);CHKERRQ(ierr);
  ierr = DSSetType(eps->ds,DSNHEP);CHKERRQ(ierr);
  if (eps->extraction==EPS_REFINED || eps->extraction==EPS_REFINED_HARMONIC) {
    ierr = DSSetRefined(eps->ds,PETSC_TRUE);CHKERRQ(ierr);
  }
  ierr = DSSetExtraRow(eps->ds,PETSC_TRUE);CHKERRQ(ierr);
  ierr = DSAllocate(eps->ds,eps->ncv+1);CHKERRQ(ierr);

  if (eps->isgeneralized && eps->ishermitian && !eps->ispositive) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Requested method does not work for indefinite problems");
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_Arnoldi(EPS eps)
{
  PetscErrorCode     ierr;
  PetscInt           k,nv,ld;
  Mat                U,Op;
  PetscScalar        *H;
  PetscReal          beta,gamma=1.0;
  PetscBool          breakdown,harmonic,refined;
  BVOrthogRefineType orthog_ref;
  EPS_ARNOLDI        *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  ierr = DSGetRefined(eps->ds,&refined);CHKERRQ(ierr);
  harmonic = (eps->extraction==EPS_HARMONIC || eps->extraction==EPS_REFINED_HARMONIC)?PETSC_TRUE:PETSC_FALSE;
  ierr = BVGetOrthogonalization(eps->V,NULL,&orthog_ref,NULL,NULL);CHKERRQ(ierr);

  /* Get the starting Arnoldi vector */
  ierr = EPSGetStartVector(eps,0,NULL);CHKERRQ(ierr);

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Arnoldi factorization */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    ierr = DSSetDimensions(eps->ds,nv,0,eps->nconv,0);CHKERRQ(ierr);
    ierr = DSGetArray(eps->ds,DS_MAT_A,&H);CHKERRQ(ierr);
    if (!arnoldi->delayed) {
      ierr = STGetOperator(eps->st,&Op);CHKERRQ(ierr);
      ierr = BVMatArnoldi(eps->V,Op,H,ld,eps->nconv,&nv,&beta,&breakdown);CHKERRQ(ierr);
      ierr = STRestoreOperator(eps->st,&Op);CHKERRQ(ierr);
    } else if (orthog_ref == BV_ORTHOG_REFINE_NEVER) {
      ierr = EPSDelayedArnoldi1(eps,H,ld,eps->nconv,&nv,&beta,&breakdown);CHKERRQ(ierr);
    } else {
      ierr = EPSDelayedArnoldi(eps,H,ld,eps->nconv,&nv,&beta,&breakdown);CHKERRQ(ierr);
    }
    ierr = DSRestoreArray(eps->ds,DS_MAT_A,&H);CHKERRQ(ierr);
    ierr = DSSetState(eps->ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(eps->V,eps->nconv,nv);CHKERRQ(ierr);

    /* Compute translation of Krylov decomposition if harmonic extraction used */
    if (harmonic) {
      ierr = DSTranslateHarmonic(eps->ds,eps->target,beta,PETSC_FALSE,NULL,&gamma);CHKERRQ(ierr);
    }

    /* Solve projected problem */
    ierr = DSSolve(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);
    ierr = DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DSUpdateExtraRow(eps->ds);CHKERRQ(ierr);
    ierr = DSSynchronize(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);

    /* Check convergence */
    ierr = EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta,0.0,gamma,&k);CHKERRQ(ierr);
    if (refined) {
      ierr = DSGetMat(eps->ds,DS_MAT_X,&U);CHKERRQ(ierr);
      ierr = BVMultInPlace(eps->V,U,eps->nconv,k+1);CHKERRQ(ierr);
      ierr = MatDestroy(&U);CHKERRQ(ierr);
      ierr = BVOrthonormalizeColumn(eps->V,k,PETSC_FALSE,NULL,NULL);CHKERRQ(ierr);
    } else {
      ierr = DSGetMat(eps->ds,DS_MAT_Q,&U);CHKERRQ(ierr);
      ierr = BVMultInPlace(eps->V,U,eps->nconv,PetscMin(k+1,nv));CHKERRQ(ierr);
      ierr = MatDestroy(&U);CHKERRQ(ierr);
    }
    ierr = (*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);
    if (eps->reason == EPS_CONVERGED_ITERATING && breakdown) {
      ierr = PetscInfo2(eps,"Breakdown in Arnoldi method (it=%D norm=%g)\n",eps->its,(double)beta);CHKERRQ(ierr);
      ierr = EPSGetStartVector(eps,k,&breakdown);CHKERRQ(ierr);
      if (breakdown) {
        eps->reason = EPS_DIVERGED_BREAKDOWN;
        ierr = PetscInfo(eps,"Unable to generate more start vectors\n");CHKERRQ(ierr);
      }
    }
    eps->nconv = k;
    ierr = EPSMonitor(eps,eps->its,eps->nconv,eps->eigr,eps->eigi,eps->errest,nv);CHKERRQ(ierr);
  }

  /* truncate Schur decomposition and change the state to raw so that
     DSVectors() computes eigenvectors from scratch */
  ierr = DSSetDimensions(eps->ds,eps->nconv,0,0,0);CHKERRQ(ierr);
  ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_Arnoldi(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      set,val;
  EPS_ARNOLDI    *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS Arnoldi Options");CHKERRQ(ierr);

    ierr = PetscOptionsBool("-eps_arnoldi_delayed","Use delayed reorthogonalization","EPSArnoldiSetDelayed",arnoldi->delayed,&val,&set);CHKERRQ(ierr);
    if (set) { ierr = EPSArnoldiSetDelayed(eps,val);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSArnoldiSetDelayed_Arnoldi(EPS eps,PetscBool delayed)
{
  EPS_ARNOLDI *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  arnoldi->delayed = delayed;
  PetscFunctionReturn(0);
}

/*@
   EPSArnoldiSetDelayed - Activates or deactivates delayed reorthogonalization
   in the Arnoldi iteration.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  delayed - boolean flag

   Options Database Key:
.  -eps_arnoldi_delayed - Activates delayed reorthogonalization in Arnoldi

   Note:
   Delayed reorthogonalization is an aggressive optimization for the Arnoldi
   eigensolver than may provide better scalability, but sometimes makes the
   solver converge less than the default algorithm.

   Level: advanced

.seealso: EPSArnoldiGetDelayed()
@*/
PetscErrorCode EPSArnoldiSetDelayed(EPS eps,PetscBool delayed)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,delayed,2);
  ierr = PetscTryMethod(eps,"EPSArnoldiSetDelayed_C",(EPS,PetscBool),(eps,delayed));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSArnoldiGetDelayed_Arnoldi(EPS eps,PetscBool *delayed)
{
  EPS_ARNOLDI *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  *delayed = arnoldi->delayed;
  PetscFunctionReturn(0);
}

/*@
   EPSArnoldiGetDelayed - Gets the type of reorthogonalization used during the Arnoldi
   iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Input Parameter:
.  delayed - boolean flag indicating if delayed reorthogonalization has been enabled

   Level: advanced

.seealso: EPSArnoldiSetDelayed()
@*/
PetscErrorCode EPSArnoldiGetDelayed(EPS eps,PetscBool *delayed)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(delayed,2);
  ierr = PetscUseMethod(eps,"EPSArnoldiGetDelayed_C",(EPS,PetscBool*),(eps,delayed));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_Arnoldi(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSArnoldiSetDelayed_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSArnoldiGetDelayed_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_Arnoldi(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isascii;
  EPS_ARNOLDI    *arnoldi = (EPS_ARNOLDI*)eps->data;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii && arnoldi->delayed) {
    ierr = PetscViewerASCIIPrintf(viewer,"  using delayed reorthogonalization\n");CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_Arnoldi(EPS eps)
{
  EPS_ARNOLDI    *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&ctx);CHKERRQ(ierr);
  eps->data = (void*)ctx;

  eps->useds = PETSC_TRUE;

  eps->ops->solve          = EPSSolve_Arnoldi;
  eps->ops->setup          = EPSSetUp_Arnoldi;
  eps->ops->setfromoptions = EPSSetFromOptions_Arnoldi;
  eps->ops->destroy        = EPSDestroy_Arnoldi;
  eps->ops->view           = EPSView_Arnoldi;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_Schur;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSArnoldiSetDelayed_C",EPSArnoldiSetDelayed_Arnoldi);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSArnoldiGetDelayed_C",EPSArnoldiGetDelayed_Arnoldi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

