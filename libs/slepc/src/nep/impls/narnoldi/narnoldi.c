/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "narnoldi"

   Method: Nonlinear Arnoldi

   Algorithm:

       Arnoldi for nonlinear eigenproblems.

   References:

       [1] H. Voss, "An Arnoldi method for nonlinear eigenvalue problems",
           BIT 44:387-401, 2004.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <../src/nep/impls/nepdefl.h>

typedef struct {
  PetscInt lag;             /* interval to rebuild preconditioner */
  KSP      ksp;             /* linear solver object */
} NEP_NARNOLDI;

PetscErrorCode NEPSetUp_NArnoldi(NEP nep)
{
  PetscErrorCode ierr;
  PetscBool      istrivial;

  PetscFunctionBegin;
  ierr = NEPSetDimensions_Default(nep,nep->nev,&nep->ncv,&nep->mpd);CHKERRQ(ierr);
  if (nep->ncv>nep->nev+nep->mpd) SETERRQ(PetscObjectComm((PetscObject)nep),1,"The value of ncv must not be larger than nev+mpd");
  if (nep->max_it==PETSC_DEFAULT) nep->max_it = nep->nev*nep->ncv;
  if (nep->which && nep->which!=NEP_TARGET_MAGNITUDE) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Wrong value of which");
  if (nep->fui!=NEP_USER_INTERFACE_SPLIT) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"NARNOLDI only available for split operator");

  ierr = RGIsTrivial(nep->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver does not support region filtering");

  ierr = NEPAllocateSolution(nep,0);CHKERRQ(ierr);
  ierr = NEPSetWorkVecs(nep,3);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSolve_NArnoldi(NEP nep)
{
  PetscErrorCode     ierr;
  NEP_NARNOLDI       *ctx = (NEP_NARNOLDI*)nep->data;
  Mat                T,H;
  Vec                f,r,u,uu;
  PetscScalar        *X,lambda=0.0,lambda2=0.0,*eigr,*Hp,*Ap,sigma;
  PetscReal          beta,resnorm=0.0,nrm,perr=0.0;
  PetscInt           n,i,j,ldds,ldh;
  PetscBool          breakdown,skip=PETSC_FALSE;
  BV                 Vext;
  DS                 ds;
  NEP_EXT_OP         extop=NULL;
  SlepcSC            sc;
  KSPConvergedReason kspreason;

  PetscFunctionBegin;
  /* get initial space and shift */
  ierr = NEPGetDefaultShift(nep,&sigma);CHKERRQ(ierr);
  if (!nep->nini) {
    ierr = BVSetRandomColumn(nep->V,0);CHKERRQ(ierr);
    ierr = BVNormColumn(nep->V,0,NORM_2,&nrm);CHKERRQ(ierr);
    ierr = BVScaleColumn(nep->V,0,1.0/nrm);CHKERRQ(ierr);
    n = 1;
  } else n = nep->nini;

  if (!ctx->ksp) { ierr = NEPNArnoldiGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
  ierr = NEPDeflationInitialize(nep,nep->V,ctx->ksp,PETSC_FALSE,nep->nev,&extop);CHKERRQ(ierr);
  ierr = NEPDeflationCreateBV(extop,nep->ncv,&Vext);CHKERRQ(ierr);

  /* prepare linear solver */
  ierr = NEPDeflationSolveSetUp(extop,sigma);CHKERRQ(ierr);

  ierr = BVGetColumn(Vext,0,&f);CHKERRQ(ierr);
  ierr = VecDuplicate(f,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(f,&u);CHKERRQ(ierr);
  ierr = BVGetColumn(nep->V,0,&uu);CHKERRQ(ierr);
  ierr = NEPDeflationCopyToExtendedVec(extop,uu,NULL,f,PETSC_FALSE);CHKERRQ(ierr);
  ierr = BVRestoreColumn(nep->V,0,&uu);CHKERRQ(ierr);
  ierr = VecCopy(f,r);CHKERRQ(ierr);
  ierr = NEPDeflationFunctionSolve(extop,r,f);CHKERRQ(ierr);
  ierr = VecNorm(f,NORM_2,&nrm);CHKERRQ(ierr);
  ierr = VecScale(f,1.0/nrm);CHKERRQ(ierr);
  ierr = BVRestoreColumn(Vext,0,&f);CHKERRQ(ierr);

  ierr = DSCreate(PetscObjectComm((PetscObject)nep),&ds);CHKERRQ(ierr);
  ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ds);CHKERRQ(ierr);
  ierr = DSSetType(ds,DSNEP);CHKERRQ(ierr);
  ierr = DSNEPSetFN(ds,nep->nt,nep->f);CHKERRQ(ierr);
  ierr = DSAllocate(ds,nep->ncv);CHKERRQ(ierr);
  ierr = DSGetSlepcSC(ds,&sc);CHKERRQ(ierr);
  sc->comparison    = nep->sc->comparison;
  sc->comparisonctx = nep->sc->comparisonctx;
  ierr = DSSetFromOptions(ds);CHKERRQ(ierr);

  /* build projected matrices for initial space */
  ierr = DSSetDimensions(ds,n,0,0,0);CHKERRQ(ierr);
  ierr = NEPDeflationProjectOperator(extop,Vext,ds,0,n);CHKERRQ(ierr);

  ierr = PetscMalloc1(nep->ncv,&eigr);CHKERRQ(ierr);

  /* Restart loop */
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    nep->its++;

    /* solve projected problem */
    ierr = DSSetDimensions(ds,n,0,0,0);CHKERRQ(ierr);
    ierr = DSSetState(ds,DS_STATE_RAW);CHKERRQ(ierr);
    ierr = DSSolve(ds,eigr,NULL);CHKERRQ(ierr);
    ierr = DSSynchronize(ds,eigr,NULL);CHKERRQ(ierr);
    if (nep->its>1) lambda2 = lambda;
    lambda = eigr[0];
    nep->eigr[nep->nconv] = lambda;

    /* compute Ritz vector, x = V*s */
    ierr = DSGetArray(ds,DS_MAT_X,&X);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(Vext,0,n);CHKERRQ(ierr);
    ierr = BVMultVec(Vext,1.0,0.0,u,X);CHKERRQ(ierr);
    ierr = DSRestoreArray(ds,DS_MAT_X,&X);CHKERRQ(ierr);

    /* compute the residual, r = T(lambda)*x */
    ierr = NEPDeflationComputeFunction(extop,lambda,&T);CHKERRQ(ierr);
    ierr = MatMult(T,u,r);CHKERRQ(ierr);

    /* convergence test */
    ierr = VecNorm(r,NORM_2,&resnorm);CHKERRQ(ierr);
    if (nep->its>1) perr=nep->errest[nep->nconv];
    ierr = (*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx);CHKERRQ(ierr);
    if (nep->errest[nep->nconv]<=nep->tol) {
      nep->nconv = nep->nconv + 1;
      ierr = NEPDeflationLocking(extop,u,lambda);CHKERRQ(ierr);
      skip = PETSC_TRUE;
    }
    ierr = (*nep->stopping)(nep,nep->its,nep->max_it,nep->nconv,nep->nev,&nep->reason,nep->stoppingctx);CHKERRQ(ierr);
    if (!skip || nep->reason>0) {
      ierr = NEPMonitor(nep,nep->its,nep->nconv,nep->eigr,nep->eigi,nep->errest,(nep->reason>0)?nep->nconv:nep->nconv+1);CHKERRQ(ierr);
    }

    if (nep->reason == NEP_CONVERGED_ITERATING) {
      if (!skip) {
        if (n>=nep->ncv) {
          nep->reason = NEP_DIVERGED_SUBSPACE_EXHAUSTED;
          break;
        }
        if (ctx->lag && !(nep->its%ctx->lag) && nep->its>=2*ctx->lag && perr && nep->errest[nep->nconv]>.5*perr) {
          ierr = NEPDeflationSolveSetUp(extop,lambda2);CHKERRQ(ierr);
        }

        /* continuation vector: f = T(sigma)\r */
        ierr = BVGetColumn(Vext,n,&f);CHKERRQ(ierr);
        ierr = NEPDeflationFunctionSolve(extop,r,f);CHKERRQ(ierr);
        ierr = BVRestoreColumn(Vext,n,&f);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(ctx->ksp,&kspreason);CHKERRQ(ierr);
        if (kspreason<0) {
          ierr = PetscInfo1(nep,"iter=%D, linear solve failed, stopping solve\n",nep->its);CHKERRQ(ierr);
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }

        /* orthonormalize */
        ierr = BVOrthonormalizeColumn(Vext,n,PETSC_FALSE,&beta,&breakdown);CHKERRQ(ierr);
        if (breakdown || beta==0.0) {
          ierr = PetscInfo1(nep,"iter=%D, orthogonalization failed, stopping solve\n",nep->its);CHKERRQ(ierr);
          nep->reason = NEP_DIVERGED_BREAKDOWN;
          break;
        }

        /* update projected matrices */
        ierr = DSSetDimensions(ds,n+1,0,0,0);CHKERRQ(ierr);
        ierr = NEPDeflationProjectOperator(extop,Vext,ds,n,n+1);CHKERRQ(ierr);
        n++;
      } else {
        nep->its--;  /* do not count this as a full iteration */
        ierr = BVGetColumn(Vext,0,&f);CHKERRQ(ierr);
        ierr = NEPDeflationSetRandomVec(extop,f);CHKERRQ(ierr);
        ierr = NEPDeflationSolveSetUp(extop,sigma);CHKERRQ(ierr);
        ierr = VecCopy(f,r);CHKERRQ(ierr);
        ierr = NEPDeflationFunctionSolve(extop,r,f);CHKERRQ(ierr);
        ierr = VecNorm(f,NORM_2,&nrm);CHKERRQ(ierr);
        ierr = VecScale(f,1.0/nrm);CHKERRQ(ierr);
        ierr = BVRestoreColumn(Vext,0,&f);CHKERRQ(ierr);
        n = 1;
        ierr = DSSetDimensions(ds,n,0,0,0);CHKERRQ(ierr);
        ierr = NEPDeflationProjectOperator(extop,Vext,ds,n-1,n);CHKERRQ(ierr);
        skip = PETSC_FALSE;
      }
    }
  }

  ierr = NEPDeflationGetInvariantPair(extop,NULL,&H);CHKERRQ(ierr);
  ierr = MatGetSize(H,NULL,&ldh);CHKERRQ(ierr);
  ierr = DSSetType(nep->ds,DSNHEP);CHKERRQ(ierr);
  ierr = DSAllocate(nep->ds,PetscMax(nep->nconv,1));CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(nep->ds,&ldds);CHKERRQ(ierr);
  ierr = MatDenseGetArray(H,&Hp);CHKERRQ(ierr);
  ierr = DSGetArray(nep->ds,DS_MAT_A,&Ap);CHKERRQ(ierr);
  for (j=0;j<nep->nconv;j++)
    for (i=0;i<nep->nconv;i++) Ap[j*ldds+i] = Hp[j*ldh+i];
  ierr = DSRestoreArray(nep->ds,DS_MAT_A,&Ap);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(H,&Hp);CHKERRQ(ierr);
  ierr = MatDestroy(&H);CHKERRQ(ierr);
  ierr = DSSetDimensions(nep->ds,nep->nconv,0,0,nep->nconv);CHKERRQ(ierr);
  ierr = DSSolve(nep->ds,nep->eigr,nep->eigi);CHKERRQ(ierr);
  ierr = NEPDeflationReset(extop);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = BVDestroy(&Vext);CHKERRQ(ierr);
  ierr = DSDestroy(&ds);CHKERRQ(ierr);
  ierr = PetscFree(eigr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNArnoldiSetLagPreconditioner_NArnoldi(NEP nep,PetscInt lag)
{
  NEP_NARNOLDI *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  if (lag<0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Lag must be non-negative");
  ctx->lag = lag;
  PetscFunctionReturn(0);
}

/*@
   NEPNArnoldiSetLagPreconditioner - Determines when the preconditioner is rebuilt in the
   nonlinear solve.

   Logically Collective on nep

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  lag - 0 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is
          computed within the nonlinear iteration, 2 means every second time
          the Jacobian is built, etc.

   Options Database Keys:
.  -nep_narnoldi_lag_preconditioner <lag>

   Notes:
   The default is 1.
   The preconditioner is ALWAYS built in the first iteration of a nonlinear solve.

   Level: intermediate

.seealso: NEPNArnoldiGetLagPreconditioner()
@*/
PetscErrorCode NEPNArnoldiSetLagPreconditioner(NEP nep,PetscInt lag)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,lag,2);
  ierr = PetscTryMethod(nep,"NEPNArnoldiSetLagPreconditioner_C",(NEP,PetscInt),(nep,lag));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNArnoldiGetLagPreconditioner_NArnoldi(NEP nep,PetscInt *lag)
{
  NEP_NARNOLDI *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  *lag = ctx->lag;
  PetscFunctionReturn(0);
}

/*@
   NEPNArnoldiGetLagPreconditioner - Indicates how often the preconditioner is rebuilt.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  lag - the lag parameter

   Level: intermediate

.seealso: NEPNArnoldiSetLagPreconditioner()
@*/
PetscErrorCode NEPNArnoldiGetLagPreconditioner(NEP nep,PetscInt *lag)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidIntPointer(lag,2);
  ierr = PetscUseMethod(nep,"NEPNArnoldiGetLagPreconditioner_C",(NEP,PetscInt*),(nep,lag));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSetFromOptions_NArnoldi(PetscOptionItems *PetscOptionsObject,NEP nep)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBool      flg;
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"NEP N-Arnoldi Options");CHKERRQ(ierr);
    i = 0;
    ierr = PetscOptionsInt("-nep_narnoldi_lag_preconditioner","Interval to rebuild preconditioner","NEPNArnoldiSetLagPreconditioner",ctx->lag,&i,&flg);CHKERRQ(ierr);
    if (flg) { ierr = NEPNArnoldiSetLagPreconditioner(nep,i);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);

  if (!ctx->ksp) { ierr = NEPNArnoldiGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
  ierr = KSPSetFromOptions(ctx->ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNArnoldiSetKSP_NArnoldi(NEP nep,KSP ksp)
{
  PetscErrorCode ierr;
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)ksp);CHKERRQ(ierr);
  ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
  ctx->ksp = ksp;
  ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->ksp);CHKERRQ(ierr);
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   NEPNArnoldiSetKSP - Associate a linear solver object (KSP) to the nonlinear
   eigenvalue solver.

   Collective on nep

   Input Parameters:
+  nep - eigenvalue solver
-  ksp - the linear solver object

   Level: advanced

.seealso: NEPNArnoldiGetKSP()
@*/
PetscErrorCode NEPNArnoldiSetKSP(NEP nep,KSP ksp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,2);
  PetscCheckSameComm(nep,1,ksp,2);
  ierr = PetscTryMethod(nep,"NEPNArnoldiSetKSP_C",(NEP,KSP),(nep,ksp));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNArnoldiGetKSP_NArnoldi(NEP nep,KSP *ksp)
{
  PetscErrorCode ierr;
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  if (!ctx->ksp) {
    ierr = KSPCreate(PetscObjectComm((PetscObject)nep),&ctx->ksp);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)ctx->ksp,(PetscObject)nep,1);CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(ctx->ksp,((PetscObject)nep)->prefix);CHKERRQ(ierr);
    ierr = KSPAppendOptionsPrefix(ctx->ksp,"nep_narnoldi_");CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->ksp);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)ctx->ksp,((PetscObject)nep)->options);CHKERRQ(ierr);
    ierr = KSPSetErrorIfNotConverged(ctx->ksp,PETSC_TRUE);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ctx->ksp,SLEPC_DEFAULT_TOL,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  }
  *ksp = ctx->ksp;
  PetscFunctionReturn(0);
}

/*@
   NEPNArnoldiGetKSP - Retrieve the linear solver object (KSP) associated with
   the nonlinear eigenvalue solver.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  ksp - the linear solver object

   Level: advanced

.seealso: NEPNArnoldiSetKSP()
@*/
PetscErrorCode NEPNArnoldiGetKSP(NEP nep,KSP *ksp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(ksp,2);
  ierr = PetscUseMethod(nep,"NEPNArnoldiGetKSP_C",(NEP,KSP*),(nep,ksp));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPView_NArnoldi(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    if (ctx->lag) {
      ierr = PetscViewerASCIIPrintf(viewer,"  updating the preconditioner every %D iterations\n",ctx->lag);CHKERRQ(ierr);
    }
    if (!ctx->ksp) { ierr = NEPNArnoldiGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = KSPView(ctx->ksp,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPReset_NArnoldi(NEP nep)
{
  PetscErrorCode ierr;
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  ierr = KSPReset(ctx->ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDestroy_NArnoldi(NEP nep)
{
  PetscErrorCode ierr;
  NEP_NARNOLDI   *ctx = (NEP_NARNOLDI*)nep->data;

  PetscFunctionBegin;
  ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
  ierr = PetscFree(nep->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiSetLagPreconditioner_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiGetLagPreconditioner_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiSetKSP_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiGetKSP_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_NArnoldi(NEP nep)
{
  PetscErrorCode ierr;
  NEP_NARNOLDI   *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(nep,&ctx);CHKERRQ(ierr);
  nep->data = (void*)ctx;
  ctx->lag  = 1;

  nep->useds = PETSC_TRUE;

  nep->ops->solve          = NEPSolve_NArnoldi;
  nep->ops->setup          = NEPSetUp_NArnoldi;
  nep->ops->setfromoptions = NEPSetFromOptions_NArnoldi;
  nep->ops->reset          = NEPReset_NArnoldi;
  nep->ops->destroy        = NEPDestroy_NArnoldi;
  nep->ops->view           = NEPView_NArnoldi;
  nep->ops->computevectors = NEPComputeVectors_Schur;

  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiSetLagPreconditioner_C",NEPNArnoldiSetLagPreconditioner_NArnoldi);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiGetLagPreconditioner_C",NEPNArnoldiGetLagPreconditioner_NArnoldi);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiSetKSP_C",NEPNArnoldiSetKSP_NArnoldi);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNArnoldiGetKSP_C",NEPNArnoldiGetKSP_NArnoldi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

