/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "rii"

   Method: Residual inverse iteration

   Algorithm:

       Simple residual inverse iteration with varying shift.

   References:

       [1] A. Neumaier, "Residual inverse iteration for the nonlinear
           eigenvalue problem", SIAM J. Numer. Anal. 22(5):914-923, 1985.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <../src/nep/impls/nepdefl.h>

typedef struct {
  PetscInt  max_inner_it;     /* maximum number of Newton iterations */
  PetscInt  lag;              /* interval to rebuild preconditioner */
  PetscBool cctol;            /* constant correction tolerance */
  PetscBool herm;             /* whether the Hermitian version of the scalar equation must be used */
  PetscReal deftol;           /* tolerance for the deflation (threshold) */
  KSP       ksp;              /* linear solver object */
} NEP_RII;

PetscErrorCode NEPSetUp_RII(NEP nep)
{
  PetscErrorCode ierr;
  PetscBool      istrivial;

  PetscFunctionBegin;
  if (nep->ncv!=PETSC_DEFAULT) { ierr = PetscInfo(nep,"Setting ncv = nev, ignoring user-provided value\n");CHKERRQ(ierr); }
  nep->ncv = nep->nev;
  if (nep->mpd!=PETSC_DEFAULT) { ierr = PetscInfo(nep,"Setting mpd = nev, ignoring user-provided value\n");CHKERRQ(ierr); }
  nep->mpd = nep->nev;
  if (nep->max_it==PETSC_DEFAULT) nep->max_it = PetscMax(5000,2*nep->n/nep->ncv);
  if (nep->which && nep->which!=NEP_TARGET_MAGNITUDE) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Wrong value of which");

  ierr = RGIsTrivial(nep->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver does not support region filtering");

  ierr = NEPAllocateSolution(nep,0);CHKERRQ(ierr);
  ierr = NEPSetWorkVecs(nep,2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSolve_RII(NEP nep)
{
  PetscErrorCode     ierr;
  NEP_RII            *ctx = (NEP_RII*)nep->data;
  Mat                T,Tp,H;
  Vec                uu,u,r,delta,t;
  PetscScalar        lambda,lambda2,sigma,a1,a2,corr,*Hp,*Ap;
  PetscReal          nrm,resnorm=1.0,ktol=0.1,perr,rtol;
  PetscBool          skip=PETSC_FALSE,lock=PETSC_FALSE;
  PetscInt           inner_its,its=0,ldh,ldds,i,j;
  NEP_EXT_OP         extop=NULL;
  KSPConvergedReason kspreason;

  PetscFunctionBegin;
  /* get initial approximation of eigenvalue and eigenvector */
  ierr = NEPGetDefaultShift(nep,&sigma);CHKERRQ(ierr);
  lambda = sigma;
  if (!nep->nini) {
    ierr = BVSetRandomColumn(nep->V,0);CHKERRQ(ierr);
    ierr = BVNormColumn(nep->V,0,NORM_2,&nrm);CHKERRQ(ierr);
    ierr = BVScaleColumn(nep->V,0,1.0/nrm);CHKERRQ(ierr);
  }
  if (!ctx->ksp) { ierr = NEPRIIGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
  ierr = NEPDeflationInitialize(nep,nep->V,ctx->ksp,PETSC_FALSE,nep->nev,&extop);CHKERRQ(ierr);
  ierr = NEPDeflationCreateVec(extop,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&delta);CHKERRQ(ierr);
  ierr = BVGetColumn(nep->V,0,&uu);CHKERRQ(ierr);
  ierr = NEPDeflationCopyToExtendedVec(extop,uu,NULL,u,PETSC_FALSE);CHKERRQ(ierr);
  ierr = BVRestoreColumn(nep->V,0,&uu);CHKERRQ(ierr);

  /* prepare linear solver */
  ierr = NEPDeflationSolveSetUp(extop,sigma);CHKERRQ(ierr);
  ierr = KSPGetTolerances(ctx->ksp,&rtol,NULL,NULL,NULL);CHKERRQ(ierr);

  ierr = VecCopy(u,r);CHKERRQ(ierr);
  ierr = NEPDeflationFunctionSolve(extop,r,u);CHKERRQ(ierr);
  ierr = VecNorm(u,NORM_2,&nrm);CHKERRQ(ierr);
  ierr = VecScale(u,1.0/nrm);CHKERRQ(ierr);

  /* Restart loop */
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    its++;

    /* Use Newton's method to compute nonlinear Rayleigh functional. Current eigenvalue
       estimate as starting value. */
    inner_its=0;
    lambda2 = lambda;
    do {
      ierr = NEPDeflationComputeFunction(extop,lambda,&T);CHKERRQ(ierr);
      ierr = MatMult(T,u,r);CHKERRQ(ierr);
      if (!ctx->herm) {
        ierr = NEPDeflationFunctionSolve(extop,r,delta);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(ctx->ksp,&kspreason);CHKERRQ(ierr);
        if (kspreason<0) {
          ierr = PetscInfo1(nep,"iter=%D, linear solve failed\n",nep->its);CHKERRQ(ierr);
        }
        t = delta;
      } else t = r;
      ierr = VecDot(t,u,&a1);CHKERRQ(ierr);
      ierr = NEPDeflationComputeJacobian(extop,lambda,&Tp);CHKERRQ(ierr);
      ierr = MatMult(Tp,u,r);CHKERRQ(ierr);
      if (!ctx->herm) {
        ierr = NEPDeflationFunctionSolve(extop,r,delta);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(ctx->ksp,&kspreason);CHKERRQ(ierr);
        if (kspreason<0) {
          ierr = PetscInfo1(nep,"iter=%D, linear solve failed\n",nep->its);CHKERRQ(ierr);
        }
        t = delta;
      } else t = r;
      ierr = VecDot(t,u,&a2);CHKERRQ(ierr);
      corr = a1/a2;
      lambda = lambda - corr;
      inner_its++;
    } while (PetscAbsScalar(corr)/PetscAbsScalar(lambda)>PETSC_SQRT_MACHINE_EPSILON && inner_its<ctx->max_inner_it);

    /* form residual,  r = T(lambda)*u */
    ierr = NEPDeflationComputeFunction(extop,lambda,&T);CHKERRQ(ierr);
    ierr = MatMult(T,u,r);CHKERRQ(ierr);

    /* convergence test */
    perr = nep->errest[nep->nconv];
    ierr = VecNorm(r,NORM_2,&resnorm);CHKERRQ(ierr);
    ierr = (*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx);CHKERRQ(ierr);
    nep->eigr[nep->nconv] = lambda;
    if (its>1 && (nep->errest[nep->nconv]<=nep->tol || nep->errest[nep->nconv]<=ctx->deftol)) {
      if (nep->errest[nep->nconv]<=ctx->deftol && !extop->ref && nep->nconv) {
        ierr = NEPDeflationExtractEigenpair(extop,nep->nconv,u,lambda,nep->ds);CHKERRQ(ierr);
        ierr = NEPDeflationSetRefine(extop,PETSC_TRUE);CHKERRQ(ierr);
        ierr = MatMult(T,u,r);CHKERRQ(ierr);
        ierr = VecNorm(r,NORM_2,&resnorm);CHKERRQ(ierr);
        ierr = (*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx);CHKERRQ(ierr);
        if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
      } else if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
    }
    if (lock) {
      ierr = NEPDeflationSetRefine(extop,PETSC_FALSE);CHKERRQ(ierr);
      nep->nconv = nep->nconv + 1;
      ierr = NEPDeflationLocking(extop,u,lambda);CHKERRQ(ierr);
      nep->its += its;
      skip = PETSC_TRUE;
      lock = PETSC_FALSE;
    }
    ierr = (*nep->stopping)(nep,nep->its+its,nep->max_it,nep->nconv,nep->nev,&nep->reason,nep->stoppingctx);CHKERRQ(ierr);
    if (!skip || nep->reason>0) {
      ierr = NEPMonitor(nep,nep->its+its,nep->nconv,nep->eigr,nep->eigi,nep->errest,(nep->reason>0)?nep->nconv:nep->nconv+1);CHKERRQ(ierr);
    }

    if (nep->reason == NEP_CONVERGED_ITERATING) {
      if (!skip) {
        /* update preconditioner and set adaptive tolerance */
        if (ctx->lag && !(its%ctx->lag) && its>=2*ctx->lag && perr && nep->errest[nep->nconv]>.5*perr) {
          ierr = NEPDeflationSolveSetUp(extop,lambda2);CHKERRQ(ierr);
        }
        if (!ctx->cctol) {
          ktol = PetscMax(ktol/2.0,rtol);
          ierr = KSPSetTolerances(ctx->ksp,ktol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
        }

        /* eigenvector correction: delta = T(sigma)\r */
        ierr = NEPDeflationFunctionSolve(extop,r,delta);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(ctx->ksp,&kspreason);CHKERRQ(ierr);
        if (kspreason<0) {
          ierr = PetscInfo1(nep,"iter=%D, linear solve failed, stopping solve\n",nep->its);CHKERRQ(ierr);
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }

        /* update eigenvector: u = u - delta */
        ierr = VecAXPY(u,-1.0,delta);CHKERRQ(ierr);

        /* normalize eigenvector */
        ierr = VecNormalize(u,NULL);CHKERRQ(ierr);
      } else {
        its = -1;
        ierr = NEPDeflationSetRandomVec(extop,u);CHKERRQ(ierr);
        ierr = NEPDeflationSolveSetUp(extop,sigma);CHKERRQ(ierr);
        ierr = VecCopy(u,r);CHKERRQ(ierr);
        ierr = NEPDeflationFunctionSolve(extop,r,u);CHKERRQ(ierr);
        ierr = VecNorm(u,NORM_2,&nrm);CHKERRQ(ierr);
        ierr = VecScale(u,1.0/nrm);CHKERRQ(ierr);
        lambda = sigma;
        skip = PETSC_FALSE;
      }
    }
  }
  ierr = NEPDeflationGetInvariantPair(extop,NULL,&H);CHKERRQ(ierr);
  ierr = MatGetSize(H,NULL,&ldh);CHKERRQ(ierr);
  ierr = DSSetType(nep->ds,DSNHEP);CHKERRQ(ierr);
  ierr = DSReset(nep->ds);CHKERRQ(ierr);
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
  ierr = VecDestroy(&delta);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSetFromOptions_RII(PetscOptionItems *PetscOptionsObject,NEP nep)
{
  PetscErrorCode ierr;
  NEP_RII        *ctx = (NEP_RII*)nep->data;
  PetscBool      flg;
  PetscInt       i;
  PetscReal      r;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"NEP RII Options");CHKERRQ(ierr);

    i = 0;
    ierr = PetscOptionsInt("-nep_rii_max_it","Maximum number of Newton iterations for updating Rayleigh functional","NEPRIISetMaximumIterations",ctx->max_inner_it,&i,&flg);CHKERRQ(ierr);
    if (flg) { ierr = NEPRIISetMaximumIterations(nep,i);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-nep_rii_const_correction_tol","Constant correction tolerance for the linear solver","NEPRIISetConstCorrectionTol",ctx->cctol,&ctx->cctol,NULL);CHKERRQ(ierr);

    ierr = PetscOptionsBool("-nep_rii_hermitian","Use Hermitian version of the scalar nonlinear equation","NEPRIISetHermitian",ctx->herm,&ctx->herm,NULL);CHKERRQ(ierr);

    i = 0;
    ierr = PetscOptionsInt("-nep_rii_lag_preconditioner","Interval to rebuild preconditioner","NEPRIISetLagPreconditioner",ctx->lag,&i,&flg);CHKERRQ(ierr);
    if (flg) { ierr = NEPRIISetLagPreconditioner(nep,i);CHKERRQ(ierr); }

    r = 0.0;
    ierr = PetscOptionsReal("-nep_rii_deflation_threshold","Tolerance used as a threshold for including deflated eigenpairs","NEPRIISetDeflationThreshold",ctx->deftol,&r,&flg);CHKERRQ(ierr);
    if (flg) { ierr = NEPRIISetDeflationThreshold(nep,r);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);

  if (!ctx->ksp) { ierr = NEPRIIGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
  ierr = KSPSetFromOptions(ctx->ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIISetMaximumIterations_RII(NEP nep,PetscInt its)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  if (its==PETSC_DEFAULT) ctx->max_inner_it = 10;
  else {
    if (its<=0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Number of iterations must be >0");
    ctx->max_inner_it = its;
  }
  PetscFunctionReturn(0);
}

/*@
   NEPRIISetMaximumIterations - Sets the maximum number of inner iterations to be
   used in the RII solver. These are the Newton iterations related to the computation
   of the nonlinear Rayleigh functional.

   Logically Collective on nep

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  its - maximum inner iterations

   Level: advanced

.seealso: NEPRIIGetMaximumIterations()
@*/
PetscErrorCode NEPRIISetMaximumIterations(NEP nep,PetscInt its)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,its,2);
  ierr = PetscTryMethod(nep,"NEPRIISetMaximumIterations_C",(NEP,PetscInt),(nep,its));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIIGetMaximumIterations_RII(NEP nep,PetscInt *its)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *its = ctx->max_inner_it;
  PetscFunctionReturn(0);
}

/*@
   NEPRIIGetMaximumIterations - Gets the maximum number of inner iterations of RII.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  its - maximum inner iterations

   Level: advanced

.seealso: NEPRIISetMaximumIterations()
@*/
PetscErrorCode NEPRIIGetMaximumIterations(NEP nep,PetscInt *its)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidIntPointer(its,2);
  ierr = PetscUseMethod(nep,"NEPRIIGetMaximumIterations_C",(NEP,PetscInt*),(nep,its));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIISetLagPreconditioner_RII(NEP nep,PetscInt lag)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  if (lag<0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Lag must be non-negative");
  ctx->lag = lag;
  PetscFunctionReturn(0);
}

/*@
   NEPRIISetLagPreconditioner - Determines when the preconditioner is rebuilt in the
   nonlinear solve.

   Logically Collective on nep

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  lag - 0 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is
          computed within the nonlinear iteration, 2 means every second time
          the Jacobian is built, etc.

   Options Database Keys:
.  -nep_rii_lag_preconditioner <lag>

   Notes:
   The default is 1.
   The preconditioner is ALWAYS built in the first iteration of a nonlinear solve.

   Level: intermediate

.seealso: NEPRIIGetLagPreconditioner()
@*/
PetscErrorCode NEPRIISetLagPreconditioner(NEP nep,PetscInt lag)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,lag,2);
  ierr = PetscTryMethod(nep,"NEPRIISetLagPreconditioner_C",(NEP,PetscInt),(nep,lag));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIIGetLagPreconditioner_RII(NEP nep,PetscInt *lag)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *lag = ctx->lag;
  PetscFunctionReturn(0);
}

/*@
   NEPRIIGetLagPreconditioner - Indicates how often the preconditioner is rebuilt.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  lag - the lag parameter

   Level: intermediate

.seealso: NEPRIISetLagPreconditioner()
@*/
PetscErrorCode NEPRIIGetLagPreconditioner(NEP nep,PetscInt *lag)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidIntPointer(lag,2);
  ierr = PetscUseMethod(nep,"NEPRIIGetLagPreconditioner_C",(NEP,PetscInt*),(nep,lag));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIISetConstCorrectionTol_RII(NEP nep,PetscBool cct)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  ctx->cctol = cct;
  PetscFunctionReturn(0);
}

/*@
   NEPRIISetConstCorrectionTol - Sets a flag to keep the tolerance used
   in the linear solver constant.

   Logically Collective on nep

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  cct - a boolean value

   Options Database Keys:
.  -nep_rii_const_correction_tol <bool> - set the boolean flag

   Notes:
   By default, an exponentially decreasing tolerance is set in the KSP used
   within the nonlinear iteration, so that each Newton iteration requests
   better accuracy than the previous one. The constant correction tolerance
   flag stops this behaviour.

   Level: intermediate

.seealso: NEPRIIGetConstCorrectionTol()
@*/
PetscErrorCode NEPRIISetConstCorrectionTol(NEP nep,PetscBool cct)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,cct,2);
  ierr = PetscTryMethod(nep,"NEPRIISetConstCorrectionTol_C",(NEP,PetscBool),(nep,cct));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIIGetConstCorrectionTol_RII(NEP nep,PetscBool *cct)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *cct = ctx->cctol;
  PetscFunctionReturn(0);
}

/*@
   NEPRIIGetConstCorrectionTol - Returns the constant tolerance flag.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  cct - the value of the constant tolerance flag

   Level: intermediate

.seealso: NEPRIISetConstCorrectionTol()
@*/
PetscErrorCode NEPRIIGetConstCorrectionTol(NEP nep,PetscBool *cct)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidBoolPointer(cct,2);
  ierr = PetscUseMethod(nep,"NEPRIIGetConstCorrectionTol_C",(NEP,PetscBool*),(nep,cct));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIISetHermitian_RII(NEP nep,PetscBool herm)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  ctx->herm = herm;
  PetscFunctionReturn(0);
}

/*@
   NEPRIISetHermitian - Sets a flag to indicate if the Hermitian version of the
   scalar nonlinear equation must be used by the solver.

   Logically Collective on nep

   Input Parameters:
+  nep  - nonlinear eigenvalue solver
-  herm - a boolean value

   Options Database Keys:
.  -nep_rii_hermitian <bool> - set the boolean flag

   Notes:
   By default, the scalar nonlinear equation x'*inv(T(sigma))*T(z)*x=0 is solved
   at each step of the nonlinear iteration. When this flag is set the simpler
   form x'*T(z)*x=0 is used, which is supposed to be valid only for Hermitian
   problems.

   Level: intermediate

.seealso: NEPRIIGetHermitian()
@*/
PetscErrorCode NEPRIISetHermitian(NEP nep,PetscBool herm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,herm,2);
  ierr = PetscTryMethod(nep,"NEPRIISetHermitian_C",(NEP,PetscBool),(nep,herm));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIIGetHermitian_RII(NEP nep,PetscBool *herm)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *herm = ctx->herm;
  PetscFunctionReturn(0);
}

/*@
   NEPRIIGetHermitian - Returns the flag about using the Hermitian version of
   the scalar nonlinear equation.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  herm - the value of the hermitian flag

   Level: intermediate

.seealso: NEPRIISetHermitian()
@*/
PetscErrorCode NEPRIIGetHermitian(NEP nep,PetscBool *herm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidBoolPointer(herm,2);
  ierr = PetscUseMethod(nep,"NEPRIIGetHermitian_C",(NEP,PetscBool*),(nep,herm));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIISetDeflationThreshold_RII(NEP nep,PetscReal deftol)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  ctx->deftol = deftol;
  PetscFunctionReturn(0);
}

/*@
   NEPRIISetDeflationThreshold - Sets the threshold value used to switch between
   deflated and non-deflated iteration.

   Logically Collective on nep

   Input Parameters:
+  nep    - nonlinear eigenvalue solver
-  deftol - the threshold value

   Options Database Keys:
.  -nep_rii_deflation_threshold <deftol> - set the threshold

   Notes:
   Normally, the solver iterates on the extended problem in order to deflate
   previously converged eigenpairs. If this threshold is set to a nonzero value,
   then once the residual error is below this threshold the solver will
   continue the iteration without deflation. The intention is to be able to
   improve the current eigenpair further, despite having previous eigenpairs
   with somewhat bad precision.

   Level: advanced

.seealso: NEPRIIGetDeflationThreshold()
@*/
PetscErrorCode NEPRIISetDeflationThreshold(NEP nep,PetscReal deftol)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,deftol,2);
  ierr = PetscTryMethod(nep,"NEPRIISetDeflationThreshold_C",(NEP,PetscReal),(nep,deftol));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIIGetDeflationThreshold_RII(NEP nep,PetscReal *deftol)
{
  NEP_RII *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  *deftol = ctx->deftol;
  PetscFunctionReturn(0);
}

/*@
   NEPRIIGetDeflationThreshold - Returns the threshold value that controls deflation.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  deftol - the threshold

   Level: advanced

.seealso: NEPRIISetDeflationThreshold()
@*/
PetscErrorCode NEPRIIGetDeflationThreshold(NEP nep,PetscReal *deftol)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidBoolPointer(deftol,2);
  ierr = PetscUseMethod(nep,"NEPRIIGetDeflationThreshold_C",(NEP,PetscReal*),(nep,deftol));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIISetKSP_RII(NEP nep,KSP ksp)
{
  PetscErrorCode ierr;
  NEP_RII        *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)ksp);CHKERRQ(ierr);
  ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
  ctx->ksp = ksp;
  ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->ksp);CHKERRQ(ierr);
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   NEPRIISetKSP - Associate a linear solver object (KSP) to the nonlinear
   eigenvalue solver.

   Collective on nep

   Input Parameters:
+  nep - eigenvalue solver
-  ksp - the linear solver object

   Level: advanced

.seealso: NEPRIIGetKSP()
@*/
PetscErrorCode NEPRIISetKSP(NEP nep,KSP ksp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,2);
  PetscCheckSameComm(nep,1,ksp,2);
  ierr = PetscTryMethod(nep,"NEPRIISetKSP_C",(NEP,KSP),(nep,ksp));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPRIIGetKSP_RII(NEP nep,KSP *ksp)
{
  PetscErrorCode ierr;
  NEP_RII        *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  if (!ctx->ksp) {
    ierr = KSPCreate(PetscObjectComm((PetscObject)nep),&ctx->ksp);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)ctx->ksp,(PetscObject)nep,1);CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(ctx->ksp,((PetscObject)nep)->prefix);CHKERRQ(ierr);
    ierr = KSPAppendOptionsPrefix(ctx->ksp,"nep_rii_");CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->ksp);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)ctx->ksp,((PetscObject)nep)->options);CHKERRQ(ierr);
    ierr = KSPSetErrorIfNotConverged(ctx->ksp,PETSC_TRUE);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ctx->ksp,SLEPC_DEFAULT_TOL,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  }
  *ksp = ctx->ksp;
  PetscFunctionReturn(0);
}

/*@
   NEPRIIGetKSP - Retrieve the linear solver object (KSP) associated with
   the nonlinear eigenvalue solver.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  ksp - the linear solver object

   Level: advanced

.seealso: NEPRIISetKSP()
@*/
PetscErrorCode NEPRIIGetKSP(NEP nep,KSP *ksp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(ksp,2);
  ierr = PetscUseMethod(nep,"NEPRIIGetKSP_C",(NEP,KSP*),(nep,ksp));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPView_RII(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  NEP_RII        *ctx = (NEP_RII*)nep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum number of inner iterations: %D\n",ctx->max_inner_it);CHKERRQ(ierr);
    if (ctx->cctol) {
      ierr = PetscViewerASCIIPrintf(viewer,"  using a constant tolerance for the linear solver\n");CHKERRQ(ierr);
    }
    if (ctx->herm) {
      ierr = PetscViewerASCIIPrintf(viewer,"  using the Hermitian version of the scalar nonlinear equation\n");CHKERRQ(ierr);
    }
    if (ctx->lag) {
      ierr = PetscViewerASCIIPrintf(viewer,"  updating the preconditioner every %D iterations\n",ctx->lag);CHKERRQ(ierr);
    }
    if (ctx->deftol) {
      ierr = PetscViewerASCIIPrintf(viewer,"  deflation threshold: %g\n",(double)ctx->deftol);CHKERRQ(ierr);
    }
    if (!ctx->ksp) { ierr = NEPRIIGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = KSPView(ctx->ksp,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPReset_RII(NEP nep)
{
  PetscErrorCode ierr;
  NEP_RII        *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  ierr = KSPReset(ctx->ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDestroy_RII(NEP nep)
{
  PetscErrorCode ierr;
  NEP_RII        *ctx = (NEP_RII*)nep->data;

  PetscFunctionBegin;
  ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
  ierr = PetscFree(nep->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetMaximumIterations_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetMaximumIterations_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetLagPreconditioner_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetLagPreconditioner_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetConstCorrectionTol_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetConstCorrectionTol_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetHermitian_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetHermitian_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetDeflationThreshold_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetDeflationThreshold_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetKSP_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetKSP_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_RII(NEP nep)
{
  PetscErrorCode ierr;
  NEP_RII        *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(nep,&ctx);CHKERRQ(ierr);
  nep->data = (void*)ctx;
  ctx->max_inner_it = 10;
  ctx->lag          = 1;
  ctx->cctol        = PETSC_FALSE;
  ctx->herm         = PETSC_FALSE;
  ctx->deftol       = 0.0;

  nep->useds = PETSC_TRUE;

  nep->ops->solve          = NEPSolve_RII;
  nep->ops->setup          = NEPSetUp_RII;
  nep->ops->setfromoptions = NEPSetFromOptions_RII;
  nep->ops->reset          = NEPReset_RII;
  nep->ops->destroy        = NEPDestroy_RII;
  nep->ops->view           = NEPView_RII;
  nep->ops->computevectors = NEPComputeVectors_Schur;

  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetMaximumIterations_C",NEPRIISetMaximumIterations_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetMaximumIterations_C",NEPRIIGetMaximumIterations_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetLagPreconditioner_C",NEPRIISetLagPreconditioner_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetLagPreconditioner_C",NEPRIIGetLagPreconditioner_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetConstCorrectionTol_C",NEPRIISetConstCorrectionTol_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetConstCorrectionTol_C",NEPRIIGetConstCorrectionTol_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetHermitian_C",NEPRIISetHermitian_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetHermitian_C",NEPRIIGetHermitian_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetDeflationThreshold_C",NEPRIISetDeflationThreshold_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetDeflationThreshold_C",NEPRIIGetDeflationThreshold_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIISetKSP_C",NEPRIISetKSP_RII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPRIIGetKSP_C",NEPRIIGetKSP_RII);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

