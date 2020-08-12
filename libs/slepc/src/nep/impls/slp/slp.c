/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "slp"

   Method: Succesive linear problems

   Algorithm:

       Newton-type iteration based on first order Taylor approximation.

   References:

       [1] A. Ruhe, "Algorithms for the nonlinear eigenvalue problem", SIAM J.
           Numer. Anal. 10(4):674-689, 1973.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <../src/nep/impls/nepdefl.h>
#include "slp.h"

typedef struct {
  NEP_EXT_OP extop;
  Vec        w;
} NEP_SLP_MATSHELL;

PetscErrorCode NEPSetUp_SLP(NEP nep)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;
  PetscBool      istrivial,flg;
  ST             st;

  PetscFunctionBegin;
  if (nep->ncv!=PETSC_DEFAULT) { ierr = PetscInfo(nep,"Setting ncv = nev, ignoring user-provided value\n");CHKERRQ(ierr); }
  nep->ncv = nep->nev;
  if (nep->mpd!=PETSC_DEFAULT) { ierr = PetscInfo(nep,"Setting mpd = nev, ignoring user-provided value\n");CHKERRQ(ierr); }
  nep->mpd = nep->nev;
  if (nep->ncv>nep->nev+nep->mpd) SETERRQ(PetscObjectComm((PetscObject)nep),1,"The value of ncv must not be larger than nev+mpd");
  if (nep->max_it==PETSC_DEFAULT) nep->max_it = PetscMax(5000,2*nep->n/nep->ncv);
  if (nep->which && nep->which!=NEP_TARGET_MAGNITUDE) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Wrong value of which");

  ierr = RGIsTrivial(nep->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver does not support region filtering");

  if (!ctx->eps) { ierr = NEPSLPGetEPS(nep,&ctx->eps);CHKERRQ(ierr); }
  ierr = EPSGetST(ctx->eps,&st);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompareAny((PetscObject)st,&flg,STSINVERT,STCAYLEY,"");CHKERRQ(ierr);
  if (flg) SETERRQ(PetscObjectComm((PetscObject)nep),1,"SLP does not support spectral transformation");
  ierr = EPSSetDimensions(ctx->eps,1,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(ctx->eps,EPS_LARGEST_MAGNITUDE);CHKERRQ(ierr);
  ierr = EPSSetTolerances(ctx->eps,nep->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL/10.0:nep->tol/10.0,nep->max_it);CHKERRQ(ierr);
  if (nep->twosided) {
    nep->ops->solve = NEPSolve_SLP_Twosided;
    nep->ops->computevectors = NULL;
    if (!ctx->epsts) { ierr = NEPSLPGetEPSLeft(nep,&ctx->epsts);CHKERRQ(ierr); }
    ierr = EPSGetST(ctx->epsts,&st);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompareAny((PetscObject)st,&flg,STSINVERT,STCAYLEY,"");CHKERRQ(ierr);
    if (flg) SETERRQ(PetscObjectComm((PetscObject)nep),1,"SLP does not support spectral transformation");
    ierr = EPSSetDimensions(ctx->epsts,1,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(ctx->epsts,EPS_LARGEST_MAGNITUDE);CHKERRQ(ierr);
    ierr = EPSSetTolerances(ctx->epsts,nep->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL/10.0:nep->tol/10.0,nep->max_it);CHKERRQ(ierr);
  } else {
    nep->ops->solve = NEPSolve_SLP;
    nep->ops->computevectors = NEPComputeVectors_Schur;
  }
  ierr = NEPAllocateSolution(nep,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_SLP(Mat M,Vec x,Vec y)
{
  PetscErrorCode   ierr;
  NEP_SLP_MATSHELL *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&ctx);CHKERRQ(ierr);
  ierr = MatMult(ctx->extop->MJ,x,ctx->w);CHKERRQ(ierr);
  ierr = NEPDeflationFunctionSolve(ctx->extop,ctx->w,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDestroy_SLP(Mat M)
{
  PetscErrorCode   ierr;
  NEP_SLP_MATSHELL *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->w);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_CUDA)
static PetscErrorCode MatCreateVecs_SLP(Mat M,Vec *left,Vec *right)
{
  PetscErrorCode   ierr;
  NEP_SLP_MATSHELL *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&ctx);CHKERRQ(ierr);
  if (right) {
    ierr = VecDuplicate(ctx->w,right);CHKERRQ(ierr);
  }
  if (left) {
    ierr = VecDuplicate(ctx->w,left);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
#endif

static PetscErrorCode NEPSLPSetUpLinearEP(NEP nep,NEP_EXT_OP extop,PetscScalar lambda,Vec u,PetscBool ini)
{
  PetscErrorCode   ierr;
  NEP_SLP          *slpctx = (NEP_SLP*)nep->data;
  Mat              Mshell;
  PetscInt         nloc,mloc;
  NEP_SLP_MATSHELL *shellctx;

  PetscFunctionBegin;
  if (ini) {
    /* Create mat shell */
    ierr = PetscNew(&shellctx);CHKERRQ(ierr);
    shellctx->extop = extop;
    ierr = NEPDeflationCreateVec(extop,&shellctx->w);CHKERRQ(ierr);
    ierr = MatGetLocalSize(nep->function,&mloc,&nloc);CHKERRQ(ierr);
    nloc += extop->szd; mloc += extop->szd;
    ierr = MatCreateShell(PetscObjectComm((PetscObject)nep),nloc,mloc,PETSC_DETERMINE,PETSC_DETERMINE,shellctx,&Mshell);CHKERRQ(ierr);
    ierr = MatShellSetOperation(Mshell,MATOP_MULT,(void(*)(void))MatMult_SLP);CHKERRQ(ierr);
    ierr = MatShellSetOperation(Mshell,MATOP_DESTROY,(void(*)(void))MatDestroy_SLP);CHKERRQ(ierr);
#if defined(PETSC_HAVE_CUDA)
    ierr = MatShellSetOperation(Mshell,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_SLP);CHKERRQ(ierr);
#endif
    ierr = EPSSetOperators(slpctx->eps,Mshell,NULL);CHKERRQ(ierr);
    ierr = MatDestroy(&Mshell);CHKERRQ(ierr);
  }
  ierr = NEPDeflationSolveSetUp(extop,lambda);CHKERRQ(ierr);
  ierr = NEPDeflationComputeJacobian(extop,lambda,NULL);CHKERRQ(ierr);
  ierr = EPSSetInitialSpace(slpctx->eps,1,&u);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSolve_SLP(NEP nep)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;
  Mat            F,H;
  Vec            uu,u,r;
  PetscScalar    sigma,lambda,mu,im,*Hp,*Ap;
  PetscReal      resnorm;
  PetscInt       nconv,ldh,ldds,i,j;
  PetscBool      skip=PETSC_FALSE,lock=PETSC_FALSE;
  NEP_EXT_OP     extop=NULL;    /* Extended operator for deflation */

  PetscFunctionBegin;
  /* get initial approximation of eigenvalue and eigenvector */
  ierr = NEPGetDefaultShift(nep,&sigma);CHKERRQ(ierr);
  if (!nep->nini) {
    ierr = BVSetRandomColumn(nep->V,0);CHKERRQ(ierr);
  }
  lambda = sigma;
  if (!ctx->ksp) { ierr = NEPSLPGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
  ierr = NEPDeflationInitialize(nep,nep->V,ctx->ksp,PETSC_TRUE,nep->nev,&extop);CHKERRQ(ierr);
  ierr = NEPDeflationCreateVec(extop,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&r);CHKERRQ(ierr);
  ierr = BVGetColumn(nep->V,0,&uu);CHKERRQ(ierr);
  ierr = NEPDeflationCopyToExtendedVec(extop,uu,NULL,u,PETSC_FALSE);CHKERRQ(ierr);
  ierr = BVRestoreColumn(nep->V,0,&uu);CHKERRQ(ierr);

  /* Restart loop */
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    nep->its++;

    /* form residual,  r = T(lambda)*u (used in convergence test only) */
    ierr = NEPDeflationComputeFunction(extop,lambda,&F);CHKERRQ(ierr);
    ierr = MatMult(F,u,r);CHKERRQ(ierr);

    /* convergence test */
    ierr = VecNorm(r,NORM_2,&resnorm);CHKERRQ(ierr);
    ierr = (*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx);CHKERRQ(ierr);
    nep->eigr[nep->nconv] = lambda;
    if (nep->errest[nep->nconv]<=nep->tol || nep->errest[nep->nconv]<=ctx->deftol) {
      if (nep->errest[nep->nconv]<=ctx->deftol && !extop->ref && nep->nconv) {
        ierr = NEPDeflationExtractEigenpair(extop,nep->nconv,u,lambda,nep->ds);CHKERRQ(ierr);
        ierr = NEPDeflationSetRefine(extop,PETSC_TRUE);CHKERRQ(ierr);
        ierr = MatMult(F,u,r);CHKERRQ(ierr);
        ierr = VecNorm(r,NORM_2,&resnorm);CHKERRQ(ierr);
        ierr = (*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx);CHKERRQ(ierr);
        if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
      } else if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
    }

    if (lock) {
      ierr = NEPDeflationSetRefine(extop,PETSC_FALSE);CHKERRQ(ierr);
      nep->nconv = nep->nconv + 1;
      skip = PETSC_TRUE;
      lock = PETSC_FALSE;
      ierr = NEPDeflationLocking(extop,u,lambda);CHKERRQ(ierr);
    }
    ierr = (*nep->stopping)(nep,nep->its,nep->max_it,nep->nconv,nep->nev,&nep->reason,nep->stoppingctx);CHKERRQ(ierr);
    if (!skip || nep->reason>0) {
      ierr = NEPMonitor(nep,nep->its,nep->nconv,nep->eigr,nep->eigi,nep->errest,(nep->reason>0)?nep->nconv:nep->nconv+1);CHKERRQ(ierr);
    }

    if (nep->reason == NEP_CONVERGED_ITERATING) {
      if (!skip) {
        /* evaluate T(lambda) and T'(lambda) */
        ierr = NEPSLPSetUpLinearEP(nep,extop,lambda,u,nep->its==1?PETSC_TRUE:PETSC_FALSE);CHKERRQ(ierr);
        /* compute new eigenvalue correction mu and eigenvector approximation u */
        ierr = EPSSolve(ctx->eps);CHKERRQ(ierr);
        ierr = EPSGetConverged(ctx->eps,&nconv);CHKERRQ(ierr);
        if (!nconv) {
          ierr = PetscInfo1(nep,"iter=%D, inner iteration failed, stopping solve\n",nep->its);CHKERRQ(ierr);
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }
        ierr = EPSGetEigenpair(ctx->eps,0,&mu,&im,u,NULL);CHKERRQ(ierr);
        mu = 1.0/mu;
        if (PetscAbsScalar(im)>PETSC_MACHINE_EPSILON) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Complex eigenvalue approximation - not implemented in real scalars");
      } else {
        nep->its--;  /* do not count this as a full iteration */
        /* use second eigenpair computed in previous iteration */
        ierr = EPSGetConverged(ctx->eps,&nconv);CHKERRQ(ierr);
        if (nconv>=2) {
          ierr = EPSGetEigenpair(ctx->eps,1,&mu,&im,u,NULL);CHKERRQ(ierr);
          mu = 1.0/mu;
        } else {
          ierr = NEPDeflationSetRandomVec(extop,u);CHKERRQ(ierr);
          mu = lambda-sigma;
        }
        skip = PETSC_FALSE;
      }
      /* correct eigenvalue */
      lambda = lambda - mu;
    }
  }
  ierr = NEPDeflationGetInvariantPair(extop,NULL,&H);CHKERRQ(ierr);
  ierr = MatGetSize(H,NULL,&ldh);CHKERRQ(ierr);
  ierr = DSReset(nep->ds);CHKERRQ(ierr);
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
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSetFromOptions_SLP(PetscOptionItems *PetscOptionsObject,NEP nep)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;
  PetscBool      flg;
  PetscReal      r;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"NEP SLP Options");CHKERRQ(ierr);

    r = 0.0;
    ierr = PetscOptionsReal("-nep_slp_deflation_threshold","Tolerance used as a threshold for including deflated eigenpairs","NEPSLPSetDeflationThreshold",ctx->deftol,&r,&flg);CHKERRQ(ierr);
    if (flg) { ierr = NEPSLPSetDeflationThreshold(nep,r);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);

  if (!ctx->eps) { ierr = NEPSLPGetEPS(nep,&ctx->eps);CHKERRQ(ierr); }
  ierr = EPSSetFromOptions(ctx->eps);CHKERRQ(ierr);
  if (nep->twosided) {
    if (!ctx->epsts) { ierr = NEPSLPGetEPSLeft(nep,&ctx->epsts);CHKERRQ(ierr); }
    ierr = EPSSetFromOptions(ctx->epsts);CHKERRQ(ierr);
  }
  if (!ctx->ksp) { ierr = NEPSLPGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
  ierr = KSPSetFromOptions(ctx->ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSLPSetDeflationThreshold_SLP(NEP nep,PetscReal deftol)
{
  NEP_SLP *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  ctx->deftol = deftol;
  PetscFunctionReturn(0);
}

/*@
   NEPSLPSetDeflationThreshold - Sets the threshold value used to switch between
   deflated and non-deflated iteration.

   Logically Collective on nep

   Input Parameters:
+  nep    - nonlinear eigenvalue solver
-  deftol - the threshold value

   Options Database Keys:
.  -nep_slp_deflation_threshold <deftol> - set the threshold

   Notes:
   Normally, the solver iterates on the extended problem in order to deflate
   previously converged eigenpairs. If this threshold is set to a nonzero value,
   then once the residual error is below this threshold the solver will
   continue the iteration without deflation. The intention is to be able to
   improve the current eigenpair further, despite having previous eigenpairs
   with somewhat bad precision.

   Level: advanced

.seealso: NEPSLPGetDeflationThreshold()
@*/
PetscErrorCode NEPSLPSetDeflationThreshold(NEP nep,PetscReal deftol)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,deftol,2);
  ierr = PetscTryMethod(nep,"NEPSLPSetDeflationThreshold_C",(NEP,PetscReal),(nep,deftol));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSLPGetDeflationThreshold_SLP(NEP nep,PetscReal *deftol)
{
  NEP_SLP *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  *deftol = ctx->deftol;
  PetscFunctionReturn(0);
}

/*@
   NEPSLPGetDeflationThreshold - Returns the threshold value that controls deflation.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  deftol - the threshold

   Level: advanced

.seealso: NEPSLPSetDeflationThreshold()
@*/
PetscErrorCode NEPSLPGetDeflationThreshold(NEP nep,PetscReal *deftol)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidBoolPointer(deftol,2);
  ierr = PetscUseMethod(nep,"NEPSLPGetDeflationThreshold_C",(NEP,PetscReal*),(nep,deftol));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSLPSetEPS_SLP(NEP nep,EPS eps)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)eps);CHKERRQ(ierr);
  ierr = EPSDestroy(&ctx->eps);CHKERRQ(ierr);
  ctx->eps = eps;
  ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->eps);CHKERRQ(ierr);
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   NEPSLPSetEPS - Associate a linear eigensolver object (EPS) to the
   nonlinear eigenvalue solver.

   Collective on nep

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  eps - the eigensolver object

   Level: advanced

.seealso: NEPSLPGetEPS()
@*/
PetscErrorCode NEPSLPSetEPS(NEP nep,EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(nep,1,eps,2);
  ierr = PetscTryMethod(nep,"NEPSLPSetEPS_C",(NEP,EPS),(nep,eps));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSLPGetEPS_SLP(NEP nep,EPS *eps)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  if (!ctx->eps) {
    ierr = EPSCreate(PetscObjectComm((PetscObject)nep),&ctx->eps);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)ctx->eps,(PetscObject)nep,1);CHKERRQ(ierr);
    ierr = EPSSetOptionsPrefix(ctx->eps,((PetscObject)nep)->prefix);CHKERRQ(ierr);
    ierr = EPSAppendOptionsPrefix(ctx->eps,"nep_slp_");CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->eps);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)ctx->eps,((PetscObject)nep)->options);CHKERRQ(ierr);
  }
  *eps = ctx->eps;
  PetscFunctionReturn(0);
}

/*@
   NEPSLPGetEPS - Retrieve the linear eigensolver object (EPS) associated
   to the nonlinear eigenvalue solver.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  eps - the eigensolver object

   Level: advanced

.seealso: NEPSLPSetEPS()
@*/
PetscErrorCode NEPSLPGetEPS(NEP nep,EPS *eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(eps,2);
  ierr = PetscUseMethod(nep,"NEPSLPGetEPS_C",(NEP,EPS*),(nep,eps));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSLPSetEPSLeft_SLP(NEP nep,EPS eps)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)eps);CHKERRQ(ierr);
  ierr = EPSDestroy(&ctx->epsts);CHKERRQ(ierr);
  ctx->epsts = eps;
  ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->epsts);CHKERRQ(ierr);
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   NEPSLPSetEPSLeft - Associate a linear eigensolver object (EPS) to the
   nonlinear eigenvalue solver, used to compute left eigenvectors in the
   two-sided variant of SLP.

   Collective on nep

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  eps - the eigensolver object

   Level: advanced

.seealso: NEPSLPGetEPSLeft(), NEPSetTwoSided()
@*/
PetscErrorCode NEPSLPSetEPSLeft(NEP nep,EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(nep,1,eps,2);
  ierr = PetscTryMethod(nep,"NEPSLPSetEPSLeft_C",(NEP,EPS),(nep,eps));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSLPGetEPSLeft_SLP(NEP nep,EPS *eps)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  if (!ctx->epsts) {
    ierr = EPSCreate(PetscObjectComm((PetscObject)nep),&ctx->epsts);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)ctx->epsts,(PetscObject)nep,1);CHKERRQ(ierr);
    ierr = EPSSetOptionsPrefix(ctx->epsts,((PetscObject)nep)->prefix);CHKERRQ(ierr);
    ierr = EPSAppendOptionsPrefix(ctx->epsts,"nep_slp_left_");CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->epsts);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)ctx->epsts,((PetscObject)nep)->options);CHKERRQ(ierr);
  }
  *eps = ctx->epsts;
  PetscFunctionReturn(0);
}

/*@
   NEPSLPGetEPSLeft - Retrieve the linear eigensolver object (EPS) associated
   to the nonlinear eigenvalue solver, used to compute left eigenvectors in the
   two-sided variant of SLP.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  eps - the eigensolver object

   Level: advanced

.seealso: NEPSLPSetEPSLeft(), NEPSetTwoSided()
@*/
PetscErrorCode NEPSLPGetEPSLeft(NEP nep,EPS *eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(eps,2);
  ierr = PetscUseMethod(nep,"NEPSLPGetEPSLeft_C",(NEP,EPS*),(nep,eps));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSLPSetKSP_SLP(NEP nep,KSP ksp)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)ksp);CHKERRQ(ierr);
  ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
  ctx->ksp = ksp;
  ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->ksp);CHKERRQ(ierr);
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   NEPSLPSetKSP - Associate a linear solver object (KSP) to the nonlinear
   eigenvalue solver.

   Collective on nep

   Input Parameters:
+  nep - eigenvalue solver
-  ksp - the linear solver object

   Level: advanced

.seealso: NEPSLPGetKSP()
@*/
PetscErrorCode NEPSLPSetKSP(NEP nep,KSP ksp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,2);
  PetscCheckSameComm(nep,1,ksp,2);
  ierr = PetscTryMethod(nep,"NEPSLPSetKSP_C",(NEP,KSP),(nep,ksp));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSLPGetKSP_SLP(NEP nep,KSP *ksp)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  if (!ctx->ksp) {
    ierr = KSPCreate(PetscObjectComm((PetscObject)nep),&ctx->ksp);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)ctx->ksp,(PetscObject)nep,1);CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(ctx->ksp,((PetscObject)nep)->prefix);CHKERRQ(ierr);
    ierr = KSPAppendOptionsPrefix(ctx->ksp,"nep_slp_");CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->ksp);CHKERRQ(ierr);
    ierr = KSPSetErrorIfNotConverged(ctx->ksp,PETSC_TRUE);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ctx->ksp,SLEPC_DEFAULT_TOL,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  }
  *ksp = ctx->ksp;
  PetscFunctionReturn(0);
}

/*@
   NEPSLPGetKSP - Retrieve the linear solver object (KSP) associated with
   the nonlinear eigenvalue solver.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  ksp - the linear solver object

   Level: advanced

.seealso: NEPSLPSetKSP()
@*/
PetscErrorCode NEPSLPGetKSP(NEP nep,KSP *ksp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(ksp,2);
  ierr = PetscUseMethod(nep,"NEPSLPGetKSP_C",(NEP,KSP*),(nep,ksp));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPView_SLP(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    if (ctx->deftol) {
      ierr = PetscViewerASCIIPrintf(viewer,"  deflation threshold: %g\n",(double)ctx->deftol);CHKERRQ(ierr);
    }
    if (!ctx->eps) { ierr = NEPSLPGetEPS(nep,&ctx->eps);CHKERRQ(ierr); }
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = EPSView(ctx->eps,viewer);CHKERRQ(ierr);
    if (nep->twosided) {
      if (!ctx->epsts) { ierr = NEPSLPGetEPSLeft(nep,&ctx->epsts);CHKERRQ(ierr); }
      ierr = EPSView(ctx->epsts,viewer);CHKERRQ(ierr);
    }
    if (!ctx->ksp) { ierr = NEPSLPGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
    ierr = KSPView(ctx->ksp,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPReset_SLP(NEP nep)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  ierr = EPSReset(ctx->eps);CHKERRQ(ierr);
  if (nep->twosided) { ierr = EPSReset(ctx->epsts);CHKERRQ(ierr); }
  ierr = KSPReset(ctx->ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDestroy_SLP(NEP nep)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;

  PetscFunctionBegin;
  ierr = KSPDestroy(&ctx->ksp);CHKERRQ(ierr);
  ierr = EPSDestroy(&ctx->eps);CHKERRQ(ierr);
  ierr = EPSDestroy(&ctx->epsts);CHKERRQ(ierr);
  ierr = PetscFree(nep->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetDeflationThreshold_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetDeflationThreshold_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetEPS_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetEPS_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetEPSLeft_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetEPSLeft_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetKSP_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetKSP_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_SLP(NEP nep)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(nep,&ctx);CHKERRQ(ierr);
  nep->data = (void*)ctx;

  nep->useds = PETSC_TRUE;
  nep->hasts = PETSC_TRUE;

  ctx->deftol = 0.0;

  nep->ops->solve          = NEPSolve_SLP;
  nep->ops->setup          = NEPSetUp_SLP;
  nep->ops->setfromoptions = NEPSetFromOptions_SLP;
  nep->ops->reset          = NEPReset_SLP;
  nep->ops->destroy        = NEPDestroy_SLP;
  nep->ops->view           = NEPView_SLP;
  nep->ops->computevectors = NEPComputeVectors_Schur;

  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetDeflationThreshold_C",NEPSLPSetDeflationThreshold_SLP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetDeflationThreshold_C",NEPSLPGetDeflationThreshold_SLP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetEPS_C",NEPSLPSetEPS_SLP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetEPS_C",NEPSLPGetEPS_SLP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetEPSLeft_C",NEPSLPSetEPSLeft_SLP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetEPSLeft_C",NEPSLPGetEPSLeft_SLP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPSetKSP_C",NEPSLPSetKSP_SLP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPSLPGetKSP_C",NEPSLPGetKSP_SLP);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

