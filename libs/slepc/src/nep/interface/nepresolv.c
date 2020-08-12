/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to resolvent T^{-1}(z) = sum_i (z-lambda_i)^{-1} x_i y_i'
*/

#include <slepc/private/nepimpl.h>       /*I "slepcnep.h" I*/

typedef struct {
  NEP              nep;
  RG               rg;
  PetscScalar      omega;
  PetscScalar      *nfactor;         /* normalization factors y_i'*T'(lambda_i)*x_i */
  PetscBool        *nfactor_avail;
  PetscScalar      *dots;            /* inner products y_i'*v */
  PetscBool        *dots_avail;
  PetscObjectId    vid;
  PetscObjectState vstate;
} NEP_RESOLVENT_MATSHELL;

static PetscErrorCode MatMult_Resolvent(Mat M,Vec v,Vec r)
{
  PetscErrorCode         ierr;
  NEP_RESOLVENT_MATSHELL *ctx;
  NEP                    nep;
  PetscInt               i,inside=1;
  PetscScalar            alpha;
  Vec                    x,y,z,w;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&ctx);CHKERRQ(ierr);
  nep = ctx->nep;
  w = nep->work[0];
  z = nep->work[1];
  if (((PetscObject)v)->id != ctx->vid || ((PetscObject)v)->state != ctx->vstate) {
    ierr = PetscArrayzero(ctx->dots_avail,ctx->nep->nconv);CHKERRQ(ierr);
    ierr = PetscObjectGetId((PetscObject)v,&ctx->vid);CHKERRQ(ierr);
    ierr = PetscObjectStateGet((PetscObject)v,&ctx->vstate);CHKERRQ(ierr);
  }
  ierr = VecSet(r,0.0);CHKERRQ(ierr);
  for (i=0;i<nep->nconv;i++) {
    if (ctx->rg) {
      ierr = RGCheckInside(ctx->rg,1,&nep->eigr[i],&nep->eigi[i],&inside);CHKERRQ(ierr);
    }
    if (inside>=0) {
      ierr = BVGetColumn(nep->V,i,&x);CHKERRQ(ierr);
      ierr = BVGetColumn(nep->W,i,&y);CHKERRQ(ierr);
      ierr = NEPApplyJacobian(nep,nep->eigr[i],x,z,w,NULL);CHKERRQ(ierr);
      if (!ctx->dots_avail[i]) {
        ierr = VecDot(v,y,&ctx->dots[i]);CHKERRQ(ierr);
        ctx->dots_avail[i] = PETSC_TRUE;
      }
      if (!ctx->nfactor_avail[i]) {
        ierr = VecDot(w,y,&ctx->nfactor[i]);CHKERRQ(ierr);
        ctx->nfactor_avail[i] = PETSC_TRUE;
      }
      alpha = ctx->dots[i]/(ctx->nfactor[i]*(ctx->omega-nep->eigr[i]));
      ierr = VecAXPY(r,alpha,x);CHKERRQ(ierr);
      ierr = BVRestoreColumn(nep->V,i,&x);CHKERRQ(ierr);
      ierr = BVRestoreColumn(nep->W,i,&y);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDestroy_Resolvent(Mat M)
{
  PetscErrorCode         ierr;
  NEP_RESOLVENT_MATSHELL *ctx;

  PetscFunctionBegin;
  if (M) {
    ierr = MatShellGetContext(M,(void**)&ctx);CHKERRQ(ierr);
    ierr = PetscFree4(ctx->nfactor,ctx->nfactor_avail,ctx->dots,ctx->dots_avail);CHKERRQ(ierr);
    ierr = PetscFree(ctx);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   NEPApplyResolvent - Applies the resolvent T^{-1}(z) to a given vector.

   Collective on nep

   Input Parameters:
+  nep   - eigensolver context obtained from NEPCreate()
.  rg    - optional region
.  omega - value where the resolvent must be evaluated
-  v     - input vector

   Output Parameter:
.  r     - result vector

   Notes:
   The resolvent T^{-1}(z) = sum_i (z-lambda_i)^{-1}*x_i*y_i' is evaluated at
   z=omega and the matrix-vector multiplication r = T^{-1}(omega)*v is computed.
   Vectors x_i and y_i are right and left eigenvectors, respectively, normalized
   so that y_i'*T'(lambda_i)*x_i=1. The sum contains only eigenvectors that have
   been previously computed with NEPSolve(), and if a region rg is given then only
   those corresponding to eigenvalues inside the region are considered.

   Level: intermediate

.seealso: NEPGetLeftEigenvector(), NEPSolve()
@*/
PetscErrorCode NEPApplyResolvent(NEP nep,RG rg,PetscScalar omega,Vec v,Vec r)
{
  PetscErrorCode         ierr;
  NEP_RESOLVENT_MATSHELL *ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(nep,omega,3);
  PetscValidHeaderSpecific(v,VEC_CLASSID,4);
  PetscValidHeaderSpecific(r,VEC_CLASSID,5);
  NEPCheckSolved(nep,1);

  ierr = PetscLogEventBegin(NEP_Resolvent,nep,0,0,0);CHKERRQ(ierr);
  if (!nep->resolvent) {
    ierr = PetscNew(&ctx);CHKERRQ(ierr);
    ctx->nep = nep;
    ierr = PetscCalloc4(nep->nconv,&ctx->nfactor,nep->nconv,&ctx->nfactor_avail,nep->nconv,&ctx->dots,nep->nconv,&ctx->dots_avail);CHKERRQ(ierr);
    ierr = MatCreateShell(PetscObjectComm((PetscObject)nep),nep->nloc,nep->nloc,nep->n,nep->n,ctx,&nep->resolvent);CHKERRQ(ierr);
    ierr = MatShellSetOperation(nep->resolvent,MATOP_MULT,(void(*)(void))MatMult_Resolvent);CHKERRQ(ierr);
    ierr = MatShellSetOperation(nep->resolvent,MATOP_DESTROY,(void(*)(void))MatDestroy_Resolvent);CHKERRQ(ierr);
  } else {
    ierr = MatShellGetContext(nep->resolvent,(void**)&ctx);CHKERRQ(ierr);
  }
  ierr = NEPComputeVectors(nep);CHKERRQ(ierr);
  ierr = NEPSetWorkVecs(nep,2);CHKERRQ(ierr);
  ctx->rg    = rg;
  ctx->omega = omega;
  ierr = MatMult(nep->resolvent,v,r);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(NEP_Resolvent,nep,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

