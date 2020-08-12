/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Two-sided variant of the NEPSLP solver.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include "slp.h"

typedef struct _n_nep_def_ctx *NEP_NEDEF_CTX;

struct _n_nep_def_ctx {
  PetscInt    n;
  PetscBool   ref;
  PetscScalar *eig;
  BV          V,W;
};

typedef struct {   /* context for two-sided solver */
  Mat         Ft;
  Mat         Jt;
  Vec         w;
} NEP_SLPTS_MATSHELL;

typedef struct {   /* context for non-equivalence deflation */
  NEP_NEDEF_CTX defctx;
  Mat           F;
  Mat           J;
  KSP           ksp;
  PetscBool     isJ;
  PetscScalar   lambda;
  Vec           w[2];
} NEP_NEDEF_MATSHELL;

static PetscErrorCode MatMult_SLPTS_Right(Mat M,Vec x,Vec y)
{
  PetscErrorCode     ierr;
  NEP_SLPTS_MATSHELL *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&ctx);CHKERRQ(ierr);
  ierr = MatMult(ctx->Jt,x,ctx->w);CHKERRQ(ierr);
  ierr = MatSolve(ctx->Ft,ctx->w,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_SLPTS_Left(Mat M,Vec x,Vec y)
{
  PetscErrorCode     ierr;
  NEP_SLPTS_MATSHELL *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&ctx);CHKERRQ(ierr);
  ierr = MatMultTranspose(ctx->Jt,x,ctx->w);CHKERRQ(ierr);
  ierr = MatSolveTranspose(ctx->Ft,ctx->w,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDestroy_SLPTS(Mat M)
{
  PetscErrorCode     ierr;
  NEP_SLPTS_MATSHELL *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->w);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_CUDA)
static PetscErrorCode MatCreateVecs_SLPTS(Mat M,Vec *left,Vec *right)
{
  PetscErrorCode     ierr;
  NEP_SLPTS_MATSHELL *ctx;

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

static PetscErrorCode NEPSLPSetUpEPSMat(NEP nep,Mat F,Mat J,PetscBool left,Mat *M)
{
  PetscErrorCode     ierr;
  Mat                Mshell;
  PetscInt           nloc,mloc;
  NEP_SLPTS_MATSHELL *shellctx;

  PetscFunctionBegin;
  /* Create mat shell */
  ierr = PetscNew(&shellctx);CHKERRQ(ierr);
  shellctx->Ft = F;
  shellctx->Jt = J;
  ierr = MatGetLocalSize(nep->function,&mloc,&nloc);CHKERRQ(ierr);
  ierr = MatCreateShell(PetscObjectComm((PetscObject)nep),nloc,mloc,PETSC_DETERMINE,PETSC_DETERMINE,shellctx,&Mshell);CHKERRQ(ierr);
  if (left) {
    ierr = MatShellSetOperation(Mshell,MATOP_MULT,(void(*)(void))MatMult_SLPTS_Left);CHKERRQ(ierr);
  } else {
    ierr = MatShellSetOperation(Mshell,MATOP_MULT,(void(*)(void))MatMult_SLPTS_Right);CHKERRQ(ierr);
  }
  ierr = MatShellSetOperation(Mshell,MATOP_DESTROY,(void(*)(void))MatDestroy_SLPTS);CHKERRQ(ierr);
#if defined(PETSC_HAVE_CUDA)
  ierr = MatShellSetOperation(Mshell,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_SLPTS);CHKERRQ(ierr);
#endif
  *M = Mshell;
  ierr = MatCreateVecs(nep->function,&shellctx->w,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Functions for deflation */
static PetscErrorCode NEPDeflationNEDestroy(NEP_NEDEF_CTX defctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!defctx) PetscFunctionReturn(0);
  ierr = PetscFree(defctx->eig);CHKERRQ(ierr);
  ierr = PetscFree(defctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationNECreate(NEP nep,BV V,BV W,PetscInt sz,NEP_NEDEF_CTX *defctx)
{
  PetscErrorCode ierr;
  NEP_NEDEF_CTX  op;

  PetscFunctionBegin;
  ierr = PetscNew(&op);CHKERRQ(ierr);
  *defctx = op;
  op->n   = 0;
  op->ref = PETSC_FALSE;
  ierr = PetscCalloc1(sz,&op->eig);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)W);CHKERRQ(ierr);
  op->V = V;
  op->W = W;
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationNEComputeFunction(NEP nep,Mat M,PetscScalar lambda)
{
  PetscErrorCode     ierr;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  if (lambda==matctx->lambda) PetscFunctionReturn(0);
  ierr = NEPComputeFunction(nep,lambda,matctx->F,matctx->F);CHKERRQ(ierr);
  if (matctx->isJ) {ierr = NEPComputeJacobian(nep,lambda,matctx->J);CHKERRQ(ierr);}
  if (matctx->ksp) {ierr = KSPSetOperators(matctx->ksp,matctx->F,matctx->F);CHKERRQ(ierr);}
  matctx->lambda = lambda;
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_NEPDeflationNE(Mat M,Vec x,Vec r)
{
  Vec                t,tt;
  PetscScalar        *h,*alpha,lambda,*eig;
  PetscInt           i,k;
  PetscErrorCode     ierr;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  if (matctx->defctx->n && !matctx->defctx->ref) {
    k = matctx->defctx->n;
    lambda = matctx->lambda;
    eig = matctx->defctx->eig;
    t = matctx->w[0];
    ierr = VecCopy(x,t);CHKERRQ(ierr);
    ierr = PetscMalloc2(k,&h,k,&alpha);CHKERRQ(ierr);
    for (i=0;i<k;i++) alpha[i] = (lambda-eig[i]-1.0)/(lambda-eig[i]);
    ierr = BVDotVec(matctx->defctx->V,t,h);CHKERRQ(ierr);
    for (i=0;i<k;i++) h[i] *= alpha[i];
    ierr = BVMultVec(matctx->defctx->W,-1.0,1.0,t,h);CHKERRQ(ierr);
    ierr = MatMult(matctx->isJ?matctx->J:matctx->F,t,r);CHKERRQ(ierr);
    if (matctx->isJ) {
      for (i=0;i<k;i++) h[i] *= (1.0/((lambda-eig[i])*(lambda-eig[i])))/alpha[i];
      tt = matctx->w[1];
      ierr = BVMultVec(matctx->defctx->W,-1.0,0.0,tt,h);CHKERRQ(ierr);
      ierr = MatMult(matctx->F,tt,t);CHKERRQ(ierr);
      ierr = VecAXPY(r,1.0,t);CHKERRQ(ierr);
    }
    ierr = PetscFree2(h,alpha);CHKERRQ(ierr);
  } else {
    ierr = MatMult(matctx->isJ?matctx->J:matctx->F,x,r);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMultTranspose_NEPDeflationNE(Mat M,Vec x,Vec r)
{
  Vec                t,tt;
  PetscScalar        *h,*alphaC,lambda,*eig;
  PetscInt           i,k;
  PetscErrorCode     ierr;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  t    = matctx->w[0];
  ierr = VecCopy(x,t);CHKERRQ(ierr);
  if (matctx->defctx->n && !matctx->defctx->ref) {
    ierr = VecConjugate(t);CHKERRQ(ierr);
    k = matctx->defctx->n;
    lambda = matctx->lambda;
    eig = matctx->defctx->eig;
    ierr = PetscMalloc2(k,&h,k,&alphaC);CHKERRQ(ierr);
    for (i=0;i<k;i++) alphaC[i] = PetscConj((lambda-eig[i]-1.0)/(lambda-eig[i]));
    ierr = BVDotVec(matctx->defctx->W,t,h);CHKERRQ(ierr);
    for (i=0;i<k;i++) h[i] *= alphaC[i];
    ierr = BVMultVec(matctx->defctx->V,-1.0,1.0,t,h);CHKERRQ(ierr);
    ierr = VecConjugate(t);CHKERRQ(ierr);
    ierr = MatMultTranspose(matctx->isJ?matctx->J:matctx->F,t,r);CHKERRQ(ierr);
    if (matctx->isJ) {
      for (i=0;i<k;i++) h[i] *= PetscConj(1.0/((lambda-eig[i])*(lambda-eig[i])))/alphaC[i];
      tt = matctx->w[1];
      ierr = BVMultVec(matctx->defctx->V,-1.0,0.0,tt,h);CHKERRQ(ierr);
      ierr = VecConjugate(tt);CHKERRQ(ierr);
      ierr = MatMultTranspose(matctx->F,tt,t);CHKERRQ(ierr);
      ierr = VecAXPY(r,1.0,t);CHKERRQ(ierr);
    }
    ierr = PetscFree2(h,alphaC);CHKERRQ(ierr);
  } else {
    ierr = MatMultTranspose(matctx->isJ?matctx->J:matctx->F,t,r);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatSolve_NEPDeflationNE(Mat M,Vec b,Vec x)
{
  PetscErrorCode     ierr;
  PetscScalar        *h,*alpha,lambda,*eig;
  PetscInt           i,k;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  if (!matctx->ksp) {
    ierr = VecCopy(b,x);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = KSPSolve(matctx->ksp,b,x);CHKERRQ(ierr);
  if (matctx->defctx->n && !matctx->defctx->ref) {
    k = matctx->defctx->n;
    lambda = matctx->lambda;
    eig = matctx->defctx->eig;
    ierr = PetscMalloc2(k,&h,k,&alpha);CHKERRQ(ierr);
    ierr = BVDotVec(matctx->defctx->V,x,h);CHKERRQ(ierr);
    for (i=0;i<k;i++) alpha[i] = (lambda-eig[i]-1.0)/(lambda-eig[i]);
    for (i=0;i<k;i++) h[i] *= alpha[i]/(1.0-alpha[i]);
    ierr = BVMultVec(matctx->defctx->W,1.0,1.0,x,h);CHKERRQ(ierr);
    ierr = PetscFree2(h,alpha);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatSolveTranspose_NEPDeflationNE(Mat M,Vec b,Vec x)
{
  PetscErrorCode     ierr;
  PetscScalar        *h,*alphaC,lambda,*eig;
  PetscInt           i,k;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  if (!matctx->ksp) {
    ierr = VecCopy(b,x);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = KSPSolveTranspose(matctx->ksp,b,x);CHKERRQ(ierr);
  if (matctx->defctx->n && !matctx->defctx->ref) {
    ierr = VecConjugate(x);CHKERRQ(ierr);
    k = matctx->defctx->n;
    lambda = matctx->lambda;
    eig = matctx->defctx->eig;
    ierr = PetscMalloc2(k,&h,k,&alphaC);CHKERRQ(ierr);
    ierr = BVDotVec(matctx->defctx->W,x,h);CHKERRQ(ierr);
    for (i=0;i<k;i++) alphaC[i] = PetscConj((lambda-eig[i]-1.0)/(lambda-eig[i]));
    for (i=0;i<k;i++) h[i] *= alphaC[i]/(1.0-alphaC[i]);
    ierr = BVMultVec(matctx->defctx->V,1.0,1.0,x,h);CHKERRQ(ierr);
    ierr = PetscFree2(h,alphaC);CHKERRQ(ierr);
    ierr = VecConjugate(x);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDestroy_NEPDeflationNE(Mat M)
{
  PetscErrorCode     ierr;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  ierr = VecDestroy(&matctx->w[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&matctx->w[1]);CHKERRQ(ierr);
  ierr = PetscFree(matctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatCreateVecs_NEPDeflationNE(Mat M,Vec *right,Vec *left)
{
  PetscErrorCode     ierr;
  NEP_NEDEF_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  ierr = MatCreateVecs(matctx->F,right,left);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationNEFunctionCreate(NEP_NEDEF_CTX defctx,NEP nep,Mat F,Mat J,KSP ksp,PetscBool isJ,Mat *Mshell)
{
  NEP_NEDEF_MATSHELL *matctx;
  PetscErrorCode     ierr;
  PetscInt           nloc,mloc;

  PetscFunctionBegin;
  /* Create mat shell */
  ierr = PetscNew(&matctx);CHKERRQ(ierr);
  ierr = MatGetLocalSize(nep->function,&mloc,&nloc);CHKERRQ(ierr);
  ierr = MatCreateShell(PetscObjectComm((PetscObject)nep),nloc,mloc,PETSC_DETERMINE,PETSC_DETERMINE,matctx,Mshell);CHKERRQ(ierr);
  matctx->F   = F;
  matctx->J   = J;
  matctx->isJ = isJ;
  matctx->ksp = ksp;
  matctx->defctx = defctx;
  matctx->lambda = PETSC_MAX_REAL;
  ierr = MatCreateVecs(F,&matctx->w[0],NULL);CHKERRQ(ierr);
  ierr = VecDuplicate(matctx->w[0],&matctx->w[1]);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Mshell,MATOP_MULT,(void(*)(void))MatMult_NEPDeflationNE);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Mshell,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMultTranspose_NEPDeflationNE);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Mshell,MATOP_SOLVE,(void(*)(void))MatSolve_NEPDeflationNE);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Mshell,MATOP_SOLVE_TRANSPOSE,(void(*)(void))MatSolveTranspose_NEPDeflationNE);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Mshell,MATOP_DESTROY,(void(*)(void))MatDestroy_NEPDeflationNE);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Mshell,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_NEPDeflationNE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationNERecoverEigenvectors(NEP_NEDEF_CTX defctx,Vec u,Vec w,PetscScalar lambda)
{
  PetscErrorCode ierr;
  PetscScalar    *h,*alpha,*eig;
  PetscInt       i,k;

  PetscFunctionBegin;
  if (w) { ierr = VecConjugate(w);CHKERRQ(ierr); }
  if (defctx->n && !defctx->ref) {
    eig = defctx->eig;
    k = defctx->n;
    ierr = PetscMalloc2(k,&h,k,&alpha);CHKERRQ(ierr);
    for (i=0;i<k;i++) alpha[i] = (lambda-eig[i]-1.0)/(lambda-eig[i]);
    ierr = BVDotVec(defctx->V,u,h);CHKERRQ(ierr);
    for (i=0;i<k;i++) h[i] *= alpha[i];
    ierr = BVMultVec(defctx->W,-1.0,1.0,u,h);CHKERRQ(ierr);
    ierr = VecNormalize(u,NULL);CHKERRQ(ierr);
    if (w) {
      ierr = BVDotVec(defctx->W,w,h);CHKERRQ(ierr);
      for (i=0;i<k;i++) alpha[i] = PetscConj((lambda-eig[i]-1.0)/(lambda-eig[i]));
      for (i=0;i<k;i++) h[i] *= alpha[i];
      ierr = BVMultVec(defctx->V,-1.0,1.0,w,h);CHKERRQ(ierr);
      ierr = VecNormalize(w,NULL);CHKERRQ(ierr);
    }
    ierr = PetscFree2(h,alpha);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationNELocking(NEP_NEDEF_CTX defctx,Vec u,Vec w,PetscScalar lambda)
{
  PetscErrorCode ierr;
  PetscInt       n;

  PetscFunctionBegin;
  n = defctx->n++;
  defctx->eig[n] = lambda;
  ierr = BVInsertVec(defctx->V,n,u);CHKERRQ(ierr);
  ierr = BVInsertVec(defctx->W,n,w);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(defctx->V,0,defctx->n);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(defctx->W,0,defctx->n);CHKERRQ(ierr);
  ierr = BVBiorthonormalizeColumn(defctx->V,defctx->W,n,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationNESetRefine(NEP_NEDEF_CTX defctx,PetscBool ref)
{
  PetscFunctionBegin;
  defctx->ref = ref;
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSolve_SLP_Twosided(NEP nep)
{
  PetscErrorCode ierr;
  NEP_SLP        *ctx = (NEP_SLP*)nep->data;
  Mat            mF,mJ,M,Mt;
  Vec            u,r,t,w;
  BV             X,Y;
  PetscScalar    sigma,lambda,mu,im=0.0,mu2,im2;
  PetscReal      resnorm,resl;
  PetscInt       nconv,nconv2,i;
  PetscBool      skip=PETSC_FALSE,lock=PETSC_FALSE;
  NEP_NEDEF_CTX  defctx=NULL;    /* Extended operator for deflation */

  PetscFunctionBegin;
  /* get initial approximation of eigenvalue and eigenvector */
  ierr = NEPGetDefaultShift(nep,&sigma);CHKERRQ(ierr);
  if (!nep->nini) {
    ierr = BVSetRandomColumn(nep->V,0);CHKERRQ(ierr);
  }
  ierr = BVSetRandomColumn(nep->W,0);CHKERRQ(ierr);
  lambda = sigma;
  if (!ctx->ksp) { ierr = NEPSLPGetKSP(nep,&ctx->ksp);CHKERRQ(ierr); }
  ierr = BVDuplicate(nep->V,&X);CHKERRQ(ierr);
  ierr = BVDuplicate(nep->W,&Y);CHKERRQ(ierr);
  ierr = NEPDeflationNECreate(nep,X,Y,nep->nev,&defctx);CHKERRQ(ierr);
  ierr = BVGetColumn(nep->V,0,&t);CHKERRQ(ierr);
  ierr = VecDuplicate(t,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(t,&w);CHKERRQ(ierr);
  ierr = BVRestoreColumn(nep->V,0,&t);CHKERRQ(ierr);
  ierr = BVCopyVec(nep->V,0,u);CHKERRQ(ierr);
  ierr = BVCopyVec(nep->W,0,w);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&r);CHKERRQ(ierr);
  ierr = NEPDeflationNEFunctionCreate(defctx,nep,nep->function,NULL,ctx->ksp,PETSC_FALSE,&mF);CHKERRQ(ierr);
  ierr = NEPDeflationNEFunctionCreate(defctx,nep,nep->function,nep->jacobian,NULL,PETSC_TRUE,&mJ);CHKERRQ(ierr);
  ierr = NEPSLPSetUpEPSMat(nep,mF,mJ,PETSC_FALSE,&M);CHKERRQ(ierr);
  ierr = NEPSLPSetUpEPSMat(nep,mF,mJ,PETSC_TRUE,&Mt);CHKERRQ(ierr);
  ierr = EPSSetOperators(ctx->eps,M,NULL);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = EPSSetOperators(ctx->epsts,Mt,NULL);CHKERRQ(ierr);
  ierr = MatDestroy(&Mt);CHKERRQ(ierr);

  /* Restart loop */
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    nep->its++;

    /* form residual,  r = T(lambda)*u (used in convergence test only) */
    ierr = NEPDeflationNEComputeFunction(nep,mF,lambda);CHKERRQ(ierr);
    ierr = MatMultTranspose(mF,w,r);CHKERRQ(ierr);
    ierr = VecNorm(r,NORM_2,&resl);CHKERRQ(ierr);
    ierr = MatMult(mF,u,r);CHKERRQ(ierr);

    /* convergence test */
    ierr = VecNorm(r,NORM_2,&resnorm);CHKERRQ(ierr);
    resnorm = PetscMax(resnorm,resl);
    ierr = (*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx);CHKERRQ(ierr);
    nep->eigr[nep->nconv] = lambda;
    if (nep->errest[nep->nconv]<=nep->tol || nep->errest[nep->nconv]<=ctx->deftol) {
      if (nep->errest[nep->nconv]<=ctx->deftol && !defctx->ref && nep->nconv) {
        ierr = NEPDeflationNERecoverEigenvectors(defctx,u,w,lambda);CHKERRQ(ierr);
        ierr = VecConjugate(w);CHKERRQ(ierr);
        ierr = NEPDeflationNESetRefine(defctx,PETSC_TRUE);CHKERRQ(ierr);
        ierr = MatMultTranspose(mF,w,r);CHKERRQ(ierr);
        ierr = VecNorm(r,NORM_2,&resl);CHKERRQ(ierr);
        ierr = MatMult(mF,u,r);CHKERRQ(ierr);
        ierr = VecNorm(r,NORM_2,&resnorm);CHKERRQ(ierr);
        resnorm = PetscMax(resnorm,resl);
        ierr = (*nep->converged)(nep,lambda,0,resnorm,&nep->errest[nep->nconv],nep->convergedctx);CHKERRQ(ierr);
        if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
      } else if (nep->errest[nep->nconv]<=nep->tol) lock = PETSC_TRUE;
    }
    if (lock) {
      lock = PETSC_FALSE;
      skip = PETSC_TRUE;
      ierr = NEPDeflationNERecoverEigenvectors(defctx,u,w,lambda);CHKERRQ(ierr);
      ierr = NEPDeflationNELocking(defctx,u,w,lambda);CHKERRQ(ierr);
      ierr = NEPDeflationNESetRefine(defctx,PETSC_FALSE);CHKERRQ(ierr);
      ierr = BVInsertVec(nep->V,nep->nconv,u);CHKERRQ(ierr);
      ierr = BVInsertVec(nep->W,nep->nconv,w);CHKERRQ(ierr);
      ierr = VecConjugate(w);CHKERRQ(ierr);
      nep->nconv = nep->nconv + 1;
    }
    ierr = (*nep->stopping)(nep,nep->its,nep->max_it,nep->nconv,nep->nev,&nep->reason,nep->stoppingctx);CHKERRQ(ierr);
    if (!skip || nep->reason>0) {
      ierr = NEPMonitor(nep,nep->its,nep->nconv,nep->eigr,nep->eigi,nep->errest,(nep->reason>0)?nep->nconv:nep->nconv+1);CHKERRQ(ierr);
    }

    if (nep->reason == NEP_CONVERGED_ITERATING) {
      if (!skip) {
        /* evaluate T(lambda) and T'(lambda) */
        ierr = NEPDeflationNEComputeFunction(nep,mF,lambda);CHKERRQ(ierr);
        ierr = NEPDeflationNEComputeFunction(nep,mJ,lambda);CHKERRQ(ierr);
        ierr = EPSSetInitialSpace(ctx->eps,1,&u);CHKERRQ(ierr);
        ierr = EPSSetInitialSpace(ctx->epsts,1,&w);CHKERRQ(ierr);

        /* compute new eigenvalue correction mu and eigenvector approximation u */
        ierr = EPSSolve(ctx->eps);CHKERRQ(ierr);
        ierr = EPSSolve(ctx->epsts);CHKERRQ(ierr);
        ierr = EPSGetConverged(ctx->eps,&nconv);CHKERRQ(ierr);
        ierr = EPSGetConverged(ctx->epsts,&nconv2);CHKERRQ(ierr);
        if (!nconv||!nconv2) {
          ierr = PetscInfo1(nep,"iter=%D, inner iteration failed, stopping solve\n",nep->its);CHKERRQ(ierr);
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }
        ierr = EPSGetEigenpair(ctx->eps,0,&mu,&im,u,NULL);CHKERRQ(ierr);
        for (i=0;i<nconv2;i++) {
          ierr = EPSGetEigenpair(ctx->epsts,i,&mu2,&im2,w,NULL);CHKERRQ(ierr);
          if (SlepcAbsEigenvalue(mu-mu2,im-im2)/SlepcAbsEigenvalue(mu,im)<nep->tol*1000) break;
        }
        if (i==nconv2) {
          ierr = PetscInfo1(nep,"iter=%D, inner iteration failed, stopping solve\n",nep->its);CHKERRQ(ierr);
          nep->reason = NEP_DIVERGED_LINEAR_SOLVE;
          break;
        }

        mu = 1.0/mu;
        if (PetscAbsScalar(im)>PETSC_MACHINE_EPSILON) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Complex eigenvalue approximation - not implemented in real scalars");
      } else {
        nep->its--;  /* do not count this as a full iteration */
        /* use second eigenpair computed in previous iteration */
        ierr = EPSGetConverged(ctx->eps,&nconv);CHKERRQ(ierr);
        if (nconv>=2 && nconv2>=2) {
          ierr = EPSGetEigenpair(ctx->eps,1,&mu,&im,u,NULL);CHKERRQ(ierr);
          ierr = EPSGetEigenpair(ctx->epsts,1,&mu2,&im2,w,NULL);CHKERRQ(ierr);
          mu = 1.0/mu;
        } else {
          ierr = BVSetRandomColumn(nep->V,nep->nconv);CHKERRQ(ierr);
          ierr = BVSetRandomColumn(nep->W,nep->nconv);CHKERRQ(ierr);
          ierr = BVCopyVec(nep->V,nep->nconv,u);CHKERRQ(ierr);
          ierr = BVCopyVec(nep->W,nep->nconv,w);CHKERRQ(ierr);
          mu = lambda-sigma;
        }
        skip = PETSC_FALSE;
      }
      /* correct eigenvalue */
      lambda = lambda - mu;
    }
  }
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = MatDestroy(&mF);CHKERRQ(ierr);
  ierr = MatDestroy(&mJ);CHKERRQ(ierr);
  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = BVDestroy(&Y);CHKERRQ(ierr);
  ierr = NEPDeflationNEDestroy(defctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

