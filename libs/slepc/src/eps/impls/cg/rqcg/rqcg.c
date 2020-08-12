/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "rqcg"

   Method: Rayleigh Quotient Conjugate Gradient

   Algorithm:

       Conjugate Gradient minimization of the Rayleigh quotient with
       periodic Rayleigh-Ritz acceleration.

   References:

       [1] L. Bergamaschi et al., "Parallel preconditioned conjugate gradient
           optimization of the Rayleigh quotient for the solution of sparse
           eigenproblems", Appl. Math. Comput. 175(2):1694-1715, 2006.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/

PetscErrorCode EPSSolve_RQCG(EPS);

typedef struct {
  PetscInt nrest;         /* user-provided reset parameter */
  PetscInt allocsize;     /* number of columns of work BV's allocated at setup */
  BV       AV,W,P,G;
} EPS_RQCG;

PetscErrorCode EPSSetUp_RQCG(EPS eps)
{
  PetscErrorCode ierr;
  PetscInt       nmat;
  EPS_RQCG       *ctx = (EPS_RQCG*)eps->data;

  PetscFunctionBegin;
  if (!eps->ishermitian || (eps->isgeneralized && !eps->ispositive)) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"RQCG only works for Hermitian problems");
  ierr = EPSSetDimensions_Default(eps,eps->nev,&eps->ncv,&eps->mpd);CHKERRQ(ierr);
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) eps->which = EPS_SMALLEST_REAL;
  if (eps->which!=EPS_SMALLEST_REAL) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Wrong value of eps->which");
  if (!eps->extraction) {
    ierr = EPSSetExtraction(eps,EPS_RITZ);CHKERRQ(ierr);
  } else if (eps->extraction!=EPS_RITZ) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");

  if (!ctx->nrest) ctx->nrest = 20;

  ierr = EPSAllocateSolution(eps,0);CHKERRQ(ierr);
  ierr = EPS_SetInnerProduct(eps);CHKERRQ(ierr);

  ierr = STGetNumMatrices(eps->st,&nmat);CHKERRQ(ierr);
  if (!ctx->allocsize) {
    ctx->allocsize = eps->mpd;
    ierr = BVDuplicateResize(eps->V,eps->mpd,&ctx->AV);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)ctx->AV);CHKERRQ(ierr);
    if (nmat>1) {
      ierr = BVDuplicate(ctx->AV,&ctx->W);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)ctx->W);CHKERRQ(ierr);
    }
    ierr = BVDuplicate(ctx->AV,&ctx->P);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)ctx->P);CHKERRQ(ierr);
    ierr = BVDuplicate(ctx->AV,&ctx->G);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)ctx->G);CHKERRQ(ierr);
  } else if (ctx->allocsize!=eps->mpd) {
    ctx->allocsize = eps->mpd;
    ierr = BVResize(ctx->AV,eps->mpd,PETSC_FALSE);CHKERRQ(ierr);
    if (nmat>1) {
      ierr = BVResize(ctx->W,eps->mpd,PETSC_FALSE);CHKERRQ(ierr);
    }
    ierr = BVResize(ctx->P,eps->mpd,PETSC_FALSE);CHKERRQ(ierr);
    ierr = BVResize(ctx->G,eps->mpd,PETSC_FALSE);CHKERRQ(ierr);
  }
  ierr = DSSetType(eps->ds,DSHEP);CHKERRQ(ierr);
  ierr = DSAllocate(eps->ds,eps->ncv);CHKERRQ(ierr);
  ierr = EPSSetWorkVecs(eps,1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   ExtractSubmatrix - Returns B = A(k+1:end,k+1:end).
*/
static PetscErrorCode ExtractSubmatrix(Mat A,PetscInt k,Mat *B)
{
  PetscErrorCode ierr;
  PetscInt       j,m,n;
  PetscScalar    *pA,*pB;

  PetscFunctionBegin;
  ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,m-k,n-k,NULL,B);CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&pA);CHKERRQ(ierr);
  ierr = MatDenseGetArray(*B,&pB);CHKERRQ(ierr);
  for (j=k;j<n;j++) {
    ierr = PetscArraycpy(pB+(j-k)*(m-k),pA+j*m+k,m-k);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(A,&pA);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(*B,&pB);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_RQCG(EPS eps)
{
  PetscErrorCode ierr;
  EPS_RQCG       *ctx = (EPS_RQCG*)eps->data;
  PetscInt       i,j,k,ld,nv,ncv = eps->ncv,kini,nmat;
  PetscScalar    *C,*gamma,g,pap,pbp,pbx,pax,nu,mu,alpha,beta;
  PetscReal      resnorm,a,b,c,d,disc,t;
  PetscBool      reset;
  Mat            A,B,Q,Q1;
  Vec            v,av,bv,p,w=eps->work[0];

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  ierr = STGetNumMatrices(eps->st,&nmat);CHKERRQ(ierr);
  ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
  if (nmat>1) { ierr = STGetMatrix(eps->st,1,&B);CHKERRQ(ierr); }
  else B = NULL;
  ierr = PetscMalloc1(eps->mpd,&gamma);CHKERRQ(ierr);

  kini = eps->nini;
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;
    nv = PetscMin(eps->nconv+eps->mpd,ncv);
    ierr = DSSetDimensions(eps->ds,nv,0,eps->nconv,0);CHKERRQ(ierr);
    for (;kini<nv;kini++) { /* Generate more initial vectors if necessary */
      ierr = BVSetRandomColumn(eps->V,kini);CHKERRQ(ierr);
      ierr = BVOrthonormalizeColumn(eps->V,kini,PETSC_TRUE,NULL,NULL);CHKERRQ(ierr);
    }
    reset = (eps->its>1 && (eps->its-1)%ctx->nrest==0)? PETSC_TRUE: PETSC_FALSE;

    if (reset) {
      /* Prevent BVDotVec below to use B-product, restored at the end */
      ierr = BVSetMatrix(eps->V,NULL,PETSC_FALSE);CHKERRQ(ierr);

      /* Compute Rayleigh quotient */
      ierr = BVSetActiveColumns(eps->V,eps->nconv,nv);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(ctx->AV,0,nv-eps->nconv);CHKERRQ(ierr);
      ierr = BVMatMult(eps->V,A,ctx->AV);CHKERRQ(ierr);
      ierr = DSGetArray(eps->ds,DS_MAT_A,&C);CHKERRQ(ierr);
      for (i=eps->nconv;i<nv;i++) {
        ierr = BVSetActiveColumns(eps->V,eps->nconv,i+1);CHKERRQ(ierr);
        ierr = BVGetColumn(ctx->AV,i-eps->nconv,&av);CHKERRQ(ierr);
        ierr = BVDotVec(eps->V,av,C+eps->nconv+i*ld);CHKERRQ(ierr);
        ierr = BVRestoreColumn(ctx->AV,i-eps->nconv,&av);CHKERRQ(ierr);
        for (j=eps->nconv;j<i-1;j++) C[i+j*ld] = PetscConj(C[j+i*ld]);
      }
      ierr = DSRestoreArray(eps->ds,DS_MAT_A,&C);CHKERRQ(ierr);
      ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);

      /* Solve projected problem */
      ierr = DSSolve(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);
      ierr = DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
      ierr = DSSynchronize(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);

      /* Update vectors V(:,idx) = V * Y(:,idx) */
      ierr = DSGetMat(eps->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
      ierr = BVMultInPlace(eps->V,Q,eps->nconv,nv);CHKERRQ(ierr);
      ierr = ExtractSubmatrix(Q,eps->nconv,&Q1);CHKERRQ(ierr);
      ierr = BVMultInPlace(ctx->AV,Q1,0,nv-eps->nconv);CHKERRQ(ierr);
      ierr = MatDestroy(&Q);CHKERRQ(ierr);
      ierr = MatDestroy(&Q1);CHKERRQ(ierr);
      if (B) { ierr = BVSetMatrix(eps->V,B,PETSC_FALSE);CHKERRQ(ierr); }
    } else {
      /* No need to do Rayleigh-Ritz, just take diag(V'*A*V) */
      for (i=eps->nconv;i<nv;i++) {
        ierr = BVGetColumn(eps->V,i,&v);CHKERRQ(ierr);
        ierr = BVGetColumn(ctx->AV,i-eps->nconv,&av);CHKERRQ(ierr);
        ierr = MatMult(A,v,av);CHKERRQ(ierr);
        ierr = VecDot(av,v,eps->eigr+i);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->V,i,&v);CHKERRQ(ierr);
        ierr = BVRestoreColumn(ctx->AV,i-eps->nconv,&av);CHKERRQ(ierr);
      }
    }

    /* Compute gradient and check convergence */
    k = -1;
    for (i=eps->nconv;i<nv;i++) {
      ierr = BVGetColumn(eps->V,i,&v);CHKERRQ(ierr);
      ierr = BVGetColumn(ctx->AV,i-eps->nconv,&av);CHKERRQ(ierr);
      ierr = BVGetColumn(ctx->G,i-eps->nconv,&p);CHKERRQ(ierr);
      if (B) {
        ierr = BVGetColumn(ctx->W,i-eps->nconv,&bv);CHKERRQ(ierr);
        ierr = MatMult(B,v,bv);CHKERRQ(ierr);
        ierr = VecWAXPY(p,-eps->eigr[i],bv,av);CHKERRQ(ierr);
        ierr = BVRestoreColumn(ctx->W,i-eps->nconv,&bv);CHKERRQ(ierr);
      } else {
        ierr = VecWAXPY(p,-eps->eigr[i],v,av);CHKERRQ(ierr);
      }
      ierr = BVRestoreColumn(eps->V,i,&v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(ctx->AV,i-eps->nconv,&av);CHKERRQ(ierr);
      ierr = VecNorm(p,NORM_2,&resnorm);CHKERRQ(ierr);
      ierr = BVRestoreColumn(ctx->G,i-eps->nconv,&p);CHKERRQ(ierr);
      ierr = (*eps->converged)(eps,eps->eigr[i],0.0,resnorm,&eps->errest[i],eps->convergedctx);CHKERRQ(ierr);
      if (k==-1 && eps->errest[i] >= eps->tol) k = i;
    }
    if (k==-1) k = nv;
    ierr = (*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);

    /* The next lines are necessary to avoid DS zeroing eigr */
    ierr = DSGetArray(eps->ds,DS_MAT_A,&C);CHKERRQ(ierr);
    for (i=eps->nconv;i<k;i++) C[i+i*ld] = eps->eigr[i];
    ierr = DSRestoreArray(eps->ds,DS_MAT_A,&C);CHKERRQ(ierr);

    if (eps->reason == EPS_CONVERGED_ITERATING) {

      /* Search direction */
      for (i=0;i<nv-eps->nconv;i++) {
        ierr = BVGetColumn(ctx->G,i,&v);CHKERRQ(ierr);
        ierr = STApply(eps->st,v,w);CHKERRQ(ierr);
        ierr = VecDot(w,v,&g);CHKERRQ(ierr);
        ierr = BVRestoreColumn(ctx->G,i,&v);CHKERRQ(ierr);
        beta = (!reset && eps->its>1)? g/gamma[i]: 0.0;
        gamma[i] = g;
        ierr = BVGetColumn(ctx->P,i,&v);CHKERRQ(ierr);
        ierr = VecAXPBY(v,1.0,beta,w);CHKERRQ(ierr);
        if (i+eps->nconv>0) {
          ierr = BVSetActiveColumns(eps->V,0,i+eps->nconv);CHKERRQ(ierr);
          ierr = BVOrthogonalizeVec(eps->V,v,NULL,NULL,NULL);CHKERRQ(ierr);
        }
        ierr = BVRestoreColumn(ctx->P,i,&v);CHKERRQ(ierr);
      }

      /* Minimization problem */
      for (i=eps->nconv;i<nv;i++) {
        ierr = BVGetColumn(eps->V,i,&v);CHKERRQ(ierr);
        ierr = BVGetColumn(ctx->AV,i-eps->nconv,&av);CHKERRQ(ierr);
        ierr = BVGetColumn(ctx->P,i-eps->nconv,&p);CHKERRQ(ierr);
        ierr = VecDot(av,v,&nu);CHKERRQ(ierr);
        ierr = VecDot(av,p,&pax);CHKERRQ(ierr);
        ierr = MatMult(A,p,w);CHKERRQ(ierr);
        ierr = VecDot(w,p,&pap);CHKERRQ(ierr);
        if (B) {
          ierr = BVGetColumn(ctx->W,i-eps->nconv,&bv);CHKERRQ(ierr);
          ierr = VecDot(bv,v,&mu);CHKERRQ(ierr);
          ierr = VecDot(bv,p,&pbx);CHKERRQ(ierr);
          ierr = BVRestoreColumn(ctx->W,i-eps->nconv,&bv);CHKERRQ(ierr);
          ierr = MatMult(B,p,w);CHKERRQ(ierr);
          ierr = VecDot(w,p,&pbp);CHKERRQ(ierr);
        } else {
          ierr = VecDot(v,v,&mu);CHKERRQ(ierr);
          ierr = VecDot(v,p,&pbx);CHKERRQ(ierr);
          ierr = VecDot(p,p,&pbp);CHKERRQ(ierr);
        }
        ierr = BVRestoreColumn(ctx->AV,i-eps->nconv,&av);CHKERRQ(ierr);
        a = PetscRealPart(pap*pbx-pax*pbp);
        b = PetscRealPart(nu*pbp-mu*pap);
        c = PetscRealPart(mu*pax-nu*pbx);
        t = PetscMax(PetscMax(PetscAbsReal(a),PetscAbsReal(b)),PetscAbsReal(c));
        if (t!=0.0) { a /= t; b /= t; c /= t; }
        disc = b*b-4.0*a*c;
        d = PetscSqrtReal(PetscAbsReal(disc));
        if (b>=0.0 && a!=0.0) alpha = (b+d)/(2.0*a);
        else if (b!=d) alpha = 2.0*c/(b-d);
        else alpha = 0;
        /* Next iterate */
        if (alpha!=0.0) {
          ierr = VecAXPY(v,alpha,p);CHKERRQ(ierr);
        }
        ierr = BVRestoreColumn(eps->V,i,&v);CHKERRQ(ierr);
        ierr = BVRestoreColumn(ctx->P,i-eps->nconv,&p);CHKERRQ(ierr);
        ierr = BVOrthonormalizeColumn(eps->V,i,PETSC_TRUE,NULL,NULL);CHKERRQ(ierr);
      }
    }

    ierr = EPSMonitor(eps,eps->its,k,eps->eigr,eps->eigi,eps->errest,nv);CHKERRQ(ierr);
    eps->nconv = k;
  }

  ierr = PetscFree(gamma);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSRQCGSetReset_RQCG(EPS eps,PetscInt nrest)
{
  EPS_RQCG *ctx = (EPS_RQCG*)eps->data;

  PetscFunctionBegin;
  if (nrest==PETSC_DEFAULT) {
    ctx->nrest = 0;
    eps->state = EPS_STATE_INITIAL;
  } else {
    if (nrest<=0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Reset parameter must be >0");
    ctx->nrest = nrest;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSRQCGSetReset - Sets the reset parameter of the RQCG iteration. Every
   nrest iterations, the solver performs a Rayleigh-Ritz projection step.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  nrest - the number of iterations between resets

   Options Database Key:
.  -eps_rqcg_reset - Sets the reset parameter

   Level: advanced

.seealso: EPSRQCGGetReset()
@*/
PetscErrorCode EPSRQCGSetReset(EPS eps,PetscInt nrest)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nrest,2);
  ierr = PetscTryMethod(eps,"EPSRQCGSetReset_C",(EPS,PetscInt),(eps,nrest));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSRQCGGetReset_RQCG(EPS eps,PetscInt *nrest)
{
  EPS_RQCG *ctx = (EPS_RQCG*)eps->data;

  PetscFunctionBegin;
  *nrest = ctx->nrest;
  PetscFunctionReturn(0);
}

/*@
   EPSRQCGGetReset - Gets the reset parameter used in the RQCG method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  nrest - the reset parameter

   Level: advanced

.seealso: EPSRQCGSetReset()
@*/
PetscErrorCode EPSRQCGGetReset(EPS eps,PetscInt *nrest)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(nrest,2);
  ierr = PetscUseMethod(eps,"EPSRQCGGetReset_C",(EPS,PetscInt*),(eps,nrest));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_RQCG(EPS eps)
{
  PetscErrorCode ierr;
  EPS_RQCG       *ctx = (EPS_RQCG*)eps->data;

  PetscFunctionBegin;
  ierr = BVDestroy(&ctx->AV);CHKERRQ(ierr);
  ierr = BVDestroy(&ctx->W);CHKERRQ(ierr);
  ierr = BVDestroy(&ctx->P);CHKERRQ(ierr);
  ierr = BVDestroy(&ctx->G);CHKERRQ(ierr);
  ctx->allocsize = 0;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_RQCG(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      flg;
  PetscInt       nrest;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS RQCG Options");CHKERRQ(ierr);

    ierr = PetscOptionsInt("-eps_rqcg_reset","Reset parameter","EPSRQCGSetReset",20,&nrest,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSRQCGSetReset(eps,nrest);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_RQCG(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSRQCGSetReset_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSRQCGGetReset_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_RQCG(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  EPS_RQCG       *ctx = (EPS_RQCG*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  reset every %D iterations\n",ctx->nrest);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_RQCG(EPS eps)
{
  EPS_RQCG       *rqcg;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&rqcg);CHKERRQ(ierr);
  eps->data = (void*)rqcg;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_RQCG;
  eps->ops->setup          = EPSSetUp_RQCG;
  eps->ops->setfromoptions = EPSSetFromOptions_RQCG;
  eps->ops->destroy        = EPSDestroy_RQCG;
  eps->ops->reset          = EPSReset_RQCG;
  eps->ops->view           = EPSView_RQCG;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_GMRES;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSRQCGSetReset_C",EPSRQCGSetReset_RQCG);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSRQCGGetReset_C",EPSRQCGGetReset_RQCG);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

