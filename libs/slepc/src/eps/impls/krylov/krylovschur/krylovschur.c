/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "krylovschur"

   Method: Krylov-Schur

   Algorithm:

       Single-vector Krylov-Schur method for non-symmetric problems,
       including harmonic extraction.

   References:

       [1] "Krylov-Schur Methods in SLEPc", SLEPc Technical Report STR-7,
           available at https://slepc.upv.es.

       [2] G.W. Stewart, "A Krylov-Schur Algorithm for Large Eigenproblems",
           SIAM J. Matrix Anal. App. 23(3):601-614, 2001.

       [3] "Practical Implementation of Harmonic Krylov-Schur", SLEPc Technical
            Report STR-9, available at https://slepc.upv.es.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include "krylovschur.h"

PetscErrorCode EPSGetArbitraryValues(EPS eps,PetscScalar *rr,PetscScalar *ri)
{
  PetscErrorCode ierr;
  PetscInt       i,newi,ld,n,l;
  Vec            xr=eps->work[0],xi=eps->work[1];
  PetscScalar    re,im,*Zr,*Zi,*X;

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  ierr = DSGetDimensions(eps->ds,&n,NULL,&l,NULL,NULL);CHKERRQ(ierr);
  for (i=l;i<n;i++) {
    re = eps->eigr[i];
    im = eps->eigi[i];
    ierr = STBackTransform(eps->st,1,&re,&im);CHKERRQ(ierr);
    newi = i;
    ierr = DSVectors(eps->ds,DS_MAT_X,&newi,NULL);CHKERRQ(ierr);
    ierr = DSGetArray(eps->ds,DS_MAT_X,&X);CHKERRQ(ierr);
    Zr = X+i*ld;
    if (newi==i+1) Zi = X+newi*ld;
    else Zi = NULL;
    ierr = EPSComputeRitzVector(eps,Zr,Zi,eps->V,xr,xi);CHKERRQ(ierr);
    ierr = DSRestoreArray(eps->ds,DS_MAT_X,&X);CHKERRQ(ierr);
    ierr = (*eps->arbitrary)(re,im,xr,xi,rr+i,ri+i,eps->arbitraryctx);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode EstimateRange(Mat A,PetscReal *left,PetscReal *right)
{
  PetscErrorCode ierr;
  PetscInt       nconv;
  PetscScalar    eig0;
  EPS            eps;

  PetscFunctionBegin;
  *left = 0.0; *right = 0.0;
  ierr = EPSCreate(PetscObjectComm((PetscObject)A),&eps);CHKERRQ(ierr);
  ierr = EPSSetOptionsPrefix(eps,"eps_filter_");CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,1e-3,50);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  if (nconv>0) {
    ierr = EPSGetEigenvalue(eps,0,&eig0,NULL);CHKERRQ(ierr);
  } else eig0 = eps->eigr[0];
  *left = PetscRealPart(eig0);
  ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  if (nconv>0) {
    ierr = EPSGetEigenvalue(eps,0,&eig0,NULL);CHKERRQ(ierr);
  } else eig0 = eps->eigr[0];
  *right = PetscRealPart(eig0);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSSetUp_KrylovSchur_Filter(EPS eps)
{
  PetscErrorCode ierr;
  SlepcSC        sc;
  PetscReal      rleft,rright;
  Mat            A;

  PetscFunctionBegin;
  if (eps->intb >= PETSC_MAX_REAL && eps->inta <= PETSC_MIN_REAL) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"The defined computational interval should have at least one of their sides bounded");
  if (!eps->ishermitian) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Polynomial filter only available for symmetric/Hermitian eigenproblems");
  if (eps->isgeneralized) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Polynomial filters not available for generalized eigenproblems");
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs cannot be used with polynomial filters");
  if (eps->tol==PETSC_DEFAULT) eps->tol = SLEPC_DEFAULT_TOL*1e-2;  /* use tighter tolerance */
  ierr = STFilterSetInterval(eps->st,eps->inta,eps->intb);CHKERRQ(ierr);
  ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
  ierr = STFilterGetRange(eps->st,&rleft,&rright);CHKERRQ(ierr);
  if (!rleft && !rright) {
    ierr = EstimateRange(A,&rleft,&rright);CHKERRQ(ierr);
    ierr = PetscInfo2(eps,"Setting eigenvalue range to [%g,%g]\n",(double)rleft,(double)rright);CHKERRQ(ierr);
    ierr = STFilterSetRange(eps->st,rleft,rright);CHKERRQ(ierr);
  }
  if (eps->ncv==PETSC_DEFAULT && eps->nev==1) eps->nev = 40;  /* user did not provide nev estimation */
  ierr = EPSSetDimensions_Default(eps,eps->nev,&eps->ncv,&eps->mpd);CHKERRQ(ierr);
  if (eps->ncv>eps->nev+eps->mpd) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv must not be larger than nev+mpd");
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);

  ierr = DSGetSlepcSC(eps->ds,&sc);CHKERRQ(ierr);
  sc->rg            = NULL;
  sc->comparison    = SlepcCompareLargestReal;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetUp_KrylovSchur(EPS eps)
{
  PetscErrorCode    ierr;
  PetscReal         eta;
  PetscBool         isfilt=PETSC_FALSE;
  BVOrthogType      otype;
  BVOrthogBlockType obtype;
  EPS_KRYLOVSCHUR   *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  enum { EPS_KS_DEFAULT,EPS_KS_SYMM,EPS_KS_SLICE,EPS_KS_FILTER,EPS_KS_INDEF,EPS_KS_TWOSIDED } variant;

  PetscFunctionBegin;
  /* spectrum slicing requires special treatment of default values */
  if (eps->which==EPS_ALL) {
    ierr = PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt);CHKERRQ(ierr);
    if (isfilt) {
      ierr = EPSSetUp_KrylovSchur_Filter(eps);CHKERRQ(ierr);
    } else {
      ierr = EPSSetUp_KrylovSchur_Slice(eps);CHKERRQ(ierr);
    }
  } else {
    ierr = EPSSetDimensions_Default(eps,eps->nev,&eps->ncv,&eps->mpd);CHKERRQ(ierr);
    if (eps->ncv>eps->nev+eps->mpd) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv must not be larger than nev+mpd");
    if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
    if (!eps->which) { ierr = EPSSetWhichEigenpairs_Default(eps);CHKERRQ(ierr); }
  }
  if (!ctx->lock && eps->mpd<eps->ncv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Should not use mpd parameter in non-locking variant");

  if (eps->isgeneralized && eps->ishermitian && !eps->ispositive && eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not implemented for indefinite problems");
  if (eps->ishermitian && eps->ispositive && (eps->which==EPS_LARGEST_IMAGINARY || eps->which==EPS_SMALLEST_IMAGINARY)) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Wrong value of eps->which");

  if (!eps->extraction) {
    ierr = EPSSetExtraction(eps,EPS_RITZ);CHKERRQ(ierr);
  } else if (eps->extraction!=EPS_RITZ && eps->extraction!=EPS_HARMONIC) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");

  if (!ctx->keep) ctx->keep = 0.5;

  ierr = EPSAllocateSolution(eps,1);CHKERRQ(ierr);
  ierr = EPS_SetInnerProduct(eps);CHKERRQ(ierr);
  if (eps->arbitrary) {
    ierr = EPSSetWorkVecs(eps,2);CHKERRQ(ierr);
  } else if (eps->ishermitian && !eps->ispositive){
    ierr = EPSSetWorkVecs(eps,1);CHKERRQ(ierr);
  }

  /* dispatch solve method */
  if (eps->ishermitian) {
    if (eps->which==EPS_ALL) {
      if (eps->isgeneralized && !eps->ispositive) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Spectrum slicing not implemented for indefinite problems");
      else variant = isfilt? EPS_KS_FILTER: EPS_KS_SLICE;
    } else if (eps->isgeneralized && !eps->ispositive) {
      variant = EPS_KS_INDEF;
    } else {
      switch (eps->extraction) {
        case EPS_RITZ:     variant = EPS_KS_SYMM; break;
        case EPS_HARMONIC: variant = EPS_KS_DEFAULT; break;
        default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");
      }
    }
  } else if (eps->twosided) {
    variant = EPS_KS_TWOSIDED;
  } else {
    switch (eps->extraction) {
      case EPS_RITZ:     variant = EPS_KS_DEFAULT; break;
      case EPS_HARMONIC: variant = EPS_KS_DEFAULT; break;
      default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");
    }
  }
  switch (variant) {
    case EPS_KS_DEFAULT:
      eps->ops->solve = EPSSolve_KrylovSchur_Default;
      eps->ops->computevectors = EPSComputeVectors_Schur;
      ierr = DSSetType(eps->ds,DSNHEP);CHKERRQ(ierr);
      ierr = DSAllocate(eps->ds,eps->ncv+1);CHKERRQ(ierr);
      break;
    case EPS_KS_SYMM:
    case EPS_KS_FILTER:
      eps->ops->solve = EPSSolve_KrylovSchur_Symm;
      eps->ops->computevectors = EPSComputeVectors_Hermitian;
      ierr = DSSetType(eps->ds,DSHEP);CHKERRQ(ierr);
      ierr = DSSetCompact(eps->ds,PETSC_TRUE);CHKERRQ(ierr);
      ierr = DSSetExtraRow(eps->ds,PETSC_TRUE);CHKERRQ(ierr);
      ierr = DSAllocate(eps->ds,eps->ncv+1);CHKERRQ(ierr);
      break;
    case EPS_KS_SLICE:
      if (eps->stopping!=EPSStoppingBasic) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Spectrum slicing does not support user-defined stopping test");
      eps->ops->solve = EPSSolve_KrylovSchur_Slice;
      eps->ops->computevectors = EPSComputeVectors_Slice;
      break;
    case EPS_KS_INDEF:
      eps->ops->solve = EPSSolve_KrylovSchur_Indefinite;
      eps->ops->computevectors = EPSComputeVectors_Indefinite;
      ierr = DSSetType(eps->ds,DSGHIEP);CHKERRQ(ierr);
      ierr = DSSetCompact(eps->ds,PETSC_TRUE);CHKERRQ(ierr);
      ierr = DSAllocate(eps->ds,eps->ncv+1);CHKERRQ(ierr);
      /* force reorthogonalization for pseudo-Lanczos */
      ierr = BVGetOrthogonalization(eps->V,&otype,NULL,&eta,&obtype);CHKERRQ(ierr);
      ierr = BVSetOrthogonalization(eps->V,otype,BV_ORTHOG_REFINE_ALWAYS,eta,obtype);CHKERRQ(ierr);
      break;
    case EPS_KS_TWOSIDED:
      eps->ops->solve = EPSSolve_KrylovSchur_TwoSided;
      eps->ops->computevectors = EPSComputeVectors_Schur;
      ierr = DSSetType(eps->ds,DSNHEP);CHKERRQ(ierr);
      ierr = DSAllocate(eps->ds,eps->ncv+1);CHKERRQ(ierr);
      ierr = DSSetType(eps->dsts,DSNHEP);CHKERRQ(ierr);
      ierr = DSAllocate(eps->dsts,eps->ncv+1);CHKERRQ(ierr);
      break;
    default: SETERRQ(PetscObjectComm((PetscObject)eps),1,"Unexpected error");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_KrylovSchur_Default(EPS eps)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,j,*pj,k,l,nv,ld,nconv;
  Mat             U,Op;
  PetscScalar     *S,*Q,*g;
  PetscReal       beta,gamma=1.0;
  PetscBool       breakdown,harmonic;

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  harmonic = (eps->extraction==EPS_HARMONIC || eps->extraction==EPS_REFINED_HARMONIC)?PETSC_TRUE:PETSC_FALSE;
  if (harmonic) { ierr = PetscMalloc1(ld,&g);CHKERRQ(ierr); }
  if (eps->arbitrary) pj = &j;
  else pj = NULL;

  /* Get the starting Arnoldi vector */
  ierr = EPSGetStartVector(eps,0,NULL);CHKERRQ(ierr);
  l = 0;

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Arnoldi factorization */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    ierr = STGetOperator(eps->st,&Op);CHKERRQ(ierr);
    ierr = DSGetArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
    ierr = BVMatArnoldi(eps->V,Op,S,ld,eps->nconv+l,&nv,&beta,&breakdown);CHKERRQ(ierr);
    ierr = DSRestoreArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
    ierr = STRestoreOperator(eps->st,&Op);CHKERRQ(ierr);
    ierr = DSSetDimensions(eps->ds,nv,0,eps->nconv,eps->nconv+l);CHKERRQ(ierr);
    if (l==0) {
      ierr = DSSetState(eps->ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
    } else {
      ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
    }
    ierr = BVSetActiveColumns(eps->V,eps->nconv,nv);CHKERRQ(ierr);

    /* Compute translation of Krylov decomposition if harmonic extraction used */
    if (harmonic) {
      ierr = DSTranslateHarmonic(eps->ds,eps->target,beta,PETSC_FALSE,g,&gamma);CHKERRQ(ierr);
    }

    /* Solve projected problem */
    ierr = DSSolve(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);
    if (eps->arbitrary) {
      ierr = EPSGetArbitraryValues(eps,eps->rr,eps->ri);CHKERRQ(ierr);
      j=1;
    }
    ierr = DSSort(eps->ds,eps->eigr,eps->eigi,eps->rr,eps->ri,pj);CHKERRQ(ierr);
    ierr = DSSynchronize(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);

    /* Check convergence */
    ierr = EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta,0.0,gamma,&k);CHKERRQ(ierr);
    ierr = (*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);
    nconv = k;

    /* Update l */
    if (eps->reason != EPS_CONVERGED_ITERATING || breakdown || k==nv) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
#if !defined(PETSC_USE_COMPLEX)
      ierr = DSGetArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
      if (S[k+l+(k+l-1)*ld] != 0.0) {
        if (k+l<nv-1) l = l+1;
        else l = l-1;
      }
      ierr = DSRestoreArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
#endif
    }
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      if (breakdown || k==nv) {
        /* Start a new Arnoldi factorization */
        ierr = PetscInfo2(eps,"Breakdown in Krylov-Schur method (it=%D norm=%g)\n",eps->its,(double)beta);CHKERRQ(ierr);
        if (k<eps->nev) {
          ierr = EPSGetStartVector(eps,k,&breakdown);CHKERRQ(ierr);
          if (breakdown) {
            eps->reason = EPS_DIVERGED_BREAKDOWN;
            ierr = PetscInfo(eps,"Unable to generate more start vectors\n");CHKERRQ(ierr);
          }
        }
      } else {
        /* Undo translation of Krylov decomposition */
        if (harmonic) {
          ierr = DSSetDimensions(eps->ds,nv,0,k,l);CHKERRQ(ierr);
          ierr = DSTranslateHarmonic(eps->ds,0.0,beta,PETSC_TRUE,g,&gamma);CHKERRQ(ierr);
          /* gamma u^ = u - U*g~ */
          ierr = BVSetActiveColumns(eps->V,0,nv);CHKERRQ(ierr);
          ierr = BVMultColumn(eps->V,-1.0,1.0,nv,g);CHKERRQ(ierr);
          ierr = BVScaleColumn(eps->V,nv,1.0/gamma);CHKERRQ(ierr);
          ierr = BVSetActiveColumns(eps->V,eps->nconv,nv);CHKERRQ(ierr);
        }
        /* Prepare the Rayleigh quotient for restart */
        ierr = DSGetArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
        ierr = DSGetArray(eps->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
        for (i=k;i<k+l;i++) {
          S[k+l+i*ld] = Q[nv-1+i*ld]*beta*gamma;
        }
        ierr = DSRestoreArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
        ierr = DSRestoreArray(eps->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
      }
    }
    /* Update the corresponding vectors V(:,idx) = V*Q(:,idx) */
    ierr = DSGetMat(eps->ds,DS_MAT_Q,&U);CHKERRQ(ierr);
    ierr = BVMultInPlace(eps->V,U,eps->nconv,k+l);CHKERRQ(ierr);
    ierr = MatDestroy(&U);CHKERRQ(ierr);

    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      ierr = BVCopyColumn(eps->V,nv,k+l);CHKERRQ(ierr);
    }
    eps->nconv = k;
    ierr = EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,nv);CHKERRQ(ierr);
  }

  if (harmonic) { ierr = PetscFree(g);CHKERRQ(ierr); }
  /* truncate Schur decomposition and change the state to raw so that
     DSVectors() computes eigenvectors from scratch */
  ierr = DSSetDimensions(eps->ds,eps->nconv,0,0,0);CHKERRQ(ierr);
  ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurSetRestart_KrylovSchur(EPS eps,PetscReal keep)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (keep==PETSC_DEFAULT) ctx->keep = 0.5;
  else {
    if (keep<0.1 || keep>0.9) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The keep argument %g must be in the range [0.1,0.9]",keep);
    ctx->keep = keep;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurSetRestart - Sets the restart parameter for the Krylov-Schur
   method, in particular the proportion of basis vectors that must be kept
   after restart.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  keep - the number of vectors to be kept at restart

   Options Database Key:
.  -eps_krylovschur_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,0.9]. The default is 0.5.

   Level: advanced

.seealso: EPSKrylovSchurGetRestart()
@*/
PetscErrorCode EPSKrylovSchurSetRestart(EPS eps,PetscReal keep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,keep,2);
  ierr = PetscTryMethod(eps,"EPSKrylovSchurSetRestart_C",(EPS,PetscReal),(eps,keep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetRestart_KrylovSchur(EPS eps,PetscReal *keep)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  *keep = ctx->keep;
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurGetRestart - Gets the restart parameter used in the
   Krylov-Schur method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  keep - the restart parameter

   Level: advanced

.seealso: EPSKrylovSchurSetRestart()
@*/
PetscErrorCode EPSKrylovSchurGetRestart(EPS eps,PetscReal *keep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidRealPointer(keep,2);
  ierr = PetscUseMethod(eps,"EPSKrylovSchurGetRestart_C",(EPS,PetscReal*),(eps,keep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurSetLocking_KrylovSchur(EPS eps,PetscBool lock)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurSetLocking - Choose between locking and non-locking variants of
   the Krylov-Schur method.

   Logically Collective on eps

   Input Parameters:
+  eps  - the eigenproblem solver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -eps_krylovschur_locking - Sets the locking flag

   Notes:
   The default is to lock converged eigenpairs when the method restarts.
   This behaviour can be changed so that all directions are kept in the
   working subspace even if already converged to working accuracy (the
   non-locking variant).

   Level: advanced

.seealso: EPSKrylovSchurGetLocking()
@*/
PetscErrorCode EPSKrylovSchurSetLocking(EPS eps,PetscBool lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,lock,2);
  ierr = PetscTryMethod(eps,"EPSKrylovSchurSetLocking_C",(EPS,PetscBool),(eps,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetLocking_KrylovSchur(EPS eps,PetscBool *lock)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurGetLocking - Gets the locking flag used in the Krylov-Schur
   method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: EPSKrylovSchurSetLocking()
@*/
PetscErrorCode EPSKrylovSchurGetLocking(EPS eps,PetscBool *lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(lock,2);
  ierr = PetscUseMethod(eps,"EPSKrylovSchurGetLocking_C",(EPS,PetscBool*),(eps,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurSetPartitions_KrylovSchur(EPS eps,PetscInt npart)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscMPIInt     size;

  PetscFunctionBegin;
  if (ctx->npart!=npart) {
    if (ctx->commset) { ierr = PetscSubcommDestroy(&ctx->subc);CHKERRQ(ierr); }
    ierr = EPSDestroy(&ctx->eps);CHKERRQ(ierr);
  }
  if (npart == PETSC_DEFAULT || npart == PETSC_DECIDE) {
    ctx->npart = 1;
  } else {
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)eps),&size);CHKERRQ(ierr);
    if (npart<1 || npart>size) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of npart");
    ctx->npart = npart;
  }
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurSetPartitions - Sets the number of partitions for the
   case of doing spectrum slicing for a computational interval with the
   communicator split in several sub-communicators.

   Logically Collective on eps

   Input Parameters:
+  eps   - the eigenproblem solver context
-  npart - number of partitions

   Options Database Key:
.  -eps_krylovschur_partitions <npart> - Sets the number of partitions

   Notes:
   By default, npart=1 so all processes in the communicator participate in
   the processing of the whole interval. If npart>1 then the interval is
   divided into npart subintervals, each of them being processed by a
   subset of processes.

   The interval is split proportionally unless the separation points are
   specified with EPSKrylovSchurSetSubintervals().

   Level: advanced

.seealso: EPSKrylovSchurSetSubintervals(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurSetPartitions(EPS eps,PetscInt npart)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,npart,2);
  ierr = PetscTryMethod(eps,"EPSKrylovSchurSetPartitions_C",(EPS,PetscInt),(eps,npart));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetPartitions_KrylovSchur(EPS eps,PetscInt *npart)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  *npart  = ctx->npart;
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurGetPartitions - Gets the number of partitions of the
   communicator in case of spectrum slicing.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  npart - number of partitions

   Level: advanced

.seealso: EPSKrylovSchurSetPartitions()
@*/
PetscErrorCode EPSKrylovSchurGetPartitions(EPS eps,PetscInt *npart)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(npart,2);
  ierr = PetscUseMethod(eps,"EPSKrylovSchurGetPartitions_C",(EPS,PetscInt*),(eps,npart));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurSetDetectZeros_KrylovSchur(EPS eps,PetscBool detect)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  ctx->detect = detect;
  eps->state  = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurSetDetectZeros - Sets a flag to enforce detection of
   zeros during the factorizations throughout the spectrum slicing computation.

   Logically Collective on eps

   Input Parameters:
+  eps    - the eigenproblem solver context
-  detect - check for zeros

   Options Database Key:
.  -eps_krylovschur_detect_zeros - Check for zeros; this takes an optional
   bool value (0/1/no/yes/true/false)

   Notes:
   A zero in the factorization indicates that a shift coincides with an eigenvalue.

   This flag is turned off by default, and may be necessary in some cases,
   especially when several partitions are being used. This feature currently
   requires an external package for factorizations with support for zero
   detection, e.g. MUMPS.

   Level: advanced

.seealso: EPSKrylovSchurSetPartitions(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurSetDetectZeros(EPS eps,PetscBool detect)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,detect,2);
  ierr = PetscTryMethod(eps,"EPSKrylovSchurSetDetectZeros_C",(EPS,PetscBool),(eps,detect));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetDetectZeros_KrylovSchur(EPS eps,PetscBool *detect)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  *detect = ctx->detect;
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurGetDetectZeros - Gets the flag that enforces zero detection
   in spectrum slicing.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  detect - whether zeros detection is enforced during factorizations

   Level: advanced

.seealso: EPSKrylovSchurSetDetectZeros()
@*/
PetscErrorCode EPSKrylovSchurGetDetectZeros(EPS eps,PetscBool *detect)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(detect,2);
  ierr = PetscUseMethod(eps,"EPSKrylovSchurGetDetectZeros_C",(EPS,PetscBool*),(eps,detect));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurSetDimensions_KrylovSchur(EPS eps,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (nev<1) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nev. Must be > 0");
  ctx->nev = nev;
  if (ncv == PETSC_DECIDE || ncv == PETSC_DEFAULT) {
    ctx->ncv = PETSC_DEFAULT;
  } else {
    if (ncv<1) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    ctx->ncv = ncv;
  }
  if (mpd == PETSC_DECIDE || mpd == PETSC_DEFAULT) {
    ctx->mpd = PETSC_DEFAULT;
  } else {
    if (mpd<1) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    ctx->mpd = mpd;
  }
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurSetDimensions - Sets the dimensions used for each subsolve
   step in case of doing spectrum slicing for a computational interval.
   The meaning of the parameters is the same as in EPSSetDimensions().

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
.  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the subsolve
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Key:
+  -eps_krylovschur_nev <nev> - Sets the number of eigenvalues
.  -eps_krylovschur_ncv <ncv> - Sets the dimension of the subspace
-  -eps_krylovschur_mpd <mpd> - Sets the maximum projected dimension

   Level: advanced

.seealso: EPSKrylovSchurGetDimensions(), EPSSetDimensions(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurSetDimensions(EPS eps,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nev,2);
  PetscValidLogicalCollectiveInt(eps,ncv,3);
  PetscValidLogicalCollectiveInt(eps,mpd,4);
  ierr = PetscTryMethod(eps,"EPSKrylovSchurSetDimensions_C",(EPS,PetscInt,PetscInt,PetscInt),(eps,nev,ncv,mpd));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetDimensions_KrylovSchur(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (nev) *nev = ctx->nev;
  if (ncv) *ncv = ctx->ncv;
  if (mpd) *mpd = ctx->mpd;
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurGetDimensions - Gets the dimensions used for each subsolve
   step in case of doing spectrum slicing for a computational interval.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the subsolve
-  mpd - the maximum dimension allowed for the projected problem

   Level: advanced

.seealso: EPSKrylovSchurSetDimensions()
@*/
PetscErrorCode EPSKrylovSchurGetDimensions(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = PetscUseMethod(eps,"EPSKrylovSchurGetDimensions_C",(EPS,PetscInt*,PetscInt*,PetscInt*),(eps,nev,ncv,mpd));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurSetSubintervals_KrylovSchur(EPS eps,PetscReal* subint)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i;

  PetscFunctionBegin;
  if (subint[0]!=eps->inta || subint[ctx->npart]!=eps->intb) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"First and last values must match the endpoints of EPSSetInterval()");
  for (i=0;i<ctx->npart;i++) if (subint[i]>subint[i+1]) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Array must contain values in ascending order");
  if (ctx->subintervals) { ierr = PetscFree(ctx->subintervals);CHKERRQ(ierr); }
  ierr = PetscMalloc1(ctx->npart+1,&ctx->subintervals);CHKERRQ(ierr);
  for (i=0;i<ctx->npart+1;i++) ctx->subintervals[i] = subint[i];
  ctx->subintset = PETSC_TRUE;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@C
   EPSKrylovSchurSetSubintervals - Sets the points that delimit the
   subintervals to be used in spectrum slicing with several partitions.

   Logically Collective on eps

   Input Parameters:
+  eps    - the eigenproblem solver context
-  subint - array of real values specifying subintervals

   Notes:
   This function must be called after EPSKrylovSchurSetPartitions(). For npart
   partitions, the argument subint must contain npart+1 real values sorted in
   ascending order: subint_0, subint_1, ..., subint_npart, where the first
   and last values must coincide with the interval endpoints set with
   EPSSetInterval().

   The subintervals are then defined by two consecutive points: [subint_0,subint_1],
   [subint_1,subint_2], and so on.

   Level: advanced

.seealso: EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubintervals(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurSetSubintervals(EPS eps,PetscReal *subint)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = PetscTryMethod(eps,"EPSKrylovSchurSetSubintervals_C",(EPS,PetscReal*),(eps,subint));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetSubintervals_KrylovSchur(EPS eps,PetscReal **subint)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i;

  PetscFunctionBegin;
  if (!ctx->subintset) {
    if (!eps->state) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
    if (!ctx->sr) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  }
  ierr = PetscMalloc1(ctx->npart+1,subint);CHKERRQ(ierr);
  for (i=0;i<=ctx->npart;i++) (*subint)[i] = ctx->subintervals[i];
  PetscFunctionReturn(0);
}

/*@C
   EPSKrylovSchurGetSubintervals - Returns the points that delimit the
   subintervals used in spectrum slicing with several partitions.

   Logically Collective on eps

   Input Parameter:
.  eps    - the eigenproblem solver context

   Output Parameter:
.  subint - array of real values specifying subintervals

   Notes:
   If the user passed values with EPSKrylovSchurSetSubintervals(), then the
   same values are returned. Otherwise, the values computed internally are
   obtained.

   This function is only available for spectrum slicing runs.

   The returned array has length npart+1 (see EPSKrylovSchurGetPartitions())
   and should be freed by the user.

   Fortran Notes:
   The calling sequence from Fortran is
.vb
   EPSKrylovSchurGetSubintervals(eps,subint,ierr)
   double precision subint(npart+1) output
.ve

   Level: advanced

.seealso: EPSKrylovSchurSetSubintervals(), EPSKrylovSchurGetPartitions(), EPSSetInterval()
@*/
PetscErrorCode EPSKrylovSchurGetSubintervals(EPS eps,PetscReal **subint)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(subint,2);
  ierr = PetscUseMethod(eps,"EPSKrylovSchurGetSubintervals_C",(EPS,PetscReal**),(eps,subint));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetInertias_KrylovSchur(EPS eps,PetscInt *n,PetscReal **shifts,PetscInt **inertias)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscInt        i,numsh;
  EPS_SR          sr = ctx->sr;

  PetscFunctionBegin;
  if (!eps->state) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
  if (!ctx->sr) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  switch (eps->state) {
  case EPS_STATE_INITIAL:
    break;
  case EPS_STATE_SETUP:
    numsh = ctx->npart+1;
    if (n) *n = numsh;
    if (shifts) {
      ierr = PetscMalloc1(numsh,shifts);CHKERRQ(ierr);
      (*shifts)[0] = eps->inta;
      if (ctx->npart==1) (*shifts)[1] = eps->intb;
      else for (i=1;i<numsh;i++) (*shifts)[i] = ctx->subintervals[i];
    }
    if (inertias) {
      ierr = PetscMalloc1(numsh,inertias);CHKERRQ(ierr);
      (*inertias)[0] = (sr->dir==1)?sr->inertia0:sr->inertia1;
      if (ctx->npart==1) (*inertias)[1] = (sr->dir==1)?sr->inertia1:sr->inertia0;
      else for (i=1;i<numsh;i++) (*inertias)[i] = (*inertias)[i-1]+ctx->nconv_loc[i-1];
    }
    break;
  case EPS_STATE_SOLVED:
  case EPS_STATE_EIGENVECTORS:
    numsh = ctx->nshifts;
    if (n) *n = numsh;
    if (shifts) {
      ierr = PetscMalloc1(numsh,shifts);CHKERRQ(ierr);
      for (i=0;i<numsh;i++) (*shifts)[i] = ctx->shifts[i];
    }
    if (inertias) {
      ierr = PetscMalloc1(numsh,inertias);CHKERRQ(ierr);
      for (i=0;i<numsh;i++) (*inertias)[i] = ctx->inertias[i];
    }
    break;
  }
  PetscFunctionReturn(0);
}

/*@C
   EPSKrylovSchurGetInertias - Gets the values of the shifts and their
   corresponding inertias in case of doing spectrum slicing for a
   computational interval.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  n        - number of shifts, including the endpoints of the interval
.  shifts   - the values of the shifts used internally in the solver
-  inertias - the values of the inertia in each shift

   Notes:
   If called after EPSSolve(), all shifts used internally by the solver are
   returned (including both endpoints and any intermediate ones). If called
   before EPSSolve() and after EPSSetUp() then only the information of the
   endpoints of subintervals is available.

   This function is only available for spectrum slicing runs.

   The returned arrays should be freed by the user. Can pass NULL in any of
   the two arrays if not required.

   Fortran Notes:
   The calling sequence from Fortran is
.vb
   EPSKrylovSchurGetInertias(eps,n,shifts,inertias,ierr)
   integer n
   double precision shifts(*)
   integer inertias(*)
.ve
   The arrays should be at least of length n. The value of n can be determined
   by an initial call
.vb
   EPSKrylovSchurGetInertias(eps,n,PETSC_NULL_REAL,PETSC_NULL_INTEGER,ierr)
.ve

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetSubintervals()
@*/
PetscErrorCode EPSKrylovSchurGetInertias(EPS eps,PetscInt *n,PetscReal **shifts,PetscInt **inertias)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(n,2);
  ierr = PetscUseMethod(eps,"EPSKrylovSchurGetInertias_C",(EPS,PetscInt*,PetscReal**,PetscInt**),(eps,n,shifts,inertias));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetSubcommInfo_KrylovSchur(EPS eps,PetscInt *k,PetscInt *n,Vec *v)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  EPS_SR          sr = ((EPS_KRYLOVSCHUR*)ctx->eps->data)->sr;

  PetscFunctionBegin;
  if (!eps->state) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
  if (!ctx->sr) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  if (k) *k = (ctx->npart==1)? 0: ctx->subc->color;
  if (n) *n = sr->numEigs;
  if (v) {
    ierr = BVCreateVec(sr->V,v);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   EPSKrylovSchurGetSubcommInfo - Gets information related to the case of
   doing spectrum slicing for a computational interval with multiple
   communicators.

   Collective on the subcommunicator (if v is given)

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  k - index of the subinterval for the calling process
.  n - number of eigenvalues found in the k-th subinterval
-  v - a vector owned by processes in the subcommunicator with dimensions
       compatible for locally computed eigenvectors (or NULL)

   Notes:
   This function is only available for spectrum slicing runs.

   The returned Vec should be destroyed by the user.

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubcommPairs()
@*/
PetscErrorCode EPSKrylovSchurGetSubcommInfo(EPS eps,PetscInt *k,PetscInt *n,Vec *v)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = PetscUseMethod(eps,"EPSKrylovSchurGetSubcommInfo_C",(EPS,PetscInt*,PetscInt*,Vec*),(eps,k,n,v));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetSubcommPairs_KrylovSchur(EPS eps,PetscInt i,PetscScalar *eig,Vec v)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  EPS_SR          sr = ((EPS_KRYLOVSCHUR*)ctx->eps->data)->sr;

  PetscFunctionBegin;
  EPSCheckSolved(eps,1);
  if (!ctx->sr) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  if (i<0 || i>=sr->numEigs) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  if (eig) *eig = sr->eigr[sr->perm[i]];
  if (v) { ierr = BVCopyVec(sr->V,sr->perm[i],v);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurGetSubcommPairs - Gets the i-th eigenpair stored
   internally in the subcommunicator to which the calling process belongs.

   Collective on the subcommunicator (if v is given)

   Input Parameter:
+  eps - the eigenproblem solver context
-  i   - index of the solution

   Output Parameters:
+  eig - the eigenvalue
-  v   - the eigenvector

   Notes:
   It is allowed to pass NULL for v if the eigenvector is not required.
   Otherwise, the caller must provide a valid Vec objects, i.e.,
   it must be created by the calling program with EPSKrylovSchurGetSubcommInfo().

   The index i should be a value between 0 and n-1, where n is the number of
   vectors in the local subinterval, see EPSKrylovSchurGetSubcommInfo().

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubcommInfo(), EPSKrylovSchurGetSubcommMats()
@*/
PetscErrorCode EPSKrylovSchurGetSubcommPairs(EPS eps,PetscInt i,PetscScalar *eig,Vec v)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (v) PetscValidLogicalCollectiveInt(v,i,2);
  ierr = PetscUseMethod(eps,"EPSKrylovSchurGetSubcommPairs_C",(EPS,PetscInt,PetscScalar*,Vec),(eps,i,eig,v));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurGetSubcommMats_KrylovSchur(EPS eps,Mat *A,Mat *B)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;

  PetscFunctionBegin;
  if (!ctx->sr) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  if (!eps->state) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
  ierr = EPSGetOperators(ctx->eps,A,B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   EPSKrylovSchurGetSubcommMats - Gets the eigenproblem matrices stored
   internally in the subcommunicator to which the calling process belongs.

   Collective on the subcommunicator

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  A  - the matrix associated with the eigensystem
-  B  - the second matrix in the case of generalized eigenproblems

   Notes:
   This is the analog of EPSGetOperators(), but returns the matrices distributed
   differently (in the subcommunicator rather than in the parent communicator).

   These matrices should not be modified by the user.

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubcommInfo()
@*/
PetscErrorCode EPSKrylovSchurGetSubcommMats(EPS eps,Mat *A,Mat *B)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = PetscTryMethod(eps,"EPSKrylovSchurGetSubcommMats_C",(EPS,Mat*,Mat*),(eps,A,B));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSKrylovSchurUpdateSubcommMats_KrylovSchur(EPS eps,PetscScalar a,PetscScalar ap,Mat Au,PetscScalar b,PetscScalar bp, Mat Bu,MatStructure str,PetscBool globalup)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data,*subctx;
  Mat             A,B=NULL,Ag,Bg=NULL;
  PetscBool       reuse=PETSC_TRUE;

  PetscFunctionBegin;
  if (!ctx->sr) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see EPSSetInterval()");
  if (!eps->state) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must call EPSSetUp() first");
  ierr = EPSGetOperators(eps,&Ag,&Bg);CHKERRQ(ierr);
  ierr = EPSGetOperators(ctx->eps,&A,&B);CHKERRQ(ierr);

  ierr = MatScale(A,a);CHKERRQ(ierr);
  if (Au) {
    ierr = MatAXPY(A,ap,Au,str);CHKERRQ(ierr);
  }
  if (B) ierr = MatScale(B,b);CHKERRQ(ierr);
  if (Bu) {
    ierr = MatAXPY(B,bp,Bu,str);CHKERRQ(ierr);
  }
  ierr = EPSSetOperators(ctx->eps,A,B);CHKERRQ(ierr);

  /* Update stored matrix state */
  subctx = (EPS_KRYLOVSCHUR*)ctx->eps->data;
  ierr = PetscObjectStateGet((PetscObject)A,&subctx->Astate);CHKERRQ(ierr);
  if (B) { ierr = PetscObjectStateGet((PetscObject)B,&subctx->Bstate);CHKERRQ(ierr); }

  /* Update matrices in the parent communicator if requested by user */
  if (globalup) {
    if (ctx->npart>1) {
      if (!ctx->isrow) {
        ierr = MatGetOwnershipIS(Ag,&ctx->isrow,&ctx->iscol);CHKERRQ(ierr);
        reuse = PETSC_FALSE;
      }
      if (str==DIFFERENT_NONZERO_PATTERN) reuse = PETSC_FALSE;
      if (ctx->submata && !reuse) {
        ierr = MatDestroyMatrices(1,&ctx->submata);CHKERRQ(ierr);
      }
      ierr = MatCreateSubMatrices(A,1,&ctx->isrow,&ctx->iscol,(reuse)?MAT_REUSE_MATRIX:MAT_INITIAL_MATRIX,&ctx->submata);CHKERRQ(ierr);
      ierr = MatCreateMPIMatConcatenateSeqMat(((PetscObject)Ag)->comm,ctx->submata[0],PETSC_DECIDE,MAT_REUSE_MATRIX,&Ag);CHKERRQ(ierr);
      if (B) {
        if (ctx->submatb && !reuse) {
          ierr = MatDestroyMatrices(1,&ctx->submatb);CHKERRQ(ierr);
        }
        ierr = MatCreateSubMatrices(B,1,&ctx->isrow,&ctx->iscol,(reuse)?MAT_REUSE_MATRIX:MAT_INITIAL_MATRIX,&ctx->submatb);CHKERRQ(ierr);
        ierr = MatCreateMPIMatConcatenateSeqMat(((PetscObject)Bg)->comm,ctx->submatb[0],PETSC_DECIDE,MAT_REUSE_MATRIX,&Bg);CHKERRQ(ierr);
      }
    }
    ierr = PetscObjectStateGet((PetscObject)Ag,&ctx->Astate);CHKERRQ(ierr);
    if (Bg) { ierr = PetscObjectStateGet((PetscObject)Bg,&ctx->Bstate);CHKERRQ(ierr); }
  }
  ierr = EPSSetOperators(eps,Ag,Bg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSKrylovSchurUpdateSubcommMats - Update the eigenproblem matrices stored
   internally in the subcommunicator to which the calling process belongs.

   Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
.  s   - scalar that multiplies the existing A matrix
.  a   - scalar used in the axpy operation on A
.  Au  - matrix used in the axpy operation on A
.  t   - scalar that multiplies the existing B matrix
.  b   - scalar used in the axpy operation on B
.  Bu  - matrix used in the axpy operation on B
.  str - structure flag
-  globalup - flag indicating if global matrices must be updated

   Notes:
   This function modifies the eigenproblem matrices at the subcommunicator level,
   and optionally updates the global matrices in the parent communicator. The updates
   are expressed as A <-- s*A + a*Au,  B <-- t*B + b*Bu.

   It is possible to update one of the matrices, or both.

   The matrices Au and Bu must be equal in all subcommunicators.

   The str flag is passed to the MatAXPY() operations to perform the updates.

   If globalup is true, communication is carried out to reconstruct the updated
   matrices in the parent communicator. The user must be warned that if global
   matrices are not in sync with subcommunicator matrices, the errors computed
   by EPSComputeError() will be wrong even if the computed solution is correct
   (the synchronization may be done only once at the end).

   Level: advanced

.seealso: EPSSetInterval(), EPSKrylovSchurSetPartitions(), EPSKrylovSchurGetSubcommMats()
@*/
PetscErrorCode EPSKrylovSchurUpdateSubcommMats(EPS eps,PetscScalar s,PetscScalar a,Mat Au,PetscScalar t,PetscScalar b,Mat Bu,MatStructure str,PetscBool globalup)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveScalar(eps,s,2);
  PetscValidLogicalCollectiveScalar(eps,a,3);
  if (Au) PetscValidHeaderSpecific(Au,MAT_CLASSID,4);
  PetscValidLogicalCollectiveScalar(eps,t,5);
  PetscValidLogicalCollectiveScalar(eps,b,6);
  if (Bu) PetscValidHeaderSpecific(Bu,MAT_CLASSID,7);
  PetscValidLogicalCollectiveEnum(eps,str,8);
  PetscValidLogicalCollectiveBool(eps,globalup,9);
  ierr = PetscTryMethod(eps,"EPSKrylovSchurUpdateSubcommMats_C",(EPS,PetscScalar,PetscScalar,Mat,PetscScalar,PetscScalar,Mat,MatStructure,PetscBool),(eps,s,a,Au,t,b,Bu,str,globalup));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_KrylovSchur(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscBool       flg,lock,b,f1,f2,f3;
  PetscReal       keep;
  PetscInt        i,j,k;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS Krylov-Schur Options");CHKERRQ(ierr);

    ierr = PetscOptionsReal("-eps_krylovschur_restart","Proportion of vectors kept after restart","EPSKrylovSchurSetRestart",0.5,&keep,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSKrylovSchurSetRestart(eps,keep);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-eps_krylovschur_locking","Choose between locking and non-locking variants","EPSKrylovSchurSetLocking",PETSC_TRUE,&lock,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSKrylovSchurSetLocking(eps,lock);CHKERRQ(ierr); }

    i = ctx->npart;
    ierr = PetscOptionsInt("-eps_krylovschur_partitions","Number of partitions of the communicator for spectrum slicing","EPSKrylovSchurSetPartitions",ctx->npart,&i,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSKrylovSchurSetPartitions(eps,i);CHKERRQ(ierr); }

    b = ctx->detect;
    ierr = PetscOptionsBool("-eps_krylovschur_detect_zeros","Check zeros during factorizations at subinterval boundaries","EPSKrylovSchurSetDetectZeros",ctx->detect,&b,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSKrylovSchurSetDetectZeros(eps,b);CHKERRQ(ierr); }

    i = 1;
    j = k = PETSC_DECIDE;
    ierr = PetscOptionsInt("-eps_krylovschur_nev","Number of eigenvalues to compute in each subsolve (only for spectrum slicing)","EPSKrylovSchurSetDimensions",40,&i,&f1);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_krylovschur_ncv","Number of basis vectors in each subsolve (only for spectrum slicing)","EPSKrylovSchurSetDimensions",80,&j,&f2);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_krylovschur_mpd","Maximum dimension of projected problem in each subsolve (only for spectrum slicing)","EPSKrylovSchurSetDimensions",80,&k,&f3);CHKERRQ(ierr);
    if (f1 || f2 || f3) { ierr = EPSKrylovSchurSetDimensions(eps,i,j,k);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_KrylovSchur(EPS eps,PetscViewer viewer)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  PetscBool       isascii,isfilt;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  %d%% of basis vectors kept after restart\n",(int)(100*ctx->keep));CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  using the %slocking variant\n",ctx->lock?"":"non-");CHKERRQ(ierr);
    if (eps->which==EPS_ALL) {
      ierr = PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt);CHKERRQ(ierr);
      if (isfilt) {
        ierr = PetscViewerASCIIPrintf(viewer,"  using filtering to extract all eigenvalues in an interval\n");CHKERRQ(ierr);
      } else {
        ierr = PetscViewerASCIIPrintf(viewer,"  doing spectrum slicing with nev=%D, ncv=%D, mpd=%D\n",ctx->nev,ctx->ncv,ctx->mpd);CHKERRQ(ierr);
        if (ctx->npart>1) {
          ierr = PetscViewerASCIIPrintf(viewer,"  multi-communicator spectrum slicing with %D partitions\n",ctx->npart);CHKERRQ(ierr);
          if (ctx->detect) { ierr = PetscViewerASCIIPrintf(viewer,"  detecting zeros when factorizing at subinterval boundaries\n");CHKERRQ(ierr); }
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_KrylovSchur(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetLocking_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetLocking_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetPartitions_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetPartitions_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetDetectZeros_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetDetectZeros_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetDimensions_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetDimensions_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetSubintervals_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubintervals_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetInertias_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommInfo_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommPairs_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommMats_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurUpdateSubcommMats_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_KrylovSchur(EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      isfilt;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilt);CHKERRQ(ierr);
  if (eps->which==EPS_ALL && !isfilt) {
    ierr = EPSReset_KrylovSchur_Slice(eps);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetDefaultST_KrylovSchur(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (eps->which==EPS_ALL) {
    if (!((PetscObject)eps->st)->type_name) {
      ierr = STSetType(eps->st,STSINVERT);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_KrylovSchur(EPS eps)
{
  EPS_KRYLOVSCHUR *ctx;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&ctx);CHKERRQ(ierr);
  eps->data   = (void*)ctx;
  ctx->lock   = PETSC_TRUE;
  ctx->nev    = 1;
  ctx->ncv    = PETSC_DEFAULT;
  ctx->mpd    = PETSC_DEFAULT;
  ctx->npart  = 1;
  ctx->detect = PETSC_FALSE;
  ctx->global = PETSC_TRUE;

  eps->useds = PETSC_TRUE;
  eps->hasts = PETSC_TRUE;

  /* solve and computevectors determined at setup */
  eps->ops->setup          = EPSSetUp_KrylovSchur;
  eps->ops->setfromoptions = EPSSetFromOptions_KrylovSchur;
  eps->ops->destroy        = EPSDestroy_KrylovSchur;
  eps->ops->reset          = EPSReset_KrylovSchur;
  eps->ops->view           = EPSView_KrylovSchur;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_KrylovSchur;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetRestart_C",EPSKrylovSchurSetRestart_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetRestart_C",EPSKrylovSchurGetRestart_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetLocking_C",EPSKrylovSchurSetLocking_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetLocking_C",EPSKrylovSchurGetLocking_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetPartitions_C",EPSKrylovSchurSetPartitions_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetPartitions_C",EPSKrylovSchurGetPartitions_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetDetectZeros_C",EPSKrylovSchurSetDetectZeros_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetDetectZeros_C",EPSKrylovSchurGetDetectZeros_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetDimensions_C",EPSKrylovSchurSetDimensions_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetDimensions_C",EPSKrylovSchurGetDimensions_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurSetSubintervals_C",EPSKrylovSchurSetSubintervals_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubintervals_C",EPSKrylovSchurGetSubintervals_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetInertias_C",EPSKrylovSchurGetInertias_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommInfo_C",EPSKrylovSchurGetSubcommInfo_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommPairs_C",EPSKrylovSchurGetSubcommPairs_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurGetSubcommMats_C",EPSKrylovSchurGetSubcommMats_KrylovSchur);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSKrylovSchurUpdateSubcommMats_C",EPSKrylovSchurUpdateSubcommMats_KrylovSchur);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

