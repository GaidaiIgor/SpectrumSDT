/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc singular value solver: "trlanczos"

   Method: Thick-restart Lanczos

   Algorithm:

       Golub-Kahan-Lanczos bidiagonalization with thick-restart.

   References:

       [1] G.H. Golub and W. Kahan, "Calculating the singular values
           and pseudo-inverse of a matrix", SIAM J. Numer. Anal. Ser.
           B 2:205-224, 1965.

       [2] V. Hernandez, J.E. Roman, and A. Tomas, "A robust and
           efficient parallel SVD solver based on restarted Lanczos
           bidiagonalization", Elec. Trans. Numer. Anal. 31:68-85,
           2008.
*/

#include <slepc/private/svdimpl.h>          /*I "slepcsvd.h" I*/

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-svd,\n"
  "   author = \"V. Hern{\\'a}ndez and J. E. Rom{\\'a}n and A. Tom{\\'a}s\",\n"
  "   title = \"A robust and efficient parallel {SVD} solver based on restarted {Lanczos} bidiagonalization\",\n"
  "   journal = \"Electron. Trans. Numer. Anal.\",\n"
  "   volume = \"31\",\n"
  "   pages = \"68--85\",\n"
  "   year = \"2008\"\n"
  "}\n";

typedef struct {
  PetscBool oneside;
} SVD_TRLANCZOS;

PetscErrorCode SVDSetUp_TRLanczos(SVD svd)
{
  PetscErrorCode ierr;
  PetscInt       N;

  PetscFunctionBegin;
  ierr = SVDMatGetSize(svd,NULL,&N);CHKERRQ(ierr);
  ierr = SVDSetDimensions_Default(svd);CHKERRQ(ierr);
  if (svd->ncv>svd->nsv+svd->mpd) SETERRQ(PetscObjectComm((PetscObject)svd),1,"The value of ncv must not be larger than nev+mpd");
  if (svd->max_it==PETSC_DEFAULT) svd->max_it = PetscMax(N/svd->ncv,100);
  svd->leftbasis = PETSC_TRUE;
  ierr = SVDAllocateSolution(svd,1);CHKERRQ(ierr);
  ierr = DSSetType(svd->ds,DSSVD);CHKERRQ(ierr);
  ierr = DSSetCompact(svd->ds,PETSC_TRUE);CHKERRQ(ierr);
  ierr = DSAllocate(svd->ds,svd->ncv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDOneSideTRLanczosMGS(SVD svd,PetscReal *alpha,PetscReal *beta,BV V,BV U,PetscInt nconv,PetscInt l,PetscInt n,PetscScalar* work)
{
  PetscErrorCode ierr;
  PetscReal      a,b;
  PetscInt       i,k=nconv+l;
  Vec            ui,ui1,vi;

  PetscFunctionBegin;
  ierr = BVGetColumn(V,k,&vi);CHKERRQ(ierr);
  ierr = BVGetColumn(U,k,&ui);CHKERRQ(ierr);
  ierr = SVDMatMult(svd,PETSC_FALSE,vi,ui);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,k,&vi);CHKERRQ(ierr);
  ierr = BVRestoreColumn(U,k,&ui);CHKERRQ(ierr);
  if (l>0) {
    for (i=0;i<l;i++) work[i]=beta[i+nconv];
    ierr = BVMultColumn(U,-1.0,1.0,k,work);CHKERRQ(ierr);
  }
  ierr = BVNormColumn(U,k,NORM_2,&a);CHKERRQ(ierr);
  ierr = BVScaleColumn(U,k,1.0/a);CHKERRQ(ierr);
  alpha[k] = a;

  for (i=k+1;i<n;i++) {
    ierr = BVGetColumn(V,i,&vi);CHKERRQ(ierr);
    ierr = BVGetColumn(U,i-1,&ui1);CHKERRQ(ierr);
    ierr = SVDMatMult(svd,PETSC_TRUE,ui1,vi);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,i,&vi);CHKERRQ(ierr);
    ierr = BVRestoreColumn(U,i-1,&ui1);CHKERRQ(ierr);
    ierr = BVOrthogonalizeColumn(V,i,NULL,&b,NULL);CHKERRQ(ierr);
    ierr = BVScaleColumn(V,i,1.0/b);CHKERRQ(ierr);
    beta[i-1] = b;

    ierr = BVGetColumn(V,i,&vi);CHKERRQ(ierr);
    ierr = BVGetColumn(U,i,&ui);CHKERRQ(ierr);
    ierr = SVDMatMult(svd,PETSC_FALSE,vi,ui);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,i,&vi);CHKERRQ(ierr);
    ierr = BVGetColumn(U,i-1,&ui1);CHKERRQ(ierr);
    ierr = VecAXPY(ui,-b,ui1);CHKERRQ(ierr);
    ierr = BVRestoreColumn(U,i-1,&ui1);CHKERRQ(ierr);
    ierr = BVRestoreColumn(U,i,&ui);CHKERRQ(ierr);
    ierr = BVNormColumn(U,i,NORM_2,&a);CHKERRQ(ierr);
    ierr = BVScaleColumn(U,i,1.0/a);CHKERRQ(ierr);
    alpha[i] = a;
  }

  ierr = BVGetColumn(V,n,&vi);CHKERRQ(ierr);
  ierr = BVGetColumn(U,n-1,&ui1);CHKERRQ(ierr);
  ierr = SVDMatMult(svd,PETSC_TRUE,ui1,vi);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,n,&vi);CHKERRQ(ierr);
  ierr = BVRestoreColumn(U,n-1,&ui1);CHKERRQ(ierr);
  ierr = BVOrthogonalizeColumn(V,n,NULL,&b,NULL);CHKERRQ(ierr);
  beta[n-1] = b;
  PetscFunctionReturn(0);
}

/*
  Custom CGS orthogonalization, preprocess after first orthogonalization
*/
static PetscErrorCode SVDOrthogonalizeCGS(BV V,PetscInt i,PetscScalar* h,PetscReal a,BVOrthogRefineType refine,PetscReal eta,PetscReal *norm)
{
  PetscErrorCode ierr;
  PetscReal      sum,onorm;
  PetscScalar    dot;
  PetscInt       j;

  PetscFunctionBegin;
  switch (refine) {
  case BV_ORTHOG_REFINE_NEVER:
    ierr = BVNormColumn(V,i,NORM_2,norm);CHKERRQ(ierr);
    break;
  case BV_ORTHOG_REFINE_ALWAYS:
    ierr = BVSetActiveColumns(V,0,i);CHKERRQ(ierr);
    ierr = BVDotColumn(V,i,h);CHKERRQ(ierr);
    ierr = BVMultColumn(V,-1.0,1.0,i,h);CHKERRQ(ierr);
    ierr = BVNormColumn(V,i,NORM_2,norm);CHKERRQ(ierr);
    break;
  case BV_ORTHOG_REFINE_IFNEEDED:
    dot = h[i];
    onorm = PetscSqrtReal(PetscRealPart(dot)) / a;
    sum = 0.0;
    for (j=0;j<i;j++) {
      sum += PetscRealPart(h[j] * PetscConj(h[j]));
    }
    *norm = PetscRealPart(dot)/(a*a) - sum;
    if (*norm>0.0) *norm = PetscSqrtReal(*norm);
    else {
      ierr = BVNormColumn(V,i,NORM_2,norm);CHKERRQ(ierr);
    }
    if (*norm < eta*onorm) {
      ierr = BVSetActiveColumns(V,0,i);CHKERRQ(ierr);
      ierr = BVDotColumn(V,i,h);CHKERRQ(ierr);
      ierr = BVMultColumn(V,-1.0,1.0,i,h);CHKERRQ(ierr);
      ierr = BVNormColumn(V,i,NORM_2,norm);CHKERRQ(ierr);
    }
    break;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDOneSideTRLanczosCGS(SVD svd,PetscReal *alpha,PetscReal *beta,BV V,BV U,PetscInt nconv,PetscInt l,PetscInt n,PetscScalar* work)
{
  PetscErrorCode     ierr;
  PetscReal          a,b,eta;
  PetscInt           i,j,k=nconv+l;
  Vec                ui,ui1,vi;
  BVOrthogRefineType refine;

  PetscFunctionBegin;
  ierr = BVGetColumn(V,k,&vi);CHKERRQ(ierr);
  ierr = BVGetColumn(U,k,&ui);CHKERRQ(ierr);
  ierr = SVDMatMult(svd,PETSC_FALSE,vi,ui);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,k,&vi);CHKERRQ(ierr);
  ierr = BVRestoreColumn(U,k,&ui);CHKERRQ(ierr);
  if (l>0) {
    for (i=0;i<l;i++) work[i]=beta[i+nconv];
    ierr = BVMultColumn(U,-1.0,1.0,k,work);CHKERRQ(ierr);
  }
  ierr = BVGetOrthogonalization(V,NULL,&refine,&eta,NULL);CHKERRQ(ierr);

  for (i=k+1;i<n;i++) {
    ierr = BVGetColumn(V,i,&vi);CHKERRQ(ierr);
    ierr = BVGetColumn(U,i-1,&ui1);CHKERRQ(ierr);
    ierr = SVDMatMult(svd,PETSC_TRUE,ui1,vi);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,i,&vi);CHKERRQ(ierr);
    ierr = BVRestoreColumn(U,i-1,&ui1);CHKERRQ(ierr);
    ierr = BVNormColumnBegin(U,i-1,NORM_2,&a);CHKERRQ(ierr);
    if (refine == BV_ORTHOG_REFINE_IFNEEDED) {
      ierr = BVSetActiveColumns(V,0,i+1);CHKERRQ(ierr);
      ierr = BVGetColumn(V,i,&vi);CHKERRQ(ierr);
      ierr = BVDotVecBegin(V,vi,work);CHKERRQ(ierr);
    } else {
      ierr = BVSetActiveColumns(V,0,i);CHKERRQ(ierr);
      ierr = BVDotColumnBegin(V,i,work);CHKERRQ(ierr);
    }
    ierr = BVNormColumnEnd(U,i-1,NORM_2,&a);CHKERRQ(ierr);
    if (refine == BV_ORTHOG_REFINE_IFNEEDED) {
      ierr = BVDotVecEnd(V,vi,work);CHKERRQ(ierr);
      ierr = BVRestoreColumn(V,i,&vi);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(V,0,i);CHKERRQ(ierr);
    } else {
      ierr = BVDotColumnEnd(V,i,work);CHKERRQ(ierr);
    }

    ierr = BVScaleColumn(U,i-1,1.0/a);CHKERRQ(ierr);
    for (j=0;j<i;j++) work[j] = work[j] / a;
    ierr = BVMultColumn(V,-1.0,1.0/a,i,work);CHKERRQ(ierr);
    ierr = SVDOrthogonalizeCGS(V,i,work,a,refine,eta,&b);CHKERRQ(ierr);
    ierr = BVScaleColumn(V,i,1.0/b);CHKERRQ(ierr);
    if (PetscAbsReal(b)<10*PETSC_MACHINE_EPSILON) SETERRQ(PetscObjectComm((PetscObject)svd),1,"Recurrence generated a zero vector; use a two-sided variant");

    ierr = BVGetColumn(V,i,&vi);CHKERRQ(ierr);
    ierr = BVGetColumn(U,i,&ui);CHKERRQ(ierr);
    ierr = BVGetColumn(U,i-1,&ui1);CHKERRQ(ierr);
    ierr = SVDMatMult(svd,PETSC_FALSE,vi,ui);CHKERRQ(ierr);
    ierr = VecAXPY(ui,-b,ui1);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,i,&vi);CHKERRQ(ierr);
    ierr = BVRestoreColumn(U,i,&ui);CHKERRQ(ierr);
    ierr = BVRestoreColumn(U,i-1,&ui1);CHKERRQ(ierr);

    alpha[i-1] = a;
    beta[i-1] = b;
  }

  ierr = BVGetColumn(V,n,&vi);CHKERRQ(ierr);
  ierr = BVGetColumn(U,n-1,&ui1);CHKERRQ(ierr);
  ierr = SVDMatMult(svd,PETSC_TRUE,ui1,vi);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,n,&vi);CHKERRQ(ierr);
  ierr = BVRestoreColumn(U,n-1,&ui1);CHKERRQ(ierr);

  ierr = BVNormColumnBegin(svd->U,n-1,NORM_2,&a);CHKERRQ(ierr);
  if (refine == BV_ORTHOG_REFINE_IFNEEDED) {
    ierr = BVSetActiveColumns(V,0,n+1);CHKERRQ(ierr);
    ierr = BVGetColumn(V,n,&vi);CHKERRQ(ierr);
    ierr = BVDotVecBegin(V,vi,work);CHKERRQ(ierr);
  } else {
    ierr = BVSetActiveColumns(V,0,n);CHKERRQ(ierr);
    ierr = BVDotColumnBegin(V,n,work);CHKERRQ(ierr);
  }
  ierr = BVNormColumnEnd(svd->U,n-1,NORM_2,&a);CHKERRQ(ierr);
  if (refine == BV_ORTHOG_REFINE_IFNEEDED) {
    ierr = BVDotVecEnd(V,vi,work);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,n,&vi);CHKERRQ(ierr);
  } else {
    ierr = BVDotColumnEnd(V,n,work);CHKERRQ(ierr);
  }

  ierr = BVScaleColumn(U,n-1,1.0/a);CHKERRQ(ierr);
  for (j=0;j<n;j++) work[j] = work[j] / a;
  ierr = BVMultColumn(V,-1.0,1.0/a,n,work);CHKERRQ(ierr);
  ierr = SVDOrthogonalizeCGS(V,n,work,a,refine,eta,&b);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(V,nconv,n);CHKERRQ(ierr);
  alpha[n-1] = a;
  beta[n-1] = b;
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSolve_TRLanczos(SVD svd)
{
  PetscErrorCode ierr;
  SVD_TRLANCZOS  *lanczos = (SVD_TRLANCZOS*)svd->data;
  PetscReal      *alpha,*beta,lastbeta,resnorm;
  PetscScalar    *Q,*swork=NULL,*w;
  PetscInt       i,k,l,nv,ld;
  Mat            U,VT;
  PetscBool      conv;
  BVOrthogType   orthog;

  PetscFunctionBegin;
  ierr = PetscCitationsRegister(citation,&cited);CHKERRQ(ierr);
  /* allocate working space */
  ierr = DSGetLeadingDimension(svd->ds,&ld);CHKERRQ(ierr);
  ierr = BVGetOrthogonalization(svd->V,&orthog,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscMalloc1(ld,&w);CHKERRQ(ierr);
  if (lanczos->oneside) {
    ierr = PetscMalloc1(svd->ncv+1,&swork);CHKERRQ(ierr);
  }

  /* normalize start vector */
  if (!svd->nini) {
    ierr = BVSetRandomColumn(svd->V,0);CHKERRQ(ierr);
    ierr = BVOrthonormalizeColumn(svd->V,0,PETSC_TRUE,NULL,NULL);CHKERRQ(ierr);
  }

  l = 0;
  while (svd->reason == SVD_CONVERGED_ITERATING) {
    svd->its++;

    /* inner loop */
    nv = PetscMin(svd->nconv+svd->mpd,svd->ncv);
    ierr = BVSetActiveColumns(svd->V,svd->nconv,nv);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(svd->U,svd->nconv,nv);CHKERRQ(ierr);
    ierr = DSGetArrayReal(svd->ds,DS_MAT_T,&alpha);CHKERRQ(ierr);
    beta = alpha + ld;
    if (lanczos->oneside) {
      if (orthog == BV_ORTHOG_MGS) {
        ierr = SVDOneSideTRLanczosMGS(svd,alpha,beta,svd->V,svd->U,svd->nconv,l,nv,swork);CHKERRQ(ierr);
      } else {
        ierr = SVDOneSideTRLanczosCGS(svd,alpha,beta,svd->V,svd->U,svd->nconv,l,nv,swork);CHKERRQ(ierr);
      }
    } else {
      ierr = SVDTwoSideLanczos(svd,alpha,beta,svd->V,svd->U,svd->nconv+l,nv);CHKERRQ(ierr);
    }
    lastbeta = beta[nv-1];
    ierr = DSRestoreArrayReal(svd->ds,DS_MAT_T,&alpha);CHKERRQ(ierr);
    ierr = BVScaleColumn(svd->V,nv,1.0/lastbeta);CHKERRQ(ierr);

    /* compute SVD of general matrix */
    ierr = DSSetDimensions(svd->ds,nv,nv,svd->nconv,svd->nconv+l);CHKERRQ(ierr);
    if (l==0) {
      ierr = DSSetState(svd->ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
    } else {
      ierr = DSSetState(svd->ds,DS_STATE_RAW);CHKERRQ(ierr);
    }
    ierr = DSSolve(svd->ds,w,NULL);CHKERRQ(ierr);
    ierr = DSSort(svd->ds,w,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DSSynchronize(svd->ds,w,NULL);CHKERRQ(ierr);

    /* compute error estimates */
    k = 0;
    conv = PETSC_TRUE;
    ierr = DSGetArray(svd->ds,DS_MAT_U,&Q);CHKERRQ(ierr);
    ierr = DSGetArrayReal(svd->ds,DS_MAT_T,&alpha);CHKERRQ(ierr);
    beta = alpha + ld;
    for (i=svd->nconv;i<nv;i++) {
      svd->sigma[i] = PetscRealPart(w[i]);
      beta[i] = PetscRealPart(Q[nv-1+i*ld])*lastbeta;
      resnorm = PetscAbsReal(beta[i]);
      ierr = (*svd->converged)(svd,svd->sigma[i],resnorm,&svd->errest[i],svd->convergedctx);CHKERRQ(ierr);
      if (conv) {
        if (svd->errest[i] < svd->tol) k++;
        else conv = PETSC_FALSE;
      }
    }
    ierr = DSRestoreArrayReal(svd->ds,DS_MAT_T,&alpha);CHKERRQ(ierr);
    ierr = DSRestoreArray(svd->ds,DS_MAT_U,&Q);CHKERRQ(ierr);

    /* check convergence and update l */
    ierr = (*svd->stopping)(svd,svd->its,svd->max_it,svd->nconv+k,svd->nsv,&svd->reason,svd->stoppingctx);CHKERRQ(ierr);
    if (svd->reason != SVD_CONVERGED_ITERATING) l = 0;
    else l = PetscMax((nv-svd->nconv-k)/2,0);

    /* compute converged singular vectors and restart vectors */
    ierr = DSGetMat(svd->ds,DS_MAT_VT,&VT);CHKERRQ(ierr);
    ierr = BVMultInPlaceTranspose(svd->V,VT,svd->nconv,svd->nconv+k+l);CHKERRQ(ierr);
    ierr = MatDestroy(&VT);CHKERRQ(ierr);
    ierr = DSGetMat(svd->ds,DS_MAT_U,&U);CHKERRQ(ierr);
    ierr = BVMultInPlace(svd->U,U,svd->nconv,svd->nconv+k+l);CHKERRQ(ierr);
    ierr = MatDestroy(&U);CHKERRQ(ierr);

    /* copy the last vector to be the next initial vector */
    if (svd->reason == SVD_CONVERGED_ITERATING) {
      ierr = BVCopyColumn(svd->V,nv,svd->nconv+k+l);CHKERRQ(ierr);
    }

    svd->nconv += k;
    ierr = SVDMonitor(svd,svd->its,svd->nconv,svd->sigma,svd->errest,nv);CHKERRQ(ierr);
  }

  /* orthonormalize U columns in one side method */
  if (lanczos->oneside) {
    for (i=0;i<svd->nconv;i++) {
      ierr = BVOrthonormalizeColumn(svd->U,i,PETSC_FALSE,NULL,NULL);CHKERRQ(ierr);
    }
  }

  /* free working space */
  ierr = PetscFree(w);CHKERRQ(ierr);
  if (swork) { ierr = PetscFree(swork);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSetFromOptions_TRLanczos(PetscOptionItems *PetscOptionsObject,SVD svd)
{
  PetscErrorCode ierr;
  PetscBool      set,val;
  SVD_TRLANCZOS  *lanczos = (SVD_TRLANCZOS*)svd->data;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"SVD TRLanczos Options");CHKERRQ(ierr);

    ierr = PetscOptionsBool("-svd_trlanczos_oneside","Use one-side reorthogonalization","SVDTRLanczosSetOneSide",lanczos->oneside,&val,&set);CHKERRQ(ierr);
    if (set) { ierr = SVDTRLanczosSetOneSide(svd,val);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDTRLanczosSetOneSide_TRLanczos(SVD svd,PetscBool oneside)
{
  SVD_TRLANCZOS *lanczos = (SVD_TRLANCZOS*)svd->data;

  PetscFunctionBegin;
  if (lanczos->oneside != oneside) {
    lanczos->oneside = oneside;
    svd->state = SVD_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   SVDTRLanczosSetOneSide - Indicate if the variant of the Lanczos method
   to be used is one-sided or two-sided.

   Logically Collective on svd

   Input Parameters:
+  svd     - singular value solver
-  oneside - boolean flag indicating if the method is one-sided or not

   Options Database Key:
.  -svd_trlanczos_oneside <boolean> - Indicates the boolean flag

   Note:
   By default, a two-sided variant is selected, which is sometimes slightly
   more robust. However, the one-sided variant is faster because it avoids
   the orthogonalization associated to left singular vectors.

   Level: advanced

.seealso: SVDLanczosSetOneSide()
@*/
PetscErrorCode SVDTRLanczosSetOneSide(SVD svd,PetscBool oneside)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,oneside,2);
  ierr = PetscTryMethod(svd,"SVDTRLanczosSetOneSide_C",(SVD,PetscBool),(svd,oneside));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDTRLanczosGetOneSide_TRLanczos(SVD svd,PetscBool *oneside)
{
  SVD_TRLANCZOS *lanczos = (SVD_TRLANCZOS*)svd->data;

  PetscFunctionBegin;
  *oneside = lanczos->oneside;
  PetscFunctionReturn(0);
}

/*@
   SVDTRLanczosGetOneSide - Gets if the variant of the Lanczos method
   to be used is one-sided or two-sided.

   Not Collective

   Input Parameters:
.  svd     - singular value solver

   Output Parameters:
.  oneside - boolean flag indicating if the method is one-sided or not

   Level: advanced

.seealso: SVDTRLanczosSetOneSide()
@*/
PetscErrorCode SVDTRLanczosGetOneSide(SVD svd,PetscBool *oneside)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidBoolPointer(oneside,2);
  ierr = PetscUseMethod(svd,"SVDTRLanczosGetOneSide_C",(SVD,PetscBool*),(svd,oneside));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDDestroy_TRLanczos(SVD svd)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(svd->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDTRLanczosSetOneSide_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDTRLanczosGetOneSide_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDView_TRLanczos(SVD svd,PetscViewer viewer)
{
  PetscErrorCode ierr;
  SVD_TRLANCZOS  *lanczos = (SVD_TRLANCZOS*)svd->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  %s-sided reorthogonalization\n",lanczos->oneside? "one": "two");CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_TRLanczos(SVD svd)
{
  PetscErrorCode ierr;
  SVD_TRLANCZOS  *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(svd,&ctx);CHKERRQ(ierr);
  svd->data = (void*)ctx;

  svd->ops->setup          = SVDSetUp_TRLanczos;
  svd->ops->solve          = SVDSolve_TRLanczos;
  svd->ops->destroy        = SVDDestroy_TRLanczos;
  svd->ops->setfromoptions = SVDSetFromOptions_TRLanczos;
  svd->ops->view           = SVDView_TRLanczos;
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDTRLanczosSetOneSide_C",SVDTRLanczosSetOneSide_TRLanczos);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDTRLanczosGetOneSide_C",SVDTRLanczosGetOneSide_TRLanczos);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

