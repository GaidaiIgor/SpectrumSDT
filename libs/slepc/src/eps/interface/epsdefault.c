/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file contains some simple default routines for common operations
*/

#include <slepc/private/epsimpl.h>   /*I "slepceps.h" I*/
#include <slepcvec.h>

PetscErrorCode EPSBackTransform_Default(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = STBackTransform(eps->st,eps->nconv,eps->eigr,eps->eigi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  EPSComputeVectors_Hermitian - Copies the Lanczos vectors as eigenvectors
  using purification for generalized eigenproblems.
 */
PetscErrorCode EPSComputeVectors_Hermitian(EPS eps)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    dot;
  PetscReal      norm;
  PetscBool      iscayley;
  Mat            B;
  Vec            w,x;

  PetscFunctionBegin;
  if (eps->purify) {
    ierr = EPS_Purify(eps,eps->nconv);CHKERRQ(ierr);
    for (i=0;i<eps->nconv;i++) {
      ierr = BVNormColumn(eps->V,i,NORM_2,&norm);CHKERRQ(ierr);
      ierr = BVScaleColumn(eps->V,i,1.0/norm);CHKERRQ(ierr);
    }
  } else {
    /* In the case of Cayley transform, eigenvectors need to be B-normalized */
    ierr = PetscObjectTypeCompare((PetscObject)eps->st,STCAYLEY,&iscayley);CHKERRQ(ierr);
    if (iscayley && eps->isgeneralized) {
      ierr = STGetMatrix(eps->st,1,&B);CHKERRQ(ierr);
      ierr = MatCreateVecs(B,NULL,&w);CHKERRQ(ierr);
      for (i=0;i<eps->nconv;i++) {
        ierr = BVGetColumn(eps->V,i,&x);CHKERRQ(ierr);
        ierr = MatMult(B,x,w);CHKERRQ(ierr);
        ierr = VecDot(w,x,&dot);CHKERRQ(ierr);
        ierr = VecScale(x,1.0/PetscSqrtScalar(dot));CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->V,i,&x);CHKERRQ(ierr);
      }
      ierr = VecDestroy(&w);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*
  EPSComputeVectors_Indefinite - similar to the Schur version but
  for indefinite problems
 */
PetscErrorCode EPSComputeVectors_Indefinite(EPS eps)
{
  PetscErrorCode ierr;
  PetscInt       n,i;
  Mat            X;
  Vec            v;
#if !defined(PETSC_USE_COMPLEX)
  Vec            v1;
#endif

  PetscFunctionBegin;
  ierr = DSGetDimensions(eps->ds,&n,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DSVectors(eps->ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);
  ierr = DSGetMat(eps->ds,DS_MAT_X,&X);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,0,n);CHKERRQ(ierr);
  ierr = BVMultInPlace(eps->V,X,0,n);CHKERRQ(ierr);
  ierr = MatDestroy(&X);CHKERRQ(ierr);

  /* purification */
  if (eps->purify) { ierr = EPS_Purify(eps,eps->nconv);CHKERRQ(ierr); }

  /* normalization */
  for (i=0;i<n;i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (eps->eigi[i] != 0.0) {
      ierr = BVGetColumn(eps->V,i,&v);CHKERRQ(ierr);
      ierr = BVGetColumn(eps->V,i+1,&v1);CHKERRQ(ierr);
      ierr = VecNormalizeComplex(v,v1,PETSC_TRUE,NULL);CHKERRQ(ierr);
      ierr = BVRestoreColumn(eps->V,i,&v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(eps->V,i+1,&v1);CHKERRQ(ierr);
      i++;
    } else
#endif
    {
      ierr = BVGetColumn(eps->V,i,&v);CHKERRQ(ierr);
      ierr = VecNormalize(v,NULL);CHKERRQ(ierr);
      ierr = BVRestoreColumn(eps->V,i,&v);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*
  EPSComputeVectors_Twosided - Adjust left eigenvectors in generalized problems: y = B^-* y.
 */
PetscErrorCode EPSComputeVectors_Twosided(EPS eps)
{
  PetscErrorCode ierr;
  PetscInt       i;
  Vec            w,y;

  PetscFunctionBegin;
  if (!eps->twosided || !eps->isgeneralized) PetscFunctionReturn(0);
  ierr = EPSSetWorkVecs(eps,1);CHKERRQ(ierr);
  w = eps->work[0];
  for (i=0;i<eps->nconv;i++) {
    ierr = BVCopyVec(eps->W,i,w);CHKERRQ(ierr);
    ierr = VecConjugate(w);CHKERRQ(ierr);
    ierr = BVGetColumn(eps->W,i,&y);CHKERRQ(ierr);
    ierr = STMatSolveTranspose(eps->st,w,y);CHKERRQ(ierr);
    ierr = VecConjugate(y);CHKERRQ(ierr);
    ierr = BVRestoreColumn(eps->W,i,&y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  EPSComputeVectors_Schur - Compute eigenvectors from the vectors
  provided by the eigensolver. This version is intended for solvers
  that provide Schur vectors. Given the partial Schur decomposition
  OP*V=V*T, the following steps are performed:
      1) compute eigenvectors of T: T*Z=Z*D
      2) compute eigenvectors of OP: X=V*Z
 */
PetscErrorCode EPSComputeVectors_Schur(EPS eps)
{
  PetscErrorCode ierr;
  PetscInt       n,i;
  Mat            Z;
  Vec            z,v;
#if !defined(PETSC_USE_COMPLEX)
  Vec            v1;
#endif

  PetscFunctionBegin;
  if (eps->ishermitian) {
    if (eps->isgeneralized && !eps->ispositive) {
      ierr =  EPSComputeVectors_Indefinite(eps);CHKERRQ(ierr);
    } else {
      ierr = EPSComputeVectors_Hermitian(eps);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }
  ierr = DSGetDimensions(eps->ds,&n,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  /* right eigenvectors */
  ierr = DSVectors(eps->ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);

  /* V = V * Z */
  ierr = DSGetMat(eps->ds,DS_MAT_X,&Z);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,0,n);CHKERRQ(ierr);
  ierr = BVMultInPlace(eps->V,Z,0,n);CHKERRQ(ierr);
  ierr = MatDestroy(&Z);CHKERRQ(ierr);

  /* Purify eigenvectors */
  if (eps->purify) { ierr = EPS_Purify(eps,eps->nconv);CHKERRQ(ierr); }

  /* Fix eigenvectors if balancing was used */
  if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
    for (i=0;i<n;i++) {
      ierr = BVGetColumn(eps->V,i,&z);CHKERRQ(ierr);
      ierr = VecPointwiseDivide(z,z,eps->D);CHKERRQ(ierr);
      ierr = BVRestoreColumn(eps->V,i,&z);CHKERRQ(ierr);
    }
  }

  /* normalize eigenvectors (when using purification or balancing) */
  if (eps->purify || (eps->balance!=EPS_BALANCE_NONE && eps->D)) {
    for (i=0;i<n;i++) {
#if !defined(PETSC_USE_COMPLEX)
      if (eps->eigi[i] != 0.0) {
        ierr = BVGetColumn(eps->V,i,&v);CHKERRQ(ierr);
        ierr = BVGetColumn(eps->V,i+1,&v1);CHKERRQ(ierr);
        ierr = VecNormalizeComplex(v,v1,PETSC_TRUE,NULL);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->V,i,&v);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->V,i+1,&v1);CHKERRQ(ierr);
        i++;
      } else
#endif
      {
        ierr = BVGetColumn(eps->V,i,&v);CHKERRQ(ierr);
        ierr = VecNormalize(v,NULL);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->V,i,&v);CHKERRQ(ierr);
      }
    }
  }

  /* left eigenvectors from eps->dsts */
  if (eps->twosided) {
    ierr = DSVectors(eps->dsts,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);
    /* W = W * Z */
    ierr = DSGetMat(eps->dsts,DS_MAT_X,&Z);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(eps->W,0,n);CHKERRQ(ierr);
    ierr = BVMultInPlace(eps->W,Z,0,n);CHKERRQ(ierr);
    ierr = MatDestroy(&Z);CHKERRQ(ierr);
    /* Fix left eigenvectors if balancing was used */
    if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
      for (i=0;i<n;i++) {
        ierr = BVGetColumn(eps->W,i,&z);CHKERRQ(ierr);
        ierr = VecPointwiseMult(z,z,eps->D);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->W,i,&z);CHKERRQ(ierr);
      }
    }
    ierr = EPSComputeVectors_Twosided(eps);CHKERRQ(ierr);
    /* normalize */
    for (i=0;i<eps->nconv;i++) {
#if !defined(PETSC_USE_COMPLEX)
      if (eps->eigi[i]!=0.0) {   /* first eigenvalue of a complex conjugate pair */
        ierr = BVGetColumn(eps->W,i,&v);CHKERRQ(ierr);
        ierr = BVGetColumn(eps->W,i+1,&v1);CHKERRQ(ierr);
        ierr = VecNormalizeComplex(v,v1,PETSC_TRUE,NULL);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->W,i,&v);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->W,i+1,&v1);CHKERRQ(ierr);
        i++;
      } else   /* real eigenvalue */
#endif
      {
        ierr = BVGetColumn(eps->W,i,&v);CHKERRQ(ierr);
        ierr = VecNormalize(v,NULL);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->W,i,&v);CHKERRQ(ierr);
      }
    }
#if !defined(PETSC_USE_COMPLEX)
    for (i=0;i<eps->nconv-1;i++) {
      if (eps->eigi[i] != 0.0) {
        if (eps->eigi[i] > 0.0) { ierr = BVScaleColumn(eps->W,i+1,-1.0);CHKERRQ(ierr); }
        i++;
      }
    }
#endif
  }
  PetscFunctionReturn(0);
}

/*@
   EPSSetWorkVecs - Sets a number of work vectors into an EPS object.

   Collective on eps

   Input Parameters:
+  eps - eigensolver context
-  nw  - number of work vectors to allocate

   Developers Note:
   This is SLEPC_EXTERN because it may be required by user plugin EPS
   implementations.

   Level: developer
@*/
PetscErrorCode EPSSetWorkVecs(EPS eps,PetscInt nw)
{
  PetscErrorCode ierr;
  Vec            t;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nw,2);
  if (nw <= 0) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"nw must be > 0: nw = %D",nw);
  if (eps->nwork < nw) {
    ierr = VecDestroyVecs(eps->nwork,&eps->work);CHKERRQ(ierr);
    eps->nwork = nw;
    ierr = BVGetColumn(eps->V,0,&t);CHKERRQ(ierr);
    ierr = VecDuplicateVecs(t,nw,&eps->work);CHKERRQ(ierr);
    ierr = BVRestoreColumn(eps->V,0,&t);CHKERRQ(ierr);
    ierr = PetscLogObjectParents(eps,nw,eps->work);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  EPSSetWhichEigenpairs_Default - Sets the default value for which,
  depending on the ST.
 */
PetscErrorCode EPSSetWhichEigenpairs_Default(EPS eps)
{
  PetscBool      target;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompareAny((PetscObject)eps->st,&target,STSINVERT,STCAYLEY,"");CHKERRQ(ierr);
  if (target) eps->which = EPS_TARGET_MAGNITUDE;
  else eps->which = EPS_LARGEST_MAGNITUDE;
  PetscFunctionReturn(0);
}

/*
  EPSConvergedRelative - Checks convergence relative to the eigenvalue.
*/
PetscErrorCode EPSConvergedRelative(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal w;

  PetscFunctionBegin;
  w = SlepcAbsEigenvalue(eigr,eigi);
  *errest = res/w;
  PetscFunctionReturn(0);
}

/*
  EPSConvergedAbsolute - Checks convergence absolutely.
*/
PetscErrorCode EPSConvergedAbsolute(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = res;
  PetscFunctionReturn(0);
}

/*
  EPSConvergedNorm - Checks convergence relative to the eigenvalue and
  the matrix norms.
*/
PetscErrorCode EPSConvergedNorm(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal w;

  PetscFunctionBegin;
  w = SlepcAbsEigenvalue(eigr,eigi);
  *errest = res / (eps->nrma + w*eps->nrmb);
  PetscFunctionReturn(0);
}

/*@C
   EPSStoppingBasic - Default routine to determine whether the outer eigensolver
   iteration must be stopped.

   Collective on eps

   Input Parameters:
+  eps    - eigensolver context obtained from EPSCreate()
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged eigenpairs
.  nev    - number of requested eigenpairs
-  ctx    - context (not used here)

   Output Parameter:
.  reason - result of the stopping test

   Notes:
   A positive value of reason indicates that the iteration has finished successfully
   (converged), and a negative value indicates an error condition (diverged). If
   the iteration needs to be continued, reason must be set to EPS_CONVERGED_ITERATING
   (zero).

   EPSStoppingBasic() will stop if all requested eigenvalues are converged, or if
   the maximum number of iterations has been reached.

   Use EPSSetStoppingTest() to provide your own test instead of using this one.

   Level: advanced

.seealso: EPSSetStoppingTest(), EPSConvergedReason, EPSGetConvergedReason()
@*/
PetscErrorCode EPSStoppingBasic(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *reason = EPS_CONVERGED_ITERATING;
  if (nconv >= nev) {
    ierr = PetscInfo2(eps,"Linear eigensolver finished successfully: %D eigenpairs converged at iteration %D\n",nconv,its);CHKERRQ(ierr);
    *reason = EPS_CONVERGED_TOL;
  } else if (its >= max_it) {
    *reason = EPS_DIVERGED_ITS;
    ierr = PetscInfo1(eps,"Linear eigensolver iteration reached maximum number of iterations (%D)\n",its);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  EPSComputeRitzVector - Computes the current Ritz vector.

  Simple case (complex scalars or real scalars with Zi=NULL):
    x = V*Zr  (V is a basis of nv vectors, Zr has length nv)

  Split case:
    x = V*Zr  y = V*Zi  (Zr and Zi have length nv)
*/
PetscErrorCode EPSComputeRitzVector(EPS eps,PetscScalar *Zr,PetscScalar *Zi,BV V,Vec x,Vec y)
{
  PetscErrorCode ierr;
  PetscInt       l,k;
  PetscReal      norm;
#if !defined(PETSC_USE_COMPLEX)
  Vec            z;
#endif

  PetscFunctionBegin;
  /* compute eigenvector */
  ierr = BVGetActiveColumns(V,&l,&k);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(V,0,k);CHKERRQ(ierr);
  ierr = BVMultVec(V,1.0,0.0,x,Zr);CHKERRQ(ierr);

  /* purify eigenvector if necessary */
  if (eps->purify) {
    ierr = STApply(eps->st,x,y);CHKERRQ(ierr);
    if (eps->ishermitian) {
      ierr = BVNormVec(eps->V,y,NORM_2,&norm);CHKERRQ(ierr);
    } else {
      ierr = VecNorm(y,NORM_2,&norm);CHKERRQ(ierr);
    }
    ierr = VecScale(y,1.0/norm);CHKERRQ(ierr);
    ierr = VecCopy(y,x);CHKERRQ(ierr);
  }
  /* fix eigenvector if balancing is used */
  if (!eps->ishermitian && eps->balance!=EPS_BALANCE_NONE && eps->D) {
    ierr = VecPointwiseDivide(x,x,eps->D);CHKERRQ(ierr);
  }
#if !defined(PETSC_USE_COMPLEX)
  /* compute imaginary part of eigenvector */
  if (Zi) {
    ierr = BVMultVec(V,1.0,0.0,y,Zi);CHKERRQ(ierr);
    if (eps->ispositive) {
      ierr = BVCreateVec(V,&z);CHKERRQ(ierr);
      ierr = STApply(eps->st,y,z);CHKERRQ(ierr);
      ierr = VecNorm(z,NORM_2,&norm);CHKERRQ(ierr);
      ierr = VecScale(z,1.0/norm);CHKERRQ(ierr);
      ierr = VecCopy(z,y);CHKERRQ(ierr);
      ierr = VecDestroy(&z);CHKERRQ(ierr);
    }
    if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
      ierr = VecPointwiseDivide(y,y,eps->D);CHKERRQ(ierr);
    }
  } else
#endif
  { ierr = VecSet(y,0.0);CHKERRQ(ierr); }

  /* normalize eigenvectors (when using balancing) */
  if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
#if !defined(PETSC_USE_COMPLEX)
    if (Zi) {
      ierr = VecNormalizeComplex(x,y,PETSC_TRUE,NULL);CHKERRQ(ierr);
    } else
#endif
    {
      ierr = VecNormalize(x,NULL);CHKERRQ(ierr);
    }
  }
  ierr = BVSetActiveColumns(V,l,k);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  EPSBuildBalance_Krylov - uses a Krylov subspace method to compute the
  diagonal matrix to be applied for balancing in non-Hermitian problems.
*/
PetscErrorCode EPSBuildBalance_Krylov(EPS eps)
{
  Vec               z,p,r;
  PetscInt          i,j;
  PetscReal         norma;
  PetscScalar       *pz,*pD;
  const PetscScalar *pr,*pp;
  PetscRandom       rand;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  ierr = EPSSetWorkVecs(eps,3);CHKERRQ(ierr);
  ierr = BVGetRandomContext(eps->V,&rand);CHKERRQ(ierr);
  r = eps->work[0];
  p = eps->work[1];
  z = eps->work[2];
  ierr = VecSet(eps->D,1.0);CHKERRQ(ierr);

  for (j=0;j<eps->balance_its;j++) {

    /* Build a random vector of +-1's */
    ierr = VecSetRandom(z,rand);CHKERRQ(ierr);
    ierr = VecGetArray(z,&pz);CHKERRQ(ierr);
    for (i=0;i<eps->nloc;i++) {
      if (PetscRealPart(pz[i])<0.5) pz[i]=-1.0;
      else pz[i]=1.0;
    }
    ierr = VecRestoreArray(z,&pz);CHKERRQ(ierr);

    /* Compute p=DA(D\z) */
    ierr = VecPointwiseDivide(r,z,eps->D);CHKERRQ(ierr);
    ierr = STApply(eps->st,r,p);CHKERRQ(ierr);
    ierr = VecPointwiseMult(p,p,eps->D);CHKERRQ(ierr);
    if (eps->balance == EPS_BALANCE_TWOSIDE) {
      if (j==0) {
        /* Estimate the matrix inf-norm */
        ierr = VecAbs(p);CHKERRQ(ierr);
        ierr = VecMax(p,NULL,&norma);CHKERRQ(ierr);
      }
      /* Compute r=D\(A'Dz) */
      ierr = VecPointwiseMult(z,z,eps->D);CHKERRQ(ierr);
      ierr = STApplyHermitianTranspose(eps->st,z,r);CHKERRQ(ierr);
      ierr = VecPointwiseDivide(r,r,eps->D);CHKERRQ(ierr);
    }

    /* Adjust values of D */
    ierr = VecGetArrayRead(r,&pr);CHKERRQ(ierr);
    ierr = VecGetArrayRead(p,&pp);CHKERRQ(ierr);
    ierr = VecGetArray(eps->D,&pD);CHKERRQ(ierr);
    for (i=0;i<eps->nloc;i++) {
      if (eps->balance == EPS_BALANCE_TWOSIDE) {
        if (PetscAbsScalar(pp[i])>eps->balance_cutoff*norma && pr[i]!=0.0)
          pD[i] *= PetscSqrtReal(PetscAbsScalar(pr[i]/pp[i]));
      } else {
        if (pp[i]!=0.0) pD[i] /= PetscAbsScalar(pp[i]);
      }
    }
    ierr = VecRestoreArrayRead(r,&pr);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(p,&pp);CHKERRQ(ierr);
    ierr = VecRestoreArray(eps->D,&pD);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

