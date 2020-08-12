/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   EPS routines related to the solution process
*/

#include <slepc/private/epsimpl.h>   /*I "slepceps.h" I*/
#include <slepc/private/bvimpl.h>    /*I "slepcbv.h" I*/
#include <petscdraw.h>

PetscErrorCode EPSComputeVectors(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  EPSCheckSolved(eps,1);
  if (eps->state==EPS_STATE_SOLVED && eps->ops->computevectors) {
    ierr = (*eps->ops->computevectors)(eps);CHKERRQ(ierr);
  }
  eps->state = EPS_STATE_EIGENVECTORS;
  PetscFunctionReturn(0);
}

#define SWAP(a,b,t) {t=a;a=b;b=t;}

static PetscErrorCode EPSComputeValues(EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      injective,iscomp,isfilter;
  PetscInt       i,n,aux,nconv0;
  Mat            A,B=NULL,G,Z;

  PetscFunctionBegin;
  switch (eps->categ) {
    case EPS_CATEGORY_KRYLOV:
    case EPS_CATEGORY_OTHER:
      ierr = STIsInjective(eps->st,&injective);CHKERRQ(ierr);
      if (injective) {
        /* one-to-one mapping: backtransform eigenvalues */
        if (eps->ops->backtransform) {
          ierr = (*eps->ops->backtransform)(eps);CHKERRQ(ierr);
        } else SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Internal error, spectral transform should have a backtransform operation");
      } else {
        /* compute eigenvalues from Rayleigh quotient */
        ierr = DSGetDimensions(eps->ds,&n,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
        if (!n) break;
        ierr = EPSGetOperators(eps,&A,&B);CHKERRQ(ierr);
        ierr = BVSetActiveColumns(eps->V,0,n);CHKERRQ(ierr);
        ierr = DSGetCompact(eps->ds,&iscomp);CHKERRQ(ierr);
        ierr = DSSetCompact(eps->ds,PETSC_FALSE);CHKERRQ(ierr);
        ierr = DSGetMat(eps->ds,DS_MAT_A,&G);CHKERRQ(ierr);
        ierr = BVMatProject(eps->V,A,eps->V,G);CHKERRQ(ierr);
        ierr = DSRestoreMat(eps->ds,DS_MAT_A,&G);CHKERRQ(ierr);
        if (B) {
          ierr = DSGetMat(eps->ds,DS_MAT_B,&G);CHKERRQ(ierr);
          ierr = BVMatProject(eps->V,B,eps->V,G);CHKERRQ(ierr);
          ierr = DSRestoreMat(eps->ds,DS_MAT_A,&G);CHKERRQ(ierr);
        }
        ierr = DSSolve(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);
        ierr = DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
        ierr = DSSynchronize(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);
        ierr = DSSetCompact(eps->ds,iscomp);CHKERRQ(ierr);
        if (eps->ishermitian && (!eps->isgeneralized || eps->ispositive)) { /* V = V * Z */
          ierr = DSVectors(eps->ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);
          ierr = DSGetMat(eps->ds,DS_MAT_X,&Z);CHKERRQ(ierr);
          ierr = BVMultInPlace(eps->V,Z,0,n);CHKERRQ(ierr);
          ierr = MatDestroy(&Z);CHKERRQ(ierr);
        }
        /* in case of STFILTER discard computed eigenvalues that lie outside the wanted interval */
        ierr = PetscObjectTypeCompare((PetscObject)eps->st,STFILTER,&isfilter);CHKERRQ(ierr);
        if (isfilter) {
          nconv0 = eps->nconv;
          for (i=0;i<eps->nconv;i++) {
            if (PetscRealPart(eps->eigr[eps->perm[i]])<eps->inta || PetscRealPart(eps->eigr[eps->perm[i]])>eps->intb) {
              eps->nconv--;
              if (i<eps->nconv) { SWAP(eps->perm[i],eps->perm[eps->nconv],aux); i--; }
            }
          }
          if (nconv0>eps->nconv) {
            ierr = PetscInfo1(eps,"Discarded %D computed eigenvalues lying outside the interval\n",nconv0-eps->nconv);CHKERRQ(ierr);
          }
        }
      }
      break;
    case EPS_CATEGORY_PRECOND:
    case EPS_CATEGORY_CONTOUR:
      /* eigenvalues already available as an output of the solver */
      break;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSSolve - Solves the eigensystem.

   Collective on eps

   Input Parameter:
.  eps - eigensolver context obtained from EPSCreate()

   Options Database Keys:
+  -eps_view - print information about the solver used
.  -eps_view_mat0 binary - save the first matrix (A) to the default binary viewer
.  -eps_view_mat1 binary - save the second matrix (B) to the default binary viewer
.  -eps_view_vectors binary - save the computed eigenvectors to the default binary viewer
.  -eps_view_values - print computed eigenvalues
.  -eps_converged_reason - print reason for convergence, and number of iterations
.  -eps_error_absolute - print absolute errors of each eigenpair
.  -eps_error_relative - print relative errors of each eigenpair
-  -eps_error_backward - print backward errors of each eigenpair

   Level: beginner

.seealso: EPSCreate(), EPSSetUp(), EPSDestroy(), EPSSetTolerances()
@*/
PetscErrorCode EPSSolve(EPS eps)
{
  PetscErrorCode ierr;
  PetscInt       i;
  STMatMode      matmode;
  Mat            A,B;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (eps->state>=EPS_STATE_SOLVED) PetscFunctionReturn(0);
  ierr = PetscLogEventBegin(EPS_Solve,eps,0,0,0);CHKERRQ(ierr);

  /* call setup */
  ierr = EPSSetUp(eps);CHKERRQ(ierr);
  eps->nconv = 0;
  eps->its   = 0;
  for (i=0;i<eps->ncv;i++) {
    eps->eigr[i]   = 0.0;
    eps->eigi[i]   = 0.0;
    eps->errest[i] = 0.0;
    eps->perm[i]   = i;
  }
  ierr = EPSViewFromOptions(eps,NULL,"-eps_view_pre");CHKERRQ(ierr);
  ierr = RGViewFromOptions(eps->rg,NULL,"-rg_view");CHKERRQ(ierr);

  /* call solver */
  ierr = (*eps->ops->solve)(eps);CHKERRQ(ierr);
  eps->state = EPS_STATE_SOLVED;

  ierr = STGetMatMode(eps->st,&matmode);CHKERRQ(ierr);
  if (matmode == ST_MATMODE_INPLACE && eps->ispositive) {
    /* Purify eigenvectors before reverting operator */
    ierr = EPSComputeVectors(eps);CHKERRQ(ierr);
  }
  ierr = STPostSolve(eps->st);CHKERRQ(ierr);

  if (!eps->reason) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_PLIB,"Internal error, solver returned without setting converged reason");

  /* Map eigenvalues back to the original problem if appropriate */
  ierr = EPSComputeValues(eps);CHKERRQ(ierr);

#if !defined(PETSC_USE_COMPLEX)
  /* reorder conjugate eigenvalues (positive imaginary first) */
  for (i=0;i<eps->nconv-1;i++) {
    if (eps->eigi[i] != 0) {
      if (eps->eigi[i] < 0) {
        eps->eigi[i] = -eps->eigi[i];
        eps->eigi[i+1] = -eps->eigi[i+1];
        /* the next correction only works with eigenvectors */
        ierr = EPSComputeVectors(eps);CHKERRQ(ierr);
        ierr = BVScaleColumn(eps->V,i+1,-1.0);CHKERRQ(ierr);
      }
      i++;
    }
  }
#endif

  /* sort eigenvalues according to eps->which parameter */
  ierr = SlepcSortEigenvalues(eps->sc,eps->nconv,eps->eigr,eps->eigi,eps->perm);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EPS_Solve,eps,0,0,0);CHKERRQ(ierr);

  /* various viewers */
  ierr = EPSViewFromOptions(eps,NULL,"-eps_view");CHKERRQ(ierr);
  ierr = EPSReasonViewFromOptions(eps);CHKERRQ(ierr);
  ierr = EPSErrorViewFromOptions(eps);CHKERRQ(ierr);
  ierr = EPSValuesViewFromOptions(eps);CHKERRQ(ierr);
  ierr = EPSVectorsViewFromOptions(eps);CHKERRQ(ierr);
  ierr = EPSGetOperators(eps,&A,&B);CHKERRQ(ierr);
  ierr = MatViewFromOptions(A,(PetscObject)eps,"-eps_view_mat0");CHKERRQ(ierr);
  if (eps->isgeneralized) {
    ierr = MatViewFromOptions(B,(PetscObject)eps,"-eps_view_mat1");CHKERRQ(ierr);
  }

  /* Remove deflation and initial subspaces */
  if (eps->nds) {
    ierr = BVSetNumConstraints(eps->V,0);CHKERRQ(ierr);
    eps->nds = 0;
  }
  eps->nini = 0;
  PetscFunctionReturn(0);
}

/*@
   EPSGetIterationNumber - Gets the current iteration number. If the
   call to EPSSolve() is complete, then it returns the number of iterations
   carried out by the solution method.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  its - number of iterations

   Note:
   During the i-th iteration this call returns i-1. If EPSSolve() is
   complete, then parameter "its" contains either the iteration number at
   which convergence was successfully reached, or failure was detected.
   Call EPSGetConvergedReason() to determine if the solver converged or
   failed and why.

   Level: intermediate

.seealso: EPSGetConvergedReason(), EPSSetTolerances()
@*/
PetscErrorCode EPSGetIterationNumber(EPS eps,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(its,2);
  *its = eps->its;
  PetscFunctionReturn(0);
}

/*@
   EPSGetConverged - Gets the number of converged eigenpairs.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  nconv - number of converged eigenpairs

   Note:
   This function should be called after EPSSolve() has finished.

   Level: beginner

.seealso: EPSSetDimensions(), EPSSolve()
@*/
PetscErrorCode EPSGetConverged(EPS eps,PetscInt *nconv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(nconv,2);
  EPSCheckSolved(eps,1);
  *nconv = eps->nconv;
  PetscFunctionReturn(0);
}

/*@
   EPSGetConvergedReason - Gets the reason why the EPSSolve() iteration was
   stopped.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  reason - negative value indicates diverged, positive value converged

   Notes:
   Possible values for reason are
+  EPS_CONVERGED_TOL - converged up to tolerance
.  EPS_CONVERGED_USER - converged due to a user-defined condition
.  EPS_DIVERGED_ITS - required more than max_it iterations to reach convergence
.  EPS_DIVERGED_BREAKDOWN - generic breakdown in method
-  EPS_DIVERGED_SYMMETRY_LOST - pseudo-Lanczos was not able to keep symmetry

   Can only be called after the call to EPSSolve() is complete.

   Level: intermediate

.seealso: EPSSetTolerances(), EPSSolve(), EPSConvergedReason
@*/
PetscErrorCode EPSGetConvergedReason(EPS eps,EPSConvergedReason *reason)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(reason,2);
  EPSCheckSolved(eps,1);
  *reason = eps->reason;
  PetscFunctionReturn(0);
}

/*@
   EPSGetInvariantSubspace - Gets an orthonormal basis of the computed invariant
   subspace.

   Not Collective, but vectors are shared by all processors that share the EPS

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  v - an array of vectors

   Notes:
   This function should be called after EPSSolve() has finished.

   The user should provide in v an array of nconv vectors, where nconv is
   the value returned by EPSGetConverged().

   The first k vectors returned in v span an invariant subspace associated
   with the first k computed eigenvalues (note that this is not true if the
   k-th eigenvalue is complex and matrix A is real; in this case the first
   k+1 vectors should be used). An invariant subspace X of A satisfies Ax
   in X for all x in X (a similar definition applies for generalized
   eigenproblems).

   Level: intermediate

.seealso: EPSGetEigenpair(), EPSGetConverged(), EPSSolve()
@*/
PetscErrorCode EPSGetInvariantSubspace(EPS eps,Vec v[])
{
  PetscErrorCode ierr;
  PetscInt       i;
  BV             V=eps->V;
  Vec            w;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(v,2);
  PetscValidHeaderSpecific(*v,VEC_CLASSID,2);
  EPSCheckSolved(eps,1);
  if (!eps->ishermitian && eps->state==EPS_STATE_EIGENVECTORS) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"EPSGetInvariantSubspace must be called before EPSGetEigenpair,EPSGetEigenvector or EPSComputeError");
  if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
    ierr = BVDuplicateResize(eps->V,eps->nconv,&V);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(eps->V,0,eps->nconv);CHKERRQ(ierr);
    ierr = BVCopy(eps->V,V);CHKERRQ(ierr);
    for (i=0;i<eps->nconv;i++) {
      ierr = BVGetColumn(V,i,&w);CHKERRQ(ierr);
      ierr = VecPointwiseDivide(w,w,eps->D);CHKERRQ(ierr);
      ierr = BVRestoreColumn(V,i,&w);CHKERRQ(ierr);
    }
    ierr = BVOrthogonalize(V,NULL);CHKERRQ(ierr);
  }
  for (i=0;i<eps->nconv;i++) {
    ierr = BVCopyVec(V,i,v[i]);CHKERRQ(ierr);
  }
  if (eps->balance!=EPS_BALANCE_NONE && eps->D) {
    ierr = BVDestroy(&V);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   EPSGetEigenpair - Gets the i-th solution of the eigenproblem as computed by
   EPSSolve(). The solution consists in both the eigenvalue and the eigenvector.

   Logically Collective on eps

   Input Parameters:
+  eps - eigensolver context
-  i   - index of the solution

   Output Parameters:
+  eigr - real part of eigenvalue
.  eigi - imaginary part of eigenvalue
.  Vr   - real part of eigenvector
-  Vi   - imaginary part of eigenvector

   Notes:
   It is allowed to pass NULL for Vr and Vi, if the eigenvector is not
   required. Otherwise, the caller must provide valid Vec objects, i.e.,
   they must be created by the calling program with e.g. MatCreateVecs().

   If the eigenvalue is real, then eigi and Vi are set to zero. If PETSc is
   configured with complex scalars the eigenvalue is stored
   directly in eigr (eigi is set to zero) and the eigenvector in Vr (Vi is
   set to zero). In both cases, the user can pass NULL in eigi and Vi.

   The index i should be a value between 0 and nconv-1 (see EPSGetConverged()).
   Eigenpairs are indexed according to the ordering criterion established
   with EPSSetWhichEigenpairs().

   The 2-norm of the eigenvector is one unless the problem is generalized
   Hermitian. In this case the eigenvector is normalized with respect to the
   norm defined by the B matrix.

   Level: beginner

.seealso: EPSGetEigenvalue(), EPSGetEigenvector(), EPSGetLeftEigenvector(), EPSSolve(),
          EPSGetConverged(), EPSSetWhichEigenpairs(), EPSGetInvariantSubspace()
@*/
PetscErrorCode EPSGetEigenpair(EPS eps,PetscInt i,PetscScalar *eigr,PetscScalar *eigi,Vec Vr,Vec Vi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);
  EPSCheckSolved(eps,1);
  if (i<0 || i>=eps->nconv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  ierr = EPSGetEigenvalue(eps,i,eigr,eigi);CHKERRQ(ierr);
  if (Vr || Vi) { ierr = EPSGetEigenvector(eps,i,Vr,Vi);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*@C
   EPSGetEigenvalue - Gets the i-th eigenvalue as computed by EPSSolve().

   Not Collective

   Input Parameters:
+  eps - eigensolver context
-  i   - index of the solution

   Output Parameters:
+  eigr - real part of eigenvalue
-  eigi - imaginary part of eigenvalue

   Notes:
   If the eigenvalue is real, then eigi is set to zero. If PETSc is
   configured with complex scalars the eigenvalue is stored
   directly in eigr (eigi is set to zero).

   The index i should be a value between 0 and nconv-1 (see EPSGetConverged()).
   Eigenpairs are indexed according to the ordering criterion established
   with EPSSetWhichEigenpairs().

   Level: beginner

.seealso: EPSSolve(), EPSGetConverged(), EPSSetWhichEigenpairs(), EPSGetEigenpair()
@*/
PetscErrorCode EPSGetEigenvalue(EPS eps,PetscInt i,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscInt k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  EPSCheckSolved(eps,1);
  if (i<0 || i>=eps->nconv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  k = eps->perm[i];
#if defined(PETSC_USE_COMPLEX)
  if (eigr) *eigr = eps->eigr[k];
  if (eigi) *eigi = 0;
#else
  if (eigr) *eigr = eps->eigr[k];
  if (eigi) *eigi = eps->eigi[k];
#endif
  PetscFunctionReturn(0);
}

/*@
   EPSGetEigenvector - Gets the i-th right eigenvector as computed by EPSSolve().

   Logically Collective on eps

   Input Parameters:
+  eps - eigensolver context
-  i   - index of the solution

   Output Parameters:
+  Vr   - real part of eigenvector
-  Vi   - imaginary part of eigenvector

   Notes:
   The caller must provide valid Vec objects, i.e., they must be created
   by the calling program with e.g. MatCreateVecs().

   If the corresponding eigenvalue is real, then Vi is set to zero. If PETSc is
   configured with complex scalars the eigenvector is stored
   directly in Vr (Vi is set to zero). In any case, the user can pass NULL in Vr
   or Vi if one of them is not required.

   The index i should be a value between 0 and nconv-1 (see EPSGetConverged()).
   Eigenpairs are indexed according to the ordering criterion established
   with EPSSetWhichEigenpairs().

   The 2-norm of the eigenvector is one unless the problem is generalized
   Hermitian. In this case the eigenvector is normalized with respect to the
   norm defined by the B matrix.

   Level: beginner

.seealso: EPSSolve(), EPSGetConverged(), EPSSetWhichEigenpairs(), EPSGetEigenpair(), EPSGetLeftEigenvector()
@*/
PetscErrorCode EPSGetEigenvector(EPS eps,PetscInt i,Vec Vr,Vec Vi)
{
  PetscErrorCode ierr;
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);
  if (Vr) { PetscValidHeaderSpecific(Vr,VEC_CLASSID,3); PetscCheckSameComm(eps,1,Vr,3); }
  if (Vi) { PetscValidHeaderSpecific(Vi,VEC_CLASSID,4); PetscCheckSameComm(eps,1,Vi,4); }
  EPSCheckSolved(eps,1);
  if (i<0 || i>=eps->nconv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  ierr = EPSComputeVectors(eps);CHKERRQ(ierr);
  k = eps->perm[i];
  ierr = BV_GetEigenvector(eps->V,k,eps->eigi[k],Vr,Vi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGetLeftEigenvector - Gets the i-th left eigenvector as computed by EPSSolve().

   Logically Collective on eps

   Input Parameters:
+  eps - eigensolver context
-  i   - index of the solution

   Output Parameters:
+  Wr   - real part of left eigenvector
-  Wi   - imaginary part of left eigenvector

   Notes:
   The caller must provide valid Vec objects, i.e., they must be created
   by the calling program with e.g. MatCreateVecs().

   If the corresponding eigenvalue is real, then Wi is set to zero. If PETSc is
   configured with complex scalars the eigenvector is stored directly in Wr
   (Wi is set to zero). In any case, the user can pass NULL in Wr or Wi if
   one of them is not required.

   The index i should be a value between 0 and nconv-1 (see EPSGetConverged()).
   Eigensolutions are indexed according to the ordering criterion established
   with EPSSetWhichEigenpairs().

   Left eigenvectors are available only if the twosided flag was set, see
   EPSSetTwoSided().

   Level: intermediate

.seealso: EPSGetEigenvector(), EPSGetConverged(), EPSSetWhichEigenpairs(), EPSSetTwoSided()
@*/
PetscErrorCode EPSGetLeftEigenvector(EPS eps,PetscInt i,Vec Wr,Vec Wi)
{
  PetscErrorCode ierr;
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);
  if (Wr) { PetscValidHeaderSpecific(Wr,VEC_CLASSID,3); PetscCheckSameComm(eps,1,Wr,3); }
  if (Wi) { PetscValidHeaderSpecific(Wi,VEC_CLASSID,4); PetscCheckSameComm(eps,1,Wi,4); }
  EPSCheckSolved(eps,1);
  if (!eps->twosided) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"Must request left vectors with EPSSetTwoSided");
  if (i<0 || i>=eps->nconv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  ierr = EPSComputeVectors(eps);CHKERRQ(ierr);
  k = eps->perm[i];
  ierr = BV_GetEigenvector(eps->W,k,eps->eigi[k],Wr,Wi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGetErrorEstimate - Returns the error estimate associated to the i-th
   computed eigenpair.

   Not Collective

   Input Parameter:
+  eps - eigensolver context
-  i   - index of eigenpair

   Output Parameter:
.  errest - the error estimate

   Notes:
   This is the error estimate used internally by the eigensolver. The actual
   error bound can be computed with EPSComputeError(). See also the users
   manual for details.

   Level: advanced

.seealso: EPSComputeError()
@*/
PetscErrorCode EPSGetErrorEstimate(EPS eps,PetscInt i,PetscReal *errest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidRealPointer(errest,3);
  EPSCheckSolved(eps,1);
  if (i<0 || i>=eps->nconv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  *errest = eps->errest[eps->perm[i]];
  PetscFunctionReturn(0);
}

/*
   EPSComputeResidualNorm_Private - Computes the norm of the residual vector
   associated with an eigenpair.

   Input Parameters:
     trans - whether A' must be used instead of A
     kr,ki - eigenvalue
     xr,xi - eigenvector
     z     - three work vectors (the second one not referenced in complex scalars)
*/
PetscErrorCode EPSComputeResidualNorm_Private(EPS eps,PetscBool trans,PetscScalar kr,PetscScalar ki,Vec xr,Vec xi,Vec *z,PetscReal *norm)
{
  PetscErrorCode ierr;
  PetscInt       nmat;
  Mat            A,B;
  Vec            u,w;
  PetscScalar    alpha;
#if !defined(PETSC_USE_COMPLEX)
  Vec            v;
  PetscReal      ni,nr;
#endif
  PetscErrorCode (*matmult)(Mat,Vec,Vec) = trans? MatMultHermitianTranspose: MatMult;

  PetscFunctionBegin;
  u = z[0]; w = z[2];
  ierr = STGetNumMatrices(eps->st,&nmat);CHKERRQ(ierr);
  ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
  if (nmat>1) { ierr = STGetMatrix(eps->st,1,&B);CHKERRQ(ierr); }

#if !defined(PETSC_USE_COMPLEX)
  v = z[1];
  if (ki == 0 || PetscAbsScalar(ki) < PetscAbsScalar(kr*PETSC_MACHINE_EPSILON)) {
#endif
    ierr = (*matmult)(A,xr,u);CHKERRQ(ierr);                          /* u=A*x */
    if (PetscAbsScalar(kr) > PETSC_MACHINE_EPSILON) {
      if (nmat>1) { ierr = (*matmult)(B,xr,w);CHKERRQ(ierr); }
      else { ierr = VecCopy(xr,w);CHKERRQ(ierr); }                    /* w=B*x */
      alpha = trans? -PetscConj(kr): -kr;
      ierr = VecAXPY(u,alpha,w);CHKERRQ(ierr);                        /* u=A*x-k*B*x */
    }
    ierr = VecNorm(u,NORM_2,norm);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  } else {
    ierr = (*matmult)(A,xr,u);CHKERRQ(ierr);                          /* u=A*xr */
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      if (nmat>1) { ierr = (*matmult)(B,xr,v);CHKERRQ(ierr); }
      else { ierr = VecCopy(xr,v);CHKERRQ(ierr); }                    /* v=B*xr */
      ierr = VecAXPY(u,-kr,v);CHKERRQ(ierr);                          /* u=A*xr-kr*B*xr */
      if (nmat>1) { ierr = (*matmult)(B,xi,w);CHKERRQ(ierr); }
      else { ierr = VecCopy(xi,w);CHKERRQ(ierr); }                    /* w=B*xi */
      ierr = VecAXPY(u,trans?-ki:ki,w);CHKERRQ(ierr);                 /* u=A*xr-kr*B*xr+ki*B*xi */
    }
    ierr = VecNorm(u,NORM_2,&nr);CHKERRQ(ierr);
    ierr = (*matmult)(A,xi,u);CHKERRQ(ierr);                          /* u=A*xi */
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      ierr = VecAXPY(u,-kr,w);CHKERRQ(ierr);                          /* u=A*xi-kr*B*xi */
      ierr = VecAXPY(u,trans?ki:-ki,v);CHKERRQ(ierr);                 /* u=A*xi-kr*B*xi-ki*B*xr */
    }
    ierr = VecNorm(u,NORM_2,&ni);CHKERRQ(ierr);
    *norm = SlepcAbsEigenvalue(nr,ni);
  }
#endif
  PetscFunctionReturn(0);
}

/*@
   EPSComputeError - Computes the error (based on the residual norm) associated
   with the i-th computed eigenpair.

   Collective on eps

   Input Parameter:
+  eps  - the eigensolver context
.  i    - the solution index
-  type - the type of error to compute

   Output Parameter:
.  error - the error

   Notes:
   The error can be computed in various ways, all of them based on the residual
   norm ||Ax-kBx||_2 where k is the eigenvalue and x is the eigenvector.

   Level: beginner

.seealso: EPSErrorType, EPSSolve(), EPSGetErrorEstimate()
@*/
PetscErrorCode EPSComputeError(EPS eps,PetscInt i,EPSErrorType type,PetscReal *error)
{
  PetscErrorCode ierr;
  Mat            A,B;
  Vec            xr,xi,w[3];
  PetscReal      t,vecnorm=1.0,errorl;
  PetscScalar    kr,ki;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);
  PetscValidLogicalCollectiveEnum(eps,type,3);
  PetscValidRealPointer(error,4);
  EPSCheckSolved(eps,1);

  /* allocate work vectors */
#if defined(PETSC_USE_COMPLEX)
  ierr = EPSSetWorkVecs(eps,3);CHKERRQ(ierr);
  xi   = NULL;
  w[1] = NULL;
#else
  ierr = EPSSetWorkVecs(eps,5);CHKERRQ(ierr);
  xi   = eps->work[3];
  w[1] = eps->work[4];
#endif
  xr   = eps->work[0];
  w[0] = eps->work[1];
  w[2] = eps->work[2];

  /* compute residual norm */
  ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
  ierr = EPSComputeResidualNorm_Private(eps,PETSC_FALSE,kr,ki,xr,xi,w,error);CHKERRQ(ierr);

  /* compute 2-norm of eigenvector */
  if (eps->problem_type==EPS_GHEP) {
    ierr = VecNorm(xr,NORM_2,&vecnorm);CHKERRQ(ierr);
  }

  /* if two-sided, compute left residual norm and take the maximum */
  if (eps->twosided) {
    ierr = EPSGetLeftEigenvector(eps,i,xr,xi);CHKERRQ(ierr);
    ierr = EPSComputeResidualNorm_Private(eps,PETSC_TRUE,kr,ki,xr,xi,w,&errorl);CHKERRQ(ierr);
    *error = PetscMax(*error,errorl);
  }

  /* compute error */
  switch (type) {
    case EPS_ERROR_ABSOLUTE:
      break;
    case EPS_ERROR_RELATIVE:
      *error /= SlepcAbsEigenvalue(kr,ki)*vecnorm;
      break;
    case EPS_ERROR_BACKWARD:
      /* initialization of matrix norms */
      if (!eps->nrma) {
        ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
        ierr = MatHasOperation(A,MATOP_NORM,&flg);CHKERRQ(ierr);
        if (!flg) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"The computation of backward errors requires a matrix norm operation");
        ierr = MatNorm(A,NORM_INFINITY,&eps->nrma);CHKERRQ(ierr);
      }
      if (eps->isgeneralized) {
        if (!eps->nrmb) {
          ierr = STGetMatrix(eps->st,1,&B);CHKERRQ(ierr);
          ierr = MatHasOperation(B,MATOP_NORM,&flg);CHKERRQ(ierr);
          if (!flg) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"The computation of backward errors requires a matrix norm operation");
          ierr = MatNorm(B,NORM_INFINITY,&eps->nrmb);CHKERRQ(ierr);
        }
      } else eps->nrmb = 1.0;
      t = SlepcAbsEigenvalue(kr,ki);
      *error /= (eps->nrma+t*eps->nrmb)*vecnorm;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid error type");
  }
  PetscFunctionReturn(0);
}

/*
   EPSGetStartVector - Generate a suitable vector to be used as the starting vector
   for the recurrence that builds the right subspace.

   Collective on eps

   Input Parameters:
+  eps - the eigensolver context
-  i   - iteration number

   Output Parameters:
.  breakdown - flag indicating that a breakdown has occurred

   Notes:
   The start vector is computed from another vector: for the first step (i=0),
   the first initial vector is used (see EPSSetInitialSpace()); otherwise a random
   vector is created. Then this vector is forced to be in the range of OP (only
   for generalized definite problems) and orthonormalized with respect to all
   V-vectors up to i-1. The resulting vector is placed in V[i].

   The flag breakdown is set to true if either i=0 and the vector belongs to the
   deflation space, or i>0 and the vector is linearly dependent with respect
   to the V-vectors.
*/
PetscErrorCode EPSGetStartVector(EPS eps,PetscInt i,PetscBool *breakdown)
{
  PetscErrorCode ierr;
  PetscReal      norm;
  PetscBool      lindep;
  Vec            w,z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);

  /* For the first step, use the first initial vector, otherwise a random one */
  if (i>0 || eps->nini==0) {
    ierr = BVSetRandomColumn(eps->V,i);CHKERRQ(ierr);
  }

  /* Force the vector to be in the range of OP for definite generalized problems */
  if (eps->ispositive || (eps->isgeneralized && eps->ishermitian)) {
    ierr = BVCreateVec(eps->V,&w);CHKERRQ(ierr);
    ierr = BVCopyVec(eps->V,i,w);CHKERRQ(ierr);
    ierr = BVGetColumn(eps->V,i,&z);CHKERRQ(ierr);
    ierr = STApply(eps->st,w,z);CHKERRQ(ierr);
    ierr = BVRestoreColumn(eps->V,i,&z);CHKERRQ(ierr);
    ierr = VecDestroy(&w);CHKERRQ(ierr);
  }

  /* Orthonormalize the vector with respect to previous vectors */
  ierr = BVOrthogonalizeColumn(eps->V,i,NULL,&norm,&lindep);CHKERRQ(ierr);
  if (breakdown) *breakdown = lindep;
  else if (lindep || norm == 0.0) {
    if (i==0) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Initial vector is zero or belongs to the deflation space");
    else SETERRQ(PetscObjectComm((PetscObject)eps),1,"Unable to generate more start vectors");
  }
  ierr = BVScaleColumn(eps->V,i,1.0/norm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   EPSGetLeftStartVector - Generate a suitable vector to be used as the left starting
   vector for the recurrence that builds the left subspace. See EPSGetStartVector().
*/
PetscErrorCode EPSGetLeftStartVector(EPS eps,PetscInt i,PetscBool *breakdown)
{
  PetscErrorCode ierr;
  PetscReal      norm;
  PetscBool      lindep;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,i,2);

  /* For the first step, use the first initial vector, otherwise a random one */
  if (i>0 || eps->ninil==0) {
    ierr = BVSetRandomColumn(eps->W,i);CHKERRQ(ierr);
  }

  /* Orthonormalize the vector with respect to previous vectors */
  ierr = BVOrthogonalizeColumn(eps->W,i,NULL,&norm,&lindep);CHKERRQ(ierr);
  if (breakdown) *breakdown = lindep;
  else if (lindep || norm == 0.0) {
    if (i==0) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Left initial vector is zero");
    else SETERRQ(PetscObjectComm((PetscObject)eps),1,"Unable to generate more left start vectors");
  }
  ierr = BVScaleColumn(eps->W,i,1.0/norm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

