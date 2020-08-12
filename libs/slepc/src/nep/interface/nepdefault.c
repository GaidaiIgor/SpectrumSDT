/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Simple default routines for common NEP operations
*/

#include <slepc/private/nepimpl.h>     /*I "slepcnep.h" I*/

/*@
   NEPSetWorkVecs - Sets a number of work vectors into a NEP object

   Collective on nep

   Input Parameters:
+  nep - nonlinear eigensolver context
-  nw  - number of work vectors to allocate

   Developers Note:
   This is SLEPC_EXTERN because it may be required by user plugin NEP
   implementations.

   Level: developer
@*/
PetscErrorCode NEPSetWorkVecs(NEP nep,PetscInt nw)
{
  PetscErrorCode ierr;
  Vec            t;

  PetscFunctionBegin;
  if (nep->nwork < nw) {
    ierr = VecDestroyVecs(nep->nwork,&nep->work);CHKERRQ(ierr);
    nep->nwork = nw;
    ierr = BVGetColumn(nep->V,0,&t);CHKERRQ(ierr);
    ierr = VecDuplicateVecs(t,nw,&nep->work);CHKERRQ(ierr);
    ierr = BVRestoreColumn(nep->V,0,&t);CHKERRQ(ierr);
    ierr = PetscLogObjectParents(nep,nw,nep->work);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  NEPGetDefaultShift - Return the value of sigma to start the nonlinear iteration.
 */
PetscErrorCode NEPGetDefaultShift(NEP nep,PetscScalar *sigma)
{
  PetscFunctionBegin;
  PetscValidScalarPointer(sigma,2);
  switch (nep->which) {
    case NEP_LARGEST_MAGNITUDE:
    case NEP_LARGEST_IMAGINARY:
    case NEP_ALL:
    case NEP_WHICH_USER:
      *sigma = 1.0;   /* arbitrary value */
      break;
    case NEP_SMALLEST_MAGNITUDE:
    case NEP_SMALLEST_IMAGINARY:
      *sigma = 0.0;
      break;
    case NEP_LARGEST_REAL:
      *sigma = PETSC_MAX_REAL;
      break;
    case NEP_SMALLEST_REAL:
      *sigma = PETSC_MIN_REAL;
      break;
    case NEP_TARGET_MAGNITUDE:
    case NEP_TARGET_REAL:
    case NEP_TARGET_IMAGINARY:
      *sigma = nep->target;
      break;
  }
  PetscFunctionReturn(0);
}

/*
  NEPConvergedRelative - Checks convergence relative to the eigenvalue.
*/
PetscErrorCode NEPConvergedRelative(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal w;

  PetscFunctionBegin;
  w = SlepcAbsEigenvalue(eigr,eigi);
  *errest = res/w;
  PetscFunctionReturn(0);
}

/*
  NEPConvergedAbsolute - Checks convergence absolutely.
*/
PetscErrorCode NEPConvergedAbsolute(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = res;
  PetscFunctionReturn(0);
}

/*
  NEPConvergedNorm - Checks convergence relative to the matrix norms.
*/
PetscErrorCode NEPConvergedNorm(NEP nep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscScalar    s;
  PetscReal      w=0.0;
  PetscInt       j;
  PetscBool      flg;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (nep->fui!=NEP_USER_INTERFACE_SPLIT) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Backward error only available in split form");
  /* initialization of matrix norms */
  if (!nep->nrma[0]) {
    for (j=0;j<nep->nt;j++) {
      ierr = MatHasOperation(nep->A[j],MATOP_NORM,&flg);CHKERRQ(ierr);
      if (!flg) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"The convergence test related to the matrix norms requires a matrix norm operation");
      ierr = MatNorm(nep->A[j],NORM_INFINITY,&nep->nrma[j]);CHKERRQ(ierr);
    }
  }
  for (j=0;j<nep->nt;j++) {
    ierr = FNEvaluateFunction(nep->f[j],eigr,&s);CHKERRQ(ierr);
    w = w + nep->nrma[j]*PetscAbsScalar(s);
  }
  *errest = res/w;
  PetscFunctionReturn(0);
}

/*@C
   NEPStoppingBasic - Default routine to determine whether the outer eigensolver
   iteration must be stopped.

   Collective on nep

   Input Parameters:
+  nep    - nonlinear eigensolver context obtained from NEPCreate()
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
   the iteration needs to be continued, reason must be set to NEP_CONVERGED_ITERATING
   (zero).

   NEPStoppingBasic() will stop if all requested eigenvalues are converged, or if
   the maximum number of iterations has been reached.

   Use NEPSetStoppingTest() to provide your own test instead of using this one.

   Level: advanced

.seealso: NEPSetStoppingTest(), NEPConvergedReason, NEPGetConvergedReason()
@*/
PetscErrorCode NEPStoppingBasic(NEP nep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,NEPConvergedReason *reason,void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *reason = NEP_CONVERGED_ITERATING;
  if (nconv >= nev) {
    ierr = PetscInfo2(nep,"Nonlinear eigensolver finished successfully: %D eigenpairs converged at iteration %D\n",nconv,its);CHKERRQ(ierr);
    *reason = NEP_CONVERGED_TOL;
  } else if (its >= max_it) {
    *reason = NEP_DIVERGED_ITS;
    ierr = PetscInfo1(nep,"Nonlinear eigensolver iteration reached maximum number of iterations (%D)\n",its);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPComputeVectors_Schur(NEP nep)
{
  PetscErrorCode ierr;
  PetscInt       n,i;
  Mat            Z;
  Vec            v;

  PetscFunctionBegin;
  ierr = DSGetDimensions(nep->ds,&n,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DSVectors(nep->ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);
  ierr = DSGetMat(nep->ds,DS_MAT_X,&Z);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(nep->V,0,n);CHKERRQ(ierr);
  ierr = BVMultInPlace(nep->V,Z,0,n);CHKERRQ(ierr);
  ierr = MatDestroy(&Z);CHKERRQ(ierr);

  /* normalization */
  for (i=0;i<n;i++) {
    ierr = BVGetColumn(nep->V,i,&v);CHKERRQ(ierr);
    ierr = VecNormalize(v,NULL);CHKERRQ(ierr);
    ierr = BVRestoreColumn(nep->V,i,&v);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

