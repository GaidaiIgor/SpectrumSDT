/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SVD routines related to the solution process
*/

#include <slepc/private/svdimpl.h>   /*I "slepcsvd.h" I*/

PetscErrorCode SVDComputeVectors(SVD svd)
{
  PetscErrorCode ierr;
  Vec            tl,uj,vj;
  PetscInt       j,oldsize;

  PetscFunctionBegin;
  SVDCheckSolved(svd,1);
  if (svd->state==SVD_STATE_SOLVED) {
    /* generate left singular vectors on U */
    if (!svd->U) { ierr = SVDGetBV(svd,NULL,&svd->U);CHKERRQ(ierr); }
    ierr = BVGetSizes(svd->U,NULL,NULL,&oldsize);CHKERRQ(ierr);
    if (!oldsize) {
      if (!((PetscObject)(svd->U))->type_name) {
        ierr = BVSetType(svd->U,BVSVEC);CHKERRQ(ierr);
      }
      ierr = SVDMatCreateVecsEmpty(svd,NULL,&tl);CHKERRQ(ierr);
      ierr = BVSetSizesFromVec(svd->U,tl,svd->ncv);CHKERRQ(ierr);
      ierr = VecDestroy(&tl);CHKERRQ(ierr);
    }
    for (j=0;j<svd->nconv;j++) {
      ierr = BVGetColumn(svd->V,j,&vj);CHKERRQ(ierr);
      ierr = BVGetColumn(svd->U,j,&uj);CHKERRQ(ierr);
      ierr = SVDMatMult(svd,PETSC_FALSE,vj,uj);CHKERRQ(ierr);
      ierr = BVRestoreColumn(svd->V,j,&vj);CHKERRQ(ierr);
      ierr = BVRestoreColumn(svd->U,j,&uj);CHKERRQ(ierr);
      ierr = BVOrthonormalizeColumn(svd->U,j,PETSC_FALSE,NULL,NULL);CHKERRQ(ierr);
    }
  }
  svd->state = SVD_STATE_VECTORS;
  PetscFunctionReturn(0);
}

/*@
   SVDSolve - Solves the singular value problem.

   Collective on svd

   Input Parameter:
.  svd - singular value solver context obtained from SVDCreate()

   Options Database Keys:
+  -svd_view - print information about the solver used
.  -svd_view_mat binary - save the matrix to the default binary viewer
.  -svd_view_vectors binary - save the computed singular vectors to the default binary viewer
.  -svd_view_values - print computed singular values
.  -svd_converged_reason - print reason for convergence, and number of iterations
.  -svd_error_absolute - print absolute errors of each singular triplet
-  -svd_error_relative - print relative errors of each singular triplet

   Level: beginner

.seealso: SVDCreate(), SVDSetUp(), SVDDestroy()
@*/
PetscErrorCode SVDSolve(SVD svd)
{
  PetscErrorCode ierr;
  PetscInt       i,*workperm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (svd->state>=SVD_STATE_SOLVED) PetscFunctionReturn(0);
  ierr = PetscLogEventBegin(SVD_Solve,svd,0,0,0);CHKERRQ(ierr);

  /* call setup */
  ierr = SVDSetUp(svd);CHKERRQ(ierr);
  svd->its = 0;
  svd->nconv = 0;
  for (i=0;i<svd->ncv;i++) {
    svd->sigma[i]  = 0.0;
    svd->errest[i] = 0.0;
    svd->perm[i]   = i;
  }
  ierr = SVDViewFromOptions(svd,NULL,"-svd_view_pre");CHKERRQ(ierr);

  ierr = (*svd->ops->solve)(svd);CHKERRQ(ierr);
  svd->state = (svd->leftbasis)? SVD_STATE_VECTORS: SVD_STATE_SOLVED;

  /* sort singular triplets */
  if (svd->which == SVD_SMALLEST) {
    ierr = PetscSortRealWithPermutation(svd->nconv,svd->sigma,svd->perm);CHKERRQ(ierr);
  } else {
    ierr = PetscMalloc1(svd->nconv,&workperm);CHKERRQ(ierr);
    for (i=0;i<svd->nconv;i++) workperm[i] = i;
    ierr = PetscSortRealWithPermutation(svd->nconv,svd->sigma,workperm);CHKERRQ(ierr);
    for (i=0;i<svd->nconv;i++) svd->perm[i] = workperm[svd->nconv-i-1];
    ierr = PetscFree(workperm);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(SVD_Solve,svd,0,0,0);CHKERRQ(ierr);

  /* various viewers */
  ierr = SVDViewFromOptions(svd,NULL,"-svd_view");CHKERRQ(ierr);
  ierr = SVDReasonViewFromOptions(svd);CHKERRQ(ierr);
  ierr = SVDErrorViewFromOptions(svd);CHKERRQ(ierr);
  ierr = SVDValuesViewFromOptions(svd);CHKERRQ(ierr);
  ierr = SVDVectorsViewFromOptions(svd);CHKERRQ(ierr);
  ierr = MatViewFromOptions(svd->OP,(PetscObject)svd,"-svd_view_mat");CHKERRQ(ierr);

  /* Remove the initial subspaces */
  svd->nini = 0;
  svd->ninil = 0;
  PetscFunctionReturn(0);
}

/*@
   SVDGetIterationNumber - Gets the current iteration number. If the
   call to SVDSolve() is complete, then it returns the number of iterations
   carried out by the solution method.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  its - number of iterations

   Note:
   During the i-th iteration this call returns i-1. If SVDSolve() is
   complete, then parameter "its" contains either the iteration number at
   which convergence was successfully reached, or failure was detected.
   Call SVDGetConvergedReason() to determine if the solver converged or
   failed and why.

   Level: intermediate

.seealso: SVDGetConvergedReason(), SVDSetTolerances()
@*/
PetscErrorCode SVDGetIterationNumber(SVD svd,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidIntPointer(its,2);
  *its = svd->its;
  PetscFunctionReturn(0);
}

/*@
   SVDGetConvergedReason - Gets the reason why the SVDSolve() iteration was
   stopped.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  reason - negative value indicates diverged, positive value converged
   (see SVDConvergedReason)

   Notes:

   Possible values for reason are
+  SVD_CONVERGED_TOL - converged up to tolerance
.  SVD_CONVERGED_USER - converged due to a user-defined condition
.  SVD_DIVERGED_ITS - required more than max_it iterations to reach convergence
-  SVD_DIVERGED_BREAKDOWN - generic breakdown in method

   Can only be called after the call to SVDSolve() is complete.

   Level: intermediate

.seealso: SVDSetTolerances(), SVDSolve(), SVDConvergedReason
@*/
PetscErrorCode SVDGetConvergedReason(SVD svd,SVDConvergedReason *reason)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidIntPointer(reason,2);
  SVDCheckSolved(svd,1);
  *reason = svd->reason;
  PetscFunctionReturn(0);
}

/*@
   SVDGetConverged - Gets the number of converged singular values.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  nconv - number of converged singular values

   Note:
   This function should be called after SVDSolve() has finished.

   Level: beginner

@*/
PetscErrorCode SVDGetConverged(SVD svd,PetscInt *nconv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidIntPointer(nconv,2);
  SVDCheckSolved(svd,1);
  *nconv = svd->nconv;
  PetscFunctionReturn(0);
}

/*@C
   SVDGetSingularTriplet - Gets the i-th triplet of the singular value decomposition
   as computed by SVDSolve(). The solution consists in the singular value and its left
   and right singular vectors.

   Not Collective, but vectors are shared by all processors that share the SVD

   Input Parameters:
+  svd - singular value solver context
-  i   - index of the solution

   Output Parameters:
+  sigma - singular value
.  u     - left singular vector
-  v     - right singular vector

   Note:
   Both U or V can be NULL if singular vectors are not required.
   Otherwise, the caller must provide valid Vec objects, i.e.,
   they must be created by the calling program with e.g. MatCreateVecs().

   The index i should be a value between 0 and nconv-1 (see SVDGetConverged()).
   Singular triplets are indexed according to the ordering criterion established
   with SVDSetWhichSingularTriplets().

   Level: beginner

.seealso: SVDSolve(), SVDGetConverged(), SVDSetWhichSingularTriplets()
@*/
PetscErrorCode SVDGetSingularTriplet(SVD svd,PetscInt i,PetscReal *sigma,Vec u,Vec v)
{
  PetscErrorCode ierr;
  PetscInt       M,N;
  Vec            w;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveInt(svd,i,2);
  SVDCheckSolved(svd,1);
  if (u) { PetscValidHeaderSpecific(u,VEC_CLASSID,4); PetscCheckSameComm(svd,1,u,4); }
  if (v) { PetscValidHeaderSpecific(v,VEC_CLASSID,5); PetscCheckSameComm(svd,1,v,5); }
  if (i<0 || i>=svd->nconv) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  if (sigma) *sigma = svd->sigma[svd->perm[i]];
  ierr = MatGetSize(svd->OP,&M,&N);CHKERRQ(ierr);
  if (M<N) { w = u; u = v; v = w; }
  if (u) {
    ierr = SVDComputeVectors(svd);CHKERRQ(ierr);
    ierr = BVCopyVec(svd->U,svd->perm[i],u);CHKERRQ(ierr);
  }
  if (v) {
    ierr = BVCopyVec(svd->V,svd->perm[i],v);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   SVDComputeResidualNorms_Private - Computes the norms of the left and
   right residuals associated with the i-th computed singular triplet.
@*/
static PetscErrorCode SVDComputeResidualNorms_Private(SVD svd,PetscInt i,PetscReal *norm1,PetscReal *norm2)
{
  PetscErrorCode ierr;
  Vec            u,v,x = NULL,y = NULL;
  PetscReal      sigma;
  PetscInt       M,N;

  PetscFunctionBegin;
  ierr = MatCreateVecs(svd->OP,&v,&u);CHKERRQ(ierr);
  ierr = SVDGetSingularTriplet(svd,i,&sigma,u,v);CHKERRQ(ierr);
  /* norm1 = ||A*v-sigma*u||_2 */
  if (norm1) {
    ierr = VecDuplicate(u,&x);CHKERRQ(ierr);
    ierr = MatMult(svd->OP,v,x);CHKERRQ(ierr);
    ierr = VecAXPY(x,-sigma,u);CHKERRQ(ierr);
    ierr = VecNorm(x,NORM_2,norm1);CHKERRQ(ierr);
  }
  /* norm2 = ||A^T*u-sigma*v||_2 */
  if (norm2) {
    ierr = VecDuplicate(v,&y);CHKERRQ(ierr);
    if (svd->A && svd->AT) {
      ierr = MatGetSize(svd->OP,&M,&N);CHKERRQ(ierr);
      if (M<N) {
        ierr = MatMult(svd->A,u,y);CHKERRQ(ierr);
      } else {
        ierr = MatMult(svd->AT,u,y);CHKERRQ(ierr);
      }
    } else {
#if defined(PETSC_USE_COMPLEX)
      ierr = MatMultHermitianTranspose(svd->OP,u,y);CHKERRQ(ierr);
#else
      ierr = MatMultTranspose(svd->OP,u,y);CHKERRQ(ierr);
#endif
    }
    ierr = VecAXPY(y,-sigma,v);CHKERRQ(ierr);
    ierr = VecNorm(y,NORM_2,norm2);CHKERRQ(ierr);
  }

  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   SVDComputeError - Computes the error (based on the residual norm) associated
   with the i-th singular triplet.

   Collective on svd

   Input Parameter:
+  svd  - the singular value solver context
.  i    - the solution index
-  type - the type of error to compute

   Output Parameter:
.  error - the error

   Notes:
   The error can be computed in various ways, all of them based on the residual
   norm obtained as sqrt(n1^2+n2^2) with n1 = ||A*v-sigma*u||_2 and
   n2 = ||A^T*u-sigma*v||_2, where sigma is the singular value, u is the left
   singular vector and v is the right singular vector.

   Level: beginner

.seealso: SVDErrorType, SVDSolve()
@*/
PetscErrorCode SVDComputeError(SVD svd,PetscInt i,SVDErrorType type,PetscReal *error)
{
  PetscErrorCode ierr;
  PetscReal      sigma,norm1,norm2;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveInt(svd,i,2);
  PetscValidLogicalCollectiveEnum(svd,type,3);
  PetscValidRealPointer(error,4);
  SVDCheckSolved(svd,1);
  ierr = SVDGetSingularTriplet(svd,i,&sigma,NULL,NULL);CHKERRQ(ierr);
  ierr = SVDComputeResidualNorms_Private(svd,i,&norm1,&norm2);CHKERRQ(ierr);
  *error = PetscSqrtReal(norm1*norm1+norm2*norm2);
  switch (type) {
    case SVD_ERROR_ABSOLUTE:
      break;
    case SVD_ERROR_RELATIVE:
      *error /= sigma;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Invalid error type");
  }
  PetscFunctionReturn(0);
}

