/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to the solution process
*/

#include <slepc/private/nepimpl.h>       /*I "slepcnep.h" I*/
#include <slepc/private/bvimpl.h>        /*I "slepcbv.h" I*/
#include <petscdraw.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-nep,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"{NEP}: a module for the parallel solution of nonlinear eigenvalue problems in {SLEPc}\",\n"
  "   journal = \"arXiv:1910.11712\",\n"
  "   url = \"https://arxiv.org/abs/1910.11712\",\n"
  "   year = \"2019\"\n"
  "}\n";

PetscErrorCode NEPComputeVectors(NEP nep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  NEPCheckSolved(nep,1);
  if (nep->state==NEP_STATE_SOLVED && nep->ops->computevectors) {
    ierr = (*nep->ops->computevectors)(nep);CHKERRQ(ierr);
  }
  nep->state = NEP_STATE_EIGENVECTORS;
  PetscFunctionReturn(0);
}

/*@
   NEPSolve - Solves the nonlinear eigensystem.

   Collective on nep

   Input Parameter:
.  nep - eigensolver context obtained from NEPCreate()

   Options Database Keys:
+  -nep_view - print information about the solver used
.  -nep_view_vectors binary - save the computed eigenvectors to the default binary viewer
.  -nep_view_values - print computed eigenvalues
.  -nep_converged_reason - print reason for convergence, and number of iterations
.  -nep_error_absolute - print absolute errors of each eigenpair
-  -nep_error_relative - print relative errors of each eigenpair

   Level: beginner

.seealso: NEPCreate(), NEPSetUp(), NEPDestroy(), NEPSetTolerances()
@*/
PetscErrorCode NEPSolve(NEP nep)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (nep->state>=NEP_STATE_SOLVED) PetscFunctionReturn(0);
  ierr = PetscCitationsRegister(citation,&cited);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(NEP_Solve,nep,0,0,0);CHKERRQ(ierr);

  /* call setup */
  ierr = NEPSetUp(nep);CHKERRQ(ierr);
  nep->nconv = 0;
  nep->its = 0;
  for (i=0;i<nep->ncv;i++) {
    nep->eigr[i]   = 0.0;
    nep->eigi[i]   = 0.0;
    nep->errest[i] = 0.0;
    nep->perm[i]   = i;
  }
  ierr = NEPViewFromOptions(nep,NULL,"-nep_view_pre");CHKERRQ(ierr);
  ierr = RGViewFromOptions(nep->rg,NULL,"-rg_view");CHKERRQ(ierr);

  ierr = (*nep->ops->solve)(nep);CHKERRQ(ierr);
  nep->state = NEP_STATE_SOLVED;

  if (!nep->reason) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_PLIB,"Internal error, solver returned without setting converged reason");

  if (nep->refine==NEP_REFINE_SIMPLE && nep->rits>0 && nep->nconv>0) {
    ierr = NEPComputeVectors(nep);CHKERRQ(ierr);
    ierr = NEPNewtonRefinementSimple(nep,&nep->rits,nep->rtol,nep->nconv);CHKERRQ(ierr);
    nep->state = NEP_STATE_EIGENVECTORS;
  }

  /* sort eigenvalues according to nep->which parameter */
  ierr = SlepcSortEigenvalues(nep->sc,nep->nconv,nep->eigr,nep->eigi,nep->perm);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(NEP_Solve,nep,0,0,0);CHKERRQ(ierr);

  /* various viewers */
  ierr = NEPViewFromOptions(nep,NULL,"-nep_view");CHKERRQ(ierr);
  ierr = NEPReasonViewFromOptions(nep);CHKERRQ(ierr);
  ierr = NEPErrorViewFromOptions(nep);CHKERRQ(ierr);
  ierr = NEPValuesViewFromOptions(nep);CHKERRQ(ierr);
  ierr = NEPVectorsViewFromOptions(nep);CHKERRQ(ierr);

  /* Remove the initial subspace */
  nep->nini = 0;

  /* Reset resolvent information */
  ierr = MatDestroy(&nep->resolvent);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   NEPProjectOperator - Computes the projection of the nonlinear operator.

   Collective on nep

   Input Parameters:
+  nep - the nonlinear eigensolver context
.  j0  - initial index
-  j1  - final index

   Notes:
   This is available for split operator only.

   The nonlinear operator T(lambda) is projected onto span(V), where V is
   an orthonormal basis built internally by the solver. The projected
   operator is equal to sum_i V'*A_i*V*f_i(lambda), so this function
   computes all matrices Ei = V'*A_i*V, and stores them in the extra
   matrices inside DS. Only rows/columns in the range [j0,j1-1] are computed,
   the previous ones are assumed to be available already.

   Level: developer

.seealso: NEPSetSplitOperator()
@*/
PetscErrorCode NEPProjectOperator(NEP nep,PetscInt j0,PetscInt j1)
{
  PetscErrorCode ierr;
  PetscInt       k;
  Mat            G;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,j0,2);
  PetscValidLogicalCollectiveInt(nep,j1,3);
  NEPCheckProblem(nep,1);
  NEPCheckSplit(nep,1);
  ierr = BVSetActiveColumns(nep->V,j0,j1);CHKERRQ(ierr);
  for (k=0;k<nep->nt;k++) {
    ierr = DSGetMat(nep->ds,DSMatExtra[k],&G);CHKERRQ(ierr);
    ierr = BVMatProject(nep->V,nep->A[k],nep->V,G);CHKERRQ(ierr);
    ierr = DSRestoreMat(nep->ds,DSMatExtra[k],&G);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   NEPApplyFunction - Applies the nonlinear function T(lambda) to a given vector.

   Collective on nep

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  lambda - scalar argument
.  x      - vector to be multiplied against
-  v      - workspace vector (used only in the case of split form)

   Output Parameters:
+  y   - result vector
.  A   - Function matrix
-  B   - optional preconditioning matrix

   Note:
   If the nonlinear operator is represented in split form, the result
   y = T(lambda)*x is computed without building T(lambda) explicitly. In
   that case, parameters A and B are not used. Otherwise, the matrix
   T(lambda) is built and the effect is the same as a call to
   NEPComputeFunction() followed by a MatMult().

   Level: developer

.seealso: NEPSetSplitOperator(), NEPComputeFunction(), NEPApplyAdjoint()
@*/
PetscErrorCode NEPApplyFunction(NEP nep,PetscScalar lambda,Vec x,Vec v,Vec y,Mat A,Mat B)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    alpha;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(nep,lambda,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  if (v) PetscValidHeaderSpecific(v,VEC_CLASSID,4);
  PetscValidHeaderSpecific(y,VEC_CLASSID,5);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,6);
  if (B) PetscValidHeaderSpecific(B,MAT_CLASSID,7);

  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    ierr = VecSet(y,0.0);CHKERRQ(ierr);
    for (i=0;i<nep->nt;i++) {
      ierr = FNEvaluateFunction(nep->f[i],lambda,&alpha);CHKERRQ(ierr);
      ierr = MatMult(nep->A[i],x,v);CHKERRQ(ierr);
      ierr = VecAXPY(y,alpha,v);CHKERRQ(ierr);
    }
  } else {
    ierr = NEPComputeFunction(nep,lambda,A,B);CHKERRQ(ierr);
    ierr = MatMult(A,x,y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   NEPApplyAdjoint - Applies the adjoint nonlinear function T(lambda)^* to a given vector.

   Collective on nep

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  lambda - scalar argument
.  x      - vector to be multiplied against
-  v      - workspace vector (used only in the case of split form)

   Output Parameters:
+  y   - result vector
.  A   - Function matrix
-  B   - optional preconditioning matrix

   Level: developer

.seealso: NEPSetSplitOperator(), NEPComputeFunction(), NEPApplyFunction()
@*/
PetscErrorCode NEPApplyAdjoint(NEP nep,PetscScalar lambda,Vec x,Vec v,Vec y,Mat A,Mat B)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    alpha;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(nep,lambda,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  if (v) PetscValidHeaderSpecific(v,VEC_CLASSID,4);
  PetscValidHeaderSpecific(y,VEC_CLASSID,5);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,6);
  if (B) PetscValidHeaderSpecific(B,MAT_CLASSID,7);

  ierr = VecConjugate(x);CHKERRQ(ierr);
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    ierr = VecSet(y,0.0);CHKERRQ(ierr);
    for (i=0;i<nep->nt;i++) {
      ierr = FNEvaluateFunction(nep->f[i],lambda,&alpha);CHKERRQ(ierr);
      ierr = MatMultTranspose(nep->A[i],x,v);CHKERRQ(ierr);
      ierr = VecAXPY(y,alpha,v);CHKERRQ(ierr);
    }
  } else {
    ierr = NEPComputeFunction(nep,lambda,A,B);CHKERRQ(ierr);
    ierr = MatMultTranspose(A,x,y);CHKERRQ(ierr);
  }
  ierr = VecConjugate(x);CHKERRQ(ierr);
  ierr = VecConjugate(y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   NEPApplyJacobian - Applies the nonlinear Jacobian T'(lambda) to a given vector.

   Collective on nep

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  lambda - scalar argument
.  x      - vector to be multiplied against
-  v      - workspace vector (used only in the case of split form)

   Output Parameters:
+  y   - result vector
-  A   - Jacobian matrix

   Note:
   If the nonlinear operator is represented in split form, the result
   y = T'(lambda)*x is computed without building T'(lambda) explicitly. In
   that case, parameter A is not used. Otherwise, the matrix
   T'(lambda) is built and the effect is the same as a call to
   NEPComputeJacobian() followed by a MatMult().

   Level: developer

.seealso: NEPSetSplitOperator(), NEPComputeJacobian()
@*/
PetscErrorCode NEPApplyJacobian(NEP nep,PetscScalar lambda,Vec x,Vec v,Vec y,Mat A)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    alpha;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveScalar(nep,lambda,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  if (v) PetscValidHeaderSpecific(v,VEC_CLASSID,4);
  PetscValidHeaderSpecific(y,VEC_CLASSID,5);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,6);

  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    ierr = VecSet(y,0.0);CHKERRQ(ierr);
    for (i=0;i<nep->nt;i++) {
      ierr = FNEvaluateDerivative(nep->f[i],lambda,&alpha);CHKERRQ(ierr);
      ierr = MatMult(nep->A[i],x,v);CHKERRQ(ierr);
      ierr = VecAXPY(y,alpha,v);CHKERRQ(ierr);
    }
  } else {
    ierr = NEPComputeJacobian(nep,lambda,A);CHKERRQ(ierr);
    ierr = MatMult(A,x,y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   NEPGetIterationNumber - Gets the current iteration number. If the
   call to NEPSolve() is complete, then it returns the number of iterations
   carried out by the solution method.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  its - number of iterations

   Note:
   During the i-th iteration this call returns i-1. If NEPSolve() is
   complete, then parameter "its" contains either the iteration number at
   which convergence was successfully reached, or failure was detected.
   Call NEPGetConvergedReason() to determine if the solver converged or
   failed and why.

   Level: intermediate

.seealso: NEPGetConvergedReason(), NEPSetTolerances()
@*/
PetscErrorCode NEPGetIterationNumber(NEP nep,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidIntPointer(its,2);
  *its = nep->its;
  PetscFunctionReturn(0);
}

/*@
   NEPGetConverged - Gets the number of converged eigenpairs.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  nconv - number of converged eigenpairs

   Note:
   This function should be called after NEPSolve() has finished.

   Level: beginner

.seealso: NEPSetDimensions(), NEPSolve()
@*/
PetscErrorCode NEPGetConverged(NEP nep,PetscInt *nconv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidIntPointer(nconv,2);
  NEPCheckSolved(nep,1);
  *nconv = nep->nconv;
  PetscFunctionReturn(0);
}

/*@
   NEPGetConvergedReason - Gets the reason why the NEPSolve() iteration was
   stopped.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  reason - negative value indicates diverged, positive value converged

   Notes:

   Possible values for reason are
+  NEP_CONVERGED_TOL - converged up to tolerance
.  NEP_CONVERGED_USER - converged due to a user-defined condition
.  NEP_DIVERGED_ITS - required more than max_it iterations to reach convergence
.  NEP_DIVERGED_BREAKDOWN - generic breakdown in method
.  NEP_DIVERGED_LINEAR_SOLVE - inner linear solve failed
-  NEP_DIVERGED_SUBSPACE_EXHAUSTED - run out of space for the basis in an
   unrestarted solver

   Can only be called after the call to NEPSolve() is complete.

   Level: intermediate

.seealso: NEPSetTolerances(), NEPSolve(), NEPConvergedReason
@*/
PetscErrorCode NEPGetConvergedReason(NEP nep,NEPConvergedReason *reason)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(reason,2);
  NEPCheckSolved(nep,1);
  *reason = nep->reason;
  PetscFunctionReturn(0);
}

/*@C
   NEPGetEigenpair - Gets the i-th solution of the eigenproblem as computed by
   NEPSolve(). The solution consists in both the eigenvalue and the eigenvector.

   Logically Collective on nep

   Input Parameters:
+  nep - nonlinear eigensolver context
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
   set to zero). In any case, the user can pass NULL in Vr or Vi if one of
   them is not required.

   The index i should be a value between 0 and nconv-1 (see NEPGetConverged()).
   Eigenpairs are indexed according to the ordering criterion established
   with NEPSetWhichEigenpairs().

   Level: beginner

.seealso: NEPSolve(), NEPGetConverged(), NEPSetWhichEigenpairs(), NEPGetLeftEigenvector()
@*/
PetscErrorCode NEPGetEigenpair(NEP nep,PetscInt i,PetscScalar *eigr,PetscScalar *eigi,Vec Vr,Vec Vi)
{
  PetscInt       k;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,i,2);
  if (Vr) { PetscValidHeaderSpecific(Vr,VEC_CLASSID,5); PetscCheckSameComm(nep,1,Vr,5); }
  if (Vi) { PetscValidHeaderSpecific(Vi,VEC_CLASSID,6); PetscCheckSameComm(nep,1,Vi,6); }
  NEPCheckSolved(nep,1);
  if (i<0 || i>=nep->nconv) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");

  ierr = NEPComputeVectors(nep);CHKERRQ(ierr);
  k = nep->perm[i];

  /* eigenvalue */
#if defined(PETSC_USE_COMPLEX)
  if (eigr) *eigr = nep->eigr[k];
  if (eigi) *eigi = 0;
#else
  if (eigr) *eigr = nep->eigr[k];
  if (eigi) *eigi = nep->eigi[k];
#endif

  /* eigenvector */
  ierr = BV_GetEigenvector(nep->V,k,nep->eigi[k],Vr,Vi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   NEPGetLeftEigenvector - Gets the i-th left eigenvector as computed by NEPSolve().

   Logically Collective on nep

   Input Parameters:
+  nep - eigensolver context
-  i   - index of the solution

   Output Parameters:
+  Wr   - real part of left eigenvector
-  Wi   - imaginary part of left eigenvector

   Notes:
   The caller must provide valid Vec objects, i.e., they must be created
   by the calling program with e.g. MatCreateVecs().

   If the corresponding eigenvalue is real, then Wi is set to zero. If PETSc is
   configured with complex scalars the eigenvector is stored directly in Wr
   (Wi is set to zero). In any case, the user can pass NULL in Wr or Wi if one of
   them is not required.

   The index i should be a value between 0 and nconv-1 (see NEPGetConverged()).
   Eigensolutions are indexed according to the ordering criterion established
   with NEPSetWhichEigenpairs().

   Left eigenvectors are available only if the twosided flag was set, see
   NEPSetTwoSided().

   Level: intermediate

.seealso: NEPGetEigenpair(), NEPGetConverged(), NEPSetWhichEigenpairs(), NEPSetTwoSided()
@*/
PetscErrorCode NEPGetLeftEigenvector(NEP nep,PetscInt i,Vec Wr,Vec Wi)
{
  PetscErrorCode ierr;
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,i,2);
  if (Wr) { PetscValidHeaderSpecific(Wr,VEC_CLASSID,3); PetscCheckSameComm(nep,1,Wr,3); }
  if (Wi) { PetscValidHeaderSpecific(Wi,VEC_CLASSID,4); PetscCheckSameComm(nep,1,Wi,4); }
  NEPCheckSolved(nep,1);
  if (!nep->twosided) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONGSTATE,"Must request left vectors with NEPSetTwoSided");
  if (i<0 || i>=nep->nconv) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  ierr = NEPComputeVectors(nep);CHKERRQ(ierr);
  k = nep->perm[i];
  ierr = BV_GetEigenvector(nep->W,k,nep->eigi[k],Wr,Wi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   NEPGetErrorEstimate - Returns the error estimate associated to the i-th
   computed eigenpair.

   Not Collective

   Input Parameter:
+  nep - nonlinear eigensolver context
-  i   - index of eigenpair

   Output Parameter:
.  errest - the error estimate

   Notes:
   This is the error estimate used internally by the eigensolver. The actual
   error bound can be computed with NEPComputeRelativeError().

   Level: advanced

.seealso: NEPComputeRelativeError()
@*/
PetscErrorCode NEPGetErrorEstimate(NEP nep,PetscInt i,PetscReal *errest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidRealPointer(errest,3);
  NEPCheckSolved(nep,1);
  if (i<0 || i>=nep->nconv) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Argument 2 out of range");
  *errest = nep->errest[nep->perm[i]];
  PetscFunctionReturn(0);
}

/*
   NEPComputeResidualNorm_Private - Computes the norm of the residual vector
   associated with an eigenpair.

   Input Parameters:
     adj    - whether the adjoint T^* must be used instead of T
     lambda - eigenvalue
     x      - eigenvector
     w      - array of work vectors (two vectors in split form, one vector otherwise)
*/
PetscErrorCode NEPComputeResidualNorm_Private(NEP nep,PetscBool adj,PetscScalar lambda,Vec x,Vec *w,PetscReal *norm)
{
  PetscErrorCode ierr;
  Vec            y,z=NULL;

  PetscFunctionBegin;
  y = w[0];
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) z = w[1];
  if (adj) {
    ierr = NEPApplyAdjoint(nep,lambda,x,z,y,nep->function,nep->function_pre);CHKERRQ(ierr);
  } else {
    ierr = NEPApplyFunction(nep,lambda,x,z,y,nep->function,nep->function_pre);CHKERRQ(ierr);
  }
  ierr = VecNorm(y,NORM_2,norm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   NEPComputeError - Computes the error (based on the residual norm) associated
   with the i-th computed eigenpair.

   Collective on nep

   Input Parameter:
+  nep  - the nonlinear eigensolver context
.  i    - the solution index
-  type - the type of error to compute

   Output Parameter:
.  error - the error

   Notes:
   The error can be computed in various ways, all of them based on the residual
   norm computed as ||T(lambda)x||_2 where lambda is the eigenvalue and x is the
   eigenvector.

   Level: beginner

.seealso: NEPErrorType, NEPSolve(), NEPGetErrorEstimate()
@*/
PetscErrorCode NEPComputeError(NEP nep,PetscInt i,NEPErrorType type,PetscReal *error)
{
  PetscErrorCode ierr;
  Vec            xr,xi=NULL;
  PetscInt       j,nwork,issplit=0;
  PetscScalar    kr,ki,s;
  PetscReal      er,z=0.0,errorl;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,i,2);
  PetscValidLogicalCollectiveEnum(nep,type,3);
  PetscValidRealPointer(error,4);
  NEPCheckSolved(nep,1);

  /* allocate work vectors */
#if defined(PETSC_USE_COMPLEX)
  nwork = 2;
#else
  nwork = 3;
#endif
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    issplit = 1;
    nwork++;  /* need an extra work vector for NEPComputeResidualNorm_Private */
  }
  ierr = NEPSetWorkVecs(nep,nwork);CHKERRQ(ierr);
  xr = nep->work[issplit+1];
#if !defined(PETSC_USE_COMPLEX)
  xi = nep->work[issplit+2];
#endif

  /* compute residual norms */
  ierr = NEPGetEigenpair(nep,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  if (ki) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Not implemented for complex eigenvalues with real scalars");
#endif
  ierr = NEPComputeResidualNorm_Private(nep,PETSC_FALSE,kr,xr,nep->work,error);CHKERRQ(ierr);
  ierr = VecNorm(xr,NORM_2,&er);CHKERRQ(ierr);

  /* if two-sided, compute left residual norm and take the maximum */
  if (nep->twosided) {
    ierr = NEPGetLeftEigenvector(nep,i,xr,xi);CHKERRQ(ierr);
    ierr = NEPComputeResidualNorm_Private(nep,PETSC_TRUE,kr,xr,nep->work,&errorl);CHKERRQ(ierr);
    *error = PetscMax(*error,errorl);
  }

  /* compute error */
  switch (type) {
    case NEP_ERROR_ABSOLUTE:
      break;
    case NEP_ERROR_RELATIVE:
      *error /= PetscAbsScalar(kr)*er;
      break;
    case NEP_ERROR_BACKWARD:
      if (nep->fui!=NEP_USER_INTERFACE_SPLIT) {
        *error = 0.0;
        ierr = PetscInfo(nep,"Backward error only available in split form\n");CHKERRQ(ierr);
        break;
      }
      /* initialization of matrix norms */
      if (!nep->nrma[0]) {
        for (j=0;j<nep->nt;j++) {
          ierr = MatHasOperation(nep->A[j],MATOP_NORM,&flg);CHKERRQ(ierr);
          if (!flg) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"The computation of backward errors requires a matrix norm operation");
          ierr = MatNorm(nep->A[j],NORM_INFINITY,&nep->nrma[j]);CHKERRQ(ierr);
        }
      }
      for (j=0;j<nep->nt;j++) {
        ierr = FNEvaluateFunction(nep->f[j],kr,&s);CHKERRQ(ierr);
        z = z + nep->nrma[j]*PetscAbsScalar(s);
      }
      *error /= z;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid error type");
  }
  PetscFunctionReturn(0);
}

/*@
   NEPComputeFunction - Computes the function matrix T(lambda) that has been
   set with NEPSetFunction().

   Collective on nep

   Input Parameters:
+  nep    - the NEP context
-  lambda - the scalar argument

   Output Parameters:
+  A   - Function matrix
-  B   - optional preconditioning matrix

   Notes:
   NEPComputeFunction() is typically used within nonlinear eigensolvers
   implementations, so most users would not generally call this routine
   themselves.

   Level: developer

.seealso: NEPSetFunction(), NEPGetFunction()
@*/
PetscErrorCode NEPComputeFunction(NEP nep,PetscScalar lambda,Mat A,Mat B)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    alpha;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckProblem(nep,1);
  switch (nep->fui) {
  case NEP_USER_INTERFACE_CALLBACK:
    if (!nep->computefunction) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_USER,"Must call NEPSetFunction() first");
    ierr = PetscLogEventBegin(NEP_FunctionEval,nep,A,B,0);CHKERRQ(ierr);
    PetscStackPush("NEP user Function function");
    ierr = (*nep->computefunction)(nep,lambda,A,B,nep->functionctx);CHKERRQ(ierr);
    PetscStackPop;
    ierr = PetscLogEventEnd(NEP_FunctionEval,nep,A,B,0);CHKERRQ(ierr);
    break;
  case NEP_USER_INTERFACE_SPLIT:
    ierr = MatZeroEntries(A);CHKERRQ(ierr);
    for (i=0;i<nep->nt;i++) {
      ierr = FNEvaluateFunction(nep->f[i],lambda,&alpha);CHKERRQ(ierr);
      ierr = MatAXPY(A,alpha,nep->A[i],nep->mstr);CHKERRQ(ierr);
    }
    if (A != B) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Not implemented");
    break;
  }
  PetscFunctionReturn(0);
}

/*@
   NEPComputeJacobian - Computes the Jacobian matrix T'(lambda) that has been
   set with NEPSetJacobian().

   Collective on nep

   Input Parameters:
+  nep    - the NEP context
-  lambda - the scalar argument

   Output Parameters:
.  A   - Jacobian matrix

   Notes:
   Most users should not need to explicitly call this routine, as it
   is used internally within the nonlinear eigensolvers.

   Level: developer

.seealso: NEPSetJacobian(), NEPGetJacobian()
@*/
PetscErrorCode NEPComputeJacobian(NEP nep,PetscScalar lambda,Mat A)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    alpha;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckProblem(nep,1);
  switch (nep->fui) {
  case NEP_USER_INTERFACE_CALLBACK:
    if (!nep->computejacobian) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_USER,"Must call NEPSetJacobian() first");
    ierr = PetscLogEventBegin(NEP_JacobianEval,nep,A,0,0);CHKERRQ(ierr);
    PetscStackPush("NEP user Jacobian function");
    ierr = (*nep->computejacobian)(nep,lambda,A,nep->jacobianctx);CHKERRQ(ierr);
    PetscStackPop;
    ierr = PetscLogEventEnd(NEP_JacobianEval,nep,A,0,0);CHKERRQ(ierr);
    break;
  case NEP_USER_INTERFACE_SPLIT:
    ierr = MatZeroEntries(A);CHKERRQ(ierr);
    for (i=0;i<nep->nt;i++) {
      ierr = FNEvaluateDerivative(nep->f[i],lambda,&alpha);CHKERRQ(ierr);
      ierr = MatAXPY(A,alpha,nep->A[i],nep->mstr);CHKERRQ(ierr);
    }
    break;
  }
  PetscFunctionReturn(0);
}

