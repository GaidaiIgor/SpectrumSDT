/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   LME routines related to the solution process
*/

#include <slepc/private/lmeimpl.h>   /*I "slepclme.h" I*/
#include <slepcblaslapack.h>

/*@
   LMESolve - Solves the linear matrix equation.

   Collective on lme

   Input Parameters:
.  lme - linear matrix equation solver context obtained from LMECreate()

   Options Database Keys:
+  -lme_view - print information about the solver used
.  -lme_view_mat binary - save the matrix to the default binary viewer
.  -lme_view_rhs binary - save right hand side to the default binary viewer
.  -lme_view_solution binary - save computed solution to the default binary viewer
-  -lme_converged_reason - print reason for convergence, and number of iterations

   Notes:
   The matrix coefficients are specified with LMESetCoefficients().
   The right-hand side is specified with LMESetRHS().
   The placeholder for the solution is specified with LMESetSolution().

   Level: beginner

.seealso: LMECreate(), LMESetUp(), LMEDestroy(), LMESetTolerances(), LMESetCoefficients(), LMESetRHS(), LMESetSolution()
@*/
PetscErrorCode LMESolve(LME lme)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);

  /* call setup */
  ierr = LMESetUp(lme);CHKERRQ(ierr);
  lme->its    = 0;
  lme->errest = 0.0;

  ierr = LMEViewFromOptions(lme,NULL,"-lme_view_pre");CHKERRQ(ierr);

  /* call solver */
  if (!lme->ops->solve[lme->problem_type]) SETERRQ1(PetscObjectComm((PetscObject)lme),PETSC_ERR_SUP,"The specified solver does not support equation type %s",LMEProblemTypes[lme->problem_type]);
  ierr = PetscLogEventBegin(LME_Solve,lme,0,0,0);CHKERRQ(ierr);
  ierr = (*lme->ops->solve[lme->problem_type])(lme);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(LME_Solve,lme,0,0,0);CHKERRQ(ierr);

  if (!lme->reason) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_PLIB,"Internal error, solver returned without setting converged reason");

  if (lme->errorifnotconverged && lme->reason < 0) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_NOT_CONVERGED,"LMESolve has not converged");

  /* various viewers */
  ierr = LMEViewFromOptions(lme,NULL,"-lme_view");CHKERRQ(ierr);
  ierr = LMEReasonViewFromOptions(lme);CHKERRQ(ierr);
  ierr = MatViewFromOptions(lme->A,(PetscObject)lme,"-lme_view_mat");CHKERRQ(ierr);
  ierr = MatViewFromOptions(lme->C,(PetscObject)lme,"-lme_view_rhs");CHKERRQ(ierr);
  ierr = MatViewFromOptions(lme->X,(PetscObject)lme,"-lme_view_solution");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   LMEGetIterationNumber - Gets the current iteration number. If the
   call to LMESolve() is complete, then it returns the number of iterations
   carried out by the solution method.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  its - number of iterations

   Note:
   During the i-th iteration this call returns i-1. If LMESolve() is
   complete, then parameter "its" contains either the iteration number at
   which convergence was successfully reached, or failure was detected.
   Call LMEGetConvergedReason() to determine if the solver converged or
   failed and why.

   Level: intermediate

.seealso: LMEGetConvergedReason(), LMESetTolerances()
@*/
PetscErrorCode LMEGetIterationNumber(LME lme,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidIntPointer(its,2);
  *its = lme->its;
  PetscFunctionReturn(0);
}

/*@
   LMEGetConvergedReason - Gets the reason why the LMESolve() iteration was
   stopped.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  reason - negative value indicates diverged, positive value converged

   Notes:

   Possible values for reason are
+  LME_CONVERGED_TOL - converged up to tolerance
.  LME_DIVERGED_ITS - required more than max_it iterations to reach convergence
-  LME_DIVERGED_BREAKDOWN - generic breakdown in method

   Can only be called after the call to LMESolve() is complete.

   Level: intermediate

.seealso: LMESetTolerances(), LMESolve(), LMEConvergedReason, LMESetErrorIfNotConverged()
@*/
PetscErrorCode LMEGetConvergedReason(LME lme,LMEConvergedReason *reason)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidIntPointer(reason,2);
  *reason = lme->reason;
  PetscFunctionReturn(0);
}

/*@
   LMEGetErrorEstimate - Returns the error estimate obtained during solve.

   Not Collective

   Input Parameter:
.  lme - linear matrix equation solver context

   Output Parameter:
.  errest - the error estimate

   Notes:
   This is the error estimated internally by the solver. The actual
   error bound can be computed with LMEComputeError(). Note that some
   solvers may not be able to provide an error estimate.

   Level: advanced

.seealso: LMEComputeError()
@*/
PetscErrorCode LMEGetErrorEstimate(LME lme,PetscReal *errest)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidRealPointer(errest,2);
  *errest = lme->errest;
  PetscFunctionReturn(0);
}

/*
   LMEComputeResidualNorm_Lyapunov - Computes the Frobenius norm of the residual matrix
   associated with the Lyapunov equation.
*/
PetscErrorCode LMEComputeResidualNorm_Lyapunov(LME lme,PetscReal *norm)
{
  PetscErrorCode    ierr;
  PetscInt          j,n,N,k,l;
  PetscBLASInt      n_,N_,k_,l_;
  PetscScalar       *Rarray,alpha=1.0,beta=0.0;
  const PetscScalar *A,*B;
  BV                W,AX,X1,C1;
  Mat               R,X1m,C1m;
  Vec               v,w;
  VecScatter        vscat;

  PetscFunctionBegin;
  ierr = MatLRCGetMats(lme->C,NULL,&C1m,NULL,NULL);CHKERRQ(ierr);
  ierr = MatLRCGetMats(lme->X,NULL,&X1m,NULL,NULL);CHKERRQ(ierr);
  ierr = BVCreateFromMat(C1m,&C1);CHKERRQ(ierr);
  ierr = BVSetFromOptions(C1);CHKERRQ(ierr);
  ierr = BVCreateFromMat(X1m,&X1);CHKERRQ(ierr);
  ierr = BVSetFromOptions(X1);CHKERRQ(ierr);
  ierr = BVGetSizes(X1,&n,&N,&k);CHKERRQ(ierr);
  ierr = BVGetSizes(C1,NULL,NULL,&l);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(N,&N_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(l,&l_);CHKERRQ(ierr);

  /* create W to store a redundant copy of a BV in each process */
  ierr = BVCreate(PETSC_COMM_SELF,&W);CHKERRQ(ierr);
  ierr = BVSetSizes(W,N,N,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(W);CHKERRQ(ierr);
  ierr = BVGetColumn(X1,0,&v);CHKERRQ(ierr);
  ierr = VecScatterCreateToAll(v,&vscat,NULL);CHKERRQ(ierr);
  ierr = BVRestoreColumn(X1,0,&v);CHKERRQ(ierr);

  /* create AX to hold the product A*X1 */
  ierr = BVDuplicate(X1,&AX);CHKERRQ(ierr);
  ierr = BVMatMult(X1,lme->A,AX);CHKERRQ(ierr);

  /* create dense matrix to hold the residual R=C1*C1'+AX*X1'+X1*AX' */
  ierr = MatCreateDense(PetscObjectComm((PetscObject)lme),n,n,N,N,NULL,&R);CHKERRQ(ierr);

  /* R=C1*C1' */
  ierr = MatDenseGetArray(R,&Rarray);CHKERRQ(ierr);
  for (j=0;j<l;j++) {
    ierr = BVGetColumn(C1,j,&v);CHKERRQ(ierr);
    ierr = BVGetColumn(W,j,&w);CHKERRQ(ierr);
    ierr = VecScatterBegin(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = BVRestoreColumn(C1,j,&v);CHKERRQ(ierr);
    ierr = BVRestoreColumn(W,j,&w);CHKERRQ(ierr);
  }
  if (n) {
    ierr = BVGetArrayRead(C1,&A);CHKERRQ(ierr);
    ierr = BVGetArrayRead(W,&B);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&n_,&N_,&l_,&alpha,(PetscScalar*)A,&n_,(PetscScalar*)B,&N_,&beta,Rarray,&n_));
    ierr = BVRestoreArrayRead(C1,&A);CHKERRQ(ierr);
    ierr = BVRestoreArrayRead(W,&B);CHKERRQ(ierr);
  }
  beta = 1.0;

  /* R+=AX*X1' */
  for (j=0;j<k;j++) {
    ierr = BVGetColumn(X1,j,&v);CHKERRQ(ierr);
    ierr = BVGetColumn(W,j,&w);CHKERRQ(ierr);
    ierr = VecScatterBegin(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = BVRestoreColumn(X1,j,&v);CHKERRQ(ierr);
    ierr = BVRestoreColumn(W,j,&w);CHKERRQ(ierr);
  }
  if (n) {
    ierr = BVGetArrayRead(AX,&A);CHKERRQ(ierr);
    ierr = BVGetArrayRead(W,&B);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&n_,&N_,&k_,&alpha,(PetscScalar*)A,&n_,(PetscScalar*)B,&N_,&beta,Rarray,&n_));
    ierr = BVRestoreArrayRead(AX,&A);CHKERRQ(ierr);
    ierr = BVRestoreArrayRead(W,&B);CHKERRQ(ierr);
  }

  /* R+=X1*AX' */
  for (j=0;j<k;j++) {
    ierr = BVGetColumn(AX,j,&v);CHKERRQ(ierr);
    ierr = BVGetColumn(W,j,&w);CHKERRQ(ierr);
    ierr = VecScatterBegin(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(vscat,v,w,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = BVRestoreColumn(AX,j,&v);CHKERRQ(ierr);
    ierr = BVRestoreColumn(W,j,&w);CHKERRQ(ierr);
  }
  if (n) {
    ierr = BVGetArrayRead(X1,&A);CHKERRQ(ierr);
    ierr = BVGetArrayRead(W,&B);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&n_,&N_,&k_,&alpha,(PetscScalar*)A,&n_,(PetscScalar*)B,&N_,&beta,Rarray,&n_));
    ierr = BVRestoreArrayRead(X1,&A);CHKERRQ(ierr);
    ierr = BVRestoreArrayRead(W,&B);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(R,&Rarray);CHKERRQ(ierr);

  /* compute ||R||_F */
  ierr = MatAssemblyBegin(R,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(R,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatNorm(R,NORM_FROBENIUS,norm);CHKERRQ(ierr);

  ierr = BVDestroy(&W);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&vscat);CHKERRQ(ierr);
  ierr = BVDestroy(&AX);CHKERRQ(ierr);
  ierr = MatDestroy(&R);CHKERRQ(ierr);
  ierr = BVDestroy(&C1);CHKERRQ(ierr);
  ierr = BVDestroy(&X1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   LMEComputeError - Computes the error (based on the residual norm) associated
   with the last equation solved.

   Collective on lme

   Input Parameter:
.  lme  - the linear matrix equation solver context

   Output Parameter:
.  error - the error

   Notes:
   This function is not scalable (in terms of memory or parallel communication),
   so it should not be called except in the case of small problem size. For
   large equations, use LMEGetErrorEstimate().

   Level: advanced

.seealso: LMESolve(), LMEGetErrorEstimate()
@*/
PetscErrorCode LMEComputeError(LME lme,PetscReal *error)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidRealPointer(error,2);

  ierr = PetscLogEventBegin(LME_ComputeError,lme,0,0,0);CHKERRQ(ierr);
  /* compute residual norm */
  switch (lme->problem_type) {
    case LME_LYAPUNOV:
      ierr = LMEComputeResidualNorm_Lyapunov(lme,error);CHKERRQ(ierr);
      break;
    default:
      SETERRQ1(PetscObjectComm((PetscObject)lme),PETSC_ERR_SUP,"Not implemented for equation type %s",LMEProblemTypes[lme->problem_type]);
  }

  /* compute error */
  /* currently we only support absolute error, so just return the norm */
  ierr = PetscLogEventEnd(LME_ComputeError,lme,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

