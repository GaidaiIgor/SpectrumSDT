/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   MFN routines related to the solution process
*/

#include <slepc/private/mfnimpl.h>   /*I "slepcmfn.h" I*/

static PetscErrorCode MFNSolve_Private(MFN mfn,Vec b,Vec x)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecSetErrorIfLocked(x,3);CHKERRQ(ierr);

  /* call setup */
  ierr = MFNSetUp(mfn);CHKERRQ(ierr);
  mfn->its = 0;

  ierr = MFNViewFromOptions(mfn,NULL,"-mfn_view_pre");CHKERRQ(ierr);

  /* check nonzero right-hand side */
  ierr = VecNorm(b,NORM_2,&mfn->bnorm);CHKERRQ(ierr);
  if (!mfn->bnorm) SETERRQ(PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_WRONG,"Cannot pass a zero b vector to MFNSolve()");

  /* call solver */
  ierr = PetscLogEventBegin(MFN_Solve,mfn,b,x,0);CHKERRQ(ierr);
  if (b!=x) { ierr = VecLockReadPush(b);CHKERRQ(ierr); }
  ierr = (*mfn->ops->solve)(mfn,b,x);CHKERRQ(ierr);
  if (b!=x) { ierr = VecLockReadPop(b);CHKERRQ(ierr); }
  ierr = PetscLogEventEnd(MFN_Solve,mfn,b,x,0);CHKERRQ(ierr);

  if (!mfn->reason) SETERRQ(PetscObjectComm((PetscObject)mfn),PETSC_ERR_PLIB,"Internal error, solver returned without setting converged reason");

  if (mfn->errorifnotconverged && mfn->reason < 0) SETERRQ(PetscObjectComm((PetscObject)mfn),PETSC_ERR_NOT_CONVERGED,"MFNSolve has not converged");

  /* various viewers */
  ierr = MFNViewFromOptions(mfn,NULL,"-mfn_view");CHKERRQ(ierr);
  ierr = MFNReasonViewFromOptions(mfn);CHKERRQ(ierr);
  ierr = MatViewFromOptions(mfn->A,(PetscObject)mfn,"-mfn_view_mat");CHKERRQ(ierr);
  ierr = VecViewFromOptions(b,(PetscObject)mfn,"-mfn_view_rhs");CHKERRQ(ierr);
  ierr = VecViewFromOptions(x,(PetscObject)mfn,"-mfn_view_solution");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   MFNSolve - Solves the matrix function problem. Given a vector b, the
   vector x = f(A)*b is returned.

   Collective on mfn

   Input Parameters:
+  mfn - matrix function context obtained from MFNCreate()
-  b   - the right hand side vector

   Output Parameter:
.  x   - the solution (this may be the same vector as b, then b will be
         overwritten with the answer)

   Options Database Keys:
+  -mfn_view - print information about the solver used
.  -mfn_view_mat binary - save the matrix to the default binary viewer
.  -mfn_view_rhs binary - save right hand side vector to the default binary viewer
.  -mfn_view_solution binary - save computed solution vector to the default binary viewer
-  -mfn_converged_reason - print reason for convergence, and number of iterations

   Notes:
   The matrix A is specified with MFNSetOperator().
   The function f is specified with MFNSetFN().

   Level: beginner

.seealso: MFNCreate(), MFNSetUp(), MFNDestroy(), MFNSetTolerances(),
          MFNSetOperator(), MFNSetFN()
@*/
PetscErrorCode MFNSolve(MFN mfn,Vec b,Vec x)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidHeaderSpecific(b,VEC_CLASSID,2);
  PetscCheckSameComm(mfn,1,b,2);
  if (b!=x) PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  if (b!=x) PetscCheckSameComm(mfn,1,x,3);
  mfn->transpose_solve = PETSC_FALSE;
  ierr = MFNSolve_Private(mfn,b,x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   MFNSolveTranspose - Solves the transpose matrix function problem. Given a vector b,
   the vector x = f(A^T)*b is returned.

   Collective on mfn

   Input Parameters:
+  mfn - matrix function context obtained from MFNCreate()
-  b   - the right hand side vector

   Output Parameter:
.  x   - the solution (this may be the same vector as b, then b will be
         overwritten with the answer)

   Note:
   See available options at MFNSolve().

   Level: beginner

.seealso: MFNSolve()
@*/
PetscErrorCode MFNSolveTranspose(MFN mfn,Vec b,Vec x)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidHeaderSpecific(b,VEC_CLASSID,2);
  PetscCheckSameComm(mfn,1,b,2);
  if (b!=x) PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  if (b!=x) PetscCheckSameComm(mfn,1,x,3);
  mfn->transpose_solve = PETSC_TRUE;
  if (!mfn->AT) { ierr = MatCreateTranspose(mfn->A,&mfn->AT);CHKERRQ(ierr); }
  ierr = MFNSolve_Private(mfn,b,x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   MFNGetIterationNumber - Gets the current iteration number. If the
   call to MFNSolve() is complete, then it returns the number of iterations
   carried out by the solution method.

   Not Collective

   Input Parameter:
.  mfn - the matrix function context

   Output Parameter:
.  its - number of iterations

   Note:
   During the i-th iteration this call returns i-1. If MFNSolve() is
   complete, then parameter "its" contains either the iteration number at
   which convergence was successfully reached, or failure was detected.
   Call MFNGetConvergedReason() to determine if the solver converged or
   failed and why.

   Level: intermediate

.seealso: MFNGetConvergedReason(), MFNSetTolerances()
@*/
PetscErrorCode MFNGetIterationNumber(MFN mfn,PetscInt *its)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidIntPointer(its,2);
  *its = mfn->its;
  PetscFunctionReturn(0);
}

/*@
   MFNGetConvergedReason - Gets the reason why the MFNSolve() iteration was
   stopped.

   Not Collective

   Input Parameter:
.  mfn - the matrix function context

   Output Parameter:
.  reason - negative value indicates diverged, positive value converged

   Notes:

   Possible values for reason are
+  MFN_CONVERGED_TOL - converged up to tolerance
.  MFN_CONVERGED_ITS - solver completed the requested number of steps
.  MFN_DIVERGED_ITS - required more than max_it iterations to reach convergence
-  MFN_DIVERGED_BREAKDOWN - generic breakdown in method

   Can only be called after the call to MFNSolve() is complete.

   Basic solvers (e.g. unrestarted Krylov iterations) cannot determine if the
   computation is accurate up to the requested tolerance. In that case, the
   converged reason is set to MFN_CONVERGED_ITS if the requested number of steps
   (for instance, the ncv value in unrestarted Krylov methods) have been
   completed successfully.

   Level: intermediate

.seealso: MFNSetTolerances(), MFNSolve(), MFNConvergedReason, MFNSetErrorIfNotConverged()
@*/
PetscErrorCode MFNGetConvergedReason(MFN mfn,MFNConvergedReason *reason)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidIntPointer(reason,2);
  *reason = mfn->reason;
  PetscFunctionReturn(0);
}

