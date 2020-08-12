/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   MFN routines related to problem setup
*/

#include <slepc/private/mfnimpl.h>       /*I "slepcmfn.h" I*/

/*@
   MFNSetUp - Sets up all the internal data structures necessary for the
   execution of the matrix function solver.

   Collective on mfn

   Input Parameter:
.  mfn   - matrix function context

   Notes:
   This function need not be called explicitly in most cases, since MFNSolve()
   calls it. It can be useful when one wants to measure the set-up time
   separately from the solve time.

   Level: developer

.seealso: MFNCreate(), MFNSolve(), MFNDestroy()
@*/
PetscErrorCode MFNSetUp(MFN mfn)
{
  PetscErrorCode ierr;
  PetscInt       N;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);

  /* reset the convergence flag from the previous solves */
  mfn->reason = MFN_CONVERGED_ITERATING;

  if (mfn->setupcalled) PetscFunctionReturn(0);
  ierr = PetscLogEventBegin(MFN_SetUp,mfn,0,0,0);CHKERRQ(ierr);

  /* Set default solver type (MFNSetFromOptions was not called) */
  if (!((PetscObject)mfn)->type_name) {
    ierr = MFNSetType(mfn,MFNKRYLOV);CHKERRQ(ierr);
  }
  if (!mfn->fn) { ierr = MFNGetFN(mfn,&mfn->fn);CHKERRQ(ierr); }
  if (!((PetscObject)mfn->fn)->type_name) {
    ierr = FNSetFromOptions(mfn->fn);CHKERRQ(ierr);
  }

  /* Check problem dimensions */
  if (!mfn->A) SETERRQ(PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_WRONGSTATE,"MFNSetOperator must be called first");
  ierr = MatGetSize(mfn->A,&N,NULL);CHKERRQ(ierr);
  if (mfn->ncv > N) mfn->ncv = N;

  /* call specific solver setup */
  ierr = (*mfn->ops->setup)(mfn);CHKERRQ(ierr);

  /* set tolerance if not yet set */
  if (mfn->tol==PETSC_DEFAULT) mfn->tol = SLEPC_DEFAULT_TOL;

  ierr = PetscLogEventEnd(MFN_SetUp,mfn,0,0,0);CHKERRQ(ierr);
  mfn->setupcalled = 1;
  PetscFunctionReturn(0);
}

/*@
   MFNSetOperator - Sets the matrix for which the matrix function is to be computed.

   Collective on mfn

   Input Parameters:
+  mfn - the matrix function context
-  A   - the problem matrix

   Notes:
   It must be called before MFNSetUp(). If it is called again after MFNSetUp() then
   the MFN object is reset.

   Level: beginner

.seealso: MFNSolve(), MFNSetUp(), MFNReset()
@*/
PetscErrorCode MFNSetOperator(MFN mfn,Mat A)
{
  PetscErrorCode ierr;
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(mfn,1,A,2);

  ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
  if (m!=n) SETERRQ(PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_WRONG,"A is a non-square matrix");
  ierr = PetscObjectReference((PetscObject)A);CHKERRQ(ierr);
  if (mfn->setupcalled) { ierr = MFNReset(mfn);CHKERRQ(ierr); }
  else { ierr = MatDestroy(&mfn->A);CHKERRQ(ierr); }
  mfn->A = A;
  mfn->setupcalled = 0;
  PetscFunctionReturn(0);
}

/*@
   MFNGetOperator - Gets the matrix associated with the MFN object.

   Collective on mfn

   Input Parameter:
.  mfn - the MFN context

   Output Parameters:
.  A  - the matrix for which the matrix function is to be computed

   Level: intermediate

.seealso: MFNSolve(), MFNSetOperator()
@*/
PetscErrorCode MFNGetOperator(MFN mfn,Mat *A)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidPointer(A,2);
  *A = mfn->A;
  PetscFunctionReturn(0);
}

/*@
   MFNAllocateSolution - Allocate memory storage for common variables such
   as the basis vectors.

   Collective on mfn

   Input Parameters:
+  mfn   - matrix function context
-  extra - number of additional positions, used for methods that require a
           working basis slightly larger than ncv

   Developers Note:
   This is SLEPC_EXTERN because it may be required by user plugin MFN
   implementations.

   Level: developer
@*/
PetscErrorCode MFNAllocateSolution(MFN mfn,PetscInt extra)
{
  PetscErrorCode ierr;
  PetscInt       oldsize,requested;
  Vec            t;

  PetscFunctionBegin;
  requested = mfn->ncv + extra;

  /* oldsize is zero if this is the first time setup is called */
  ierr = BVGetSizes(mfn->V,NULL,NULL,&oldsize);CHKERRQ(ierr);

  /* allocate basis vectors */
  if (!mfn->V) { ierr = MFNGetBV(mfn,&mfn->V);CHKERRQ(ierr); }
  if (!oldsize) {
    if (!((PetscObject)(mfn->V))->type_name) {
      ierr = BVSetType(mfn->V,BVSVEC);CHKERRQ(ierr);
    }
    ierr = MatCreateVecsEmpty(mfn->A,&t,NULL);CHKERRQ(ierr);
    ierr = BVSetSizesFromVec(mfn->V,t,requested);CHKERRQ(ierr);
    ierr = VecDestroy(&t);CHKERRQ(ierr);
  } else {
    ierr = BVResize(mfn->V,requested,PETSC_FALSE);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

