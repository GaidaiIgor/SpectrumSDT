/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   LME routines related to problem setup
*/

#include <slepc/private/lmeimpl.h>       /*I "slepclme.h" I*/

PETSC_STATIC_INLINE PetscErrorCode LMESetUp_Lyapunov(LME lme)
{
  PetscErrorCode ierr;
  Mat            C1,C2,X1,X2;
  Vec            dc,dx;

  PetscFunctionBegin;
  ierr = MatLRCGetMats(lme->C,NULL,&C1,&dc,&C2);CHKERRQ(ierr);
  if (C1!=C2) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"Lyapunov matrix equation requires symmetric right-hand side C");
  if (dc) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"Lyapunov solvers currently require positive-definite right-hand side C");
  if (lme->X) {
    ierr = MatLRCGetMats(lme->X,NULL,&X1,&dx,&X2);CHKERRQ(ierr);
    if (X1!=X2) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"Lyapunov matrix equation requires symmetric solution X");
    if (dx) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"Lyapunov solvers currently assume a positive-definite solution X");
  }
  PetscFunctionReturn(0);
}

/*@
   LMESetUp - Sets up all the internal data structures necessary for the
   execution of the linear matrix equation solver.

   Collective on lme

   Input Parameter:
.  lme   - linear matrix equation solver context

   Notes:
   This function need not be called explicitly in most cases, since LMESolve()
   calls it. It can be useful when one wants to measure the set-up time
   separately from the solve time.

   Level: developer

.seealso: LMECreate(), LMESolve(), LMEDestroy()
@*/
PetscErrorCode LMESetUp(LME lme)
{
  PetscErrorCode ierr;
  PetscInt       N;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);

  /* reset the convergence flag from the previous solves */
  lme->reason = LME_CONVERGED_ITERATING;

  if (lme->setupcalled) PetscFunctionReturn(0);
  ierr = PetscLogEventBegin(LME_SetUp,lme,0,0,0);CHKERRQ(ierr);

  /* Set default solver type (LMESetFromOptions was not called) */
  if (!((PetscObject)lme)->type_name) {
    ierr = LMESetType(lme,LMEKRYLOV);CHKERRQ(ierr);
  }

  /* Check problem dimensions */
  if (!lme->A) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONGSTATE,"LMESetCoefficients must be called first");
  ierr = MatGetSize(lme->A,&N,NULL);CHKERRQ(ierr);
  if (lme->ncv > N) lme->ncv = N;

  /* setup options for the particular equation type */
  switch (lme->problem_type) {
    case LME_LYAPUNOV:
      ierr = LMESetUp_Lyapunov(lme);CHKERRQ(ierr);
      break;
    case LME_SYLVESTER:
      LMECheckCoeff(lme,lme->B,"B","Sylvester");
      break;
    case LME_GEN_LYAPUNOV:
      LMECheckCoeff(lme,lme->D,"D","Generalized Lyapunov");
      break;
    case LME_GEN_SYLVESTER:
      LMECheckCoeff(lme,lme->B,"B","Generalized Sylvester");
      LMECheckCoeff(lme,lme->D,"D","Generalized Sylvester");
      LMECheckCoeff(lme,lme->E,"E","Generalized Sylvester");
      break;
    case LME_DT_LYAPUNOV:
      break;
    case LME_STEIN:
      LMECheckCoeff(lme,lme->D,"D","Stein");
      break;
  }
  if (lme->problem_type!=LME_LYAPUNOV) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_SUP,"There is no solver yet for this matrix equation type");

  /* call specific solver setup */
  ierr = (*lme->ops->setup)(lme);CHKERRQ(ierr);

  /* set tolerance if not yet set */
  if (lme->tol==PETSC_DEFAULT) lme->tol = SLEPC_DEFAULT_TOL;

  ierr = PetscLogEventEnd(LME_SetUp,lme,0,0,0);CHKERRQ(ierr);
  lme->setupcalled = 1;
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode LMESetCoefficients_Private(LME lme,Mat A,Mat *lmeA)
{
  PetscErrorCode ierr;
  PetscInt       m,n;

  PetscFunctionBegin;
  ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
  if (m!=n) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONG,"Matrix is non-square");
  if (!lme->setupcalled) { ierr = MatDestroy(lmeA);CHKERRQ(ierr); }
  ierr = PetscObjectReference((PetscObject)A);CHKERRQ(ierr);
  *lmeA = A;
  PetscFunctionReturn(0);
}

/*@
   LMESetCoefficients - Sets the coefficient matrices that define the linear matrix
   equation to be solved.

   Collective on lme

   Input Parameters:
+  lme - the matrix function context
.  A   - first coefficient matrix
.  B   - second coefficient matrix
.  D   - third coefficient matrix
-  E   - fourth coefficient matrix

   Notes:
   The matrix equation takes the general form A*X*E+D*X*B=C, where matrix C is not
   provided here but with LMESetRHS(). Not all four matrices must be passed, some
   can be NULL instead, see LMESetProblemType() for details.

   It must be called before LMESetUp(). If it is called again after LMESetUp() then
   the LME object is reset.

   In order to delete a previously set matrix, pass a NULL in the corresponding
   argument.

   Level: beginner

.seealso: LMESolve(), LMESetUp(), LMESetRHS(), LMESetProblemType()
@*/
PetscErrorCode LMESetCoefficients(LME lme,Mat A,Mat B,Mat D,Mat E)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(lme,1,A,2);
  if (B) {
    PetscValidHeaderSpecific(B,MAT_CLASSID,3);
    PetscCheckSameComm(lme,1,B,3);
  }
  if (D) {
    PetscValidHeaderSpecific(D,MAT_CLASSID,4);
    PetscCheckSameComm(lme,1,D,4);
  }
  if (E) {
    PetscValidHeaderSpecific(E,MAT_CLASSID,5);
    PetscCheckSameComm(lme,1,E,5);
  }

  if (lme->setupcalled) { ierr = LMEReset(lme);CHKERRQ(ierr); }

  ierr = LMESetCoefficients_Private(lme,A,&lme->A);CHKERRQ(ierr);
  if (B) { ierr = LMESetCoefficients_Private(lme,B,&lme->B);CHKERRQ(ierr); }
  else if (!lme->setupcalled) { ierr = MatDestroy(&lme->B);CHKERRQ(ierr); }
  if (D) { ierr = LMESetCoefficients_Private(lme,D,&lme->D);CHKERRQ(ierr); }
  else if (!lme->setupcalled) { ierr = MatDestroy(&lme->D);CHKERRQ(ierr); }
  if (E) { ierr = LMESetCoefficients_Private(lme,E,&lme->E);CHKERRQ(ierr); }
  else if (!lme->setupcalled) { ierr = MatDestroy(&lme->E);CHKERRQ(ierr); }

  lme->setupcalled = 0;
  PetscFunctionReturn(0);
}

/*@
   LMEGetCoefficients - Gets the coefficient matrices of the matrix equation.

   Collective on lme

   Input Parameter:
.  lme - the LME context

   Output Parameters:
+  A   - first coefficient matrix
.  B   - second coefficient matrix
.  D   - third coefficient matrix
-  E   - fourth coefficient matrix

   Level: intermediate

.seealso: LMESolve(), LMESetCoefficients()
@*/
PetscErrorCode LMEGetCoefficients(LME lme,Mat *A,Mat *B,Mat *D,Mat *E)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (A) *A = lme->A;
  if (B) *B = lme->B;
  if (D) *D = lme->D;
  if (E) *E = lme->E;
  PetscFunctionReturn(0);
}

/*@
   LMESetRHS - Sets the right-hand side of the matrix equation, as a low-rank
   matrix.

   Collective on lme

   Input Parameters:
+  lme - the matrix function context
-  C   - the right-hand side matrix

   Notes:
   The matrix equation takes the general form A*X*E+D*X*B=C, where matrix C is
   given with this function. C must be a low-rank matrix of type MATLRC, that is,
   C = U*D*V' where D is diagonal of order k, and U, V are dense tall-skinny
   matrices with k columns. No sparse matrix must be provided when creating the
   MATLRC matrix.

   In equation types that require C to be symmetric, such as Lyapunov, C must be
   created with V=U (or V=NULL).

   Level: beginner

.seealso: LMESetSolution(), LMESetProblemType()
@*/
PetscErrorCode LMESetRHS(LME lme,Mat C)
{
  PetscErrorCode ierr;
  PetscBool      match;
  Mat            A;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidHeaderSpecific(C,MAT_CLASSID,2);
  PetscCheckSameComm(lme,1,C,2);

  ierr = PetscObjectTypeCompare((PetscObject)C,MATLRC,&match);CHKERRQ(ierr);
  if (!match) SETERRQ(PetscObjectComm((PetscObject)C),PETSC_ERR_SUP,"Mat argument must have been created with MatCreateLRC");
  ierr = MatLRCGetMats(C,&A,NULL,NULL,NULL);CHKERRQ(ierr);
  if (A) SETERRQ(PetscObjectComm((PetscObject)C),PETSC_ERR_SUP,"The MatLRC must not have a sparse matrix term");

  ierr = PetscObjectReference((PetscObject)C);CHKERRQ(ierr);
  ierr = MatDestroy(&lme->C);CHKERRQ(ierr);
  lme->C = C;
  PetscFunctionReturn(0);
}

/*@
   LMEGetRHS - Gets the right-hand side of the matrix equation.

   Collective on lme

   Input Parameter:
.  lme - the LME context

   Output Parameters:
.  C   - the low-rank matrix

   Level: intermediate

.seealso: LMESolve(), LMESetRHS()
@*/
PetscErrorCode LMEGetRHS(LME lme,Mat *C)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidPointer(C,2);
  *C = lme->C;
  PetscFunctionReturn(0);
}

/*@
   LMESetSolution - Sets the placeholder for the solution of the matrix
   equation, as a low-rank matrix.

   Collective on lme

   Input Parameters:
+  lme - the matrix function context
-  X   - the solution matrix

   Notes:
   The matrix equation takes the general form A*X*E+D*X*B=C, where the solution
   matrix is of low rank and is written in factored form X = U*D*V'. This function
   provides a Mat object of type MATLRC that stores U, V and (optionally) D.
   These factors will be computed during LMESolve().

   In equation types whose solution X is symmetric, such as Lyapunov, X must be
   created with V=U (or V=NULL).

   If the user provides X with this function, then the solver will
   return a solution with rank at most the number of columns of U. Alternatively,
   it is possible to let the solver choose the rank of the solution, by
   setting X to NULL and then calling LMEGetSolution() after LMESolve().

   Level: intermediate

.seealso: LMEGetSolution(), LMESetRHS(), LMESetProblemType(), LMESolve()
@*/
PetscErrorCode LMESetSolution(LME lme,Mat X)
{
  PetscErrorCode ierr;
  PetscBool      match;
  Mat            A;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (X) {
    PetscValidHeaderSpecific(X,MAT_CLASSID,2);
    PetscCheckSameComm(lme,1,X,2);
    ierr = PetscObjectTypeCompare((PetscObject)X,MATLRC,&match);CHKERRQ(ierr);
    if (!match) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_SUP,"Mat argument must have been created with MatCreateLRC");
    ierr = MatLRCGetMats(X,&A,NULL,NULL,NULL);CHKERRQ(ierr);
    if (A) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_SUP,"The MatLRC must not have a sparse matrix term");
    ierr = PetscObjectReference((PetscObject)X);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&lme->X);CHKERRQ(ierr);
  lme->X = X;
  PetscFunctionReturn(0);
}

/*@
   LMEGetSolution - Gets the solution of the matrix equation.

   Collective on lme

   Input Parameter:
.  lme - the LME context

   Output Parameters:
.  X   - the low-rank matrix

   Level: intermediate

.seealso: LMESolve(), LMESetSolution()
@*/
PetscErrorCode LMEGetSolution(LME lme,Mat *X)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidPointer(X,2);
  *X = lme->X;
  PetscFunctionReturn(0);
}

/*@
   LMEAllocateSolution - Allocate memory storage for common variables such
   as the basis vectors.

   Collective on lme

   Input Parameters:
+  lme   - linear matrix equation solver context
-  extra - number of additional positions, used for methods that require a
           working basis slightly larger than ncv

   Developers Note:
   This is SLEPC_EXTERN because it may be required by user plugin LME
   implementations.

   Level: developer
@*/
PetscErrorCode LMEAllocateSolution(LME lme,PetscInt extra)
{
  PetscErrorCode ierr;
  PetscInt       oldsize,requested;
  Vec            t;

  PetscFunctionBegin;
  requested = lme->ncv + extra;

  /* oldsize is zero if this is the first time setup is called */
  ierr = BVGetSizes(lme->V,NULL,NULL,&oldsize);CHKERRQ(ierr);

  /* allocate basis vectors */
  if (!lme->V) { ierr = LMEGetBV(lme,&lme->V);CHKERRQ(ierr); }
  if (!oldsize) {
    if (!((PetscObject)(lme->V))->type_name) {
      ierr = BVSetType(lme->V,BVSVEC);CHKERRQ(ierr);
    }
    ierr = MatCreateVecsEmpty(lme->A,&t,NULL);CHKERRQ(ierr);
    ierr = BVSetSizesFromVec(lme->V,t,requested);CHKERRQ(ierr);
    ierr = VecDestroy(&t);CHKERRQ(ierr);
  } else {
    ierr = BVResize(lme->V,requested,PETSC_FALSE);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

