/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates the use of shell spectral transformations. "
  "The problem to be solved is the same as ex1.c and"
  "corresponds to the Laplacian operator in 1 dimension.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

#include <slepceps.h>

/* Define context for user-provided spectral transformation */
typedef struct {
  KSP        ksp;
} SampleShellST;

/* Declare routines for user-provided spectral transformation */
PetscErrorCode STCreate_User(SampleShellST**);
PetscErrorCode STSetUp_User(SampleShellST*,ST);
PetscErrorCode STApply_User(ST,Vec,Vec);
PetscErrorCode STApplyTranspose_User(ST,Vec,Vec);
PetscErrorCode STBackTransform_User(ST,PetscInt,PetscScalar*,PetscScalar*);
PetscErrorCode STDestroy_User(SampleShellST*);

int main (int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  ST             st;              /* spectral transformation context */
  SampleShellST  *shell;          /* user-defined spectral transform context */
  EPSType        type;
  PetscInt       n=30,i,Istart,Iend,nev;
  PetscBool      isShell,terse;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem (shell-enabled), n=%D\n\n",n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(A,i,i-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A,i,i+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,i,i,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /*
     Initialize shell spectral transformation if selected by user
  */
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)st,STSHELL,&isShell);CHKERRQ(ierr);
  if (isShell) {
    /* Change sorting criterion since this ST example computes values
       closest to 0 */
    ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);

    /* (Required) Create a context for the user-defined spectral transform;
       this context can be defined to contain any application-specific data. */
    ierr = STCreate_User(&shell);CHKERRQ(ierr);
    ierr = STShellSetContext(st,shell);CHKERRQ(ierr);

    /* (Required) Set the user-defined routine for applying the operator */
    ierr = STShellSetApply(st,STApply_User);CHKERRQ(ierr);

    /* (Optional) Set the user-defined routine for applying the transposed operator */
    ierr = STShellSetApplyTranspose(st,STApplyTranspose_User);CHKERRQ(ierr);

    /* (Optional) Set the user-defined routine for back-transformation */
    ierr = STShellSetBackTransform(st,STBackTransform_User);CHKERRQ(ierr);

    /* (Optional) Set a name for the transformation, used for STView() */
    ierr = PetscObjectSetName((PetscObject)st,"MyTransformation");CHKERRQ(ierr);

    /* (Optional) Do any setup required for the new transformation */
    ierr = STSetUp_User(shell,st);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  if (isShell) {
    ierr = STDestroy_User(shell);CHKERRQ(ierr);
  }
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/***********************************************************************/
/*     Routines for a user-defined shell spectral transformation       */
/***********************************************************************/

/*
   STCreate_User - This routine creates a user-defined
   spectral transformation context.

   Output Parameter:
.  shell - user-defined spectral transformation context
*/
PetscErrorCode STCreate_User(SampleShellST **shell)
{
  SampleShellST  *newctx;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr   = PetscNew(&newctx);CHKERRQ(ierr);
  ierr   = KSPCreate(PETSC_COMM_WORLD,&newctx->ksp);CHKERRQ(ierr);
  ierr   = KSPAppendOptionsPrefix(newctx->ksp,"st_");CHKERRQ(ierr);
  *shell = newctx;
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
/*
   STSetUp_User - This routine sets up a user-defined
   spectral transformation context.

   Input Parameters:
+  shell - user-defined spectral transformation context
-  st    - spectral transformation context containing the operator matrices

   Output Parameter:
.  shell - fully set up user-defined transformation context

   Notes:
   In this example, the user-defined transformation is simply OP=A^-1.
   Therefore, the eigenpairs converge in reversed order. The KSP object
   used for the solution of linear systems with A is handled via the
   user-defined context SampleShellST.
*/
PetscErrorCode STSetUp_User(SampleShellST *shell,ST st)
{
  Mat            A;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = STGetMatrix(st,0,&A);CHKERRQ(ierr);
  ierr = KSPSetOperators(shell->ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(shell->ksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
/*
   STApply_User - This routine demonstrates the use of a
   user-provided spectral transformation.

   Input Parameters:
+  st - spectral transformation context
-  x - input vector

   Output Parameter:
.  y - output vector

   Notes:
   The transformation implemented in this code is just OP=A^-1 and
   therefore it is of little use, merely as an example of working with
   a STSHELL.
*/
PetscErrorCode STApply_User(ST st,Vec x,Vec y)
{
  SampleShellST  *shell;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = STShellGetContext(st,(void**)&shell);CHKERRQ(ierr);
  ierr = KSPSolve(shell->ksp,x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
/*
   STApplyTranspose_User - This is not required unless using a two-sided
   eigensolver.

   Input Parameters:
+  st - spectral transformation context
-  x - input vector

   Output Parameter:
.  y - output vector
*/
PetscErrorCode STApplyTranspose_User(ST st,Vec x,Vec y)
{
  SampleShellST  *shell;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = STShellGetContext(st,(void**)&shell);CHKERRQ(ierr);
  ierr = KSPSolveTranspose(shell->ksp,x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
/*
   STBackTransform_User - This routine demonstrates the use of a
   user-provided spectral transformation.

   Input Parameters:
+  st - spectral transformation context
-  n  - number of eigenvalues to transform

   Input/Output Parameters:
+  eigr - pointer to real part of eigenvalues
-  eigi - pointer to imaginary part of eigenvalues

   Notes:
   This code implements the back transformation of eigenvalues in
   order to retrieve the eigenvalues of the original problem. In this
   example, simply set k_i = 1/k_i.
*/
PetscErrorCode STBackTransform_User(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscInt j;

  PetscFunctionBeginUser;
  for (j=0;j<n;j++) {
    eigr[j] = 1.0 / eigr[j];
  }
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */
/*
   STDestroy_User - This routine destroys a user-defined
   spectral transformation context.

   Input Parameter:
.  shell - user-defined spectral transformation context
*/
PetscErrorCode STDestroy_User(SampleShellST *shell)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = KSPDestroy(&shell->ksp);CHKERRQ(ierr);
  ierr = PetscFree(shell);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   testset:
      args: -eps_nev 5 -eps_non_hermitian -terse
      requires: !single
      output_file: output/ex10_1.out
      test:
         suffix: 1_sinvert
         args: -st_type sinvert
      test:
         suffix: 1_sinvert_twoside
         args: -st_type sinvert -eps_balance twoside
      test:
         suffix: 1_shell
         args: -st_type shell
      test:
         suffix: 1_shell_twoside
         args: -st_type shell -eps_balance twoside

TEST*/
