/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Shell spectral transformations with a non-injective mapping. "
  "Implements spectrum folding for the 2-D Laplacian, as in ex24.c.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n";

#include <slepceps.h>

/* Context for spectrum folding spectral transformation */
typedef struct {
  Mat         A;
  Vec         w;
  PetscScalar target;
} FoldShellST;

/* Routines for shell spectral transformation */
PetscErrorCode STCreate_Fold(Mat,PetscScalar,FoldShellST**);
PetscErrorCode STApply_Fold(ST,Vec,Vec);
PetscErrorCode STDestroy_Fold(FoldShellST*);

int main (int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  ST             st;              /* spectral transformation context */
  FoldShellST    *fold;           /* user-defined spectral transform context */
  EPSType        type;
  PetscInt       N,n=10,m,i,j,II,Istart,Iend,nev;
  PetscBool      isShell,terse,flag;
  PetscScalar    target=1.1;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag);CHKERRQ(ierr);
  if (!flag) m = n;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-target",&target,NULL);CHKERRQ(ierr);
  N = n*m;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSpectrum Folding via shell ST, N=%D (%Dx%D grid) target=%3.2f\n\n",N,n,m,(double)PetscRealPart(target));CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the 5-point stencil Laplacian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) { ierr = MatSetValue(A,II,II-n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<m-1) { ierr = MatSetValue(A,II,II+n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j>0) { ierr = MatSetValue(A,II,II-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(A,II,II+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,II,II,4.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  ierr = EPSSetTarget(eps,target);CHKERRQ(ierr);
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = STSetType(st,STSHELL);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /*
     Initialize shell spectral transformation
  */
  ierr = PetscObjectTypeCompare((PetscObject)st,STSHELL,&isShell);CHKERRQ(ierr);
  if (isShell) {
    /* Change sorting criterion since this shell ST computes eigenvalues
       of the transformed operator closest to 0 */
    ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);

    /* Create the context for the user-defined spectral transform */
    ierr = STCreate_Fold(A,target,&fold);CHKERRQ(ierr);
    ierr = STShellSetContext(st,fold);CHKERRQ(ierr);

    /* Set callback function for applying the operator (in this case we do not
       provide a back-transformation callback since the mapping is not one-to-one) */
    ierr = STShellSetApply(st,STApply_Fold);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)st,"STFOLD");CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);
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
    ierr = STDestroy_Fold(fold);CHKERRQ(ierr);
  }
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*
   STCreate_Fold - Creates the spectrum folding ST context.

   Input Parameter:
+  A - problem matrix
-  target - target value

   Output Parameter:
.  fold - user-defined spectral transformation context
*/
PetscErrorCode STCreate_Fold(Mat A,PetscScalar target,FoldShellST **fold)
{
  FoldShellST    *newctx;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscNew(&newctx);CHKERRQ(ierr);
  newctx->A = A;
  ierr = PetscObjectReference((PetscObject)A);CHKERRQ(ierr);
  newctx->target = target;
  ierr = MatCreateVecs(A,&newctx->w,NULL);CHKERRQ(ierr);
  *fold = newctx;
  PetscFunctionReturn(0);
}

/*
   STApply_Fold - Applies the operator (A-target*I)^2 to a given vector.

   Input Parameters:
+  st - spectral transformation context
-  x  - input vector

   Output Parameter:
.  y - output vector
*/
PetscErrorCode STApply_Fold(ST st,Vec x,Vec y)
{
  FoldShellST    *fold;
  PetscScalar    sigma;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = STShellGetContext(st,(void**)&fold);CHKERRQ(ierr);
  sigma = -fold->target;
  ierr = MatMult(fold->A,x,fold->w);CHKERRQ(ierr);
  ierr = VecAXPY(fold->w,sigma,x);CHKERRQ(ierr);
  ierr = MatMult(fold->A,fold->w,y);CHKERRQ(ierr);
  ierr = VecAXPY(y,sigma,fold->w);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   STDestroy_Fold - This routine destroys the shell ST context.

   Input Parameter:
.  fold - user-defined spectral transformation context
*/
PetscErrorCode STDestroy_Fold(FoldShellST *fold)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = MatDestroy(&fold->A);CHKERRQ(ierr);
  ierr = VecDestroy(&fold->w);CHKERRQ(ierr);
  ierr = PetscFree(fold);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   testset:
      args: -m 11 -eps_nev 4 -terse
      output_file: output/ex35_1.out
      test:
         suffix: 1
         requires: !single
      test:
         suffix: 1_single
         args: -eps_tol 1e-5
         requires: single

TEST*/
