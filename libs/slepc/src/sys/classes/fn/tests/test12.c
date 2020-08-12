/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test matrix function evaluation via diagonalization.\n\n";

#include <slepcfn.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  FN             fn;
  Mat            A,F,G;
  PetscInt       i,j,n=10;
  PetscReal      nrm;
  PetscScalar    *As,alpha,beta;
  PetscViewer    viewer;
  PetscBool      verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix function of symmetric/Hermitian matrix, n=%D.\n",n);CHKERRQ(ierr);

  /* Create function object */
  ierr = FNCreate(PETSC_COMM_WORLD,&fn);CHKERRQ(ierr);
  ierr = FNSetType(fn,FNEXP);CHKERRQ(ierr);   /* default to exponential */
#if defined(PETSC_USE_COMPLEX)
  alpha = PetscCMPLX(0.3,0.8);
  beta  = PetscCMPLX(1.1,-0.1);
#else
  alpha = 0.3;
  beta  = 1.1;
#endif
  ierr = FNSetScale(fn,alpha,beta);CHKERRQ(ierr);
  ierr = FNSetFromOptions(fn);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Create a symmetric/Hermitian Toeplitz matrix */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A,"A");CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  for (i=0;i<n;i++) As[i+i*n]=2.0;
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) {
#if defined(PETSC_USE_COMPLEX)
      As[i+(i+j)*n]=PetscCMPLX(1.0,0.1); As[(i+j)+i*n]=PetscCMPLX(1.0,-0.1);
#else
      As[i+(i+j)*n]=0.5; As[(i+j)+i*n]=0.5;
#endif
    }
  }
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(A,viewer);CHKERRQ(ierr);
  }

  /* compute matrix function */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&F);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)F,"F");CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMat(fn,A,F);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed f(A) - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(F,viewer);CHKERRQ(ierr);
  }

  /* Repeat with MAT_HERMITIAN flag set */
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&G);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)G,"G");CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMat(fn,A,G);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed f(A) symm - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(G,viewer);CHKERRQ(ierr);
  }

  /* compare the two results */
  ierr = MatAXPY(F,-1.0,G,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatNorm(F,NORM_FROBENIUS,&nrm);CHKERRQ(ierr);
  if (nrm>100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of F-G is %g\n",(double)nrm);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed results match.\n");CHKERRQ(ierr);
  }

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&F);CHKERRQ(ierr);
  ierr = MatDestroy(&G);CHKERRQ(ierr);
  ierr = FNDestroy(&fn);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -fn_type {{exp sqrt}shared output}
      output_file: output/test12_1.out

   test:
      suffix: 1_rational
      nsize: 1
      args: -fn_type rational -fn_rational_numerator 2,-1.5 -fn_rational_denominator 1,0.8
      output_file: output/test12_1.out

TEST*/
