/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This example implements one of the problems found at
       NLEVP: A Collection of Nonlinear Eigenvalue Problems,
       The University of Manchester.
   The details of the collection can be found at:
       [1] T. Betcke et al., "NLEVP: A Collection of Nonlinear Eigenvalue
           Problems", ACM Trans. Math. Software 39(2), Article 7, 2013.

   The sleeper problem is a proportionally damped QEP describing the
   oscillations of a rail track resting on sleepers.
*/

static char help[] = "Oscillations of a rail track resting on sleepers.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = dimension of the matrices.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            M,C,K,A[3];      /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscInt       n=10,Istart,Iend,i;
  PetscBool      terse;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nRailtrack resting on sleepers, n=%D\n\n",n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is a pentadiagonal */
  ierr = MatCreate(PETSC_COMM_WORLD,&K);CHKERRQ(ierr);
  ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(K);CHKERRQ(ierr);
  ierr = MatSetUp(K);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(K,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i==0) {
      ierr = MatSetValue(K,i,n-1,-3.0,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatSetValue(K,i,n-2,1.0,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (i==1) { ierr = MatSetValue(K,i,n-1,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i>0) { ierr = MatSetValue(K,i,i-1,-3.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i>1) { ierr = MatSetValue(K,i,i-2,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(K,i,i,5.0,INSERT_VALUES);CHKERRQ(ierr);
    if (i==n-1) {
      ierr = MatSetValue(K,i,0,-3.0,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatSetValue(K,i,1,1.0,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (i==n-2) { ierr = MatSetValue(K,i,0,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(K,i,i+1,-3.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-2) { ierr = MatSetValue(K,i,i+2,1.0,INSERT_VALUES);CHKERRQ(ierr); }
  }

  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* C is a circulant matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&C);CHKERRQ(ierr);
  ierr = MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(C);CHKERRQ(ierr);
  ierr = MatSetUp(C);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(C,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i==0) {
      ierr = MatSetValue(C,i,n-1,-4.0,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatSetValue(C,i,n-2,1.0,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (i==1) { ierr = MatSetValue(C,i,n-1,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i>0) { ierr = MatSetValue(C,i,i-1,-4.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i>1) { ierr = MatSetValue(C,i,i-2,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(C,i,i,7.0,INSERT_VALUES);CHKERRQ(ierr);
    if (i==n-1) {
      ierr = MatSetValue(C,i,0,-4.0,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatSetValue(C,i,1,1.0,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (i==n-2) { ierr = MatSetValue(C,i,0,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(C,i,i+1,-4.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-2) { ierr = MatSetValue(C,i,i+2,1.0,INSERT_VALUES);CHKERRQ(ierr); }
  }

  ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* M is the identity matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&M);CHKERRQ(ierr);
  ierr = MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(M);CHKERRQ(ierr);
  ierr = MatSetUp(M);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(M,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(M,i,i,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  A[0] = K; A[1] = C; A[2] = M;
  ierr = PEPSetOperators(pep,3,A);CHKERRQ(ierr);
  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);
  ierr = PEPSolve(pep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = PEPReasonView(pep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -n 100 -pep_nev 4 -pep_ncv 24 -st_type sinvert -terse
      output_file: output/sleeper_1.out
      test:
         suffix: 1
         args: -pep_type {{toar linear}} -pep_ncv 20
      test:
         suffix: 1_qarnoldi
         args: -pep_type qarnoldi -pep_qarnoldi_restart 0.4

   testset:
      args: -n 24 -pep_nev 4 -pep_ncv 9 -pep_target -.62 -terse
      output_file: output/sleeper_2.out
      test:
         suffix: 2_toar
         args: -pep_type toar -pep_toar_restart .3 -st_type sinvert
         requires: !single
      test:
         suffix: 2_jd
         args: -pep_type jd -pep_jd_restart .3 -pep_jd_projection orthogonal

   test:
      suffix: 3
      args: -n 275 -pep_type stoar -pep_hermitian -st_type sinvert -pep_nev 2 -pep_target -.89 -terse
      requires: !single

   test:
      suffix: 4
      args: -n 270 -pep_type stoar -pep_hermitian -pep_interval -3,-2.51 -st_type sinvert -st_pc_type cholesky -terse
      requires: !single

TEST*/
