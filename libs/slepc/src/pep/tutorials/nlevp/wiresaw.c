/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This example implements two of the problems found at
       NLEVP: A Collection of Nonlinear Eigenvalue Problems,
       The University of Manchester.
   The details of the collection can be found at:
       [1] T. Betcke et al., "NLEVP: A Collection of Nonlinear Eigenvalue
           Problems", ACM Trans. Math. Software 39(2), Article 7, 2013.

   WIRESAW1 is a gyroscopic QEP from vibration analysis of a wiresaw,
   where the parameter V represents the speed of the wire. When the
   parameter eta is nonzero, then it turns into the WIRESAW2 problem
   (with added viscous damping, e.g. eta=0.8).

       [2] S. Wei and I. Kao, "Vibration analysis of wire and frequency
           response in the modern wiresaw manufacturing process", J. Sound
           Vib. 213(5):1383-1395, 2000.
*/

static char help[] = "Vibration analysis of a wiresaw.\n\n"
  "The command line options are:\n"
  "  -n <n> ... dimension of the matrices (default 10).\n"
  "  -v <value> ... velocity of the wire (default 0.01).\n"
  "  -eta <value> ... viscous damping (default 0.0).\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            M,D,K,A[3];      /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscInt       n=10,Istart,Iend,j,k;
  PetscReal      v=0.01,eta=0.0;
  PetscBool      terse;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-v",&v,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-eta",&eta,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nVibration analysis of a wiresaw, n=%D v=%g eta=%g\n\n",n,(double)v,(double)eta);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*D+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is a diagonal matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&K);CHKERRQ(ierr);
  ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(K);CHKERRQ(ierr);
  ierr = MatSetUp(K);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(K,&Istart,&Iend);CHKERRQ(ierr);
  for (j=Istart;j<Iend;j++) {
    ierr = MatSetValue(K,j,j,(j+1)*(j+1)*PETSC_PI*PETSC_PI*(1.0-v*v),INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatScale(K,0.5);CHKERRQ(ierr);

  /* D is a tridiagonal */
  ierr = MatCreate(PETSC_COMM_WORLD,&D);CHKERRQ(ierr);
  ierr = MatSetSizes(D,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(D);CHKERRQ(ierr);
  ierr = MatSetUp(D);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(D,&Istart,&Iend);CHKERRQ(ierr);
  for (j=Istart;j<Iend;j++) {
    for (k=0;k<n;k++) {
      if ((j+k)%2) {
        ierr = MatSetValue(D,j,k,8.0*(j+1)*(k+1)*v/((j+1)*(j+1)-(k+1)*(k+1)),INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  ierr = MatAssemblyBegin(D,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(D,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatScale(D,0.5);CHKERRQ(ierr);

  /* M is a diagonal matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&M);CHKERRQ(ierr);
  ierr = MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(M);CHKERRQ(ierr);
  ierr = MatSetUp(M);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(M,&Istart,&Iend);CHKERRQ(ierr);
  for (j=Istart;j<Iend;j++) {
    ierr = MatSetValue(M,j,j,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatScale(M,0.5);CHKERRQ(ierr);

  /* add damping */
  if (eta>0.0) {
    ierr = MatAXPY(K,eta,D,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr); /* K = K + eta*D */
    ierr = MatShift(D,eta);CHKERRQ(ierr); /* D = D + eta*eye(n) */
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  A[0] = K; A[1] = D; A[2] = M;
  ierr = PEPSetOperators(pep,3,A);CHKERRQ(ierr);
  if (eta==0.0) {
    ierr = PEPSetProblemType(pep,PEP_GYROSCOPIC);CHKERRQ(ierr);
  } else {
    ierr = PEPSetProblemType(pep,PEP_GENERAL);CHKERRQ(ierr);
  }
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
  ierr = MatDestroy(&D);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -pep_nev 4 -terse
      requires: !single
      output_file: output/wiresaw_1.out
      test:
         suffix: 1
         args: -pep_type {{toar qarnoldi}}
      test:
         suffix: 1_linear_h1
         args: -pep_type linear -pep_linear_explicitmatrix -pep_linear_linearization 1,0 -pep_linear_st_ksp_type bcgs -pep_linear_st_pc_type kaczmarz
      test:
         suffix: 1_linear_h2
         args: -pep_type linear -pep_linear_explicitmatrix -pep_linear_linearization 0,1

TEST*/
