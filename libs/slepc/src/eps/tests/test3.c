/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests multiple calls to EPSSolve with different matrix.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A1,A2;       /* problem matrices */
  EPS            eps;         /* eigenproblem solver context */
  PetscReal      tol=1000*PETSC_MACHINE_EPSILON,v;
  Vec            d;
  PetscInt       n=30,i,Istart,Iend;
  PetscRandom    myrand;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nTridiagonal with random diagonal, n=%D\n\n",n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           Create matrix tridiag([-1 0 -1])
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatCreate(PETSC_COMM_WORLD,&A1);CHKERRQ(ierr);
  ierr = MatSetSizes(A1,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A1);CHKERRQ(ierr);
  ierr = MatSetUp(A1);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A1,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(A1,i,i-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A1,i,i+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
  }
  ierr = MatAssemblyBegin(A1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create two matrices by filling the diagonal with rand values
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatDuplicate(A1,MAT_COPY_VALUES,&A2);CHKERRQ(ierr);
  ierr = MatCreateVecs(A1,NULL,&d);CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&myrand);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(myrand);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(myrand,0.0,1.0);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = PetscRandomGetValueReal(myrand,&v);CHKERRQ(ierr);
    ierr = VecSetValue(d,i,v,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(d);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(d);CHKERRQ(ierr);
  ierr = MatDiagonalSet(A1,d,INSERT_VALUES);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = PetscRandomGetValueReal(myrand,&v);CHKERRQ(ierr);
    ierr = VecSetValue(d,i,v,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(d);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(d);CHKERRQ(ierr);
  ierr = MatDiagonalSet(A2,d,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecDestroy(&d);CHKERRQ(ierr);
  ierr = PetscRandomDestroy(&myrand);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Create the eigensolver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,tol,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A1,NULL);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve first eigenproblem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," - - - First matrix - - -\n");CHKERRQ(ierr);
  ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve second eigenproblem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSetOperators(eps,A2,NULL);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," - - - Second matrix - - -\n");CHKERRQ(ierr);
  ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A1);CHKERRQ(ierr);
  ierr = MatDestroy(&A2);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -eps_nev 4
      requires: !single
      output_file: output/test3_1.out
      test:
         suffix: 1
         args: -eps_type {{krylovschur subspace arnoldi lanczos lapack}}
      test:
         suffix: 1_power
         args: -eps_type power -eps_max_it 20000
      test:
         suffix: 1_jd
         args: -eps_type jd -eps_jd_initial_size 7
      test:
         suffix: 1_gd
         args: -eps_type gd -eps_gd_initial_size 7
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion
      test:
         suffix: 1_arpack
         args: -eps_type arpack
         requires: arpack
      test:
         suffix: 1_primme
         args: -eps_type primme -eps_conv_abs -eps_primme_blocksize 4
         requires: primme
      test:
         suffix: 1_trlan
         args: -eps_type trlan
         requires: trlan
      test:
         suffix: 1_blzpack
         args: -eps_type blzpack -st_type sinvert -eps_target 2.7
         requires: blzpack

   testset:
      args: -eps_nev 4 -eps_smallest_real -eps_max_it 500
      requires: !single
      output_file: output/test3_2.out
      test:
         suffix: 2_rqcg
         args: -eps_type rqcg -eps_rqcg_reset 5 -eps_ncv 32
      test:
         suffix: 2_lobpcg
         args: -eps_type lobpcg -eps_lobpcg_blocksize 6 -st_pc_type none
      test:
         suffix: 2_lanczos
         args: -eps_type lanczos
      test:
         suffix: 2_lanczos_delayed
         args: -eps_type lanczos -eps_lanczos_reorthog delayed -eps_tol 1e-8
      test:
         suffix: 2_trlan
         args: -eps_type trlan
         requires: trlan
      test:
         suffix: 2_blopex
         args: -eps_type blopex -eps_conv_abs -st_shift -2
         requires: blopex

TEST*/
