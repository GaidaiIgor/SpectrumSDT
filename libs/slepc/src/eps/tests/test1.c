/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests B-orthonormality of eigenvectors in a GHEP problem.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat               A,B;        /* matrices */
  EPS               eps;        /* eigenproblem solver context */
  ST                st;
  Vec               *X,v;
  PetscReal         lev,tol=1000*PETSC_MACHINE_EPSILON;
  PetscInt          N,n=45,m,Istart,Iend,II,i,j,nconv;
  PetscBool         flag;
  EPSPowerShiftType variant;
  PetscErrorCode    ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag);CHKERRQ(ierr);
  if (!flag) m=n;
  N = n*m;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized Symmetric Eigenproblem, N=%D (%Dx%D grid)\n\n",N,n,m);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(B);CHKERRQ(ierr);
  ierr = MatSetUp(B);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) { ierr = MatSetValue(A,II,II-n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<m-1) { ierr = MatSetValue(A,II,II+n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j>0) { ierr = MatSetValue(A,II,II-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(A,II,II+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,II,II,4.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(B,II,II,2.0/PetscLogScalar(II+2),INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateVecs(B,&v,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,tol,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = EPSSetConvergenceTest(eps,EPS_CONV_NORM);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* illustrate how to extract parameters from specific solver types */
  ierr = PetscObjectTypeCompare((PetscObject)eps,EPSPOWER,&flag);CHKERRQ(ierr);
  if (flag) {
    ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)st,STSHIFT,&flag);CHKERRQ(ierr);
    if (flag) {
      ierr = EPSPowerGetShiftType(eps,&variant);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Type of shifts used during power iteration: %s\n",EPSPowerShiftTypes[variant]);CHKERRQ(ierr);
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSGetTolerances(eps,&tol,NULL);CHKERRQ(ierr);
  ierr = EPSErrorView(eps,EPS_ERROR_BACKWARD,NULL);CHKERRQ(ierr);
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  if (nconv>1) {
    ierr = VecDuplicateVecs(v,nconv,&X);CHKERRQ(ierr);
    for (i=0;i<nconv;i++) {
      ierr = EPSGetEigenvector(eps,i,X[i],NULL);CHKERRQ(ierr);
    }
    ierr = VecCheckOrthogonality(X,nconv,NULL,nconv,B,NULL,&lev);CHKERRQ(ierr);
    if (lev<10*tol) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality below the tolerance\n");CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)lev);CHKERRQ(ierr);
    }
    ierr = VecDestroyVecs(nconv,&X);CHKERRQ(ierr);
  }

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -n 18 -eps_nev 4 -eps_max_it 1500
      requires: !single
      output_file: output/test1_1.out
      test:
         suffix: 1
         args: -eps_type {{krylovschur subspace arnoldi gd jd lapack}}
      test:
         suffix: 1_ks_nopurify
         args: -eps_purify 0
      test:
         suffix: 1_ks_trueres
         args: -eps_true_residual
      test:
         suffix: 1_ks_sinvert
         args: -st_type sinvert -eps_target 22
      test:
         suffix: 1_ks_cayley
         args: -st_type cayley -eps_target 22
      test:
         suffix: 1_lanczos
         args: -eps_type lanczos -eps_lanczos_reorthog full
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion
      test:
         suffix: 1_gd_borth
         args: -eps_type gd -eps_gd_borth
      test:
         suffix: 1_jd_borth
         args: -eps_type jd -eps_jd_borth
      test:
         suffix: 1_ciss
         args: -eps_type ciss -rg_interval_endpoints 20.8,22 -eps_largest_real
      test:
         suffix: 1_ciss_trapezoidal
         args: -eps_type ciss -rg_interval_endpoints 20.8,22 -eps_largest_real -eps_ciss_quadrule trapezoidal -eps_ciss_integration_points 24 -eps_ciss_extraction hankel -eps_ciss_delta 1e-10 -eps_tol 5e-11
      test:
         suffix: 1_lobpcg
         args: -eps_type lobpcg -st_shift 22 -eps_largest_real
      test:
         suffix: 1_cholesky
         args: -mat_type sbaij

   testset:
      requires: !single
      args: -eps_tol 1e-10 -st_type sinvert -st_ksp_type preonly -st_pc_type cholesky
      test:
         suffix: 2
         args: -eps_interval .1,1.1
      test:
         suffix: 2_open
         args: -eps_interval -inf,1.1
      test:
         suffix: 2_parallel
         requires: mumps
         nsize: 3
         args: -eps_interval .1,1.1 -eps_krylovschur_partitions 2 -st_pc_factor_mat_solver_type mumps -mat_mumps_icntl_13 1
         output_file: output/test1_2.out

   test:
      suffix: 3
      requires: !single
      args: -n 18 -eps_type power -eps_nev 3

   test:
      suffix: 4
      requires: !single
      args: -n 18 -eps_type power -eps_nev 3 -st_type sinvert -eps_target 1.149 -eps_power_shift_type {{constant rayleigh wilkinson}}

   testset:
      args: -n 18 -eps_nev 3 -eps_smallest_real -eps_max_it 500 -st_pc_type icc
      output_file: output/test1_5.out
      test:
         suffix: 5_rqcg
         args: -eps_type rqcg
      test:
         suffix: 5_lobpcg
         args: -eps_type lobpcg -eps_lobpcg_blocksize 3
         requires: !single
      test:
         suffix: 5_blopex
         args: -eps_type blopex -eps_conv_abs -st_shift 0.1
         requires: blopex

   testset:
      args: -n 18 -eps_nev 12 -eps_mpd 8 -eps_max_it 3000
      requires: !single
      output_file: output/test1_6.out
      test:
         suffix: 6
         args: -eps_type {{krylovschur subspace arnoldi gd}}
      test:
         suffix: 6_lanczos
         args: -eps_type lanczos -eps_lanczos_reorthog full

   testset:
      args: -n 18 -eps_nev 4 -eps_max_it 1500 -mat_type aijcusparse
      requires: cuda !single
      output_file: output/test1_1.out
      test:
         suffix: 7
         args: -eps_type {{krylovschur subspace arnoldi gd jd}}
      test:
         suffix: 7_ks_sinvert
         args: -st_type sinvert -eps_target 22
      test:
         suffix: 7_lanczos
         args: -eps_type lanczos -eps_lanczos_reorthog full
      test:
         suffix: 7_ciss
         args: -eps_type ciss -rg_interval_endpoints 20.8,22 -eps_largest_real

   testset:
      args: -n 18 -eps_nev 3 -eps_smallest_real -eps_max_it 500 -st_pc_type sor -mat_type aijcusparse
      requires: cuda
      output_file: output/test1_5.out
      test:
         suffix: 8_rqcg
         args: -eps_type rqcg
      test:
         suffix: 8_lobpcg
         args: -eps_type lobpcg -eps_lobpcg_blocksize 3

   testset:
      nsize: 2
      args: -n 18 -eps_nev 7 -eps_ncv 32 -ds_parallel synchronized
      requires: !single
      filter: grep -v "orthogonality"
      output_file: output/test1_9.out
      test:
         suffix: 9_ks_ghep
         args: -eps_gen_hermitian -st_pc_type redundant -st_type sinvert
      test:
         suffix: 9_ks_gnhep
         args: -eps_gen_non_hermitian -st_pc_type redundant -st_type sinvert
      test:
         suffix: 9_ks_ghiep
         args: -eps_gen_indefinite -st_pc_type redundant -st_type sinvert
      test:
         suffix: 9_lobpcg_ghep
         args: -eps_gen_hermitian -eps_type lobpcg -eps_max_it 200 -eps_lobpcg_blocksize 6
         timeoutfactor: 2
      test:
         suffix: 9_jd_gnhep
         args: -eps_gen_non_hermitian -eps_type jd -eps_target 0 -eps_ncv 64
         timeoutfactor: 2

TEST*/
