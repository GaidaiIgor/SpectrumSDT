/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Example based on spring problem in NLEVP collection [1]. See the parameters
   meaning at Example 2 in [2].

   [1] T. Betcke, N. J. Higham, V. Mehrmann, C. Schroder, and F. Tisseur,
       NLEVP: A Collection of Nonlinear Eigenvalue Problems, MIMS EPrint
       2010.98, November 2010.
   [2] F. Tisseur, Backward error and condition of polynomial eigenvalue
       problems, Linear Algebra and its Applications, 309 (2000), pp. 339--361,
       April 2000.
*/

static char help[] = "Test the solution of a PEP from a finite element model of "
  "damped mass-spring system (problem from NLEVP collection).\n\n"
  "The command line options are:\n"
  "  -n <n> ... number of grid subdivisions.\n"
  "  -mu <value> ... mass (default 1).\n"
  "  -tau <value> ... damping constant of the dampers (default 10).\n"
  "  -kappa <value> ... damping constant of the springs (default 5).\n"
  "  -initv ... set an initial vector.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            M,C,K,A[3];      /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscErrorCode ierr;
  PetscInt       n=30,Istart,Iend,i,nev;
  PetscScalar    mu=1.0,tau=10.0,kappa=5.0;
  PetscBool      initv=PETSC_FALSE;
  Vec            IV[2];

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-mu",&mu,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-tau",&tau,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-kappa",&kappa,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-initv",&initv,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is a tridiagonal */
  ierr = MatCreate(PETSC_COMM_WORLD,&K);CHKERRQ(ierr);
  ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(K);CHKERRQ(ierr);
  ierr = MatSetUp(K);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(K,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) {
      ierr = MatSetValue(K,i,i-1,-kappa,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatSetValue(K,i,i,kappa*3.0,INSERT_VALUES);CHKERRQ(ierr);
    if (i<n-1) {
      ierr = MatSetValue(K,i,i+1,-kappa,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* C is a tridiagonal */
  ierr = MatCreate(PETSC_COMM_WORLD,&C);CHKERRQ(ierr);
  ierr = MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(C);CHKERRQ(ierr);
  ierr = MatSetUp(C);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(C,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) {
      ierr = MatSetValue(C,i,i-1,-tau,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatSetValue(C,i,i,tau*3.0,INSERT_VALUES);CHKERRQ(ierr);
    if (i<n-1) {
      ierr = MatSetValue(C,i,i+1,-tau,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* M is a diagonal matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&M);CHKERRQ(ierr);
  ierr = MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(M);CHKERRQ(ierr);
  ierr = MatSetUp(M);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(M,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(M,i,i,mu,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  A[0] = K; A[1] = C; A[2] = M;
  ierr = PEPSetOperators(pep,3,A);CHKERRQ(ierr);
  ierr = PEPSetProblemType(pep,PEP_GENERAL);CHKERRQ(ierr);
  ierr = PEPSetTolerances(pep,PETSC_SMALL,PETSC_DEFAULT);CHKERRQ(ierr);
  if (initv) { /* initial vector */
    ierr = MatCreateVecs(K,&IV[0],NULL);CHKERRQ(ierr);
    ierr = VecSetValue(IV[0],0,-1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(IV[0],1,0.5,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(IV[0]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(IV[0]);CHKERRQ(ierr);
    ierr = MatCreateVecs(K,&IV[1],NULL);CHKERRQ(ierr);
    ierr = VecSetValue(IV[1],0,4.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(IV[1],2,1.5,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(IV[1]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(IV[1]);CHKERRQ(ierr);
    ierr = PEPSetInitialSpace(pep,2,IV);CHKERRQ(ierr);
    ierr = VecDestroy(&IV[0]);CHKERRQ(ierr);
    ierr = VecDestroy(&IV[1]);CHKERRQ(ierr);
  }
  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPSolve(pep);CHKERRQ(ierr);
  ierr = PEPGetDimensions(pep,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL);CHKERRQ(ierr);
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -pep_nev 4 -initv
      requires: !single
      output_file: output/test2_1.out
      test:
         suffix: 1
         args: -pep_type {{toar linear}}
      test:
         suffix: 1_toar_mgs
         args: -pep_type toar -bv_orthog_type mgs
      test:
         suffix: 1_qarnoldi
         args: -pep_type qarnoldi -bv_orthog_refine never
      test:
         suffix: 1_linear_gd
         args: -pep_type linear -pep_linear_eps_type gd -pep_linear_explicitmatrix

   testset:
      args: -pep_target -0.43 -pep_nev 4 -pep_ncv 20 -st_type sinvert
      output_file: output/test2_2.out
      test:
         suffix: 2
         args: -pep_type {{toar linear}}
      test:
         suffix: 2_toar_scaleboth
         args: -pep_type toar -pep_scale both
      test:
         suffix: 2_toar_transform
         args: -pep_type toar -st_transform
      test:
         suffix: 2_qarnoldi
         args: -pep_type qarnoldi -bv_orthog_refine always
      test:
         suffix: 2_linear_explicit
         args: -pep_type linear -pep_linear_explicitmatrix -pep_linear_linearization 0,1
      test:
         suffix: 2_linear_explicit_her
         args: -pep_type linear -pep_linear_explicitmatrix -pep_hermitian -pep_linear_linearization 0,1
      test:
         suffix: 2_stoar
         args: -pep_type stoar -pep_hermitian
         requires: !single
      test:
         suffix: 2_jd
         args: -pep_type jd -st_type precond -pep_max_it 200
         requires: !single

   test:
      suffix: 3
      args: -pep_nev 12 -pep_extract {{none norm residual structured}} -pep_monitor_cancel
      requires: !single

   testset:
      args: -st_type sinvert -pep_target -0.43 -pep_nev 4
      output_file: output/test2_2.out
      test:
         suffix: 4_schur
         args: -pep_refine simple -pep_refine_scheme schur
      test:
         suffix: 4_mbe
         args: -pep_refine simple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
      test:
         suffix: 4_explicit
         args: -pep_refine simple -pep_refine_scheme explicit
      test:
         suffix: 4_multiple_schur
         args: -pep_refine multiple -pep_refine_scheme schur
      test:
         suffix: 4_multiple_mbe
         args: -pep_refine multiple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu -pep_refine_pc_factor_shift_type nonzero
      test:
         suffix: 4_multiple_explicit
         args: -pep_refine multiple -pep_refine_scheme explicit

   test:
      suffix: 5
      nsize: 2
      args: -pep_type linear -pep_linear_explicitmatrix -pep_general -pep_target -0.43 -pep_nev 4 -pep_ncv 20 -st_type sinvert -pep_linear_st_ksp_type bcgs -pep_linear_st_pc_type bjacobi
      output_file: output/test2_2.out

   test:
      suffix: 6
      args: -pep_type linear -pep_nev 12 -pep_extract {{none norm residual}}
      requires: !single
      output_file: output/test2_3.out

   test:
      suffix: 7
      args: -pep_nev 12 -pep_extract {{none norm residual structured}} -pep_refine multiple
      requires: !single
      output_file: output/test2_3.out

   testset:
      args: -st_type sinvert -pep_target -0.43 -pep_nev 4 -st_transform
      output_file: output/test2_2.out
      test:
         suffix: 8_schur
         args: -pep_refine simple -pep_refine_scheme schur
      test:
         suffix: 8_mbe
         args: -pep_refine simple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
      test:
         suffix: 8_explicit
         args: -pep_refine simple -pep_refine_scheme explicit
      test:
         suffix: 8_multiple_schur
         args: -pep_refine multiple -pep_refine_scheme schur
      test:
         suffix: 8_multiple_mbe
         args: -pep_refine multiple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
      test:
         suffix: 8_multiple_explicit
         args: -pep_refine multiple -pep_refine_scheme explicit

   testset:
      nsize: 2
      args: -st_type sinvert -pep_target -0.49 -pep_nev 4 -pep_refine_partitions 2 -st_ksp_type bcgs -st_pc_type bjacobi -pep_scale diagonal -pep_scale_its 4
      output_file: output/test2_2.out
      requires: !single
      test:
         suffix: 9_mbe
         args: -pep_refine simple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
      test:
         suffix: 9_explicit
         args: -pep_refine simple -pep_refine_scheme explicit
      test:
         suffix: 9_multiple_mbe
         args: -pep_refine multiple -pep_refine_scheme mbe -pep_refine_ksp_type preonly -pep_refine_pc_type lu
      test:
         suffix: 9_multiple_explicit
         args: -pep_refine multiple -pep_refine_scheme explicit

   test:
      suffix: 10
      nsize: 4
      args: -st_type sinvert -pep_target -0.43 -pep_nev 4 -pep_refine simple -pep_refine_scheme explicit -pep_refine_partitions 2 -st_ksp_type bcgs -st_pc_type bjacobi -pep_scale diagonal -pep_scale_its 4
      requires: double
      output_file: output/test2_2.out

   testset:
      args: -pep_nev 4 -initv -mat_type aijcusparse
      output_file: output/test2_1.out
      requires: cuda !single
      test:
         suffix: 11_cuda
         args: -pep_type {{toar linear}}
      test:
         suffix: 11_cuda_qarnoldi
         args: -pep_type qarnoldi -bv_orthog_refine never
      test:
         suffix: 11_cuda_linear_gd
         args: -pep_type linear -pep_linear_eps_type gd -pep_linear_explicitmatrix

   test:
      suffix: 12
      nsize: 2
      args: -pep_type jd -ds_parallel synchronized -pep_target -0.43 -pep_nev 4 -pep_ncv 20
      requires: !single

   test:
      suffix: 13
      args: -pep_nev 12 -pep_view_values draw -pep_monitor_lg
      requires: x !single
      output_file: output/test2_3.out

TEST*/
