/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test SVD with user-provided initial vectors.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = row dimension.\n"
  "  -m <m>, where <m> = column dimension.\n\n";

#include <slepcsvd.h>

/*
   This example computes the singular values of a rectangular nxm Grcar matrix:

              |  1  1  1  1               |
              | -1  1  1  1  1            |
              |    -1  1  1  1  1         |
          A = |       .  .  .  .  .       |
              |          .  .  .  .  .    |
              |            -1  1  1  1  1 |
              |               -1  1  1  1 |

 */

int main(int argc,char **argv)
{
  Mat            A;               /* Grcar matrix */
  SVD            svd;             /* singular value solver context */
  Vec            v0,w0;           /* initial vectors */
  PetscInt       N=35,M=30,Istart,Iend,i,col[5];
  PetscScalar    value[] = { -1, 1, 1, 1, 1 };
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&M,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSVD of a rectangular Grcar matrix, %Dx%D\n\n",N,M);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,M);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    col[0]=i-1; col[1]=i; col[2]=i+1; col[3]=i+2; col[4]=i+3;
    if (i==0) {
      ierr = MatSetValues(A,1,&i,PetscMin(4,M-i+1),col+1,value+1,INSERT_VALUES);CHKERRQ(ierr);
    } else {
      ierr = MatSetValues(A,1,&i,PetscMin(5,M-i+1),col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Create the SVD context and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
  ierr = SVDSetOperator(svd,A);CHKERRQ(ierr);
  ierr = SVDSetFromOptions(svd);CHKERRQ(ierr);

  /*
     Set the initial vectors. This is optional, if not done the initial
     vectors are set to random values
  */
  ierr = MatCreateVecs(A,&v0,&w0);CHKERRQ(ierr);
  ierr = VecSet(v0,1.0);CHKERRQ(ierr);
  ierr = VecSet(w0,1.0);CHKERRQ(ierr);
  ierr = SVDSetInitialSpaces(svd,1,&v0,1,&w0);CHKERRQ(ierr);

  /*
     Compute solution
  */
  ierr = SVDSolve(svd);CHKERRQ(ierr);
  ierr = SVDErrorView(svd,SVD_ERROR_RELATIVE,NULL);CHKERRQ(ierr);

  /*
     Free work space
  */
  ierr = VecDestroy(&v0);CHKERRQ(ierr);
  ierr = VecDestroy(&w0);CHKERRQ(ierr);
  ierr = SVDDestroy(&svd);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -svd_nsv 4
      requires: !single
      output_file: output/test3_1.out
      test:
         suffix: 1_lanczos
         args: -svd_type lanczos
      test:
         suffix: 1_lanczos_one
         args: -svd_type lanczos -svd_lanczos_oneside
      test:
         suffix: 1_trlanczos
         args: -svd_type trlanczos
      test:
         suffix: 1_trlanczos_one
         args: -svd_type trlanczos -svd_trlanczos_oneside
      test:
         suffix: 1_trlanczos_one_mgs
         args: -svd_type trlanczos -svd_trlanczos_oneside -bv_orthog_type mgs
      test:
         suffix: 1_trlanczos_one_always
         args: -svd_type trlanczos -svd_trlanczos_oneside -bv_orthog_refine always
      test:
         suffix: 1_cross
         args: -svd_type cross
      test:
         suffix: 1_cross_exp
         args: -svd_type cross -svd_cross_explicitmatrix
      test:
         suffix: 1_cyclic
         args: -svd_type cyclic
      test:
         suffix: 1_cyclic_exp
         args: -svd_type cyclic -svd_cyclic_explicitmatrix
      test:
         suffix: 1_lapack
         args: -svd_type lapack
      test:
         suffix: 1_primme
         args: -svd_type primme
         requires: primme

   testset:
      args: -svd_implicittranspose -svd_nsv 4 -svd_tol 1e-5
      requires: !single
      output_file: output/test3_1.out
      test:
         suffix: 2_lanczos
         args: -svd_type lanczos
      test:
         suffix: 2_lanczos_one
         args: -svd_type lanczos -svd_lanczos_oneside
      test:
         suffix: 2_trlanczos
         args: -svd_type trlanczos
      test:
         suffix: 2_trlanczos_one
         args: -svd_type trlanczos -svd_trlanczos_oneside
      test:
         suffix: 2_trlanczos_one_mgs
         args: -svd_type trlanczos -svd_trlanczos_oneside -bv_orthog_type mgs
      test:
         suffix: 2_trlanczos_one_always
         args: -svd_type trlanczos -svd_trlanczos_oneside -bv_orthog_refine always
      test:
         suffix: 2_cross
         args: -svd_type cross
      test:
         suffix: 2_cross_exp
         args: -svd_type cross -svd_cross_explicitmatrix
         requires: !complex
      test:
         suffix: 2_cyclic
         args: -svd_type cyclic
      test:
         suffix: 2_lapack
         args: -svd_type lapack

   testset:
      args: -svd_nsv 4 -mat_type aijcusparse
      requires: cuda !single
      output_file: output/test3_1.out
      test:
         suffix: 3_cuda_lanczos
         args: -svd_type lanczos
      test:
         suffix: 3_cuda_lanczos_one
         args: -svd_type lanczos -svd_lanczos_oneside
      test:
         suffix: 3_cuda_trlanczos
         args: -svd_type trlanczos
      test:
         suffix: 3_cuda_trlanczos_one
         args: -svd_type trlanczos -svd_trlanczos_oneside
      test:
         suffix: 3_cuda_trlanczos_one_mgs
         args: -svd_type trlanczos -svd_trlanczos_oneside -bv_orthog_type mgs
      test:
         suffix: 3_cuda_trlanczos_one_always
         args: -svd_type trlanczos -svd_trlanczos_oneside -bv_orthog_refine always
      test:
         suffix: 3_cuda_cross
         args: -svd_type cross
      test:
         suffix: 3_cuda_cyclic
         args: -svd_type cyclic
      test:
         suffix: 3_cuda_cyclic_exp
         args: -svd_type cyclic -svd_cyclic_explicitmatrix

   test:
      suffix: 4
      requires: !single
      args: -svd_type lapack -svd_nsv 4
      output_file: output/test3_1.out
      nsize: 2

   test:
      suffix: 5
      args: -svd_nsv 4 -svd_view_values draw -svd_monitor_lg
      requires: x !single
      output_file: output/test3_1.out

TEST*/
