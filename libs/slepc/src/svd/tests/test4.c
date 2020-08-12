/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test an SVD problem with more columns than rows.\n\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = matrix rows.\n"
  "  -n <n>, where <n> = matrix columns (defaults to m+2).\n\n";

#include <slepcsvd.h>

/*
   This example computes the singular values of a rectangular bidiagonal matrix

              |  1  2                     |
              |     1  2                  |
              |        1  2               |
          A = |          .  .             |
              |             .  .          |
              |                1  2       |
              |                   1  2    |
 */

int main(int argc,char **argv)
{
  Mat                  A,B;
  SVD                  svd;
  SVDConv              conv;
  SVDStop              stop;
  SVDWhich             which;
  SVDConvergedReason   reason;
  PetscInt             m=20,n,Istart,Iend,i,col[2],its;
  PetscScalar          value[] = { 1, 2 };
  PetscBool            flg,tmode;
  PetscErrorCode       ierr;
  PetscViewerAndFormat *vf;
  const char           *ctest[] = { "absolute", "relative to the singular value", "user-defined" };
  const char           *stest[] = { "basic", "user-defined" };

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,&flg);CHKERRQ(ierr);
  if (!flg) n=m+2;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nRectangular bidiagonal matrix, m=%D n=%D\n\n",m,n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    col[0]=i; col[1]=i+1;
    if (i<n-1) {
      ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    } else if (i==n-1) {
      ierr = MatSetValue(A,i,col[0],value[0],INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Compute singular values
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
  ierr = SVDSetOperator(svd,A);CHKERRQ(ierr);

  /* test some interface functions */
  ierr = SVDGetOperator(svd,&B);CHKERRQ(ierr);
  ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = SVDSetConvergenceTest(svd,SVD_CONV_ABS);CHKERRQ(ierr);
  ierr = SVDSetStoppingTest(svd,SVD_STOP_BASIC);CHKERRQ(ierr);
  /* test monitors */
  ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);CHKERRQ(ierr);
  ierr = SVDMonitorSet(svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))SVDMonitorFirst,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);
  /* ierr = SVDMonitorCancel(svd);CHKERRQ(ierr); */
  ierr = SVDSetFromOptions(svd);CHKERRQ(ierr);

  /* query properties and print them */
  ierr = SVDGetImplicitTranspose(svd,&tmode);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Transpose mode is %s\n",tmode?"implicit":"explicit");CHKERRQ(ierr);
  ierr = SVDGetConvergenceTest(svd,&conv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Convergence test is %s\n",ctest[conv]);CHKERRQ(ierr);
  ierr = SVDGetStoppingTest(svd,&stop);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping test is %s\n",stest[stop]);CHKERRQ(ierr);
  ierr = SVDGetWhichSingularTriplets(svd,&which);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Which = %s\n",which?"largest":"smallest");CHKERRQ(ierr);

  /* call the solver */
  ierr = SVDSolve(svd);CHKERRQ(ierr);
  ierr = SVDGetConvergedReason(svd,&reason);CHKERRQ(ierr);
  ierr = SVDGetIterationNumber(svd,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Finished - converged reason = %d\n",(int)reason);CHKERRQ(ierr);
  /* ierr = PetscPrintf(PETSC_COMM_WORLD," its = %D\n",its);CHKERRQ(ierr); */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SVDErrorView(svd,SVD_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = SVDDestroy(&svd);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -svd_monitor_cancel
      filter: grep -v "Transpose mode"
      output_file: output/test4_1.out
      test:
         suffix: 1_lanczos
         args: -svd_type lanczos
      test:
         suffix: 1_trlanczos
         args: -svd_type trlanczos -svd_ncv 12
      test:
         suffix: 1_cross
         args: -svd_type cross
      test:
         suffix: 1_cross_exp
         args: -svd_type cross -svd_cross_explicitmatrix
      test:
         suffix: 1_cross_exp_imp
         args: -svd_type cross -svd_cross_explicitmatrix -svd_implicittranspose
         requires: !complex
      test:
         suffix: 1_cyclic
         args: -svd_type cyclic
      test:
         suffix: 1_cyclic_imp
         args: -svd_type cyclic -svd_implicittranspose
      test:
         suffix: 1_cyclic_exp
         args: -svd_type cyclic -svd_cyclic_explicitmatrix
      test:
         suffix: 1_lapack
         args: -svd_type lapack
         requires: !single

   testset:
      args: -svd_monitor_cancel  -mat_type aijcusparse
      requires: cuda !single
      filter: grep -v "Transpose mode" | sed -e "s/seqaijcusparse/seqaij/"
      output_file: output/test4_1.out
      test:
         suffix: 2_cuda_lanczos
         args: -svd_type lanczos
      test:
         suffix: 2_cuda_trlanczos
         args: -svd_type trlanczos -svd_ncv 12
      test:
         suffix: 2_cuda_cross
         args: -svd_type cross

   test:
      suffix: 3
      nsize: 2
      args: -svd_type trlanczos -svd_ncv 14 -svd_monitor_cancel -ds_parallel synchronized

TEST*/
