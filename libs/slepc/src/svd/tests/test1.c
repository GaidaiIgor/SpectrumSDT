/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test the solution of a SVD without calling SVDSetFromOptions (based on ex8.c).\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = matrix dimension.\n"
  "  -type <svd_type> = svd type to test.\n\n";

#include <slepcsvd.h>

/*
   This example computes the singular values of an nxn Grcar matrix,
   which is a nonsymmetric Toeplitz matrix:

              |  1  1  1  1               |
              | -1  1  1  1  1            |
              |    -1  1  1  1  1         |
              |       .  .  .  .  .       |
          A = |          .  .  .  .  .    |
              |            -1  1  1  1  1 |
              |               -1  1  1  1 |
              |                  -1  1  1 |
              |                     -1  1 |

 */

int main(int argc,char **argv)
{
  Mat            A;               /* Grcar matrix */
  SVD            svd;             /* singular value solver context */
  PetscInt       N=30,Istart,Iend,i,col[5],nconv1,nconv2;
  PetscScalar    value[] = { -1, 1, 1, 1, 1 };
  PetscReal      sigma_1,sigma_n;
  char           svdtype[30] = "cross",epstype[30] = "";
  PetscBool      flg;
  EPS            eps;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-type",svdtype,30,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-epstype",epstype,30,&flg);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nEstimate the condition number of a Grcar matrix, n=%D",N);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    col[0]=i-1; col[1]=i; col[2]=i+1; col[3]=i+2; col[4]=i+3;
    if (i==0) {
      ierr = MatSetValues(A,1,&i,4,col+1,value+1,INSERT_VALUES);CHKERRQ(ierr);
    } else {
      ierr = MatSetValues(A,1,&i,PetscMin(5,N-i+1),col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Create the singular value solver and set the solution method
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create singular value context
  */
  ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);

  /*
     Set operator
  */
  ierr = SVDSetOperator(svd,A);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = SVDSetType(svd,svdtype);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscObjectTypeCompare((PetscObject)svd,SVDCROSS,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = SVDCrossGetEPS(svd,&eps);CHKERRQ(ierr);
      ierr = EPSSetType(eps,epstype);CHKERRQ(ierr);
    }
    ierr = PetscObjectTypeCompare((PetscObject)svd,SVDCYCLIC,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = SVDCyclicGetEPS(svd,&eps);CHKERRQ(ierr);
      ierr = EPSSetType(eps,epstype);CHKERRQ(ierr);
    }
  }
  ierr = SVDSetDimensions(svd,1,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = SVDSetTolerances(svd,1e-6,1000);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Compute the singular values
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     First request the largest singular value
  */
  ierr = SVDSetWhichSingularTriplets(svd,SVD_LARGEST);CHKERRQ(ierr);
  ierr = SVDSolve(svd);CHKERRQ(ierr);
  /*
     Get number of converged singular values
  */
  ierr = SVDGetConverged(svd,&nconv1);CHKERRQ(ierr);
  /*
     Get converged singular values: largest singular value is stored in sigma_1.
     In this example, we are not interested in the singular vectors
  */
  if (nconv1 > 0) {
    ierr = SVDGetSingularTriplet(svd,0,&sigma_1,NULL,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD," Unable to compute large singular value!\n\n");CHKERRQ(ierr);
  }

  /*
     Request the smallest singular value
  */
  ierr = SVDSetWhichSingularTriplets(svd,SVD_SMALLEST);CHKERRQ(ierr);
  ierr = SVDSolve(svd);CHKERRQ(ierr);
  /*
     Get number of converged triplets
  */
  ierr = SVDGetConverged(svd,&nconv2);CHKERRQ(ierr);
  /*
     Get converged singular values: smallest singular value is stored in sigma_n.
     As before, we are not interested in the singular vectors
  */
  if (nconv2 > 0) {
    ierr = SVDGetSingularTriplet(svd,0,&sigma_n,NULL,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD," Unable to compute small singular value!\n\n");CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (nconv1 > 0 && nconv2 > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD," Computed singular values: sigma_1=%.4f, sigma_n=%.4f\n",(double)sigma_1,(double)sigma_n);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Estimated condition number: sigma_1/sigma_n=%.4f\n\n",(double)(sigma_1/sigma_n));CHKERRQ(ierr);
  }

  /*
     Free work space
  */
  ierr = SVDDestroy(&svd);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -type {{lanczos trlanczos cross cyclic lapack}}

   test:
      suffix: 1_cross_gd
      args: -type cross -epstype gd
      output_file: output/test1_1.out

   test:
      suffix: 1_cyclic_gd
      args: -type cyclic -epstype gd
      output_file: output/test1_1.out

   test:
      suffix: 1_primme
      args: -type primme
      requires: primme
      output_file: output/test1_1.out

TEST*/
