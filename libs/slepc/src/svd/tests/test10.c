/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Lanczos SVD. Also illustrates the use of SVDSetBV().\n\n"
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
  Mat            A;
  SVD            svd;
  PetscInt       m=20,n,Istart,Iend,i,k=6,col[2];
  PetscScalar    value[] = { 1, 2 };
  PetscBool      flg,oneside=PETSC_FALSE;
  const char     *prefix;
  BV             U,V;
  Vec            u,v;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,&flg);CHKERRQ(ierr);
  if (!flg) n=m+2;
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
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
  ierr = MatCreateVecs(A,&v,&u);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Create standalone BV objects to illustrate use of SVDSetBV()
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = BVCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)U,"U");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(U,u,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(U);CHKERRQ(ierr);
  ierr = BVCreate(PETSC_COMM_WORLD,&V);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)V,"V");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(V,v,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(V);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Compute singular values
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
  ierr = SVDSetBV(svd,V,U);CHKERRQ(ierr);
  ierr = SVDSetOptionsPrefix(svd,"check_");CHKERRQ(ierr);
  ierr = SVDAppendOptionsPrefix(svd,"myprefix_");CHKERRQ(ierr);
  ierr = SVDGetOptionsPrefix(svd,&prefix);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"SVD prefix is currently: %s\n\n",prefix);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)svd,"SVD_solver");CHKERRQ(ierr);

  ierr = SVDSetOperator(svd,A);CHKERRQ(ierr);
  ierr = SVDSetType(svd,SVDLANCZOS);CHKERRQ(ierr);
  ierr = SVDSetFromOptions(svd);CHKERRQ(ierr);

  ierr = PetscObjectTypeCompare((PetscObject)svd,SVDLANCZOS,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = SVDLanczosGetOneSide(svd,&oneside);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Running Lanczos %s\n\n",oneside?"(onesided)":"");CHKERRQ(ierr);
  }
  ierr = PetscObjectTypeCompare((PetscObject)svd,SVDTRLANCZOS,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = SVDTRLanczosGetOneSide(svd,&oneside);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Running thick-restart Lanczos %s\n\n",oneside?"(onesided)":"");CHKERRQ(ierr);
  }

  ierr = SVDSolve(svd);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SVDErrorView(svd,SVD_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = BVDestroy(&U);CHKERRQ(ierr);
  ierr = BVDestroy(&V);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = SVDDestroy(&svd);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -check_myprefix_svd_nsv 3
      requires: double
      test:
         suffix: 1
         args: -check_myprefix_svd_view_vectors ::ascii_info
      test:
         suffix: 2
         args: -check_myprefix_svd_type trlanczos -check_myprefix_svd_monitor -check_myprefix_svd_view_values ::ascii_matlab
         filter: sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      test:
         suffix: 3
         args: -m 22 -n 20

TEST*/
