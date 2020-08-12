/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "SVD problem with user-defined stopping test.\n\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = matrix rows.\n"
  "  -n <n>, where <n> = matrix columns (defaults to m+2).\n\n";

#include <slepcsvd.h>
#include <petsctime.h>

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

/*
    Function for user-defined stopping test.

    Checks that the computing time has not exceeded the deadline.
*/
PetscErrorCode MyStoppingTest(SVD svd,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,SVDConvergedReason *reason,void *ctx)
{
  PetscErrorCode ierr;
  PetscLogDouble now,deadline = *(PetscLogDouble*)ctx;

  PetscFunctionBeginUser;
  /* check if usual termination conditions are met */
  ierr = SVDStoppingBasic(svd,its,max_it,nconv,nev,reason,NULL);CHKERRQ(ierr);
  if (*reason==SVD_CONVERGED_ITERATING) {
    /* check if deadline has expired */
    ierr = PetscTime(&now);CHKERRQ(ierr);
    if (now>deadline) *reason = SVD_CONVERGED_USER;
  }
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  Mat                A;
  SVD                svd;
  SVDConvergedReason reason;
  PetscInt           m=20,n,Istart,Iend,i,col[2],nconv;
  PetscReal          seconds=2.5;     /* maximum time allowed for computation */
  PetscLogDouble     deadline;        /* time to abort computation */
  PetscScalar        value[] = { 1, 2 };
  PetscBool          terse,flg;
  PetscViewer        viewer;
  PetscErrorCode     ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,&flg);CHKERRQ(ierr);
  if (!flg) n=m+2;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nRectangular bidiagonal matrix, m=%D n=%D\n\n",m,n);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-seconds",&seconds,NULL);CHKERRQ(ierr);
  deadline = seconds;
  ierr = PetscTimeAdd(&deadline);CHKERRQ(ierr);

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
  ierr = SVDSetWhichSingularTriplets(svd,SVD_SMALLEST);CHKERRQ(ierr);
  ierr = SVDSetTolerances(svd,PETSC_DEFAULT,1000);CHKERRQ(ierr);
  ierr = SVDSetType(svd,SVDTRLANCZOS);CHKERRQ(ierr);
  ierr = SVDSetStoppingTestFunction(svd,MyStoppingTest,&deadline,NULL);CHKERRQ(ierr);
  ierr = SVDSetFromOptions(svd);CHKERRQ(ierr);

  /* call the solver */
  ierr = SVDSolve(svd);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = SVDErrorView(svd,SVD_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = SVDGetConvergedReason(svd,&reason);CHKERRQ(ierr);
    if (reason!=SVD_CONVERGED_USER) {
      ierr = SVDReasonView(svd,viewer);CHKERRQ(ierr);
      ierr = SVDErrorView(svd,SVD_ERROR_RELATIVE,viewer);CHKERRQ(ierr);
    } else {
      ierr = SVDGetConverged(svd,&nconv);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"SVD solve finished with %D converged eigenpairs; reason=%s\n",nconv,SVDConvergedReasons[reason]);CHKERRQ(ierr);
    }
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  }
  ierr = SVDDestroy(&svd);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -m 750 -seconds 0.1 -svd_max_it 10000

TEST*/
