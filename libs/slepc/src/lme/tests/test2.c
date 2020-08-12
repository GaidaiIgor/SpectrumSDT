/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test dense LME functions.\n\n";

#include <slepclme.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  LME            lme;
  Mat            A,B,C,X;
  PetscInt       i,j,n=10,k=2;
  PetscScalar    *As,*Bs,*Cs,*Xs;
  PetscViewer    viewer;
  PetscBool      verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Dense matrix equations, n=%D.\n",n);CHKERRQ(ierr);

  /* Create LME object */
  ierr = LMECreate(PETSC_COMM_WORLD,&lme);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Create matrices */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A,"A");CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&B);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)B,"B");CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,k,NULL,&C);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)C,"C");CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)X,"X");CHKERRQ(ierr);

  /* Fill A with an upper Hessenberg Toeplitz matrix */
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  for (i=0;i<n;i++) As[i+i*n]=3.0-(PetscReal)n/2;
  for (i=0;i<n-1;i++) As[i+1+i*n]=0.5;
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) As[i+(i+j)*n]=1.0;
  }
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);

  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(A,viewer);CHKERRQ(ierr);
  }

  /* Fill B with the 1-D Laplacian matrix */
  ierr = MatDenseGetArray(B,&Bs);CHKERRQ(ierr);
  for (i=0;i<n;i++) Bs[i+i*n]=2.0;
  for (i=0;i<n-1;i++) { Bs[i+1+i*n]=-1; Bs[i+(i+1)*n]=-1; }
  ierr = MatDenseRestoreArray(B,&Bs);CHKERRQ(ierr);

  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix B - - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(B,viewer);CHKERRQ(ierr);
  }

  /* Solve Lyapunov equation A*X+X*A'= -B */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Solving Lyapunov equation for B\n");CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  ierr = MatDenseGetArray(B,&Bs);CHKERRQ(ierr);
  ierr = MatDenseGetArray(X,&Xs);CHKERRQ(ierr);
  ierr = LMEDenseLyapunov(lme,n,As,n,Bs,n,Xs,n);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(B,&Bs);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(X,&Xs);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Solution X - - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(X,viewer);CHKERRQ(ierr);
  }

  /* Fill C with a full-rank nx2 matrix */
  ierr = MatDenseGetArray(C,&Cs);CHKERRQ(ierr);
  for (i=0;i<k;i++) Cs[i+i*n] = (i%2)? -1.0: 1.0;
  ierr = MatDenseRestoreArray(C,&Cs);CHKERRQ(ierr);

  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix C - - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(C,viewer);CHKERRQ(ierr);
  }

  /* Solve Lyapunov equation A*X+X*A'= -C*C' */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Solving Lyapunov equation for C (Cholesky)\n");CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  ierr = MatDenseGetArray(C,&Cs);CHKERRQ(ierr);
  ierr = MatDenseGetArray(X,&Xs);CHKERRQ(ierr);
  ierr = LMEDenseHessLyapunovChol(lme,n,As,n,2,Cs,n,Xs,n,NULL);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(C,&Cs);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(X,&Xs);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Solution X - - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(X,viewer);CHKERRQ(ierr);
  }

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = MatDestroy(&X);CHKERRQ(ierr);
  ierr = LMEDestroy(&lme);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      args: -info :lme
      requires: double
      filter: sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/1e-\\1/g"

TEST*/
