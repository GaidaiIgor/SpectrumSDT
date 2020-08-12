/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves a Lypunov equation with the shifted 2-D Laplacian.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepclme.h>

int main(int argc,char **argv)
{
  Mat                A;           /* problem matrix */
  Mat                C,C1;        /* right-hand side */
  Mat                X,X1;        /* solution */
  LME                lme;
  PetscReal          tol,errest,error;
  PetscScalar        *u,sigma=0.0;
  PetscInt           N,n=10,m,Istart,Iend,II,maxit,its,ncv,i,j,rank=0;
  PetscErrorCode     ierr;
  PetscBool          flag;
  LMEConvergedReason reason;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag);CHKERRQ(ierr);
  if (!flag) m=n;
  N = n*m;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-sigma",&sigma,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-rank",&rank,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nLyapunov equation, N=%D (%Dx%D grid)\n\n",N,n,m);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Create the 2-D Laplacian, A
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) { ierr = MatSetValue(A,II,II-n,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<m-1) { ierr = MatSetValue(A,II,II+n,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j>0) { ierr = MatSetValue(A,II,II-1,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(A,II,II+1,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,II,II,-4.0-sigma,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create a low-rank Mat to store the right-hand side C = C1*C1'
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&C1);CHKERRQ(ierr);
  ierr = MatSetSizes(C1,PETSC_DECIDE,PETSC_DECIDE,N,2);CHKERRQ(ierr);
  ierr = MatSetType(C1,MATDENSE);CHKERRQ(ierr);
  ierr = MatSetUp(C1);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(C1,&Istart,&Iend);CHKERRQ(ierr);
  ierr = MatDenseGetArray(C1,&u);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i<N/2) u[i-Istart] = 1.0;
    if (i==0) u[i+Iend-2*Istart] = -2.0;
    if (i==1) u[i+Iend-2*Istart] = -1.0;
    if (i==2) u[i+Iend-2*Istart] = -1.0;
  }
  ierr = MatDenseRestoreArray(C1,&u);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(C1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateLRC(NULL,C1,NULL,NULL,&C);CHKERRQ(ierr);
  ierr = MatDestroy(&C1);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create the matrix equation solver context
  */
  ierr = LMECreate(PETSC_COMM_WORLD,&lme);CHKERRQ(ierr);

  /*
     Set the type of equation
  */
  ierr = LMESetProblemType(lme,LME_LYAPUNOV);CHKERRQ(ierr);

  /*
     Set the matrix coefficients, the right-hand side, and the solution.
     In this case, it is a Lyapunov equation A*X+X*A'=-C where both
     C and X are symmetric and low-rank, C=C1*C1', X=X1*X1'
  */
  ierr = LMESetCoefficients(lme,A,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = LMESetRHS(lme,C);CHKERRQ(ierr);

  if (rank) {  /* Create X only if the user has specified a nonzero value of rank */
    ierr = PetscPrintf(PETSC_COMM_WORLD," Computing a solution with prescribed rank=%d\n",rank);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&X1);CHKERRQ(ierr);
    ierr = MatSetSizes(X1,PETSC_DECIDE,PETSC_DECIDE,N,rank);CHKERRQ(ierr);
    ierr = MatSetType(X1,MATDENSE);CHKERRQ(ierr);
    ierr = MatSetUp(X1);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(X1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(X1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatCreateLRC(NULL,X1,NULL,NULL,&X);CHKERRQ(ierr);
    ierr = MatDestroy(&X1);CHKERRQ(ierr);
    ierr = LMESetSolution(lme,X);CHKERRQ(ierr);
    ierr = MatDestroy(&X);CHKERRQ(ierr);
  }

  /*
     (Optional) Set other solver options
  */
  ierr = LMESetTolerances(lme,1e-07,PETSC_DEFAULT);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = LMESetFromOptions(lme);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                   Solve the matrix equation, A*X+X*A'=-C
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = LMESolve(lme);CHKERRQ(ierr);
  ierr = LMEGetConvergedReason(lme,&reason);CHKERRQ(ierr);
  if (reason<0) SETERRQ(PETSC_COMM_WORLD,1,"Solver did not converge");

  if (!rank) {  /* X1 was created by the solver, so extract it and see how many columns it has */
    ierr = LMEGetSolution(lme,&X);CHKERRQ(ierr);
    ierr = MatLRCGetMats(X,NULL,&X1,NULL,NULL);CHKERRQ(ierr);
    ierr = MatGetSize(X1,NULL,&rank);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," The solver has computed a solution with rank=%d\n",rank);CHKERRQ(ierr);
  }

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = LMEGetIterationNumber(lme,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = LMEGetDimensions(lme,&ncv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Subspace dimension: %D\n",ncv);CHKERRQ(ierr);
  ierr = LMEGetTolerances(lme,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Compute residual error
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = LMEGetErrorEstimate(lme,&errest);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Error estimate reported by the solver: %.4g\n",(double)errest);CHKERRQ(ierr);
  if (n<=150) {
    ierr = LMEComputeError(lme,&error);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Computed residual norm: %.4g\n\n",(double)error);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD," Matrix too large to compute residual norm\n\n");CHKERRQ(ierr);
  }

  /*
     Free work space
  */
  ierr = LMEDestroy(&lme);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      requires: !single

   test:
      suffix: 2
      args: -rank 40
      requires: !single

TEST*/
