/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test LME interface functions, based on ex32.c.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepclme.h>

int main(int argc,char **argv)
{
  Mat                  A,B,C,C1,D;
  LME                  lme;
  PetscReal            tol,errest,error;
  PetscScalar          *u;
  PetscInt             N,n=10,m,Istart,Iend,II,maxit,ncv,i,j;
  PetscErrorCode       ierr;
  PetscBool            flg,testprefix=PETSC_FALSE,viewmatrices=PETSC_FALSE;
  const char           *prefix;
  LMEType              type;
  LMEProblemType       ptype;
  PetscViewerAndFormat *vf;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&flg);CHKERRQ(ierr);
  if (!flg) m=n;
  N = n*m;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nLyapunov equation, N=%D (%Dx%D grid)\n\n",N,n,m);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-test_prefix",&testprefix,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-view_matrices",&viewmatrices,NULL);CHKERRQ(ierr);

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
    ierr = MatSetValue(A,II,II,-4.0,INSERT_VALUES);CHKERRQ(ierr);
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
  ierr = LMECreate(PETSC_COMM_WORLD,&lme);CHKERRQ(ierr);
  ierr = LMESetProblemType(lme,LME_SYLVESTER);CHKERRQ(ierr);
  ierr = LMEGetProblemType(lme,&ptype);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Equation type set to %D\n",ptype);CHKERRQ(ierr);
  ierr = LMESetProblemType(lme,LME_LYAPUNOV);CHKERRQ(ierr);
  ierr = LMEGetProblemType(lme,&ptype);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Equation type changed to %D\n",ptype);CHKERRQ(ierr);
  ierr = LMESetCoefficients(lme,A,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = LMESetRHS(lme,C);CHKERRQ(ierr);

  /* test prefix usage */
  if (testprefix) {
    ierr = LMESetOptionsPrefix(lme,"check_");CHKERRQ(ierr);
    ierr = LMEAppendOptionsPrefix(lme,"myprefix_");CHKERRQ(ierr);
    ierr = LMEGetOptionsPrefix(lme,&prefix);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," LME prefix is currently: %s\n",prefix);CHKERRQ(ierr);
  }

  /* test some interface functions */
  ierr = LMEGetCoefficients(lme,&B,NULL,NULL,NULL);CHKERRQ(ierr);
  if (viewmatrices) { ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); }
  ierr = LMEGetRHS(lme,&D);CHKERRQ(ierr);
  if (viewmatrices) { ierr = MatView(D,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); }
  ierr = LMESetTolerances(lme,PETSC_DEFAULT,100);CHKERRQ(ierr);
  ierr = LMESetDimensions(lme,21);CHKERRQ(ierr);
  ierr = LMESetErrorIfNotConverged(lme,PETSC_TRUE);CHKERRQ(ierr);
  /* test monitors */
  ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);CHKERRQ(ierr);
  ierr = LMEMonitorSet(lme,(PetscErrorCode (*)(LME,PetscInt,PetscReal,void*))LMEMonitorDefault,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);
  /* ierr = LMEMonitorCancel(lme);CHKERRQ(ierr); */
  ierr = LMESetFromOptions(lme);CHKERRQ(ierr);

  ierr = LMEGetType(lme,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solver being used: %s\n",type);CHKERRQ(ierr);

  /* query properties and print them */
  ierr = LMEGetTolerances(lme,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Tolerance: %g, max iterations: %D\n",(double)tol,maxit);CHKERRQ(ierr);
  ierr = LMEGetDimensions(lme,&ncv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Subspace dimension: %D\n",ncv);CHKERRQ(ierr);
  ierr = LMEGetErrorIfNotConverged(lme,&flg);CHKERRQ(ierr);
  if (flg) { ierr = PetscPrintf(PETSC_COMM_WORLD," Erroring out if convergence fails\n");CHKERRQ(ierr); }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the matrix equation and compute residual error
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = LMESolve(lme);CHKERRQ(ierr);
  ierr = LMEGetErrorEstimate(lme,&errest);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Error estimate reported by the solver: %.4g\n",(double)errest);CHKERRQ(ierr);
  ierr = LMEComputeError(lme,&error);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Computed residual norm: %.4g\n\n",(double)error);CHKERRQ(ierr);

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
      args: -lme_monitor_cancel -lme_converged_reason -lme_view -view_matrices -log_exclude lme,bv
      requires: double
      filter: sed -e "s/4.0[0-9]*e-10/4.03e-10/"

   test:
      suffix: 2
      args: -test_prefix -check_myprefix_lme_monitor
      requires: double
      filter: sed -e "s/estimate [0-9]\.[0-9]*e[+-]\([0-9]*\)/estimate (removed)/g" | sed -e "s/4.0[0-9]*e-10/4.03e-10/"

   test:
      suffix: 3
      args: -lme_monitor_cancel -info
      requires: double
      filter: sed -e "s/equation = [0-9]\.[0-9]*e[+-]\([0-9]*\)/equation = (removed)/g" | sed -e "s/4.0[0-9]*e-10/4.03e-10/" | grep -v Comm | grep -v machine | grep -v PetscGetHostName | grep -v OpenMP | grep -v "Rank of the Cholesky factor" | grep -v "potrf failed"

TEST*/
