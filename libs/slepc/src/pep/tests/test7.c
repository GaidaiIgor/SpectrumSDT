/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test interface functions of spectrum-slicing STOAR.\n\n"
  "This is based on ex38.c. The command line options are:\n"
  "  -n <n> ... dimension of the matrices.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            M,C,K,A[3]; /* problem matrices */
  PEP            pep;        /* polynomial eigenproblem solver context */
  ST             st;         /* spectral transformation context */
  KSP            ksp;
  PC             pc;
  PetscBool      showinertia=PETSC_TRUE,lock,detect,checket;
  PetscInt       n=100,Istart,Iend,i,*inertias,ns,nev,ncv,mpd;
  PetscReal      mu=1,tau=10,kappa=5,int0,int1,*shifts;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-showinertia",&showinertia,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSpectrum slicing on PEP, n=%D\n\n",n);CHKERRQ(ierr);

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
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  A[0] = K; A[1] = C; A[2] = M;
  ierr = PEPSetOperators(pep,3,A);CHKERRQ(ierr);
  ierr = PEPSetProblemType(pep,PEP_HYPERBOLIC);CHKERRQ(ierr);
  ierr = PEPSetType(pep,PEPSTOAR);CHKERRQ(ierr);

  /*
     Set interval and other settings for spectrum slicing
  */
  int0 = -11.3;
  int1 = -9.5;
  ierr = PEPSetInterval(pep,int0,int1);CHKERRQ(ierr);
  ierr = PEPSetWhichEigenpairs(pep,PEP_ALL);CHKERRQ(ierr);
  ierr = PEPGetST(pep,&st);CHKERRQ(ierr);
  ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr);

  /*
     Test interface functions of STOAR solver
  */
  ierr = PEPSTOARGetDetectZeros(pep,&detect);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Detect zeros before changing = %d",(int)detect);CHKERRQ(ierr);
  ierr = PEPSTOARSetDetectZeros(pep,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PEPSTOARGetDetectZeros(pep,&detect);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)detect);CHKERRQ(ierr);

  ierr = PEPSTOARGetLocking(pep,&lock);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Locking flag before changing = %d",(int)lock);CHKERRQ(ierr);
  ierr = PEPSTOARSetLocking(pep,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PEPSTOARGetLocking(pep,&lock);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)lock);CHKERRQ(ierr);

  ierr = PEPSTOARGetCheckEigenvalueType(pep,&checket);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Check eigenvalue type flag before changing = %d",(int)checket);CHKERRQ(ierr);
  ierr = PEPSTOARSetCheckEigenvalueType(pep,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PEPSTOARGetCheckEigenvalueType(pep,&checket);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)checket);CHKERRQ(ierr);

  ierr = PEPSTOARGetDimensions(pep,&nev,&ncv,&mpd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Sub-solve dimensions before changing = [%D,%D,%D]",nev,ncv,mpd);CHKERRQ(ierr);
  ierr = PEPSTOARSetDimensions(pep,30,60,60);CHKERRQ(ierr);
  ierr = PEPSTOARGetDimensions(pep,&nev,&ncv,&mpd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to [%D,%D,%D]\n",nev,ncv,mpd);CHKERRQ(ierr);

  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Compute all eigenvalues in interval and display info
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPSetUp(pep);CHKERRQ(ierr);
  ierr = PEPSTOARGetInertias(pep,&ns,&shifts,&inertias);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Inertias (after setup):\n");CHKERRQ(ierr);
  for (i=0;i<ns;i++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD," .. %g (%D)\n",(double)shifts[i],inertias[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(shifts);CHKERRQ(ierr);
  ierr = PetscFree(inertias);CHKERRQ(ierr);

  ierr = PEPSolve(pep);CHKERRQ(ierr);
  ierr = PEPGetDimensions(pep,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PEPGetInterval(pep,&int0,&int1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Found %D eigenvalues in interval [%g,%g]\n",nev,(double)int0,(double)int1);CHKERRQ(ierr);

  if (showinertia) {
    ierr = PEPSTOARGetInertias(pep,&ns,&shifts,&inertias);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Used %D shifts (inertia):\n",ns);CHKERRQ(ierr);
    for (i=0;i<ns;i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD," .. %g (%D)\n",(double)shifts[i],inertias[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree(shifts);CHKERRQ(ierr);
    ierr = PetscFree(inertias);CHKERRQ(ierr);
  }

  ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      requires: !single
      args: -showinertia 0

TEST*/
