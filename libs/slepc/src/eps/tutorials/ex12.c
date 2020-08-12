/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Compute all eigenvalues in an interval of a symmetric-definite problem.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A,B;         /* matrices */
  EPS            eps;         /* eigenproblem solver context */
  ST             st;          /* spectral transformation context */
  KSP            ksp;
  PC             pc;
  PetscInt       N,n=35,m,Istart,Iend,II,nev,i,j,k,*inertias;
  PetscBool      flag,showinertia=PETSC_TRUE;
  PetscReal      int0,int1,*shifts;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetBool(NULL,NULL,"-showinertia",&showinertia,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag);CHKERRQ(ierr);
  if (!flag) m=n;
  N = n*m;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSymmetric-definite problem with two intervals, N=%D (%Dx%D grid)\n\n",N,n,m);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(B);CHKERRQ(ierr);
  ierr = MatSetUp(B);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) { ierr = MatSetValue(A,II,II-n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<m-1) { ierr = MatSetValue(A,II,II+n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j>0) { ierr = MatSetValue(A,II,II-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(A,II,II+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,II,II,4.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(B,II,II,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  if (Istart==0) {
    ierr = MatSetValue(B,0,0,6.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(B,0,1,-1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(B,1,0,-1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(B,1,1,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);

  /*
     Set first interval and other settings for spectrum slicing
  */
  ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRQ(ierr);
  ierr = EPSSetInterval(eps,1.1,1.3);CHKERRQ(ierr);
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Solve for first interval and display info
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = EPSGetInterval(eps,&int0,&int1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Found %D eigenvalues in interval [%g,%g]\n",nev,(double)int0,(double)int1);CHKERRQ(ierr);
  if (showinertia) {
    ierr = EPSKrylovSchurGetInertias(eps,&k,&shifts,&inertias);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Used %D shifts (inertia):\n",k);CHKERRQ(ierr);
    for (i=0;i<k;i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD," .. %g (%D)\n",(double)shifts[i],inertias[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree(shifts);CHKERRQ(ierr);
    ierr = PetscFree(inertias);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Solve for second interval and display info
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSetInterval(eps,1.499,1.6);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = EPSGetInterval(eps,&int0,&int1);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Found %D eigenvalues in interval [%g,%g]\n",nev,(double)int0,(double)int1);CHKERRQ(ierr);
  if (showinertia) {
    ierr = EPSKrylovSchurGetInertias(eps,&k,&shifts,&inertias);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Used %D shifts (inertia):\n",k);CHKERRQ(ierr);
    for (i=0;i<k;i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD," .. %g (%D)\n",(double)shifts[i],inertias[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree(shifts);CHKERRQ(ierr);
    ierr = PetscFree(inertias);CHKERRQ(ierr);
  }

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -showinertia 0 -eps_error_relative
      requires: !single

TEST*/
