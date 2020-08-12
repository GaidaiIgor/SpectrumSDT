/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates region filtering. "
  "Based on ex5.\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = number of grid subdivisions in each dimension.\n\n";

#include <slepceps.h>

/*
   User-defined routines
*/
PetscErrorCode MatMarkovModel(PetscInt m,Mat A);

int main(int argc,char **argv)
{
  Mat            A;
  EPS            eps;
  ST             st;
  RG             rg;
  PetscReal      radius,tol=1000*PETSC_MACHINE_EPSILON;
  PetscScalar    target=0.5;
  PetscInt       N,m=15,nev;
  PetscErrorCode ierr;
  char           str[50];

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  N = m*(m+1)/2;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nMarkov Model, N=%D (m=%D)\n",N,m);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-target",&target,NULL);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,target,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Searching closest eigenvalues to the right of %s.\n\n",str);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = MatMarkovModel(m,A);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create a standalone spectral transformation
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = STCreate(PETSC_COMM_WORLD,&st);CHKERRQ(ierr);
  ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
  ierr = STSetShift(st,target);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Create a region for filtering
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = RGCreate(PETSC_COMM_WORLD,&rg);CHKERRQ(ierr);
  ierr = RGSetType(rg,RGELLIPSE);CHKERRQ(ierr);
  radius = (1.0-PetscRealPart(target))/2.0;
  ierr = RGEllipseSetParameters(rg,target+radius,radius,0.1);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetST(eps,st);CHKERRQ(ierr);
  ierr = EPSSetRG(eps,rg);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,tol,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = EPSSetTarget(eps,target);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = STDestroy(&st);CHKERRQ(ierr);
  ierr = RGDestroy(&rg);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

PetscErrorCode MatMarkovModel(PetscInt m,Mat A)
{
  const PetscReal cst = 0.5/(PetscReal)(m-1);
  PetscReal       pd,pu;
  PetscInt        Istart,Iend,i,j,jmax,ix=0;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=1;i<=m;i++) {
    jmax = m-i+1;
    for (j=1;j<=jmax;j++) {
      ix = ix + 1;
      if (ix-1<Istart || ix>Iend) continue;  /* compute only owned rows */
      if (j!=jmax) {
        pd = cst*(PetscReal)(i+j-1);
        /* north */
        if (i==1) {
          ierr = MatSetValue(A,ix-1,ix,2*pd,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          ierr = MatSetValue(A,ix-1,ix,pd,INSERT_VALUES);CHKERRQ(ierr);
        }
        /* east */
        if (j==1) {
          ierr = MatSetValue(A,ix-1,ix+jmax-1,2*pd,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          ierr = MatSetValue(A,ix-1,ix+jmax-1,pd,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
      /* south */
      pu = 0.5 - cst*(PetscReal)(i+j-3);
      if (j>1) {
        ierr = MatSetValue(A,ix-1,ix-2,pu,INSERT_VALUES);CHKERRQ(ierr);
      }
      /* west */
      if (i>1) {
        ierr = MatSetValue(A,ix-1,ix-jmax-2,pu,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      args: -eps_nev 4 -eps_ncv 20
      output_file: output/test11_1.out
      requires: !single

TEST*/
