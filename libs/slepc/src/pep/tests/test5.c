/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test PEP view and monitor functionality.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            A[3];
  PEP            pep;
  PetscInt       n=6,Istart,Iend,i;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPEP of diagonal problem, n=%D\n\n",n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatCreate(PETSC_COMM_WORLD,&A[0]);CHKERRQ(ierr);
  ierr = MatSetSizes(A[0],PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A[0]);CHKERRQ(ierr);
  ierr = MatSetUp(A[0]);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A[0],&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(A[0],i,i,i+1,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A[0],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A[0],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A[1]);CHKERRQ(ierr);
  ierr = MatSetSizes(A[1],PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A[1]);CHKERRQ(ierr);
  ierr = MatSetUp(A[1]);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(A[1],i,i,-1.5,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A[2]);CHKERRQ(ierr);
  ierr = MatSetSizes(A[2],PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A[2]);CHKERRQ(ierr);
  ierr = MatSetUp(A[2]);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(A[2],i,i,-1.0/(i+1),INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A[2],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A[2],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                     Create the PEP solver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)pep,"pep");CHKERRQ(ierr);
  ierr = PEPSetOperators(pep,3,A);CHKERRQ(ierr);
  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the eigensystem and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PEPSolve(pep);CHKERRQ(ierr);
  ierr = PEPErrorView(pep,PEP_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  ierr = MatDestroy(&A[0]);CHKERRQ(ierr);
  ierr = MatDestroy(&A[1]);CHKERRQ(ierr);
  ierr = MatDestroy(&A[2]);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -pep_error_backward ::ascii_info_detail -pep_largest_real -pep_view_values -pep_monitor_conv -pep_error_absolute ::ascii_matlab -pep_monitor_all -pep_converged_reason -pep_view
      requires: !single
      filter: grep -v "tolerance" | grep -v "problem type" | sed -e "s/[+-]0.000000i//" -e "s/\([0-9]\.[5]*\)[+-][0-9]\.[0-9]*e-[0-9]*i/\\1/g" -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g"

   test:
      suffix: 2
      args: -n 12 -pep_largest_real -pep_monitor -pep_view_values ::ascii_matlab
      requires: double
      filter: sed -e "s/[+-][0-9]\.[0-9]*e-[0-9]*i//" -e "s/[0-9]\.[0-9]*e-\([0-9]*\)/removed/g" -e "s/5\.\([49]\)999999[0-9]*e+00/5.\\1999999999999999e+00/"

TEST*/
