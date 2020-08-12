/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test NEP view and monitor functionality.\n\n";

#include <slepcnep.h>

int main(int argc,char **argv)
{
  Mat                A[3];
  FN                 f[3];
  NEP                nep;
  Vec                xr,xi;
  PetscScalar        kr,ki,coeffs[3];
  PetscInt           n=6,i,Istart,Iend,nconv,its;
  PetscErrorCode     ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nDiagonal Nonlinear Eigenproblem, n=%D\n\n",n);CHKERRQ(ierr);

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

  /*
     Functions: f0=1.0, f1=lambda, f2=lambda^2
  */
  ierr = FNCreate(PETSC_COMM_WORLD,&f[0]);CHKERRQ(ierr);
  ierr = FNSetType(f[0],FNRATIONAL);CHKERRQ(ierr);
  coeffs[0] = 1.0;
  ierr = FNRationalSetNumerator(f[0],1,coeffs);CHKERRQ(ierr);

  ierr = FNCreate(PETSC_COMM_WORLD,&f[1]);CHKERRQ(ierr);
  ierr = FNSetType(f[1],FNRATIONAL);CHKERRQ(ierr);
  coeffs[0] = 1.0; coeffs[1] = 0.0;
  ierr = FNRationalSetNumerator(f[1],2,coeffs);CHKERRQ(ierr);

  ierr = FNCreate(PETSC_COMM_WORLD,&f[2]);CHKERRQ(ierr);
  ierr = FNSetType(f[2],FNRATIONAL);CHKERRQ(ierr);
  coeffs[0] = 1.0; coeffs[1] = 0.0; coeffs[2] = 0.0;
  ierr = FNRationalSetNumerator(f[2],3,coeffs);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Create the NEP solver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = NEPCreate(PETSC_COMM_WORLD,&nep);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)nep,"nep");CHKERRQ(ierr);
  ierr = NEPSetSplitOperator(nep,3,A,f,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = NEPSetTarget(nep,1.1);CHKERRQ(ierr);
  ierr = NEPSetWhichEigenpairs(nep,NEP_TARGET_MAGNITUDE);CHKERRQ(ierr);
  ierr = NEPSetFromOptions(nep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the eigensystem and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = NEPSolve(nep);CHKERRQ(ierr);
  ierr = NEPGetConverged(nep,&nconv);CHKERRQ(ierr);
  ierr = NEPGetIterationNumber(nep,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," %D converged eigenpairs after %D iterations\n",nconv,its);CHKERRQ(ierr);
  if (nconv>0) {
    ierr = MatCreateVecs(A[0],&xr,&xi);CHKERRQ(ierr);
    ierr = NEPGetEigenpair(nep,0,&kr,&ki,xr,xi);CHKERRQ(ierr);
    ierr = VecDestroy(&xr);CHKERRQ(ierr);
    ierr = VecDestroy(&xi);CHKERRQ(ierr);
  }

  ierr = NEPErrorView(nep,NEP_ERROR_BACKWARD,NULL);CHKERRQ(ierr);
  ierr = NEPDestroy(&nep);CHKERRQ(ierr);
  ierr = MatDestroy(&A[0]);CHKERRQ(ierr);
  ierr = MatDestroy(&A[1]);CHKERRQ(ierr);
  ierr = MatDestroy(&A[2]);CHKERRQ(ierr);
  ierr = FNDestroy(&f[0]);CHKERRQ(ierr);
  ierr = FNDestroy(&f[1]);CHKERRQ(ierr);
  ierr = FNDestroy(&f[2]);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -nep_type slp -nep_target -.5 -nep_error_backward ::ascii_info_detail -nep_view_values -nep_error_absolute ::ascii_matlab -nep_monitor_all -nep_converged_reason -nep_view
      filter: grep -v "tolerance" | grep -v "problem type" | sed -e "s/[+-]0.000000i//" -e "s/+0i//" -e "s/[+-][0-9]\.[0-9]*e-[0-9]*i//g" -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double

   test:
      suffix: 2
      args: -nep_type rii -nep_target -.5 -nep_rii_hermitian -nep_monitor -nep_view_values ::ascii_matlab
      filter: sed -e "s/[+-][0-9]\.[0-9]*e-[0-9]*i//" -e "s/([0-9]\.[0-9]*e[+-]\([0-9]*\))/(removed)/g"
      requires: double

TEST*/
