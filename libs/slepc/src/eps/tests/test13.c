/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test EPSSetArbitrarySelection.\n\n";

#include <slepceps.h>

PetscErrorCode MyArbitrarySelection(PetscScalar eigr,PetscScalar eigi,Vec xr,Vec xi,PetscScalar *rr,PetscScalar *ri,void *ctx)
{
  PetscErrorCode  ierr;
  Vec             xref = *(Vec*)ctx;

  PetscFunctionBeginUser;
  ierr = VecDot(xr,xref,rr);CHKERRQ(ierr);
  *rr = PetscAbsScalar(*rr);
  if (ri) *ri = 0.0;
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrices */
  EPS            eps;         /* eigenproblem solver context */
  PetscScalar    seigr,seigi;
  PetscReal      tol=1000*PETSC_MACHINE_EPSILON;
  Vec            sxr,sxi;
  PetscInt       n=30,i,Istart,Iend,nconv;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nTridiagonal with zero diagonal, n=%D\n\n",n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           Create matrix tridiag([-1 0 -1])
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(A,i,i-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A,i,i+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Create the eigensolver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,tol,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve eigenproblem and store some solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,&sxr,NULL);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,&sxi,NULL);CHKERRQ(ierr);
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  if (nconv>0) {
    ierr = EPSGetEigenpair(eps,0,&seigr,&seigi,sxr,sxi);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Solve eigenproblem using an arbitrary selection
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = EPSSetArbitrarySelection(eps,MyArbitrarySelection,&sxr);CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE);CHKERRQ(ierr);
    ierr = EPSSolve(eps);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem: no eigenpairs converged.\n");CHKERRQ(ierr);
  }

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = VecDestroy(&sxr);CHKERRQ(ierr);
  ierr = VecDestroy(&sxi);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -eps_max_it 5000 -st_pc_type jacobi
      output_file: output/test13_1.out
      test:
         suffix: 1
         args: -eps_type {{krylovschur gd jd}}
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion
      test:
         suffix: 2
         args: -eps_non_hermitian -eps_type {{krylovschur gd jd}}
      test:
         suffix: 2_gd2
         args: -eps_non_hermitian -eps_type gd -eps_gd_double_expansion
         timeoutfactor: 2

TEST*/
