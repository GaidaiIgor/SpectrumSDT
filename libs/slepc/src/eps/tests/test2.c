/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests multiple calls to EPSSolve with the same matrix.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  ST             st;
  PetscReal      tol=PetscMax(1000*PETSC_MACHINE_EPSILON,1e-9);
  PetscInt       n=30,i,Istart,Iend;
  PetscBool      flg;
  PetscErrorCode ierr;
  EPSLanczosReorthogType reorth;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem, n=%D\n\n",n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(A,i,i-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A,i,i+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,i,i,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Create the eigensolver
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,tol,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* illustrate how to extract parameters from specific solver types */
  ierr = PetscObjectTypeCompare((PetscObject)eps,EPSLANCZOS,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EPSLanczosGetReorthog(eps,&reorth);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reorthogonalization type used in Lanczos: %s\n",EPSLanczosReorthogTypes[reorth]);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Solve for largest eigenvalues
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," - - - Largest eigenvalues - - -\n");CHKERRQ(ierr);
  ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Solve for smallest eigenvalues
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," - - - Smallest eigenvalues - - -\n");CHKERRQ(ierr);
  ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Solve for interior eigenvalues (target=2.1)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);CHKERRQ(ierr);
  ierr = EPSSetTarget(eps,2.1);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)eps,EPSLANCZOS,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
    ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
  } else {
    ierr = PetscObjectTypeCompare((PetscObject)eps,EPSKRYLOVSCHUR,&flg);CHKERRQ(ierr);
    if (!flg) {
      ierr = PetscObjectTypeCompare((PetscObject)eps,EPSARNOLDI,&flg);CHKERRQ(ierr);
    }
    if (flg) {
      ierr = EPSSetExtraction(eps,EPS_HARMONIC);CHKERRQ(ierr);
    }
  }
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," - - - Interior eigenvalues - - -\n");CHKERRQ(ierr);
  ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -eps_nev 4
      requires: !single
      output_file: output/test2_1.out
      test:
         suffix: 1
         args: -eps_type {{arnoldi gd jd lapack}}
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion
         timeoutfactor: 2
      test:
         suffix: 1_krylovschur
         args: -eps_type krylovschur -eps_krylovschur_locking {{0 1}}

   testset:
      args: -eps_type lanczos -eps_nev 4
      requires: !single
      filter: grep -v "Lanczos"
      output_file: output/test2_1.out
      test:
         suffix: 2
         args: -eps_lanczos_reorthog {{local full selective periodic partial}}

   testset:
      args: -n 32 -eps_nev 4
      requires: !single
      output_file: output/test2_3.out
      test:
         nsize: 2
         suffix: 3
         args: -eps_type {{krylovschur lapack}}
      test:
         nsize: 2
         suffix: 3_gd
         args: -eps_type gd -eps_gd_krylov_start
         timeoutfactor: 2
      test:
         suffix: 3_jd
         args: -eps_type jd -eps_jd_krylov_start -eps_ncv 18

   testset:
      args: -eps_nev 4 -mat_type aijcusparse
      requires: cuda !single
      output_file: output/test2_1.out
      test:
         suffix: 4_cuda
         args: -eps_type {{krylovschur arnoldi gd jd}}

TEST*/
