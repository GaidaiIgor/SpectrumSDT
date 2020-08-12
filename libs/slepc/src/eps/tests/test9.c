/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Eigenvalue problem associated with a Markov model of a random walk on a triangular grid. "
  "It is a standard nonsymmetric eigenproblem with real eigenvalues and the rightmost eigenvalue is known to be 1.\n"
  "This example illustrates how the user can set the initial vector.\n\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = number of grid subdivisions in each dimension.\n\n";

#include <slepceps.h>

/*
   User-defined routines
*/
PetscErrorCode MatMarkovModel(PetscInt m,Mat A);
PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx);

int main(int argc,char **argv)
{
  Vec            v0;              /* initial vector */
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  PetscReal      tol=1000*PETSC_MACHINE_EPSILON;
  PetscInt       N,m=15,nev;
  PetscScalar    origin=0.0;
  PetscBool      flg,delay;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  N = m*(m+1)/2;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nMarkov Model, N=%D (m=%D)\n\n",N,m);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = MatMarkovModel(m,A);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,tol,PETSC_DEFAULT);CHKERRQ(ierr);

  /*
     Set the custom comparing routine in order to obtain the eigenvalues
     closest to the target on the right only
  */
  ierr = EPSSetEigenvalueComparison(eps,MyEigenSort,&origin);CHKERRQ(ierr);


  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)eps,EPSARNOLDI,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EPSArnoldiGetDelayed(eps,&delay);CHKERRQ(ierr);
    if (delay) {
      ierr = PetscPrintf(PETSC_COMM_WORLD," Warning: delayed reorthogonalization may be unstable\n");CHKERRQ(ierr);
    }
  }

  /*
     Set the initial vector. This is optional, if not done the initial
     vector is set to random values
  */
  ierr = MatCreateVecs(A,&v0,NULL);CHKERRQ(ierr);
  ierr = VecSetValue(v0,0,-1.5,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecSetValue(v0,1,2.1,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(v0);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(v0);CHKERRQ(ierr);
  ierr = EPSSetInitialSpace(eps,1,&v0);CHKERRQ(ierr);

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
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&v0);CHKERRQ(ierr);
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

/*
    Function for user-defined eigenvalue ordering criterion.

    Given two eigenvalues ar+i*ai and br+i*bi, the subroutine must choose
    one of them as the preferred one according to the criterion.
    In this example, the preferred value is the one furthest to the origin.
*/
PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  PetscScalar origin = *(PetscScalar*)ctx;
  PetscReal   d;

  PetscFunctionBeginUser;
  d = (SlepcAbsEigenvalue(br-origin,bi) - SlepcAbsEigenvalue(ar-origin,ai))/PetscMax(SlepcAbsEigenvalue(ar-origin,ai),SlepcAbsEigenvalue(br-origin,bi));
  *r = d > PETSC_SQRT_MACHINE_EPSILON ? 1 : (d < -PETSC_SQRT_MACHINE_EPSILON ? -1 : PetscSign(PetscRealPart(br)));
  PetscFunctionReturn(0);
}

/*TEST

   testset:
      args: -eps_nev 4
      requires: !single
      output_file: output/test9_1.out
      test:
         suffix: 1
         args: -eps_type {{krylovschur arnoldi lapack}} -eps_ncv 7 -eps_max_it 300
      test:
         suffix: 1_gd
         args: -eps_type gd -st_pc_type none
      test:
         suffix: 1_gd2
         args: -eps_type gd -eps_gd_double_expansion -st_pc_type none

   test:
      suffix: 2
      args: -eps_balance {{none oneside twoside}} -eps_krylovschur_locking {{0 1}} -eps_nev 4 -eps_ncv 7 -eps_max_it 500
      requires: double
      output_file: output/test9_1.out

   test:
      suffix: 3
      nsize: 2
      args: -eps_type arnoldi -eps_arnoldi_delayed -eps_largest_real -eps_nev 3 -eps_tol 1e-7 -bv_orthog_refine {{never ifneeded}}
      requires: !single
      output_file: output/test9_3.out

   test:
      suffix: 4
      args: -eps_nev 4 -eps_true_residual
      requires: !single
      output_file: output/test9_1.out

   test:
      suffix: 5
      args: -eps_type jd -eps_nev 3 -eps_target .5 -eps_harmonic -st_ksp_type bicg -st_pc_type lu -eps_jd_minv 2
      requires: !single

   test:
      suffix: 5_arpack
      args: -eps_nev 3 -st_type sinvert -eps_target .5 -eps_type arpack -eps_ncv 10
      requires: arpack
      output_file: output/test9_5.out

   testset:
      args: -eps_type ciss -eps_tol 1e-9 -rg_type ellipse -rg_ellipse_center 0.55 -rg_ellipse_radius 0.05 -rg_ellipse_vscale 0.1 -eps_ciss_usest 0 -eps_all
      requires: !single
      output_file: output/test9_6.out
      test:
         suffix: 6
      test:
         suffix: 6_hankel
         args: -eps_ciss_extraction hankel -eps_ciss_spurious_threshold 1e-6
      test:
         suffix: 6_cheby
         args: -eps_ciss_quadrule chebyshev
      test:
         suffix: 6_hankel_cheby
         args: -eps_ciss_extraction hankel -eps_ciss_quadrule chebyshev
      test:
         suffix: 6_refine
         args: -eps_ciss_refine_inner 1 -eps_ciss_refine_blocksize 1
      test:
         suffix: 6_bcgs
         args: -eps_ciss_realmats -eps_ciss_ksp_type bcgs -eps_ciss_pc_type sor -eps_ciss_integration_points 12

   test:
      suffix: 7
      args: -eps_nev 4 -eps_two_sided -eps_view_vectors ::ascii_info -eps_view_values
      requires: !single
      filter: sed -e "s/\(0x[0-9a-fA-F]*\)/objectid/"

   test:
      suffix: 8
      args: -eps_nev 4 -eps_view_values draw -eps_monitor_lg
      requires: x !single
      output_file: output/test9_1.out

   test:
      suffix: 9
      args: -eps_nev 3 -st_type sinvert -eps_target .5 -st_ksp_type hpddm -st_ksp_hpddm_type gcrodr -eps_ncv 10
      requires: hpddm
      output_file: output/test9_5.out

   test:
      suffix: 10
      args: -eps_nev 3 -st_type sinvert -eps_target .5 -st_pc_type hpddm -st_pc_hpddm_coarse_pc_type lu -eps_ncv 10
      requires: hpddm
      output_file: output/test9_5.out

TEST*/
