/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Spectrum slicing on generalized symmetric eigenproblem.\n\n"
  "The problem is similar to ex13.c.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A,B;         /* matrices */
  EPS            eps;         /* eigenproblem solver context */
  ST             st;          /* spectral transformation context */
  KSP            ksp;
  PC             pc;
  EPSType        type;
  PetscInt       N,n=10,m,Istart,Iend,II,nev,i,j,*inertias,ns;
  PetscReal      inta,intb,*shifts;
  PetscBool      flag,show=PETSC_FALSE,terse;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-show_inertias",&show,NULL);CHKERRQ(ierr);
  if (!flag) m=n;
  N = n*m;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSpectrum slicing on GHEP, N=%D (%Dx%D grid)\n\n",N,n,m);CHKERRQ(ierr);

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
    ierr = MatSetValue(B,II,II,4.0,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators and set problem type
  */
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);

  /*
     Set interval for spectrum slicing
  */
  inta = 0.1;
  intb = 0.2;
  ierr = EPSSetInterval(eps,inta,intb);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRQ(ierr);

  /*
     Spectrum slicing requires Krylov-Schur
  */
  ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRQ(ierr);

  /*
     Set shift-and-invert with Cholesky; select MUMPS if available
  */

  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);

  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr);

#if defined(PETSC_HAVE_MUMPS)
#if defined(PETSC_USE_COMPLEX)
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Spectrum slicing with MUMPS is not available for complex scalars");
#endif
  ierr = EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE);CHKERRQ(ierr);  /* enforce zero detection */
  ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);CHKERRQ(ierr);
  /*
     Add several MUMPS options (see ex43.c for a better way of setting them in program):
     '-mat_mumps_icntl_13 1': turn off ScaLAPACK for matrix inertia
     '-mat_mumps_icntl_24 1': detect null pivots in factorization (for the case that a shift is equal to an eigenvalue)
     '-mat_mumps_cntl_3 <tol>': a tolerance used for null pivot detection (must be larger than machine epsilon)

     Note: depending on the interval, it may be necessary also to increase the workspace:
     '-mat_mumps_icntl_14 <percentage>': increase workspace with a percentage (50, 100 or more)
  */
  ierr = PetscOptionsInsertString(NULL,"-mat_mumps_icntl_13 1 -mat_mumps_icntl_24 1 -mat_mumps_cntl_3 1e-12");CHKERRQ(ierr);
#endif

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSetUp(eps);CHKERRQ(ierr);
  if (show) {
    ierr = EPSKrylovSchurGetInertias(eps,&ns,&shifts,&inertias);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Subintervals (after setup):\n");CHKERRQ(ierr);
    for (i=0;i<ns;i++) { ierr = PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %D \n",(double)shifts[i],inertias[i]);CHKERRQ(ierr); }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    ierr = PetscFree(shifts);CHKERRQ(ierr);
    ierr = PetscFree(inertias);CHKERRQ(ierr);
  }
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  if (show) {
    ierr = EPSKrylovSchurGetInertias(eps,&ns,&shifts,&inertias);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"All shifts (after solve):\n");CHKERRQ(ierr);
    for (i=0;i<ns;i++) { ierr = PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %D \n",(double)shifts[i],inertias[i]);CHKERRQ(ierr); }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    ierr = PetscFree(shifts);CHKERRQ(ierr);
    ierr = PetscFree(inertias);CHKERRQ(ierr);
  }

  /*
     Show eigenvalues in interval and print solution
  */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = EPSGetInterval(eps,&inta,&intb);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," %D eigenvalues found in [%g, %g]\n",nev,(double)inta,(double)intb);CHKERRQ(ierr);

  /*
     Show detailed info unless -terse option is given by user
   */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -terse
      output_file: output/ex25_1.out
      test:
         suffix: 1
         requires: !single
      test:
         suffix: 1_single
         args: -eps_tol 1e-5
         requires: single

TEST*/
