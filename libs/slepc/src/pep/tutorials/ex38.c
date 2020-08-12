/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Spectrum slicing on quadratic symmetric eigenproblem.\n\n"
  "The command line options are:\n"
  "  -n <n> ... dimension of the matrices.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat            M,C,K,A[3]; /* problem matrices */
  PEP            pep;        /* polynomial eigenproblem solver context */
  ST             st;         /* spectral transformation context */
  KSP            ksp;
  PC             pc;
  PEPType        type;
  PetscBool      show=PETSC_FALSE,terse;
  PetscInt       n=100,Istart,Iend,nev,i,*inertias,ns;
  PetscReal      mu=1,tau=10,kappa=5,inta,intb,*shifts;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-show_inertias",&show,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
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

  /*
     Create eigensolver context
  */
  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);

  /*
     Set operators and set problem type
  */
  A[0] = K; A[1] = C; A[2] = M;
  ierr = PEPSetOperators(pep,3,A);CHKERRQ(ierr);
  ierr = PEPSetProblemType(pep,PEP_HYPERBOLIC);CHKERRQ(ierr);

  /*
     Set interval for spectrum slicing
  */
  inta = -11.3;
  intb = -9.5;
  ierr = PEPSetInterval(pep,inta,intb);CHKERRQ(ierr);
  ierr = PEPSetWhichEigenpairs(pep,PEP_ALL);CHKERRQ(ierr);

  /*
     Spectrum slicing requires STOAR
  */
  ierr = PEPSetType(pep,PEPSTOAR);CHKERRQ(ierr);

  /*
     Set shift-and-invert with Cholesky; select MUMPS if available
  */
  ierr = PEPGetST(pep,&st);CHKERRQ(ierr);
  ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);

  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr);

#if defined(PETSC_HAVE_MUMPS)
#if defined(PETSC_USE_COMPLEX)
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Spectrum slicing with MUMPS is not available for complex scalars");
#endif
  ierr = PEPSTOARSetDetectZeros(pep,PETSC_TRUE);CHKERRQ(ierr);  /* enforce zero detection */
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
  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PEPSetUp(pep);CHKERRQ(ierr);
  if (show) {
    ierr = PEPSTOARGetInertias(pep,&ns,&shifts,&inertias);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Subintervals (after setup):\n");CHKERRQ(ierr);
    for (i=0;i<ns;i++) { ierr = PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %D \n",(double)shifts[i],inertias[i]);CHKERRQ(ierr); }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    ierr = PetscFree(shifts);CHKERRQ(ierr);
    ierr = PetscFree(inertias);CHKERRQ(ierr);
  }
  ierr = PEPSolve(pep);CHKERRQ(ierr);
  if (show && !terse) {
    ierr = PEPSTOARGetInertias(pep,&ns,&shifts,&inertias);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"All shifts (after solve):\n");CHKERRQ(ierr);
    for (i=0;i<ns;i++) { ierr = PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %D \n",(double)shifts[i],inertias[i]);CHKERRQ(ierr); }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    ierr = PetscFree(shifts);CHKERRQ(ierr);
    ierr = PetscFree(inertias);CHKERRQ(ierr);
  }

  /*
     Show eigenvalues in interval and print solution
  */
  ierr = PEPGetType(pep,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = PEPGetDimensions(pep,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PEPGetInterval(pep,&inta,&intb);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," %D eigenvalues found in [%g, %g]\n",nev,(double)inta,(double)intb);CHKERRQ(ierr);

  /*
     Show detailed info unless -terse option is given by user
   */
  if (terse) {
    ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = PEPReasonView(pep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

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

   testset:
      requires: !single
      args: -show_inertias -terse
      output_file: output/ex38_1.out
      test:
         suffix: 1
         requires: !complex
      test:
         suffix: 1_complex
         requires: complex !mumps

TEST*/
