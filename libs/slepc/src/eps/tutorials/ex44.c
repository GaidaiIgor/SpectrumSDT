/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Compute rightmost eigenvalues with Lyapunov inverse iteration.\n\n"
  "Loads matrix from a file or builds the same problem as ex36.c (with fixed parameters).\n\n"
  "The command line options are:\n"
  "  -file <filename>, where <filename> = matrix file in PETSc binary form.\n"
  "  -shift <sigma>, shift to make the matrix stable.\n"
  "  -n <n>, block dimension of the 2x2 block matrix (if matrix is generated).\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscScalar    alpha,beta,tau1,tau2,delta1,delta2,L,h,sigma=0.0;
  PetscInt       n=30,i,Istart,Iend,nev;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;
  PetscBool      flg,terse;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetString(NULL,NULL,"-file",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (flg) {

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Load the matrix from file
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nEigenproblem stored in file.\n\n");CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrix from a binary file...\n");CHKERRQ(ierr);
#else
    ierr = PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrix from a binary file...\n");CHKERRQ(ierr);
#endif
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatLoad(A,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  } else {

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Generate Brusselator matrix
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nBrusselator wave model, n=%D\n\n",n);CHKERRQ(ierr);

    alpha  = 2.0;
    beta   = 5.45;
    delta1 = 0.008;
    delta2 = 0.004;
    L      = 0.51302;

    h = 1.0 / (PetscReal)(n+1);
    tau1 = delta1 / ((h*L)*(h*L));
    tau2 = delta2 / ((h*L)*(h*L));

    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,2*n,2*n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
    for (i=Istart;i<Iend;i++) {
      if (i<n) {  /* upper blocks */
        if (i>0) { ierr = MatSetValue(A,i,i-1,tau1,INSERT_VALUES);CHKERRQ(ierr); }
        if (i<n-1) { ierr = MatSetValue(A,i,i+1,tau1,INSERT_VALUES);CHKERRQ(ierr); }
        ierr = MatSetValue(A,i,i,-2.0*tau1+beta-1.0,INSERT_VALUES);CHKERRQ(ierr);
        ierr = MatSetValue(A,i,i+n,alpha*alpha,INSERT_VALUES);CHKERRQ(ierr);
      } else {  /* lower blocks */
        if (i>n) { ierr = MatSetValue(A,i,i-1,tau2,INSERT_VALUES);CHKERRQ(ierr); }
        if (i<2*n-1) { ierr = MatSetValue(A,i,i+1,tau2,INSERT_VALUES);CHKERRQ(ierr); }
        ierr = MatSetValue(A,i,i,-2.0*tau2-alpha*alpha,INSERT_VALUES);CHKERRQ(ierr);
        ierr = MatSetValue(A,i,i-n,-beta,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  /* Shift the matrix to make it stable, A-sigma*I */
  ierr = PetscOptionsGetScalar(NULL,NULL,"-shift",&sigma,NULL);CHKERRQ(ierr);
  ierr = MatShift(A,-sigma);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -eps_nev 6 -shift 0.1 -eps_type {{krylovschur lyapii}} -eps_tol 1e-7 -terse
      requires: double
      filter: grep -v method | sed -e "s/-0.09981-2.13938i, -0.09981+2.13938i/-0.09981+2.13938i, -0.09981-2.13938i/" | sed -e "s/-0.77192-2.52712i, -0.77192+2.52712i/-0.77192+2.52712i, -0.77192-2.52712i/" | sed -e "s/-1.88445-3.02666i, -1.88445+3.02666i/-1.88445+3.02666i, -1.88445-3.02666i/"
      output_file: output/ex44_1.out
      test:
         suffix: 1
      test:
         suffix: 2
         args: -eps_lyapii_ranks 8,20 -options_left no

TEST*/
