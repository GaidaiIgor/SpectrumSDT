/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This example implements one of the problems found at
       NLEVP: A Collection of Nonlinear Eigenvalue Problems,
       The University of Manchester.
   The details of the collection can be found at:
       [1] T. Betcke et al., "NLEVP: A Collection of Nonlinear Eigenvalue
           Problems", ACM Trans. Math. Software 39(2), Article 7, 2013.

   The pdde_stability problem is a complex-symmetric QEP from the stability
   analysis of a discretized partial delay-differential equation. It requires
   complex scalars.
*/

static char help[] = "Stability analysis of a discretized partial delay-differential equation.\n\n"
  "The command line options are:\n"
  "  -m <m>, grid size, the matrices have dimension n=m*m.\n"
  "  -c <a0,b0,a1,b1,a2,b2,phi1>, comma-separated list of 7 real parameters.\n\n";

#include <slepcpep.h>

#define NMAT 3

/*
    Function for user-defined eigenvalue ordering criterion.

    Given two eigenvalues ar+i*ai and br+i*bi, the subroutine must choose
    one of them as the preferred one according to the criterion.
    In this example, the preferred value is the one with absolute value closest to 1.
*/
PetscErrorCode MyEigenSort(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *r,void *ctx)
{
  PetscReal aa,ab;

  PetscFunctionBeginUser;
  aa = PetscAbsReal(SlepcAbsEigenvalue(ar,ai)-1.0);
  ab = PetscAbsReal(SlepcAbsEigenvalue(br,bi)-1.0);
  *r = aa > ab ? 1 : (aa < ab ? -1 : 0);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  Mat            A[NMAT];         /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscInt       m=15,n,II,Istart,Iend,i,j,k;
  PetscReal      h,xi,xj,c[7] = { 2, .3, -2, .2, -2, -.3, -PETSC_PI/2 };
  PetscScalar    alpha,beta,gamma;
  PetscBool      flg,terse;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
#if !defined(PETSC_USE_COMPLEX)
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP, "This example requires complex scalars");
#endif

  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  n = m*m;
  h = PETSC_PI/(m+1);
  gamma = PetscExpScalar(PETSC_i*c[6]);
  gamma = gamma/PetscAbsScalar(gamma);
  k = 7;
  ierr = PetscOptionsGetRealArray(NULL,NULL,"-c",c,&k,&flg);CHKERRQ(ierr);
  if (flg && k!=7) SETERRQ1(PETSC_COMM_WORLD,1,"The number of parameters -c should be 7, you provided %D",k);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPDDE stability, n=%D (m=%D)\n\n",n,m);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                     Compute the polynomial matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* initialize matrices */
  for (i=0;i<NMAT;i++) {
    ierr = MatCreate(PETSC_COMM_WORLD,&A[i]);CHKERRQ(ierr);
    ierr = MatSetSizes(A[i],PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A[i]);CHKERRQ(ierr);
    ierr = MatSetUp(A[i]);CHKERRQ(ierr);
  }
  ierr = MatGetOwnershipRange(A[0],&Istart,&Iend);CHKERRQ(ierr);

  /* A[1] has a pattern similar to the 2D Laplacian */
  for (II=Istart;II<Iend;II++) {
    i = II/m; j = II-i*m;
    xi = (i+1)*h; xj = (j+1)*h;
    alpha = c[0]+c[1]*PetscSinReal(xi)+gamma*(c[2]+c[3]*xi*(1.0-PetscExpReal(xi-PETSC_PI)));
    beta = c[0]+c[1]*PetscSinReal(xj)-gamma*(c[2]+c[3]*xj*(1.0-PetscExpReal(xj-PETSC_PI)));
    ierr = MatSetValue(A[1],II,II,alpha+beta-4.0/(h*h),INSERT_VALUES);CHKERRQ(ierr);
    if (j>0) { ierr = MatSetValue(A[1],II,II-1,1.0/(h*h),INSERT_VALUES);CHKERRQ(ierr); }
    if (j<m-1) { ierr = MatSetValue(A[1],II,II+1,1.0/(h*h),INSERT_VALUES);CHKERRQ(ierr); }
    if (i>0) { ierr = MatSetValue(A[1],II,II-m,1.0/(h*h),INSERT_VALUES);CHKERRQ(ierr); }
    if (i<m-1) { ierr = MatSetValue(A[1],II,II+m,1.0/(h*h),INSERT_VALUES);CHKERRQ(ierr); }
  }

  /* A[0] and A[2] are diagonal */
  for (II=Istart;II<Iend;II++) {
    i = II/m; j = II-i*m;
    xi = (i+1)*h; xj = (j+1)*h;
    alpha = c[4]+c[5]*xi*(PETSC_PI-xi);
    beta = c[4]+c[5]*xj*(PETSC_PI-xj);
    ierr = MatSetValue(A[0],II,II,alpha,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(A[2],II,II,beta,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* assemble matrices */
  for (i=0;i<NMAT;i++) {
    ierr = MatAssemblyBegin(A[i],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  for (i=0;i<NMAT;i++) {
    ierr = MatAssemblyEnd(A[i],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  ierr = PEPSetOperators(pep,NMAT,A);CHKERRQ(ierr);
  ierr = PEPSetEigenvalueComparison(pep,MyEigenSort,NULL);CHKERRQ(ierr);
  ierr = PEPSetDimensions(pep,4,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);
  ierr = PEPSolve(pep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = PEPReasonView(pep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  for (i=0;i<NMAT;i++) {
    ierr = MatDestroy(&A[i]);CHKERRQ(ierr);
  }
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   build:
      requires: complex

   test:
      suffix: 1
      args: -pep_type {{toar qarnoldi linear}} -pep_ncv 25 -terse
      requires: complex

TEST*/
