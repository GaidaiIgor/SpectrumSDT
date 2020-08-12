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

   The loaded_string problem is a rational eigenvalue problem for the
   finite element model of a loaded vibrating string.
*/

static char help[] = "Illustrates computation of left eigenvectors and resolvent.\n\n"
  "This is based on loaded_string from the NLEVP collection.\n"
  "The command line options are:\n"
  "  -n <n>, dimension of the matrices.\n"
  "  -kappa <kappa>, stiffness of elastic spring.\n"
  "  -mass <m>, mass of the attached load.\n\n";

#include <slepcnep.h>

#define NMAT 3

int main(int argc,char **argv)
{
  Mat            A[NMAT];         /* problem matrices */
  FN             f[NMAT];         /* functions to define the nonlinear operator */
  NEP            nep;             /* nonlinear eigensolver context */
  RG             rg;
  Vec            v,r,z,w;
  PetscInt       n=100,Istart,Iend,i,nconv;
  PetscReal      kappa=1.0,m=1.0,nrm,tol;
  PetscScalar    lambda,sigma,numer[2],denom[2],omega1,omega2;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-kappa",&kappa,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-mass",&m,NULL);CHKERRQ(ierr);
  sigma = kappa/m;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Loaded vibrating string, n=%D kappa=%g m=%g\n\n",n,(double)kappa,(double)m);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Build the problem matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* initialize matrices */
  for (i=0;i<NMAT;i++) {
    ierr = MatCreate(PETSC_COMM_WORLD,&A[i]);CHKERRQ(ierr);
    ierr = MatSetSizes(A[i],PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A[i]);CHKERRQ(ierr);
    ierr = MatSetUp(A[i]);CHKERRQ(ierr);
  }
  ierr = MatGetOwnershipRange(A[0],&Istart,&Iend);CHKERRQ(ierr);

  /* A0 */
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(A[0],i,i,(i==n-1)?1.0*n:2.0*n,INSERT_VALUES);CHKERRQ(ierr);
    if (i>0) { ierr = MatSetValue(A[0],i,i-1,-1.0*n,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A[0],i,i+1,-1.0*n,INSERT_VALUES);CHKERRQ(ierr); }
  }

  /* A1 */
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(A[1],i,i,(i==n-1)?2.0/(6.0*n):4.0/(6.0*n),INSERT_VALUES);CHKERRQ(ierr);
    if (i>0) { ierr = MatSetValue(A[1],i,i-1,1.0/(6.0*n),INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A[1],i,i+1,1.0/(6.0*n),INSERT_VALUES);CHKERRQ(ierr); }
  }

  /* A2 */
  if (Istart<=n-1 && n-1<Iend) {
    ierr = MatSetValue(A[2],n-1,n-1,kappa,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* assemble matrices */
  for (i=0;i<NMAT;i++) {
    ierr = MatAssemblyBegin(A[i],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  for (i=0;i<NMAT;i++) {
    ierr = MatAssemblyEnd(A[i],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Create the problem functions
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* f1=1 */
  ierr = FNCreate(PETSC_COMM_WORLD,&f[0]);CHKERRQ(ierr);
  ierr = FNSetType(f[0],FNRATIONAL);CHKERRQ(ierr);
  numer[0] = 1.0;
  ierr = FNRationalSetNumerator(f[0],1,numer);CHKERRQ(ierr);

  /* f2=-lambda */
  ierr = FNCreate(PETSC_COMM_WORLD,&f[1]);CHKERRQ(ierr);
  ierr = FNSetType(f[1],FNRATIONAL);CHKERRQ(ierr);
  numer[0] = -1.0; numer[1] = 0.0;
  ierr = FNRationalSetNumerator(f[1],2,numer);CHKERRQ(ierr);

  /* f3=lambda/(lambda-sigma) */
  ierr = FNCreate(PETSC_COMM_WORLD,&f[2]);CHKERRQ(ierr);
  ierr = FNSetType(f[2],FNRATIONAL);CHKERRQ(ierr);
  numer[0] = 1.0; numer[1] = 0.0;
  denom[0] = 1.0; denom[1] = -sigma;
  ierr = FNRationalSetNumerator(f[2],2,numer);CHKERRQ(ierr);
  ierr = FNRationalSetDenominator(f[2],2,denom);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = NEPCreate(PETSC_COMM_WORLD,&nep);CHKERRQ(ierr);
  ierr = NEPSetSplitOperator(nep,3,A,f,SUBSET_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = NEPSetProblemType(nep,NEP_RATIONAL);CHKERRQ(ierr);
  ierr = NEPSetDimensions(nep,8,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  /* set two-sided NLEIGS solver */
  ierr = NEPSetType(nep,NEPNLEIGS);CHKERRQ(ierr);
  ierr = NEPNLEIGSSetFullBasis(nep,PETSC_TRUE);CHKERRQ(ierr);
  ierr = NEPSetTwoSided(nep,PETSC_TRUE);CHKERRQ(ierr);
  ierr = NEPGetRG(nep,&rg);CHKERRQ(ierr);
  ierr = RGSetType(rg,RGINTERVAL);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = RGIntervalSetEndpoints(rg,4.0,700.0,-0.001,0.001);CHKERRQ(ierr);
#else
  ierr = RGIntervalSetEndpoints(rg,4.0,700.0,0,0);CHKERRQ(ierr);
#endif
  ierr = NEPSetTarget(nep,5.0);CHKERRQ(ierr);

  ierr = NEPSetFromOptions(nep);CHKERRQ(ierr);
  ierr = NEPSolve(nep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Check left residual
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatCreateVecs(A[0],&v,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&w);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&z);CHKERRQ(ierr);
  ierr = NEPGetConverged(nep,&nconv);CHKERRQ(ierr);
  ierr = NEPGetTolerances(nep,&tol,NULL);CHKERRQ(ierr);
  for (i=0;i<nconv;i++) {
    ierr = NEPGetEigenpair(nep,i,&lambda,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = NEPGetLeftEigenvector(nep,i,v,NULL);CHKERRQ(ierr);
    ierr = NEPApplyAdjoint(nep,lambda,v,w,r,NULL,NULL);CHKERRQ(ierr);
    ierr = VecNorm(r,NORM_2,&nrm);CHKERRQ(ierr);
    if (nrm>tol*PetscAbsScalar(lambda)) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Left residual i=%D is above tolerance --> %g\n",i,nrm/PetscAbsScalar(lambda));
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Operate with resolvent
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  omega1 = 20.0;
  omega2 = 150.0;
  ierr = VecSet(v,0.0);CHKERRQ(ierr);
  ierr = VecSetValue(v,0,-1.0,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecSetValue(v,1,3.0,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
  ierr = NEPApplyResolvent(nep,NULL,omega1,v,r);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&nrm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"resolvent, omega=%g: norm of computed vector=%g\n",(double)PetscRealPart(omega1),(double)nrm);CHKERRQ(ierr);
  ierr = NEPApplyResolvent(nep,NULL,omega2,v,r);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&nrm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"resolvent, omega=%g: norm of computed vector=%g\n",(double)PetscRealPart(omega2),(double)nrm);CHKERRQ(ierr);
  ierr = VecSet(v,1.0);CHKERRQ(ierr);
  ierr = NEPApplyResolvent(nep,NULL,omega1,v,r);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&nrm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"resolvent, omega=%g: norm of computed vector=%g\n",(double)PetscRealPart(omega1),(double)nrm);CHKERRQ(ierr);
  ierr = NEPApplyResolvent(nep,NULL,omega2,v,r);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&nrm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"resolvent, omega=%g: norm of computed vector=%g\n",(double)PetscRealPart(omega2),(double)nrm);CHKERRQ(ierr);

  /* clean up */
  ierr = NEPDestroy(&nep);CHKERRQ(ierr);
  for (i=0;i<NMAT;i++) {
    ierr = MatDestroy(&A[i]);CHKERRQ(ierr);
    ierr = FNDestroy(&f[i]);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = VecDestroy(&z);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      requires: !single

TEST*/
