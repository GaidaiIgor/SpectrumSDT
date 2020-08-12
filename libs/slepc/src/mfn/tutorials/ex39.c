/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This example illustrates the use of Phi functions in exponential integrators.
   In particular, it implements the Norsett-Euler scheme of stiff order 1.

   The problem is the 1-D heat equation with source term

             y_t = y_xx + 1/(1+u^2) + psi

   where psi is chosen so that the exact solution is yex = x*(1-x)*exp(tend).
   The space domain is [0,1] and the time interval is [0,tend].

       [1] M. Hochbruck and A. Ostermann, "Explicit exponential Runge-Kutta
           methods for semilinear parabolic problems", SIAM J. Numer. Anal. 43(3),
           1069-1090, 2005.
*/

static char help[] = "Exponential integrator for the heat equation with source term.\n\n"
  "The command line options are:\n"
  "  -n <idim>, where <idim> = dimension of the spatial discretization.\n"
  "  -tend <rval>, where <rval> = real value that corresponding to the final time.\n"
  "  -deltat <rval>, where <rval> = real value for the time increment.\n"
  "  -combine <bool>, to represent the phi function with FNCOMBINE instead of FNPHI.\n\n";

#include <slepcmfn.h>

/*
   BuildFNPhi: builds an FNCOMBINE object representing the phi_1 function

        f(x) = (exp(x)-1)/x

   with the following tree:

            f(x)                  f(x)              (combined by division)
           /    \                 p(x) = x          (polynomial)
        a(x)    p(x)              a(x)              (combined by addition)
       /    \                     e(x) = exp(x)     (exponential)
     e(x)   c(x)                  c(x) = -1         (constant)
*/
PetscErrorCode BuildFNPhi(FN fphi)
{
  PetscErrorCode ierr;
  FN             fexp,faux,fconst,fpol;
  PetscScalar    coeffs[2];

  PetscFunctionBeginUser;
  ierr = FNCreate(PETSC_COMM_WORLD,&fexp);CHKERRQ(ierr);
  ierr = FNCreate(PETSC_COMM_WORLD,&fconst);CHKERRQ(ierr);
  ierr = FNCreate(PETSC_COMM_WORLD,&faux);CHKERRQ(ierr);
  ierr = FNCreate(PETSC_COMM_WORLD,&fpol);CHKERRQ(ierr);

  ierr = FNSetType(fexp,FNEXP);CHKERRQ(ierr);

  ierr = FNSetType(fconst,FNRATIONAL);CHKERRQ(ierr);
  coeffs[0] = -1.0;
  ierr = FNRationalSetNumerator(fconst,1,coeffs);CHKERRQ(ierr);

  ierr = FNSetType(faux,FNCOMBINE);CHKERRQ(ierr);
  ierr = FNCombineSetChildren(faux,FN_COMBINE_ADD,fexp,fconst);CHKERRQ(ierr);

  ierr = FNSetType(fpol,FNRATIONAL);CHKERRQ(ierr);
  coeffs[0] = 1.0; coeffs[1] = 0.0;
  ierr = FNRationalSetNumerator(fpol,2,coeffs);CHKERRQ(ierr);

  ierr = FNSetType(fphi,FNCOMBINE);CHKERRQ(ierr);
  ierr = FNCombineSetChildren(fphi,FN_COMBINE_DIVIDE,faux,fpol);CHKERRQ(ierr);

  ierr = FNDestroy(&faux);CHKERRQ(ierr);
  ierr = FNDestroy(&fpol);CHKERRQ(ierr);
  ierr = FNDestroy(&fconst);CHKERRQ(ierr);
  ierr = FNDestroy(&fexp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  Mat               L;
  Vec               u,w,z,yex;
  MFN               mfnexp,mfnphi;
  FN                fexp,fphi;
  PetscBool         combine=PETSC_FALSE;
  PetscInt          i,k,Istart,Iend,n=199,steps;
  PetscReal         t,tend=1.0,deltat=0.01,nrmd,nrmu,x,h;
  const PetscReal   half=0.5;
  PetscScalar       value,c,uval,*warray;
  const PetscScalar *uarray;
  PetscErrorCode    ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-tend",&tend,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-deltat",&deltat,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-combine",&combine,NULL);CHKERRQ(ierr);
  h = 1.0/(n+1.0);
  c = (n+1)*(n+1);

  steps = (PetscInt)(tend/deltat);
  if (PetscAbsReal(tend-steps*deltat)>10*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_WORLD,1,"This example requires tend being a multiple of deltat");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nHeat equation via phi functions, n=%D, tend=%g, deltat=%g%s\n\n",n,(double)tend,(double)deltat,combine?" (combine)":"");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Build the 1-D Laplacian and various vectors
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatCreate(PETSC_COMM_WORLD,&L);CHKERRQ(ierr);
  ierr = MatSetSizes(L,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(L);CHKERRQ(ierr);
  ierr = MatSetUp(L);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(L,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(L,i,i-1,c,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(L,i,i+1,c,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(L,i,i,-2.0*c,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(L,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(L,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateVecs(L,NULL,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&yex);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&w);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&z);CHKERRQ(ierr);

  /*
     Compute various vectors:
     - the exact solution yex = x*(1-x)*exp(tend)
     - the initial condition u = abs(x-0.5)-0.5
  */
  for (i=Istart;i<Iend;i++) {
    x = (i+1)*h;
    value = x*(1.0-x)*PetscExpReal(tend);
    ierr = VecSetValue(yex,i,value,INSERT_VALUES);CHKERRQ(ierr);
    value = PetscAbsReal(x-half)-half;
    ierr = VecSetValue(u,i,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(yex);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(yex);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
  ierr = VecViewFromOptions(yex,NULL,"-exact_sol");CHKERRQ(ierr);
  ierr = VecViewFromOptions(u,NULL,"-initial_cond");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              Create two MFN solvers, for exp() and phi_1()
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MFNCreate(PETSC_COMM_WORLD,&mfnexp);CHKERRQ(ierr);
  ierr = MFNSetOperator(mfnexp,L);CHKERRQ(ierr);
  ierr = MFNGetFN(mfnexp,&fexp);CHKERRQ(ierr);
  ierr = FNSetType(fexp,FNEXP);CHKERRQ(ierr);
  ierr = FNSetScale(fexp,deltat,1.0);CHKERRQ(ierr);
  ierr = MFNSetErrorIfNotConverged(mfnexp,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MFNSetFromOptions(mfnexp);CHKERRQ(ierr);

  ierr = MFNCreate(PETSC_COMM_WORLD,&mfnphi);CHKERRQ(ierr);
  ierr = MFNSetOperator(mfnphi,L);CHKERRQ(ierr);
  ierr = MFNGetFN(mfnphi,&fphi);CHKERRQ(ierr);
  if (combine) {
    ierr = BuildFNPhi(fphi);CHKERRQ(ierr);
  } else {
    ierr = FNSetType(fphi,FNPHI);CHKERRQ(ierr);
    ierr = FNPhiSetIndex(fphi,1);CHKERRQ(ierr);
  }
  ierr = FNSetScale(fphi,deltat,1.0);CHKERRQ(ierr);
  ierr = MFNSetErrorIfNotConverged(mfnphi,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MFNSetFromOptions(mfnphi);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Solve the problem with the Norsett-Euler scheme
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  t = 0.0;
  for (k=0;k<steps;k++) {

    /* evaluate nonlinear part */
    ierr = VecGetArrayRead(u,&uarray);CHKERRQ(ierr);
    ierr = VecGetArray(w,&warray);CHKERRQ(ierr);
    for (i=Istart;i<Iend;i++) {
      x = (i+1)*h;
      uval = uarray[i-Istart];
      value = x*(1.0-x)*PetscExpReal(t);
      value = value + 2.0*PetscExpReal(t) - 1.0/(1.0+value*value);
      value = value + 1.0/(1.0+uval*uval);
      warray[i-Istart] = deltat*value;
    }
    ierr = VecRestoreArrayRead(u,&uarray);CHKERRQ(ierr);
    ierr = VecRestoreArray(w,&warray);CHKERRQ(ierr);
    ierr = MFNSolve(mfnphi,w,z);CHKERRQ(ierr);

    /* evaluate linear part */
    ierr = MFNSolve(mfnexp,u,u);CHKERRQ(ierr);
    ierr = VecAXPY(u,1.0,z);CHKERRQ(ierr);
    t = t + deltat;

  }
  ierr = VecViewFromOptions(u,NULL,"-computed_sol");CHKERRQ(ierr);

  /*
     Compare with exact solution and show error norm
  */
  ierr = VecCopy(u,z);CHKERRQ(ierr);
  ierr = VecAXPY(z,-1.0,yex);CHKERRQ(ierr);
  ierr = VecNorm(z,NORM_2,&nrmd);CHKERRQ(ierr);
  ierr = VecNorm(u,NORM_2,&nrmu);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," The relative error at t=%g is %.4f\n\n",(double)t,(double)(nrmd/nrmu));CHKERRQ(ierr);

  /*
     Free work space
  */
  ierr = MFNDestroy(&mfnexp);CHKERRQ(ierr);
  ierr = MFNDestroy(&mfnphi);CHKERRQ(ierr);
  ierr = MatDestroy(&L);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&yex);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = VecDestroy(&z);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   build:
      requires: c99

   test:
      suffix: 1
      args: -n 127 -tend 0.125 -mfn_tol 1e-3 -deltat 0.025
      timeoutfactor: 2

   test:
      suffix: 2
      args: -n 127 -tend 0.125 -mfn_tol 1e-3 -deltat 0.025 -combine
      filter: sed -e "s/ (combine)//"
      requires: !single
      output_file: output/ex39_1.out
      timeoutfactor: 2

TEST*/
