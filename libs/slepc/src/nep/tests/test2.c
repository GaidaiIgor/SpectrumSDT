/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test NEP interface functions.\n\n";

#include <slepcnep.h>

int main(int argc,char **argv)
{
  Mat                  A[3],B;      /* problem matrices */
  FN                   f[3],g;      /* problem functions */
  NEP                  nep;         /* eigenproblem solver context */
  DS                   ds;
  RG                   rg;
  PetscReal            tol;
  PetscScalar          coeffs[2],target;
  PetscInt             n=20,i,its,nev,ncv,mpd,Istart,Iend,nterm;
  PetscBool            twoside;
  NEPWhich             which;
  NEPConvergedReason   reason;
  NEPType              type;
  NEPRefine            refine;
  NEPRefineScheme      rscheme;
  NEPConv              conv;
  NEPStop              stop;
  NEPProblemType       ptype;
  MatStructure         mstr;
  PetscErrorCode       ierr;
  PetscViewerAndFormat *vf;
  const char           *str=NULL;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nDiagonal Nonlinear Eigenproblem, n=%D\n\n",n);CHKERRQ(ierr);

  /*
     Matrices
  */
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
  ierr = MatGetOwnershipRange(A[1],&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(A[1],i,i,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A[2]);CHKERRQ(ierr);
  ierr = MatSetSizes(A[2],PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A[2]);CHKERRQ(ierr);
  ierr = MatSetUp(A[2]);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A[1],&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(A[2],i,i,n/(PetscReal)(i+1),INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A[2],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A[2],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /*
     Functions: f0=-lambda, f1=1.0, f2=sqrt(lambda)
  */
  ierr = FNCreate(PETSC_COMM_WORLD,&f[0]);CHKERRQ(ierr);
  ierr = FNSetType(f[0],FNRATIONAL);CHKERRQ(ierr);
  coeffs[0] = -1.0; coeffs[1] = 0.0;
  ierr = FNRationalSetNumerator(f[0],2,coeffs);CHKERRQ(ierr);

  ierr = FNCreate(PETSC_COMM_WORLD,&f[1]);CHKERRQ(ierr);
  ierr = FNSetType(f[1],FNRATIONAL);CHKERRQ(ierr);
  coeffs[0] = 1.0;
  ierr = FNRationalSetNumerator(f[1],1,coeffs);CHKERRQ(ierr);

  ierr = FNCreate(PETSC_COMM_WORLD,&f[2]);CHKERRQ(ierr);
  ierr = FNSetType(f[2],FNSQRT);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Create eigensolver and test interface functions
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = NEPCreate(PETSC_COMM_WORLD,&nep);CHKERRQ(ierr);
  ierr = NEPSetSplitOperator(nep,3,A,f,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = NEPGetSplitOperatorInfo(nep,&nterm,&mstr);CHKERRQ(ierr);
  switch (mstr) {
    case DIFFERENT_NONZERO_PATTERN: str = "different"; break;
    case SUBSET_NONZERO_PATTERN:    str = "subset"; break;
    case SAME_NONZERO_PATTERN:      str = "same"; break;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD," Nonlinear function with %d terms, with %s nonzero pattern\n",(int)nterm,str);CHKERRQ(ierr);
  ierr = NEPGetSplitOperatorTerm(nep,0,&B,&g);CHKERRQ(ierr);
  ierr = MatView(B,NULL);CHKERRQ(ierr);
  ierr = FNView(g,NULL);CHKERRQ(ierr);

  ierr = NEPSetType(nep,NEPRII);CHKERRQ(ierr);
  ierr = NEPGetType(nep,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Type set to %s\n",type);CHKERRQ(ierr);
  ierr = NEPGetTwoSided(nep,&twoside);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Two-sided flag = %s\n",twoside?"true":"false");CHKERRQ(ierr);

  ierr = NEPGetProblemType(nep,&ptype);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Problem type before changing = %d",(int)ptype);CHKERRQ(ierr);
  ierr = NEPSetProblemType(nep,NEP_RATIONAL);CHKERRQ(ierr);
  ierr = NEPGetProblemType(nep,&ptype);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d.\n",(int)ptype);CHKERRQ(ierr);

  ierr = NEPSetRefine(nep,NEP_REFINE_SIMPLE,1,1e-9,2,NEP_REFINE_SCHEME_EXPLICIT);CHKERRQ(ierr);
  ierr = NEPGetRefine(nep,&refine,NULL,&tol,&its,&rscheme);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Refinement: %s, tol=%g, its=%D, scheme=%s\n",NEPRefineTypes[refine],(double)tol,its,NEPRefineSchemes[rscheme]);CHKERRQ(ierr);

  ierr = NEPSetTarget(nep,1.1);CHKERRQ(ierr);
  ierr = NEPGetTarget(nep,&target);CHKERRQ(ierr);
  ierr = NEPSetWhichEigenpairs(nep,NEP_TARGET_MAGNITUDE);CHKERRQ(ierr);
  ierr = NEPGetWhichEigenpairs(nep,&which);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Which = %d, target = %g\n",(int)which,(double)PetscRealPart(target));CHKERRQ(ierr);

  ierr = NEPSetDimensions(nep,1,12,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = NEPGetDimensions(nep,&nev,&ncv,&mpd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Dimensions: nev=%D, ncv=%D, mpd=%D\n",nev,ncv,mpd);CHKERRQ(ierr);

  ierr = NEPSetTolerances(nep,1.0e-6,200);CHKERRQ(ierr);
  ierr = NEPGetTolerances(nep,&tol,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Tolerance = %.6f, max_its = %D\n",(double)tol,its);CHKERRQ(ierr);

  ierr = NEPSetConvergenceTest(nep,NEP_CONV_ABS);CHKERRQ(ierr);
  ierr = NEPGetConvergenceTest(nep,&conv);CHKERRQ(ierr);
  ierr = NEPSetStoppingTest(nep,NEP_STOP_BASIC);CHKERRQ(ierr);
  ierr = NEPGetStoppingTest(nep,&stop);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Convergence test = %d, stopping test = %d\n",(int)conv,(int)stop);CHKERRQ(ierr);

  ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);CHKERRQ(ierr);
  ierr = NEPMonitorSet(nep,(PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))NEPMonitorFirst,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);
  ierr = NEPMonitorCancel(nep);CHKERRQ(ierr);

  ierr = NEPGetDS(nep,&ds);CHKERRQ(ierr);
  ierr = DSView(ds,NULL);CHKERRQ(ierr);
  ierr = NEPSetFromOptions(nep);CHKERRQ(ierr);

  ierr = NEPGetRG(nep,&rg);CHKERRQ(ierr);
  ierr = RGView(rg,NULL);CHKERRQ(ierr);

  ierr = NEPSolve(nep);CHKERRQ(ierr);
  ierr = NEPGetConvergedReason(nep,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Finished - converged reason = %d\n",(int)reason);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = NEPErrorView(nep,NEP_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
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
      requires: !single

TEST*/
