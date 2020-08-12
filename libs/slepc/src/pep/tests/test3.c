/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test PEP interface functions.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat                A[3],B;      /* problem matrices */
  PEP                pep;         /* eigenproblem solver context */
  ST                 st;
  KSP                ksp;
  DS                 ds;
  PetscReal          tol,alpha;
  PetscScalar        target;
  PetscInt           n=20,i,its,nev,ncv,mpd,Istart,Iend,nmat;
  PEPWhich           which;
  PEPConvergedReason reason;
  PEPType            type;
  PEPExtract         extr;
  PEPBasis           basis;
  PEPScale           scale;
  PEPRefine          refine;
  PEPRefineScheme    rscheme;
  PEPConv            conv;
  PEPStop            stop;
  PEPProblemType     ptype;
  PetscErrorCode     ierr;
  PetscViewerAndFormat *vf;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nDiagonal Quadratic Eigenproblem, n=%D\n\n",n);CHKERRQ(ierr);

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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Create eigensolver and test interface functions
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  ierr = PEPSetOperators(pep,3,A);CHKERRQ(ierr);
  ierr = PEPGetNumMatrices(pep,&nmat);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Polynomial of degree %d\n",(int)nmat-1);CHKERRQ(ierr);
  ierr = PEPGetOperators(pep,0,&B);CHKERRQ(ierr);
  ierr = MatView(B,NULL);CHKERRQ(ierr);

  ierr = PEPSetType(pep,PEPTOAR);CHKERRQ(ierr);
  ierr = PEPGetType(pep,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Type set to %s\n",type);CHKERRQ(ierr);

  ierr = PEPGetProblemType(pep,&ptype);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Problem type before changing = %d",(int)ptype);CHKERRQ(ierr);
  ierr = PEPSetProblemType(pep,PEP_HERMITIAN);CHKERRQ(ierr);
  ierr = PEPGetProblemType(pep,&ptype);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d.\n",(int)ptype);CHKERRQ(ierr);

  ierr = PEPGetExtract(pep,&extr);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Extraction before changing = %d",(int)extr);CHKERRQ(ierr);
  ierr = PEPSetExtract(pep,PEP_EXTRACT_STRUCTURED);CHKERRQ(ierr);
  ierr = PEPGetExtract(pep,&extr);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)extr);CHKERRQ(ierr);

  ierr = PEPSetScale(pep,PEP_SCALE_SCALAR,.1,NULL,NULL,5,1.0);CHKERRQ(ierr);
  ierr = PEPGetScale(pep,&scale,&alpha,NULL,NULL,&its,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Scaling: %s, alpha=%g, its=%D\n",PEPScaleTypes[scale],(double)alpha,its);CHKERRQ(ierr);

  ierr = PEPSetBasis(pep,PEP_BASIS_CHEBYSHEV1);CHKERRQ(ierr);
  ierr = PEPGetBasis(pep,&basis);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Polynomial basis: %s\n",PEPBasisTypes[basis]);CHKERRQ(ierr);

  ierr = PEPSetRefine(pep,PEP_REFINE_SIMPLE,1,1e-9,2,PEP_REFINE_SCHEME_SCHUR);CHKERRQ(ierr);
  ierr = PEPGetRefine(pep,&refine,NULL,&tol,&its,&rscheme);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Refinement: %s, tol=%g, its=%D, scheme=%s\n",PEPRefineTypes[refine],(double)tol,its,PEPRefineSchemes[rscheme]);CHKERRQ(ierr);

  ierr = PEPSetTarget(pep,4.8);CHKERRQ(ierr);
  ierr = PEPGetTarget(pep,&target);CHKERRQ(ierr);
  ierr = PEPSetWhichEigenpairs(pep,PEP_TARGET_MAGNITUDE);CHKERRQ(ierr);
  ierr = PEPGetWhichEigenpairs(pep,&which);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Which = %d, target = %g\n",(int)which,(double)PetscRealPart(target));CHKERRQ(ierr);

  ierr = PEPSetDimensions(pep,4,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = PEPGetDimensions(pep,&nev,&ncv,&mpd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Dimensions: nev=%D, ncv=%D, mpd=%D\n",nev,ncv,mpd);CHKERRQ(ierr);

  ierr = PEPSetTolerances(pep,2.2e-4,200);CHKERRQ(ierr);
  ierr = PEPGetTolerances(pep,&tol,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Tolerance = %.5f, max_its = %D\n",(double)tol,its);CHKERRQ(ierr);

  ierr = PEPSetConvergenceTest(pep,PEP_CONV_ABS);CHKERRQ(ierr);
  ierr = PEPGetConvergenceTest(pep,&conv);CHKERRQ(ierr);
  ierr = PEPSetStoppingTest(pep,PEP_STOP_BASIC);CHKERRQ(ierr);
  ierr = PEPGetStoppingTest(pep,&stop);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Convergence test = %d, stopping test = %d\n",(int)conv,(int)stop);CHKERRQ(ierr);

  ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);CHKERRQ(ierr);
  ierr = PEPMonitorSet(pep,(PetscErrorCode (*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))PEPMonitorFirst,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);
  ierr = PEPMonitorCancel(pep);CHKERRQ(ierr);

  ierr = PEPGetST(pep,&st);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1e-8,1e-35,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = STView(st,NULL);CHKERRQ(ierr);
  ierr = PEPGetDS(pep,&ds);CHKERRQ(ierr);
  ierr = DSView(ds,NULL);CHKERRQ(ierr);

  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);
  ierr = PEPSolve(pep);CHKERRQ(ierr);
  ierr = PEPGetConvergedReason(pep,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Finished - converged reason = %d\n",(int)reason);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PEPErrorView(pep,PEP_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  ierr = MatDestroy(&A[0]);CHKERRQ(ierr);
  ierr = MatDestroy(&A[1]);CHKERRQ(ierr);
  ierr = MatDestroy(&A[2]);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -pep_tol 1e-6 -pep_ncv 22

TEST*/
