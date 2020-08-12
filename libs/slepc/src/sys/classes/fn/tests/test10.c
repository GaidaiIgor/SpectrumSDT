/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test Phi functions.\n\n";

#include <slepcfn.h>

/*
   Evaluates phi_k function on a scalar and on a matrix
 */
PetscErrorCode TestPhiFunction(FN fn,PetscScalar x,Mat A,PetscBool verbose)
{
  PetscErrorCode ierr;
  PetscScalar    y,yp;
  char           strx[50],str[50];
  Vec            v,f;
  PetscReal      nrm;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  ierr = FNView(fn,NULL);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(strx,50,x,PETSC_FALSE);CHKERRQ(ierr);
  ierr = FNEvaluateFunction(fn,x,&y);CHKERRQ(ierr);
  ierr = FNEvaluateDerivative(fn,x,&yp);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,y,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nf(%s)=%s\n",strx,str);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,yp,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"f'(%s)=%s\n",strx,str);CHKERRQ(ierr);
  /* compute phi_k(A)*e_1 */
  ierr = MatCreateVecs(A,&v,&f);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMatVec(fn,A,f);CHKERRQ(ierr);  /* reference result by diagonalization */
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMatVec(fn,A,v);CHKERRQ(ierr);
  ierr = VecAXPY(v,-1.0,f);CHKERRQ(ierr);
  ierr = VecNorm(v,NORM_2,&nrm);CHKERRQ(ierr);
  if (nrm>100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-ref is %g\n",(double)nrm);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"f(A)*e_1 =\n");CHKERRQ(ierr);
    ierr = VecView(v,NULL);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&f);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  FN             phi0,phi1,phik,phicopy;
  Mat            A;
  PetscInt       i,j,n=8,k;
  PetscScalar    tau,eta,*As;
  PetscBool      verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test Phi functions, n=%D.\n",n);CHKERRQ(ierr);

  /* Create matrix, fill it with 1-D Laplacian */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A,"A");CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  for (i=0;i<n;i++) As[i+i*n]=2.0;
  j=1;
  for (i=0;i<n-j;i++) { As[i+(i+j)*n]=-1.0; As[(i+j)+i*n]=-1.0; }
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);

  /* phi_0(x) = exp(x) */
  ierr = FNCreate(PETSC_COMM_WORLD,&phi0);CHKERRQ(ierr);
  ierr = FNSetType(phi0,FNPHI);CHKERRQ(ierr);
  ierr = FNPhiSetIndex(phi0,0);CHKERRQ(ierr);
  ierr = TestPhiFunction(phi0,2.2,A,verbose);CHKERRQ(ierr);

  /* phi_1(x) = (exp(x)-1)/x with scaling factors eta*phi_1(tau*x) */
  ierr = FNCreate(PETSC_COMM_WORLD,&phi1);CHKERRQ(ierr);
  ierr = FNSetType(phi1,FNPHI);CHKERRQ(ierr);  /* default index should be 1 */
  tau = 0.2;
  eta = 1.3;
  ierr = FNSetScale(phi1,tau,eta);CHKERRQ(ierr);
  ierr = TestPhiFunction(phi1,2.2,A,verbose);CHKERRQ(ierr);

  /* phi_k(x) with index set from command-line arguments */
  ierr = FNCreate(PETSC_COMM_WORLD,&phik);CHKERRQ(ierr);
  ierr = FNSetType(phik,FNPHI);CHKERRQ(ierr);
  ierr = FNSetFromOptions(phik);CHKERRQ(ierr);

  ierr = FNDuplicate(phik,PETSC_COMM_WORLD,&phicopy);CHKERRQ(ierr);
  ierr = FNPhiGetIndex(phicopy,&k);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Index of phi function is %D\n",k);CHKERRQ(ierr);
  ierr = TestPhiFunction(phicopy,2.2,A,verbose);CHKERRQ(ierr);

  ierr = FNDestroy(&phi0);CHKERRQ(ierr);
  ierr = FNDestroy(&phi1);CHKERRQ(ierr);
  ierr = FNDestroy(&phik);CHKERRQ(ierr);
  ierr = FNDestroy(&phicopy);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -fn_phi_index 3

TEST*/
