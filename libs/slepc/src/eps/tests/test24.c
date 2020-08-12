/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Eigenproblem for the 1-D Laplacian with constraints. "
  "Based on ex1.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

#include <slepceps.h>

int main(int argc,char **argv)
{
  Mat            A;
  EPS            eps;
  EPSType        type;
  Vec            *vi=NULL,*vc=NULL,t;
  PetscInt       n=30,nev=4,i,j,Istart,Iend,nini=0,ncon=0,bs;
  PetscReal      alpha,beta,restart;
  PetscBool      flg,lock;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-nini",&nini,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ncon",&ncon,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem, n=%D nini=%D ncon=%D\n\n",n,nini,ncon);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(A,i,i-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A,i,i+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,i,i,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  ierr = EPSSetType(eps,EPSLOBPCG);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
  ierr = EPSSetConvergenceTest(eps,EPS_CONV_ABS);CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps,nev,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = EPSLOBPCGSetBlockSize(eps,nev);CHKERRQ(ierr);
  ierr = EPSLOBPCGSetRestart(eps,0.7);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,1e-8,1200);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  ierr = MatCreateVecs(A,&t,NULL);CHKERRQ(ierr);
  if (nini) {
    ierr = VecDuplicateVecs(t,nini,&vi);CHKERRQ(ierr);
    for (i=0;i<nini;i++) {
      ierr = VecSetRandom(vi[i],NULL);CHKERRQ(ierr);
    }
    ierr = EPSSetInitialSpace(eps,nini,vi);CHKERRQ(ierr);
  }
  if (ncon) {   /* constraints are exact eigenvectors of lowest eigenvalues */
    alpha = PETSC_PI/(n+1);
    beta  = PetscSqrtReal(2.0/(n+1));
    ierr = VecDuplicateVecs(t,ncon,&vc);CHKERRQ(ierr);
    for (i=0;i<ncon;i++) {
      for (j=0;j<n;j++) {
        ierr = VecSetValue(vc[i],j,PetscSinReal(alpha*(j+1)*(i+1))*beta,INSERT_VALUES);CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(vc[i]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(vc[i]);CHKERRQ(ierr);
    }
    ierr = EPSSetDeflationSpace(eps,ncon,vc);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)eps,EPSLOBPCG,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EPSLOBPCGGetLocking(eps,&lock);CHKERRQ(ierr);
    if (lock) { ierr = PetscPrintf(PETSC_COMM_WORLD," Using soft locking\n");CHKERRQ(ierr); }
    ierr = EPSLOBPCGGetRestart(eps,&restart);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," LOBPCG Restart parameter=%.4g\n",(double)restart);CHKERRQ(ierr);
    ierr = EPSLOBPCGGetBlockSize(eps,&bs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," LOBPCG Block size=%D\n",bs);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroyVecs(nini,&vi);CHKERRQ(ierr);
  ierr = VecDestroyVecs(ncon,&vc);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -ncon 2
      output_file: output/test24_1.out
      test:
         suffix: 1
         requires: !single
      test:
         suffix: 2_cuda
         args: -mat_type aijcusparse
         requires: cuda !single

TEST*/
