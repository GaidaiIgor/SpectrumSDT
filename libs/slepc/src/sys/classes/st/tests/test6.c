/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test STPRECOND operations.\n\n";

#include <slepcst.h>

int main(int argc,char **argv)
{
  Mat            A,P,mat[1];
  ST             st;
  KSP            ksp;
  Vec            v,w;
  PetscScalar    sigma;
  PetscInt       n=10,i,Istart,Iend;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPreconditioner for 1-D Laplacian, n=%D\n\n",n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Compute the operator matrix for the 1-D Laplacian
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
  ierr = MatCreateVecs(A,&v,&w);CHKERRQ(ierr);
  ierr = VecSet(v,1.0);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the spectral transformation object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = STCreate(PETSC_COMM_WORLD,&st);CHKERRQ(ierr);
  mat[0] = A;
  ierr = STSetMatrices(st,1,mat);CHKERRQ(ierr);
  ierr = STSetType(st,STPRECOND);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = STSetFromOptions(st);CHKERRQ(ierr);

  /* set up */
  /* - the transform flag is necessary so that A-sigma*I is built explicitly */
  ierr = STSetTransform(st,PETSC_TRUE);CHKERRQ(ierr);
  /* - the ksphasmat flag is necessary when using STApply(), otherwise can only use PCApply() */
  ierr = STPrecondSetKSPHasMat(st,PETSC_TRUE);CHKERRQ(ierr);
  /* no need to call STSetUp() explicitly */
  ierr = STSetUp(st);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                       Apply the preconditioner
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* default shift */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With default shift\n");CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* change shift */
  sigma = 0.1;
  ierr = STSetShift(st,sigma);CHKERRQ(ierr);
  ierr = STGetShift(st,&sigma);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma));CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);
  ierr = STPostSolve(st);CHKERRQ(ierr);   /* undo changes if inplace */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Test a user-provided preconditioner matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&P);CHKERRQ(ierr);
  ierr = MatSetSizes(P,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(P);CHKERRQ(ierr);
  ierr = MatSetUp(P);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(P,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(P,i,i,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  if (Istart==0) {
    ierr = MatSetValue(P,1,0,-1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(P,0,1,-1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* apply new preconditioner */
  ierr = STPrecondSetMatForPC(st,P);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With user-provided matrix\n");CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                             Clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = STDestroy(&st);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&P);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -st_matmode {{copy inplace}}
      requires: !single

   test:
      suffix: 2
      args: -st_matmode shell

TEST*/
