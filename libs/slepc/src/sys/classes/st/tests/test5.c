/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test various ST interface functions.\n\n";

#include <slepceps.h>

int main (int argc,char **argv)
{
  ST             st;
  KSP            ksp;
  PC             pc;
  Mat            A,mat[1],Op;
  Vec            v,w;
  PetscInt       N,n=4,i,j,II,Istart,Iend;
  PetscScalar    d;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  N = n*n;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create non-symmetric matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    d = 0.0;
    if (i>0) { ierr = MatSetValue(A,II,II-n,1.0,INSERT_VALUES);CHKERRQ(ierr); d=d+1.0; }
    if (i<n-1) { ierr = MatSetValue(A,II,II+n,-1.0,INSERT_VALUES);CHKERRQ(ierr); d=d+1.0; }
    if (j>0) { ierr = MatSetValue(A,II,II-1,1.0,INSERT_VALUES);CHKERRQ(ierr); d=d+1.0; }
    if (j<n-1) { ierr = MatSetValue(A,II,II+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); d=d+1.0; }
    ierr = MatSetValue(A,II,II,d,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,&v,&w);CHKERRQ(ierr);
  ierr = VecSetValue(v,0,-.5,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecSetValue(v,1,1.5,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecSetValue(v,2,2,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
  ierr = VecView(v,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create the spectral transformation object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = STCreate(PETSC_COMM_WORLD,&st);CHKERRQ(ierr);
  mat[0] = A;
  ierr = STSetMatrices(st,1,mat);CHKERRQ(ierr);
  ierr = STSetType(st,STCAYLEY);CHKERRQ(ierr);
  ierr = STSetShift(st,2.0);CHKERRQ(ierr);
  ierr = STCayleySetAntishift(st,1.0);CHKERRQ(ierr);
  ierr = STSetTransform(st,PETSC_TRUE);CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,100*PETSC_MACHINE_EPSILON,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = STSetKSP(st,ksp);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  ierr = STSetFromOptions(st);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Apply the operator
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  ierr = STApplyTranspose(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  ierr = STMatMult(st,1,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  ierr = STMatMultTranspose(st,1,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  ierr = STMatSolve(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  ierr = STMatSolveTranspose(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Get the operator matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = STGetOperator(st,&Op);CHKERRQ(ierr);
  ierr = MatMult(Op,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);
  ierr = MatMultTranspose(Op,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);
  ierr = STRestoreOperator(st,&Op);CHKERRQ(ierr);

  ierr = STDestroy(&st);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      output_file: output/test5_1.out
      requires: !single
      test:
         args: -st_matmode {{copy inplace}}
      test:
         args: -st_matmode shell -ksp_type bcgs -pc_type jacobi

TEST*/
