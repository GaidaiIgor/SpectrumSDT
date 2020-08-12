/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV matrix projection.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  Vec            t,v;
  Mat            B,G,H0,H1;
  BV             X,Y,Z;
  PetscInt       i,j,n=20,kx=6,lx=3,ky=5,ly=2,Istart,Iend,col[5];
  PetscScalar    alpha,value[] = { -1, 1, 1, 1, 1 };
  PetscViewer    view;
  PetscReal      norm;
  PetscBool      verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-kx",&kx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-lx",&lx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ky",&ky,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ly",&ly,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test BV projection (n=%D).\n",n);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"X has %D active columns (%D leading columns).\n",kx,lx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Y has %D active columns (%D leading columns).\n",ky,ly);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Create non-symmetric matrix G (Toeplitz) */
  ierr = MatCreate(PETSC_COMM_WORLD,&G);CHKERRQ(ierr);
  ierr = MatSetSizes(G,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(G);CHKERRQ(ierr);
  ierr = MatSetUp(G);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)G,"G");CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(G,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    col[0]=i-1; col[1]=i; col[2]=i+1; col[3]=i+2; col[4]=i+3;
    if (i==0) {
      ierr = MatSetValues(G,1,&i,PetscMin(4,n-i),col+1,value+1,INSERT_VALUES);CHKERRQ(ierr);
    } else {
      ierr = MatSetValues(G,1,&i,PetscMin(5,n-i+1),col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(G,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(G,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (verbose) {
    ierr = MatView(G,view);CHKERRQ(ierr);
  }

  /* Create symmetric matrix B (1-D Laplacian) */
  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(B);CHKERRQ(ierr);
  ierr = MatSetUp(B);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)B,"B");CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(B,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(B,i,i-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(B,i,i+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(B,i,i,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateVecs(B,&t,NULL);CHKERRQ(ierr);
  if (verbose) {
    ierr = MatView(B,view);CHKERRQ(ierr);
  }

  /* Create BV object X */
  ierr = BVCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)X,"X");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(X,t,kx+2);CHKERRQ(ierr);  /* two extra columns to test active columns */
  ierr = BVSetFromOptions(X);CHKERRQ(ierr);

  /* Fill X entries */
  for (j=0;j<kx+2;j++) {
    ierr = BVGetColumn(X,j,&v);CHKERRQ(ierr);
    ierr = VecSet(v,0.0);CHKERRQ(ierr);
    for (i=0;i<4;i++) {
      if (i+j<n) {
#if defined(PETSC_USE_COMPLEX)
        alpha = PetscCMPLX((PetscReal)(3*i+j-2),(PetscReal)(2*i));
#else
        alpha = (PetscReal)(3*i+j-2);
#endif
        ierr = VecSetValue(v,i+j,alpha,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
    ierr = BVRestoreColumn(X,j,&v);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Duplicate BV object and store Z=G*X */
  ierr = BVDuplicate(X,&Z);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)Z,"Z");CHKERRQ(ierr);
  ierr = BVSetActiveColumns(X,0,kx);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(Z,0,kx);CHKERRQ(ierr);
  ierr = BVMatMult(X,G,Z);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(X,lx,kx);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(Z,lx,kx);CHKERRQ(ierr);

  /* Create BV object Y */
  ierr = BVCreate(PETSC_COMM_WORLD,&Y);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)Y,"Y");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(Y,t,ky+1);CHKERRQ(ierr);
  ierr = BVSetFromOptions(Y);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(Y,ly,ky);CHKERRQ(ierr);

  /* Fill Y entries */
  for (j=0;j<ky+1;j++) {
    ierr = BVGetColumn(Y,j,&v);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    alpha = PetscCMPLX((PetscReal)(j+1)/4.0,-(PetscReal)j);
#else
    alpha = (PetscReal)(j+1)/4.0;
#endif
    ierr = VecSet(v,(PetscScalar)(j+1)/4.0);CHKERRQ(ierr);
    ierr = BVRestoreColumn(Y,j,&v);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = BVView(Y,view);CHKERRQ(ierr);
  }

  /* Test BVMatProject for non-symmetric matrix G */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,ky,kx,NULL,&H0);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)H0,"H0");CHKERRQ(ierr);
  ierr = BVMatProject(X,G,Y,H0);CHKERRQ(ierr);
  if (verbose) {
    ierr = MatView(H0,view);CHKERRQ(ierr);
  }

  /* Test BVMatProject with previously stored G*X */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,ky,kx,NULL,&H1);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)H1,"H1");CHKERRQ(ierr);
  ierr = BVMatProject(Z,NULL,Y,H1);CHKERRQ(ierr);
  if (verbose) {
    ierr = MatView(H1,view);CHKERRQ(ierr);
  }

  /* Check that H0 and H1 are equal */
  ierr = MatAXPY(H0,-1.0,H1,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatNorm(H0,NORM_1,&norm);CHKERRQ(ierr);
  if (norm<10*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||H0-H1|| < 10*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||H0-H1||=%g\n",(double)norm);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&H0);CHKERRQ(ierr);
  ierr = MatDestroy(&H1);CHKERRQ(ierr);

  /* Test BVMatProject for symmetric matrix B with orthogonal projection */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,kx,kx,NULL,&H0);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)H0,"H0");CHKERRQ(ierr);
  ierr = BVMatProject(X,B,X,H0);CHKERRQ(ierr);
  if (verbose) {
    ierr = MatView(H0,view);CHKERRQ(ierr);
  }

  /* Repeat previous test with symmetry flag set */
  ierr = MatSetOption(B,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,kx,kx,NULL,&H1);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)H1,"H1");CHKERRQ(ierr);
  ierr = BVMatProject(X,B,X,H1);CHKERRQ(ierr);
  if (verbose) {
    ierr = MatView(H1,view);CHKERRQ(ierr);
  }

  /* Check that H0 and H1 are equal */
  ierr = MatAXPY(H0,-1.0,H1,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatNorm(H0,NORM_1,&norm);CHKERRQ(ierr);
  if (norm<10*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||H0-H1|| < 10*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||H0-H1||=%g\n",(double)norm);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&H0);CHKERRQ(ierr);
  ierr = MatDestroy(&H1);CHKERRQ(ierr);

  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = BVDestroy(&Y);CHKERRQ(ierr);
  ierr = BVDestroy(&Z);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = MatDestroy(&G);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test9_1.out

   test:
      suffix: 1_svec_vecs
      nsize: 1
      args: -bv_type svec -bv_matmult vecs
      output_file: output/test9_1.out

   test:
      suffix: 1_cuda
      nsize: 1
      args: -bv_type svec -mat_type aijcusparse
      requires: cuda
      output_file: output/test9_1.out

   test:
      suffix: 2
      nsize: 2
      args: -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test9_1.out

   test:
      suffix: 2_svec_vecs
      nsize: 2
      args: -bv_type svec -bv_matmult vecs
      output_file: output/test9_1.out

TEST*/
