/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV operations with non-standard inner product.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  Vec            t,v;
  Mat            B,M;
  BV             X;
  PetscInt       i,j,n=10,k=5,Istart,Iend;
  PetscScalar    alpha;
  PetscReal      nrm;
  PetscViewer    view;
  PetscBool      verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test BV with non-standard inner product (n=%D, k=%D).\n",n,k);CHKERRQ(ierr);

  /* Create inner product matrix */
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

  /* Create BV object X */
  ierr = BVCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)X,"X");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(X,t,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(X);CHKERRQ(ierr);
  ierr = BVSetMatrix(X,B,PETSC_FALSE);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Fill X entries */
  for (j=0;j<k;j++) {
    ierr = BVGetColumn(X,j,&v);CHKERRQ(ierr);
    ierr = VecSet(v,0.0);CHKERRQ(ierr);
    for (i=0;i<4;i++) {
      if (i+j<n) {
        ierr = VecSetValue(v,i+j,(PetscScalar)(3*i+j-2),INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
    ierr = BVRestoreColumn(X,j,&v);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = MatView(B,view);CHKERRQ(ierr);
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Test BVNormColumn */
  ierr = BVNormColumn(X,0,NORM_2,&nrm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"B-Norm of X[0] = %g\n",(double)nrm);CHKERRQ(ierr);

  /* Test BVOrthogonalizeColumn */
  for (j=0;j<k;j++) {
    ierr = BVOrthogonalizeColumn(X,j,NULL,&nrm,NULL);CHKERRQ(ierr);
    alpha = 1.0/nrm;
    ierr = BVScaleColumn(X,j,alpha);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Check orthogonality */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M);CHKERRQ(ierr);
  ierr = BVDot(X,X,M);CHKERRQ(ierr);
  ierr = MatShift(M,-1.0);CHKERRQ(ierr);
  ierr = MatNorm(M,NORM_1,&nrm);CHKERRQ(ierr);
  if (nrm<100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality < 100*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)nrm);CHKERRQ(ierr);
  }

  /* Test BVNormVecBegin/End */
  ierr = BVGetColumn(X,0,&v);CHKERRQ(ierr);
  ierr = BVNormVecBegin(X,v,NORM_1,&nrm);CHKERRQ(ierr);
  ierr = BVNormVecEnd(X,v,NORM_1,&nrm);CHKERRQ(ierr);
  ierr = BVRestoreColumn(X,0,&v);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"B-Norm of X[0] = %g\n",(double)nrm);CHKERRQ(ierr);

  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test3_1.out

   test:
      suffix: 1_svec_vecs
      nsize: 1
      args: -bv_type svec -bv_matmult vecs
      output_file: output/test3_1.out

   test:
      suffix: 1_cuda
      nsize: 1
      args: -bv_type svec -mat_type aijcusparse
      requires: cuda
      output_file: output/test3_1.out

   test:
      suffix: 2
      nsize: 2
      args: -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test3_1.out

   test:
      suffix: 3
      nsize: 2
      args: -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_type mgs
      output_file: output/test3_1.out

TEST*/
