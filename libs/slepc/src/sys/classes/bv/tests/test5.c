/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV operations with indefinite inner product.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  Vec            t,v,w,omega;
  Mat            B,M;
  BV             X,Y;
  PetscInt       i,j,n=10,k=5,l,Istart,Iend;
  PetscScalar    alpha;
  PetscReal      nrm;
  PetscViewer    view;
  PetscBool      verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test BV with indefinite inner product (n=%D, k=%D).\n",n,k);CHKERRQ(ierr);

  /* Create inner product matrix (standard involutionary permutation) */
  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(B);CHKERRQ(ierr);
  ierr = MatSetUp(B);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)B,"B");CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(B,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(B,i,n-i-1,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateVecs(B,&t,NULL);CHKERRQ(ierr);

  /* Create BV object X */
  ierr = BVCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)X,"X");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(X,t,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(X);CHKERRQ(ierr);
  ierr = BVSetMatrix(X,B,PETSC_TRUE);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Fill X entries */
  l = -3;
  for (j=0;j<k;j++) {
    ierr = BVGetColumn(X,j,&v);CHKERRQ(ierr);
    ierr = VecSet(v,-1.0);CHKERRQ(ierr);
    for (i=0;i<n/2;i++) {
      if (i+j<n) {
        l = (l + 3*i+j-2) % n;
        ierr = VecSetValue(v,i+j,(PetscScalar)l,INSERT_VALUES);CHKERRQ(ierr);
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

  /* Create a copy on Y */
  ierr = BVDuplicate(X,&Y);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)Y,"Y");CHKERRQ(ierr);
  ierr = BVCopy(X,Y);CHKERRQ(ierr);

  /* Check orthogonality */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M);CHKERRQ(ierr);
  ierr = BVDot(Y,Y,M);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF,k,&omega);CHKERRQ(ierr);
  ierr = BVGetSignature(Y,omega);CHKERRQ(ierr);
  ierr = VecScale(omega,-1.0);CHKERRQ(ierr);
  ierr = MatDiagonalSet(M,omega,ADD_VALUES);CHKERRQ(ierr);
  ierr = MatNorm(M,NORM_1,&nrm);CHKERRQ(ierr);
  if (nrm<100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality < 100*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)nrm);CHKERRQ(ierr);
  }

  /* Test BVSetSignature */
  ierr = VecScale(omega,-1.0);CHKERRQ(ierr);
  ierr = BVSetSignature(Y,omega);CHKERRQ(ierr);
  ierr = VecDestroy(&omega);CHKERRQ(ierr);

  /* Test BVApplyMatrix */
  ierr = VecDuplicate(t,&w);CHKERRQ(ierr);
  ierr = BVGetColumn(X,0,&v);CHKERRQ(ierr);
  ierr = BVApplyMatrix(X,v,w);CHKERRQ(ierr);
  ierr = BVApplyMatrix(X,w,t);CHKERRQ(ierr);
  ierr = VecAXPY(t,-1.0,v);CHKERRQ(ierr);
  ierr = BVRestoreColumn(X,0,&v);CHKERRQ(ierr);
  ierr = VecNorm(t,NORM_2,&nrm);CHKERRQ(ierr);
  if (PetscAbsReal(nrm)>10*PETSC_MACHINE_EPSILON) SETERRQ1(PETSC_COMM_WORLD,1,"Wrong value, nrm = %g\n",(double)nrm);

  ierr = BVApplyMatrixBV(X,Y);CHKERRQ(ierr);
  ierr = BVGetColumn(Y,0,&v);CHKERRQ(ierr);
  ierr = VecAXPY(w,-1.0,v);CHKERRQ(ierr);
  ierr = BVRestoreColumn(Y,0,&v);CHKERRQ(ierr);
  ierr = VecNorm(w,NORM_2,&nrm);CHKERRQ(ierr);
  if (PetscAbsReal(nrm)>10*PETSC_MACHINE_EPSILON) SETERRQ1(PETSC_COMM_WORLD,1,"Wrong value, nrm = %g\n",(double)nrm);

  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = BVDestroy(&Y);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -bv_orthog_refine always -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test5_1.out

   test:
      suffix: 1_cuda
      nsize: 1
      args: -bv_orthog_refine always -bv_type svec -mat_type aijcusparse
      requires: cuda
      output_file: output/test5_1.out

   test:
      suffix: 2
      nsize: 1
      args: -bv_orthog_refine always -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_type mgs
      output_file: output/test5_1.out

   test:
      suffix: 2_cuda
      nsize: 1
      args: -bv_orthog_refine always -bv_type svec -mat_type aijcusparse -bv_orthog_type mgs
      requires: cuda
      output_file: output/test5_1.out


TEST*/
