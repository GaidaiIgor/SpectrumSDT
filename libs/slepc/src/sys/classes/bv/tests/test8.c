/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV orthogonalization with selected columns.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  BV             X;
  Vec            v,t,z;
  PetscInt       i,j,n=20,k=8;
  PetscViewer    view;
  PetscBool      verbose,*which;
  PetscReal      norm;
  PetscScalar    alpha,*pz;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test BV orthogonalization with selected columns of length %D.\n",n);CHKERRQ(ierr);

  /* Create template vector */
  ierr = VecCreate(PETSC_COMM_WORLD,&t);CHKERRQ(ierr);
  ierr = VecSetSizes(t,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(t);CHKERRQ(ierr);

  /* Create BV object X */
  ierr = BVCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)X,"X");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(X,t,k);CHKERRQ(ierr);
  ierr = BVSetOrthogonalization(X,BV_ORTHOG_MGS,BV_ORTHOG_REFINE_IFNEEDED,PETSC_DEFAULT,BV_ORTHOG_BLOCK_GS);CHKERRQ(ierr);
  ierr = BVSetFromOptions(X);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Fill X entries */
  for (j=0;j<k;j++) {
    ierr = BVGetColumn(X,j,&v);CHKERRQ(ierr);
    ierr = VecSet(v,0.0);CHKERRQ(ierr);
    for (i=0;i<=n/2;i++) {
      if (i+j<n) {
        alpha = (3.0*i+j-2)/(2*(i+j+1));
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

  /* Orthonormalize first k-1 columns */
  for (j=0;j<k-1;j++) {
    ierr = BVOrthogonalizeColumn(X,j,NULL,&norm,NULL);CHKERRQ(ierr);
    alpha = 1.0/norm;
    ierr = BVScaleColumn(X,j,alpha);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Select odd columns and orthogonalize last column against those only */
  ierr = PetscMalloc1(k,&which);CHKERRQ(ierr);
  for (i=0;i<k;i++) which[i] = (i%2)? PETSC_TRUE: PETSC_FALSE;
  ierr = BVOrthogonalizeSomeColumn(X,k-1,which,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscFree(which);CHKERRQ(ierr);
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Check orthogonality */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Orthogonalization coefficients:\n");CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF,k-1,&z);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)z,"z");CHKERRQ(ierr);
  ierr = VecGetArray(z,&pz);CHKERRQ(ierr);
  ierr = BVDotColumn(X,k-1,pz);CHKERRQ(ierr);
  for (i=0;i<k-1;i++) {
    if (PetscAbsScalar(pz[i])<5.0*PETSC_MACHINE_EPSILON) pz[i]=0.0;
  }
  ierr = VecRestoreArray(z,&pz);CHKERRQ(ierr);
  ierr = VecView(z,view);CHKERRQ(ierr);
  ierr = VecDestroy(&z);CHKERRQ(ierr);

  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test8_1.out

   test:
      suffix: 1_cuda
      nsize: 1
      args: -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test8_1.out

   test:
      suffix: 2
      nsize: 1
      args: -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_refine never
      requires: !single
      output_file: output/test8_1.out

   test:
      suffix: 3
      nsize: 1
      args: -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_refine always
      output_file: output/test8_1.out

   test:
      suffix: 4
      nsize: 1
      args: -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_type mgs
      output_file: output/test8_1.out

   test:
      suffix: 4_cuda
      nsize: 1
      args: -bv_type svec -vec_type cuda -bv_orthog_type mgs
      requires: cuda
      output_file: output/test8_1.out

TEST*/
