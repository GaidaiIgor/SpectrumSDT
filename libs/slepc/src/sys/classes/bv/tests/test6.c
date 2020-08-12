/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV orthogonalization functions with constraints.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  BV             X;
  Mat            M;
  Vec            v,t,*C;
  PetscInt       i,j,n=20,k=8,nc=2,nloc;
  PetscViewer    view;
  PetscBool      verbose;
  PetscReal      norm;
  PetscScalar    alpha;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-nc",&nc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test BV orthogonalization with %D columns + %D constraints, of length %D.\n",k,nc,n);CHKERRQ(ierr);

  /* Create template vector */
  ierr = VecCreate(PETSC_COMM_WORLD,&t);CHKERRQ(ierr);
  ierr = VecSetSizes(t,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(t);CHKERRQ(ierr);
  ierr = VecGetLocalSize(t,&nloc);CHKERRQ(ierr);

  /* Create BV object X */
  ierr = BVCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)X,"X");CHKERRQ(ierr);
  ierr = BVSetSizes(X,nloc,n,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(X);CHKERRQ(ierr);

  /* Generate constraints and attach them to X */
  if (nc>0) {
    ierr = VecDuplicateVecs(t,nc,&C);CHKERRQ(ierr);
    for (j=0;j<nc;j++) {
      for (i=0;i<=j;i++) {
        ierr = VecSetValue(C[j],i,1.0,INSERT_VALUES);CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(C[j]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(C[j]);CHKERRQ(ierr);
    }
    ierr = BVInsertConstraints(X,&nc,C);CHKERRQ(ierr);
    ierr = VecDestroyVecs(nc,&C);CHKERRQ(ierr);
  }

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

  /* Test BVOrthogonalizeColumn */
  for (j=0;j<k;j++) {
    ierr = BVOrthogonalizeColumn(X,j,NULL,&norm,NULL);CHKERRQ(ierr);
    alpha = 1.0/norm;
    ierr = BVScaleColumn(X,j,alpha);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Check orthogonality */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M);CHKERRQ(ierr);
  ierr = BVDot(X,X,M);CHKERRQ(ierr);
  ierr = MatShift(M,-1.0);CHKERRQ(ierr);
  ierr = MatNorm(M,NORM_1,&norm);CHKERRQ(ierr);
  if (norm<100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality < 100*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)norm);CHKERRQ(ierr);
  }

  ierr = MatDestroy(&M);CHKERRQ(ierr);
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
      output_file: output/test6_1.out

   test:
      suffix: 1_cuda
      nsize: 1
      args: -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test6_1.out

   test:
      suffix: 2
      nsize: 2
      args: -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test6_1.out

   test:
      suffix: 3
      nsize: 1
      args: -bv_type {{vecs contiguous svec mat}shared output} -bv_orthog_type mgs
      output_file: output/test6_1.out

   test:
      suffix: 3_cuda
      nsize: 1
      args: -bv_type svec -vec_type cuda -bv_orthog_type mgs
      requires: cuda
      output_file: output/test6_1.out

TEST*/
