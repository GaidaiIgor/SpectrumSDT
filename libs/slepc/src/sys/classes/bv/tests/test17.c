/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV bi-orthogonalization functions.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  BV             X,Y;
  Mat            M;
  Vec            v,t;
  PetscInt       i,j,n=20,k=7;
  PetscViewer    view;
  PetscBool      verbose;
  PetscReal      norm,condn=1.0;
  PetscScalar    alpha;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-condn",&condn,NULL);CHKERRQ(ierr);
  if (condn<1.0) SETERRQ(PETSC_COMM_WORLD,1,"The condition number must be > 1");
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test BV bi-orthogonalization with %D columns of length %D.\n",k,n);CHKERRQ(ierr);
  if (condn>1.0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD," - Using random BVs with condition number = %g\n",(double)condn);CHKERRQ(ierr);
  }

  /* Create template vector */
  ierr = VecCreate(PETSC_COMM_WORLD,&t);CHKERRQ(ierr);
  ierr = VecSetSizes(t,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(t);CHKERRQ(ierr);

  /* Create BV object X */
  ierr = BVCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)X,"X");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(X,t,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(X);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Fill X entries */
  if (condn==1.0) {
    for (j=0;j<k;j++) {
      ierr = BVGetColumn(X,j,&v);CHKERRQ(ierr);
      ierr = VecSet(v,0.0);CHKERRQ(ierr);
      for (i=0;i<=n/2;i++) {
        if (i+j<n) {
          alpha = (3.0*i+j-2)/(2*(i+j+1));
#if defined(PETSC_USE_COMPLEX)
          alpha += (1.2*i+j-2)/(0.1*(i+j+1))*PETSC_i;
#endif
          ierr = VecSetValue(v,i+j,alpha,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
      ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(X,j,&v);CHKERRQ(ierr);
    }
  } else {
    ierr = BVSetRandomCond(X,condn);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Create Y and fill its entries */
  ierr = BVDuplicate(X,&Y);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)Y,"Y");CHKERRQ(ierr);
  if (condn==1.0) {
    for (j=0;j<k;j++) {
      ierr = BVGetColumn(Y,j,&v);CHKERRQ(ierr);
      ierr = VecSet(v,0.0);CHKERRQ(ierr);
      for (i=PetscMin(n,2+(2*j)%6);i<PetscMin(n,6+(3*j)%9);i++) {
        if (i%5 && i!=j) {
          alpha = (1.5*i+j)/(2.2*(i-j));
          ierr = VecSetValue(v,i+j,alpha,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
      ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(Y,j,&v);CHKERRQ(ierr);
    }
  } else {
    ierr = BVSetRandomCond(Y,condn);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = BVView(Y,view);CHKERRQ(ierr);
  }

  /* Test BVBiorthonormalizeColumn */
  for (j=0;j<k;j++) {
    ierr = BVBiorthonormalizeColumn(X,Y,j,NULL);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
    ierr = BVView(Y,view);CHKERRQ(ierr);
  }

  /* Check orthogonality */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)M,"M");CHKERRQ(ierr);
  ierr = BVDot(X,Y,M);CHKERRQ(ierr);
  if (verbose) {
    ierr = MatView(M,view);CHKERRQ(ierr);
  }
  ierr = MatShift(M,-1.0);CHKERRQ(ierr);
  ierr = MatNorm(M,NORM_1,&norm);CHKERRQ(ierr);
  if (norm<200*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of bi-orthogonality < 200*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of bi-orthogonality: %g\n",(double)norm);CHKERRQ(ierr);
  }

  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = BVDestroy(&Y);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -bv_type {{vecs contiguous svec mat}} -bv_orthog_type cgs
      requires: double
      output_file: output/test17_1.out

   test:
      suffix: 1_cuda
      nsize: 1
      args: -bv_type svec -vec_type cuda -bv_orthog_type cgs
      requires: cuda
      output_file: output/test17_1.out

   test:
      suffix: 2
      nsize: 1
      args: -bv_type {{vecs contiguous svec mat}} -bv_orthog_type mgs
      output_file: output/test17_1.out

   test:
      suffix: 2_cuda
      nsize: 1
      args: -bv_type svec -vec_type cuda -bv_orthog_type mgs
      requires: cuda
      output_file: output/test17_1.out

TEST*/
