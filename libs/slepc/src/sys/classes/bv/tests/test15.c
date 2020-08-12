/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BVGetSplit().\n\n";

#include <slepcbv.h>

/*
   Print the first row of a BV
 */
PetscErrorCode PrintFirstRow(BV X)
{
  PetscErrorCode    ierr;
  PetscMPIInt       rank;
  PetscInt          i,nloc,k,nc;
  const PetscScalar *pX;
  const char        *name;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)X),&rank);CHKERRQ(ierr);
  if (!rank) {
    ierr = BVGetActiveColumns(X,NULL,&k);CHKERRQ(ierr);
    ierr = BVGetSizes(X,&nloc,NULL,NULL);CHKERRQ(ierr);
    ierr = BVGetNumConstraints(X,&nc);CHKERRQ(ierr);
    ierr = PetscObjectGetName((PetscObject)X,&name);CHKERRQ(ierr);
    ierr = PetscPrintf(PetscObjectComm((PetscObject)X),"First row of %s =\n",name);CHKERRQ(ierr);
    ierr = BVGetArrayRead(X,&pX);CHKERRQ(ierr);
    for (i=0;i<nc+k;i++) {
      ierr = PetscPrintf(PetscObjectComm((PetscObject)X),"%g ",(double)PetscRealPart(pX[i*nloc]));CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PetscObjectComm((PetscObject)X),"\n");CHKERRQ(ierr);
    ierr = BVRestoreArrayRead(X,&pX);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  Vec            t,v,*C;
  BV             X,L,R;
  PetscInt       i,j,n=10,k=5,l=3,nc=0,nloc;
  PetscReal      norm;
  PetscScalar    alpha;
  PetscViewer    view;
  PetscBool      verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-nc",&nc,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test BVGetSplit (length %D, l=%D, k=%D, nc=%D).\n",n,l,k,nc);CHKERRQ(ierr);

  /* Create template vector */
  ierr = VecCreate(PETSC_COMM_WORLD,&t);CHKERRQ(ierr);
  ierr = VecSetSizes(t,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(t);CHKERRQ(ierr);
  ierr = VecGetLocalSize(t,&nloc);CHKERRQ(ierr);

  /* Create BV object X */
  ierr = BVCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)X,"X");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(X,t,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(X);CHKERRQ(ierr);

  /* Generate constraints and attach them to X */
  if (nc>0) {
    ierr = VecDuplicateVecs(t,nc,&C);CHKERRQ(ierr);
    for (j=0;j<nc;j++) {
      for (i=0;i<=j;i++) {
        ierr = VecSetValue(C[j],nc-i+1,1.0,INSERT_VALUES);CHKERRQ(ierr);
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
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Get split BVs */
  ierr = BVSetActiveColumns(X,l,k);CHKERRQ(ierr);
  ierr = BVGetSplit(X,&L,&R);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)L,"L");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)R,"R");CHKERRQ(ierr);

  if (verbose) {
    ierr = BVView(L,view);CHKERRQ(ierr);
    ierr = BVView(R,view);CHKERRQ(ierr);
  }

  /* Modify first column of R */
  ierr = BVGetColumn(R,0,&v);CHKERRQ(ierr);
  ierr = VecSet(v,-1.0);CHKERRQ(ierr);
  ierr = BVRestoreColumn(R,0,&v);CHKERRQ(ierr);

  /* Finished using the split BVs */
  ierr = BVRestoreSplit(X,&L,&R);CHKERRQ(ierr);
  ierr = PrintFirstRow(X);CHKERRQ(ierr);
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Get the left split BV only */
  ierr = BVGetSplit(X,&L,NULL);CHKERRQ(ierr);
  for (j=0;j<l;j++) {
    ierr = BVOrthogonalizeColumn(L,j,NULL,&norm,NULL);CHKERRQ(ierr);
    alpha = 1.0/norm;
    ierr = BVScaleColumn(L,j,alpha);CHKERRQ(ierr);
  }
  ierr = BVRestoreSplit(X,&L,NULL);CHKERRQ(ierr);
  ierr = PrintFirstRow(X);CHKERRQ(ierr);
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  /* Now get the right split BV after changing the number of leading columns */
  ierr = BVSetActiveColumns(X,l-1,k);CHKERRQ(ierr);
  ierr = BVGetSplit(X,NULL,&R);CHKERRQ(ierr);
  ierr = BVGetColumn(R,0,&v);CHKERRQ(ierr);
  ierr = BVInsertVec(X,0,v);CHKERRQ(ierr);
  ierr = BVRestoreColumn(R,0,&v);CHKERRQ(ierr);
  ierr = BVRestoreSplit(X,NULL,&R);CHKERRQ(ierr);
  ierr = PrintFirstRow(X);CHKERRQ(ierr);
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 2
      args: -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test15_1.out

   test:
      suffix: 1_cuda
      nsize: 2
      args: -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test15_1.out

   test:
      suffix: 2
      nsize: 2
      args: -nc 2 -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test15_2.out

   test:
      suffix: 2_cuda
      nsize: 2
      args: -nc 2 -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test15_2.out

TEST*/
