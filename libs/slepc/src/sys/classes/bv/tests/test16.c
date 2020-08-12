/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test tensor BV.\n\n";

#include <slepcbv.h>

int main(int argc,char **argv)
{
  PetscErrorCode    ierr;
  Vec               t,v;
  Mat               S,M,Q;
  BV                U,V,UU;
  PetscInt          i,ii,j,jj,n=10,k=6,l=3,d=3,deg,id,lds;
  PetscScalar       *pS,*q;
  PetscReal         norm;
  PetscViewer       view;
  PetscBool         verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-d",&d,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test tensor BV of degree %D with %D columns of dimension %D*d.\n",d,k,n);CHKERRQ(ierr);

  /* Create template vector */
  ierr = VecCreate(PETSC_COMM_WORLD,&t);CHKERRQ(ierr);
  ierr = VecSetSizes(t,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(t);CHKERRQ(ierr);

  /* Create BV object U */
  ierr = BVCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)U,"U");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(U,t,k+d-1);CHKERRQ(ierr);
  ierr = BVSetFromOptions(U);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)U,"U");CHKERRQ(ierr);

  /* Fill first d columns of U */
  for (j=0;j<d;j++) {
    ierr = BVGetColumn(U,j,&v);CHKERRQ(ierr);
    ierr = VecSet(v,0.0);CHKERRQ(ierr);
    for (i=0;i<4;i++) {
      if (i+j<n) {
        ierr = VecSetValue(v,i+j,(PetscScalar)(3*i+j-2),INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
    ierr = BVRestoreColumn(U,j,&v);CHKERRQ(ierr);
  }

  /* Create tensor BV */
  ierr = BVCreateTensor(U,d,&V);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)V,"V");CHKERRQ(ierr);
  ierr = BVTensorGetDegree(V,&deg);CHKERRQ(ierr);
  if (deg!=d) SETERRQ(PETSC_COMM_WORLD,1,"Wrong degree");

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
  ierr = BVView(V,view);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(view);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = BVView(V,view);CHKERRQ(ierr);
  }

  /* Build first column from previously introduced coefficients */
  ierr = BVTensorBuildFirstColumn(V,d);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"After building the first column - - - - -\n");CHKERRQ(ierr);
    ierr = BVView(V,view);CHKERRQ(ierr);
  }

  /* Test orthogonalization */
  ierr = BVTensorGetFactors(V,&UU,&S);CHKERRQ(ierr);
  ierr = BVGetActiveColumns(UU,NULL,&j);CHKERRQ(ierr);
  ierr = BVGetSizes(UU,NULL,NULL,&id);CHKERRQ(ierr);
  if (id!=k+d-1) SETERRQ(PETSC_COMM_WORLD,1,"Wrong dimensions");
  lds = id*d;
  for (jj=1;jj<k;jj++) {
    /* set new orthogonal column in U */
    ierr = BVGetColumn(UU,j,&v);CHKERRQ(ierr);
    ierr = VecSet(v,0.0);CHKERRQ(ierr);
    for (i=0;i<4;i++) {
      if (i+j<n) {
        ierr = VecSetValue(v,i+j,(PetscScalar)(3*i+j-2),INSERT_VALUES);CHKERRQ(ierr);
      }
    }
    ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
    ierr = BVRestoreColumn(UU,j,&v);CHKERRQ(ierr);
    ierr = BVOrthonormalizeColumn(UU,j,PETSC_TRUE,NULL,NULL);CHKERRQ(ierr);
    j++;
    ierr = BVSetActiveColumns(UU,0,j);CHKERRQ(ierr);
    /* set new column of S */
    ierr = MatDenseGetArray(S,&pS);CHKERRQ(ierr);
    for (ii=0;ii<d;ii++) {
      for (i=0;i<ii+jj+1;i++) {
        pS[i+ii*id+jj*lds] = (PetscScalar)(2*ii+i+0.5*jj);
      }
    }
    ierr = MatDenseRestoreArray(S,&pS);CHKERRQ(ierr);
    ierr = BVOrthonormalizeColumn(V,jj,PETSC_TRUE,NULL,NULL);CHKERRQ(ierr);
  }
  ierr = BVTensorRestoreFactors(V,&UU,&S);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"After orthogonalization - - - - -\n");CHKERRQ(ierr);
    ierr = BVView(V,view);CHKERRQ(ierr);
  }

  /* Check orthogonality */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M);CHKERRQ(ierr);
  ierr = BVDot(V,V,M);CHKERRQ(ierr);
  ierr = MatShift(M,-1.0);CHKERRQ(ierr);
  ierr = MatNorm(M,NORM_1,&norm);CHKERRQ(ierr);
  if (norm<100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality < 100*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)norm);CHKERRQ(ierr);
  }

  /* Test BVTensorCompress */
  ierr = BVSetActiveColumns(V,0,l);CHKERRQ(ierr);
  ierr = BVTensorCompress(V,0);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"After BVTensorCompress - - - - -\n");CHKERRQ(ierr);
    ierr = BVView(V,view);CHKERRQ(ierr);
  }

  /* Create Mat */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,l,NULL,&Q);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)Q,"Q");CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  for (i=0;i<k;i++)
    for (j=0;j<l;j++)
      q[i+j*k] = (i<j)? 2.0: -0.5;
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  if (verbose) {
    ierr = MatView(Q,NULL);CHKERRQ(ierr);
  }

  /* Test BVMultInPlace */
  ierr = BVMultInPlace(V,Q,1,l);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"After BVMultInPlace - - - - -\n");CHKERRQ(ierr);
    ierr = BVView(V,view);CHKERRQ(ierr);
  }

  /* Test BVNorm */
  ierr = BVNorm(V,NORM_1,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm: %g\n",(double)norm);CHKERRQ(ierr);

  ierr = BVDestroy(&U);CHKERRQ(ierr);
  ierr = BVDestroy(&V);CHKERRQ(ierr);
  ierr = MatDestroy(&Q);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 2
      args: -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test16_1.out
      filter: grep -v "doing matmult"

TEST*/
