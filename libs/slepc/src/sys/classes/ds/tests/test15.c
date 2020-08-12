/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSPEP.\n\n";

#include <slepcds.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DS             ds;
  SlepcSC        sc;
  Mat            X;
  Vec            x0;
  PetscScalar    *K,*C,*M,*wr,*wi,z=1.0;
  PetscReal      re,im,nrm,*pbc;
  PetscInt       i,j,n=10,d=2,ld;
  PetscViewer    viewer;
  PetscBool      verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Solve a Dense System of type PEP - n=%D.\n",n);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);

  /* Create DS object */
  ierr = DSCreate(PETSC_COMM_WORLD,&ds);CHKERRQ(ierr);
  ierr = DSSetType(ds,DSPEP);CHKERRQ(ierr);
  ierr = DSSetFromOptions(ds);CHKERRQ(ierr);
  ierr = DSPEPSetDegree(ds,d);CHKERRQ(ierr);

  /* Set dimensions */
  ld = n+2;  /* test leading dimension larger than n */
  ierr = DSAllocate(ds,ld);CHKERRQ(ierr);
  ierr = DSSetDimensions(ds,n,0,0,0);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
  ierr = DSView(ds,viewer);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Fill matrices */
  ierr = DSGetArray(ds,DS_MAT_E0,&K);CHKERRQ(ierr);
  for (i=0;i<n-1;i++) K[i+i*ld] = 2.0*n;
  K[n-1+(n-1)*ld] = 1.0*n;
  for (i=1;i<n;i++) {
    K[i+(i-1)*ld] = -1.0*n;
    K[(i-1)+i*ld] = -1.0*n;
  }
  ierr = DSRestoreArray(ds,DS_MAT_E0,&K);CHKERRQ(ierr);
  ierr = DSGetArray(ds,DS_MAT_E1,&C);CHKERRQ(ierr);
  C[n-1+(n-1)*ld] = 2.0*PETSC_PI/z;
  ierr = DSRestoreArray(ds,DS_MAT_E1,&C);CHKERRQ(ierr);
  ierr = DSGetArray(ds,DS_MAT_E2,&M);CHKERRQ(ierr);
  for (i=0;i<n-1;i++) M[i+i*ld] = -4.0*PETSC_PI*PETSC_PI/n;
  M[i+i*ld] = -2.0*PETSC_PI*PETSC_PI/n;
  ierr = DSRestoreArray(ds,DS_MAT_E2,&M);CHKERRQ(ierr);

  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Initial - - - - - - - - -\n");CHKERRQ(ierr);
    ierr = DSView(ds,viewer);CHKERRQ(ierr);
  }

  /* Solve */
  ierr = PetscMalloc2(d*n,&wr,d*n,&wi);CHKERRQ(ierr);
  ierr = DSGetSlepcSC(ds,&sc);CHKERRQ(ierr);
  sc->comparison    = SlepcCompareLargestReal;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  ierr = DSSolve(ds,wr,wi);CHKERRQ(ierr);
  ierr = DSSort(ds,wr,wi,NULL,NULL,NULL);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"After solve - - - - - - - - -\n");CHKERRQ(ierr);
    ierr = DSView(ds,viewer);CHKERRQ(ierr);
  }

  /* Print polynomial coefficients */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Polynomial coefficients (alpha,beta,gamma) =\n");CHKERRQ(ierr);
  ierr = DSPEPGetCoefficients(ds,&pbc);CHKERRQ(ierr);
  for (j=0;j<3;j++) {
    for (i=0;i<d+1;i++) {
      ierr = PetscViewerASCIIPrintf(viewer,"  %.5f",(double)pbc[j+3*i]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  ierr = PetscFree(pbc);CHKERRQ(ierr);

  /* Print eigenvalues */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed eigenvalues =\n");CHKERRQ(ierr);
  for (i=0;i<d*n;i++) {
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(wr[i]);
    im = PetscImaginaryPart(wr[i]);
#else
    re = wr[i];
    im = wi[i];
#endif
    if (PetscAbs(im)<1e-10) {
      ierr = PetscViewerASCIIPrintf(viewer,"  %.5f\n",(double)re);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"  %.5f%+.5fi\n",(double)re,(double)im);CHKERRQ(ierr);
    }
  }

  /* Eigenvectors */
  ierr = DSVectors(ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);  /* all eigenvectors */
  ierr = DSGetMat(ds,DS_MAT_X,&X);CHKERRQ(ierr);
  ierr = MatCreateVecs(X,NULL,&x0);CHKERRQ(ierr);
  ierr = MatGetColumnVector(X,x0,0);CHKERRQ(ierr);
  ierr = VecNorm(x0,NORM_2,&nrm);CHKERRQ(ierr);
  ierr = MatDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&x0);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of 1st vector = %.3f\n",(double)nrm);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"After vectors - - - - - - - - -\n");CHKERRQ(ierr);
    ierr = DSView(ds,viewer);CHKERRQ(ierr);
  }

  ierr = PetscFree2(wr,wi);CHKERRQ(ierr);
  ierr = DSDestroy(&ds);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      requires: !single

TEST*/
