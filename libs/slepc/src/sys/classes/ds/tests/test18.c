/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test DSSynchronize() on a NHEP.\n\n";

#include <slepcds.h>

PetscErrorCode CheckArray(PetscScalar *A,const char *label,PetscInt k)
{
  PetscErrorCode ierr;
  PetscInt       j;
  PetscMPIInt    p,size,rank;
  PetscScalar    dif,*buf;
  PetscReal      error;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  if (rank) {
    ierr = MPI_Send(A,k,MPIU_SCALAR,0,111,PETSC_COMM_WORLD);CHKERRQ(ierr);
  } else {
    ierr = PetscMalloc1(k,&buf);CHKERRQ(ierr);
    for (p=1;p<size;p++) {
      ierr = MPI_Recv(buf,k,MPIU_SCALAR,p,111,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);CHKERRQ(ierr);
      dif = 0.0;
      for (j=0;j<k;j++) dif += A[j]-buf[j];
      error = PetscAbsScalar(dif);
      if (error>10*PETSC_MACHINE_EPSILON) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," Array %s differs in proc %d: %g\n",label,(int)p,(double)error);CHKERRQ(ierr);
      }
    }
    ierr = PetscFree(buf);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DS             ds;
  SlepcSC        sc;
  PetscScalar    *A,*Q,*wr,*wi;
  PetscReal      re,im;
  PetscInt       i,j,n=10;
  PetscMPIInt    size;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Solve a Dense System of type NHEP - dimension %D.\n",n);CHKERRQ(ierr);

  /* Create DS object */
  ierr = DSCreate(PETSC_COMM_WORLD,&ds);CHKERRQ(ierr);
  ierr = DSSetType(ds,DSNHEP);CHKERRQ(ierr);
  ierr = DSSetFromOptions(ds);CHKERRQ(ierr);
  ierr = DSAllocate(ds,n);CHKERRQ(ierr);
  ierr = DSSetDimensions(ds,n,0,0,0);CHKERRQ(ierr);

  /* Fill with Grcar matrix */
  ierr = DSGetArray(ds,DS_MAT_A,&A);CHKERRQ(ierr);
  for (i=1;i<n;i++) A[i+(i-1)*n]=-1.0;
  for (j=0;j<4;j++) {
    for (i=0;i<n-j;i++) A[i+(i+j)*n]=1.0;
  }
  ierr = DSRestoreArray(ds,DS_MAT_A,&A);CHKERRQ(ierr);
  ierr = DSSetState(ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);

  /* Solve */
  ierr = PetscMalloc2(n,&wr,n,&wi);CHKERRQ(ierr);
  ierr = DSGetSlepcSC(ds,&sc);CHKERRQ(ierr);
  sc->comparison    = SlepcCompareLargestMagnitude;
  sc->comparisonctx = NULL;
  sc->map           = NULL;
  sc->mapobj        = NULL;
  ierr = DSSolve(ds,wr,wi);CHKERRQ(ierr);
  ierr = DSSort(ds,wr,wi,NULL,NULL,NULL);CHKERRQ(ierr);

  /* Print eigenvalues */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed eigenvalues =\n");CHKERRQ(ierr);
  for (i=0;i<n;i++) {
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(wr[i]);
    im = PetscImaginaryPart(wr[i]);
#else
    re = wr[i];
    im = wi[i];
#endif
    if (PetscAbs(im)<1e-10) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  %.5f\n",(double)re);CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  %.5f%+.5fi\n",(double)re,(double)im);CHKERRQ(ierr);
    }
  }

  /* Synchronize data and check */
  ierr = DSSynchronize(ds,wr,wi);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size>1) {
    ierr = CheckArray(wr,"wr",n);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    ierr = CheckArray(wi,"wi",n);CHKERRQ(ierr);
#endif
    ierr = DSGetArray(ds,DS_MAT_A,&A);CHKERRQ(ierr);
    ierr = CheckArray(A,"A",n*n);CHKERRQ(ierr);
    ierr = DSRestoreArray(ds,DS_MAT_A,&A);CHKERRQ(ierr);
    ierr = DSGetArray(ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
    ierr = CheckArray(Q,"Q",n*n);CHKERRQ(ierr);
    ierr = DSRestoreArray(ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
  }

  ierr = PetscFree2(wr,wi);CHKERRQ(ierr);
  ierr = DSDestroy(&ds);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: {{1 2 3}}
      args: -ds_parallel redundant
      filter: sed -e "s/[+-]\([0-9]\.[0-9]*i\)/+-\\1/"

   test:
      suffix: 2
      nsize: {{1 2 3}}
      args: -ds_parallel synchronized
      filter: sed -e "s/[+-]\([0-9]\.[0-9]*i\)/+-\\1/"
      output_file: output/test18_1.out

TEST*/
