/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc singular value solver: "cyclic" (CUDA implementation)
*/
#include <slepc/private/svdimpl.h>
#include "../src/svd/impls/cyclic/cyclic.h"

PetscErrorCode MatMult_Cyclic_CUDA(Mat B,Vec x,Vec y)
{
  PetscErrorCode ierr;
  SVD            svd;
  SVD_CYCLIC     *cyclic;
  PetscScalar    *d_px,*d_py;
  PetscInt       m;

  PetscFunctionBegin;
  ierr = MatShellGetContext(B,(void**)&svd);CHKERRQ(ierr);
  cyclic = (SVD_CYCLIC*)svd->data;
  ierr = SVDMatGetLocalSize(svd,&m,NULL);CHKERRQ(ierr);
  ierr = VecCUDAGetArrayRead(x,(const PetscScalar**)&d_px);CHKERRQ(ierr);
  ierr = VecCUDAGetArray(y,&d_py);CHKERRQ(ierr);
  ierr = VecCUDAPlaceArray(cyclic->x1,d_px);CHKERRQ(ierr);
  ierr = VecCUDAPlaceArray(cyclic->x2,d_px+m);CHKERRQ(ierr);
  ierr = VecCUDAPlaceArray(cyclic->y1,d_py);CHKERRQ(ierr);
  ierr = VecCUDAPlaceArray(cyclic->y2,d_py+m);CHKERRQ(ierr);
  ierr = SVDMatMult(svd,PETSC_FALSE,cyclic->x2,cyclic->y1);CHKERRQ(ierr);
  ierr = SVDMatMult(svd,PETSC_TRUE,cyclic->x1,cyclic->y2);CHKERRQ(ierr);
  ierr = VecCUDAResetArray(cyclic->x1);CHKERRQ(ierr);
  ierr = VecCUDAResetArray(cyclic->x2);CHKERRQ(ierr);
  ierr = VecCUDAResetArray(cyclic->y1);CHKERRQ(ierr);
  ierr = VecCUDAResetArray(cyclic->y2);CHKERRQ(ierr);
  ierr = VecCUDARestoreArrayRead(x,(const PetscScalar**)&d_px);CHKERRQ(ierr);
  ierr = VecCUDARestoreArray(y,&d_py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatCreateVecs_Cyclic_CUDA(Mat B,Vec *right,Vec *left)
{
  PetscErrorCode ierr;
  SVD            svd;
  SVD_CYCLIC     *cyclic;
  PetscInt       M,N,m,n;
  PetscMPIInt    size;

  PetscFunctionBegin;
  ierr = MatShellGetContext(B,(void**)&svd);CHKERRQ(ierr);
  cyclic = (SVD_CYCLIC*)svd->data;
  ierr = SVDMatGetSize(svd,&M,&N);CHKERRQ(ierr);
  ierr = SVDMatGetLocalSize(svd,&m,&n);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)cyclic->mat),&size);CHKERRQ(ierr);
  if (right) {
    ierr = VecCreate(PetscObjectComm((PetscObject)cyclic->mat),right);CHKERRQ(ierr);
    ierr = VecSetSizes(*right,m+n,M+N);CHKERRQ(ierr);
    if (size>1) {
      ierr = VecSetType(*right,VECMPICUDA);CHKERRQ(ierr);
    } else {
      ierr = VecSetType(*right,VECSEQCUDA);CHKERRQ(ierr);
    }
  }
  if (left) {
    ierr = VecCreate(PetscObjectComm((PetscObject)cyclic->mat),left);CHKERRQ(ierr);
    ierr = VecSetSizes(*left,m+n,M+N);CHKERRQ(ierr);
    if (size>1) {
      ierr = VecSetType(*left,VECMPICUDA);CHKERRQ(ierr);
    } else {
      ierr = VecSetType(*left,VECSEQCUDA);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

