/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV orthogonalization routines (CUDA)
*/

#include <slepc/private/bvimpl.h>          /*I   "slepcbv.h"   I*/
#include <slepcblaslapack.h>
#include <petsccublas.h>

/*
   BV_CleanCoefficients_CUDA - Sets to zero all entries of column j of the bv buffer
*/
PetscErrorCode BV_CleanCoefficients_CUDA(BV bv,PetscInt j,PetscScalar *h)
{
  PetscErrorCode ierr;
  PetscScalar    *d_hh,*d_a;
  PetscInt       i;
  cudaError_t    cerr;

  PetscFunctionBegin;
  if (!h) {
    ierr = VecCUDAGetArray(bv->buffer,&d_a);CHKERRQ(ierr);
    d_hh = d_a + j*(bv->nc+bv->m);
    cerr = cudaMemset(d_hh,0,(bv->nc+j)*sizeof(PetscScalar));CHKERRCUDA(cerr);
    cerr = WaitForGPU();CHKERRCUDA(cerr);
    ierr = VecCUDARestoreArray(bv->buffer,&d_a);CHKERRQ(ierr);
  } else { /* cpu memory */
    for (i=0;i<bv->nc+j;i++) h[i] = 0.0;
  }
  PetscFunctionReturn(0);
}

/*
   BV_AddCoefficients_CUDA - Add the contents of the scratch (0-th column) of the bv buffer
   into column j of the bv buffer
 */
PetscErrorCode BV_AddCoefficients_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscScalar *c)
{
  PetscErrorCode ierr;
  PetscScalar    *d_h,*d_c,sone=1.0;
  PetscInt       i;
  PetscBLASInt   idx,one=1;
  cublasStatus_t cberr;
  cublasHandle_t cublasv2handle;

  PetscFunctionBegin;
  if (!h) {
    ierr = PetscCUBLASGetHandle(&cublasv2handle);CHKERRQ(ierr);
    ierr = VecCUDAGetArray(bv->buffer,&d_c);CHKERRQ(ierr);
    d_h = d_c + j*(bv->nc+bv->m);
    ierr = PetscBLASIntCast(bv->nc+j,&idx);CHKERRQ(ierr);
    ierr = PetscLogGpuTimeBegin();CHKERRQ(ierr);
    cberr = cublasXaxpy(cublasv2handle,idx,&sone,d_c,one,d_h,one);CHKERRCUBLAS(cberr);
    ierr = PetscLogGpuTimeEnd();CHKERRQ(ierr);
    ierr = PetscLogGpuFlops(1.0*bv->nc+j);CHKERRQ(ierr);
    ierr = VecCUDARestoreArray(bv->buffer,&d_c);CHKERRQ(ierr);
  } else { /* cpu memory */
    for (i=0;i<bv->nc+j;i++) h[i] += c[i];
    ierr = PetscLogFlops(1.0*bv->nc+j);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   BV_SetValue_CUDA - Sets value in row j (counted after the constraints) of column k
   of the coefficients array
*/
PetscErrorCode BV_SetValue_CUDA(BV bv,PetscInt j,PetscInt k,PetscScalar *h,PetscScalar value)
{
  PetscErrorCode ierr;
  PetscScalar    *d_h,*a;
  cudaError_t    cerr;

  PetscFunctionBegin;
  if (!h) {
    ierr = VecCUDAGetArray(bv->buffer,&a);CHKERRQ(ierr);
    d_h = a + k*(bv->nc+bv->m) + bv->nc+j;
    cerr = cudaMemcpy(d_h,&value,sizeof(PetscScalar),cudaMemcpyHostToDevice);CHKERRCUDA(cerr);
    ierr = PetscLogCpuToGpu(sizeof(PetscScalar));CHKERRQ(ierr);
    cerr = WaitForGPU();CHKERRCUDA(cerr);
    ierr = VecCUDARestoreArray(bv->buffer,&a);CHKERRQ(ierr);
  } else { /* cpu memory */
    h[bv->nc+j] = value;
  }
  PetscFunctionReturn(0);
}

/*
   BV_SquareSum_CUDA - Returns the value h'*h, where h represents the contents of the
   coefficients array (up to position j)
*/
PetscErrorCode BV_SquareSum_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscReal *sum)
{
  PetscErrorCode    ierr;
  const PetscScalar *d_h;
  PetscScalar       dot;
  PetscInt          i;
  PetscBLASInt      idx,one=1;
  cublasStatus_t    cberr;
  cublasHandle_t    cublasv2handle;

  PetscFunctionBegin;
  if (!h) {
    ierr = PetscCUBLASGetHandle(&cublasv2handle);CHKERRQ(ierr);
    ierr = VecCUDAGetArrayRead(bv->buffer,&d_h);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(bv->nc+j,&idx);CHKERRQ(ierr);
    ierr = PetscLogGpuTimeBegin();CHKERRQ(ierr);
    cberr = cublasXdotc(cublasv2handle,idx,d_h,one,d_h,one,&dot);CHKERRCUBLAS(cberr);
    ierr = PetscLogGpuTimeEnd();CHKERRQ(ierr);
    ierr = PetscLogGpuFlops(2.0*bv->nc+j);CHKERRQ(ierr);
    *sum = PetscRealPart(dot);
    ierr = VecCUDARestoreArrayRead(bv->buffer,&d_h);CHKERRQ(ierr);
  } else { /* cpu memory */
    *sum = 0.0;
    for (i=0;i<bv->nc+j;i++) *sum += PetscRealPart(h[i]*PetscConj(h[i]));
    ierr = PetscLogFlops(2.0*bv->nc+j);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#define X_AXIS        0
#define BLOCK_SIZE_X 64
#define TILE_SIZE_X  16 /* work to be done by any thread on axis x */

/*
   Set the kernels grid dimensions
   xcount: number of kernel calls needed for the requested size
 */
PetscErrorCode SetGrid1D(PetscInt n, dim3 *dimGrid, dim3 *dimBlock,PetscInt *xcount)
{
  PetscInt              one=1;
  PetscBLASInt          card;
  struct cudaDeviceProp devprop;
  cudaError_t           cerr;

  PetscFunctionBegin;
  *xcount = 1;
  if (n>BLOCK_SIZE_X) {
    dimBlock->x = BLOCK_SIZE_X;
    dimGrid->x = (n+BLOCK_SIZE_X*TILE_SIZE_X-one)/BLOCK_SIZE_X*TILE_SIZE_X;
  } else {
    dimBlock->x = (n+TILE_SIZE_X-one)/TILE_SIZE_X;
    dimGrid->x = one;
  }
  cerr = cudaGetDevice(&card);CHKERRCUDA(cerr);
  cerr = cudaGetDeviceProperties(&devprop,card);CHKERRCUDA(cerr);
  if (dimGrid->x>(unsigned)devprop.maxGridSize[X_AXIS]) {
    *xcount = (dimGrid->x+devprop.maxGridSize[X_AXIS]-one)/devprop.maxGridSize[X_AXIS];
    dimGrid->x = devprop.maxGridSize[X_AXIS];
  }
  PetscFunctionReturn(0);
}

/* pointwise multiplication */
__global__ void PointwiseMult_kernel(PetscInt xcount,PetscScalar *a,const PetscScalar *b,PetscInt n)
{
  PetscInt i,x;

  x = xcount*gridDim.x*blockDim.x+blockIdx.x*blockDim.x*TILE_SIZE_X+threadIdx.x*TILE_SIZE_X;
  for (i=x;i<x+TILE_SIZE_X&&i<n;i++) {
    a[i] *= PetscRealPart(b[i]);
  }
}

/* pointwise division */
__global__ void PointwiseDiv_kernel(PetscInt xcount,PetscScalar *a,const PetscScalar *b,PetscInt n)
{
  PetscInt i,x;

  x = xcount*gridDim.x*blockDim.x+blockIdx.x*blockDim.x*TILE_SIZE_X+threadIdx.x*TILE_SIZE_X;
  for (i=x;i<x+TILE_SIZE_X&&i<n;i++) {
    a[i] /= PetscRealPart(b[i]);
  }
}

/*
   BV_ApplySignature_CUDA - Computes the pointwise product h*omega, where h represents
   the contents of the coefficients array (up to position j) and omega is the signature;
   if inverse=TRUE then the operation is h/omega
*/
PetscErrorCode BV_ApplySignature_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscBool inverse)
{
  PetscErrorCode    ierr;
  PetscScalar       *d_h;
  const PetscScalar *d_omega,*omega;
  PetscInt          i,xcount;
  dim3              blocks3d, threads3d;
  cudaError_t       cerr;

  PetscFunctionBegin;
  if (!(bv->nc+j)) PetscFunctionReturn(0);
  if (!h) {
    ierr = VecCUDAGetArray(bv->buffer,&d_h);CHKERRQ(ierr);
    ierr = VecCUDAGetArrayRead(bv->omega,&d_omega);CHKERRQ(ierr);
    ierr = SetGrid1D(bv->nc+j,&blocks3d,&threads3d,&xcount);CHKERRQ(ierr);
    ierr = PetscLogGpuTimeBegin();CHKERRQ(ierr);
    if (inverse) {
      for (i=0;i<xcount;i++) {
        PointwiseDiv_kernel<<<blocks3d,threads3d>>>(i,d_h,d_omega,bv->nc+j);
      }
    } else {
      for (i=0;i<xcount;i++) {
        PointwiseMult_kernel<<<blocks3d,threads3d>>>(i,d_h,d_omega,bv->nc+j);
      }
    }
    cerr = cudaGetLastError();CHKERRCUDA(cerr);
    ierr = PetscLogGpuTimeEnd();CHKERRQ(ierr);
    ierr = PetscLogGpuFlops(1.0*bv->nc+j);CHKERRQ(ierr);
    cerr = WaitForGPU();CHKERRCUDA(cerr);
    ierr = VecCUDARestoreArrayRead(bv->omega,&d_omega);CHKERRQ(ierr);
    ierr = VecCUDARestoreArray(bv->buffer,&d_h);CHKERRQ(ierr);
  } else {
    ierr = VecGetArrayRead(bv->omega,&omega);CHKERRQ(ierr);
    if (inverse) for (i=0;i<bv->nc+j;i++) h[i] /= PetscRealPart(omega[i]);
    else for (i=0;i<bv->nc+j;i++) h[i] *= PetscRealPart(omega[i]);
    ierr = VecRestoreArrayRead(bv->omega,&omega);CHKERRQ(ierr);
    ierr = PetscLogFlops(1.0*bv->nc+j);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   BV_SquareRoot_CUDA - Returns the square root of position j (counted after the constraints)
   of the coefficients array
*/
PetscErrorCode BV_SquareRoot_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscReal *beta)
{
  PetscErrorCode    ierr;
  const PetscScalar *d_h;
  PetscScalar       hh;
  cudaError_t       cerr;

  PetscFunctionBegin;
  if (!h) {
    ierr = VecCUDAGetArrayRead(bv->buffer,&d_h);CHKERRQ(ierr);
    cerr = cudaMemcpy(&hh,d_h+bv->nc+j,sizeof(PetscScalar),cudaMemcpyDeviceToHost);CHKERRCUDA(cerr);
    ierr = PetscLogGpuToCpu(sizeof(PetscScalar));CHKERRQ(ierr);
    cerr = WaitForGPU();CHKERRCUDA(cerr);
    ierr = BV_SafeSqrt(bv,hh,beta);CHKERRQ(ierr);
    ierr = VecCUDARestoreArrayRead(bv->buffer,&d_h);CHKERRQ(ierr);
  } else {
    ierr = BV_SafeSqrt(bv,h[bv->nc+j],beta);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   BV_StoreCoefficients_CUDA - Copy the contents of the coefficients array to an array dest
   provided by the caller (only values from l to j are copied)
*/
PetscErrorCode BV_StoreCoefficients_CUDA(BV bv,PetscInt j,PetscScalar *h,PetscScalar *dest)
{
  PetscErrorCode    ierr;
  const PetscScalar *d_h,*d_a;
  PetscInt          i;
  cudaError_t       cerr;

  PetscFunctionBegin;
  if (!h) {
    ierr = VecCUDAGetArrayRead(bv->buffer,&d_a);CHKERRQ(ierr);
    d_h = d_a + j*(bv->nc+bv->m)+bv->nc;
    cerr = cudaMemcpy(dest-bv->l,d_h,(j-bv->l)*sizeof(PetscScalar),cudaMemcpyDeviceToHost);CHKERRCUDA(cerr);
    ierr = PetscLogGpuToCpu((j-bv->l)*sizeof(PetscScalar));CHKERRQ(ierr);
    cerr = WaitForGPU();CHKERRCUDA(cerr);
    ierr = VecCUDARestoreArrayRead(bv->buffer,&d_a);CHKERRQ(ierr);
  } else {
    for (i=bv->l;i<j;i++) dest[i-bv->l] = h[bv->nc+i];
  }
  PetscFunctionReturn(0);
}

