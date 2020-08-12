/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV private kernels that use the BLAS
*/

#include <slepc/private/bvimpl.h>
#include <slepcblaslapack.h>

#define BLOCKSIZE 64

/*
    C := alpha*A*B + beta*C

    A is mxk (ld=m), B is kxn (ld=ldb), C is mxn (ld=m)
*/
PetscErrorCode BVMult_BLAS_Private(BV bv,PetscInt m_,PetscInt n_,PetscInt k_,PetscInt ldb_,PetscScalar alpha,const PetscScalar *A,const PetscScalar *B,PetscScalar beta,PetscScalar *C)
{
  PetscErrorCode ierr;
  PetscBLASInt   m,n,k,ldb;
#if defined(PETSC_HAVE_FBLASLAPACK) || defined(PETSC_HAVE_F2CBLASLAPACK)
  PetscBLASInt   l,bs=BLOCKSIZE;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(m_,&m);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k_,&k);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldb_,&ldb);CHKERRQ(ierr);
#if defined(PETSC_HAVE_FBLASLAPACK) || defined(PETSC_HAVE_F2CBLASLAPACK)
  l = m % bs;
  if (l) PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&l,&n,&k,&alpha,(PetscScalar*)A,&m,(PetscScalar*)B,&ldb,&beta,C,&m));
  for (;l<m;l+=bs) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&bs,&n,&k,&alpha,(PetscScalar*)A+l,&m,(PetscScalar*)B,&ldb,&beta,C+l,&m));
  }
#else
  if (m) PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&m,&n,&k,&alpha,(PetscScalar*)A,&m,(PetscScalar*)B,&ldb,&beta,C,&m));
#endif
  ierr = PetscLogFlops(2.0*m*n*k);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    y := alpha*A*x + beta*y

    A is nxk (ld=n)
*/
PetscErrorCode BVMultVec_BLAS_Private(BV bv,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *A,const PetscScalar *x,PetscScalar beta,PetscScalar *y)
{
  PetscErrorCode ierr;
  PetscBLASInt   n,k,one=1;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k_,&k);CHKERRQ(ierr);
  if (n) PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n,&k,&alpha,A,&n,x,&one,&beta,y,&one));
  ierr = PetscLogFlops(2.0*n*k);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    A(:,s:e-1) := A*B(:,s:e-1)

    A is mxk (ld=m), B is kxn (ld=ldb)  n=e-s
*/
PetscErrorCode BVMultInPlace_BLAS_Private(BV bv,PetscInt m_,PetscInt k_,PetscInt ldb_,PetscInt s,PetscInt e,PetscScalar *A,const PetscScalar *B,PetscBool btrans)
{
  PetscErrorCode ierr;
  PetscScalar    *pb,zero=0.0,one=1.0;
  PetscBLASInt   m,n,k,l,ldb,bs=BLOCKSIZE;
  PetscInt       j,n_=e-s;
  const char     *bt;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(m_,&m);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k_,&k);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldb_,&ldb);CHKERRQ(ierr);
  ierr = BVAllocateWork_Private(bv,BLOCKSIZE*n_);CHKERRQ(ierr);
  if (btrans) {
    pb = (PetscScalar*)B+s;
    bt = "C";
  } else {
    pb = (PetscScalar*)B+s*ldb;
    bt = "N";
  }
  l = m % bs;
  if (l) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N",bt,&l,&n,&k,&one,A,&m,pb,&ldb,&zero,bv->work,&l));
    for (j=0;j<n;j++) {
      ierr = PetscArraycpy(A+(s+j)*m,bv->work+j*l,l);CHKERRQ(ierr);
    }
  }
  for (;l<m;l+=bs) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N",bt,&bs,&n,&k,&one,A+l,&m,pb,&ldb,&zero,bv->work,&bs));
    for (j=0;j<n;j++) {
      ierr = PetscArraycpy(A+(s+j)*m+l,bv->work+j*bs,bs);CHKERRQ(ierr);
    }
  }
  ierr = PetscLogFlops(2.0*m*n*k);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    V := V*B

    V is mxn (ld=m), B is nxn (ld=k)
*/
PetscErrorCode BVMultInPlace_Vecs_Private(BV bv,PetscInt m_,PetscInt n_,PetscInt k_,Vec *V,const PetscScalar *B,PetscBool btrans)
{
  PetscErrorCode    ierr;
  PetscScalar       zero=0.0,one=1.0,*out,*pout;
  const PetscScalar *pin;
  PetscBLASInt      m,n,k,l,bs=BLOCKSIZE;
  PetscInt          j;
  const char        *bt;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(m_,&m);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k_,&k);CHKERRQ(ierr);
  ierr = BVAllocateWork_Private(bv,2*BLOCKSIZE*n_);CHKERRQ(ierr);
  out = bv->work+BLOCKSIZE*n_;
  if (btrans) bt = "C";
  else bt = "N";
  l = m % bs;
  if (l) {
    for (j=0;j<n;j++) {
      ierr = VecGetArrayRead(V[j],&pin);CHKERRQ(ierr);
      ierr = PetscArraycpy(bv->work+j*l,pin,l);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(V[j],&pin);CHKERRQ(ierr);
    }
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N",bt,&l,&n,&n,&one,bv->work,&l,(PetscScalar*)B,&k,&zero,out,&l));
    for (j=0;j<n;j++) {
      ierr = VecGetArray(V[j],&pout);CHKERRQ(ierr);
      ierr = PetscArraycpy(pout,out+j*l,l);CHKERRQ(ierr);
      ierr = VecRestoreArray(V[j],&pout);CHKERRQ(ierr);
    }
  }
  for (;l<m;l+=bs) {
    for (j=0;j<n;j++) {
      ierr = VecGetArrayRead(V[j],&pin);CHKERRQ(ierr);
      ierr = PetscArraycpy(bv->work+j*bs,pin+l,bs);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(V[j],&pin);CHKERRQ(ierr);
    }
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N",bt,&bs,&n,&n,&one,bv->work,&bs,(PetscScalar*)B,&k,&zero,out,&bs));
    for (j=0;j<n;j++) {
      ierr = VecGetArray(V[j],&pout);CHKERRQ(ierr);
      ierr = PetscArraycpy(pout+l,out+j*bs,bs);CHKERRQ(ierr);
      ierr = VecRestoreArray(V[j],&pout);CHKERRQ(ierr);
    }
  }
  ierr = PetscLogFlops(2.0*n*n*k);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    B := alpha*A + beta*B

    A,B are nxk (ld=n)
*/
PetscErrorCode BVAXPY_BLAS_Private(BV bv,PetscInt n_,PetscInt k_,PetscScalar alpha,const PetscScalar *A,PetscScalar beta,PetscScalar *B)
{
  PetscErrorCode ierr;
  PetscBLASInt   m,one=1;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(n_*k_,&m);CHKERRQ(ierr);
  if (beta!=(PetscScalar)1.0) {
    PetscStackCallBLAS("BLASscal",BLASscal_(&m,&beta,B,&one));
    ierr = PetscLogFlops(m);CHKERRQ(ierr);
  }
  PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&m,&alpha,A,&one,B,&one));
  ierr = PetscLogFlops(2.0*m);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    C := A'*B

    A' is mxk (ld=k), B is kxn (ld=k), C is mxn (ld=ldc)
*/
PetscErrorCode BVDot_BLAS_Private(BV bv,PetscInt m_,PetscInt n_,PetscInt k_,PetscInt ldc_,const PetscScalar *A,const PetscScalar *B,PetscScalar *C,PetscBool mpi)
{
  PetscErrorCode ierr;
  PetscScalar    zero=0.0,one=1.0,*CC;
  PetscBLASInt   m,n,k,ldc,j;
  PetscMPIInt    len;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(m_,&m);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k_,&k);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldc_,&ldc);CHKERRQ(ierr);
  if (mpi) {
    if (ldc==m) {
      ierr = BVAllocateWork_Private(bv,m*n);CHKERRQ(ierr);
      if (k) PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&m,&n,&k,&one,(PetscScalar*)A,&k,(PetscScalar*)B,&k,&zero,bv->work,&ldc));
      else { ierr = PetscArrayzero(bv->work,m*n);CHKERRQ(ierr); }
      ierr = PetscMPIIntCast(m*n,&len);CHKERRQ(ierr);
      ierr = MPI_Allreduce(bv->work,C,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv));CHKERRQ(ierr);
    } else {
      ierr = BVAllocateWork_Private(bv,2*m*n);CHKERRQ(ierr);
      CC = bv->work+m*n;
      if (k) PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&m,&n,&k,&one,(PetscScalar*)A,&k,(PetscScalar*)B,&k,&zero,bv->work,&m));
      else { ierr = PetscArrayzero(bv->work,m*n);CHKERRQ(ierr); }
      ierr = PetscMPIIntCast(m*n,&len);CHKERRQ(ierr);
      ierr = MPI_Allreduce(bv->work,CC,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv));CHKERRQ(ierr);
      for (j=0;j<n;j++) {
        ierr = PetscArraycpy(C+j*ldc,CC+j*m,m);CHKERRQ(ierr);
      }
    }
  } else {
    if (k) PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&m,&n,&k,&one,(PetscScalar*)A,&k,(PetscScalar*)B,&k,&zero,C,&ldc));
  }
  ierr = PetscLogFlops(2.0*m*n*k);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    y := A'*x

    A is nxk (ld=n)
*/
PetscErrorCode BVDotVec_BLAS_Private(BV bv,PetscInt n_,PetscInt k_,const PetscScalar *A,const PetscScalar *x,PetscScalar *y,PetscBool mpi)
{
  PetscErrorCode ierr;
  PetscScalar    zero=0.0,done=1.0;
  PetscBLASInt   n,k,one=1;
  PetscMPIInt    len;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k_,&k);CHKERRQ(ierr);
  if (mpi) {
    ierr = BVAllocateWork_Private(bv,k);CHKERRQ(ierr);
    if (n) {
      PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&n,&k,&done,A,&n,x,&one,&zero,bv->work,&one));
    } else {
      ierr = PetscArrayzero(bv->work,k);CHKERRQ(ierr);
    }
    ierr = PetscMPIIntCast(k,&len);CHKERRQ(ierr);
    ierr = MPI_Allreduce(bv->work,y,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)bv));CHKERRQ(ierr);
  } else {
    if (n) PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&n,&k,&done,A,&n,x,&one,&zero,y,&one));
  }
  ierr = PetscLogFlops(2.0*n*k);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    Scale n scalars
*/
PetscErrorCode BVScale_BLAS_Private(BV bv,PetscInt n_,PetscScalar *A,PetscScalar alpha)
{
  PetscErrorCode ierr;
  PetscBLASInt   n,one=1;

  PetscFunctionBegin;
  if (alpha == (PetscScalar)0.0) {
    ierr = PetscArrayzero(A,n_);CHKERRQ(ierr);
  } else if (alpha!=(PetscScalar)1.0) {
    ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASscal",BLASscal_(&n,&alpha,A,&one));
    ierr = PetscLogFlops(n);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

