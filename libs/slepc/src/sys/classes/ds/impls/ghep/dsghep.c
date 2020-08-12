/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/dsimpl.h>
#include <slepcblaslapack.h>

PetscErrorCode DSAllocate_GHEP(DS ds,PetscInt ld)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DSAllocateMat_Private(ds,DS_MAT_A);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_B);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_Q);CHKERRQ(ierr);
  ierr = PetscFree(ds->perm);CHKERRQ(ierr);
  ierr = PetscMalloc1(ld,&ds->perm);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)ds,ld*sizeof(PetscInt));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DSView_GHEP(DS ds,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  PetscViewerFormat format;

  PetscFunctionBegin;
  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
  ierr = DSViewMat(ds,viewer,DS_MAT_A);CHKERRQ(ierr);
  ierr = DSViewMat(ds,viewer,DS_MAT_B);CHKERRQ(ierr);
  if (ds->state>DS_STATE_INTERMEDIATE) { ierr = DSViewMat(ds,viewer,DS_MAT_Q);CHKERRQ(ierr); }
  if (ds->mat[DS_MAT_X]) { ierr = DSViewMat(ds,viewer,DS_MAT_X);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode DSVectors_GHEP(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscScalar    *Q = ds->mat[DS_MAT_Q];
  PetscInt       ld = ds->ld,i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (rnorm) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
  switch (mat) {
    case DS_MAT_X:
    case DS_MAT_Y:
      if (j) {
        if (ds->state>=DS_STATE_CONDENSED) {
          ierr = PetscArraycpy(ds->mat[mat]+(*j)*ld,Q+(*j)*ld,ld);CHKERRQ(ierr);
        } else {
          ierr = PetscArrayzero(ds->mat[mat]+(*j)*ld,ld);CHKERRQ(ierr);
          *(ds->mat[mat]+(*j)+(*j)*ld) = 1.0;
        }
      } else {
        if (ds->state>=DS_STATE_CONDENSED) {
          ierr = PetscArraycpy(ds->mat[mat],Q,ld*ld);CHKERRQ(ierr);
        } else {
          ierr = PetscArrayzero(ds->mat[mat],ld*ld);CHKERRQ(ierr);
          for (i=0;i<ds->n;i++) *(ds->mat[mat]+i+i*ld) = 1.0;
        }
      }
      break;
    case DS_MAT_U:
    case DS_MAT_VT:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSort_GHEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscErrorCode ierr;
  PetscInt       n,l,i,*perm,ld=ds->ld;
  PetscScalar    *A;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(0);
  n = ds->n;
  l = ds->l;
  A  = ds->mat[DS_MAT_A];
  perm = ds->perm;
  for (i=l;i<n;i++) wr[i] = A[i+i*ld];
  if (rr) {
    ierr = DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_FALSE);CHKERRQ(ierr);
  } else {
    ierr = DSSortEigenvalues_Private(ds,wr,NULL,perm,PETSC_FALSE);CHKERRQ(ierr);
  }
  for (i=l;i<n;i++) A[i+i*ld] = wr[perm[i]];
  for (i=l;i<n;i++) wr[i] = A[i+i*ld];
  ierr = DSPermuteColumns_Private(ds,l,n,DS_MAT_Q,perm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DSSolve_GHEP(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscScalar    *work,*A,*B,*Q;
  PetscBLASInt   itype = 1,*iwork,info,n1,liwork,ld,lrwork=0,lwork;
  PetscInt       off,i;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork,*rr;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n-ds->l,&n1);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(5*ds->n+3,&liwork);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscBLASIntCast(ds->n*ds->n+2*ds->n,&lwork);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(2*ds->n*ds->n+5*ds->n+1+n1,&lrwork);CHKERRQ(ierr);
#else
  ierr = PetscBLASIntCast(2*ds->n*ds->n+6*ds->n+1,&lwork);CHKERRQ(ierr);
#endif
  ierr = DSAllocateWork_Private(ds,lwork,lrwork,liwork);CHKERRQ(ierr);
  work = ds->work;
  iwork = ds->iwork;
  off = ds->l+ds->l*ld;
  A = ds->mat[DS_MAT_A];
  B = ds->mat[DS_MAT_B];
  Q = ds->mat[DS_MAT_Q];
#if defined(PETSC_USE_COMPLEX)
  rr = ds->rwork;
  rwork = ds->rwork + n1;
  lrwork = ds->lrwork - n1;
  PetscStackCallBLAS("LAPACKsygvd",LAPACKsygvd_(&itype,"V","U",&n1,A+off,&ld,B+off,&ld,rr,work,&lwork,rwork,&lrwork,iwork,&liwork,&info));
  for (i=0;i<n1;i++) wr[ds->l+i] = rr[i];
#else
  PetscStackCallBLAS("LAPACKsygvd",LAPACKsygvd_(&itype,"V","U",&n1,A+off,&ld,B+off,&ld,wr+ds->l,work,&lwork,iwork,&liwork,&info));
#endif
  SlepcCheckLapackInfo("sygvd",info);
  ierr = PetscArrayzero(Q+ds->l*ld,n1*ld);CHKERRQ(ierr);
  for (i=ds->l;i<ds->n;i++) {
    ierr = PetscArraycpy(Q+ds->l+i*ld,A+ds->l+i*ld,n1);CHKERRQ(ierr);
  }
  ierr = PetscArrayzero(B+ds->l*ld,n1*ld);CHKERRQ(ierr);
  ierr = PetscArrayzero(A+ds->l*ld,n1*ld);CHKERRQ(ierr);
  for (i=ds->l;i<ds->n;i++) {
    if (wi) wi[i] = 0.0;
    B[i+i*ld] = 1.0;
    A[i+i*ld] = wr[i];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSynchronize_GHEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscErrorCode ierr;
  PetscInt       ld=ds->ld,l=ds->l,k;
  PetscMPIInt    n,rank,off=0,size,ldn;

  PetscFunctionBegin;
  k = 2*(ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += (ds->n-l)*ld;
  if (eigr) k += (ds->n-l);
  ierr = DSAllocateWork_Private(ds,k,0,0);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(k*sizeof(PetscScalar),&size);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ds->n-l,&n);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ld*(ds->n-l),&ldn);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank);CHKERRQ(ierr);
  if (!rank) {
    ierr = MPI_Pack(ds->mat[DS_MAT_A]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    ierr = MPI_Pack(ds->mat[DS_MAT_B]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    if (ds->state>DS_STATE_RAW) {
      ierr = MPI_Pack(ds->mat[DS_MAT_Q]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigr) {
      ierr = MPI_Pack(eigr+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
  }
  ierr = MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
  if (rank) {
    ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_A]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_B]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    if (ds->state>DS_STATE_RAW) {
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_Q]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigr) {
      ierr = MPI_Unpack(ds->work,size,&off,eigr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSHermitian_GHEP(DS ds,DSMatType m,PetscBool *flg)
{
  PetscFunctionBegin;
  if (m==DS_MAT_A || m==DS_MAT_B) *flg = PETSC_TRUE;
  else *flg = PETSC_FALSE;
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode DSCreate_GHEP(DS ds)
{
  PetscFunctionBegin;
  ds->ops->allocate      = DSAllocate_GHEP;
  ds->ops->view          = DSView_GHEP;
  ds->ops->vectors       = DSVectors_GHEP;
  ds->ops->solve[0]      = DSSolve_GHEP;
  ds->ops->sort          = DSSort_GHEP;
  ds->ops->synchronize   = DSSynchronize_GHEP;
  ds->ops->hermitian     = DSHermitian_GHEP;
  PetscFunctionReturn(0);
}

