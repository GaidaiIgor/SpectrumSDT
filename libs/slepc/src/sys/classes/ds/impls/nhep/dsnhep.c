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

PetscErrorCode DSAllocate_NHEP(DS ds,PetscInt ld)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DSAllocateMat_Private(ds,DS_MAT_A);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_Q);CHKERRQ(ierr);
  ierr = PetscFree(ds->perm);CHKERRQ(ierr);
  ierr = PetscMalloc1(ld,&ds->perm);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)ds,ld*sizeof(PetscInt));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DSView_NHEP(DS ds,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  PetscViewerFormat format;

  PetscFunctionBegin;
  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
  ierr = DSViewMat(ds,viewer,DS_MAT_A);CHKERRQ(ierr);
  if (ds->state>DS_STATE_INTERMEDIATE) { ierr = DSViewMat(ds,viewer,DS_MAT_Q);CHKERRQ(ierr); }
  if (ds->mat[DS_MAT_X]) { ierr = DSViewMat(ds,viewer,DS_MAT_X);CHKERRQ(ierr); }
  if (ds->mat[DS_MAT_Y]) { ierr = DSViewMat(ds,viewer,DS_MAT_Y);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSVectors_NHEP_Refined_Some(DS ds,PetscInt *k,PetscReal *rnorm,PetscBool left)
{
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscBLASInt   info,ld,n,n1,lwork,inc=1;
  PetscScalar    sdummy,done=1.0,zero=0.0;
  PetscReal      *sigma;
  PetscBool      iscomplex = PETSC_FALSE;
  PetscScalar    *A = ds->mat[DS_MAT_A];
  PetscScalar    *Q = ds->mat[DS_MAT_Q];
  PetscScalar    *X = ds->mat[left?DS_MAT_Y:DS_MAT_X];
  PetscScalar    *W;

  PetscFunctionBegin;
  if (left) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented for left vectors");
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  n1 = n+1;
  if ((*k)<n-1 && A[(*k)+1+(*k)*ld]!=0.0) iscomplex = PETSC_TRUE;
  if (iscomplex) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Not implemented for complex eigenvalues yet");
  ierr = DSAllocateWork_Private(ds,5*ld,6*ld,0);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_W);CHKERRQ(ierr);
  W = ds->mat[DS_MAT_W];
  lwork = 5*ld;
  sigma = ds->rwork+5*ld;

  /* build A-w*I in W */
  for (j=0;j<n;j++)
    for (i=0;i<=n;i++)
      W[i+j*ld] = A[i+j*ld];
  for (i=0;i<n;i++)
    W[i+i*ld] -= A[(*k)+(*k)*ld];

  /* compute SVD of W */
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","O",&n1,&n,W,&ld,sigma,&sdummy,&ld,&sdummy,&ld,ds->work,&lwork,&info));
#else
  PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","O",&n1,&n,W,&ld,sigma,&sdummy,&ld,&sdummy,&ld,ds->work,&lwork,ds->rwork,&info));
#endif
  SlepcCheckLapackInfo("gesvd",info);

  /* the smallest singular value is the new error estimate */
  if (rnorm) *rnorm = sigma[n-1];

  /* update vector with right singular vector associated to smallest singular value,
     accumulating the transformation matrix Q */
  PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&done,Q,&ld,W+n-1,&ld,&zero,X+(*k)*ld,&inc));
  PetscFunctionReturn(0);
}

static PetscErrorCode DSVectors_NHEP_Refined_All(DS ds,PetscBool left)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  for (i=0;i<ds->n;i++) {
    ierr = DSVectors_NHEP_Refined_Some(ds,&i,NULL,left);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSVectors_NHEP_Eigen_Some(DS ds,PetscInt *k,PetscReal *rnorm,PetscBool left)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBLASInt   mm=1,mout,info,ld,n,*select,inc = 1;
  PetscScalar    tmp,done=1.0,zero=0.0;
  PetscReal      norm;
  PetscBool      iscomplex = PETSC_FALSE;
  PetscScalar    *A = ds->mat[DS_MAT_A];
  PetscScalar    *Q = ds->mat[DS_MAT_Q];
  PetscScalar    *X = ds->mat[left?DS_MAT_Y:DS_MAT_X];
  PetscScalar    *Y;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  ierr = DSAllocateWork_Private(ds,0,0,ld);CHKERRQ(ierr);
  select = ds->iwork;
  for (i=0;i<n;i++) select[i] = (PetscBLASInt)PETSC_FALSE;

  /* compute k-th eigenvector Y of A */
  Y = X+(*k)*ld;
  select[*k] = (PetscBLASInt)PETSC_TRUE;
#if !defined(PETSC_USE_COMPLEX)
  if ((*k)<n-1 && A[(*k)+1+(*k)*ld]!=0.0) iscomplex = PETSC_TRUE;
  mm = iscomplex? 2: 1;
  if (iscomplex) select[(*k)+1] = (PetscBLASInt)PETSC_TRUE;
  ierr = DSAllocateWork_Private(ds,3*ld,0,0);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKtrevc",LAPACKtrevc_(left?"L":"R","S",select,&n,A,&ld,Y,&ld,Y,&ld,&mm,&mout,ds->work,&info));
#else
  ierr = DSAllocateWork_Private(ds,2*ld,ld,0);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKtrevc",LAPACKtrevc_(left?"L":"R","S",select,&n,A,&ld,Y,&ld,Y,&ld,&mm,&mout,ds->work,ds->rwork,&info));
#endif
  SlepcCheckLapackInfo("trevc",info);
  if (mout != mm) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Inconsistent arguments");

  /* accumulate and normalize eigenvectors */
  if (ds->state>=DS_STATE_CONDENSED) {
    ierr = PetscArraycpy(ds->work,Y,mout*ld);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&done,Q,&ld,ds->work,&inc,&zero,Y,&inc));
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&done,Q,&ld,ds->work+ld,&inc,&zero,Y+ld,&inc));
#endif
    norm = BLASnrm2_(&n,Y,&inc);
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) {
      tmp = BLASnrm2_(&n,Y+ld,&inc);
      norm = SlepcAbsEigenvalue(norm,tmp);
    }
#endif
    tmp = 1.0 / norm;
    PetscStackCallBLAS("BLASscal",BLASscal_(&n,&tmp,Y,&inc));
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) PetscStackCallBLAS("BLASscal",BLASscal_(&n,&tmp,Y+ld,&inc));
#endif
  }

  /* set output arguments */
  if (iscomplex) (*k)++;
  if (rnorm) {
    if (iscomplex) *rnorm = SlepcAbsEigenvalue(Y[n-1],Y[n-1+ld]);
    else *rnorm = PetscAbsScalar(Y[n-1]);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSVectors_NHEP_Eigen_All(DS ds,PetscBool left)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBLASInt   n,ld,mout,info,inc = 1;
  PetscBool      iscomplex;
  PetscScalar    *X,*Y,*Z,*A = ds->mat[DS_MAT_A],tmp;
  PetscReal      norm;
  const char     *side,*back;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  if (left) {
    X = NULL;
    Y = ds->mat[DS_MAT_Y];
    side = "L";
  } else {
    X = ds->mat[DS_MAT_X];
    Y = NULL;
    side = "R";
  }
  Z = left? Y: X;
  if (ds->state>=DS_STATE_CONDENSED) {
    /* DSSolve() has been called, backtransform with matrix Q */
    back = "B";
    ierr = PetscArraycpy(Z,ds->mat[DS_MAT_Q],ld*ld);CHKERRQ(ierr);
  } else back = "A";
#if !defined(PETSC_USE_COMPLEX)
  ierr = DSAllocateWork_Private(ds,3*ld,0,0);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKtrevc",LAPACKtrevc_(side,back,NULL,&n,A,&ld,Y,&ld,X,&ld,&n,&mout,ds->work,&info));
#else
  ierr = DSAllocateWork_Private(ds,2*ld,ld,0);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKtrevc",LAPACKtrevc_(side,back,NULL,&n,A,&ld,Y,&ld,X,&ld,&n,&mout,ds->work,ds->rwork,&info));
#endif
  SlepcCheckLapackInfo("trevc",info);

  /* normalize eigenvectors */
  for (i=0;i<n;i++) {
    iscomplex = (i<n-1 && A[i+1+i*ld]!=0.0)? PETSC_TRUE: PETSC_FALSE;
    norm = BLASnrm2_(&n,Z+i*ld,&inc);
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) {
      tmp = BLASnrm2_(&n,Z+(i+1)*ld,&inc);
      norm = SlepcAbsEigenvalue(norm,tmp);
    }
#endif
    tmp = 1.0 / norm;
    PetscStackCallBLAS("BLASscal",BLASscal_(&n,&tmp,Z+i*ld,&inc));
#if !defined(PETSC_USE_COMPLEX)
    if (iscomplex) PetscStackCallBLAS("BLASscal",BLASscal_(&n,&tmp,Z+(i+1)*ld,&inc));
#endif
    if (iscomplex) i++;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSVectors_NHEP(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_X:
      if (ds->refined) {
        if (!ds->extrarow) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Refined vectors require activating the extra row");
        if (j) {
          ierr = DSVectors_NHEP_Refined_Some(ds,j,rnorm,PETSC_FALSE);CHKERRQ(ierr);
        } else {
          ierr = DSVectors_NHEP_Refined_All(ds,PETSC_FALSE);CHKERRQ(ierr);
        }
      } else {
        if (j) {
          ierr = DSVectors_NHEP_Eigen_Some(ds,j,rnorm,PETSC_FALSE);CHKERRQ(ierr);
        } else {
          ierr = DSVectors_NHEP_Eigen_All(ds,PETSC_FALSE);CHKERRQ(ierr);
        }
      }
      break;
    case DS_MAT_Y:
      if (ds->refined) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
      if (j) {
        ierr = DSVectors_NHEP_Eigen_Some(ds,j,rnorm,PETSC_TRUE);CHKERRQ(ierr);
      } else {
        ierr = DSVectors_NHEP_Eigen_All(ds,PETSC_TRUE);CHKERRQ(ierr);
      }
      break;
    case DS_MAT_U:
    case DS_MAT_VT:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  if (ds->state < DS_STATE_CONDENSED) {
    ierr = DSSetState(ds,DS_STATE_CONDENSED);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSSort_NHEP_Arbitrary(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBLASInt   info,n,ld,mout,lwork,*selection;
  PetscScalar    *T = ds->mat[DS_MAT_A],*Q = ds->mat[DS_MAT_Q],*work;
#if !defined(PETSC_USE_COMPLEX)
  PetscBLASInt   *iwork,liwork;
#endif

  PetscFunctionBegin;
  if (!k) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"Must supply argument k");
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  lwork = n;
  liwork = 1;
  ierr = DSAllocateWork_Private(ds,lwork,0,liwork+n);CHKERRQ(ierr);
  work = ds->work;
  lwork = ds->lwork;
  selection = ds->iwork;
  iwork = ds->iwork + n;
  liwork = ds->liwork - n;
#else
  lwork = 1;
  ierr = DSAllocateWork_Private(ds,lwork,0,n);CHKERRQ(ierr);
  work = ds->work;
  selection = ds->iwork;
#endif
  /* Compute the selected eigenvalue to be in the leading position */
  ierr = DSSortEigenvalues_Private(ds,rr,ri,ds->perm,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscArrayzero(selection,n);CHKERRQ(ierr);
  for (i=0;i<*k;i++) selection[ds->perm[i]] = 1;
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKtrsen",LAPACKtrsen_("N","V",selection,&n,T,&ld,Q,&ld,wr,wi,&mout,NULL,NULL,work,&lwork,iwork,&liwork,&info));
#else
  PetscStackCallBLAS("LAPACKtrsen",LAPACKtrsen_("N","V",selection,&n,T,&ld,Q,&ld,wr,&mout,NULL,NULL,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("trsen",info);
  *k = mout;
  PetscFunctionReturn(0);
}

static PetscErrorCode DSSort_NHEP_Total(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscScalar    re;
  PetscInt       i,j,pos,result;
  PetscBLASInt   ifst,ilst,info,n,ld;
  PetscScalar    *T = ds->mat[DS_MAT_A];
  PetscScalar    *Q = ds->mat[DS_MAT_Q];
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    *work,im;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  ierr = DSAllocateWork_Private(ds,ld,0,0);CHKERRQ(ierr);
  work = ds->work;
#endif
  /* selection sort */
  for (i=ds->l;i<n-1;i++) {
    re = wr[i];
#if !defined(PETSC_USE_COMPLEX)
    im = wi[i];
#endif
    pos = 0;
    j=i+1; /* j points to the next eigenvalue */
#if !defined(PETSC_USE_COMPLEX)
    if (im != 0) j=i+2;
#endif
    /* find minimum eigenvalue */
    for (;j<n;j++) {
#if !defined(PETSC_USE_COMPLEX)
      ierr = SlepcSCCompare(ds->sc,re,im,wr[j],wi[j],&result);CHKERRQ(ierr);
#else
      ierr = SlepcSCCompare(ds->sc,re,0.0,wr[j],0.0,&result);CHKERRQ(ierr);
#endif
      if (result > 0) {
        re = wr[j];
#if !defined(PETSC_USE_COMPLEX)
        im = wi[j];
#endif
        pos = j;
      }
#if !defined(PETSC_USE_COMPLEX)
      if (wi[j] != 0) j++;
#endif
    }
    if (pos) {
      /* interchange blocks */
      ierr = PetscBLASIntCast(pos+1,&ifst);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(i+1,&ilst);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
      PetscStackCallBLAS("LAPACKtrexc",LAPACKtrexc_("V",&n,T,&ld,Q,&ld,&ifst,&ilst,work,&info));
#else
      PetscStackCallBLAS("LAPACKtrexc",LAPACKtrexc_("V",&n,T,&ld,Q,&ld,&ifst,&ilst,&info));
#endif
      SlepcCheckLapackInfo("trexc",info);
      /* recover original eigenvalues from T matrix */
      for (j=i;j<n;j++) {
        wr[j] = T[j+j*ld];
#if !defined(PETSC_USE_COMPLEX)
        if (j<n-1 && T[j+1+j*ld] != 0.0) {
          /* complex conjugate eigenvalue */
          wi[j] = PetscSqrtReal(PetscAbsReal(T[j+1+j*ld])) *
                  PetscSqrtReal(PetscAbsReal(T[j+(j+1)*ld]));
          wr[j+1] = wr[j];
          wi[j+1] = -wi[j];
          j++;
        } else {
          wi[j] = 0.0;
        }
#endif
      }
    }
#if !defined(PETSC_USE_COMPLEX)
    if (wi[i] != 0) i++;
#endif
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSort_NHEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!rr || wr == rr) {
    ierr = DSSort_NHEP_Total(ds,wr,wi);CHKERRQ(ierr);
  } else {
    ierr = DSSort_NHEP_Arbitrary(ds,wr,wi,rr,ri,k);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSSortWithPermutation_NHEP(DS ds,PetscInt *perm,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscInt       i,j,pos,inc=1;
  PetscBLASInt   ifst,ilst,info,n,ld;
  PetscScalar    *T = ds->mat[DS_MAT_A];
  PetscScalar    *Q = ds->mat[DS_MAT_Q];
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    *work;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  ierr = DSAllocateWork_Private(ds,ld,0,0);CHKERRQ(ierr);
  work = ds->work;
#endif
  for (i=ds->l;i<n-1;i++) {
    pos = perm[i];
#if !defined(PETSC_USE_COMPLEX)
    inc = (pos<n-1 && T[pos+1+pos*ld] != 0.0)? 2: 1;
#endif
    if (pos!=i) {
#if !defined(PETSC_USE_COMPLEX)
      if ((T[pos+(pos-1)*ld] != 0.0 && perm[i+1]!=pos-1) || (T[pos+1+pos*ld] != 0.0 && perm[i+1]!=pos+1))
 SETERRQ1(PETSC_COMM_SELF,1,"Invalid permutation due to a 2x2 block at position %D",pos);
#endif

      /* interchange blocks */
      ierr = PetscBLASIntCast(pos+1,&ifst);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(i+1,&ilst);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
      PetscStackCallBLAS("LAPACKtrexc",LAPACKtrexc_("V",&n,T,&ld,Q,&ld,&ifst,&ilst,work,&info));
#else
      PetscStackCallBLAS("LAPACKtrexc",LAPACKtrexc_("V",&n,T,&ld,Q,&ld,&ifst,&ilst,&info));
#endif
      SlepcCheckLapackInfo("trexc",info);
      for (j=i+1;j<n;j++) {
        if (perm[j]>=i && perm[j]<pos) perm[j]+=inc;
      }
      perm[i] = i;
      if (inc==2) perm[i+1] = i+1;
    }
    if (inc==2) i++;
  }
  /* recover original eigenvalues from T matrix */
  for (j=ds->l;j<n;j++) {
    wr[j] = T[j+j*ld];
#if !defined(PETSC_USE_COMPLEX)
    if (j<n-1 && T[j+1+j*ld] != 0.0) {
      /* complex conjugate eigenvalue */
      wi[j] = PetscSqrtReal(PetscAbsReal(T[j+1+j*ld])) * PetscSqrtReal(PetscAbsReal(T[j+(j+1)*ld]));
      wr[j+1] = wr[j];
      wi[j+1] = -wi[j];
      j++;
    } else {
      wi[j] = 0.0;
    }
#endif
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSUpdateExtraRow_NHEP(DS ds)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBLASInt   n,ld,incx=1;
  PetscScalar    *A,*Q,*x,*y,one=1.0,zero=0.0;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  A  = ds->mat[DS_MAT_A];
  Q  = ds->mat[DS_MAT_Q];
  ierr = DSAllocateWork_Private(ds,2*ld,0,0);CHKERRQ(ierr);
  x = ds->work;
  y = ds->work+ld;
  for (i=0;i<n;i++) x[i] = PetscConj(A[n+i*ld]);
  PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&one,Q,&ld,x,&incx,&zero,y,&incx));
  for (i=0;i<n;i++) A[n+i*ld] = PetscConj(y[i]);
  ds->k = n;
  PetscFunctionReturn(0);
}

/*
   Reduce a matrix A to upper Hessenberg form Q'*A*Q, where Q is an orthogonal matrix.
   The result overwrites A. Matrix A has the form

         [ S | * ]
     A = [-------]
         [ R | H ]

  where S is an upper (quasi-)triangular matrix of order k, H is an upper Hessenberg
  matrix of order n-k, and R is all zeros except the first row (the arrow).
  The algorithm uses elementary reflectors to annihilate entries in the arrow, and
  then proceeds upwards.
  If ilo>1, then it is assumed that the first ilo-1 entries of the arrow are zero, and
  hence the first ilo-1 rows and columns of Q are set to the identity matrix.

  Required workspace is 2*n.
*/
static PetscErrorCode ArrowHessenberg(PetscBLASInt n,PetscBLASInt k,PetscBLASInt ilo,PetscScalar *A,PetscBLASInt lda,PetscScalar *Q,PetscBLASInt ldq,PetscScalar *work)
{
  PetscBLASInt i,j,n0,m,inc=1,incn=-1;
  PetscScalar  t,*v=work,*w=work+n,tau,tauc;

  PetscFunctionBegin;
  m = n-ilo+1;
  for (i=k;i>ilo;i--) {
    for (j=0;j<=i-ilo;j++) v[j] = A[i+(i-j-1)*lda];  /* _larfg does not allow negative inc, so use buffer */
    n0 = i-ilo+1;
    PetscStackCallBLAS("LAPACKlarfg",LAPACKlarfg_(&n0,v,v+1,&inc,&tau));
    for (j=1;j<=i-ilo;j++) v[j] = PetscConj(v[j]);
    tauc = PetscConj(tau);
    A[i+(i-1)*lda] = v[0];
    v[0] = 1.0;
    PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("R",&i,&n0,v,&incn,&tauc,A+(ilo-1)*lda,&lda,w));
    for (j=1;j<=i-ilo;j++) A[i+(i-j-1)*lda] = 0.0;
    PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n0,&m,v,&incn,&tau,A+ilo-1+(ilo-1)*lda,&lda,w));
    PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n0,&m,v,&incn,&tau,Q+ilo-1+(ilo-1)*ldq,&ldq,w));
  }
  /* trivial in-place transposition of Q */
  for (j=ilo-1;j<k;j++) {
    for (i=j;i<k;i++) {
      t = Q[i+j*ldq];
      if (i!=j) Q[i+j*ldq] = PetscConj(Q[j+i*ldq]);
      Q[j+i*ldq] = PetscConj(t);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSolve_NHEP(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscScalar    *work,*tau;
  PetscInt       i,j;
  PetscBLASInt   ilo,lwork,info,n,k,ld;
  PetscScalar    *A = ds->mat[DS_MAT_A];
  PetscScalar    *Q = ds->mat[DS_MAT_Q];

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscValidScalarPointer(wi,3);
#endif
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->l+1,&ilo);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->k,&k);CHKERRQ(ierr);
  ierr = DSAllocateWork_Private(ds,ld+6*ld,0,0);CHKERRQ(ierr);
  tau  = ds->work;
  work = ds->work+ld;
  lwork = 6*ld;

  /* initialize orthogonal matrix */
  ierr = PetscArrayzero(Q,ld*ld);CHKERRQ(ierr);
  for (i=0;i<n;i++) Q[i+i*ld] = 1.0;
  if (n==1) { /* quick return */
    wr[0] = A[0];
    wi[0] = 0.0;
    PetscFunctionReturn(0);
  }

  /* reduce to upper Hessenberg form */
  if (ds->state<DS_STATE_INTERMEDIATE) {
    if (PETSC_FALSE && k>0) {
      ierr = ArrowHessenberg(n,k,ilo,A,ld,Q,ld,work);CHKERRQ(ierr);
    } else {
      PetscStackCallBLAS("LAPACKgehrd",LAPACKgehrd_(&n,&ilo,&n,A,&ld,tau,work,&lwork,&info));
      SlepcCheckLapackInfo("gehrd",info);
      for (j=0;j<n-1;j++) {
        for (i=j+2;i<n;i++) {
          Q[i+j*ld] = A[i+j*ld];
          A[i+j*ld] = 0.0;
        }
      }
      PetscStackCallBLAS("LAPACKorghr",LAPACKorghr_(&n,&ilo,&n,Q,&ld,tau,work,&lwork,&info));
      SlepcCheckLapackInfo("orghr",info);
    }
  }

  /* compute the (real) Schur form */
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKhseqr",LAPACKhseqr_("S","V",&n,&ilo,&n,A,&ld,wr,wi,Q,&ld,work,&lwork,&info));
  for (j=0;j<ds->l;j++) {
    if (j==n-1 || A[j+1+j*ld] == 0.0) {
      /* real eigenvalue */
      wr[j] = A[j+j*ld];
      wi[j] = 0.0;
    } else {
      /* complex eigenvalue */
      wr[j] = A[j+j*ld];
      wr[j+1] = A[j+j*ld];
      wi[j] = PetscSqrtReal(PetscAbsReal(A[j+1+j*ld])) *
              PetscSqrtReal(PetscAbsReal(A[j+(j+1)*ld]));
      wi[j+1] = -wi[j];
      j++;
    }
  }
#else
  PetscStackCallBLAS("LAPACKhseqr",LAPACKhseqr_("S","V",&n,&ilo,&n,A,&ld,wr,Q,&ld,work,&lwork,&info));
  if (wi) for (i=ds->l;i<n;i++) wi[i] = 0.0;
#endif
  SlepcCheckLapackInfo("hseqr",info);
  PetscFunctionReturn(0);
}

PetscErrorCode DSSynchronize_NHEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscErrorCode ierr;
  PetscInt       ld=ds->ld,l=ds->l,k;
  PetscMPIInt    n,rank,off=0,size,ldn;

  PetscFunctionBegin;
  k = (ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += (ds->n-l)*ld;
  if (eigr) k += ds->n-l;
  if (eigi) k += ds->n-l;
  ierr = DSAllocateWork_Private(ds,k,0,0);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(k*sizeof(PetscScalar),&size);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ds->n-l,&n);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ld*(ds->n-l),&ldn);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank);CHKERRQ(ierr);
  if (!rank) {
    ierr = MPI_Pack(ds->mat[DS_MAT_A]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    if (ds->state>DS_STATE_RAW) {
      ierr = MPI_Pack(ds->mat[DS_MAT_Q]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigr) {
      ierr = MPI_Pack(eigr+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigi) {
      ierr = MPI_Pack(eigi+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
  }
  ierr = MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
  if (rank) {
    ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_A]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    if (ds->state>DS_STATE_RAW) {
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_Q]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigr) {
      ierr = MPI_Unpack(ds->work,size,&off,eigr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigi) {
      ierr = MPI_Unpack(ds->work,size,&off,eigi+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSTruncate_NHEP(DS ds,PetscInt n)
{
  PetscInt    i,newn,ld=ds->ld,l=ds->l;
  PetscScalar *A;

  PetscFunctionBegin;
  if (ds->state==DS_STATE_CONDENSED) ds->t = ds->n;
  A = ds->mat[DS_MAT_A];
  /* be careful not to break a diagonal 2x2 block */
  if (A[n+(n-1)*ld]==0.0) newn = n;
  else {
    if (n<ds->n-1) newn = n+1;
    else newn = n-1;
  }
  if (ds->extrarow && ds->k==ds->n) {
    /* copy entries of extra row to the new position, then clean last row */
    for (i=l;i<newn;i++) A[newn+i*ld] = A[ds->n+i*ld];
    for (i=l;i<ds->n;i++) A[ds->n+i*ld] = 0.0;
  }
  ds->k = 0;
  ds->n = newn;
  PetscFunctionReturn(0);
}

PetscErrorCode DSCond_NHEP(DS ds,PetscReal *cond)
{
  PetscErrorCode ierr;
  PetscScalar    *work;
  PetscReal      *rwork;
  PetscBLASInt   *ipiv;
  PetscBLASInt   lwork,info,n,ld;
  PetscReal      hn,hin;
  PetscScalar    *A;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  lwork = 8*ld;
  ierr = DSAllocateWork_Private(ds,lwork,ld,ld);CHKERRQ(ierr);
  work  = ds->work;
  rwork = ds->rwork;
  ipiv  = ds->iwork;

  /* use workspace matrix W to avoid overwriting A */
  ierr = DSAllocateMat_Private(ds,DS_MAT_W);CHKERRQ(ierr);
  A = ds->mat[DS_MAT_W];
  ierr = PetscArraycpy(A,ds->mat[DS_MAT_A],ds->ld*ds->ld);CHKERRQ(ierr);

  /* norm of A */
  if (ds->state<DS_STATE_INTERMEDIATE) hn = LAPACKlange_("I",&n,&n,A,&ld,rwork);
  else hn = LAPACKlanhs_("I",&n,A,&ld,rwork);

  /* norm of inv(A) */
  PetscStackCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n,&n,A,&ld,ipiv,&info));
  SlepcCheckLapackInfo("getrf",info);
  PetscStackCallBLAS("LAPACKgetri",LAPACKgetri_(&n,A,&ld,ipiv,work,&lwork,&info));
  SlepcCheckLapackInfo("getri",info);
  hin = LAPACKlange_("I",&n,&n,A,&ld,rwork);

  *cond = hn*hin;
  PetscFunctionReturn(0);
}

PetscErrorCode DSTranslateHarmonic_NHEP(DS ds,PetscScalar tau,PetscReal beta,PetscBool recover,PetscScalar *gin,PetscReal *gamma)
{
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscBLASInt   *ipiv,info,n,ld,one=1,ncol;
  PetscScalar    *A,*B,*Q,*g=gin,*ghat;
  PetscScalar    done=1.0,dmone=-1.0,dzero=0.0;
  PetscReal      gnorm;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  A  = ds->mat[DS_MAT_A];

  if (!recover) {

    ierr = DSAllocateWork_Private(ds,0,0,ld);CHKERRQ(ierr);
    ipiv = ds->iwork;
    if (!g) {
      ierr = DSAllocateWork_Private(ds,ld,0,0);CHKERRQ(ierr);
      g = ds->work;
    }
    /* use workspace matrix W to factor A-tau*eye(n) */
    ierr = DSAllocateMat_Private(ds,DS_MAT_W);CHKERRQ(ierr);
    B = ds->mat[DS_MAT_W];
    ierr = PetscArraycpy(B,A,ld*ld);CHKERRQ(ierr);

    /* Vector g initialy stores b = beta*e_n^T */
    ierr = PetscArrayzero(g,n);CHKERRQ(ierr);
    g[n-1] = beta;

    /* g = (A-tau*eye(n))'\b */
    for (i=0;i<n;i++) B[i+i*ld] -= tau;
    PetscStackCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n,&n,B,&ld,ipiv,&info));
    SlepcCheckLapackInfo("getrf",info);
    ierr = PetscLogFlops(2.0*n*n*n/3.0);CHKERRQ(ierr);
    PetscStackCallBLAS("LAPACKgetrs",LAPACKgetrs_("C",&n,&one,B,&ld,ipiv,g,&ld,&info));
    SlepcCheckLapackInfo("getrs",info);
    ierr = PetscLogFlops(2.0*n*n-n);CHKERRQ(ierr);

    /* A = A + g*b' */
    for (i=0;i<n;i++) A[i+(n-1)*ld] += g[i]*beta;

  } else { /* recover */

    ierr = DSAllocateWork_Private(ds,ld,0,0);CHKERRQ(ierr);
    ghat = ds->work;
    Q    = ds->mat[DS_MAT_Q];

    /* g^ = -Q(:,idx)'*g */
    ierr = PetscBLASIntCast(ds->l+ds->k,&ncol);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&n,&ncol,&dmone,Q,&ld,g,&one,&dzero,ghat,&one));

    /* A = A + g^*b' */
    for (i=0;i<ds->l+ds->k;i++)
      for (j=ds->l;j<ds->l+ds->k;j++)
        A[i+j*ld] += ghat[i]*Q[n-1+j*ld]*beta;

    /* g~ = (I-Q(:,idx)*Q(:,idx)')*g = g+Q(:,idx)*g^ */
    PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n,&ncol,&done,Q,&ld,ghat,&one,&done,g,&one));
  }

  /* Compute gamma factor */
  if (gamma) {
    gnorm = 0.0;
    for (i=0;i<n;i++) gnorm = gnorm + PetscRealPart(g[i]*PetscConj(g[i]));
    *gamma = PetscSqrtReal(1.0+gnorm);
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode DSCreate_NHEP(DS ds)
{
  PetscFunctionBegin;
  ds->ops->allocate      = DSAllocate_NHEP;
  ds->ops->view          = DSView_NHEP;
  ds->ops->vectors       = DSVectors_NHEP;
  ds->ops->solve[0]      = DSSolve_NHEP;
  ds->ops->sort          = DSSort_NHEP;
  ds->ops->sortperm      = DSSortWithPermutation_NHEP;
  ds->ops->synchronize   = DSSynchronize_NHEP;
  ds->ops->truncate      = DSTruncate_NHEP;
  ds->ops->update        = DSUpdateExtraRow_NHEP;
  ds->ops->cond          = DSCond_NHEP;
  ds->ops->transharm     = DSTranslateHarmonic_NHEP;
  PetscFunctionReturn(0);
}

