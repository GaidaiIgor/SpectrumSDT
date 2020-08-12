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

PetscErrorCode DSAllocate_SVD(DS ds,PetscInt ld)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DSAllocateMat_Private(ds,DS_MAT_A);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_U);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_VT);CHKERRQ(ierr);
  ierr = DSAllocateMatReal_Private(ds,DS_MAT_T);CHKERRQ(ierr);
  ierr = PetscFree(ds->perm);CHKERRQ(ierr);
  ierr = PetscMalloc1(ld,&ds->perm);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)ds,ld*sizeof(PetscInt));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*   0       l           k                 n-1
    -----------------------------------------
    |*       .           .                  |
    |  *     .           .                  |
    |    *   .           .                  |
    |      * .           .                  |
    |        o           o                  |
    |          o         o                  |
    |            o       o                  |
    |              o     o                  |
    |                o   o                  |
    |                  o o                  |
    |                    o x                |
    |                      x x              |
    |                        x x            |
    |                          x x          |
    |                            x x        |
    |                              x x      |
    |                                x x    |
    |                                  x x  |
    |                                    x x|
    |                                      x|
    -----------------------------------------
*/

static PetscErrorCode DSSwitchFormat_SVD(DS ds)
{
  PetscErrorCode ierr;
  PetscReal      *T = ds->rmat[DS_MAT_T];
  PetscScalar    *A = ds->mat[DS_MAT_A];
  PetscInt       i,m=ds->m,k=ds->k,ld=ds->ld;

  PetscFunctionBegin;
  if (!m) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"m was not set");
  /* switch from compact (arrow) to dense storage */
  ierr = PetscArrayzero(A,ld*ld);CHKERRQ(ierr);
  for (i=0;i<k;i++) {
    A[i+i*ld] = T[i];
    A[i+k*ld] = T[i+ld];
  }
  A[k+k*ld] = T[k];
  for (i=k+1;i<m;i++) {
    A[i+i*ld]   = T[i];
    A[i-1+i*ld] = T[i-1+ld];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSView_SVD(DS ds,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  PetscViewerFormat format;
  PetscInt          i,j,r,c;
  PetscReal         value;

  PetscFunctionBegin;
  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
  if (ds->compact) {
    if (!ds->m) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"m was not set");
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      ierr = PetscViewerASCIIPrintf(viewer,"%% Size = %D %D\n",ds->n,ds->m);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"zzz = zeros(%D,3);\n",2*ds->n);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"zzz = [\n");CHKERRQ(ierr);
      for (i=0;i<PetscMin(ds->n,ds->m);i++) {
        ierr = PetscViewerASCIIPrintf(viewer,"%D %D  %18.16e\n",i+1,i+1,(double)*(ds->rmat[DS_MAT_T]+i));CHKERRQ(ierr);
      }
      for (i=0;i<PetscMin(ds->n,ds->m)-1;i++) {
        r = PetscMax(i+2,ds->k+1);
        c = i+1;
        ierr = PetscViewerASCIIPrintf(viewer,"%D %D  %18.16e\n",c,r,(double)*(ds->rmat[DS_MAT_T]+ds->ld+i));CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer,"];\n%s = spconvert(zzz);\n",DSMatName[DS_MAT_T]);CHKERRQ(ierr);
    } else {
      for (i=0;i<ds->n;i++) {
        for (j=0;j<ds->m;j++) {
          if (i==j) value = *(ds->rmat[DS_MAT_T]+i);
          else if (i<ds->k && j==ds->k) value = *(ds->rmat[DS_MAT_T]+ds->ld+PetscMin(i,j));
          else if (i==j+1 && i>ds->k) value = *(ds->rmat[DS_MAT_T]+ds->ld+i-1);
          else value = 0.0;
          ierr = PetscViewerASCIIPrintf(viewer," %18.16e ",(double)value);CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
      }
    }
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  } else {
    ierr = DSViewMat(ds,viewer,DS_MAT_A);CHKERRQ(ierr);
  }
  if (ds->state>DS_STATE_INTERMEDIATE) {
    ierr = DSViewMat(ds,viewer,DS_MAT_U);CHKERRQ(ierr);
    ierr = DSViewMat(ds,viewer,DS_MAT_VT);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSVectors_SVD(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_U:
    case DS_MAT_VT:
      if (rnorm) *rnorm = 0.0;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSort_SVD(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscErrorCode ierr;
  PetscInt       n,l,i,*perm,ld=ds->ld;
  PetscScalar    *A;
  PetscReal      *d;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(0);
  l = ds->l;
  n = PetscMin(ds->n,ds->m);
  A = ds->mat[DS_MAT_A];
  d = ds->rmat[DS_MAT_T];
  perm = ds->perm;
  if (!rr) {
    ierr = DSSortEigenvaluesReal_Private(ds,d,perm);CHKERRQ(ierr);
  } else {
    ierr = DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_FALSE);CHKERRQ(ierr);
  }
  for (i=l;i<n;i++) wr[i] = d[perm[i]];
  ierr = DSPermuteBoth_Private(ds,l,n,DS_MAT_U,DS_MAT_VT,perm);CHKERRQ(ierr);
  for (i=l;i<n;i++) d[i] = PetscRealPart(wr[i]);
  if (!ds->compact) {
    for (i=l;i<n;i++) A[i+i*ld] = wr[i];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSolve_SVD_DC(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBLASInt   n1,n2,n3,m2,m3,info,l,n,m,nm,ld,off,lwork;
  PetscScalar    *A,*U,*VT,qwork;
  PetscReal      *d,*e,*Ur,*VTr;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->m,&m);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->l,&l);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->k-l+1,&n1);CHKERRQ(ierr); /* size of leading block, excl. locked */
  ierr = PetscBLASIntCast(n-ds->k-1,&n2);CHKERRQ(ierr); /* size of trailing block */
  ierr = PetscBLASIntCast(m-ds->k-1,&m2);CHKERRQ(ierr);
  n3 = n1+n2;
  m3 = n1+m2;
  off = l+l*ld;
  A  = ds->mat[DS_MAT_A];
  U  = ds->mat[DS_MAT_U];
  VT = ds->mat[DS_MAT_VT];
  d  = ds->rmat[DS_MAT_T];
  e  = ds->rmat[DS_MAT_T]+ld;
  ierr = PetscArrayzero(U,ld*ld);CHKERRQ(ierr);
  for (i=0;i<l;i++) U[i+i*ld] = 1.0;
  ierr = PetscArrayzero(VT,ld*ld);CHKERRQ(ierr);
  for (i=0;i<l;i++) VT[i+i*ld] = 1.0;

  if (ds->state>DS_STATE_RAW) {
    /* solve bidiagonal SVD problem */
    for (i=0;i<l;i++) wr[i] = d[i];
    ierr = DSAllocateWork_Private(ds,0,3*ld*ld+4*ld,8*ld);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = DSAllocateMatReal_Private(ds,DS_MAT_U);CHKERRQ(ierr);
    ierr = DSAllocateMatReal_Private(ds,DS_MAT_VT);CHKERRQ(ierr);
    Ur  = ds->rmat[DS_MAT_U];
    VTr = ds->rmat[DS_MAT_VT];
#else
    Ur  = U;
    VTr = VT;
#endif
    PetscStackCallBLAS("LAPACKbdsdc",LAPACKbdsdc_("U","I",&n3,d+l,e+l,Ur+off,&ld,VTr+off,&ld,NULL,NULL,ds->rwork,ds->iwork,&info));
    SlepcCheckLapackInfo("bdsdc",info);
#if defined(PETSC_USE_COMPLEX)
    for (i=l;i<n;i++) {
      for (j=0;j<n;j++) {
        U[i+j*ld] = Ur[i+j*ld];
        VT[i+j*ld] = VTr[i+j*ld];
      }
    }
#endif
  } else {
    /* solve general rectangular SVD problem */
    if (ds->compact) { ierr = DSSwitchFormat_SVD(ds);CHKERRQ(ierr); }
    for (i=0;i<l;i++) wr[i] = d[i];
    nm = PetscMin(n,m);
    ierr = DSAllocateWork_Private(ds,0,0,8*nm);CHKERRQ(ierr);
    lwork = -1;
#if defined(PETSC_USE_COMPLEX)
    ierr = DSAllocateWork_Private(ds,0,5*nm*nm+7*nm,0);CHKERRQ(ierr);
    PetscStackCallBLAS("LAPACKgesdd",LAPACKgesdd_("A",&n3,&m3,A+off,&ld,d+l,U+off,&ld,VT+off,&ld,&qwork,&lwork,ds->rwork,ds->iwork,&info));
#else
    PetscStackCallBLAS("LAPACKgesdd",LAPACKgesdd_("A",&n3,&m3,A+off,&ld,d+l,U+off,&ld,VT+off,&ld,&qwork,&lwork,ds->iwork,&info));
#endif
    SlepcCheckLapackInfo("gesdd",info);
    ierr = PetscBLASIntCast((PetscInt)PetscRealPart(qwork),&lwork);CHKERRQ(ierr);
    ierr = DSAllocateWork_Private(ds,lwork,0,0);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    PetscStackCallBLAS("LAPACKgesdd",LAPACKgesdd_("A",&n3,&m3,A+off,&ld,d+l,U+off,&ld,VT+off,&ld,ds->work,&lwork,ds->rwork,ds->iwork,&info));
#else
    PetscStackCallBLAS("LAPACKgesdd",LAPACKgesdd_("A",&n3,&m3,A+off,&ld,d+l,U+off,&ld,VT+off,&ld,ds->work,&lwork,ds->iwork,&info));
#endif
    SlepcCheckLapackInfo("gesdd",info);
  }
  for (i=l;i<PetscMin(ds->n,ds->m);i++) wr[i] = d[i];

  /* create diagonal matrix as a result */
  if (ds->compact) {
    ierr = PetscArrayzero(e,n-1);CHKERRQ(ierr);
  } else {
    for (i=l;i<n;i++) {
      ierr = PetscArrayzero(A+l+i*ld,n-l);CHKERRQ(ierr);
    }
    for (i=l;i<n;i++) A[i+i*ld] = d[i];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSynchronize_SVD(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscErrorCode ierr;
  PetscInt       ld=ds->ld,l=ds->l,k=0,kr=0;
  PetscMPIInt    n,rank,off=0,size,ldn,ld3;

  PetscFunctionBegin;
  if (ds->compact) kr = 3*ld;
  else k = (ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += 2*(ds->n-l)*ld;
  if (eigr) k += ds->n-l;
  ierr = DSAllocateWork_Private(ds,k+kr,0,0);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(k*sizeof(PetscScalar)+kr*sizeof(PetscReal),&size);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ds->n-l,&n);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ld*(ds->n-l),&ldn);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(3*ld,&ld3);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank);CHKERRQ(ierr);
  if (!rank) {
    if (ds->compact) {
      ierr = MPI_Pack(ds->rmat[DS_MAT_T],ld3,MPIU_REAL,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    } else {
      ierr = MPI_Pack(ds->mat[DS_MAT_A]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (ds->state>DS_STATE_RAW) {
      ierr = MPI_Pack(ds->mat[DS_MAT_U]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
      ierr = MPI_Pack(ds->mat[DS_MAT_VT]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigr) {
      ierr = MPI_Pack(eigr+l,n,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
  }
  ierr = MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
  if (rank) {
    if (ds->compact) {
      ierr = MPI_Unpack(ds->work,size,&off,ds->rmat[DS_MAT_T],ld3,MPIU_REAL,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    } else {
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_A]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (ds->state>DS_STATE_RAW) {
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_U]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_VT]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigr) {
      ierr = MPI_Unpack(ds->work,size,&off,eigr+l,n,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSMatGetSize_SVD(DS ds,DSMatType t,PetscInt *rows,PetscInt *cols)
{
  PetscFunctionBegin;
  switch (t) {
    case DS_MAT_A:
    case DS_MAT_T:
      *rows = ds->n;
      *cols = ds->m;
      break;
    case DS_MAT_U:
      *rows = ds->n;
      *cols = ds->n;
      break;
    case DS_MAT_VT:
      *rows = ds->m;
      *cols = ds->m;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid t parameter");
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode DSCreate_SVD(DS ds)
{
  PetscFunctionBegin;
  ds->ops->allocate      = DSAllocate_SVD;
  ds->ops->view          = DSView_SVD;
  ds->ops->vectors       = DSVectors_SVD;
  ds->ops->solve[0]      = DSSolve_SVD_DC;
  ds->ops->sort          = DSSort_SVD;
  ds->ops->synchronize   = DSSynchronize_SVD;
  ds->ops->matgetsize    = DSMatGetSize_SVD;
  PetscFunctionReturn(0);
}

