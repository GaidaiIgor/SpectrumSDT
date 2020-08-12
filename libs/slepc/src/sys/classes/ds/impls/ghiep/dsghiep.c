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

PetscErrorCode DSAllocate_GHIEP(DS ds,PetscInt ld)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DSAllocateMat_Private(ds,DS_MAT_A);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_B);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_Q);CHKERRQ(ierr);
  ierr = DSAllocateMatReal_Private(ds,DS_MAT_T);CHKERRQ(ierr);
  ierr = DSAllocateMatReal_Private(ds,DS_MAT_D);CHKERRQ(ierr);
  ierr = PetscFree(ds->perm);CHKERRQ(ierr);
  ierr = PetscMalloc1(ld,&ds->perm);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)ds,ld*sizeof(PetscInt));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DSSwitchFormat_GHIEP(DS ds,PetscBool tocompact)
{
  PetscErrorCode ierr;
  PetscReal      *T,*S;
  PetscScalar    *A,*B;
  PetscInt       i,n,ld;

  PetscFunctionBegin;
  A = ds->mat[DS_MAT_A];
  B = ds->mat[DS_MAT_B];
  T = ds->rmat[DS_MAT_T];
  S = ds->rmat[DS_MAT_D];
  n = ds->n;
  ld = ds->ld;
  if (tocompact) { /* switch from dense (arrow) to compact storage */
    ierr = PetscArrayzero(T,3*ld);CHKERRQ(ierr);
    ierr = PetscArrayzero(S,ld);CHKERRQ(ierr);
    for (i=0;i<n-1;i++) {
      T[i]    = PetscRealPart(A[i+i*ld]);
      T[ld+i] = PetscRealPart(A[i+1+i*ld]);
      S[i]    = PetscRealPart(B[i+i*ld]);
    }
    T[n-1] = PetscRealPart(A[n-1+(n-1)*ld]);
    S[n-1] = PetscRealPart(B[n-1+(n-1)*ld]);
    for (i=ds->l;i<ds->k;i++) T[2*ld+i] = PetscRealPart(A[ds->k+i*ld]);
  } else { /* switch from compact (arrow) to dense storage */
    ierr = PetscArrayzero(A,ld*ld);CHKERRQ(ierr);
    ierr = PetscArrayzero(B,ld*ld);CHKERRQ(ierr);
    for (i=0;i<n-1;i++) {
      A[i+i*ld]     = T[i];
      A[i+1+i*ld]   = T[ld+i];
      A[i+(i+1)*ld] = T[ld+i];
      B[i+i*ld]     = S[i];
    }
    A[n-1+(n-1)*ld] = T[n-1];
    B[n-1+(n-1)*ld] = S[n-1];
    for (i=ds->l;i<ds->k;i++) {
      A[ds->k+i*ld] = T[2*ld+i];
      A[i+ds->k*ld] = T[2*ld+i];
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSView_GHIEP(DS ds,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  PetscViewerFormat format;
  PetscInt          i,j;
  PetscReal         value;
  const char        *methodname[] = {
                     "QR + Inverse Iteration",
                     "HZ method",
                     "QR"
  };
  const int         nmeth=sizeof(methodname)/sizeof(methodname[0]);

  PetscFunctionBegin;
  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
    if (ds->method<nmeth) {
      ierr = PetscViewerASCIIPrintf(viewer,"solving the problem with: %s\n",methodname[ds->method]);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }
  if (ds->compact) {
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      ierr = PetscViewerASCIIPrintf(viewer,"%% Size = %D %D\n",ds->n,ds->n);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"zzz = zeros(%D,3);\n",3*ds->n);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"zzz = [\n");CHKERRQ(ierr);
      for (i=0;i<ds->n;i++) {
        ierr = PetscViewerASCIIPrintf(viewer,"%D %D  %18.16e\n",i+1,i+1,(double)*(ds->rmat[DS_MAT_T]+i));CHKERRQ(ierr);
      }
      for (i=0;i<ds->n-1;i++) {
        if (*(ds->rmat[DS_MAT_T]+ds->ld+i) !=0 && i!=ds->k-1) {
          ierr = PetscViewerASCIIPrintf(viewer,"%D %D  %18.16e\n",i+2,i+1,(double)*(ds->rmat[DS_MAT_T]+ds->ld+i));CHKERRQ(ierr);
          ierr = PetscViewerASCIIPrintf(viewer,"%D %D  %18.16e\n",i+1,i+2,(double)*(ds->rmat[DS_MAT_T]+ds->ld+i));CHKERRQ(ierr);
        }
      }
      for (i = ds->l;i<ds->k;i++) {
        if (*(ds->rmat[DS_MAT_T]+2*ds->ld+i)) {
          ierr = PetscViewerASCIIPrintf(viewer,"%D %D  %18.16e\n",ds->k+1,i+1,(double)*(ds->rmat[DS_MAT_T]+2*ds->ld+i));CHKERRQ(ierr);
          ierr = PetscViewerASCIIPrintf(viewer,"%D %D  %18.16e\n",i+1,ds->k+1,(double)*(ds->rmat[DS_MAT_T]+2*ds->ld+i));CHKERRQ(ierr);
        }
      }
      ierr = PetscViewerASCIIPrintf(viewer,"];\n%s = spconvert(zzz);\n",DSMatName[DS_MAT_A]);CHKERRQ(ierr);

      ierr = PetscViewerASCIIPrintf(viewer,"%% Size = %D %D\n",ds->n,ds->n);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"omega = zeros(%D,3);\n",3*ds->n);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"omega = [\n");CHKERRQ(ierr);
      for (i=0;i<ds->n;i++) {
        ierr = PetscViewerASCIIPrintf(viewer,"%D %D  %18.16e\n",i+1,i+1,(double)*(ds->rmat[DS_MAT_D]+i));CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer,"];\n%s = spconvert(omega);\n",DSMatName[DS_MAT_B]);CHKERRQ(ierr);

    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"T\n");CHKERRQ(ierr);
      for (i=0;i<ds->n;i++) {
        for (j=0;j<ds->n;j++) {
          if (i==j) value = *(ds->rmat[DS_MAT_T]+i);
          else if (i==j+1 || j==i+1) value = *(ds->rmat[DS_MAT_T]+ds->ld+PetscMin(i,j));
          else if ((i<ds->k && j==ds->k) || (i==ds->k && j<ds->k)) value = *(ds->rmat[DS_MAT_T]+2*ds->ld+PetscMin(i,j));
          else value = 0.0;
          ierr = PetscViewerASCIIPrintf(viewer," %18.16e ",(double)value);CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer,"omega\n");CHKERRQ(ierr);
      for (i=0;i<ds->n;i++) {
        for (j=0;j<ds->n;j++) {
          if (i==j) value = *(ds->rmat[DS_MAT_D]+i);
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
    ierr = DSViewMat(ds,viewer,DS_MAT_B);CHKERRQ(ierr);
  }
  if (ds->state>DS_STATE_INTERMEDIATE) { ierr = DSViewMat(ds,viewer,DS_MAT_Q);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSVectors_GHIEP_Eigen_Some(DS ds,PetscInt *idx,PetscReal *rnorm)
{
  PetscErrorCode ierr;
  PetscReal      b[4],M[4],d1,d2,s1,s2,e;
  PetscReal      scal1,scal2,wr1,wr2,wi,ep,norm;
  PetscScalar    *Q,*X,Y[4],alpha,zeroS = 0.0;
  PetscInt       k;
  PetscBLASInt   two = 2,n_,ld,one=1;
#if !defined(PETSC_USE_COMPLEX)
  PetscBLASInt   four=4;
#endif

  PetscFunctionBegin;
  X = ds->mat[DS_MAT_X];
  Q = ds->mat[DS_MAT_Q];
  k = *idx;
  ierr = PetscBLASIntCast(ds->n,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  if (k < ds->n-1) e = (ds->compact)?*(ds->rmat[DS_MAT_T]+ld+k):PetscRealPart(*(ds->mat[DS_MAT_A]+(k+1)+ld*k));
  else e = 0.0;
  if (e == 0.0) { /* Real */
    if (ds->state>=DS_STATE_CONDENSED) {
      ierr = PetscArraycpy(X+k*ld,Q+k*ld,ld);CHKERRQ(ierr);
    } else {
      ierr = PetscArrayzero(X+k*ds->ld,ds->ld);CHKERRQ(ierr);
      X[k+k*ds->ld] = 1.0;
    }
    if (rnorm) *rnorm = PetscAbsScalar(X[ds->n-1+k*ld]);
  } else { /* 2x2 block */
    if (ds->compact) {
      s1 = *(ds->rmat[DS_MAT_D]+k);
      d1 = *(ds->rmat[DS_MAT_T]+k);
      s2 = *(ds->rmat[DS_MAT_D]+k+1);
      d2 = *(ds->rmat[DS_MAT_T]+k+1);
    } else {
      s1 = PetscRealPart(*(ds->mat[DS_MAT_B]+k*ld+k));
      d1 = PetscRealPart(*(ds->mat[DS_MAT_A]+k+k*ld));
      s2 = PetscRealPart(*(ds->mat[DS_MAT_B]+(k+1)*ld+k+1));
      d2 = PetscRealPart(*(ds->mat[DS_MAT_A]+k+1+(k+1)*ld));
    }
    M[0] = d1; M[1] = e; M[2] = e; M[3]= d2;
    b[0] = s1; b[1] = 0.0; b[2] = 0.0; b[3] = s2;
    ep = LAPACKlamch_("S");
    /* Compute eigenvalues of the block */
    PetscStackCallBLAS("LAPACKlag2",LAPACKlag2_(M,&two,b,&two,&ep,&scal1,&scal2,&wr1,&wr2,&wi));
    if (wi==0.0) SETERRQ(PETSC_COMM_SELF,1,"Real block in DSVectors_GHIEP");
    else { /* Complex eigenvalues */
      if (scal1<ep) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"Nearly infinite eigenvalue");
      wr1 /= scal1;
      wi  /= scal1;
#if !defined(PETSC_USE_COMPLEX)
      if (SlepcAbs(s1*d1-wr1,wi)<SlepcAbs(s2*d2-wr1,wi)) {
        Y[0] = wr1-s2*d2; Y[1] = s2*e; Y[2] = wi; Y[3] = 0.0;
      } else {
        Y[0] = s1*e; Y[1] = wr1-s1*d1; Y[2] = 0.0; Y[3] = wi;
      }
      norm = BLASnrm2_(&four,Y,&one);
      norm = 1.0/norm;
      if (ds->state >= DS_STATE_CONDENSED) {
        alpha = norm;
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&two,&two,&alpha,ds->mat[DS_MAT_Q]+k*ld,&ld,Y,&two,&zeroS,X+k*ld,&ld));
        if (rnorm) *rnorm = SlepcAbsEigenvalue(X[ds->n-1+k*ld],X[ds->n-1+(k+1)*ld]);
      } else {
        ierr = PetscArrayzero(X+k*ld,2*ld);CHKERRQ(ierr);
        X[k*ld+k]       = Y[0]*norm;
        X[k*ld+k+1]     = Y[1]*norm;
        X[(k+1)*ld+k]   = Y[2]*norm;
        X[(k+1)*ld+k+1] = Y[3]*norm;
      }
#else
      if (SlepcAbs(s1*d1-wr1,wi)<SlepcAbs(s2*d2-wr1,wi)) {
        Y[0] = PetscCMPLX(wr1-s2*d2,wi);
        Y[1] = s2*e;
      } else {
        Y[0] = s1*e;
        Y[1] = PetscCMPLX(wr1-s1*d1,wi);
      }
      norm = BLASnrm2_(&two,Y,&one);
      norm = 1.0/norm;
      if (ds->state >= DS_STATE_CONDENSED) {
        alpha = norm;
        PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&two,&alpha,ds->mat[DS_MAT_Q]+k*ld,&ld,Y,&one,&zeroS,X+k*ld,&one));
        if (rnorm) *rnorm = PetscAbsScalar(X[ds->n-1+k*ld]);
      } else {
        ierr = PetscArrayzero(X+k*ld,2*ld);CHKERRQ(ierr);
        X[k*ld+k]   = Y[0]*norm;
        X[k*ld+k+1] = Y[1]*norm;
      }
      X[(k+1)*ld+k]   = PetscConj(X[k*ld+k]);
      X[(k+1)*ld+k+1] = PetscConj(X[k*ld+k+1]);
#endif
      (*idx)++;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSVectors_GHIEP(DS ds,DSMatType mat,PetscInt *k,PetscReal *rnorm)
{
  PetscInt       i;
  PetscReal      e;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_X:
    case DS_MAT_Y:
      if (k) {
        ierr = DSVectors_GHIEP_Eigen_Some(ds,k,rnorm);CHKERRQ(ierr);
      } else {
        for (i=0; i<ds->n; i++) {
          e = (ds->compact)?*(ds->rmat[DS_MAT_T]+ds->ld+i):PetscRealPart(*(ds->mat[DS_MAT_A]+(i+1)+ds->ld*i));
          if (e == 0.0) { /* real */
            if (ds->state >= DS_STATE_CONDENSED) {
              ierr = PetscArraycpy(ds->mat[mat]+i*ds->ld,ds->mat[DS_MAT_Q]+i*ds->ld,ds->ld);CHKERRQ(ierr);
            } else {
              ierr = PetscArrayzero(ds->mat[mat]+i*ds->ld,ds->ld);CHKERRQ(ierr);
              *(ds->mat[mat]+i+i*ds->ld) = 1.0;
            }
          } else {
            ierr = DSVectors_GHIEP_Eigen_Some(ds,&i,rnorm);CHKERRQ(ierr);
          }
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

/*
  Extract the eigenvalues contained in the block-diagonal of the indefinite problem.
  Only the index range n0..n1 is processed.
*/
PetscErrorCode DSGHIEPComplexEigs(DS ds,PetscInt n0,PetscInt n1,PetscScalar *wr,PetscScalar *wi)
{
  PetscInt     k,ld;
  PetscBLASInt two=2;
  PetscScalar  *A,*B;
  PetscReal    *D,*T;
  PetscReal    b[4],M[4],d1,d2,s1,s2,e;
  PetscReal    scal1,scal2,ep,wr1,wr2,wi1;

  PetscFunctionBegin;
  ld = ds->ld;
  A = ds->mat[DS_MAT_A];
  B = ds->mat[DS_MAT_B];
  D = ds->rmat[DS_MAT_D];
  T = ds->rmat[DS_MAT_T];
  for (k=n0;k<n1;k++) {
    if (k < n1-1) e = (ds->compact)?T[ld+k]:PetscRealPart(A[(k+1)+ld*k]);
    else e = 0.0;
    if (e==0.0) { /* real eigenvalue */
      wr[k] = (ds->compact)?T[k]/D[k]:A[k+k*ld]/B[k+k*ld];
#if !defined(PETSC_USE_COMPLEX)
      wi[k] = 0.0 ;
#endif
    } else { /* diagonal block */
      if (ds->compact) {
        s1 = D[k];
        d1 = T[k];
        s2 = D[k+1];
        d2 = T[k+1];
      } else {
        s1 = PetscRealPart(B[k*ld+k]);
        d1 = PetscRealPart(A[k+k*ld]);
        s2 = PetscRealPart(B[(k+1)*ld+k+1]);
        d2 = PetscRealPart(A[k+1+(k+1)*ld]);
      }
      M[0] = d1; M[1] = e; M[2] = e; M[3]= d2;
      b[0] = s1; b[1] = 0.0; b[2] = 0.0; b[3] = s2;
      ep = LAPACKlamch_("S");
      /* Compute eigenvalues of the block */
      PetscStackCallBLAS("LAPACKlag2",LAPACKlag2_(M,&two,b,&two,&ep,&scal1,&scal2,&wr1,&wr2,&wi1));
      if (scal1<ep) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"Nearly infinite eigenvalue");
      if (wi1==0.0) { /* Real eigenvalues */
        if (scal2<ep) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"Nearly infinite eigenvalue");
        wr[k] = wr1/scal1; wr[k+1] = wr2/scal2;
#if !defined(PETSC_USE_COMPLEX)
        wi[k] = wi[k+1] = 0.0;
#endif
      } else { /* Complex eigenvalues */
#if !defined(PETSC_USE_COMPLEX)
        wr[k]   = wr1/scal1;
        wr[k+1] = wr[k];
        wi[k]   = wi1/scal1;
        wi[k+1] = -wi[k];
#else
        wr[k]   = PetscCMPLX(wr1,wi1)/scal1;
        wr[k+1] = PetscConj(wr[k]);
#endif
      }
      k++;
    }
  }
#if defined(PETSC_USE_COMPLEX)
  if (wi) {
    for (k=n0;k<n1;k++) wi[k] = 0.0;
  }
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode DSSort_GHIEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscErrorCode ierr;
  PetscInt       n,i,*perm;
  PetscReal      *d,*e,*s;

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscValidScalarPointer(wi,3);
#endif
  n = ds->n;
  d = ds->rmat[DS_MAT_T];
  e = d + ds->ld;
  s = ds->rmat[DS_MAT_D];
  ierr = DSAllocateWork_Private(ds,ds->ld,ds->ld,0);CHKERRQ(ierr);
  perm = ds->perm;
  if (!rr) {
    rr = wr;
    ri = wi;
  }
  ierr = DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_TRUE);CHKERRQ(ierr);
  if (!ds->compact) { ierr = DSSwitchFormat_GHIEP(ds,PETSC_TRUE);CHKERRQ(ierr); }
  ierr = PetscArraycpy(ds->work,wr,n);CHKERRQ(ierr);
  for (i=ds->l;i<n;i++) wr[i] = *(ds->work+perm[i]);
#if !defined(PETSC_USE_COMPLEX)
  ierr = PetscArraycpy(ds->work,wi,n);CHKERRQ(ierr);
  for (i=ds->l;i<n;i++) wi[i] = *(ds->work+perm[i]);
#endif
  ierr = PetscArraycpy(ds->rwork,s,n);CHKERRQ(ierr);
  for (i=ds->l;i<n;i++) s[i] = *(ds->rwork+perm[i]);
  ierr = PetscArraycpy(ds->rwork,d,n);CHKERRQ(ierr);
  for (i=ds->l;i<n;i++) d[i] = *(ds->rwork+perm[i]);
  ierr = PetscArraycpy(ds->rwork,e,n-1);CHKERRQ(ierr);
  ierr = PetscArrayzero(e+ds->l,n-1-ds->l);CHKERRQ(ierr);
  for (i=ds->l;i<n-1;i++) {
    if (perm[i]<n-1) e[i] = *(ds->rwork+perm[i]);
  }
  if (!ds->compact) { ierr = DSSwitchFormat_GHIEP(ds,PETSC_FALSE);CHKERRQ(ierr); }
  ierr = DSPermuteColumns_Private(ds,ds->l,n,DS_MAT_Q,perm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Get eigenvectors with inverse iteration.
  The system matrix is in Hessenberg form.
*/
PetscErrorCode DSGHIEPInverseIteration(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscInt       i,off;
  PetscBLASInt   *select,*infoC,ld,n1,mout,info;
  PetscScalar    *A,*B,*H,*X;
  PetscReal      *s,*d,*e;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->n-ds->l,&n1);CHKERRQ(ierr);
  ierr = DSAllocateWork_Private(ds,ld*ld+2*ld,ld,2*ld);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_W);CHKERRQ(ierr);
  A = ds->mat[DS_MAT_A];
  B = ds->mat[DS_MAT_B];
  H = ds->mat[DS_MAT_W];
  s = ds->rmat[DS_MAT_D];
  d = ds->rmat[DS_MAT_T];
  e = d + ld;
  select = ds->iwork;
  infoC = ds->iwork + ld;
  off = ds->l+ds->l*ld;
  if (ds->compact) {
    H[off] = d[ds->l]*s[ds->l];
    H[off+ld] = e[ds->l]*s[ds->l];
    for (i=ds->l+1;i<ds->n-1;i++) {
      H[i+(i-1)*ld] = e[i-1]*s[i];
      H[i+i*ld] = d[i]*s[i];
      H[i+(i+1)*ld] = e[i]*s[i];
    }
    H[ds->n-1+(ds->n-2)*ld] = e[ds->n-2]*s[ds->n-1];
    H[ds->n-1+(ds->n-1)*ld] = d[ds->n-1]*s[ds->n-1];
  } else {
    s[ds->l]  = PetscRealPart(B[off]);
    H[off]    = A[off]*s[ds->l];
    H[off+ld] = A[off+ld]*s[ds->l];
    for (i=ds->l+1;i<ds->n-1;i++) {
      s[i] = PetscRealPart(B[i+i*ld]);
      H[i+(i-1)*ld] = A[i+(i-1)*ld]*s[i];
      H[i+i*ld]     = A[i+i*ld]*s[i];
      H[i+(i+1)*ld] = A[i+(i+1)*ld]*s[i];
    }
    s[ds->n-1] = PetscRealPart(B[ds->n-1+(ds->n-1)*ld]);
    H[ds->n-1+(ds->n-2)*ld] = A[ds->n-1+(ds->n-2)*ld]*s[ds->n-1];
    H[ds->n-1+(ds->n-1)*ld] = A[ds->n-1+(ds->n-1)*ld]*s[ds->n-1];
  }
  ierr = DSAllocateMat_Private(ds,DS_MAT_X);CHKERRQ(ierr);
  X = ds->mat[DS_MAT_X];
  for (i=0;i<n1;i++) select[i] = 1;
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKhsein",LAPACKhsein_("R","N","N",select,&n1,H+off,&ld,wr+ds->l,wi+ds->l,NULL,&ld,X+off,&ld,&n1,&mout,ds->work,NULL,infoC,&info));
#else
  PetscStackCallBLAS("LAPACKhsein",LAPACKhsein_("R","N","N",select,&n1,H+off,&ld,wr+ds->l,NULL,&ld,X+off,&ld,&n1,&mout,ds->work,ds->rwork,NULL,infoC,&info));

  /* Separate real and imaginary part of complex eigenvectors */
  for (j=ds->l;j<ds->n;j++) {
    if (PetscAbsReal(PetscImaginaryPart(wr[j])) > PetscAbsScalar(wr[j])*PETSC_SQRT_MACHINE_EPSILON) {
      for (i=ds->l;i<ds->n;i++) {
        X[i+(j+1)*ds->ld] = PetscImaginaryPart(X[i+j*ds->ld]);
        X[i+j*ds->ld] = PetscRealPart(X[i+j*ds->ld]);
      }
      j++;
    }
  }
#endif
  SlepcCheckLapackInfo("hsein",info);
  ierr = DSGHIEPOrthogEigenv(ds,DS_MAT_X,wr,wi,PETSC_TRUE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Undo 2x2 blocks that have real eigenvalues.
*/
PetscErrorCode DSGHIEPRealBlocks(DS ds)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      e,d1,d2,s1,s2,ss1,ss2,t,dd,ss;
  PetscReal      maxy,ep,scal1,scal2,snorm;
  PetscReal      *T,*D,b[4],M[4],wr1,wr2,wi;
  PetscScalar    *A,*B,Y[4],oneS = 1.0,zeroS = 0.0;
  PetscBLASInt   m,two=2,ld;
  PetscBool      isreal;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->n-ds->l,&m);CHKERRQ(ierr);
  A = ds->mat[DS_MAT_A];
  B = ds->mat[DS_MAT_B];
  T = ds->rmat[DS_MAT_T];
  D = ds->rmat[DS_MAT_D];
  ierr = DSAllocateWork_Private(ds,2*m,0,0);CHKERRQ(ierr);
  for (i=ds->l;i<ds->n-1;i++) {
    e = (ds->compact)?T[ld+i]:PetscRealPart(A[(i+1)+ld*i]);
    if (e != 0.0) { /* 2x2 block */
      if (ds->compact) {
        s1 = D[i];
        d1 = T[i];
        s2 = D[i+1];
        d2 = T[i+1];
      } else {
        s1 = PetscRealPart(B[i*ld+i]);
        d1 = PetscRealPart(A[i*ld+i]);
        s2 = PetscRealPart(B[(i+1)*ld+i+1]);
        d2 = PetscRealPart(A[(i+1)*ld+i+1]);
      }
      isreal = PETSC_FALSE;
      if (s1==s2) { /* apply a Jacobi rotation to compute the eigendecomposition */
        dd = d1-d2;
        if (2*PetscAbsReal(e) <= dd) {
          t = 2*e/dd;
          t = t/(1 + PetscSqrtReal(1+t*t));
        } else {
          t = dd/(2*e);
          ss = (t>=0)?1.0:-1.0;
          t = ss/(PetscAbsReal(t)+PetscSqrtReal(1+t*t));
        }
        Y[0] = 1/PetscSqrtReal(1 + t*t); Y[3] = Y[0]; /* c */
        Y[1] = Y[0]*t; Y[2] = -Y[1]; /* s */
        wr1 = d1+t*e; wr2 = d2-t*e;
        ss1 = s1; ss2 = s2;
        isreal = PETSC_TRUE;
      } else {
        ss1 = 1.0; ss2 = 1.0,
        M[0] = d1; M[1] = e; M[2] = e; M[3]= d2;
        b[0] = s1; b[1] = 0.0; b[2] = 0.0; b[3] = s2;
        ep = LAPACKlamch_("S");

        /* Compute eigenvalues of the block */
        PetscStackCallBLAS("LAPACKlag2",LAPACKlag2_(M,&two,b,&two,&ep,&scal1,&scal2,&wr1,&wr2,&wi));
        if (wi==0.0) { /* Real eigenvalues */
          isreal = PETSC_TRUE;
          if (scal1<ep||scal2<ep) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"Nearly infinite eigenvalue");
          wr1 /= scal1;
          wr2 /= scal2;
          if (PetscAbsReal(s1*d1-wr1)<PetscAbsReal(s2*d2-wr1)) {
            Y[0] = wr1-s2*d2;
            Y[1] = s2*e;
          } else {
            Y[0] = s1*e;
            Y[1] = wr1-s1*d1;
          }
          /* normalize with a signature*/
          maxy = PetscMax(PetscAbsScalar(Y[0]),PetscAbsScalar(Y[1]));
          scal1 = PetscRealPart(Y[0])/maxy;
          scal2 = PetscRealPart(Y[1])/maxy;
          snorm = scal1*scal1*s1 + scal2*scal2*s2;
          if (snorm<0) { ss1 = -1.0; snorm = -snorm; }
          snorm = maxy*PetscSqrtReal(snorm);
          Y[0] = Y[0]/snorm;
          Y[1] = Y[1]/snorm;
          if (PetscAbsReal(s1*d1-wr2)<PetscAbsReal(s2*d2-wr2)) {
            Y[2] = wr2-s2*d2;
            Y[3] = s2*e;
          } else {
            Y[2] = s1*e;
            Y[3] = wr2-s1*d1;
          }
          maxy = PetscMax(PetscAbsScalar(Y[2]),PetscAbsScalar(Y[3]));
          scal1 = PetscRealPart(Y[2])/maxy;
          scal2 = PetscRealPart(Y[3])/maxy;
          snorm = scal1*scal1*s1 + scal2*scal2*s2;
          if (snorm<0) { ss2 = -1.0; snorm = -snorm; }
          snorm = maxy*PetscSqrtReal(snorm); Y[2] = Y[2]/snorm; Y[3] = Y[3]/snorm;
        }
        wr1 *= ss1; wr2 *= ss2;
      }
      if (isreal) {
        if (ds->compact) {
          D[i]    = ss1;
          T[i]    = wr1;
          D[i+1]  = ss2;
          T[i+1]  = wr2;
          T[ld+i] = 0.0;
        } else {
          B[i*ld+i]       = ss1;
          A[i*ld+i]       = wr1;
          B[(i+1)*ld+i+1] = ss2;
          A[(i+1)*ld+i+1] = wr2;
          A[(i+1)+ld*i]   = 0.0;
          A[i+ld*(i+1)]   = 0.0;
        }
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&m,&two,&two,&oneS,ds->mat[DS_MAT_Q]+ds->l+i*ld,&ld,Y,&two,&zeroS,ds->work,&m));
        ierr = PetscArraycpy(ds->mat[DS_MAT_Q]+ds->l+i*ld,ds->work,m);CHKERRQ(ierr);
        ierr = PetscArraycpy(ds->mat[DS_MAT_Q]+ds->l+(i+1)*ld,ds->work+m,m);CHKERRQ(ierr);
      }
      i++;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSolve_GHIEP_QR_II(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscInt       i,off;
  PetscBLASInt   n1,ld,one,info,lwork;
  PetscScalar    *H,*A,*B,*Q;
  PetscReal      *d,*e,*s;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscValidPointer(wi,3);
#endif
  one = 1;
  ierr = PetscBLASIntCast(ds->n-ds->l,&n1);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  off = ds->l + ds->l*ld;
  A = ds->mat[DS_MAT_A];
  B = ds->mat[DS_MAT_B];
  Q = ds->mat[DS_MAT_Q];
  d = ds->rmat[DS_MAT_T];
  e = ds->rmat[DS_MAT_T] + ld;
  s = ds->rmat[DS_MAT_D];
  ierr = DSAllocateWork_Private(ds,ld*ld,2*ld,ld*2);CHKERRQ(ierr);
  lwork = ld*ld;

  /* Quick return if possible */
  if (n1 == 1) {
    for (i=0;i<=ds->l;i++) Q[i+i*ld] = 1.0;
    ierr = DSGHIEPComplexEigs(ds,0,ds->l,wr,wi);CHKERRQ(ierr);
    if (!ds->compact) {
      d[ds->l] = PetscRealPart(A[off]);
      s[ds->l] = PetscRealPart(B[off]);
    }
    wr[ds->l] = d[ds->l]/s[ds->l];
    if (wi) wi[ds->l] = 0.0;
    PetscFunctionReturn(0);
  }
  /* Reduce to pseudotriadiagonal form */
  ierr = DSIntermediate_GHIEP(ds);CHKERRQ(ierr);

  /* Compute Eigenvalues (QR)*/
  ierr = DSAllocateMat_Private(ds,DS_MAT_W);CHKERRQ(ierr);
  H = ds->mat[DS_MAT_W];
  if (ds->compact) {
    H[off]    = d[ds->l]*s[ds->l];
    H[off+ld] = e[ds->l]*s[ds->l];
    for (i=ds->l+1;i<ds->n-1;i++) {
      H[i+(i-1)*ld] = e[i-1]*s[i];
      H[i+i*ld]     = d[i]*s[i];
      H[i+(i+1)*ld] = e[i]*s[i];
    }
    H[ds->n-1+(ds->n-2)*ld] = e[ds->n-2]*s[ds->n-1];
    H[ds->n-1+(ds->n-1)*ld] = d[ds->n-1]*s[ds->n-1];
  } else {
    s[ds->l]  = PetscRealPart(B[off]);
    H[off]    = A[off]*s[ds->l];
    H[off+ld] = A[off+ld]*s[ds->l];
    for (i=ds->l+1;i<ds->n-1;i++) {
      s[i] = PetscRealPart(B[i+i*ld]);
      H[i+(i-1)*ld] = A[i+(i-1)*ld]*s[i];
      H[i+i*ld]     = A[i+i*ld]*s[i];
      H[i+(i+1)*ld] = A[i+(i+1)*ld]*s[i];
    }
    s[ds->n-1] = PetscRealPart(B[ds->n-1+(ds->n-1)*ld]);
    H[ds->n-1+(ds->n-2)*ld] = A[ds->n-1+(ds->n-2)*ld]*s[ds->n-1];
    H[ds->n-1+(ds->n-1)*ld] = A[ds->n-1+(ds->n-1)*ld]*s[ds->n-1];
  }

#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKhseqr",LAPACKhseqr_("E","N",&n1,&one,&n1,H+off,&ld,wr+ds->l,wi+ds->l,NULL,&ld,ds->work,&lwork,&info));
#else
  PetscStackCallBLAS("LAPACKhseqr",LAPACKhseqr_("E","N",&n1,&one,&n1,H+off,&ld,wr+ds->l,NULL,&ld,ds->work,&lwork,&info));
  for (i=ds->l;i<ds->n;i++) if (PetscAbsReal(PetscImaginaryPart(wr[i]))<10*PETSC_MACHINE_EPSILON) wr[i] = PetscRealPart(wr[i]);
  /* Sort to have consecutive conjugate pairs */
  for (i=ds->l;i<ds->n;i++) {
      j=i+1;
      while (j<ds->n && (PetscAbsScalar(wr[i]-PetscConj(wr[j]))>PetscAbsScalar(wr[i])*PETSC_SQRT_MACHINE_EPSILON)) j++;
      if (j==ds->n) {
        if (PetscAbsReal(PetscImaginaryPart(wr[i]))<PetscAbsScalar(wr[i])*PETSC_SQRT_MACHINE_EPSILON) wr[i]=PetscRealPart(wr[i]);
        else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"Found complex without conjugate pair");
      } else { /* complex eigenvalue */
        wr[j] = wr[i+1];
        if (PetscImaginaryPart(wr[i])<0) wr[i] = PetscConj(wr[i]);
        wr[i+1] = PetscConj(wr[i]);
        i++;
      }
  }
#endif
  SlepcCheckLapackInfo("hseqr",info);
  /* Compute Eigenvectors with Inverse Iteration */
  ierr = DSGHIEPInverseIteration(ds,wr,wi);CHKERRQ(ierr);

  /* Recover eigenvalues from diagonal */
  ierr = DSGHIEPComplexEigs(ds,0,ds->l,wr,wi);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  if (wi) {
    for (i=ds->l;i<ds->n;i++) wi[i] = 0.0;
  }
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode DSSolve_GHIEP_QR(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscInt       i,j,off,nwu=0,n,lw,lwr,nwru=0;
  PetscBLASInt   n_,ld,info,lwork,ilo,ihi;
  PetscScalar    *H,*A,*B,*Q,*X;
  PetscReal      *d,*s,*scale,nrm,*rcde,*rcdv;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       k;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscValidPointer(wi,3);
#endif
  n = ds->n-ds->l;
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  off = ds->l + ds->l*ld;
  A = ds->mat[DS_MAT_A];
  B = ds->mat[DS_MAT_B];
  Q = ds->mat[DS_MAT_Q];
  d = ds->rmat[DS_MAT_T];
  s = ds->rmat[DS_MAT_D];
  lw = 14*ld+ld*ld;
  lwr = 7*ld;
  ierr = DSAllocateWork_Private(ds,lw,lwr,0);CHKERRQ(ierr);
  scale = ds->rwork+nwru;
  nwru += ld;
  rcde = ds->rwork+nwru;
  nwru += ld;
  rcdv = ds->rwork+nwru;
  nwru += ld;
  /* Quick return if possible */
  if (n_ == 1) {
    for (i=0;i<=ds->l;i++) Q[i+i*ld] = 1.0;
    ierr = DSGHIEPComplexEigs(ds,0,ds->l,wr,wi);CHKERRQ(ierr);
    if (!ds->compact) {
      d[ds->l] = PetscRealPart(A[off]);
      s[ds->l] = PetscRealPart(B[off]);
    }
    wr[ds->l] = d[ds->l]/s[ds->l];
    if (wi) wi[ds->l] = 0.0;
    PetscFunctionReturn(0);
  }

  /* Form pseudo-symmetric matrix */
  H =  ds->work+nwu;
  nwu += n*n;
  ierr = PetscArrayzero(H,n*n);CHKERRQ(ierr);
  if (ds->compact) {
    for (i=0;i<n-1;i++) {
      H[i+i*n]     = s[ds->l+i]*d[ds->l+i];
      H[i+1+i*n]   = s[ds->l+i+1]*d[ld+ds->l+i];
      H[i+(i+1)*n] = s[ds->l+i]*d[ld+ds->l+i];
    }
    H[n-1+(n-1)*n] = s[ds->l+n-1]*d[ds->l+n-1];
    for (i=0;i<ds->k-ds->l;i++) {
      H[ds->k-ds->l+i*n] = s[ds->k]*d[2*ld+ds->l+i];
      H[i+(ds->k-ds->l)*n] = s[i+ds->l]*d[2*ld+ds->l+i];
    }
  } else {
    for (j=0;j<n;j++) {
      for (i=0;i<n;i++) H[i+j*n] = B[off+i+i*ld]*A[off+i+j*ld];
    }
  }

  /* Compute eigenpairs */
  ierr = PetscBLASIntCast(lw-nwu,&lwork);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_X);CHKERRQ(ierr);
  X = ds->mat[DS_MAT_X];
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgeevx",LAPACKgeevx_("B","N","V","N",&n_,H,&n_,wr+ds->l,wi+ds->l,NULL,&ld,X+off,&ld,&ilo,&ihi,scale,&nrm,rcde,rcdv,ds->work+nwu,&lwork,NULL,&info));
#else
  PetscStackCallBLAS("LAPACKgeevx",LAPACKgeevx_("B","N","V","N",&n_,H,&n_,wr+ds->l,NULL,&ld,X+off,&ld,&ilo,&ihi,scale,&nrm,rcde,rcdv,ds->work+nwu,&lwork,ds->rwork+nwru,&info));

  /* Sort to have consecutive conjugate pairs
     Separate real and imaginary part of complex eigenvectors*/
  for (i=ds->l;i<ds->n;i++) {
    j=i+1;
    while (j<ds->n && (PetscAbsScalar(wr[i]-PetscConj(wr[j]))>PetscAbsScalar(wr[i])*PETSC_SQRT_MACHINE_EPSILON)) j++;
    if (j==ds->n) {
      if (PetscAbsReal(PetscImaginaryPart(wr[i]))<PetscAbsScalar(wr[i])*PETSC_SQRT_MACHINE_EPSILON) {
        wr[i]=PetscRealPart(wr[i]); /* real eigenvalue */
        for (k=ds->l;k<ds->n;k++) {
          X[k+i*ds->ld] = PetscRealPart(X[k+i*ds->ld]);
        }
      } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"Found complex without conjugate pair");
    } else { /* complex eigenvalue */
      if (j!=i+1) {
        wr[j] = wr[i+1];
        ierr = PetscArraycpy(X+j*ds->ld,X+(i+1)*ds->ld,ds->ld);CHKERRQ(ierr);
      }
      if (PetscImaginaryPart(wr[i])<0) {
        wr[i] = PetscConj(wr[i]);
        for (k=ds->l;k<ds->n;k++) {
          X[k+(i+1)*ds->ld] = -PetscImaginaryPart(X[k+i*ds->ld]);
          X[k+i*ds->ld] = PetscRealPart(X[k+i*ds->ld]);
        }
      } else {
        for (k=ds->l;k<ds->n;k++) {
          X[k+(i+1)*ds->ld] = PetscImaginaryPart(X[k+i*ds->ld]);
          X[k+i*ds->ld] = PetscRealPart(X[k+i*ds->ld]);
        }
      }
      wr[i+1] = PetscConj(wr[i]);
      i++;
    }
  }
#endif
  SlepcCheckLapackInfo("geevx",info);

  /* Compute real s-orthonormal basis */
  ierr = DSGHIEPOrthogEigenv(ds,DS_MAT_X,wr,wi,PETSC_FALSE);CHKERRQ(ierr);

  /* Recover eigenvalues from diagonal */
  ierr = DSGHIEPComplexEigs(ds,0,ds->l,wr,wi);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  if (wi) {
    for (i=ds->l;i<ds->n;i++) wi[i] = 0.0;
  }
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode DSSynchronize_GHIEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscErrorCode ierr;
  PetscInt       ld=ds->ld,l=ds->l,k=0,kr=0;
  PetscMPIInt    n,rank,off=0,size,ldn,ld3,ld_;

  PetscFunctionBegin;
  if (ds->compact) kr = 4*ld;
  else k = 2*(ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += (ds->n-l)*ld;
  if (eigr) k += (ds->n-l);
  if (eigi) k += (ds->n-l);
  ierr = DSAllocateWork_Private(ds,k+kr,0,0);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(k*sizeof(PetscScalar)+kr*sizeof(PetscReal),&size);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ds->n-l,&n);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ld*(ds->n-l),&ldn);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ld*3,&ld3);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ld,&ld_);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank);CHKERRQ(ierr);
  if (!rank) {
    if (ds->compact) {
      ierr = MPI_Pack(ds->rmat[DS_MAT_T],ld3,MPIU_REAL,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
      ierr = MPI_Pack(ds->rmat[DS_MAT_D],ld_,MPIU_REAL,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    } else {
      ierr = MPI_Pack(ds->mat[DS_MAT_A]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
      ierr = MPI_Pack(ds->mat[DS_MAT_B]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
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
    if (ds->compact) {
      ierr = MPI_Unpack(ds->work,size,&off,ds->rmat[DS_MAT_T],ld3,MPIU_REAL,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
      ierr = MPI_Unpack(ds->work,size,&off,ds->rmat[DS_MAT_D],ld_,MPIU_REAL,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    } else {
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_A]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_B]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
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

PetscErrorCode DSHermitian_GHIEP(DS ds,DSMatType m,PetscBool *flg)
{
  PetscFunctionBegin;
  if (m==DS_MAT_A || m==DS_MAT_B) *flg = PETSC_TRUE;
  else *flg = PETSC_FALSE;
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode DSCreate_GHIEP(DS ds)
{
  PetscFunctionBegin;
  ds->ops->allocate      = DSAllocate_GHIEP;
  ds->ops->view          = DSView_GHIEP;
  ds->ops->vectors       = DSVectors_GHIEP;
  ds->ops->solve[0]      = DSSolve_GHIEP_QR_II;
  ds->ops->solve[1]      = DSSolve_GHIEP_HZ;
  ds->ops->solve[2]      = DSSolve_GHIEP_QR;
  ds->ops->sort          = DSSort_GHIEP;
  ds->ops->synchronize   = DSSynchronize_GHIEP;
  ds->ops->hermitian     = DSHermitian_GHIEP;
  PetscFunctionReturn(0);
}

