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

/*
  1) Patterns of A and B
      DS_STATE_RAW:       DS_STATE_INTERM/CONDENSED
       0       n-1              0       n-1
      -------------            -------------
    0 |* * * * * *|          0 |* * * * * *|
      |* * * * * *|            |  * * * * *|
      |* * * * * *|            |    * * * *|
      |* * * * * *|            |    * * * *|
      |* * * * * *|            |        * *|
  n-1 |* * * * * *|        n-1 |          *|
      -------------            -------------

  2) Moreover, P and Q are assumed to be the identity in DS_STATE_INTERMEDIATE.
*/


static PetscErrorCode CleanDenseSchur(PetscInt n,PetscInt k,PetscScalar *S,PetscInt ldS,PetscScalar *T,PetscInt ldT,PetscScalar *X,PetscInt ldX,PetscScalar *Y,PetscInt ldY);

PetscErrorCode DSAllocate_GNHEP(DS ds,PetscInt ld)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DSAllocateMat_Private(ds,DS_MAT_A);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_B);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_Z);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_Q);CHKERRQ(ierr);
  ierr = PetscFree(ds->perm);CHKERRQ(ierr);
  ierr = PetscMalloc1(ld,&ds->perm);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)ds,ld*sizeof(PetscInt));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DSView_GNHEP(DS ds,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  PetscViewerFormat format;

  PetscFunctionBegin;
  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
  ierr = DSViewMat(ds,viewer,DS_MAT_A);CHKERRQ(ierr);
  ierr = DSViewMat(ds,viewer,DS_MAT_B);CHKERRQ(ierr);
  if (ds->state>DS_STATE_INTERMEDIATE) {
    ierr = DSViewMat(ds,viewer,DS_MAT_Z);CHKERRQ(ierr);
    ierr = DSViewMat(ds,viewer,DS_MAT_Q);CHKERRQ(ierr);
  }
  if (ds->mat[DS_MAT_X]) { ierr = DSViewMat(ds,viewer,DS_MAT_X);CHKERRQ(ierr); }
  if (ds->mat[DS_MAT_Y]) { ierr = DSViewMat(ds,viewer,DS_MAT_Y);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSVectors_GNHEP_Eigen_Some(DS ds,PetscInt *k,PetscReal *rnorm,PetscBool left)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBLASInt   n,ld,mout,info,*select,mm,inc = 1;
  PetscScalar    *X,*Y,*Z,*A = ds->mat[DS_MAT_A],*B = ds->mat[DS_MAT_B],tmp,fone=1.0,fzero=0.0;
  PetscReal      norm;
  PetscBool      iscomplex = PETSC_FALSE;
  const char     *side;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  if (left) {
    X = NULL;
    Y = &ds->mat[DS_MAT_Y][ld*(*k)];
    side = "L";
  } else {
    X = &ds->mat[DS_MAT_X][ld*(*k)];
    Y = NULL;
    side = "R";
  }
  Z = left? Y: X;
  ierr = DSAllocateWork_Private(ds,0,0,ld);CHKERRQ(ierr);
  select = ds->iwork;
  for (i=0;i<n;i++) select[i] = (PetscBLASInt)PETSC_FALSE;
  if (ds->state <= DS_STATE_INTERMEDIATE) {
    ierr = DSSetIdentity(ds,DS_MAT_Q);CHKERRQ(ierr);
    ierr = DSSetIdentity(ds,DS_MAT_Z);CHKERRQ(ierr);
  }
  ierr = CleanDenseSchur(n,0,A,ld,B,ld,ds->mat[DS_MAT_Q],ld,ds->mat[DS_MAT_Z],ld);CHKERRQ(ierr);
  if (ds->state < DS_STATE_CONDENSED) { ierr = DSSetState(ds,DS_STATE_CONDENSED);CHKERRQ(ierr); }

  /* compute k-th eigenvector */
  select[*k] = (PetscBLASInt)PETSC_TRUE;
#if defined(PETSC_USE_COMPLEX)
  mm = 1;
  ierr = DSAllocateWork_Private(ds,2*ld,2*ld,0);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKtgevc",LAPACKtgevc_(side,"S",select,&n,A,&ld,B,&ld,Y,&ld,X,&ld,&mm,&mout,ds->work,ds->rwork,&info));
#else
  if ((*k)<n-1 && (A[ld*(*k)+(*k)+1] != 0.0 || B[ld*(*k)+(*k)+1] != 0.0)) iscomplex = PETSC_TRUE;
  mm = iscomplex? 2: 1;
  if (iscomplex) select[(*k)+1] = (PetscBLASInt)PETSC_TRUE;
  ierr = DSAllocateWork_Private(ds,6*ld,0,0);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKtgevc",LAPACKtgevc_(side,"S",select,&n,A,&ld,B,&ld,Y,&ld,X,&ld,&mm,&mout,ds->work,&info));
#endif
  SlepcCheckLapackInfo("tgevc",info);
  if (!select[*k] || mout != mm) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong arguments in call to Lapack xTGEVC");

  /* accumulate and normalize eigenvectors */
  ierr = PetscArraycpy(ds->work,Z,mm*ld);CHKERRQ(ierr);
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&mm,&n,&fone,ds->mat[left?DS_MAT_Z:DS_MAT_Q],&ld,ds->work,&ld,&fzero,Z,&ld));
  norm = BLASnrm2_(&n,Z,&inc);
#if !defined(PETSC_USE_COMPLEX)
  if (iscomplex) {
    tmp = BLASnrm2_(&n,Z+ld,&inc);
    norm = SlepcAbsEigenvalue(norm,tmp);
  }
#endif
  tmp = 1.0 / norm;
  PetscStackCallBLAS("BLASscal",BLASscal_(&n,&tmp,Z,&inc));
#if !defined(PETSC_USE_COMPLEX)
  if (iscomplex) PetscStackCallBLAS("BLASscal",BLASscal_(&n,&tmp,Z+ld,&inc));
#endif

  /* set output arguments */
  if (iscomplex) (*k)++;
  if (rnorm) {
    if (iscomplex) *rnorm = SlepcAbsEigenvalue(Z[n-1],Z[n-1+ld]);
    else *rnorm = PetscAbsScalar(Z[n-1]);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSVectors_GNHEP_Eigen_All(DS ds,PetscBool left)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBLASInt   n,ld,mout,info,inc = 1;
  PetscBool      iscomplex;
  PetscScalar    *X,*Y,*Z,*A = ds->mat[DS_MAT_A],*B = ds->mat[DS_MAT_B],tmp;
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
  if (ds->state <= DS_STATE_INTERMEDIATE) {
    ierr = DSSetIdentity(ds,DS_MAT_Q);CHKERRQ(ierr);
    ierr = DSSetIdentity(ds,DS_MAT_Z);CHKERRQ(ierr);
  }
  ierr = CleanDenseSchur(n,0,A,ld,B,ld,ds->mat[DS_MAT_Q],ld,ds->mat[DS_MAT_Z],ld);CHKERRQ(ierr);
  if (ds->state>=DS_STATE_CONDENSED) {
    /* DSSolve() has been called, backtransform with matrix Q */
    back = "B";
    ierr = PetscArraycpy(left?Y:X,ds->mat[left?DS_MAT_Z:DS_MAT_Q],ld*ld);CHKERRQ(ierr);
  } else {
    back = "A";
    ierr = DSSetState(ds,DS_STATE_CONDENSED);CHKERRQ(ierr);
  }
#if defined(PETSC_USE_COMPLEX)
  ierr = DSAllocateWork_Private(ds,2*ld,2*ld,0);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKtgevc",LAPACKtgevc_(side,back,NULL,&n,A,&ld,B,&ld,Y,&ld,X,&ld,&n,&mout,ds->work,ds->rwork,&info));
#else
  ierr = DSAllocateWork_Private(ds,6*ld,0,0);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKtgevc",LAPACKtgevc_(side,back,NULL,&n,A,&ld,B,&ld,Y,&ld,X,&ld,&n,&mout,ds->work,&info));
#endif
  SlepcCheckLapackInfo("tgevc",info);

  /* normalize eigenvectors */
  for (i=0;i<n;i++) {
    iscomplex = (i<n-1 && (A[i+1+i*ld]!=0.0 || B[i+1+i*ld]!=0.0))? PETSC_TRUE: PETSC_FALSE;
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

PetscErrorCode DSVectors_GNHEP(DS ds,DSMatType mat,PetscInt *k,PetscReal *rnorm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  switch (mat) {
    case DS_MAT_X:
    case DS_MAT_Y:
      if (k) {
        ierr = DSVectors_GNHEP_Eigen_Some(ds,k,rnorm,mat == DS_MAT_Y?PETSC_TRUE:PETSC_FALSE);CHKERRQ(ierr);
      } else {
        ierr = DSVectors_GNHEP_Eigen_All(ds,mat == DS_MAT_Y?PETSC_TRUE:PETSC_FALSE);CHKERRQ(ierr);
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSSort_GNHEP_Arbitrary(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBLASInt   info,n,ld,mout,lwork,liwork,*iwork,*selection,zero_=0,true_=1;
  PetscScalar    *S = ds->mat[DS_MAT_A],*T = ds->mat[DS_MAT_B],*Q = ds->mat[DS_MAT_Q],*Z = ds->mat[DS_MAT_Z],*work,*beta;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(0);
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  lwork = 4*n+16;
#else
  lwork = 1;
#endif
  liwork = 1;
  ierr = DSAllocateWork_Private(ds,lwork+2*n,0,liwork+n);CHKERRQ(ierr);
  beta      = ds->work;
  work      = ds->work + n;
  lwork     = ds->lwork - n;
  selection = ds->iwork;
  iwork     = ds->iwork + n;
  liwork    = ds->liwork - n;
  /* Compute the selected eigenvalue to be in the leading position */
  ierr = DSSortEigenvalues_Private(ds,rr,ri,ds->perm,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscArrayzero(selection,n);CHKERRQ(ierr);
  for (i=0; i<*k; i++) selection[ds->perm[i]] = 1;
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKtgsen",LAPACKtgsen_(&zero_,&true_,&true_,selection,&n,S,&ld,T,&ld,wr,wi,beta,Z,&ld,Q,&ld,&mout,NULL,NULL,NULL,work,&lwork,iwork,&liwork,&info));
#else
  PetscStackCallBLAS("LAPACKtgsen",LAPACKtgsen_(&zero_,&true_,&true_,selection,&n,S,&ld,T,&ld,wr,beta,Z,&ld,Q,&ld,&mout,NULL,NULL,NULL,work,&lwork,iwork,&liwork,&info));
#endif
  SlepcCheckLapackInfo("tgsen",info);
  *k = mout;
  for (i=0;i<n;i++) {
    if (beta[i]==0.0) wr[i] = (PetscRealPart(wr[i])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wr[i] /= beta[i];
#if !defined(PETSC_USE_COMPLEX)
    if (beta[i]==0.0) wi[i] = (wi[i]>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wi[i] /= beta[i];
#endif
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSSort_GNHEP_Total(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscScalar    re;
  PetscInt       i,j,pos,result;
  PetscBLASInt   ifst,ilst,info,n,ld,one=1;
  PetscScalar    *S = ds->mat[DS_MAT_A],*T = ds->mat[DS_MAT_B],*Z = ds->mat[DS_MAT_Z],*Q = ds->mat[DS_MAT_Q];
#if !defined(PETSC_USE_COMPLEX)
  PetscBLASInt   lwork;
  PetscScalar    *work,a,safmin,scale1,scale2,im;
#endif

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(0);
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  lwork = -1;
  PetscStackCallBLAS("LAPACKtgexc",LAPACKtgexc_(&one,&one,&ld,NULL,&ld,NULL,&ld,NULL,&ld,NULL,&ld,&one,&one,&a,&lwork,&info));
  SlepcCheckLapackInfo("tgexc",info);
  safmin = LAPACKlamch_("S");
  ierr = PetscBLASIntCast((PetscInt)a,&lwork);CHKERRQ(ierr);
  ierr = DSAllocateWork_Private(ds,lwork,0,0);CHKERRQ(ierr);
  work = ds->work;
#endif
  /* selection sort */
  for (i=ds->l;i<n-1;i++) {
    re = wr[i];
#if !defined(PETSC_USE_COMPLEX)
    im = wi[i];
#endif
    pos = 0;
    j = i+1; /* j points to the next eigenvalue */
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
      PetscStackCallBLAS("LAPACKtgexc",LAPACKtgexc_(&one,&one,&n,S,&ld,T,&ld,Z,&ld,Q,&ld,&ifst,&ilst,work,&lwork,&info));
#else
      PetscStackCallBLAS("LAPACKtgexc",LAPACKtgexc_(&one,&one,&n,S,&ld,T,&ld,Z,&ld,Q,&ld,&ifst,&ilst,&info));
#endif
      SlepcCheckLapackInfo("tgexc",info);
      /* recover original eigenvalues from T and S matrices */
      for (j=i;j<n;j++) {
#if !defined(PETSC_USE_COMPLEX)
        if (j<n-1 && S[j*ld+j+1] != 0.0) {
          /* complex conjugate eigenvalue */
          PetscStackCallBLAS("LAPACKlag2",LAPACKlag2_(S+j*ld+j,&ld,T+j*ld+j,&ld,&safmin,&scale1,&scale2,&re,&a,&im));
          wr[j] = re / scale1;
          wi[j] = im / scale1;
          wr[j+1] = a / scale2;
          wi[j+1] = -wi[j];
          j++;
        } else
#endif
        {
          if (T[j*ld+j] == 0.0) wr[j] = (PetscRealPart(S[j*ld+j])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
          else wr[j] = S[j*ld+j] / T[j*ld+j];
#if !defined(PETSC_USE_COMPLEX)
          wi[j] = 0.0;
#endif
        }
      }
    }
#if !defined(PETSC_USE_COMPLEX)
    if (wi[i] != 0.0) i++;
#endif
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSort_GNHEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *k)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!rr || wr == rr) {
    ierr = DSSort_GNHEP_Total(ds,wr,wi);CHKERRQ(ierr);
  } else {
    ierr = DSSort_GNHEP_Arbitrary(ds,wr,wi,rr,ri,k);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   Write zeros from the column k to n in the lower triangular part of the
   matrices S and T, and inside 2-by-2 diagonal blocks of T in order to
   make (S,T) a valid Schur decompositon.
*/
static PetscErrorCode CleanDenseSchur(PetscInt n,PetscInt k,PetscScalar *S,PetscInt ldS,PetscScalar *T,PetscInt ldT,PetscScalar *X,PetscInt ldX,PetscScalar *Y,PetscInt ldY)
{
  PetscInt       i;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
  PetscScalar    s;
#else
  PetscErrorCode ierr;
  PetscBLASInt   ldS_,ldT_,n_i,n_i_2,one=1,n_,i_2,i_;
  PetscScalar    b11,b22,sr,cr,sl,cl;
#endif

  PetscFunctionBegin;
#if defined(PETSC_USE_COMPLEX)
  for (i=k; i<n; i++) {
    /* Some functions need the diagonal elements in T be real */
    if (T && PetscImaginaryPart(T[ldT*i+i]) != 0.0) {
      s = PetscConj(T[ldT*i+i])/PetscAbsScalar(T[ldT*i+i]);
      for (j=0;j<=i;j++) {
        T[ldT*i+j] *= s;
        S[ldS*i+j] *= s;
      }
      T[ldT*i+i] = PetscRealPart(T[ldT*i+i]);
      if (X) for (j=0;j<n;j++) X[ldX*i+j] *= s;
    }
    j = i+1;
    if (j<n) {
      S[ldS*i+j] = 0.0;
      if (T) T[ldT*i+j] = 0.0;
    }
  }
#else
  ierr = PetscBLASIntCast(ldS,&ldS_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldT,&ldT_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  for (i=k;i<n-1;i++) {
    if (S[ldS*i+i+1] != 0.0) {
      /* Check if T(i+1,i) and T(i,i+1) are zero */
      if (T[ldT*(i+1)+i] != 0.0 || T[ldT*i+i+1] != 0.0) {
        /* Check if T(i+1,i) and T(i,i+1) are negligible */
        if (PetscAbs(T[ldT*(i+1)+i])+PetscAbs(T[ldT*i+i+1]) < (PetscAbs(T[ldT*i+i])+PetscAbs(T[ldT*(i+1)+i+1]))*PETSC_MACHINE_EPSILON) {
          T[ldT*i+i+1] = 0.0;
          T[ldT*(i+1)+i] = 0.0;
        } else {
          /* If one of T(i+1,i) or T(i,i+1) is negligible, we make zero the other element */
          if (PetscAbs(T[ldT*i+i+1]) < (PetscAbs(T[ldT*i+i])+PetscAbs(T[ldT*(i+1)+i+1])+PetscAbs(T[ldT*(i+1)+i]))*PETSC_MACHINE_EPSILON) {
            PetscStackCallBLAS("LAPACKlasv2",LAPACKlasv2_(&T[ldT*i+i],&T[ldT*(i+1)+i],&T[ldT*(i+1)+i+1],&b22,&b11,&sl,&cl,&sr,&cr));
          } else if (PetscAbs(T[ldT*(i+1)+i]) < (PetscAbs(T[ldT*i+i])+PetscAbs(T[ldT*(i+1)+i+1])+PetscAbs(T[ldT*i+i+1]))*PETSC_MACHINE_EPSILON) {
            PetscStackCallBLAS("LAPACKlasv2",LAPACKlasv2_(&T[ldT*i+i],&T[ldT*i+i+1],&T[ldT*(i+1)+i+1],&b22,&b11,&sr,&cr,&sl,&cl));
          } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unsupported format. Call DSSolve before this function");
          ierr = PetscBLASIntCast(n-i,&n_i);CHKERRQ(ierr);
          n_i_2 = n_i - 2;
          ierr = PetscBLASIntCast(i+2,&i_2);CHKERRQ(ierr);
          ierr = PetscBLASIntCast(i,&i_);CHKERRQ(ierr);
          if (b11 < 0.0) {
            cr = -cr; sr = -sr;
            b11 = -b11; b22 = -b22;
          }
          PetscStackCallBLAS("BLASrot",BLASrot_(&n_i,&S[ldS*i+i],&ldS_,&S[ldS*i+i+1],&ldS_,&cl,&sl));
          PetscStackCallBLAS("BLASrot",BLASrot_(&i_2,&S[ldS*i],&one,&S[ldS*(i+1)],&one,&cr,&sr));
          PetscStackCallBLAS("BLASrot",BLASrot_(&n_i_2,&T[ldT*(i+2)+i],&ldT_,&T[ldT*(i+2)+i+1],&ldT_,&cl,&sl));
          PetscStackCallBLAS("BLASrot",BLASrot_(&i_,&T[ldT*i],&one,&T[ldT*(i+1)],&one,&cr,&sr));
          if (X) PetscStackCallBLAS("BLASrot",BLASrot_(&n_,&X[ldX*i],&one,&X[ldX*(i+1)],&one,&cr,&sr));
          if (Y) PetscStackCallBLAS("BLASrot",BLASrot_(&n_,&Y[ldY*i],&one,&Y[ldY*(i+1)],&one,&cl,&sl));
          T[ldT*i+i] = b11; T[ldT*i+i+1] = 0.0;
          T[ldT*(i+1)+i] = 0.0; T[ldT*(i+1)+i+1] = b22;
        }
      }
      i++;
    }
  }
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode DSSolve_GNHEP(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  PetscScalar    *work,*beta,a;
  PetscInt       i;
  PetscBLASInt   lwork,info,n,ld,iaux;
  PetscScalar    *A = ds->mat[DS_MAT_A],*B = ds->mat[DS_MAT_B],*Z = ds->mat[DS_MAT_Z],*Q = ds->mat[DS_MAT_Q];

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  PetscValidScalarPointer(wi,3);
#endif
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld,&ld);CHKERRQ(ierr);
  lwork = -1;
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgges",LAPACKgges_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,wi,NULL,Z,&ld,Q,&ld,&a,&lwork,NULL,&info));
  ierr = PetscBLASIntCast((PetscInt)a,&lwork);CHKERRQ(ierr);
  ierr = DSAllocateWork_Private(ds,lwork+ld,0,0);CHKERRQ(ierr);
  beta = ds->work;
  work = beta+ds->n;
  ierr = PetscBLASIntCast(ds->lwork-ds->n,&lwork);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgges",LAPACKgges_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,wi,beta,Z,&ld,Q,&ld,work,&lwork,NULL,&info));
#else
  PetscStackCallBLAS("LAPACKgges",LAPACKgges_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,NULL,Z,&ld,Q,&ld,&a,&lwork,NULL,NULL,&info));
  ierr = PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork);CHKERRQ(ierr);
  ierr = DSAllocateWork_Private(ds,lwork+ld,8*ld,0);CHKERRQ(ierr);
  beta = ds->work;
  work = beta+ds->n;
  ierr = PetscBLASIntCast(ds->lwork-ds->n,&lwork);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgges",LAPACKgges_("V","V","N",NULL,&n,A,&ld,B,&ld,&iaux,wr,beta,Z,&ld,Q,&ld,work,&lwork,ds->rwork,NULL,&info));
#endif
  SlepcCheckLapackInfo("gges",info);
  for (i=0;i<n;i++) {
    if (beta[i]==0.0) wr[i] = (PetscRealPart(wr[i])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wr[i] /= beta[i];
#if !defined(PETSC_USE_COMPLEX)
    if (beta[i]==0.0) wi[i] = (wi[i]>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wi[i] /= beta[i];
#else
    if (wi) wi[i] = 0.0;
#endif
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSynchronize_GNHEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscErrorCode ierr;
  PetscInt       ld=ds->ld,l=ds->l,k;
  PetscMPIInt    n,rank,off=0,size,ldn;

  PetscFunctionBegin;
  k = 2*(ds->n-l)*ld;
  if (ds->state>DS_STATE_RAW) k += 2*(ds->n-l)*ld;
  if (eigr) k += (ds->n-l);
  if (eigi) k += (ds->n-l);
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
      ierr = MPI_Pack(ds->mat[DS_MAT_Z]+l*ld,ldn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
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
    ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_B]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    if (ds->state>DS_STATE_RAW) {
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_Q]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_Z]+l*ld,ldn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
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

SLEPC_EXTERN PetscErrorCode DSCreate_GNHEP(DS ds)
{
  PetscFunctionBegin;
  ds->ops->allocate      = DSAllocate_GNHEP;
  ds->ops->view          = DSView_GNHEP;
  ds->ops->vectors       = DSVectors_GNHEP;
  ds->ops->solve[0]      = DSSolve_GNHEP;
  ds->ops->sort          = DSSort_GNHEP;
  ds->ops->synchronize   = DSSynchronize_GNHEP;
  PetscFunctionReturn(0);
}

