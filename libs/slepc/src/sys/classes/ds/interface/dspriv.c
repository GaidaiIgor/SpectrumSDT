/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Private DS routines
*/

#include <slepc/private/dsimpl.h>      /*I "slepcds.h" I*/
#include <slepcblaslapack.h>

PetscErrorCode DSAllocateMatrix_Private(DS ds,DSMatType m,PetscBool isreal)
{
  size_t         sz;
  PetscInt       n,d,nelem;
  PetscBool      ispep;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)ds,DSPEP,&ispep);CHKERRQ(ierr);
  if (ispep) {
    ierr = DSPEPGetDegree(ds,&d);CHKERRQ(ierr);
  }
  if (ispep && (m==DS_MAT_A || m==DS_MAT_B || m==DS_MAT_W || m==DS_MAT_U || m==DS_MAT_X || m==DS_MAT_Y)) n = d*ds->ld;
  else n = ds->ld;
  switch (m) {
    case DS_MAT_T:
      nelem = 3*ds->ld;
      break;
    case DS_MAT_D:
      nelem = ds->ld;
      break;
    case DS_MAT_X:
      nelem = ds->ld*n;
      break;
    case DS_MAT_Y:
      nelem = ds->ld*n;
      break;
    default:
      nelem = n*n;
  }
  if (isreal) {
    sz = nelem*sizeof(PetscReal);
    if (ds->rmat[m]) {
      ierr = PetscFree(ds->rmat[m]);CHKERRQ(ierr);
    } else {
      ierr = PetscLogObjectMemory((PetscObject)ds,sz);CHKERRQ(ierr);
    }
    ierr = PetscCalloc1(nelem,&ds->rmat[m]);CHKERRQ(ierr);
  } else {
    sz = nelem*sizeof(PetscScalar);
    if (ds->mat[m]) {
      ierr = PetscFree(ds->mat[m]);CHKERRQ(ierr);
    } else {
      ierr = PetscLogObjectMemory((PetscObject)ds,sz);CHKERRQ(ierr);
    }
    ierr = PetscCalloc1(nelem,&ds->mat[m]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSAllocateWork_Private(DS ds,PetscInt s,PetscInt r,PetscInt i)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (s>ds->lwork) {
    ierr = PetscFree(ds->work);CHKERRQ(ierr);
    ierr = PetscMalloc1(s,&ds->work);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)ds,(s-ds->lwork)*sizeof(PetscScalar));CHKERRQ(ierr);
    ds->lwork = s;
  }
  if (r>ds->lrwork) {
    ierr = PetscFree(ds->rwork);CHKERRQ(ierr);
    ierr = PetscMalloc1(r,&ds->rwork);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)ds,(r-ds->lrwork)*sizeof(PetscReal));CHKERRQ(ierr);
    ds->lrwork = r;
  }
  if (i>ds->liwork) {
    ierr = PetscFree(ds->iwork);CHKERRQ(ierr);
    ierr = PetscMalloc1(i,&ds->iwork);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)ds,(i-ds->liwork)*sizeof(PetscBLASInt));CHKERRQ(ierr);
    ds->liwork = i;
  }
  PetscFunctionReturn(0);
}

/*@C
   DSViewMat - Prints one of the internal DS matrices.

   Collective on ds

   Input Parameters:
+  ds     - the direct solver context
.  viewer - visualization context
-  m      - matrix to display

   Note:
   Works only for ascii viewers. Set the viewer in Matlab format if
   want to paste into Matlab.

   Level: developer

.seealso: DSView()
@*/
PetscErrorCode DSViewMat(DS ds,PetscViewer viewer,DSMatType m)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,rows,cols;
  PetscScalar       *v;
  PetscViewerFormat format;
#if defined(PETSC_USE_COMPLEX)
  PetscBool         allreal = PETSC_TRUE;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(ds,m,3);
  DSCheckValidMat(ds,m,3);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ds),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(ds,1,viewer,2);
  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
  ierr = DSMatGetSize(ds,m,&rows,&cols);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  /* determine if matrix has all real values */
  v = ds->mat[m];
  for (i=0;i<rows;i++)
    for (j=0;j<cols;j++)
      if (PetscImaginaryPart(v[i+j*ds->ld])) { allreal = PETSC_FALSE; break; }
#endif
  if (format == PETSC_VIEWER_ASCII_MATLAB) {
    ierr = PetscViewerASCIIPrintf(viewer,"%% Size = %D %D\n",rows,cols);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%s = [\n",DSMatName[m]);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerASCIIPrintf(viewer,"Matrix %s =\n",DSMatName[m]);CHKERRQ(ierr);
  }

  for (i=0;i<rows;i++) {
    v = ds->mat[m]+i;
    for (j=0;j<cols;j++) {
#if defined(PETSC_USE_COMPLEX)
      if (allreal) {
        ierr = PetscViewerASCIIPrintf(viewer,"%18.16e ",(double)PetscRealPart(*v));CHKERRQ(ierr);
      } else {
        ierr = PetscViewerASCIIPrintf(viewer,"%18.16e%+18.16ei ",(double)PetscRealPart(*v),(double)PetscImaginaryPart(*v));CHKERRQ(ierr);
      }
#else
      ierr = PetscViewerASCIIPrintf(viewer,"%18.16e ",(double)*v);CHKERRQ(ierr);
#endif
      v += ds->ld;
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }

  if (format == PETSC_VIEWER_ASCII_MATLAB) {
    ierr = PetscViewerASCIIPrintf(viewer,"];\n");CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DSSortEigenvalues_Private(DS ds,PetscScalar *wr,PetscScalar *wi,PetscInt *perm,PetscBool isghiep)
{
  PetscErrorCode ierr;
  PetscScalar    re,im,wi0;
  PetscInt       n,i,j,result,tmp1,tmp2=0,d=1;

  PetscFunctionBegin;
  n = ds->t;   /* sort only first t pairs if truncated */
  /* insertion sort */
  i=ds->l+1;
#if !defined(PETSC_USE_COMPLEX)
  if (wi && wi[perm[i-1]]!=0.0) i++; /* initial value is complex */
#else
  if (isghiep && PetscImaginaryPart(wr[perm[i-1]])!=0.0) i++;
#endif
  for (;i<n;i+=d) {
    re = wr[perm[i]];
    if (wi) im = wi[perm[i]];
    else im = 0.0;
    tmp1 = perm[i];
#if !defined(PETSC_USE_COMPLEX)
    if (im!=0.0) { d = 2; tmp2 = perm[i+1]; }
    else d = 1;
#else
    if (isghiep && PetscImaginaryPart(re)!=0.0) { d = 2; tmp2 = perm[i+1]; }
    else d = 1;
#endif
    j = i-1;
    if (wi) wi0 = wi[perm[j]];
    else wi0 = 0.0;
    ierr = SlepcSCCompare(ds->sc,re,im,wr[perm[j]],wi0,&result);CHKERRQ(ierr);
    while (result<0 && j>=ds->l) {
      perm[j+d] = perm[j];
      j--;
#if !defined(PETSC_USE_COMPLEX)
      if (wi && wi[perm[j+1]]!=0)
#else
      if (isghiep && PetscImaginaryPart(wr[perm[j+1]])!=0)
#endif
        { perm[j+d] = perm[j]; j--; }

      if (j>=ds->l) {
        if (wi) wi0 = wi[perm[j]];
        else wi0 = 0.0;
        ierr = SlepcSCCompare(ds->sc,re,im,wr[perm[j]],wi0,&result);CHKERRQ(ierr);
      }
    }
    perm[j+1] = tmp1;
    if (d==2) perm[j+2] = tmp2;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSortEigenvaluesReal_Private(DS ds,PetscReal *eig,PetscInt *perm)
{
  PetscErrorCode ierr;
  PetscScalar    re;
  PetscInt       i,j,result,tmp,l,n;

  PetscFunctionBegin;
  n = ds->t;   /* sort only first t pairs if truncated */
  l = ds->l;
  /* insertion sort */
  for (i=l+1;i<n;i++) {
    re = eig[perm[i]];
    j = i-1;
    ierr = SlepcSCCompare(ds->sc,re,0.0,eig[perm[j]],0.0,&result);CHKERRQ(ierr);
    while (result<0 && j>=l) {
      tmp = perm[j]; perm[j] = perm[j+1]; perm[j+1] = tmp; j--;
      if (j>=l) {
        ierr = SlepcSCCompare(ds->sc,re,0.0,eig[perm[j]],0.0,&result);CHKERRQ(ierr);
      }
    }
  }
  PetscFunctionReturn(0);
}

/*
  DSCopyMatrix_Private - Copies the trailing block of a matrix (from
  rows/columns l to n).
*/
PetscErrorCode DSCopyMatrix_Private(DS ds,DSMatType dst,DSMatType src)
{
  PetscErrorCode ierr;
  PetscInt    j,m,off,ld;
  PetscScalar *S,*D;

  PetscFunctionBegin;
  ld  = ds->ld;
  m   = ds->n-ds->l;
  off = ds->l+ds->l*ld;
  S   = ds->mat[src];
  D   = ds->mat[dst];
  for (j=0;j<m;j++) {
    ierr = PetscArraycpy(D+off+j*ld,S+off+j*ld,m);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSPermuteColumns_Private(DS ds,PetscInt l,PetscInt n,DSMatType mat,PetscInt *perm)
{
  PetscInt    i,j,k,p,ld;
  PetscScalar *Q,rtmp;

  PetscFunctionBegin;
  ld = ds->ld;
  Q  = ds->mat[mat];
  for (i=l;i<n;i++) {
    p = perm[i];
    if (p != i) {
      j = i + 1;
      while (perm[j] != i) j++;
      perm[j] = p; perm[i] = i;
      /* swap columns i and j */
      for (k=0;k<n;k++) {
        rtmp = Q[k+p*ld]; Q[k+p*ld] = Q[k+i*ld]; Q[k+i*ld] = rtmp;
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSPermuteRows_Private(DS ds,PetscInt l,PetscInt n,DSMatType mat,PetscInt *perm)
{
  PetscInt    i,j,m=ds->m,k,p,ld;
  PetscScalar *Q,rtmp;

  PetscFunctionBegin;
  if (!m) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"m was not set");
  ld = ds->ld;
  Q  = ds->mat[mat];
  for (i=l;i<n;i++) {
    p = perm[i];
    if (p != i) {
      j = i + 1;
      while (perm[j] != i) j++;
      perm[j] = p; perm[i] = i;
      /* swap rows i and j */
      for (k=0;k<m;k++) {
        rtmp = Q[p+k*ld]; Q[p+k*ld] = Q[i+k*ld]; Q[i+k*ld] = rtmp;
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSPermuteBoth_Private(DS ds,PetscInt l,PetscInt n,DSMatType mat1,DSMatType mat2,PetscInt *perm)
{
  PetscInt    i,j,m=ds->m,k,p,ld;
  PetscScalar *U,*VT,rtmp;

  PetscFunctionBegin;
  if (!m) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"m was not set");
  ld = ds->ld;
  U  = ds->mat[mat1];
  VT = ds->mat[mat2];
  for (i=l;i<n;i++) {
    p = perm[i];
    if (p != i) {
      j = i + 1;
      while (perm[j] != i) j++;
      perm[j] = p; perm[i] = i;
      /* swap columns i and j of U */
      for (k=0;k<n;k++) {
        rtmp = U[k+p*ld]; U[k+p*ld] = U[k+i*ld]; U[k+i*ld] = rtmp;
      }
      /* swap rows i and j of VT */
      for (k=0;k<m;k++) {
        rtmp = VT[p+k*ld]; VT[p+k*ld] = VT[i+k*ld]; VT[i+k*ld] = rtmp;
      }
    }
  }
  PetscFunctionReturn(0);
}

/*@
   DSSetIdentity - Copy the identity (a diagonal matrix with ones) on the
   active part of a matrix.

   Logically Collective on ds

   Input Parameters:
+  ds  - the direct solver context
-  mat - the matrix to modify

   Level: intermediate
@*/
PetscErrorCode DSSetIdentity(DS ds,DSMatType mat)
{
  PetscErrorCode ierr;
  PetscScalar    *x;
  PetscInt       i,ld,n,l;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(ds,mat,2);
  DSCheckValidMat(ds,mat,2);

  ierr = DSGetDimensions(ds,&n,NULL,&l,NULL,NULL);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(ds,&ld);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(DS_Other,ds,0,0,0);CHKERRQ(ierr);
  ierr = DSGetArray(ds,mat,&x);CHKERRQ(ierr);
  ierr = PetscArrayzero(&x[ld*l],ld*(n-l));CHKERRQ(ierr);
  for (i=l;i<n;i++) x[i+i*ld] = 1.0;
  ierr = DSRestoreArray(ds,mat,&x);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(DS_Other,ds,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   DSOrthogonalize - Orthogonalize the columns of a matrix.

   Logically Collective on ds

   Input Parameters:
+  ds   - the direct solver context
.  mat  - a matrix
-  cols - number of columns to orthogonalize (starting from column zero)

   Output Parameter:
.  lindcols - (optional) number of linearly independent columns of the matrix

   Level: developer

.seealso: DSPseudoOrthogonalize()
@*/
PetscErrorCode DSOrthogonalize(DS ds,DSMatType mat,PetscInt cols,PetscInt *lindcols)
{
  PetscErrorCode ierr;
  PetscInt       n,l,ld;
  PetscBLASInt   ld_,rA,cA,info,ltau,lw;
  PetscScalar    *A,*tau,*w,saux,dummy;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  PetscValidLogicalCollectiveEnum(ds,mat,2);
  DSCheckValidMat(ds,mat,2);
  PetscValidLogicalCollectiveInt(ds,cols,3);

  ierr = DSGetDimensions(ds,&n,NULL,&l,NULL,NULL);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(ds,&ld);CHKERRQ(ierr);
  n = n - l;
  if (cols > n) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"Invalid number of columns");
  if (n == 0 || cols == 0) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(DS_Other,ds,0,0,0);CHKERRQ(ierr);
  ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
  ierr = DSGetArray(ds,mat,&A);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(PetscMin(cols,n),&ltau);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ld,&ld_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&rA);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(cols,&cA);CHKERRQ(ierr);
  lw = -1;
  PetscStackCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&rA,&cA,A,&ld_,&dummy,&saux,&lw,&info));
  SlepcCheckLapackInfo("geqrf",info);
  lw = (PetscBLASInt)PetscRealPart(saux);
  ierr = DSAllocateWork_Private(ds,lw+ltau,0,0);CHKERRQ(ierr);
  tau = ds->work;
  w = &tau[ltau];
  PetscStackCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&rA,&cA,&A[ld*l+l],&ld_,tau,w,&lw,&info));
  SlepcCheckLapackInfo("geqrf",info);
  PetscStackCallBLAS("LAPACKorgqr",LAPACKorgqr_(&rA,&ltau,&ltau,&A[ld*l+l],&ld_,tau,w,&lw,&info));
  SlepcCheckLapackInfo("orgqr",info);
  if (lindcols) *lindcols = ltau;

  ierr = PetscFPTrapPop();CHKERRQ(ierr);
  ierr = PetscLogEventEnd(DS_Other,ds,0,0,0);CHKERRQ(ierr);
  ierr = DSRestoreArray(ds,mat,&A);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)ds);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Compute C <- a*A*B + b*C, where
    ldC, the leading dimension of C,
    ldA, the leading dimension of A,
    rA, cA, rows and columns of A,
    At, if true use the transpose of A instead,
    ldB, the leading dimension of B,
    rB, cB, rows and columns of B,
    Bt, if true use the transpose of B instead
*/
static PetscErrorCode SlepcMatDenseMult(PetscScalar *C,PetscInt _ldC,PetscScalar b,PetscScalar a,const PetscScalar *A,PetscInt _ldA,PetscInt rA,PetscInt cA,PetscBool At,const PetscScalar *B,PetscInt _ldB,PetscInt rB,PetscInt cB,PetscBool Bt)
{
  PetscErrorCode ierr;
  PetscInt       tmp;
  PetscBLASInt   m, n, k, ldA = _ldA, ldB = _ldB, ldC = _ldC;
  const char     *N = "N", *T = "C", *qA = N, *qB = N;

  PetscFunctionBegin;
  if ((rA == 0) || (cB == 0)) PetscFunctionReturn(0);
  PetscValidScalarPointer(C,1);
  PetscValidScalarPointer(A,5);
  PetscValidScalarPointer(B,10);

  /* Transpose if needed */
  if (At) tmp = rA, rA = cA, cA = tmp, qA = T;
  if (Bt) tmp = rB, rB = cB, cB = tmp, qB = T;

  /* Check size */
  if (cA != rB) SETERRQ(PETSC_COMM_SELF,1,"Matrix dimensions do not match");

  /* Do stub */
  if ((rA == 1) && (cA == 1) && (cB == 1)) {
    if (!At && !Bt) *C = *A * *B;
    else if (At && !Bt) *C = PetscConj(*A) * *B;
    else if (!At && Bt) *C = *A * PetscConj(*B);
    else *C = PetscConj(*A) * PetscConj(*B);
    m = n = k = 1;
  } else {
    m = rA; n = cB; k = cA;
    PetscStackCallBLAS("BLASgemm",BLASgemm_(qA,qB,&m,&n,&k,&a,(PetscScalar*)A,&ldA,(PetscScalar*)B,&ldB,&b,C,&ldC));
  }

  ierr = PetscLogFlops(2.0*m*n*k);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   DSPseudoOrthogonalize - Orthogonalize the columns of a matrix with Modified
   Gram-Schmidt in an indefinite inner product space defined by a signature.

   Logically Collective on ds

   Input Parameters:
+  ds   - the direct solver context
.  mat  - the matrix
.  cols - number of columns to orthogonalize (starting from column zero)
-  s    - the signature that defines the inner product

   Output Parameters:
+  lindcols - (optional) linearly independent columns of the matrix
-  ns   - (optional) the new signature of the vectors

   Note:
   After the call the matrix satisfies A'*s*A = ns.

   Level: developer

.seealso: DSOrthogonalize()
@*/
PetscErrorCode DSPseudoOrthogonalize(DS ds,DSMatType mat,PetscInt cols,PetscReal *s,PetscInt *lindcols,PetscReal *ns)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,l,n,ld;
  PetscBLASInt   one=1,rA_;
  PetscScalar    alpha,*A,*A_,*m,*h,nr0;
  PetscReal      nr_o,nr,*ns_;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  DSCheckAlloc(ds,1);
  PetscValidLogicalCollectiveEnum(ds,mat,2);
  DSCheckValidMat(ds,mat,2);
  PetscValidLogicalCollectiveInt(ds,cols,3);
  PetscValidRealPointer(s,4);
  ierr = DSGetDimensions(ds,&n,NULL,&l,NULL,NULL);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(ds,&ld);CHKERRQ(ierr);
  n = n - l;
  if (cols > n) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONG,"Invalid number of columns");
  if (n == 0 || cols == 0) PetscFunctionReturn(0);
  ierr = PetscBLASIntCast(n,&rA_);CHKERRQ(ierr);
  ierr = DSGetArray(ds,mat,&A_);CHKERRQ(ierr);
  A = &A_[ld*l+l];
  ierr = DSAllocateWork_Private(ds,n+cols,ns?0:cols,0);CHKERRQ(ierr);
  m = ds->work;
  h = &m[n];
  ns_ = ns ? ns : ds->rwork;
  ierr = PetscLogEventBegin(DS_Other,ds,0,0,0);CHKERRQ(ierr);
  for (i=0; i<cols; i++) {
    /* m <- diag(s)*A[i] */
    for (k=0; k<n; k++) m[k] = s[k]*A[k+i*ld];
    /* nr_o <- mynorm(A[i]'*m), mynorm(x) = sign(x)*sqrt(|x|) */
    ierr = SlepcMatDenseMult(&nr0,1,0.0,1.0,&A[ld*i],ld,n,1,PETSC_TRUE,m,n,n,1,PETSC_FALSE);CHKERRQ(ierr);
    nr = nr_o = PetscSign(PetscRealPart(nr0))*PetscSqrtReal(PetscAbsScalar(nr0));
    for (j=0; j<3 && i>0; j++) {
      /* h <- A[0:i-1]'*m */
      ierr = SlepcMatDenseMult(h,i,0.0,1.0,A,ld,n,i,PETSC_TRUE,m,n,n,1,PETSC_FALSE);CHKERRQ(ierr);
      /* h <- diag(ns)*h */
      for (k=0; k<i; k++) h[k] *= ns_[k];
      /* A[i] <- A[i] - A[0:i-1]*h */
      ierr = SlepcMatDenseMult(&A[ld*i],ld,1.0,-1.0,A,ld,n,i,PETSC_FALSE,h,i,i,1,PETSC_FALSE);CHKERRQ(ierr);
      /* m <- diag(s)*A[i] */
      for (k=0; k<n; k++) m[k] = s[k]*A[k+i*ld];
      /* nr_o <- mynorm(A[i]'*m) */
      ierr = SlepcMatDenseMult(&nr0,1,0.0,1.0,&A[ld*i],ld,n,1,PETSC_TRUE,m,n,n,1,PETSC_FALSE);CHKERRQ(ierr);
      nr = PetscSign(PetscRealPart(nr0))*PetscSqrtReal(PetscAbsScalar(nr0));
      if (PetscAbs(nr) < PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_SELF,1,"Linear dependency detected");
      if (PetscAbs(nr) > 0.7*PetscAbs(nr_o)) break;
      nr_o = nr;
    }
    ns_[i] = PetscSign(nr);
    /* A[i] <- A[i]/|nr| */
    alpha = 1.0/PetscAbs(nr);
    PetscStackCallBLAS("BLASscal",BLASscal_(&rA_,&alpha,&A[i*ld],&one));
  }
  ierr = PetscLogEventEnd(DS_Other,ds,0,0,0);CHKERRQ(ierr);
  ierr = DSRestoreArray(ds,mat,&A_);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)ds);CHKERRQ(ierr);
  if (lindcols) *lindcols = cols;
  PetscFunctionReturn(0);
}
