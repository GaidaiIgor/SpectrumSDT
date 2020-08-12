/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/dsimpl.h>       /*I "slepcds.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  PetscInt  d;              /* polynomial degree */
  PetscReal *pbc;           /* polynomial basis coefficients */
} DS_PEP;

PetscErrorCode DSAllocate_PEP(DS ds,PetscInt ld)
{
  PetscErrorCode ierr;
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (!ctx->d) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONGSTATE,"DSPEP requires specifying the polynomial degree via DSPEPSetDegree()");
  ierr = DSAllocateMat_Private(ds,DS_MAT_X);CHKERRQ(ierr);
  ierr = DSAllocateMat_Private(ds,DS_MAT_Y);CHKERRQ(ierr);
  for (i=0;i<=ctx->d;i++) {
    ierr = DSAllocateMat_Private(ds,DSMatExtra[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(ds->perm);CHKERRQ(ierr);
  ierr = PetscMalloc1(ld*ctx->d,&ds->perm);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)ds,ld*ctx->d*sizeof(PetscInt));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DSView_PEP(DS ds,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  DS_PEP            *ctx = (DS_PEP*)ds->data;
  PetscViewerFormat format;
  PetscInt          i;

  PetscFunctionBegin;
  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"polynomial degree: %D\n",ctx->d);CHKERRQ(ierr);
  if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
  for (i=0;i<=ctx->d;i++) {
    ierr = DSViewMat(ds,viewer,DSMatExtra[i]);CHKERRQ(ierr);
  }
  if (ds->state>DS_STATE_INTERMEDIATE) { ierr = DSViewMat(ds,viewer,DS_MAT_X);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode DSVectors_PEP(DS ds,DSMatType mat,PetscInt *j,PetscReal *rnorm)
{
  PetscFunctionBegin;
  if (rnorm) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_SUP,"Not implemented yet");
  switch (mat) {
    case DS_MAT_X:
      break;
    case DS_MAT_Y:
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mat parameter");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSort_PEP(DS ds,PetscScalar *wr,PetscScalar *wi,PetscScalar *rr,PetscScalar *ri,PetscInt *kout)
{
  PetscErrorCode ierr;
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       n,i,j,k,p,*perm,told,ld;
  PetscScalar    *A,*X,*Y,rtmp,rtmp2;

  PetscFunctionBegin;
  if (!ds->sc) PetscFunctionReturn(0);
  n = ds->n*ctx->d;
  A  = ds->mat[DS_MAT_A];
  perm = ds->perm;
  for (i=0;i<n;i++) perm[i] = i;
  told = ds->t;
  ds->t = n;  /* force the sorting routines to consider d*n eigenvalues */
  if (rr) {
    ierr = DSSortEigenvalues_Private(ds,rr,ri,perm,PETSC_FALSE);CHKERRQ(ierr);
  } else {
    ierr = DSSortEigenvalues_Private(ds,wr,wi,perm,PETSC_FALSE);CHKERRQ(ierr);
  }
  ds->t = told;  /* restore value of t */
  for (i=0;i<n;i++) A[i]  = wr[perm[i]];
  for (i=0;i<n;i++) wr[i] = A[i];
  for (i=0;i<n;i++) A[i]  = wi[perm[i]];
  for (i=0;i<n;i++) wi[i] = A[i];
  /* cannot use DSPermuteColumns_Private() since matrix is not square */
  ld = ds->ld;
  X  = ds->mat[DS_MAT_X];
  Y  = ds->mat[DS_MAT_Y];
  for (i=0;i<n;i++) {
    p = perm[i];
    if (p != i) {
      j = i + 1;
      while (perm[j] != i) j++;
      perm[j] = p; perm[i] = i;
      /* swap columns i and j */
      for (k=0;k<ds->n;k++) {
        rtmp  = X[k+p*ld]; X[k+p*ld] = X[k+i*ld]; X[k+i*ld] = rtmp;
        rtmp2 = Y[k+p*ld]; Y[k+p*ld] = Y[k+i*ld]; Y[k+i*ld] = rtmp2;
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSolve_PEP_QZ(DS ds,PetscScalar *wr,PetscScalar *wi)
{
  PetscErrorCode ierr;
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       i,j,k,off;
  PetscScalar    *A,*B,*W,*X,*U,*Y,*E,*work,*beta,norm;
  PetscReal      *ca,*cb,*cg;
  PetscBLASInt   info,n,ldd,nd,lrwork=0,lwork,one=1;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork;
#else
  PetscScalar    norm0;
#endif

  PetscFunctionBegin;
  if (!ds->mat[DS_MAT_A]) {
    ierr = DSAllocateMat_Private(ds,DS_MAT_A);CHKERRQ(ierr);
  }
  if (!ds->mat[DS_MAT_B]) {
    ierr = DSAllocateMat_Private(ds,DS_MAT_B);CHKERRQ(ierr);
  }
  if (!ds->mat[DS_MAT_W]) {
    ierr = DSAllocateMat_Private(ds,DS_MAT_W);CHKERRQ(ierr);
  }
  if (!ds->mat[DS_MAT_U]) {
    ierr = DSAllocateMat_Private(ds,DS_MAT_U);CHKERRQ(ierr);
  }
  ierr = PetscBLASIntCast(ds->n*ctx->d,&nd);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->n,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ds->ld*ctx->d,&ldd);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscBLASIntCast(nd+2*nd,&lwork);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(8*nd,&lrwork);CHKERRQ(ierr);
#else
  ierr = PetscBLASIntCast(nd+8*nd,&lwork);CHKERRQ(ierr);
#endif
  ierr = DSAllocateWork_Private(ds,lwork,lrwork,0);CHKERRQ(ierr);
  beta = ds->work;
  work = ds->work + nd;
  lwork -= nd;
  A = ds->mat[DS_MAT_A];
  B = ds->mat[DS_MAT_B];
  W = ds->mat[DS_MAT_W];
  U = ds->mat[DS_MAT_U];
  X = ds->mat[DS_MAT_X];
  Y = ds->mat[DS_MAT_Y];
  E = ds->mat[DSMatExtra[ctx->d]];

  /* build matrices A and B of the linearization */
  ierr = PetscArrayzero(A,ldd*ldd);CHKERRQ(ierr);
  if (!ctx->pbc) { /* monomial basis */
    for (i=0;i<nd-ds->n;i++) A[i+(i+ds->n)*ldd] = 1.0;
    for (i=0;i<ctx->d;i++) {
      off = i*ds->n*ldd+(ctx->d-1)*ds->n;
      for (j=0;j<ds->n;j++) {
        ierr = PetscArraycpy(A+off+j*ldd,ds->mat[DSMatExtra[i]]+j*ds->ld,ds->n);CHKERRQ(ierr);
      }
    }
  } else {
    ca = ctx->pbc;
    cb = ca+ctx->d+1;
    cg = cb+ctx->d+1;
    for (i=0;i<ds->n;i++) {
      A[i+(i+ds->n)*ldd] = ca[0];
      A[i+i*ldd] = cb[0];
    }
    for (;i<nd-ds->n;i++) {
      j = i/ds->n;
      A[i+(i+ds->n)*ldd] = ca[j];
      A[i+i*ldd] = cb[j];
      A[i+(i-ds->n)*ldd] = cg[j];
    }
    for (i=0;i<ctx->d-2;i++) {
      off = i*ds->n*ldd+(ctx->d-1)*ds->n;
      for (j=0;j<ds->n;j++)
        for (k=0;k<ds->n;k++)
          *(A+off+j*ldd+k) = *(ds->mat[DSMatExtra[i]]+j*ds->ld+k)*ca[ctx->d-1];
    }
    off = i*ds->n*ldd+(ctx->d-1)*ds->n;
    for (j=0;j<ds->n;j++)
      for (k=0;k<ds->n;k++)
        *(A+off+j*ldd+k) = *(ds->mat[DSMatExtra[i]]+j*ds->ld+k)*ca[ctx->d-1]-E[j*ds->ld+k]*cg[ctx->d-1];
    off = (++i)*ds->n*ldd+(ctx->d-1)*ds->n;
    for (j=0;j<ds->n;j++)
      for (k=0;k<ds->n;k++)
        *(A+off+j*ldd+k) = *(ds->mat[DSMatExtra[i]]+j*ds->ld+k)*ca[ctx->d-1]-E[j*ds->ld+k]*cb[ctx->d-1];
  }
  ierr = PetscArrayzero(B,ldd*ldd);CHKERRQ(ierr);
  for (i=0;i<nd-ds->n;i++) B[i+i*ldd] = 1.0;
  off = (ctx->d-1)*ds->n*(ldd+1);
  for (j=0;j<ds->n;j++) {
    for (i=0;i<ds->n;i++) B[off+i+j*ldd] = -E[i+j*ds->ld];
  }

  /* solve generalized eigenproblem */
#if defined(PETSC_USE_COMPLEX)
  rwork = ds->rwork;
  PetscStackCallBLAS("LAPACKggev",LAPACKggev_("V","V",&nd,A,&ldd,B,&ldd,wr,beta,U,&ldd,W,&ldd,work,&lwork,rwork,&info));
#else
  PetscStackCallBLAS("LAPACKggev",LAPACKggev_("V","V",&nd,A,&ldd,B,&ldd,wr,wi,beta,U,&ldd,W,&ldd,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("ggev",info);

  /* copy eigenvalues */
  for (i=0;i<nd;i++) {
    if (beta[i]==0.0) wr[i] = (PetscRealPart(wr[i])>0.0)? PETSC_MAX_REAL: PETSC_MIN_REAL;
    else wr[i] /= beta[i];
#if !defined(PETSC_USE_COMPLEX)
    if (beta[i]==0.0) wi[i] = 0.0;
    else wi[i] /= beta[i];
#else
    if (wi) wi[i] = 0.0;
#endif
  }

  /* copy and normalize eigenvectors */
  for (j=0;j<nd;j++) {
    ierr = PetscArraycpy(X+j*ds->ld,W+j*ldd,ds->n);CHKERRQ(ierr);
    ierr = PetscArraycpy(Y+j*ds->ld,U+ds->n*(ctx->d-1)+j*ldd,ds->n);CHKERRQ(ierr);
  }
  for (j=0;j<nd;j++) {
#if !defined(PETSC_USE_COMPLEX)
    if (wi[j] != 0.0) {
      norm = BLASnrm2_(&n,X+j*ds->ld,&one);
      norm0 = BLASnrm2_(&n,X+(j+1)*ds->ld,&one);
      norm = 1.0/SlepcAbsEigenvalue(norm,norm0);
      PetscStackCallBLAS("BLASscal",BLASscal_(&n,&norm,X+j*ds->ld,&one));
      PetscStackCallBLAS("BLASscal",BLASscal_(&n,&norm,X+(j+1)*ds->ld,&one));
      norm = BLASnrm2_(&n,Y+j*ds->ld,&one);
      norm0 = BLASnrm2_(&n,Y+(j+1)*ds->ld,&one);
      norm = 1.0/SlepcAbsEigenvalue(norm,norm0);
      PetscStackCallBLAS("BLASscal",BLASscal_(&n,&norm,Y+j*ds->ld,&one));
      PetscStackCallBLAS("BLASscal",BLASscal_(&n,&norm,Y+(j+1)*ds->ld,&one));
      j++;
    } else
#endif
    {
      norm = 1.0/BLASnrm2_(&n,X+j*ds->ld,&one);
      PetscStackCallBLAS("BLASscal",BLASscal_(&n,&norm,X+j*ds->ld,&one));
      norm = 1.0/BLASnrm2_(&n,Y+j*ds->ld,&one);
      PetscStackCallBLAS("BLASscal",BLASscal_(&n,&norm,Y+j*ds->ld,&one));
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSSynchronize_PEP(DS ds,PetscScalar eigr[],PetscScalar eigi[])
{
  PetscErrorCode ierr;
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       ld=ds->ld,k=0;
  PetscMPIInt    ldnd,rank,off=0,size,dn;

  PetscFunctionBegin;
  if (ds->state>=DS_STATE_CONDENSED) k += 2*ctx->d*ds->n*ld;
  if (eigr) k += ctx->d*ds->n;
  if (eigi) k += ctx->d*ds->n;
  ierr = DSAllocateWork_Private(ds,k,0,0);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(k*sizeof(PetscScalar),&size);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ds->n*ctx->d*ld,&ldnd);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(ctx->d*ds->n,&dn);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)ds),&rank);CHKERRQ(ierr);
  if (!rank) {
    if (ds->state>=DS_STATE_CONDENSED) {
      ierr = MPI_Pack(ds->mat[DS_MAT_X],ldnd,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
      ierr = MPI_Pack(ds->mat[DS_MAT_Y],ldnd,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigr) {
      ierr = MPI_Pack(eigr,dn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigi) {
      ierr = MPI_Pack(eigi,dn,MPIU_SCALAR,ds->work,size,&off,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
  }
  ierr = MPI_Bcast(ds->work,size,MPI_BYTE,0,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
  if (rank) {
    if (ds->state>=DS_STATE_CONDENSED) {
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_X],ldnd,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
      ierr = MPI_Unpack(ds->work,size,&off,ds->mat[DS_MAT_Y],ldnd,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigr) {
      ierr = MPI_Unpack(ds->work,size,&off,eigr,dn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
    if (eigi) {
      ierr = MPI_Unpack(ds->work,size,&off,eigi,dn,MPIU_SCALAR,PetscObjectComm((PetscObject)ds));CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DSPEPSetDegree_PEP(DS ds,PetscInt d)
{
  DS_PEP *ctx = (DS_PEP*)ds->data;

  PetscFunctionBegin;
  if (d<0) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"The degree must be a non-negative integer");
  if (d>=DS_NUM_EXTRA) SETERRQ1(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_OUTOFRANGE,"Only implemented for polynomials of degree at most %D",DS_NUM_EXTRA-1);
  ctx->d = d;
  PetscFunctionReturn(0);
}

/*@
   DSPEPSetDegree - Sets the polynomial degree for a DSPEP.

   Logically Collective on ds

   Input Parameters:
+  ds - the direct solver context
-  d  - the degree

   Level: intermediate

.seealso: DSPEPGetDegree()
@*/
PetscErrorCode DSPEPSetDegree(DS ds,PetscInt d)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidLogicalCollectiveInt(ds,d,2);
  ierr = PetscTryMethod(ds,"DSPEPSetDegree_C",(DS,PetscInt),(ds,d));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DSPEPGetDegree_PEP(DS ds,PetscInt *d)
{
  DS_PEP *ctx = (DS_PEP*)ds->data;

  PetscFunctionBegin;
  *d = ctx->d;
  PetscFunctionReturn(0);
}

/*@
   DSPEPGetDegree - Returns the polynomial degree for a DSPEP.

   Not collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
.  d - the degree

   Level: intermediate

.seealso: DSPEPSetDegree()
@*/
PetscErrorCode DSPEPGetDegree(DS ds,PetscInt *d)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidIntPointer(d,2);
  ierr = PetscUseMethod(ds,"DSPEPGetDegree_C",(DS,PetscInt*),(ds,d));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DSPEPSetCoefficients_PEP(DS ds,PetscReal *pbc)
{
  PetscErrorCode ierr;
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (!ctx->d) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONGSTATE,"Must first specify the polynomial degree via DSPEPSetDegree()");
  if (ctx->pbc) { ierr = PetscFree(ctx->pbc);CHKERRQ(ierr); }
  ierr = PetscMalloc1(3*(ctx->d+1),&ctx->pbc);CHKERRQ(ierr);
  for (i=0;i<3*(ctx->d+1);i++) ctx->pbc[i] = pbc[i];
  ds->state = DS_STATE_RAW;
  PetscFunctionReturn(0);
}

/*@C
   DSPEPSetCoefficients - Sets the polynomial basis coefficients for a DSPEP.

   Logically Collective on ds

   Input Parameters:
+  ds  - the direct solver context
-  pbc - the polynomial basis coefficients

   Notes:
   This function is required only in the case of a polynomial specified in a
   non-monomial basis, to provide the coefficients that will be used
   during the linearization, multiplying the identity blocks on the three main
   diagonal blocks. Depending on the polynomial basis (Chebyshev, Legendre, ...)
   the coefficients must be different.

   There must be a total of 3*(d+1) coefficients, where d is the degree of the
   polynomial. The coefficients are arranged in three groups: alpha, beta, and
   gamma, according to the definition of the three-term recurrence. In the case
   of the monomial basis, alpha=1 and beta=gamma=0, in which case it is not
   necessary to invoke this function.

   Level: advanced

.seealso: DSPEPGetCoefficients(), DSPEPSetDegree()
@*/
PetscErrorCode DSPEPSetCoefficients(DS ds,PetscReal *pbc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  ierr = PetscTryMethod(ds,"DSPEPSetCoefficients_C",(DS,PetscReal*),(ds,pbc));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DSPEPGetCoefficients_PEP(DS ds,PetscReal **pbc)
{
  PetscErrorCode ierr;
  DS_PEP         *ctx = (DS_PEP*)ds->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (!ctx->d) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONGSTATE,"Must first specify the polynomial degree via DSPEPSetDegree()");
  ierr = PetscCalloc1(3*(ctx->d+1),pbc);CHKERRQ(ierr);
  if (ctx->pbc) for (i=0;i<3*(ctx->d+1);i++) (*pbc)[i] = ctx->pbc[i];
  else for (i=0;i<ctx->d+1;i++) (*pbc)[i] = 1.0;
  PetscFunctionReturn(0);
}

/*@C
   DSPEPGetCoefficients - Returns the polynomial basis coefficients for a DSPEP.

   Not collective

   Input Parameter:
.  ds - the direct solver context

   Output Parameters:
.  pbc - the polynomial basis coefficients

   Note:
   The returned array has length 3*(d+1) and should be freed by the user.

   Fortran Note:
   The calling sequence from Fortran is
.vb
   DSPEPGetCoefficients(eps,pbc,ierr)
   double precision pbc(d+1) output
.ve

   Level: advanced

.seealso: DSPEPSetCoefficients()
@*/
PetscErrorCode DSPEPGetCoefficients(DS ds,PetscReal **pbc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(ds,DS_CLASSID,1);
  PetscValidPointer(pbc,2);
  ierr = PetscUseMethod(ds,"DSPEPGetCoefficients_C",(DS,PetscReal**),(ds,pbc));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DSDestroy_PEP(DS ds)
{
  PetscErrorCode ierr;
  DS_PEP         *ctx = (DS_PEP*)ds->data;

  PetscFunctionBegin;
  if (ctx->pbc) { ierr = PetscFree(ctx->pbc);CHKERRQ(ierr); }
  ierr = PetscFree(ds->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ds,"DSPEPSetDegree_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ds,"DSPEPGetDegree_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ds,"DSPEPSetCoefficients_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ds,"DSPEPGetCoefficients_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DSMatGetSize_PEP(DS ds,DSMatType t,PetscInt *rows,PetscInt *cols)
{
  DS_PEP *ctx = (DS_PEP*)ds->data;

  PetscFunctionBegin;
  if (!ctx->d) SETERRQ(PetscObjectComm((PetscObject)ds),PETSC_ERR_ARG_WRONGSTATE,"DSPEP requires specifying the polynomial degree via DSPEPSetDegree()");
  *rows = ds->n;
  if (t==DS_MAT_A || t==DS_MAT_B || t==DS_MAT_W || t==DS_MAT_U) *rows *= ctx->d;
  *cols = ds->n;
  if (t==DS_MAT_A || t==DS_MAT_B || t==DS_MAT_W || t==DS_MAT_U || t==DS_MAT_X || t==DS_MAT_Y) *cols *= ctx->d;
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode DSCreate_PEP(DS ds)
{
  DS_PEP         *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(ds,&ctx);CHKERRQ(ierr);
  ds->data = (void*)ctx;

  ds->ops->allocate      = DSAllocate_PEP;
  ds->ops->view          = DSView_PEP;
  ds->ops->vectors       = DSVectors_PEP;
  ds->ops->solve[0]      = DSSolve_PEP_QZ;
  ds->ops->sort          = DSSort_PEP;
  ds->ops->synchronize   = DSSynchronize_PEP;
  ds->ops->destroy       = DSDestroy_PEP;
  ds->ops->matgetsize    = DSMatGetSize_PEP;
  ierr = PetscObjectComposeFunction((PetscObject)ds,"DSPEPSetDegree_C",DSPEPSetDegree_PEP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ds,"DSPEPGetDegree_C",DSPEPGetDegree_PEP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ds,"DSPEPSetCoefficients_C",DSPEPSetCoefficients_PEP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)ds,"DSPEPGetCoefficients_C",DSPEPGetCoefficients_PEP);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

