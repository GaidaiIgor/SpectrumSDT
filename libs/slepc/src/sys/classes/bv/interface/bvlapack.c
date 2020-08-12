/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV private kernels that use the LAPACK
*/

#include <slepc/private/bvimpl.h>
#include <slepcblaslapack.h>

/*
    Compute ||A|| for an mxn matrix
*/
PetscErrorCode BVNorm_LAPACK_Private(BV bv,PetscInt m_,PetscInt n_,const PetscScalar *A,NormType type,PetscReal *nrm,PetscBool mpi)
{
  PetscErrorCode ierr;
  PetscBLASInt   m,n,i,j;
  PetscMPIInt    len;
  PetscReal      lnrm,*rwork=NULL,*rwork2=NULL;

  PetscFunctionBegin;
  ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m_,&m);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  if (type==NORM_FROBENIUS || type==NORM_2) {
    lnrm = LAPACKlange_("F",&m,&n,(PetscScalar*)A,&m,rwork);
    if (mpi) {
      lnrm = lnrm*lnrm;
      ierr = MPI_Allreduce(&lnrm,nrm,1,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)bv));CHKERRQ(ierr);
      *nrm = PetscSqrtReal(*nrm);
    } else *nrm = lnrm;
    ierr = PetscLogFlops(2.0*m*n);CHKERRQ(ierr);
  } else if (type==NORM_1) {
    if (mpi) {
      ierr = BVAllocateWork_Private(bv,2*n_);CHKERRQ(ierr);
      rwork = (PetscReal*)bv->work;
      rwork2 = rwork+n_;
      ierr = PetscArrayzero(rwork,n_);CHKERRQ(ierr);
      ierr = PetscArrayzero(rwork2,n_);CHKERRQ(ierr);
      for (j=0;j<n_;j++) {
        for (i=0;i<m_;i++) {
          rwork[j] += PetscAbsScalar(A[i+j*m_]);
        }
      }
      ierr = PetscMPIIntCast(n_,&len);CHKERRQ(ierr);
      ierr = MPI_Allreduce(rwork,rwork2,len,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)bv));CHKERRQ(ierr);
      *nrm = 0.0;
      for (j=0;j<n_;j++) if (rwork2[j] > *nrm) *nrm = rwork2[j];
    } else {
      *nrm = LAPACKlange_("O",&m,&n,(PetscScalar*)A,&m,rwork);
    }
    ierr = PetscLogFlops(1.0*m*n);CHKERRQ(ierr);
  } else if (type==NORM_INFINITY) {
    ierr = BVAllocateWork_Private(bv,m_);CHKERRQ(ierr);
    rwork = (PetscReal*)bv->work;
    lnrm = LAPACKlange_("I",&m,&n,(PetscScalar*)A,&m,rwork);
    if (mpi) {
      ierr = MPI_Allreduce(&lnrm,nrm,1,MPIU_REAL,MPIU_MAX,PetscObjectComm((PetscObject)bv));CHKERRQ(ierr);
    } else *nrm = lnrm;
    ierr = PetscLogFlops(1.0*m*n);CHKERRQ(ierr);
  }
  ierr = PetscFPTrapPop();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Compute the upper Cholesky factor in R and its inverse in S.
   If S == R then the inverse overwrites the Cholesky factor.
 */
PetscErrorCode BVMatCholInv_LAPACK_Private(BV bv,Mat R,Mat S)
{
  PetscErrorCode ierr;
  PetscInt       i,k,l,n,m,ld,lds;
  PetscScalar    *pR,*pS;
  PetscBLASInt   info,n_,m_,ld_,lds_;

  PetscFunctionBegin;
  l = bv->l;
  k = bv->k;
  ierr = MatGetSize(R,&m,NULL);CHKERRQ(ierr);
  n = k-l;
  ierr = PetscBLASIntCast(m,&m_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ld  = m;
  ld_ = m_;
  ierr = MatDenseGetArray(R,&pR);CHKERRQ(ierr);

  if (S==R) {
    ierr = BVAllocateWork_Private(bv,m*k);CHKERRQ(ierr);
    pS = bv->work;
    lds = ld;
    lds_ = ld_;
  } else {
    ierr = MatDenseGetArray(S,&pS);CHKERRQ(ierr);
    ierr = MatGetSize(S,&lds,NULL);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  }

  /* save a copy of matrix in S */
  for (i=l;i<k;i++) {
    ierr = PetscArraycpy(pS+i*lds+l,pR+i*ld+l,n);CHKERRQ(ierr);
  }

  /* compute upper Cholesky factor in R */
  PetscStackCallBLAS("LAPACKpotrf",LAPACKpotrf_("U",&n_,pR+l*ld+l,&ld_,&info));
  ierr = PetscLogFlops((1.0*n*n*n)/3.0);CHKERRQ(ierr);

  if (info) {  /* LAPACKpotrf failed, retry on diagonally perturbed matrix */
    for (i=l;i<k;i++) {
      ierr = PetscArraycpy(pR+i*ld+l,pS+i*lds+l,n);CHKERRQ(ierr);
      pR[i+i*ld] += 50.0*PETSC_MACHINE_EPSILON;
    }
    PetscStackCallBLAS("LAPACKpotrf",LAPACKpotrf_("U",&n_,pR+l*ld+l,&ld_,&info));
    SlepcCheckLapackInfo("potrf",info);
    ierr = PetscLogFlops((1.0*n*n*n)/3.0);CHKERRQ(ierr);
  }

  /* compute S = inv(R) */
  if (S==R) {
    PetscStackCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,pR+l*ld+l,&ld_,&info));
  } else {
    ierr = PetscArrayzero(pS+l*lds,(k-l)*k);CHKERRQ(ierr);
    for (i=l;i<k;i++) {
      ierr = PetscArraycpy(pS+i*lds+l,pR+i*ld+l,n);CHKERRQ(ierr);
    }
    PetscStackCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,pS+l*lds+l,&lds_,&info));
  }
  SlepcCheckLapackInfo("trtri",info);
  ierr = PetscLogFlops(0.33*n*n*n);CHKERRQ(ierr);

  /* Zero out entries below the diagonal */
  for (i=l;i<k-1;i++) {
    ierr = PetscArrayzero(pR+i*ld+i+1,(k-i-1));CHKERRQ(ierr);
    if (S!=R) { ierr = PetscArrayzero(pS+i*lds+i+1,(k-i-1));CHKERRQ(ierr); }
  }
  ierr = MatDenseRestoreArray(R,&pR);CHKERRQ(ierr);
  if (S!=R) { ierr = MatDenseRestoreArray(S,&pS);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   Compute the inverse of an upper triangular matrix R, store it in S.
   If S == R then the inverse overwrites R.
 */
PetscErrorCode BVMatTriInv_LAPACK_Private(BV bv,Mat R,Mat S)
{
  PetscErrorCode ierr;
  PetscInt       i,k,l,n,m,ld,lds;
  PetscScalar    *pR,*pS;
  PetscBLASInt   info,n_,m_,ld_,lds_;

  PetscFunctionBegin;
  l = bv->l;
  k = bv->k;
  ierr = MatGetSize(R,&m,NULL);CHKERRQ(ierr);
  n = k-l;
  ierr = PetscBLASIntCast(m,&m_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ld  = m;
  ld_ = m_;
  ierr = MatDenseGetArray(R,&pR);CHKERRQ(ierr);

  if (S==R) {
    ierr = BVAllocateWork_Private(bv,m*k);CHKERRQ(ierr);
    pS = bv->work;
    lds = ld;
    lds_ = ld_;
  } else {
    ierr = MatDenseGetArray(S,&pS);CHKERRQ(ierr);
    ierr = MatGetSize(S,&lds,NULL);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  }

  /* compute S = inv(R) */
  if (S==R) {
    PetscStackCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,pR+l*ld+l,&ld_,&info));
  } else {
    ierr = PetscArrayzero(pS+l*lds,(k-l)*k);CHKERRQ(ierr);
    for (i=l;i<k;i++) {
      ierr = PetscArraycpy(pS+i*lds+l,pR+i*ld+l,n);CHKERRQ(ierr);
    }
    PetscStackCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,pS+l*lds+l,&lds_,&info));
  }
  SlepcCheckLapackInfo("trtri",info);
  ierr = PetscLogFlops(0.33*n*n*n);CHKERRQ(ierr);

  ierr = MatDenseRestoreArray(R,&pR);CHKERRQ(ierr);
  if (S!=R) { ierr = MatDenseRestoreArray(S,&pS);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   Compute the matrix to be used for post-multiplying the basis in the SVQB
   block orthogonalization method.
   On input R = V'*V, on output S = D*U*Lambda^{-1/2} where (U,Lambda) is
   the eigendecomposition of D*R*D with D=diag(R)^{-1/2}.
   If S == R then the result overwrites R.
 */
PetscErrorCode BVMatSVQB_LAPACK_Private(BV bv,Mat R,Mat S)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,l,n,m,ld,lds;
  PetscScalar    *pR,*pS,*D,*work,a;
  PetscReal      *eig,dummy;
  PetscBLASInt   info,lwork,n_,m_,ld_,lds_;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork,rdummy;
#endif

  PetscFunctionBegin;
  l = bv->l;
  k = bv->k;
  ierr = MatGetSize(R,&m,NULL);CHKERRQ(ierr);
  n = k-l;
  ierr = PetscBLASIntCast(m,&m_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ld  = m;
  ld_ = m_;
  ierr = MatDenseGetArray(R,&pR);CHKERRQ(ierr);

  if (S==R) {
    pS = pR;
    lds = ld;
    lds_ = ld_;
  } else {
    ierr = MatDenseGetArray(S,&pS);CHKERRQ(ierr);
    ierr = MatGetSize(S,&lds,NULL);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  }

  /* workspace query and memory allocation */
  lwork = -1;
#if defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n_,pS,&lds_,&dummy,&a,&lwork,&rdummy,&info));
  ierr = PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork);CHKERRQ(ierr);
  ierr = PetscMalloc4(n,&eig,n,&D,lwork,&work,PetscMax(1,3*n-2),&rwork);CHKERRQ(ierr);
#else
  PetscStackCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n_,pS,&lds_,&dummy,&a,&lwork,&info));
  ierr = PetscBLASIntCast((PetscInt)PetscRealPart(a),&lwork);CHKERRQ(ierr);
  ierr = PetscMalloc3(n,&eig,n,&D,lwork,&work);CHKERRQ(ierr);
#endif

  /* copy and scale matrix */
  for (i=l;i<k;i++) D[i-l] = 1.0/PetscSqrtReal(PetscRealPart(pR[i+i*ld]));
  for (i=l;i<k;i++) for (j=l;j<k;j++) pS[i+j*lds] = pR[i+j*ld]*D[i-l];
  for (j=l;j<k;j++) for (i=l;i<k;i++) pS[i+j*lds] *= D[j-l];

  /* compute eigendecomposition */
#if defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n_,pS+l*lds+l,&lds_,eig,work,&lwork,rwork,&info));
#else
  PetscStackCallBLAS("LAPACKsyev",LAPACKsyev_("V","L",&n_,pS+l*lds+l,&lds_,eig,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("syev",info);

  if (S!=R) {   /* R = U' */
    for (i=l;i<k;i++) for (j=l;j<k;j++) pR[i+j*ld] = pS[j+i*lds];
  }

  /* compute S = D*U*Lambda^{-1/2} */
  for (i=l;i<k;i++) for (j=l;j<k;j++) pS[i+j*lds] *= D[i-l];
  for (j=l;j<k;j++) for (i=l;i<k;i++) pS[i+j*lds] /= PetscSqrtReal(eig[j-l]);

  if (S!=R) {   /* compute R = inv(S) = Lambda^{1/2}*U'/D */
    for (i=l;i<k;i++) for (j=l;j<k;j++) pR[i+j*ld] *= PetscSqrtReal(eig[i-l]);
    for (j=l;j<k;j++) for (i=l;i<k;i++) pR[i+j*ld] /= D[j-l];
  }

#if defined(PETSC_USE_COMPLEX)
  ierr = PetscFree4(eig,D,work,rwork);CHKERRQ(ierr);
#else
  ierr = PetscFree3(eig,D,work);CHKERRQ(ierr);
#endif
  ierr = PetscLogFlops(9.0*n*n*n);CHKERRQ(ierr);

  ierr = MatDenseRestoreArray(R,&pR);CHKERRQ(ierr);
  if (S!=R) { ierr = MatDenseRestoreArray(S,&pS);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
    QR factorization of an mxn matrix via parallel TSQR
*/
PetscErrorCode BVOrthogonalize_LAPACK_TSQR(BV bv,PetscInt m_,PetscInt n_,PetscScalar *Q,PetscScalar *R,PetscInt ldr)
{
  PetscErrorCode ierr;
  PetscInt       level,plevel,nlevels,powtwo,lda,worklen;
  PetscBLASInt   m,n,i,j,k,l,s,nb,sz,lwork,info;
  PetscScalar    *tau,*work,*A=NULL,*QQ=NULL,*Qhalf,*C=NULL,one=1.0,zero=0.0;
  PetscMPIInt    rank,size,count,stride;
  MPI_Datatype   tmat;

  PetscFunctionBegin;
  ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m_,&m);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  k  = PetscMin(m,n);
  nb = 16;
  lda = 2*n;
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)bv),&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)bv),&rank);CHKERRQ(ierr);
  nlevels = (PetscInt)PetscCeilReal(PetscLog2Real((PetscReal)size));
  powtwo  = PetscPowInt(2,(PetscInt)PetscFloorReal(PetscLog2Real((PetscReal)size)));
  worklen = n+n*nb;
  if (nlevels) worklen += n*lda+n*lda*nlevels+n*lda;
  ierr = BVAllocateWork_Private(bv,worklen);CHKERRQ(ierr);
  tau  = bv->work;
  work = bv->work+n;
  ierr = PetscBLASIntCast(n*nb,&lwork);CHKERRQ(ierr);
  if (nlevels) {
    A  = bv->work+n+n*nb;
    QQ = bv->work+n+n*nb+n*lda;
    C  = bv->work+n+n*nb+n*lda+n*lda*nlevels;
  }

  /* Compute QR */
  PetscStackCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&m,&n,Q,&m,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("geqrf",info);

  /* Extract R */
  if (R || nlevels) {
    for (j=0;j<n;j++) {
      for (i=0;i<=PetscMin(j,m-1);i++) {
        if (nlevels) A[i+j*lda] = Q[i+j*m];
        else R[i+j*ldr] = Q[i+j*m];
      }
      for (i=PetscMin(j,m-1)+1;i<n;i++) {
        if (nlevels) A[i+j*lda] = 0.0;
        else R[i+j*ldr] = 0.0;
      }
    }
  }

  /* Compute orthogonal matrix in Q */
  PetscStackCallBLAS("LAPACKorgqr",LAPACKorgqr_(&m,&k,&k,Q,&m,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("orgqr",info);

  if (nlevels) {

    ierr = PetscMPIIntCast(n,&count);CHKERRQ(ierr);
    ierr = PetscMPIIntCast(lda,&stride);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(lda,&l);CHKERRQ(ierr);
    ierr = MPI_Type_vector(count,count,stride,MPIU_SCALAR,&tmat);CHKERRQ(ierr);
    ierr = MPI_Type_commit(&tmat);CHKERRQ(ierr);

    for (level=nlevels;level>=1;level--) {

      plevel = PetscPowInt(2,level);
      ierr = PetscBLASIntCast(plevel*PetscFloorReal(rank/(PetscReal)plevel)+(rank+PetscPowInt(2,level-1))%plevel,&s);CHKERRQ(ierr);

      /* Stack triangular matrices */
      if (rank<s && s<size) {  /* send top part, receive bottom part */
        ierr = MPI_Sendrecv(A,1,tmat,s,111,A+n,1,tmat,s,111,PetscObjectComm((PetscObject)bv),MPI_STATUS_IGNORE);CHKERRQ(ierr);
      } else if (s<size) {  /* copy top to bottom, receive top part */
        ierr = MPI_Sendrecv(A,1,tmat,rank,111,A+n,1,tmat,rank,111,PetscObjectComm((PetscObject)bv),MPI_STATUS_IGNORE);CHKERRQ(ierr);
        ierr = MPI_Sendrecv(A+n,1,tmat,s,111,A,1,tmat,s,111,PetscObjectComm((PetscObject)bv),MPI_STATUS_IGNORE);CHKERRQ(ierr);
      }
      if (level<nlevels && size!=powtwo) {  /* for cases when size is not a power of 2 */
        if (rank<size-powtwo) {  /* send bottom part */
          ierr = MPI_Send(A+n,1,tmat,rank+powtwo,111,PetscObjectComm((PetscObject)bv));CHKERRQ(ierr);
        } else if (rank>=powtwo) {  /* receive bottom part */
          ierr = MPI_Recv(A+n,1,tmat,rank-powtwo,111,PetscObjectComm((PetscObject)bv),MPI_STATUS_IGNORE);CHKERRQ(ierr);
        }
      }
      /* Compute QR and build orthogonal matrix */
      if (level<nlevels || (level==nlevels && s<size)) {
        PetscStackCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&l,&n,A,&l,tau,work,&lwork,&info));
        SlepcCheckLapackInfo("geqrf",info);
        ierr = PetscArraycpy(QQ+(level-1)*n*lda,A,n*lda);CHKERRQ(ierr);
        PetscStackCallBLAS("LAPACKorgqr",LAPACKorgqr_(&l,&n,&n,QQ+(level-1)*n*lda,&l,tau,work,&lwork,&info));
        SlepcCheckLapackInfo("orgqr",info);
        for (j=0;j<n;j++) {
          for (i=j+1;i<n;i++) A[i+j*lda] = 0.0;
        }
      } else if (level==nlevels) {  /* only one triangular matrix, set Q=I */
        ierr = PetscArrayzero(QQ+(level-1)*n*lda,n*lda);CHKERRQ(ierr);
        for (j=0;j<n;j++) QQ[j+j*lda+(level-1)*n*lda] = 1.0;
      }
    }

    /* Extract R */
    if (R) {
      for (j=0;j<n;j++) {
        for (i=0;i<=j;i++) R[i+j*ldr] = A[i+j*lda];
        for (i=j+1;i<n;i++) R[i+j*ldr] = 0.0;
      }
    }

    /* Accumulate orthogonal matrices */
    for (level=1;level<=nlevels;level++) {
      plevel = PetscPowInt(2,level);
      ierr = PetscBLASIntCast(plevel*PetscFloorReal(rank/(PetscReal)plevel)+(rank+PetscPowInt(2,level-1))%plevel,&s);CHKERRQ(ierr);
      Qhalf = (rank<s)? QQ+(level-1)*n*lda: QQ+(level-1)*n*lda+n;
      if (level<nlevels) {
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&l,&n,&n,&one,QQ+level*n*lda,&l,Qhalf,&l,&zero,C,&l));
        ierr = PetscArraycpy(QQ+level*n*lda,C,n*lda);CHKERRQ(ierr);
      } else {
        for (i=0;i<m/l;i++) {
          PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&l,&n,&n,&one,Q+i*l,&m,Qhalf,&l,&zero,C,&l));
          for (j=0;j<n;j++) { ierr = PetscArraycpy(Q+i*l+j*m,C+j*l,l);CHKERRQ(ierr); }
        }
        sz = m%l;
        if (sz) {
          PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&sz,&n,&n,&one,Q+(m/l)*l,&m,Qhalf,&l,&zero,C,&l));
          for (j=0;j<n;j++) { ierr = PetscArraycpy(Q+(m/l)*l+j*m,C+j*l,sz);CHKERRQ(ierr); }
        }
      }
    }

    ierr = MPI_Type_free(&tmat);CHKERRQ(ierr);
  }

  ierr = PetscLogFlops(3.0*m*n*n);CHKERRQ(ierr);
  ierr = PetscFPTrapPop();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    Reduction operation to compute [~,Rout]=qr([Rin1;Rin2]) in the TSQR algorithm;
    all matrices are upper triangular stored in packed format
*/
SLEPC_EXTERN void MPIAPI SlepcGivensPacked(void *in,void *inout,PetscMPIInt *len,MPI_Datatype *datatype)
{
  PetscErrorCode ierr;
  PetscBLASInt   n,i,j,k,one=1;
  PetscMPIInt    tsize;
  PetscScalar    v,s,*R2=(PetscScalar*)in,*R1=(PetscScalar*)inout;
  PetscReal      c;

  PetscFunctionBegin;
  ierr = MPI_Type_size(*datatype,&tsize);CHKERRABORT(PETSC_COMM_SELF,ierr);  /* we assume len=1 */
  tsize /= sizeof(PetscScalar);
  n = (-1+(PetscBLASInt)PetscSqrtReal(1+8*tsize))/2;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      LAPACKlartg_(R1+(2*n-j-1)*j/2+j,R2+(2*n-i-1)*i/2+j,&c,&s,&v);
      R1[(2*n-j-1)*j/2+j] = v;
      k = n-j-1;
      if (k) BLASrot_(&k,R1+(2*n-j-1)*j/2+j+1,&one,R2+(2*n-i-1)*i/2+j+1,&one,&c,&s);
    }
  }
  PetscFunctionReturnVoid();
}

/*
    Computes the R factor of the QR factorization of an mxn matrix via parallel TSQR
*/
PetscErrorCode BVOrthogonalize_LAPACK_TSQR_OnlyR(BV bv,PetscInt m_,PetscInt n_,PetscScalar *Q,PetscScalar *R,PetscInt ldr)
{
  PetscErrorCode ierr;
  PetscInt       worklen;
  PetscBLASInt   m,n,i,j,s,nb,lwork,info;
  PetscScalar    *tau,*work,*A=NULL,*R1=NULL,*R2=NULL;
  PetscMPIInt    size,count;
  MPI_Datatype   tmat;

  PetscFunctionBegin;
  ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m_,&m);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  nb = 16;
  s  = n+n*(n-1)/2;  /* length of packed triangular matrix */
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)bv),&size);CHKERRQ(ierr);
  worklen = n+n*nb+2*s+m*n;
  ierr = BVAllocateWork_Private(bv,worklen);CHKERRQ(ierr);
  tau  = bv->work;
  work = bv->work+n;
  R1   = bv->work+n+n*nb;
  R2   = bv->work+n+n*nb+s;
  A    = bv->work+n+n*nb+2*s;
  ierr = PetscBLASIntCast(n*nb,&lwork);CHKERRQ(ierr);
  ierr = PetscArraycpy(A,Q,m*n);CHKERRQ(ierr);

  /* Compute QR */
  PetscStackCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&m,&n,A,&m,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("geqrf",info);

  if (size==1) {
    /* Extract R */
    for (j=0;j<n;j++) {
      for (i=0;i<=PetscMin(j,m-1);i++) R[i+j*ldr] = A[i+j*m];
      for (i=PetscMin(j,m-1)+1;i<n;i++) R[i+j*ldr] = 0.0;
    }
  } else {
    /* Use MPI reduction operation to obtain global R */
    ierr = PetscMPIIntCast(s,&count);CHKERRQ(ierr);
    ierr = MPI_Type_contiguous(count,MPIU_SCALAR,&tmat);CHKERRQ(ierr);
    ierr = MPI_Type_commit(&tmat);CHKERRQ(ierr);
    for (i=0;i<n;i++) {
      for (j=i;j<n;j++) R1[(2*n-i-1)*i/2+j] = (i<m)?A[i+j*m]:0.0;
    }
    ierr = MPI_Allreduce(R1,R2,1,tmat,MPIU_TSQR,PetscObjectComm((PetscObject)bv));CHKERRQ(ierr);
    for (i=0;i<n;i++) {
      for (j=0;j<i;j++) R[i+j*ldr] = 0.0;
      for (j=i;j<n;j++) R[i+j*ldr] = R2[(2*n-i-1)*i/2+j];
    }
    ierr = MPI_Type_free(&tmat);CHKERRQ(ierr);
  }

  ierr = PetscLogFlops(3.0*m*n*n);CHKERRQ(ierr);
  ierr = PetscFPTrapPop();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

