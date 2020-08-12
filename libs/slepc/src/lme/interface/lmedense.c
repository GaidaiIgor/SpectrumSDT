/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Routines for solving dense matrix equations, in some cases calling SLICOT
*/

#include <slepc/private/lmeimpl.h>     /*I "slepclme.h" I*/
#include <slepcblaslapack.h>

/*
   LMEDenseRankSVD - given a square matrix A, compute its SVD U*S*V', and determine the
   numerical rank. On exit, U contains U*S and A is overwritten with V'
*/
PetscErrorCode LMEDenseRankSVD(LME lme,PetscInt n,PetscScalar *A,PetscInt lda,PetscScalar *U,PetscInt ldu,PetscInt *rank)
{
  PetscErrorCode ierr;
  PetscInt       i,j,rk=0;
  PetscScalar    *work;
  PetscReal      tol,*sg,*rwork;
  PetscBLASInt   n_,lda_,ldu_,info,lw_;

  PetscFunctionBegin;
  ierr = PetscCalloc3(n,&sg,10*n,&work,5*n,&rwork);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lda,&lda_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldu,&ldu_);CHKERRQ(ierr);
  lw_ = 10*n_;
#if !defined (PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","O",&n_,&n_,A,&lda_,sg,U,&ldu_,NULL,&n_,work,&lw_,&info));
#else
  PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","O",&n_,&n_,A,&lda_,sg,U,&ldu_,NULL,&n_,work,&lw_,rwork,&info));
#endif
  SlepcCheckLapackInfo("gesvd",info);
  tol = 10*PETSC_MACHINE_EPSILON*n*sg[0];
  for (j=0;j<n;j++) {
    if (sg[j]>tol) {
      for (i=0;i<n;i++) U[i+j*n] *= sg[j];
      rk++;
    } else break;
  }
  *rank = rk;
  ierr = PetscFree3(sg,work,rwork);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_USE_INFO)
/*
   LyapunovCholResidual - compute the residual norm ||A*U'*U+U'*U*A'+B*B'||
*/
static PetscErrorCode LyapunovCholResidual(PetscInt m,PetscScalar *A,PetscInt lda,PetscInt k,PetscScalar *B,PetscInt ldb,PetscScalar *U,PetscInt ldu,PetscReal *res)
{
  PetscErrorCode ierr;
  PetscBLASInt   n,kk,la,lb,lu;
  PetscScalar    *M,*R,zero=0.0,done=1.0;

  PetscFunctionBegin;
  *res = 0;
  ierr = PetscBLASIntCast(lda,&la);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldb,&lb);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldu,&lu);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&kk);CHKERRQ(ierr);
  ierr = PetscMalloc2(m*m,&M,m*m,&R);CHKERRQ(ierr);

  /* R = B*B' */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&kk,&done,B,&lb,B,&lb,&zero,R,&n));
  /* M = A*U' */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&n,&done,A,&la,U,&lu,&zero,M,&n));
  /* R = R+M*U */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&done,M,&n,U,&lu,&done,R,&n));
  /* R = R+U'*M' */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("C","C",&n,&n,&n,&done,U,&lu,M,&n,&done,R,&n));

  *res = LAPACKlange_("F",&n,&n,R,&n,NULL);
  ierr = PetscFree2(M,R);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   LyapunovResidual - compute the residual norm ||A*X+X*A'+B||
*/
static PetscErrorCode LyapunovResidual(PetscInt m,PetscScalar *A,PetscInt lda,PetscScalar *B,PetscInt ldb,PetscScalar *X,PetscInt ldx,PetscReal *res)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBLASInt   n,la,lb,lx;
  PetscScalar    *R,done=1.0;

  PetscFunctionBegin;
  *res = 0;
  ierr = PetscBLASIntCast(lda,&la);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldb,&lb);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldx,&lx);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ierr = PetscMalloc1(m*m,&R);CHKERRQ(ierr);

  /* R = B+A*X */
  for (i=0;i<m;i++) {
    ierr = PetscArraycpy(R+i*m,B+i*ldb,m);CHKERRQ(ierr);
  }
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&done,A,&la,X,&lx,&done,R,&n));
  /* R = R+X*A' */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&n,&done,X,&lx,A,&la,&done,R,&n));

  *res = LAPACKlange_("F",&n,&n,R,&n,NULL);
  ierr = PetscFree(R);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

#if defined(SLEPC_HAVE_SLICOT)
/*
   HessLyapunovChol_SLICOT - implementation used when SLICOT is available
*/
static PetscErrorCode HessLyapunovChol_SLICOT(PetscInt m,PetscScalar *H,PetscInt ldh,PetscInt k,PetscScalar *B,PetscInt ldb,PetscScalar *U,PetscInt ldu,PetscReal *res)
{
  PetscErrorCode ierr;
  PetscBLASInt   lwork,info,n,kk,lu,ione=1,sdim;
  PetscInt       i,j;
  PetscReal      scal;
  PetscScalar    *Q,*W,*wr,*wi,*work;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ldu,&lu);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&kk);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(6*m,&lwork);CHKERRQ(ierr);

  /* transpose W = H' */
  ierr = PetscMalloc5(m*m,&W,m*m,&Q,m,&wr,m,&wi,lwork,&work);CHKERRQ(ierr);
  for (j=0;j<m;j++) {
    for (i=0;i<m;i++) W[i+j*m] = H[j+i*ldh];
  }

  /* compute the real Schur form of W */
  PetscStackCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,W,&n,&sdim,wr,wi,Q,&n,work,&lwork,NULL,&info));
  SlepcCheckLapackInfo("gees",info);
#if defined(PETSC_USE_DEBUG)
  for (i=0;i<m;i++) if (PetscRealPart(wr[i])>=0.0) SETERRQ(PETSC_COMM_SELF,1,"Eigenvalue with non-negative real part, the coefficient matrix is not stable");
#endif

  /* copy B' into first rows of U */
  for (i=0;i<k;i++) {
    for (j=0;j<m;j++) U[i+j*ldu] = B[j+i*ldb];
  }

  /* solve Lyapunov equation (Hammarling) */
  PetscStackCallBLAS("SLICOTsb03od",SLICOTsb03od_("C","F","N",&n,&kk,W,&n,Q,&n,U,&lu,&scal,wr,wi,work,&lwork,&info));
  if (info) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in SLICOT subroutine SB03OD: info=%d",(int)info);
  if (scal!=1.0) SETERRQ1(PETSC_COMM_SELF,1,"Current implementation cannot handle scale factor %g",scal);

  /* resnorm = norm(H(m+1,:)*U'*U), use Q(:,1) = U'*U(:,m) */
  if (res) {
    for (j=0;j<m;j++) Q[j] = U[j+(m-1)*ldu];
    PetscStackCallBLAS("BLAStrmv",BLAStrmv_("U","C","N",&n,U,&lu,Q,&ione));
    if (k!=1) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"Residual error is intended for k=1 only, but you set k=%d",(int)k);
    *res *= BLASnrm2_(&n,Q,&ione);
  }

  ierr = PetscFree5(W,Q,wr,wi,work);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#else

/*
   Compute the upper Cholesky factor of A
 */
static PetscErrorCode CholeskyFactor(PetscInt m,PetscScalar *A,PetscInt lda)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    *S;
  PetscBLASInt   info,n,ld;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lda,&ld);CHKERRQ(ierr);
  ierr = PetscMalloc1(m*m,&S);CHKERRQ(ierr);

  /* save a copy of matrix in S */
  for (i=0;i<m;i++) {
    ierr = PetscArraycpy(S+i*m,A+i*lda,m);CHKERRQ(ierr);
  }

  /* compute upper Cholesky factor in R */
  PetscStackCallBLAS("LAPACKpotrf",LAPACKpotrf_("U",&n,A,&ld,&info));
  ierr = PetscLogFlops((1.0*n*n*n)/3.0);CHKERRQ(ierr);

  if (info) {
    ierr = PetscInfo(NULL,"potrf failed, retry on diagonally perturbed matrix\n");CHKERRQ(ierr);
    for (i=0;i<m;i++) {
      ierr = PetscArraycpy(A+i*lda,S+i*m,m);CHKERRQ(ierr);
      A[i+i*lda] += 50.0*PETSC_MACHINE_EPSILON;
    }
    PetscStackCallBLAS("LAPACKpotrf",LAPACKpotrf_("U",&n,A,&ld,&info));
    SlepcCheckLapackInfo("potrf",info);
    ierr = PetscLogFlops((1.0*n*n*n)/3.0);CHKERRQ(ierr);
  }

  /* Zero out entries below the diagonal */
  for (i=0;i<m-1;i++) {
    ierr = PetscArrayzero(A+i*lda+i+1,m-i-1);CHKERRQ(ierr);
  }
  ierr = PetscFree(S);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   HessLyapunovChol_LAPACK - alternative implementation when SLICOT is not available
*/
static PetscErrorCode HessLyapunovChol_LAPACK(PetscInt m,PetscScalar *H,PetscInt ldh,PetscInt k,PetscScalar *B,PetscInt ldb,PetscScalar *U,PetscInt ldu,PetscReal *res)
{
  PetscErrorCode ierr;
  PetscBLASInt   ilo=1,lwork,info,n,kk,lu,lb,ione=1;
  PetscInt       i,j;
  PetscReal      scal;
  PetscScalar    *Q,*C,*W,*Z,*wr,*work,zero=0.0,done=1.0,dmone=-1.0;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    *wi;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ldb,&lb);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldu,&lu);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&kk);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(6*m,&lwork);CHKERRQ(ierr);
  C = U;

#if !defined(PETSC_USE_COMPLEX)
  ierr = PetscMalloc6(m*m,&Q,m*m,&W,m*k,&Z,m,&wr,m,&wi,lwork,&work);CHKERRQ(ierr);
#else
  ierr = PetscMalloc5(m*m,&Q,m*m,&W,m*k,&Z,m,&wr,lwork,&work);CHKERRQ(ierr);
#endif

  /* save a copy W = H */
  for (j=0;j<m;j++) {
    for (i=0;i<m;i++) W[i+j*m] = H[i+j*ldh];
  }

  /* compute the (real) Schur form of W */
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKhseqr",LAPACKhseqr_("S","I",&n,&ilo,&n,W,&n,wr,wi,Q,&n,work,&lwork,&info));
#else
  PetscStackCallBLAS("LAPACKhseqr",LAPACKhseqr_("S","I",&n,&ilo,&n,W,&n,wr,Q,&n,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("hseqr",info);
#if defined(PETSC_USE_DEBUG)
  for (i=0;i<m;i++) if (PetscRealPart(wr[i])>=0.0) SETERRQ1(PETSC_COMM_SELF,1,"Eigenvalue with non-negative real part %g, the coefficient matrix is not stable",PetscRealPart(wr[i]));
#endif

  /* C = -Z*Z', Z = Q'*B */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&n,&kk,&n,&done,Q,&n,B,&lb,&zero,Z,&n));
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&kk,&dmone,Z,&n,Z,&n,&zero,C,&lu));

  /* solve triangular Sylvester equation */
  PetscStackCallBLAS("LAPACKtrsyl",LAPACKtrsyl_("N","C",&ione,&n,&n,W,&n,W,&n,C,&lu,&scal,&info));
  SlepcCheckLapackInfo("trsyl",info);
  if (scal!=1.0) SETERRQ1(PETSC_COMM_SELF,1,"Current implementation cannot handle scale factor %g",scal);

  /* back-transform C = Q*C*Q' */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&done,Q,&n,C,&n,&zero,W,&n));
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&n,&done,W,&n,Q,&n,&zero,C,&lu));

  /* resnorm = norm(H(m+1,:)*Y) */
  if (res) {
    if (k!=1) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"Residual error is intended for k=1 only, but you set k=%d",(int)k);
    *res *= BLASnrm2_(&n,C+m-1,&n);
  }

  /* U = chol(C) */
  ierr = CholeskyFactor(m,C,ldu);CHKERRQ(ierr);

#if !defined(PETSC_USE_COMPLEX)
  ierr = PetscFree6(Q,W,Z,wr,wi,work);CHKERRQ(ierr);
#else
  ierr = PetscFree5(Q,W,Z,wr,work);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

#endif /* SLEPC_HAVE_SLICOT */

/*@C
   LMEDenseHessLyapunovChol - Computes the Cholesky factor of the solution of a
   dense Lyapunov equation with an upper Hessenberg coefficient matrix.

   Logically Collective on lme

   Input Parameters:
+  lme - linear matrix equation solver context
.  m   - number of rows and columns of H
.  H   - coefficient matrix
.  ldh - leading dimension of H
.  k   - number of columns of B
.  B   - right-hand side matrix
.  ldb - leading dimension of B
-  ldu - leading dimension of U

   Output Parameter:
.  U   - Cholesky factor of the solution

   Input/Output Parameter:
.  res - (optional) residual norm, on input it should contain H(m+1,m)

   Note:
   The Lyapunov equation has the form H*X + X*H' = -B*B', where H is an mxm
   upper Hessenberg matrix, B is an mxk matrix and the solution is expressed
   as X = U'*U, where U is upper triangular. H is supposed to be stable.

   When k=1 and the res argument is provided, the last row of X is used to
   compute the residual norm of a Lyapunov equation projected via Arnoldi.

   Level: developer

.seealso: LMEDenseLyapunov(), LMESolve()
@*/
PetscErrorCode LMEDenseHessLyapunovChol(LME lme,PetscInt m,PetscScalar *H,PetscInt ldh,PetscInt k,PetscScalar *B,PetscInt ldb,PetscScalar *U,PetscInt ldu,PetscReal *res)
{
  PetscErrorCode ierr;
#if defined(PETSC_USE_INFO)
  PetscReal      error;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveInt(lme,m,2);
  PetscValidPointer(H,3);
  PetscValidLogicalCollectiveInt(lme,ldh,4);
  PetscValidLogicalCollectiveInt(lme,k,5);
  PetscValidPointer(B,6);
  PetscValidLogicalCollectiveInt(lme,ldb,7);
  PetscValidPointer(U,8);
  PetscValidLogicalCollectiveInt(lme,ldu,9);
  if (res) PetscValidLogicalCollectiveInt(lme,*res,10);

#if defined(SLEPC_HAVE_SLICOT)
  ierr = HessLyapunovChol_SLICOT(m,H,ldh,k,B,ldb,U,ldu,res);CHKERRQ(ierr);
#else
  ierr = HessLyapunovChol_LAPACK(m,H,ldh,k,B,ldb,U,ldu,res);CHKERRQ(ierr);
#endif

#if defined(PETSC_USE_INFO)
  if (PetscLogPrintInfo) {
    ierr = LyapunovCholResidual(m,H,ldh,k,B,ldb,U,ldu,&error);CHKERRQ(ierr);
    ierr = PetscInfo1(lme,"Residual norm of dense Lyapunov equation = %g\n",error);CHKERRQ(ierr);
  }
#endif
  PetscFunctionReturn(0);
}

#if defined(SLEPC_HAVE_SLICOT)
/*
   Lyapunov_SLICOT - implementation used when SLICOT is available
*/
static PetscErrorCode Lyapunov_SLICOT(PetscInt m,PetscScalar *H,PetscInt ldh,PetscScalar *B,PetscInt ldb,PetscScalar *X,PetscInt ldx)
{
  PetscErrorCode ierr;
  PetscBLASInt   sdim,lwork,info,n,lx,*iwork;
  PetscInt       i,j;
  PetscReal      scal,sep,ferr,*work;
  PetscScalar    *Q,*W,*wr,*wi;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ldx,&lx);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(PetscMax(20,m*m),&lwork);CHKERRQ(ierr);

  /* transpose W = H' */
  ierr = PetscMalloc6(m*m,&W,m*m,&Q,m,&wr,m,&wi,m*m,&iwork,lwork,&work);CHKERRQ(ierr);
  for (j=0;j<m;j++) {
    for (i=0;i<m;i++) W[i+j*m] = H[j+i*ldh];
  }

  /* compute the real Schur form of W */
  PetscStackCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,W,&n,&sdim,wr,wi,Q,&n,work,&lwork,NULL,&info));
  SlepcCheckLapackInfo("gees",info);

  /* copy -B into X */
  for (i=0;i<m;i++) {
    for (j=0;j<m;j++) X[i+j*ldx] = -B[i+j*ldb];
  }

  /* solve Lyapunov equation (Hammarling) */
  PetscStackCallBLAS("SLICOTsb03od",SLICOTsb03md_("C","X","F","N",&n,W,&n,Q,&n,X,&lx,&scal,&sep,&ferr,wr,wi,iwork,work,&lwork,&info));
  if (info) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in SLICOT subroutine SB03OD: info=%d",(int)info);
  if (scal!=1.0) SETERRQ1(PETSC_COMM_SELF,1,"Current implementation cannot handle scale factor %g",scal);

  ierr = PetscFree6(W,Q,wr,wi,iwork,work);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#else

/*
   Lyapunov_LAPACK - alternative implementation when SLICOT is not available
*/
static PetscErrorCode Lyapunov_LAPACK(PetscInt m,PetscScalar *A,PetscInt lda,PetscScalar *B,PetscInt ldb,PetscScalar *X,PetscInt ldx)
{
  PetscErrorCode ierr;
  PetscBLASInt   sdim,lwork,info,n,lx,lb,ione=1;
  PetscInt       i,j;
  PetscReal      scal;
  PetscScalar    *Q,*W,*Z,*wr,*work,zero=0.0,done=1.0,dmone=-1.0;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork;
#else
  PetscScalar    *wi;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ldb,&lb);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldx,&lx);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(6*m,&lwork);CHKERRQ(ierr);

#if !defined(PETSC_USE_COMPLEX)
  ierr = PetscMalloc6(m*m,&Q,m*m,&W,m*m,&Z,m,&wr,m,&wi,lwork,&work);CHKERRQ(ierr);
#else
  ierr = PetscMalloc6(m*m,&Q,m*m,&W,m*m,&Z,m,&wr,lwork,&work,m,&rwork);CHKERRQ(ierr);
#endif

  /* save a copy W = A */
  for (j=0;j<m;j++) {
    for (i=0;i<m;i++) W[i+j*m] = A[i+j*lda];
  }

  /* compute the (real) Schur form of W */
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,W,&n,&sdim,wr,wi,Q,&n,work,&lwork,NULL,&info));
#else
  PetscStackCallBLAS("LAPACKgees",LAPACKgees_("V","N",NULL,&n,W,&n,&sdim,wr,Q,&n,work,&lwork,rwork,NULL,&info));
#endif
  SlepcCheckLapackInfo("gees",info);

  /* X = -Q'*B*Q */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&n,&n,&n,&done,Q,&n,B,&lb,&zero,Z,&n));
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&dmone,Z,&n,Q,&n,&zero,X,&lx));

  /* solve triangular Sylvester equation */
  PetscStackCallBLAS("LAPACKtrsyl",LAPACKtrsyl_("N","C",&ione,&n,&n,W,&n,W,&n,X,&lx,&scal,&info));
  SlepcCheckLapackInfo("trsyl",info);
  if (scal!=1.0) SETERRQ1(PETSC_COMM_SELF,1,"Current implementation cannot handle scale factor %g",scal);

  /* back-transform X = Q*X*Q' */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&done,Q,&n,X,&n,&zero,W,&n));
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&n,&n,&n,&done,W,&n,Q,&n,&zero,X,&lx));

#if !defined(PETSC_USE_COMPLEX)
  ierr = PetscFree6(Q,W,Z,wr,wi,work);CHKERRQ(ierr);
#else
  ierr = PetscFree6(Q,W,Z,wr,work,rwork);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

#endif /* SLEPC_HAVE_SLICOT */

/*@C
   LMEDenseLyapunov - Computes the solution of a dense continuous-time Lyapunov
   equation.

   Logically Collective on lme

   Input Parameters:
+  lme - linear matrix equation solver context
.  m   - number of rows and columns of A
.  A   - coefficient matrix
.  lda - leading dimension of A
.  B   - right-hand side matrix
.  ldb - leading dimension of B
-  ldx - leading dimension of X

   Output Parameter:
.  X   - the solution

   Note:
   The Lyapunov equation has the form A*X + X*A' = -B, where all are mxm
   matrices, and B is symmetric.

   Level: developer

.seealso: LMEDenseHessLyapunovChol(), LMESolve()
@*/
PetscErrorCode LMEDenseLyapunov(LME lme,PetscInt m,PetscScalar *A,PetscInt lda,PetscScalar *B,PetscInt ldb,PetscScalar *X,PetscInt ldx)
{
  PetscErrorCode ierr;
#if defined(PETSC_USE_INFO)
  PetscReal      error;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveInt(lme,m,2);
  PetscValidPointer(A,3);
  PetscValidLogicalCollectiveInt(lme,lda,4);
  PetscValidPointer(B,5);
  PetscValidLogicalCollectiveInt(lme,ldb,6);
  PetscValidPointer(X,7);
  PetscValidLogicalCollectiveInt(lme,ldx,8);

#if defined(SLEPC_HAVE_SLICOT)
  ierr = Lyapunov_SLICOT(m,A,lda,B,ldb,X,ldx);CHKERRQ(ierr);
#else
  ierr = Lyapunov_LAPACK(m,A,lda,B,ldb,X,ldx);CHKERRQ(ierr);
#endif

#if defined(PETSC_USE_INFO)
  if (PetscLogPrintInfo) {
    ierr = LyapunovResidual(m,A,lda,B,ldb,X,ldx,&error);CHKERRQ(ierr);
    ierr = PetscInfo1(lme,"Residual norm of dense Lyapunov equation = %g\n",error);CHKERRQ(ierr);
  }
#endif
  PetscFunctionReturn(0);
}

