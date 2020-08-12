/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Tensor BV that is represented in compact form as V = (I otimes U) S
*/

#include <slepc/private/bvimpl.h>      /*I "slepcbv.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  BV          U;        /* first factor */
  Mat         S;        /* second factor */
  PetscScalar *qB;      /* auxiliary matrix used in non-standard inner products */
  PetscScalar *sw;      /* work space */
  PetscInt    d;        /* degree of the tensor BV */
  PetscInt    ld;       /* leading dimension of a single block in S */
  PetscInt    puk;      /* copy of the k value */
  Vec         u;        /* auxiliary work vector */
} BV_TENSOR;

PetscErrorCode BVMultInPlace_Tensor(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscScalar    *pS,*q;
  PetscInt       ldq,lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(ctx->S,&pS);CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  ierr = BVMultInPlace_BLAS_Private(V,lds,V->k-V->l,ldq,s-V->l,e-V->l,pS+(V->nc+V->l)*lds,q+V->l*ldq+V->l,PETSC_FALSE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(ctx->S,&pS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVMultInPlaceTranspose_Tensor(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscScalar    *pS,*q;
  PetscInt       ldq,lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(ctx->S,&pS);CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  ierr = BVMultInPlace_BLAS_Private(V,lds,V->k-V->l,ldq,s-V->l,e-V->l,pS+(V->nc+V->l)*lds,q+V->l*ldq+V->l,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(ctx->S,&pS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDot_Tensor(BV X,BV Y,Mat M)
{
  PetscErrorCode ierr;
  BV_TENSOR      *x = (BV_TENSOR*)X->data,*y = (BV_TENSOR*)Y->data;
  PetscScalar    *m,*px,*py;
  PetscInt       ldm,lds = x->ld*x->d;

  PetscFunctionBegin;
  if (x->U!=y->U) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_SUP,"BVDot() in BVTENSOR requires that both operands have the same U factor");
  if (lds!=y->ld*y->d) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Mismatching dimensions ld*d %D %D",lds,y->ld*y->d);
  ierr = MatGetSize(M,&ldm,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(x->S,&px);CHKERRQ(ierr);
  ierr = MatDenseGetArray(y->S,&py);CHKERRQ(ierr);
  ierr = MatDenseGetArray(M,&m);CHKERRQ(ierr);
  ierr = BVDot_BLAS_Private(X,Y->k-Y->l,X->k-X->l,lds,ldm,py+(Y->nc+Y->l)*lds,px+(X->nc+X->l)*lds,m+X->l*ldm+Y->l,PETSC_FALSE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(M,&m);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(x->S,&px);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(y->S,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVScale_Tensor(BV bv,PetscInt j,PetscScalar alpha)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)bv->data;
  PetscScalar    *pS;
  PetscInt       lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(ctx->S,&pS);CHKERRQ(ierr);
  if (j<0) {
    ierr = BVScale_BLAS_Private(bv,(bv->k-bv->l)*lds,pS+(bv->nc+bv->l)*lds,alpha);CHKERRQ(ierr);
  } else {
    ierr = BVScale_BLAS_Private(bv,lds,pS+(bv->nc+j)*lds,alpha);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(ctx->S,&pS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_Tensor(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)bv->data;
  PetscScalar    *pS;
  PetscInt       lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(ctx->S,&pS);CHKERRQ(ierr);
  if (j<0) {
    ierr = BVNorm_LAPACK_Private(bv,lds,bv->k-bv->l,pS+(bv->nc+bv->l)*lds,type,val,PETSC_FALSE);CHKERRQ(ierr);
  } else {
    ierr = BVNorm_LAPACK_Private(bv,lds,1,pS+(bv->nc+j)*lds,type,val,PETSC_FALSE);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(ctx->S,&pS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVCopyColumn_Tensor(BV V,PetscInt j,PetscInt i)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscScalar    *pS;
  PetscInt       lds = ctx->ld*ctx->d;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(ctx->S,&pS);CHKERRQ(ierr);
  ierr = PetscArraycpy(pS+(V->nc+i)*lds,pS+(V->nc+j)*lds,lds);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(ctx->S,&pS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode BVTensorNormColumn(BV bv,PetscInt j,PetscReal *norm)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)bv->data;
  PetscBLASInt   one=1,lds_;
  PetscScalar    sone=1.0,szero=0.0,*S,*x,dot;
  PetscReal      alpha=1.0,scale=0.0,aval;
  PetscInt       i,lds,ld=ctx->ld;

  PetscFunctionBegin;
  lds = ld*ctx->d;
  ierr = MatDenseGetArray(ctx->S,&S);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  if (ctx->qB) {
    x = ctx->sw;
    PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&lds_,&sone,ctx->qB,&lds_,S+j*lds,&one,&szero,x,&one));
    dot = PetscRealPart(BLASdot_(&lds_,S+j*lds,&one,x,&one));
    ierr = BV_SafeSqrt(bv,dot,norm);CHKERRQ(ierr);
  } else {
    /* Compute *norm = BLASnrm2_(&lds_,S+j*lds,&one); */
    if (lds==1) *norm = PetscAbsScalar(S[j*lds]);
    else {
      for (i=0;i<lds;i++) {
        aval = PetscAbsScalar(S[i+j*lds]);
        if (aval!=0.0) {
          if (scale<aval) {
            alpha = 1.0 + alpha*PetscSqr(scale/aval);
            scale = aval;
          } else alpha += PetscSqr(aval/scale);
        }
      }
      *norm = scale*PetscSqrtReal(alpha);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVOrthogonalizeGS1_Tensor(BV bv,PetscInt k,Vec v,PetscBool *which,PetscScalar *h,PetscScalar *c,PetscReal *onorm,PetscReal *norm)
{
  PetscErrorCode    ierr;
  BV_TENSOR         *ctx = (BV_TENSOR*)bv->data;
  PetscScalar       *pS,*cc,*x,dot,sonem=-1.0,sone=1.0,szero=0.0;
  PetscInt          i,lds = ctx->ld*ctx->d;
  PetscBLASInt      lds_,k_,one=1;
  const PetscScalar *omega;

  PetscFunctionBegin;
  if (v) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Orthogonalization against an external vector is not allowed in BVTENSOR");
  ierr = MatDenseGetArray(ctx->S,&pS);CHKERRQ(ierr);
  if (!c) {
    ierr = VecGetArray(bv->buffer,&cc);CHKERRQ(ierr);
  } else cc = c;
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);

  if (onorm) { ierr = BVTensorNormColumn(bv,k,onorm);CHKERRQ(ierr); }

  if (ctx->qB) x = ctx->sw;
  else x = pS+k*lds;

  if (bv->orthog_type==BV_ORTHOG_MGS) {  /* modified Gram-Schmidt */

    if (bv->indef) { /* signature */
      ierr = VecGetArrayRead(bv->omega,&omega);CHKERRQ(ierr);
    }
    for (i=-bv->nc;i<k;i++) {
      if (which && i>=0 && !which[i]) continue;
      if (ctx->qB) PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&lds_,&sone,ctx->qB,&lds_,pS+k*lds,&one,&szero,x,&one));
      /* c_i = (s_k, s_i) */
      dot = PetscRealPart(BLASdot_(&lds_,pS+i*lds,&one,x,&one));
      if (bv->indef) dot /= PetscRealPart(omega[i]);
      ierr = BV_SetValue(bv,i,0,cc,dot);CHKERRQ(ierr);
      /* s_k = s_k - c_i s_i */
      dot = -dot;
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&lds_,&dot,pS+i*lds,&one,pS+k*lds,&one));
    }
    if (bv->indef) {
      ierr = VecRestoreArrayRead(bv->omega,&omega);CHKERRQ(ierr);
    }

  } else {  /* classical Gram-Schmidt */
    if (ctx->qB) PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&lds_,&sone,ctx->qB,&lds_,pS+k*lds,&one,&szero,x,&one));

    /* cc = S_{0:k-1}^* s_k */
    PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&lds_,&k_,&sone,pS,&lds_,x,&one,&szero,cc,&one));

    /* s_k = s_k - S_{0:k-1} cc */
    if (bv->indef) { ierr = BV_ApplySignature(bv,k,cc,PETSC_TRUE);CHKERRQ(ierr); }
     PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&k_,&sonem,pS,&lds_,cc,&one,&sone,pS+k*lds,&one));
    if (bv->indef) { ierr = BV_ApplySignature(bv,k,cc,PETSC_FALSE);CHKERRQ(ierr); }
  }

  if (norm) { ierr = BVTensorNormColumn(bv,k,norm);CHKERRQ(ierr); }
  ierr = BV_AddCoefficients(bv,k,h,cc);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(ctx->S,&pS);CHKERRQ(ierr);
  ierr = VecRestoreArray(bv->buffer,&cc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVView_Tensor(BV bv,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  BV_TENSOR         *ctx = (BV_TENSOR*)bv->data;
  PetscViewerFormat format;
  PetscBool         isascii;
  const char        *bvname,*uname,*sname;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) {
      ierr = PetscViewerASCIIPrintf(viewer,"number of tensor blocks (degree): %D\n",ctx->d);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"number of columns of U factor: %D\n",ctx->ld);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    ierr = BVView(ctx->U,viewer);CHKERRQ(ierr);
    ierr = MatView(ctx->S,viewer);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      ierr = PetscObjectGetName((PetscObject)bv,&bvname);CHKERRQ(ierr);
      ierr = PetscObjectGetName((PetscObject)ctx->U,&uname);CHKERRQ(ierr);
      ierr = PetscObjectGetName((PetscObject)ctx->S,&sname);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%s=kron(eye(%D),%s)*%s(:,1:%D);\n",bvname,ctx->d,uname,sname,bv->k);CHKERRQ(ierr);
    }
  } else {
    ierr = BVView(ctx->U,viewer);CHKERRQ(ierr);
    ierr = MatView(ctx->S,viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode BVTensorUpdateMatrix(BV V,PetscInt ini,PetscInt end)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscInt       i,j,r,c,l,k,ld=ctx->ld,lds=ctx->d*ctx->ld;
  PetscScalar    *qB,*sqB;
  Vec            u;
  Mat            A;

  PetscFunctionBegin;
  if (!V->matrix) PetscFunctionReturn(0);
  l = ctx->U->l; k = ctx->U->k;
  /* update inner product matrix */
  if (!ctx->qB) {
    ierr = PetscCalloc2(lds*lds,&ctx->qB,lds,&ctx->sw);CHKERRQ(ierr);
    ierr = VecDuplicate(ctx->U->t,&ctx->u);CHKERRQ(ierr);
  }
  ctx->U->l = 0;
  for (r=0;r<ctx->d;r++) {
    for (c=0;c<=r;c++) {
      ierr = MatNestGetSubMat(V->matrix,r,c,&A);CHKERRQ(ierr);
      if (A) {
        qB = ctx->qB+c*ld*lds+r*ld;
        for (i=ini;i<end;i++) {
          ierr = BVGetColumn(ctx->U,i,&u);CHKERRQ(ierr);
          ierr = MatMult(A,u,ctx->u);CHKERRQ(ierr);
          ctx->U->k = i+1;
          ierr = BVDotVec(ctx->U,ctx->u,qB+i*lds);CHKERRQ(ierr);
          ierr = BVRestoreColumn(ctx->U,i,&u);CHKERRQ(ierr);
          for (j=0;j<i;j++) qB[i+j*lds] = PetscConj(qB[j+i*lds]);
          qB[i*lds+i] = PetscRealPart(qB[i+i*lds]);
        }
        if (c!=r) {
          sqB = ctx->qB+r*ld*lds+c*ld;
          for (i=ini;i<end;i++) for (j=0;j<=i;j++) {
            sqB[i+j*lds] = PetscConj(qB[j+i*lds]);
            sqB[j+i*lds] = qB[j+i*lds];
          }
        }
      }
    }
  }
  ctx->U->l = l; ctx->U->k = k;
  PetscFunctionReturn(0);
}

static PetscErrorCode BVTensorBuildFirstColumn_Tensor(BV V,PetscInt k)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscInt       i,nq=0;
  PetscScalar    *pS,*omega;
  PetscReal      norm;
  PetscBool      breakdown=PETSC_FALSE;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(ctx->S,&pS);CHKERRQ(ierr);
  for (i=0;i<ctx->d;i++) {
    if (i>=k) {
      ierr = BVSetRandomColumn(ctx->U,nq);CHKERRQ(ierr);
    } else {
      ierr = BVCopyColumn(ctx->U,i,nq);CHKERRQ(ierr);
    }
    ierr = BVOrthogonalizeColumn(ctx->U,nq,pS+i*ctx->ld,&norm,&breakdown);CHKERRQ(ierr);
    if (!breakdown) {
      ierr = BVScaleColumn(ctx->U,nq,1.0/norm);CHKERRQ(ierr);
      pS[nq+i*ctx->ld] = norm;
      nq++;
    }
  }
  ierr = MatDenseRestoreArray(ctx->S,&pS);CHKERRQ(ierr);
  if (!nq) SETERRQ1(PetscObjectComm((PetscObject)V),1,"Cannot build first column of tensor BV; U should contain k=%D nonzero columns",k);
  ierr = BVTensorUpdateMatrix(V,0,nq);CHKERRQ(ierr);
  ierr = BVTensorNormColumn(V,0,&norm);CHKERRQ(ierr);
  ierr = BVScale_Tensor(V,0,1.0/norm);CHKERRQ(ierr);
  if (V->indef) {
    ierr = BV_AllocateSignature(V);CHKERRQ(ierr);
    ierr = VecGetArray(V->omega,&omega);CHKERRQ(ierr);
    omega[0] = (norm<0.0)? -1.0: 1.0;
    ierr = VecRestoreArray(V->omega,&omega);CHKERRQ(ierr);
  }
  /* set active columns */
  ctx->U->l = 0;
  ctx->U->k = nq;
  PetscFunctionReturn(0);
}

/*@
   BVTensorBuildFirstColumn - Builds the first column of the tensor basis vectors
   V from the data contained in the first k columns of U.

   Collective on V

   Input Parameters:
+  V - the basis vectors context
-  k - the number of columns of U with relevant information

   Notes:
   At most d columns are considered, where d is the degree of the tensor BV.
   Given V = (I otimes U) S, this function computes the first column of V, that
   is, it computes the coefficients of the first column of S by orthogonalizing
   the first d columns of U. If k is less than d (or linearly dependent columns
   are found) then additional random columns are used.

   The computed column has unit norm.

   Level: advanced

.seealso: BVCreateTensor()
@*/
PetscErrorCode BVTensorBuildFirstColumn(BV V,PetscInt k)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(V,k,2);
  ierr = PetscUseMethod(V,"BVTensorBuildFirstColumn_C",(BV,PetscInt),(V,k));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode BVTensorCompress_Tensor(BV V,PetscInt newc)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;
  PetscInt       nwu=0,nnc,nrow,lwa,r,c;
  PetscInt       i,j,k,n,lds=ctx->ld*ctx->d,deg=ctx->d,lock,cs1=V->k,rs1=ctx->U->k,rk=0,offu;
  PetscScalar    *S,*M,*Z,*pQ,*SS,*SS2,t,sone=1.0,zero=0.0,mone=-1.0,*p,*tau,*work,*qB,*sqB;
  PetscReal      *sg,tol,*rwork;
  PetscBLASInt   ld_,cs1_,rs1_,cs1tdeg,n_,info,lw_,newc_,newctdeg,nnc_,nrow_,nnctdeg,lds_,rk_;
  Mat            Q,A;

  PetscFunctionBegin;
  if (!cs1) PetscFunctionReturn(0);
  lwa = 6*ctx->ld*lds+2*cs1;
  n = PetscMin(rs1,deg*cs1);
  lock = ctx->U->l;
  nnc = cs1-lock-newc;
  nrow = rs1-lock;
  ierr = PetscCalloc6(deg*newc*nnc,&SS,newc*nnc,&SS2,(rs1+lock+newc)*n,&pQ,deg*rs1,&tau,lwa,&work,6*n,&rwork);CHKERRQ(ierr);
  offu = lock*(rs1+1);
  M = work+nwu;
  nwu += rs1*cs1*deg;
  sg = rwork;
  Z = work+nwu;
  nwu += deg*cs1*n;
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(nnc,&nnc_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(cs1,&cs1_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(rs1,&rs1_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(newc,&newc_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(newc*deg,&newctdeg);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(nnc*deg,&nnctdeg);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(cs1*deg,&cs1tdeg);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lwa-nwu,&lw_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(nrow,&nrow_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = MatDenseGetArray(ctx->S,&S);CHKERRQ(ierr);

  if (newc>0) {
    /* truncate columns associated with new converged eigenpairs */
    for (j=0;j<deg;j++) {
      for (i=lock;i<lock+newc;i++) {
        ierr = PetscArraycpy(M+(i-lock+j*newc)*nrow,S+i*lds+j*ctx->ld+lock,nrow);CHKERRQ(ierr);
      }
    }
#if !defined (PETSC_USE_COMPLEX)
    PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&nrow_,&newctdeg,M,&nrow_,sg,pQ+offu,&rs1_,Z,&n_,work+nwu,&lw_,&info));
#else
    PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&nrow_,&newctdeg,M,&nrow_,sg,pQ+offu,&rs1_,Z,&n_,work+nwu,&lw_,rwork+n,&info));
#endif
    SlepcCheckLapackInfo("gesvd",info);
    /* SVD has rank min(newc,nrow) */
    rk = PetscMin(newc,nrow);
    for (i=0;i<rk;i++) {
      t = sg[i];
      PetscStackCallBLAS("BLASscal",BLASscal_(&newctdeg,&t,Z+i,&n_));
    }
    for (i=0;i<deg;i++) {
      for (j=lock;j<lock+newc;j++) {
        ierr = PetscArraycpy(S+j*lds+i*ctx->ld+lock,Z+(newc*i+j-lock)*n,rk);CHKERRQ(ierr);
        ierr = PetscArrayzero(S+j*lds+i*ctx->ld+lock+rk,(ctx->ld-lock-rk));CHKERRQ(ierr);
      }
    }
    /*
      update columns associated with non-converged vectors, orthogonalize
      against pQ so that next M has rank nnc+d-1 insted of nrow+d-1
    */
    for (i=0;i<deg;i++) {
      PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&newc_,&nnc_,&nrow_,&sone,pQ+offu,&rs1_,S+(lock+newc)*lds+i*ctx->ld+lock,&lds_,&zero,SS+i*newc*nnc,&newc_));
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&nrow_,&nnc_,&newc_,&mone,pQ+offu,&rs1_,SS+i*newc*nnc,&newc_,&sone,S+(lock+newc)*lds+i*ctx->ld+lock,&lds_));
      /* repeat orthogonalization step */
      PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&newc_,&nnc_,&nrow_,&sone,pQ+offu,&rs1_,S+(lock+newc)*lds+i*ctx->ld+lock,&lds_,&zero,SS2,&newc_));
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&nrow_,&nnc_,&newc_,&mone,pQ+offu,&rs1_,SS2,&newc_,&sone,S+(lock+newc)*lds+i*ctx->ld+lock,&lds_));
      for (j=0;j<newc*nnc;j++) *(SS+i*newc*nnc+j) += SS2[j];
    }
  }

  /* truncate columns associated with non-converged eigenpairs */
  for (j=0;j<deg;j++) {
    for (i=lock+newc;i<cs1;i++) {
      ierr = PetscArraycpy(M+(i-lock-newc+j*nnc)*nrow,S+i*lds+j*ctx->ld+lock,nrow);CHKERRQ(ierr);
    }
  }
#if !defined (PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&nrow_,&nnctdeg,M,&nrow_,sg,pQ+offu+newc*rs1,&rs1_,Z,&n_,work+nwu,&lw_,&info));
#else
  PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&nrow_,&nnctdeg,M,&nrow_,sg,pQ+offu+newc*rs1,&rs1_,Z,&n_,work+nwu,&lw_,rwork+n,&info));
#endif
  SlepcCheckLapackInfo("gesvd",info);
  tol = PetscMax(rs1,deg*cs1)*PETSC_MACHINE_EPSILON*sg[0];
  rk = 0;
  for (i=0;i<PetscMin(nrow,nnctdeg);i++) if (sg[i]>tol) rk++;
  rk = PetscMin(nnc+deg-1,rk);
  /* the SVD has rank (at most) nnc+deg-1 */
  for (i=0;i<rk;i++) {
    t = sg[i];
    PetscStackCallBLAS("BLASscal",BLASscal_(&nnctdeg,&t,Z+i,&n_));
  }
  /* update S */
  ierr = PetscArrayzero(S+cs1*lds,(V->m-cs1)*lds);CHKERRQ(ierr);
  k = ctx->ld-lock-newc-rk;
  for (i=0;i<deg;i++) {
    for (j=lock+newc;j<cs1;j++) {
      ierr = PetscArraycpy(S+j*lds+i*ctx->ld+lock+newc,Z+(nnc*i+j-lock-newc)*n,rk);CHKERRQ(ierr);
      ierr = PetscArrayzero(S+j*lds+i*ctx->ld+lock+newc+rk,k);CHKERRQ(ierr);
    }
  }
  if (newc>0) {
    for (i=0;i<deg;i++) {
      p = SS+nnc*newc*i;
      for (j=lock+newc;j<cs1;j++) {
        for (k=0;k<newc;k++) S[j*lds+i*ctx->ld+lock+k] = *(p++);
      }
    }
  }

  /* orthogonalize pQ */
  rk = rk+newc;
  ierr = PetscBLASIntCast(rk,&rk_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(cs1-lock,&nnc_);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&nrow_,&rk_,pQ+offu,&rs1_,tau,work+nwu,&lw_,&info));
  SlepcCheckLapackInfo("geqrf",info);
  for (i=0;i<deg;i++) {
    PetscStackCallBLAS("BLAStrmm",BLAStrmm_("L","U","N","N",&rk_,&nnc_,&sone,pQ+offu,&rs1_,S+lock*lds+lock+i*ctx->ld,&lds_));
  }
  PetscStackCallBLAS("LAPACKorgqr",LAPACKorgqr_(&nrow_,&rk_,&rk_,pQ+offu,&rs1_,tau,work+nwu,&lw_,&info));
  SlepcCheckLapackInfo("orgqr",info);

  /* update vectors U(:,idx) = U*Q(:,idx) */
  rk = rk+lock;
  for (i=0;i<lock;i++) pQ[i*(1+rs1)] = 1.0;
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,rs1,rk,pQ,&Q);CHKERRQ(ierr);
  ctx->U->k = rs1;
  ierr = BVMultInPlace(ctx->U,Q,lock,rk);CHKERRQ(ierr);
  ierr = MatDestroy(&Q);CHKERRQ(ierr);

  if (ctx->qB) {
   /* update matrix qB */
    ierr = PetscBLASIntCast(ctx->ld,&ld_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(rk,&rk_);CHKERRQ(ierr);
    for (r=0;r<ctx->d;r++) {
      for (c=0;c<=r;c++) {
        ierr = MatNestGetSubMat(V->matrix,r,c,&A);CHKERRQ(ierr);
        if (A) {
          qB = ctx->qB+r*ctx->ld+c*ctx->ld*lds;
          PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&rs1_,&rk_,&rs1_,&sone,qB,&lds_,pQ,&rs1_,&zero,work+nwu,&rs1_));
          PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&rk_,&rk_,&rs1_,&sone,pQ,&rs1_,work+nwu,&rs1_,&zero,qB,&lds_));
          for (i=0;i<rk;i++) {
            for (j=0;j<i;j++) qB[i+j*lds] = PetscConj(qB[j+i*lds]);
            qB[i+i*lds] = PetscRealPart(qB[i+i*lds]);
          }
          for (i=rk;i<ctx->ld;i++) {
            ierr = PetscArrayzero(qB+i*lds,ctx->ld);CHKERRQ(ierr);
          }
          for (i=0;i<rk;i++) {
            ierr = PetscArrayzero(qB+i*lds+rk,(ctx->ld-rk));CHKERRQ(ierr);
          }
          if (c!=r) {
            sqB = ctx->qB+r*ctx->ld*lds+c*ctx->ld;
            for (i=0;i<ctx->ld;i++) for (j=0;j<ctx->ld;j++) sqB[i+j*lds] = PetscConj(qB[j+i*lds]);
          }
        }
      }
    }
  }

  /* free work space */
  ierr = PetscFree6(SS,SS2,pQ,tau,work,rwork);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(ctx->S,&S);CHKERRQ(ierr);

  /* set active columns */
  if (newc) ctx->U->l += newc;
  ctx->U->k = rk;
  PetscFunctionReturn(0);
}

/*@
   BVTensorCompress - Updates the U and S factors of the tensor basis vectors
   object V by means of an SVD, removing redundant information.

   Collective on V

   Input Parameters:
+  V - the tensor basis vectors context
-  newc - additional columns to be locked

   Notes:
   This function is typically used when restarting Krylov solvers. Truncating a
   tensor BV V = (I otimes U) S to its leading columns amounts to keeping the
   leading columns of S. However, to effectively reduce the size of the
   decomposition, it is necessary to compress it in a way that fewer columns of
   U are employed. This can be achieved by means of an update that involves the
   SVD of the low-rank matrix [S_0 S_1 ... S_{d-1}], where S_i are the pieces of S.

   If newc is nonzero, then newc columns are added to the leading columns of V.
   This means that the corresponding columns of the U and S factors will remain
   invariant in subsequent operations.

   Level: advanced

.seealso: BVCreateTensor(), BVSetActiveColumns()
@*/
PetscErrorCode BVTensorCompress(BV V,PetscInt newc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(V,newc,2);
  ierr = PetscUseMethod(V,"BVTensorCompress_C",(BV,PetscInt),(V,newc));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode BVTensorGetDegree_Tensor(BV bv,PetscInt *d)
{
  BV_TENSOR *ctx = (BV_TENSOR*)bv->data;

  PetscFunctionBegin;
  *d = ctx->d;
  PetscFunctionReturn(0);
}

/*@
   BVTensorGetDegree - Returns the number of blocks (degree) of the tensor BV.

   Not collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  d - the degree

   Level: advanced

.seealso: BVCreateTensor()
@*/
PetscErrorCode BVTensorGetDegree(BV bv,PetscInt *d)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidIntPointer(d,2);
  ierr = PetscUseMethod(bv,"BVTensorGetDegree_C",(BV,PetscInt*),(bv,d));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode BVTensorGetFactors_Tensor(BV V,BV *U,Mat *S)
{
  BV_TENSOR *ctx = (BV_TENSOR*)V->data;

  PetscFunctionBegin;
  if (ctx->puk>-1) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ORDER,"Previous call to BVTensonGetFactors without a BVTensorRestoreFactors call");
  ctx->puk = ctx->U->k;
  if (U) *U = ctx->U;
  if (S) *S = ctx->S;
  PetscFunctionReturn(0);
}

/*@C
   BVTensorGetFactors - Returns the two factors involved in the definition of the
   tensor basis vectors object, V = (I otimes U) S.

   Logically Collective on V

   Input Parameter:
.  V - the basis vectors context

   Output Parameters:
+  U - the BV factor
-  S - the Mat factor

   Notes:
   The returned factors are references (not copies) of the internal factors,
   so modifying them will change the tensor BV as well. Some operations of the
   tensor BV assume that U has orthonormal columns, so if the user modifies U
   this restriction must be taken into account.

   The returned factors must not be destroyed. BVTensorRestoreFactors() must
   be called when they are no longer needed.

   Pass a NULL vector for any of the arguments that is not needed.

   Level: advanced

.seealso: BVTensorRestoreFactors()
@*/
PetscErrorCode BVTensorGetFactors(BV V,BV *U,Mat *S)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  ierr = PetscUseMethod(V,"BVTensorGetFactors_C",(BV,BV*,Mat*),(V,U,S));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode BVTensorRestoreFactors_Tensor(BV V,BV *U,Mat *S)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)V->data;

  PetscFunctionBegin;
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  if (U) *U = NULL;
  if (S) *S = NULL;
  ierr = BVTensorUpdateMatrix(V,ctx->puk,ctx->U->k);CHKERRQ(ierr);
  ctx->puk = -1;
  PetscFunctionReturn(0);
}

/*@C
   BVTensorRestoreFactors - Restore the two factors that were obtained with
   BVTensorGetFactors().

   Logically Collective on V

   Input Parameters:
+  V - the basis vectors context
.  U - the BV factor (or NULL)
-  S - the Mat factor (or NULL)

   Nots:
   The arguments must match the corresponding call to BVTensorGetFactors().

   Level: advanced

.seealso: BVTensorGetFactors()
@*/
PetscErrorCode BVTensorRestoreFactors(BV V,BV *U,Mat *S)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  if (U) PetscValidHeaderSpecific(*U,BV_CLASSID,2);
  if (S) PetscValidHeaderSpecific(*S,MAT_CLASSID,3);
  ierr = PetscUseMethod(V,"BVTensorRestoreFactors_C",(BV,BV*,Mat*),(V,U,S));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDestroy_Tensor(BV bv)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx = (BV_TENSOR*)bv->data;

  PetscFunctionBegin;
  ierr = BVDestroy(&ctx->U);CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->S);CHKERRQ(ierr);
  if (ctx->u) {
    ierr = PetscFree2(ctx->qB,ctx->sw);CHKERRQ(ierr);
    ierr = VecDestroy(&ctx->u);CHKERRQ(ierr);
  }
  ierr = PetscFree(bv->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorBuildFirstColumn_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorCompress_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorGetDegree_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorGetFactors_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorRestoreFactors_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode BVCreate_Tensor(BV bv)
{
  PetscErrorCode ierr;
  BV_TENSOR      *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(bv,&ctx);CHKERRQ(ierr);
  bv->data = (void*)ctx;
  ctx->puk = -1;

  bv->ops->multinplace      = BVMultInPlace_Tensor;
  bv->ops->multinplacetrans = BVMultInPlaceTranspose_Tensor;
  bv->ops->dot              = BVDot_Tensor;
  bv->ops->scale            = BVScale_Tensor;
  bv->ops->norm             = BVNorm_Tensor;
  bv->ops->copycolumn       = BVCopyColumn_Tensor;
  bv->ops->gramschmidt      = BVOrthogonalizeGS1_Tensor;
  bv->ops->destroy          = BVDestroy_Tensor;
  bv->ops->view             = BVView_Tensor;

  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorBuildFirstColumn_C",BVTensorBuildFirstColumn_Tensor);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorCompress_C",BVTensorCompress_Tensor);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorGetDegree_C",BVTensorGetDegree_Tensor);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorGetFactors_C",BVTensorGetFactors_Tensor);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)bv,"BVTensorRestoreFactors_C",BVTensorRestoreFactors_Tensor);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVCreateTensor - Creates a tensor BV that is represented in compact form
   as V = (I otimes U) S, where U has orthonormal columns.

   Collective on U

   Input Parameters:
+  U - a basis vectors object
-  d - the number of blocks (degree) of the tensor BV

   Output Parameter:
.  V - the new basis vectors context

   Notes:
   The new basis vectors object is V = (I otimes U) S, where otimes denotes
   the Kronecker product, I is the identity matrix of order d, and S is a
   sequential matrix allocated internally. This compact representation is
   used e.g. to represent the Krylov basis generated with the linearization
   of a matrix polynomial of degree d.

   The size of V (number of rows) is equal to d times n, where n is the size
   of U. The dimensions of S are d times m rows and m-d+1 columns, where m is
   the number of columns of U, so m should be at least d.

   The communicator of V will be the same as U.

   On input, the content of U is irrelevant. Alternatively, it may contain
   some nonzero columns that will be used by BVTensorBuildFirstColumn().

   Level: advanced

.seealso: BVTensorGetDegree(), BVTensorGetFactors(), BVTensorBuildFirstColumn()
@*/
PetscErrorCode BVCreateTensor(BV U,PetscInt d,BV *V)
{
  PetscErrorCode ierr;
  PetscBool      match;
  PetscInt       n,N,m;
  BV_TENSOR      *ctx;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(U,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(U,d,2);
  ierr = PetscObjectTypeCompare((PetscObject)U,BVTENSOR,&match);CHKERRQ(ierr);
  if (match) SETERRQ(PetscObjectComm((PetscObject)U),PETSC_ERR_SUP,"U cannot be of type tensor");

  ierr = BVCreate(PetscObjectComm((PetscObject)U),V);CHKERRQ(ierr);
  ierr = BVGetSizes(U,&n,&N,&m);CHKERRQ(ierr);
  if (m<d) SETERRQ2(PetscObjectComm((PetscObject)U),PETSC_ERR_ARG_SIZ,"U has %D columns, it should have at least d=%D",m,d);
  ierr = BVSetSizes(*V,d*n,d*N,m-d+1);CHKERRQ(ierr);
  ierr = PetscObjectChangeTypeName((PetscObject)*V,BVTENSOR);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(BV_Create,*V,0,0,0);CHKERRQ(ierr);
  ierr = BVCreate_Tensor(*V);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_Create,*V,0,0,0);CHKERRQ(ierr);

  ctx = (BV_TENSOR*)(*V)->data;
  ctx->U  = U;
  ctx->d  = d;
  ctx->ld = m;
  ierr = PetscObjectReference((PetscObject)U);CHKERRQ(ierr);
  ierr = PetscLogObjectParent((PetscObject)*V,(PetscObject)U);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,d*m,m-d+1,NULL,&ctx->S);CHKERRQ(ierr);
  ierr = PetscLogObjectParent((PetscObject)*V,(PetscObject)ctx->S);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)ctx->S,"S");CHKERRQ(ierr);

  /* Copy user-provided attributes of U */
  (*V)->orthog_type  = U->orthog_type;
  (*V)->orthog_ref   = U->orthog_ref;
  (*V)->orthog_eta   = U->orthog_eta;
  (*V)->orthog_block = U->orthog_block;
  (*V)->vmm          = U->vmm;
  (*V)->rrandom      = U->rrandom;
  PetscFunctionReturn(0);
}

