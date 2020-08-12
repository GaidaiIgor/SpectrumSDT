/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV implemented as a single Vec
*/

#include <slepc/private/bvimpl.h>
#include "svec.h"

PetscErrorCode BVMult_Svec(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  PetscErrorCode    ierr;
  BV_SVEC           *y = (BV_SVEC*)Y->data,*x = (BV_SVEC*)X->data;
  const PetscScalar *px;
  PetscScalar       *py,*q;
  PetscInt          ldq;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(x->v,&px);CHKERRQ(ierr);
  ierr = VecGetArray(y->v,&py);CHKERRQ(ierr);
  if (Q) {
    ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
    ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
    ierr = BVMult_BLAS_Private(Y,Y->n,Y->k-Y->l,X->k-X->l,ldq,alpha,px+(X->nc+X->l)*X->n,q+Y->l*ldq+X->l,beta,py+(Y->nc+Y->l)*Y->n);CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  } else {
    ierr = BVAXPY_BLAS_Private(Y,Y->n,Y->k-Y->l,alpha,px+(X->nc+X->l)*X->n,beta,py+(Y->nc+Y->l)*Y->n);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(x->v,&px);CHKERRQ(ierr);
  ierr = VecRestoreArray(y->v,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVMultVec_Svec(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar *q)
{
  PetscErrorCode ierr;
  BV_SVEC        *x = (BV_SVEC*)X->data;
  PetscScalar    *px,*py,*qq=q;

  PetscFunctionBegin;
  ierr = VecGetArray(x->v,&px);CHKERRQ(ierr);
  ierr = VecGetArray(y,&py);CHKERRQ(ierr);
  if (!q) { ierr = VecGetArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = BVMultVec_BLAS_Private(X,X->n,X->k-X->l,alpha,px+(X->nc+X->l)*X->n,qq,beta,py);CHKERRQ(ierr);
  if (!q) { ierr = VecRestoreArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = VecRestoreArray(x->v,&px);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVMultInPlace_Svec(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)V->data;
  PetscScalar    *pv,*q;
  PetscInt       ldq;

  PetscFunctionBegin;
  ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
  ierr = VecGetArray(ctx->v,&pv);CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  ierr = BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,ldq,s-V->l,e-V->l,pv+(V->nc+V->l)*V->n,q+V->l*ldq+V->l,PETSC_FALSE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  ierr = VecRestoreArray(ctx->v,&pv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVMultInPlaceTranspose_Svec(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)V->data;
  PetscScalar    *pv,*q;
  PetscInt       ldq;

  PetscFunctionBegin;
  ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
  ierr = VecGetArray(ctx->v,&pv);CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  ierr = BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,ldq,s-V->l,e-V->l,pv+(V->nc+V->l)*V->n,q+V->l*ldq+V->l,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  ierr = VecRestoreArray(ctx->v,&pv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDot_Svec(BV X,BV Y,Mat M)
{
  PetscErrorCode    ierr;
  BV_SVEC           *x = (BV_SVEC*)X->data,*y = (BV_SVEC*)Y->data;
  const PetscScalar *px,*py;
  PetscScalar       *m;
  PetscInt          ldm;

  PetscFunctionBegin;
  ierr = MatGetSize(M,&ldm,NULL);CHKERRQ(ierr);
  ierr = VecGetArrayRead(x->v,&px);CHKERRQ(ierr);
  ierr = VecGetArrayRead(y->v,&py);CHKERRQ(ierr);
  ierr = MatDenseGetArray(M,&m);CHKERRQ(ierr);
  ierr = BVDot_BLAS_Private(X,Y->k-Y->l,X->k-X->l,X->n,ldm,py+(Y->nc+Y->l)*Y->n,px+(X->nc+X->l)*X->n,m+X->l*ldm+Y->l,x->mpi);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(M,&m);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(x->v,&px);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(y->v,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_Svec(BV X,Vec y,PetscScalar *q)
{
  PetscErrorCode    ierr;
  BV_SVEC           *x = (BV_SVEC*)X->data;
  const PetscScalar *px,*py;
  PetscScalar       *qq=q;
  Vec               z = y;

  PetscFunctionBegin;
  if (X->matrix) {
    ierr = BV_IPMatMult(X,y);CHKERRQ(ierr);
    z = X->Bx;
  }
  ierr = VecGetArrayRead(x->v,&px);CHKERRQ(ierr);
  ierr = VecGetArrayRead(z,&py);CHKERRQ(ierr);
  if (!q) { ierr = VecGetArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = BVDotVec_BLAS_Private(X,X->n,X->k-X->l,px+(X->nc+X->l)*X->n,py,qq,x->mpi);CHKERRQ(ierr);
  if (!q) { ierr = VecRestoreArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = VecRestoreArrayRead(z,&py);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(x->v,&px);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_Local_Svec(BV X,Vec y,PetscScalar *m)
{
  PetscErrorCode ierr;
  BV_SVEC        *x = (BV_SVEC*)X->data;
  PetscScalar    *px,*py;
  Vec            z = y;

  PetscFunctionBegin;
  if (X->matrix) {
    ierr = BV_IPMatMult(X,y);CHKERRQ(ierr);
    z = X->Bx;
  }
  ierr = VecGetArray(x->v,&px);CHKERRQ(ierr);
  ierr = VecGetArray(z,&py);CHKERRQ(ierr);
  ierr = BVDotVec_BLAS_Private(X,X->n,X->k-X->l,px+(X->nc+X->l)*X->n,py,m,PETSC_FALSE);CHKERRQ(ierr);
  ierr = VecRestoreArray(z,&py);CHKERRQ(ierr);
  ierr = VecRestoreArray(x->v,&px);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVScale_Svec(BV bv,PetscInt j,PetscScalar alpha)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *array;

  PetscFunctionBegin;
  ierr = VecGetArray(ctx->v,&array);CHKERRQ(ierr);
  if (j<0) {
    ierr = BVScale_BLAS_Private(bv,(bv->k-bv->l)*bv->n,array+(bv->nc+bv->l)*bv->n,alpha);CHKERRQ(ierr);
  } else {
    ierr = BVScale_BLAS_Private(bv,bv->n,array+(bv->nc+j)*bv->n,alpha);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(ctx->v,&array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_Svec(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *array;

  PetscFunctionBegin;
  ierr = VecGetArray(ctx->v,&array);CHKERRQ(ierr);
  if (j<0) {
    ierr = BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->n,type,val,ctx->mpi);CHKERRQ(ierr);
  } else {
    ierr = BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->n,type,val,ctx->mpi);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(ctx->v,&array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_Local_Svec(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *array;

  PetscFunctionBegin;
  ierr = VecGetArray(ctx->v,&array);CHKERRQ(ierr);
  if (j<0) {
    ierr = BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->n,type,val,PETSC_FALSE);CHKERRQ(ierr);
  } else {
    ierr = BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->n,type,val,PETSC_FALSE);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(ctx->v,&array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVMatMult_Svec(BV V,Mat A,BV W)
{
  PetscErrorCode ierr;
  PetscInt       j;
  PetscBool      flg;
  Mat            Vmat,Wmat;
  Vec            vv,ww;

  PetscFunctionBegin;
  ierr = MatHasOperation(A,MATOP_MAT_MULT,&flg);CHKERRQ(ierr);
  if (V->vmm && flg) {
    ierr = BVGetMat(V,&Vmat);CHKERRQ(ierr);
    ierr = BVGetMat(W,&Wmat);CHKERRQ(ierr);
    ierr = MatProductCreateWithMat(A,Vmat,NULL,Wmat);CHKERRQ(ierr);
    ierr = MatProductSetType(Wmat,MATPRODUCT_AB);CHKERRQ(ierr);
    ierr = MatProductSetFromOptions(Wmat);CHKERRQ(ierr);
    ierr = MatProductSymbolic(Wmat);CHKERRQ(ierr);
    ierr = MatProductNumeric(Wmat);CHKERRQ(ierr);
    ierr = MatProductClear(Wmat);CHKERRQ(ierr);
    ierr = BVRestoreMat(V,&Vmat);CHKERRQ(ierr);
    ierr = BVRestoreMat(W,&Wmat);CHKERRQ(ierr);
  } else {
    for (j=0;j<V->k-V->l;j++) {
      ierr = BVGetColumn(V,V->l+j,&vv);CHKERRQ(ierr);
      ierr = BVGetColumn(W,W->l+j,&ww);CHKERRQ(ierr);
      ierr = MatMult(A,vv,ww);CHKERRQ(ierr);
      ierr = BVRestoreColumn(V,V->l+j,&vv);CHKERRQ(ierr);
      ierr = BVRestoreColumn(W,W->l+j,&ww);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVCopy_Svec(BV V,BV W)
{
  PetscErrorCode ierr;
  BV_SVEC        *v = (BV_SVEC*)V->data,*w = (BV_SVEC*)W->data;
  PetscScalar    *pv,*pw,*pvc,*pwc;

  PetscFunctionBegin;
  ierr = VecGetArray(v->v,&pv);CHKERRQ(ierr);
  ierr = VecGetArray(w->v,&pw);CHKERRQ(ierr);
  pvc = pv+(V->nc+V->l)*V->n;
  pwc = pw+(W->nc+W->l)*W->n;
  ierr = PetscArraycpy(pwc,pvc,(V->k-V->l)*V->n);CHKERRQ(ierr);
  ierr = VecRestoreArray(v->v,&pv);CHKERRQ(ierr);
  ierr = VecRestoreArray(w->v,&pw);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVCopyColumn_Svec(BV V,PetscInt j,PetscInt i)
{
  PetscErrorCode ierr;
  BV_SVEC        *v = (BV_SVEC*)V->data;
  PetscScalar    *pv;

  PetscFunctionBegin;
  ierr = VecGetArray(v->v,&pv);CHKERRQ(ierr);
  ierr = PetscArraycpy(pv+(V->nc+i)*V->n,pv+(V->nc+j)*V->n,V->n);CHKERRQ(ierr);
  ierr = VecRestoreArray(v->v,&pv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVResize_Svec(BV bv,PetscInt m,PetscBool copy)
{
  PetscErrorCode    ierr;
  BV_SVEC           *ctx = (BV_SVEC*)bv->data;
  PetscScalar       *pnew;
  const PetscScalar *pv;
  PetscInt          bs,lsplit;
  Vec               vnew,vpar;
  char              str[50];
  BV                parent;

  PetscFunctionBegin;
  if (bv->issplit==2) {
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    vpar = ((BV_SVEC*)parent->data)->v;
    ierr = VecGetArrayRead(vpar,&pv);CHKERRQ(ierr);
    ierr = VecResetArray(ctx->v);CHKERRQ(ierr);
    ierr = VecPlaceArray(ctx->v,pv+lsplit*bv->n);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(vpar,&pv);CHKERRQ(ierr);
  } else if (!bv->issplit) {
    ierr = VecGetBlockSize(bv->t,&bs);CHKERRQ(ierr);
    ierr = VecCreate(PetscObjectComm((PetscObject)bv->t),&vnew);CHKERRQ(ierr);
    ierr = VecSetType(vnew,((PetscObject)bv->t)->type_name);CHKERRQ(ierr);
    ierr = VecSetSizes(vnew,m*bv->n,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = VecSetBlockSize(vnew,bs);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)vnew);CHKERRQ(ierr);
    if (((PetscObject)bv)->name) {
      ierr = PetscSNPrintf(str,50,"%s_0",((PetscObject)bv)->name);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)vnew,str);CHKERRQ(ierr);
    }
    if (copy) {
      ierr = VecGetArrayRead(ctx->v,&pv);CHKERRQ(ierr);
      ierr = VecGetArray(vnew,&pnew);CHKERRQ(ierr);
      ierr = PetscArraycpy(pnew,pv,PetscMin(m,bv->m)*bv->n);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(ctx->v,&pv);CHKERRQ(ierr);
      ierr = VecRestoreArray(vnew,&pnew);CHKERRQ(ierr);
    }
    ierr = VecDestroy(&ctx->v);CHKERRQ(ierr);
    ctx->v = vnew;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetColumn_Svec(BV bv,PetscInt j,Vec *v)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscScalar    *pv;
  PetscInt       l;

  PetscFunctionBegin;
  l = BVAvailableVec;
  ierr = VecGetArray(ctx->v,&pv);CHKERRQ(ierr);
  ierr = VecPlaceArray(bv->cv[l],pv+(bv->nc+j)*bv->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreColumn_Svec(BV bv,PetscInt j,Vec *v)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;
  PetscInt       l;

  PetscFunctionBegin;
  l = (j==bv->ci[0])? 0: 1;
  ierr = VecResetArray(bv->cv[l]);CHKERRQ(ierr);
  ierr = VecRestoreArray(ctx->v,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetArray_Svec(BV bv,PetscScalar **a)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  ierr = VecGetArray(ctx->v,a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreArray_Svec(BV bv,PetscScalar **a)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  ierr = VecRestoreArray(ctx->v,a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetArrayRead_Svec(BV bv,const PetscScalar **a)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(ctx->v,a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreArrayRead_Svec(BV bv,const PetscScalar **a)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  ierr = VecRestoreArrayRead(ctx->v,a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVView_Svec(BV bv,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  BV_SVEC           *ctx = (BV_SVEC*)bv->data;
  PetscViewerFormat format;
  PetscBool         isascii;
  const char        *bvname,*name;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
    ierr = VecView(ctx->v,viewer);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      ierr = PetscObjectGetName((PetscObject)bv,&bvname);CHKERRQ(ierr);
      ierr = PetscObjectGetName((PetscObject)ctx->v,&name);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%s=reshape(%s,%D,%D);clear %s\n",bvname,name,bv->N,bv->nc+bv->m,name);CHKERRQ(ierr);
      if (bv->nc) {
        ierr = PetscViewerASCIIPrintf(viewer,"%s=%s(:,%D:end);\n",bvname,bvname,bv->nc+1);CHKERRQ(ierr);
      }
    }
  } else {
    ierr = VecView(ctx->v,viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVDestroy_Svec(BV bv)
{
  PetscErrorCode ierr;
  BV_SVEC        *ctx = (BV_SVEC*)bv->data;

  PetscFunctionBegin;
  ierr = VecDestroy(&ctx->v);CHKERRQ(ierr);
  ierr = VecDestroy(&bv->cv[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&bv->cv[1]);CHKERRQ(ierr);
  ierr = PetscFree(bv->data);CHKERRQ(ierr);
  bv->cuda = PETSC_FALSE;
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode BVCreate_Svec(BV bv)
{
  PetscErrorCode    ierr;
  BV_SVEC           *ctx;
  PetscInt          nloc,N,bs,tglobal=0,tlocal,lsplit;
  PetscBool         seq;
  PetscScalar       *aa,*vv;
  const PetscScalar *array,*ptr;
  char              str[50];
  BV                parent;
  Vec               vpar;
#if defined(PETSC_HAVE_CUDA)
  PetscScalar       *gpuarray,*gptr;
#endif

  PetscFunctionBegin;
  ierr = PetscNewLog(bv,&ctx);CHKERRQ(ierr);
  bv->data = (void*)ctx;

  ierr = PetscObjectTypeCompareAny((PetscObject)bv->t,&bv->cuda,VECSEQCUDA,VECMPICUDA,"");CHKERRQ(ierr);
  ierr = PetscObjectTypeCompareAny((PetscObject)bv->t,&ctx->mpi,VECMPI,VECMPICUDA,"");CHKERRQ(ierr);

  ierr = PetscObjectTypeCompare((PetscObject)bv->t,VECSEQ,&seq);CHKERRQ(ierr);
  if (!seq && !ctx->mpi && !bv->cuda) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"BVSVEC does not support the type of the provided template vector");

  ierr = VecGetLocalSize(bv->t,&nloc);CHKERRQ(ierr);
  ierr = VecGetSize(bv->t,&N);CHKERRQ(ierr);
  ierr = VecGetBlockSize(bv->t,&bs);CHKERRQ(ierr);
  tlocal  = bv->m*nloc;
  ierr = PetscIntMultError(bv->m,N,&tglobal);CHKERRQ(ierr);

  if (bv->issplit) {
    /* split BV: create Vec sharing the memory of the parent BV */
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    vpar = ((BV_SVEC*)parent->data)->v;
    if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
      ierr = VecCUDAGetArray(vpar,&gpuarray);CHKERRQ(ierr);
      gptr = (bv->issplit==1)? gpuarray: gpuarray+lsplit*nloc;
      ierr = VecCUDARestoreArray(vpar,&gpuarray);CHKERRQ(ierr);
      if (ctx->mpi) {
        ierr = VecCreateMPICUDAWithArray(PetscObjectComm((PetscObject)bv->t),bs,tlocal,bv->m*N,NULL,&ctx->v);CHKERRQ(ierr);
      } else {
        ierr = VecCreateSeqCUDAWithArray(PetscObjectComm((PetscObject)bv->t),bs,tlocal,NULL,&ctx->v);CHKERRQ(ierr);
      }
      ierr = VecCUDAPlaceArray(ctx->v,gptr);CHKERRQ(ierr);
#endif
    } else {
      ierr = VecGetArrayRead(vpar,&array);CHKERRQ(ierr);
      ptr = (bv->issplit==1)? array: array+lsplit*nloc;
      ierr = VecRestoreArrayRead(vpar,&array);CHKERRQ(ierr);
      if (ctx->mpi) {
        ierr = VecCreateMPIWithArray(PetscObjectComm((PetscObject)bv->t),bs,tlocal,bv->m*N,NULL,&ctx->v);CHKERRQ(ierr);
      } else {
        ierr = VecCreateSeqWithArray(PetscObjectComm((PetscObject)bv->t),bs,tlocal,NULL,&ctx->v);CHKERRQ(ierr);
      }
      ierr = VecPlaceArray(ctx->v,ptr);CHKERRQ(ierr);
    }
  } else {
    /* regular BV: create Vec to store the BV entries */
    ierr = VecCreate(PetscObjectComm((PetscObject)bv->t),&ctx->v);CHKERRQ(ierr);
    ierr = VecSetType(ctx->v,((PetscObject)bv->t)->type_name);CHKERRQ(ierr);
    ierr = VecSetSizes(ctx->v,tlocal,tglobal);CHKERRQ(ierr);
    ierr = VecSetBlockSize(ctx->v,bs);CHKERRQ(ierr);
  }
  ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)ctx->v);CHKERRQ(ierr);
  if (((PetscObject)bv)->name) {
    ierr = PetscSNPrintf(str,50,"%s_0",((PetscObject)bv)->name);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)ctx->v,str);CHKERRQ(ierr);
  }

  if (bv->Acreate) {
    ierr = MatDenseGetArray(bv->Acreate,&aa);CHKERRQ(ierr);
    ierr = VecGetArray(ctx->v,&vv);CHKERRQ(ierr);
    ierr = PetscArraycpy(vv,aa,tlocal);CHKERRQ(ierr);
    ierr = VecRestoreArray(ctx->v,&vv);CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(bv->Acreate,&aa);CHKERRQ(ierr);
    ierr = MatDestroy(&bv->Acreate);CHKERRQ(ierr);
  }

  ierr = VecDuplicateEmpty(bv->t,&bv->cv[0]);CHKERRQ(ierr);
  ierr = VecDuplicateEmpty(bv->t,&bv->cv[1]);CHKERRQ(ierr);

  if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
    bv->ops->mult             = BVMult_Svec_CUDA;
    bv->ops->multvec          = BVMultVec_Svec_CUDA;
    bv->ops->multinplace      = BVMultInPlace_Svec_CUDA;
    bv->ops->multinplacetrans = BVMultInPlaceTranspose_Svec_CUDA;
    bv->ops->dot              = BVDot_Svec_CUDA;
    bv->ops->dotvec           = BVDotVec_Svec_CUDA;
    bv->ops->dotvec_local     = BVDotVec_Local_Svec_CUDA;
    bv->ops->scale            = BVScale_Svec_CUDA;
    bv->ops->matmult          = BVMatMult_Svec_CUDA;
    bv->ops->copy             = BVCopy_Svec_CUDA;
    bv->ops->copycolumn       = BVCopyColumn_Svec_CUDA;
    bv->ops->resize           = BVResize_Svec_CUDA;
    bv->ops->getcolumn        = BVGetColumn_Svec_CUDA;
    bv->ops->restorecolumn    = BVRestoreColumn_Svec_CUDA;
    bv->ops->restoresplit     = BVRestoreSplit_Svec_CUDA;
#endif
  } else {
    bv->ops->mult             = BVMult_Svec;
    bv->ops->multvec          = BVMultVec_Svec;
    bv->ops->multinplace      = BVMultInPlace_Svec;
    bv->ops->multinplacetrans = BVMultInPlaceTranspose_Svec;
    bv->ops->dot              = BVDot_Svec;
    bv->ops->dotvec           = BVDotVec_Svec;
    bv->ops->dotvec_local     = BVDotVec_Local_Svec;
    bv->ops->scale            = BVScale_Svec;
    bv->ops->matmult          = BVMatMult_Svec;
    bv->ops->copy             = BVCopy_Svec;
    bv->ops->copycolumn       = BVCopyColumn_Svec;
    bv->ops->resize           = BVResize_Svec;
    bv->ops->getcolumn        = BVGetColumn_Svec;
    bv->ops->restorecolumn    = BVRestoreColumn_Svec;
  }
  bv->ops->norm             = BVNorm_Svec;
  bv->ops->norm_local       = BVNorm_Local_Svec;
  bv->ops->getarray         = BVGetArray_Svec;
  bv->ops->restorearray     = BVRestoreArray_Svec;
  bv->ops->getarrayread     = BVGetArrayRead_Svec;
  bv->ops->restorearrayread = BVRestoreArrayRead_Svec;
  bv->ops->destroy          = BVDestroy_Svec;
  if (!ctx->mpi) bv->ops->view = BVView_Svec;
  PetscFunctionReturn(0);
}

