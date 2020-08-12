/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV implemented with a dense Mat
*/

#include <slepc/private/bvimpl.h>

typedef struct {
  Mat       A;
  PetscBool mpi;
} BV_MAT;

PetscErrorCode BVMult_Mat(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  PetscErrorCode ierr;
  BV_MAT         *y = (BV_MAT*)Y->data,*x = (BV_MAT*)X->data;
  PetscScalar    *px,*py,*q;
  PetscInt       ldq;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(x->A,&px);CHKERRQ(ierr);
  ierr = MatDenseGetArray(y->A,&py);CHKERRQ(ierr);
  if (Q) {
    ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
    ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
    ierr = BVMult_BLAS_Private(Y,Y->n,Y->k-Y->l,X->k-X->l,ldq,alpha,px+(X->nc+X->l)*X->n,q+Y->l*ldq+X->l,beta,py+(Y->nc+Y->l)*Y->n);CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  } else {
    ierr = BVAXPY_BLAS_Private(Y,Y->n,Y->k-Y->l,alpha,px+(X->nc+X->l)*X->n,beta,py+(Y->nc+Y->l)*Y->n);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(x->A,&px);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(y->A,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVMultVec_Mat(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar *q)
{
  PetscErrorCode ierr;
  BV_MAT         *x = (BV_MAT*)X->data;
  PetscScalar    *px,*py,*qq=q;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(x->A,&px);CHKERRQ(ierr);
  ierr = VecGetArray(y,&py);CHKERRQ(ierr);
  if (!q) { ierr = VecGetArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = BVMultVec_BLAS_Private(X,X->n,X->k-X->l,alpha,px+(X->nc+X->l)*X->n,qq,beta,py);CHKERRQ(ierr);
  if (!q) { ierr = VecRestoreArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = MatDenseRestoreArray(x->A,&px);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVMultInPlace_Mat(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)V->data;
  PetscScalar    *pv,*q;
  PetscInt       ldq;

  PetscFunctionBegin;
  ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(ctx->A,&pv);CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  ierr = BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,ldq,s-V->l,e-V->l,pv+(V->nc+V->l)*V->n,q+V->l*ldq+V->l,PETSC_FALSE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(ctx->A,&pv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVMultInPlaceTranspose_Mat(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)V->data;
  PetscScalar    *pv,*q;
  PetscInt       ldq;

  PetscFunctionBegin;
  ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(ctx->A,&pv);CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  ierr = BVMultInPlace_BLAS_Private(V,V->n,V->k-V->l,ldq,s-V->l,e-V->l,pv+(V->nc+V->l)*V->n,q+V->l*ldq+V->l,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(ctx->A,&pv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDot_Mat(BV X,BV Y,Mat M)
{
  PetscErrorCode ierr;
  BV_MAT         *x = (BV_MAT*)X->data,*y = (BV_MAT*)Y->data;
  PetscScalar    *px,*py,*m;
  PetscInt       ldm;

  PetscFunctionBegin;
  ierr = MatGetSize(M,&ldm,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(x->A,&px);CHKERRQ(ierr);
  ierr = MatDenseGetArray(y->A,&py);CHKERRQ(ierr);
  ierr = MatDenseGetArray(M,&m);CHKERRQ(ierr);
  ierr = BVDot_BLAS_Private(X,Y->k-Y->l,X->k-X->l,X->n,ldm,py+(Y->nc+Y->l)*Y->n,px+(X->nc+X->l)*X->n,m+X->l*ldm+Y->l,x->mpi);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(M,&m);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(x->A,&px);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(y->A,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_Mat(BV X,Vec y,PetscScalar *q)
{
  PetscErrorCode    ierr;
  BV_MAT            *x = (BV_MAT*)X->data;
  PetscScalar       *px,*qq=q;
  const PetscScalar *py;
  Vec               z = y;

  PetscFunctionBegin;
  if (X->matrix) {
    ierr = BV_IPMatMult(X,y);CHKERRQ(ierr);
    z = X->Bx;
  }
  ierr = MatDenseGetArray(x->A,&px);CHKERRQ(ierr);
  ierr = VecGetArrayRead(z,&py);CHKERRQ(ierr);
  if (!q) { ierr = VecGetArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = BVDotVec_BLAS_Private(X,X->n,X->k-X->l,px+(X->nc+X->l)*X->n,py,qq,x->mpi);CHKERRQ(ierr);
  if (!q) { ierr = VecRestoreArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = VecRestoreArrayRead(z,&py);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(x->A,&px);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_Local_Mat(BV X,Vec y,PetscScalar *m)
{
  PetscErrorCode ierr;
  BV_MAT         *x = (BV_MAT*)X->data;
  PetscScalar    *px,*py;
  Vec            z = y;

  PetscFunctionBegin;
  if (X->matrix) {
    ierr = BV_IPMatMult(X,y);CHKERRQ(ierr);
    z = X->Bx;
  }
  ierr = MatDenseGetArray(x->A,&px);CHKERRQ(ierr);
  ierr = VecGetArray(z,&py);CHKERRQ(ierr);
  ierr = BVDotVec_BLAS_Private(X,X->n,X->k-X->l,px+(X->nc+X->l)*X->n,py,m,PETSC_FALSE);CHKERRQ(ierr);
  ierr = VecRestoreArray(z,&py);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(x->A,&px);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVScale_Mat(BV bv,PetscInt j,PetscScalar alpha)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *array;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(ctx->A,&array);CHKERRQ(ierr);
  if (j<0) {
    ierr = BVScale_BLAS_Private(bv,(bv->k-bv->l)*bv->n,array+(bv->nc+bv->l)*bv->n,alpha);CHKERRQ(ierr);
  } else {
    ierr = BVScale_BLAS_Private(bv,bv->n,array+(bv->nc+j)*bv->n,alpha);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(ctx->A,&array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_Mat(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *array;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(ctx->A,&array);CHKERRQ(ierr);
  if (j<0) {
    ierr = BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->n,type,val,ctx->mpi);CHKERRQ(ierr);
  } else {
    ierr = BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->n,type,val,ctx->mpi);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(ctx->A,&array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_Local_Mat(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *array;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(ctx->A,&array);CHKERRQ(ierr);
  if (j<0) {
    ierr = BVNorm_LAPACK_Private(bv,bv->n,bv->k-bv->l,array+(bv->nc+bv->l)*bv->n,type,val,PETSC_FALSE);CHKERRQ(ierr);
  } else {
    ierr = BVNorm_LAPACK_Private(bv,bv->n,1,array+(bv->nc+j)*bv->n,type,val,PETSC_FALSE);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(ctx->A,&array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVMatMult_Mat(BV V,Mat A,BV W)
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

PetscErrorCode BVCopy_Mat(BV V,BV W)
{
  PetscErrorCode ierr;
  BV_MAT         *v = (BV_MAT*)V->data,*w = (BV_MAT*)W->data;
  PetscScalar    *pv,*pw,*pvc,*pwc;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(v->A,&pv);CHKERRQ(ierr);
  ierr = MatDenseGetArray(w->A,&pw);CHKERRQ(ierr);
  pvc = pv+(V->nc+V->l)*V->n;
  pwc = pw+(W->nc+W->l)*W->n;
  ierr = PetscArraycpy(pwc,pvc,(V->k-V->l)*V->n);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(v->A,&pv);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(w->A,&pw);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVCopyColumn_Mat(BV V,PetscInt j,PetscInt i)
{
  PetscErrorCode ierr;
  BV_MAT         *v = (BV_MAT*)V->data;
  PetscScalar    *pv;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(v->A,&pv);CHKERRQ(ierr);
  ierr = PetscArraycpy(pv+(V->nc+i)*V->n,pv+(V->nc+j)*V->n,V->n);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(v->A,&pv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVResize_Mat(BV bv,PetscInt m,PetscBool copy)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *pA,*pnew;
  PetscInt       lsplit;
  Mat            A;
  char           str[50];
  PetscScalar    *array;
  BV             parent;

  PetscFunctionBegin;
  if (bv->issplit==2) {
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    ierr = MatDenseGetArray(((BV_MAT*)parent->data)->A,&array);CHKERRQ(ierr);
    ierr = MatDensePlaceArray(ctx->A,array+lsplit*bv->n);CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(((BV_MAT*)parent->data)->A,&array);CHKERRQ(ierr);
  } else if (!bv->issplit) {
    ierr = MatCreateDense(PetscObjectComm((PetscObject)bv->t),bv->n,PETSC_DECIDE,PETSC_DECIDE,m,NULL,&A);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)A);CHKERRQ(ierr);
    if (((PetscObject)bv)->name) {
      ierr = PetscSNPrintf(str,50,"%s_0",((PetscObject)bv)->name);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)A,str);CHKERRQ(ierr);
    }
    if (copy) {
      ierr = MatDenseGetArray(ctx->A,&pA);CHKERRQ(ierr);
      ierr = MatDenseGetArray(A,&pnew);CHKERRQ(ierr);
      ierr = PetscArraycpy(pnew,pA,PetscMin(m,bv->m)*bv->n);CHKERRQ(ierr);
      ierr = MatDenseRestoreArray(ctx->A,&pA);CHKERRQ(ierr);
      ierr = MatDenseRestoreArray(A,&pnew);CHKERRQ(ierr);
    }
    ierr = MatDestroy(&ctx->A);CHKERRQ(ierr);
    ctx->A = A;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetColumn_Mat(BV bv,PetscInt j,Vec *v)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *pA;
  PetscInt       l;

  PetscFunctionBegin;
  l = BVAvailableVec;
  ierr = MatDenseGetArray(ctx->A,&pA);CHKERRQ(ierr);
  ierr = VecPlaceArray(bv->cv[l],pA+(bv->nc+j)*bv->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreColumn_Mat(BV bv,PetscInt j,Vec *v)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;
  PetscScalar    *pA;
  PetscInt       l;

  PetscFunctionBegin;
  l = (j==bv->ci[0])? 0: 1;
  ierr = VecResetArray(bv->cv[l]);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(ctx->A,&pA);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetArray_Mat(BV bv,PetscScalar **a)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(ctx->A,a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreArray_Mat(BV bv,PetscScalar **a)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  if (a) { ierr = MatDenseRestoreArray(ctx->A,a);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetArrayRead_Mat(BV bv,const PetscScalar **a)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(ctx->A,(PetscScalar**)a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreArrayRead_Mat(BV bv,const PetscScalar **a)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  if (a) { ierr = MatDenseRestoreArray(ctx->A,(PetscScalar**)a);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode BVView_Mat(BV bv,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  BV_MAT            *ctx = (BV_MAT*)bv->data;
  PetscViewerFormat format;
  PetscBool         isascii;
  const char        *bvname,*name;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
    ierr = MatView(ctx->A,viewer);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_MATLAB) {
      ierr = PetscObjectGetName((PetscObject)bv,&bvname);CHKERRQ(ierr);
      ierr = PetscObjectGetName((PetscObject)ctx->A,&name);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%s=%s;clear %s\n",bvname,name,name);CHKERRQ(ierr);
      if (bv->nc) {
        ierr = PetscViewerASCIIPrintf(viewer,"%s=%s(:,%D:end);\n",bvname,bvname,bv->nc+1);CHKERRQ(ierr);
      }
    }
  } else {
    ierr = MatView(ctx->A,viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVDestroy_Mat(BV bv)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx = (BV_MAT*)bv->data;

  PetscFunctionBegin;
  ierr = MatDestroy(&ctx->A);CHKERRQ(ierr);
  ierr = VecDestroy(&bv->cv[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&bv->cv[1]);CHKERRQ(ierr);
  ierr = PetscFree(bv->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode BVCreate_Mat(BV bv)
{
  PetscErrorCode ierr;
  BV_MAT         *ctx;
  PetscInt       nloc,bs,lsplit;
  PetscBool      seq;
  char           str[50];
  PetscScalar    *array,*ptr;
  BV             parent;

  PetscFunctionBegin;
  ierr = PetscNewLog(bv,&ctx);CHKERRQ(ierr);
  bv->data = (void*)ctx;

  ierr = PetscObjectTypeCompare((PetscObject)bv->t,VECMPI,&ctx->mpi);CHKERRQ(ierr);
  if (!ctx->mpi) {
    ierr = PetscObjectTypeCompare((PetscObject)bv->t,VECSEQ,&seq);CHKERRQ(ierr);
    if (!seq) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Cannot create a BVMAT from a non-standard template vector");
  }

  ierr = VecGetLocalSize(bv->t,&nloc);CHKERRQ(ierr);
  ierr = VecGetBlockSize(bv->t,&bs);CHKERRQ(ierr);

  if (bv->issplit) {
    /* split BV: share the memory of the parent BV */
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    ierr = MatDenseGetArray(((BV_MAT*)parent->data)->A,&array);CHKERRQ(ierr);
    ptr = (bv->issplit==1)? array: array+lsplit*nloc;
    ierr = MatDenseRestoreArray(((BV_MAT*)parent->data)->A,&array);CHKERRQ(ierr);
  } else {
    /* regular BV: allocate memory for the BV entries */
    ptr = NULL;
  }
  ierr = MatCreateDense(PetscObjectComm((PetscObject)bv->t),nloc,PETSC_DECIDE,PETSC_DECIDE,bv->m,ptr,&ctx->A);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(ctx->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(ctx->A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)ctx->A);CHKERRQ(ierr);
  if (((PetscObject)bv)->name) {
    ierr = PetscSNPrintf(str,50,"%s_0",((PetscObject)bv)->name);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)ctx->A,str);CHKERRQ(ierr);
  }

  if (bv->Acreate) {
    ierr = MatCopy(bv->Acreate,ctx->A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = MatDestroy(&bv->Acreate);CHKERRQ(ierr);
  }

  ierr = VecDuplicateEmpty(bv->t,&bv->cv[0]);CHKERRQ(ierr);
  ierr = VecDuplicateEmpty(bv->t,&bv->cv[1]);CHKERRQ(ierr);

  bv->ops->mult             = BVMult_Mat;
  bv->ops->multvec          = BVMultVec_Mat;
  bv->ops->multinplace      = BVMultInPlace_Mat;
  bv->ops->multinplacetrans = BVMultInPlaceTranspose_Mat;
  bv->ops->dot              = BVDot_Mat;
  bv->ops->dotvec           = BVDotVec_Mat;
  bv->ops->dotvec_local     = BVDotVec_Local_Mat;
  bv->ops->scale            = BVScale_Mat;
  bv->ops->norm             = BVNorm_Mat;
  bv->ops->norm_local       = BVNorm_Local_Mat;
  bv->ops->matmult          = BVMatMult_Mat;
  bv->ops->copy             = BVCopy_Mat;
  bv->ops->copycolumn       = BVCopyColumn_Mat;
  bv->ops->resize           = BVResize_Mat;
  bv->ops->getcolumn        = BVGetColumn_Mat;
  bv->ops->restorecolumn    = BVRestoreColumn_Mat;
  bv->ops->getarray         = BVGetArray_Mat;
  bv->ops->restorearray     = BVRestoreArray_Mat;
  bv->ops->getarrayread     = BVGetArrayRead_Mat;
  bv->ops->restorearrayread = BVRestoreArrayRead_Mat;
  bv->ops->destroy          = BVDestroy_Mat;
  if (!ctx->mpi) bv->ops->view = BVView_Mat;
  PetscFunctionReturn(0);
}

