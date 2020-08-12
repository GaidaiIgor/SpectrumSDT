/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV implemented as an array of independent Vecs
*/

#include <slepc/private/bvimpl.h>

typedef struct {
  Vec      *V;
  PetscInt vmip;   /* Version of BVMultInPlace:
       0: memory-efficient version, uses VecGetArray (default in CPU)
       1: version that allocates (e-s) work vectors in every call (default in GPU) */
} BV_VECS;

PetscErrorCode BVMult_Vecs(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  PetscErrorCode ierr;
  BV_VECS        *y = (BV_VECS*)Y->data,*x = (BV_VECS*)X->data;
  PetscScalar    *q,*s=NULL;
  PetscInt       i,j,ldq;

  PetscFunctionBegin;
  if (Q) {
    ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
    if (alpha!=1.0) {
      ierr = BVAllocateWork_Private(Y,X->k-X->l);CHKERRQ(ierr);
      s = Y->work;
    }
    ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
    for (j=Y->l;j<Y->k;j++) {
      ierr = VecScale(y->V[Y->nc+j],beta);CHKERRQ(ierr);
      if (alpha!=1.0) {
        for (i=X->l;i<X->k;i++) s[i-X->l] = alpha*q[i+j*ldq];
      } else s = q+j*ldq+X->l;
      ierr = VecMAXPY(y->V[Y->nc+j],X->k-X->l,s,x->V+X->nc+X->l);CHKERRQ(ierr);
    }
    ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  } else {
    for (j=0;j<Y->k-Y->l;j++) {
      ierr = VecScale(y->V[Y->nc+Y->l+j],beta);CHKERRQ(ierr);
      ierr = VecAXPY(y->V[Y->nc+Y->l+j],alpha,x->V[X->nc+X->l+j]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVMultVec_Vecs(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar *q)
{
  PetscErrorCode ierr;
  BV_VECS        *x = (BV_VECS*)X->data;
  PetscScalar    *s=NULL,*qq=q;
  PetscInt       i;

  PetscFunctionBegin;
  if (alpha!=1.0) {
    ierr = BVAllocateWork_Private(X,X->k-X->l);CHKERRQ(ierr);
    s = X->work;
  }
  if (!q) { ierr = VecGetArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = VecScale(y,beta);CHKERRQ(ierr);
  if (alpha!=1.0) {
    for (i=0;i<X->k-X->l;i++) s[i] = alpha*qq[i];
  } else s = qq;
  ierr = VecMAXPY(y,X->k-X->l,s,x->V+X->nc+X->l);CHKERRQ(ierr);
  if (!q) { ierr = VecRestoreArray(X->buffer,&qq);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   BVMultInPlace_Vecs_ME - V(:,s:e-1) = V*Q(:,s:e-1) for regular vectors.

   Memory-efficient version, uses VecGetArray (default in CPU)

   Writing V = [ V1 V2 V3 ] and Q(:,s:e-1) = [ Q1 Q2 Q3 ]', where V2
   corresponds to the columns s:e-1, the computation is done as
                  V2 := V2*Q2 + V1*Q1 + V3*Q3
*/
PetscErrorCode BVMultInPlace_Vecs_ME(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)V->data;
  PetscScalar    *q;
  PetscInt       i,ldq;

  PetscFunctionBegin;
  ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  /* V2 := V2*Q2 */
  ierr = BVMultInPlace_Vecs_Private(V,V->n,e-s,ldq,ctx->V+V->nc+s,q+s*ldq+s,PETSC_FALSE);CHKERRQ(ierr);
  /* V2 += V1*Q1 + V3*Q3 */
  for (i=s;i<e;i++) {
    if (s>V->l) {
      ierr = VecMAXPY(ctx->V[V->nc+i],s-V->l,q+i*ldq+V->l,ctx->V+V->nc+V->l);CHKERRQ(ierr);
    }
    if (V->k>e) {
      ierr = VecMAXPY(ctx->V[V->nc+i],V->k-e,q+i*ldq+e,ctx->V+V->nc+e);CHKERRQ(ierr);
    }
  }
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   BVMultInPlace_Vecs_Alloc - V(:,s:e-1) = V*Q(:,s:e-1) for regular vectors.

   Version that allocates (e-s) work vectors in every call (default in GPU)
*/
PetscErrorCode BVMultInPlace_Vecs_Alloc(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)V->data;
  PetscScalar    *q;
  PetscInt       i,ldq;
  Vec            *W;

  PetscFunctionBegin;
  ierr = MatGetSize(Q,&ldq,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(V->t,e-s,&W);CHKERRQ(ierr);
  for (i=s;i<e;i++) {
    ierr = VecMAXPY(W[i-s],V->k-V->l,q+i*ldq+V->l,ctx->V+V->nc+V->l);CHKERRQ(ierr);
  }
  for (i=s;i<e;i++) {
    ierr = VecCopy(W[i-s],ctx->V[V->nc+i]);CHKERRQ(ierr);
  }
  ierr = VecDestroyVecs(e-s,&W);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   BVMultInPlaceTranspose_Vecs - V(:,s:e-1) = V*Q'(:,s:e-1) for regular vectors.
*/
PetscErrorCode BVMultInPlaceTranspose_Vecs(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)V->data;
  PetscScalar    *q;
  PetscInt       i,j,ldq,n;

  PetscFunctionBegin;
  ierr = MatGetSize(Q,&ldq,&n);CHKERRQ(ierr);
  ierr = MatDenseGetArray(Q,&q);CHKERRQ(ierr);
  /* V2 := V2*Q2' */
  ierr = BVMultInPlace_Vecs_Private(V,V->n,e-s,ldq,ctx->V+V->nc+s,q+s*ldq+s,PETSC_TRUE);CHKERRQ(ierr);
  /* V2 += V1*Q1' + V3*Q3' */
  for (i=s;i<e;i++) {
    for (j=V->l;j<s;j++) {
      ierr = VecAXPY(ctx->V[V->nc+i],q[i+j*ldq],ctx->V[V->nc+j]);CHKERRQ(ierr);
    }
    for (j=e;j<n;j++) {
      ierr = VecAXPY(ctx->V[V->nc+i],q[i+j*ldq],ctx->V[V->nc+j]);CHKERRQ(ierr);
    }
  }
  ierr = MatDenseRestoreArray(Q,&q);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDot_Vecs(BV X,BV Y,Mat M)
{
  PetscErrorCode ierr;
  BV_VECS        *x = (BV_VECS*)X->data,*y = (BV_VECS*)Y->data;
  PetscScalar    *m;
  PetscInt       j,ldm;

  PetscFunctionBegin;
  ierr = MatGetSize(M,&ldm,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(M,&m);CHKERRQ(ierr);
  for (j=X->l;j<X->k;j++) {
    ierr = VecMDot(x->V[X->nc+j],Y->k-Y->l,y->V+Y->nc+Y->l,m+j*ldm+Y->l);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(M,&m);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_Vecs(BV X,Vec y,PetscScalar *q)
{
  PetscErrorCode ierr;
  BV_VECS        *x = (BV_VECS*)X->data;
  Vec            z = y;
  PetscScalar    *qq=q;

  PetscFunctionBegin;
  if (X->matrix) {
    ierr = BV_IPMatMult(X,y);CHKERRQ(ierr);
    z = X->Bx;
  }
  if (!q) { ierr = VecGetArray(X->buffer,&qq);CHKERRQ(ierr); }
  ierr = VecMDot(z,X->k-X->l,x->V+X->nc+X->l,qq);CHKERRQ(ierr);
  if (!q) { ierr = VecRestoreArray(X->buffer,&qq);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_Begin_Vecs(BV X,Vec y,PetscScalar *m)
{
  PetscErrorCode ierr;
  BV_VECS        *x = (BV_VECS*)X->data;
  Vec            z = y;

  PetscFunctionBegin;
  if (X->matrix) {
    ierr = BV_IPMatMult(X,y);CHKERRQ(ierr);
    z = X->Bx;
  }
  ierr = VecMDotBegin(z,X->k-X->l,x->V+X->nc+X->l,m);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDotVec_End_Vecs(BV X,Vec y,PetscScalar *m)
{
  PetscErrorCode ierr;
  BV_VECS        *x = (BV_VECS*)X->data;

  PetscFunctionBegin;
  ierr = VecMDotEnd(y,X->k-X->l,x->V+X->nc+X->l,m);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVScale_Vecs(BV bv,PetscInt j,PetscScalar alpha)
{
  PetscErrorCode ierr;
  PetscInt       i;
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  if (j<0) {
    for (i=bv->l;i<bv->k;i++) {
      ierr = VecScale(ctx->V[bv->nc+i],alpha);CHKERRQ(ierr);
    }
  } else {
    ierr = VecScale(ctx->V[bv->nc+j],alpha);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_Vecs(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      nrm;
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  if (j<0) {
    if (type==NORM_FROBENIUS) {
      *val = 0.0;
      for (i=bv->l;i<bv->k;i++) {
        ierr = VecNorm(ctx->V[bv->nc+i],NORM_2,&nrm);CHKERRQ(ierr);
        *val += nrm*nrm;
      }
      *val = PetscSqrtReal(*val);
    } else SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not implemented in BVVECS");
  } else {
    ierr = VecNorm(ctx->V[bv->nc+j],type,val);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_Begin_Vecs(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  if (j<0) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not implemented in BVVECS");
  else {
    ierr = VecNormBegin(ctx->V[bv->nc+j],type,val);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVNorm_End_Vecs(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  if (j<0) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not implemented in BVVECS");
  else {
    ierr = VecNormEnd(ctx->V[bv->nc+j],type,val);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVMatMult_Vecs(BV V,Mat A,BV W)
{
  PetscErrorCode ierr;
  BV_VECS        *v = (BV_VECS*)V->data,*w = (BV_VECS*)W->data;
  PetscInt       j;
  PetscBool      flg;
  Mat            Vmat,Wmat;

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
      ierr = MatMult(A,v->V[V->nc+V->l+j],w->V[W->nc+W->l+j]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVCopy_Vecs(BV V,BV W)
{
  PetscErrorCode ierr;
  BV_VECS        *v = (BV_VECS*)V->data,*w = (BV_VECS*)W->data;
  PetscInt       j;

  PetscFunctionBegin;
  for (j=0;j<V->k-V->l;j++) {
    ierr = VecCopy(v->V[V->nc+V->l+j],w->V[W->nc+W->l+j]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVCopyColumn_Vecs(BV V,PetscInt j,PetscInt i)
{
  PetscErrorCode ierr;
  BV_VECS        *v = (BV_VECS*)V->data;

  PetscFunctionBegin;
  ierr = VecCopy(v->V[V->nc+j],v->V[V->nc+i]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVResize_Vecs(BV bv,PetscInt m,PetscBool copy)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)bv->data;
  Vec            *newV;
  PetscInt       j,lsplit;
  char           str[50];
  BV             parent;

  PetscFunctionBegin;
  if (bv->issplit==2) {
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    ctx->V = ((BV_VECS*)parent->data)->V+lsplit;
  } else if (!bv->issplit) {
    ierr = VecDuplicateVecs(bv->t,m,&newV);CHKERRQ(ierr);
    ierr = PetscLogObjectParents(bv,m,newV);CHKERRQ(ierr);
    if (((PetscObject)bv)->name) {
      for (j=0;j<m;j++) {
        ierr = PetscSNPrintf(str,50,"%s_%D",((PetscObject)bv)->name,j);CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)newV[j],str);CHKERRQ(ierr);
      }
    }
    if (copy) {
      for (j=0;j<PetscMin(m,bv->m);j++) {
        ierr = VecCopy(ctx->V[j],newV[j]);CHKERRQ(ierr);
      }
    }
    ierr = VecDestroyVecs(bv->m,&ctx->V);CHKERRQ(ierr);
    ctx->V = newV;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetColumn_Vecs(BV bv,PetscInt j,Vec *v)
{
  BV_VECS  *ctx = (BV_VECS*)bv->data;
  PetscInt l;

  PetscFunctionBegin;
  l = BVAvailableVec;
  bv->cv[l] = ctx->V[bv->nc+j];
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetArray_Vecs(BV bv,PetscScalar **a)
{
  PetscErrorCode    ierr;
  BV_VECS           *ctx = (BV_VECS*)bv->data;
  PetscInt          j;
  const PetscScalar *p;

  PetscFunctionBegin;
  ierr = PetscMalloc1((bv->nc+bv->m)*bv->n,a);CHKERRQ(ierr);
  for (j=0;j<bv->nc+bv->m;j++) {
    ierr = VecGetArrayRead(ctx->V[j],&p);CHKERRQ(ierr);
    ierr = PetscArraycpy(*a+j*bv->n,p,bv->n);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(ctx->V[j],&p);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreArray_Vecs(BV bv,PetscScalar **a)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)bv->data;
  PetscInt       j;
  PetscScalar    *p;

  PetscFunctionBegin;
  for (j=0;j<bv->nc+bv->m;j++) {
    ierr = VecGetArray(ctx->V[j],&p);CHKERRQ(ierr);
    ierr = PetscArraycpy(p,*a+j*bv->n,bv->n);CHKERRQ(ierr);
    ierr = VecRestoreArray(ctx->V[j],&p);CHKERRQ(ierr);
  }
  ierr = PetscFree(*a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVGetArrayRead_Vecs(BV bv,const PetscScalar **a)
{
  PetscErrorCode    ierr;
  BV_VECS           *ctx = (BV_VECS*)bv->data;
  PetscInt          j;
  const PetscScalar *p;

  PetscFunctionBegin;
  ierr = PetscMalloc1((bv->nc+bv->m)*bv->n,(PetscScalar**)a);CHKERRQ(ierr);
  for (j=0;j<bv->nc+bv->m;j++) {
    ierr = VecGetArrayRead(ctx->V[j],&p);CHKERRQ(ierr);
    ierr = PetscArraycpy((PetscScalar*)*a+j*bv->n,p,bv->n);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(ctx->V[j],&p);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVRestoreArrayRead_Vecs(BV bv,const PetscScalar **a)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(*a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Sets the value of vmip flag and resets ops->multinplace accordingly
 */
PETSC_STATIC_INLINE PetscErrorCode BVVecsSetVmip(BV bv,PetscInt vmip)
{
  typedef PetscErrorCode (*fmultinplace)(BV,Mat,PetscInt,PetscInt);
  fmultinplace multinplace[2] = {BVMultInPlace_Vecs_ME, BVMultInPlace_Vecs_Alloc};
  BV_VECS      *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  ctx->vmip            = vmip;
  bv->ops->multinplace = multinplace[vmip];
  PetscFunctionReturn(0);
}

PetscErrorCode BVSetFromOptions_Vecs(PetscOptionItems *PetscOptionsObject,BV bv)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"BV Vecs Options");CHKERRQ(ierr);

    ierr = PetscOptionsRangeInt("-bv_vecs_vmip","Version of BVMultInPlace operation","",ctx->vmip,&ctx->vmip,NULL,0,1);CHKERRQ(ierr);
    ierr = BVVecsSetVmip(bv,ctx->vmip);CHKERRQ(ierr);

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVView_Vecs(BV bv,PetscViewer viewer)
{
  PetscErrorCode    ierr;
  BV_VECS           *ctx = (BV_VECS*)bv->data;
  PetscInt          j;
  PetscViewerFormat format;
  PetscBool         isascii,ismatlab=PETSC_FALSE;
  const char        *bvname,*name;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    if (format == PETSC_VIEWER_ASCII_INFO || format == PETSC_VIEWER_ASCII_INFO_DETAIL) PetscFunctionReturn(0);
    if (format == PETSC_VIEWER_ASCII_MATLAB) ismatlab = PETSC_TRUE;
  }
  if (ismatlab) {
    ierr = PetscObjectGetName((PetscObject)bv,&bvname);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%s=[];\n",bvname);CHKERRQ(ierr);
  }
  for (j=bv->nc;j<bv->nc+bv->m;j++) {
    ierr = VecView(ctx->V[j],viewer);CHKERRQ(ierr);
    if (ismatlab) {
      ierr = PetscObjectGetName((PetscObject)ctx->V[j],&name);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%s=[%s,%s];clear %s\n",bvname,bvname,name,name);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BVDestroy_Vecs(BV bv)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)bv->data;

  PetscFunctionBegin;
  if (!bv->issplit) { ierr = VecDestroyVecs(bv->nc+bv->m,&ctx->V);CHKERRQ(ierr); }
  ierr = PetscFree(bv->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode BVDuplicate_Vecs(BV V,BV W)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx = (BV_VECS*)V->data;

  PetscFunctionBegin;
  ierr = BVVecsSetVmip(W,ctx->vmip);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode BVCreate_Vecs(BV bv)
{
  PetscErrorCode ierr;
  BV_VECS        *ctx;
  PetscInt       j,lsplit;
  PetscBool      isgpu;
  char           str[50];
  BV             parent;
  Vec            *Vpar;

  PetscFunctionBegin;
  ierr = PetscNewLog(bv,&ctx);CHKERRQ(ierr);
  bv->data = (void*)ctx;

  if (bv->issplit) {
    /* split BV: share the Vecs of the parent BV */
    parent = bv->splitparent;
    lsplit = parent->lsplit;
    Vpar   = ((BV_VECS*)parent->data)->V;
    ctx->V = (bv->issplit==1)? Vpar: Vpar+lsplit;
  } else {
    /* regular BV: create array of Vecs to store the BV columns */
    ierr = VecDuplicateVecs(bv->t,bv->m,&ctx->V);CHKERRQ(ierr);
    ierr = PetscLogObjectParents(bv,bv->m,ctx->V);CHKERRQ(ierr);
    if (((PetscObject)bv)->name) {
      for (j=0;j<bv->m;j++) {
        ierr = PetscSNPrintf(str,50,"%s_%D",((PetscObject)bv)->name,j);CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)ctx->V[j],str);CHKERRQ(ierr);
      }
    }
  }

  if (bv->Acreate) {
    for (j=0;j<bv->m;j++) {
      ierr = MatGetColumnVector(bv->Acreate,ctx->V[j],j);CHKERRQ(ierr);
    }
    ierr = MatDestroy(&bv->Acreate);CHKERRQ(ierr);
  }

  /* Default version of BVMultInPlace */
  ierr = PetscObjectTypeCompareAny((PetscObject)bv->t,&isgpu,VECSEQCUDA,VECMPICUDA,"");CHKERRQ(ierr);
  ctx->vmip = isgpu? 1: 0;

  /* Default BVMatMult method */
  bv->vmm = BV_MATMULT_VECS;

  /* Deferred call to setfromoptions */
  if (bv->defersfo) {
    ierr = PetscObjectOptionsBegin((PetscObject)bv);CHKERRQ(ierr);
    ierr = BVSetFromOptions_Vecs(PetscOptionsObject,bv);CHKERRQ(ierr);
    ierr = PetscOptionsEnd();CHKERRQ(ierr);
  }
  ierr = BVVecsSetVmip(bv,ctx->vmip);CHKERRQ(ierr);

  bv->ops->mult             = BVMult_Vecs;
  bv->ops->multvec          = BVMultVec_Vecs;
  bv->ops->multinplacetrans = BVMultInPlaceTranspose_Vecs;
  bv->ops->dot              = BVDot_Vecs;
  bv->ops->dotvec           = BVDotVec_Vecs;
  bv->ops->dotvec_begin     = BVDotVec_Begin_Vecs;
  bv->ops->dotvec_end       = BVDotVec_End_Vecs;
  bv->ops->scale            = BVScale_Vecs;
  bv->ops->norm             = BVNorm_Vecs;
  bv->ops->norm_begin       = BVNorm_Begin_Vecs;
  bv->ops->norm_end         = BVNorm_End_Vecs;
  bv->ops->matmult          = BVMatMult_Vecs;
  bv->ops->copy             = BVCopy_Vecs;
  bv->ops->copycolumn       = BVCopyColumn_Vecs;
  bv->ops->resize           = BVResize_Vecs;
  bv->ops->getcolumn        = BVGetColumn_Vecs;
  bv->ops->getarray         = BVGetArray_Vecs;
  bv->ops->restorearray     = BVRestoreArray_Vecs;
  bv->ops->getarrayread     = BVGetArrayRead_Vecs;
  bv->ops->restorearrayread = BVRestoreArrayRead_Vecs;
  bv->ops->destroy          = BVDestroy_Vecs;
  bv->ops->duplicate        = BVDuplicate_Vecs;
  bv->ops->setfromoptions   = BVSetFromOptions_Vecs;
  bv->ops->view             = BVView_Vecs;
  PetscFunctionReturn(0);
}

