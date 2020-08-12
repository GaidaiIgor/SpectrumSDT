/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <slepcblaslapack.h>
#include "nepdefl.h"

PetscErrorCode NEPDeflationGetInvariantPair(NEP_EXT_OP extop,BV *X,Mat *H)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (X) *X = extop->X;
  if (H) {
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,extop->szd+1,extop->szd+1,extop->H,H);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationExtendInvariantPair(NEP_EXT_OP extop,Vec u,PetscScalar lambda,PetscInt k)
{
  PetscErrorCode ierr;
  Vec            uu;
  PetscInt       ld,i;
  PetscMPIInt    np;
  PetscReal      norm;

  PetscFunctionBegin;
  ierr = BVGetColumn(extop->X,k,&uu);CHKERRQ(ierr);
  ld = extop->szd+1;
  ierr = NEPDeflationCopyToExtendedVec(extop,uu,extop->H+k*ld,u,PETSC_TRUE);CHKERRQ(ierr);
  ierr = BVRestoreColumn(extop->X,k,&uu);CHKERRQ(ierr);
  ierr = BVNormColumn(extop->X,k,NORM_2,&norm);CHKERRQ(ierr);
  ierr = BVScaleColumn(extop->X,k,1.0/norm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)u),&np);CHKERRQ(ierr);
  for (i=0;i<k;i++) extop->H[k*ld+i] *= PetscSqrtReal(np)/norm;
  extop->H[k*(ld+1)] = lambda;
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationExtractEigenpair(NEP_EXT_OP extop,PetscInt k,Vec u,PetscScalar lambda,DS ds)
{
  PetscErrorCode ierr;
  PetscScalar    *Ap;
  PetscInt       ldh=extop->szd+1,ldds,i,j,k1=k+1;
  PetscScalar    *eigr,*eigi,*t,*Z;
  Vec            x;

  PetscFunctionBegin;
  ierr = NEPDeflationExtendInvariantPair(extop,u,lambda,k);CHKERRQ(ierr);
  ierr = PetscCalloc3(k1,&eigr,k1,&eigi,extop->szd,&t);CHKERRQ(ierr);
  ierr = DSReset(ds);CHKERRQ(ierr);
  ierr = DSSetType(ds,DSNHEP);CHKERRQ(ierr);
  ierr = DSAllocate(ds,ldh);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(ds,&ldds);CHKERRQ(ierr);
  ierr = DSGetArray(ds,DS_MAT_A,&Ap);CHKERRQ(ierr);
  for (j=0;j<k1;j++)
    for (i=0;i<k1;i++) Ap[j*ldds+i] = extop->H[j*ldh+i];
  ierr = DSRestoreArray(ds,DS_MAT_A,&Ap);CHKERRQ(ierr);
  ierr = DSSetDimensions(ds,k1,0,0,k1);CHKERRQ(ierr);
  ierr = DSSolve(ds,eigr,eigi);CHKERRQ(ierr);
  ierr = DSVectors(ds,DS_MAT_X,&k,NULL);CHKERRQ(ierr);
  ierr = DSGetArray(ds,DS_MAT_X,&Z);CHKERRQ(ierr);
  ierr = BVMultColumn(extop->X,1.0,Z[k*ldds+k],k,Z+k*ldds);CHKERRQ(ierr);
  ierr = DSRestoreArray(ds,DS_MAT_X,&Z);CHKERRQ(ierr);
  ierr = BVGetColumn(extop->X,k,&x);CHKERRQ(ierr);
  ierr = NEPDeflationCopyToExtendedVec(extop,x,t,u,PETSC_FALSE);CHKERRQ(ierr);
  ierr = BVRestoreColumn(extop->X,k,&x);CHKERRQ(ierr);
  ierr = PetscFree3(eigr,eigi,t);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationCopyToExtendedVec(NEP_EXT_OP extop,Vec v,PetscScalar *a,Vec vex,PetscBool back)
{
  PetscErrorCode ierr;
  PetscMPIInt    np,rk,count;
  PetscScalar    *array1,*array2;
  PetscInt       nloc;

  PetscFunctionBegin;
  if (extop->szd) {
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)vex),&rk);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)vex),&np);CHKERRQ(ierr);
    ierr = BVGetSizes(extop->nep->V,&nloc,NULL,NULL);CHKERRQ(ierr);
    if (v) {
      ierr = VecGetArray(v,&array1);CHKERRQ(ierr);
      ierr = VecGetArray(vex,&array2);CHKERRQ(ierr);
      if (back) {
        ierr = PetscArraycpy(array1,array2,nloc);CHKERRQ(ierr);
      } else {
        ierr = PetscArraycpy(array2,array1,nloc);CHKERRQ(ierr);
      }
      ierr = VecRestoreArray(v,&array1);CHKERRQ(ierr);
      ierr = VecRestoreArray(vex,&array2);CHKERRQ(ierr);
    }
    if (a) {
      ierr = VecGetArray(vex,&array2);CHKERRQ(ierr);
      if (back) {
        ierr = PetscArraycpy(a,array2+nloc,extop->szd);CHKERRQ(ierr);
        ierr = PetscMPIIntCast(extop->szd,&count);CHKERRQ(ierr);
        ierr = MPI_Bcast(a,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)vex));CHKERRQ(ierr);
      } else {
        ierr = PetscArraycpy(array2+nloc,a,extop->szd);CHKERRQ(ierr);
        ierr = PetscMPIIntCast(extop->szd,&count);CHKERRQ(ierr);
        ierr = MPI_Bcast(array2+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)vex));CHKERRQ(ierr);
      }
      ierr = VecRestoreArray(vex,&array2);CHKERRQ(ierr);
    }
  } else {
    if (back) {ierr = VecCopy(vex,v);CHKERRQ(ierr);}
    else {ierr = VecCopy(v,vex);CHKERRQ(ierr);}
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationCreateVec(NEP_EXT_OP extop,Vec *v)
{
  PetscErrorCode ierr;
  PetscInt       nloc;
  Vec            u;
  VecType        type;

  PetscFunctionBegin;
  if (extop->szd) {
    ierr = BVGetColumn(extop->nep->V,0,&u);CHKERRQ(ierr);
    ierr = VecGetType(u,&type);CHKERRQ(ierr);
    ierr = BVRestoreColumn(extop->nep->V,0,&u);CHKERRQ(ierr);
    ierr = VecCreate(PetscObjectComm((PetscObject)extop->nep),v);CHKERRQ(ierr);
    ierr = VecSetType(*v,type);CHKERRQ(ierr);
    ierr = BVGetSizes(extop->nep->V,&nloc,NULL,NULL);CHKERRQ(ierr);
    nloc += extop->szd;
    ierr = VecSetSizes(*v,nloc,PETSC_DECIDE);CHKERRQ(ierr);
  } else {
    ierr = BVCreateVec(extop->nep->V,v);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationCreateBV(NEP_EXT_OP extop,PetscInt sz,BV *V)
{
  PetscErrorCode     ierr;
  PetscInt           nloc;
  BVType             type;
  BVOrthogType       otype;
  BVOrthogRefineType oref;
  PetscReal          oeta;
  BVOrthogBlockType  oblock;
  NEP                nep=extop->nep;

  PetscFunctionBegin;
  if (extop->szd) {
    ierr = BVGetSizes(nep->V,&nloc,NULL,NULL);CHKERRQ(ierr);
    ierr = BVCreate(PetscObjectComm((PetscObject)nep),V);CHKERRQ(ierr);
    ierr = BVSetSizes(*V,nloc+extop->szd,PETSC_DECIDE,sz);CHKERRQ(ierr);
    ierr = BVGetType(nep->V,&type);CHKERRQ(ierr);
    ierr = BVSetType(*V,type);CHKERRQ(ierr);
    ierr = BVGetOrthogonalization(nep->V,&otype,&oref,&oeta,&oblock);CHKERRQ(ierr);
    ierr = BVSetOrthogonalization(*V,otype,oref,oeta,oblock);CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease((PetscObject)*V);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)*V);CHKERRQ(ierr);
  } else {
    ierr = BVDuplicateResize(nep->V,sz,V);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationSetRandomVec(NEP_EXT_OP extop,Vec v)
{
  PetscErrorCode ierr;
  PetscInt       n,next,i;
  PetscRandom    rand;
  PetscScalar    *array;
  PetscMPIInt    nn,np;

  PetscFunctionBegin;
  ierr = BVGetRandomContext(extop->nep->V,&rand);CHKERRQ(ierr);
  ierr = VecSetRandom(v,rand);CHKERRQ(ierr);
  if (extop->szd) {
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)v),&np);CHKERRQ(ierr);
    ierr = BVGetSizes(extop->nep->V,&n,NULL,NULL);CHKERRQ(ierr);
    ierr = VecGetLocalSize(v,&next);CHKERRQ(ierr);
    ierr = VecGetArray(v,&array);CHKERRQ(ierr);
    for (i=n+extop->n;i<next;i++) array[i] = 0.0;
    for (i=n;i<n+extop->n;i++) array[i] /= PetscSqrtReal(np);
    ierr = PetscMPIIntCast(extop->n,&nn);CHKERRQ(ierr);
    ierr = MPI_Bcast(array+n,nn,MPIU_SCALAR,0,PetscObjectComm((PetscObject)v));CHKERRQ(ierr);
    ierr = VecRestoreArray(v,&array);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationEvaluateBasisMat(NEP_EXT_OP extop,PetscInt idx,PetscBool hat,PetscScalar *bval,PetscScalar *Hj,PetscScalar *Hjprev)
{
  PetscErrorCode ierr;
  PetscInt       i,k,n=extop->n,ldhj=extop->szd,ldh=extop->szd+1;
  PetscScalar    sone=1.0,zero=0.0;
  PetscBLASInt   ldh_,ldhj_,n_;

  PetscFunctionBegin;
  i = (idx<0)?extop->szd*extop->szd*(-idx):extop->szd*extop->szd;
  ierr = PetscArrayzero(Hj,i);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldhj+1,&ldh_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldhj,&ldhj_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  if (idx<1) {
    if (!hat) for (i=0;i<extop->n;i++) Hj[i+i*ldhj] = 1.0;
    else for (i=0;i<extop->n;i++) Hj[i+i*ldhj] = 0.0;
  } else {
      for (i=0;i<n;i++) extop->H[i*ldh+i] -= extop->bc[idx-1];
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,extop->H,&ldh_,Hjprev,&ldhj_,&zero,Hj,&ldhj_));
      for (i=0;i<n;i++) extop->H[i*ldh+i] += extop->bc[idx-1];
      if (hat) for (i=0;i<n;i++) Hj[i*(ldhj+1)] += bval[idx-1];
  }
  if (idx<0) {
    idx = -idx;
    for (k=1;k<idx;k++) {
      for (i=0;i<n;i++) extop->H[i*ldh+i] -= extop->bc[k-1];
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,extop->H,&ldh_,Hj+(k-1)*ldhj*ldhj,&ldhj_,&zero,Hj+k*ldhj*ldhj,&ldhj_));
      for (i=0;i<n;i++) extop->H[i*ldh+i] += extop->bc[k-1];
      if (hat) for (i=0;i<n;i++) Hj[i*(ldhj+1)] += bval[k-1];
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationLocking(NEP_EXT_OP extop,Vec u,PetscScalar lambda)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = NEPDeflationExtendInvariantPair(extop,u,lambda,extop->n);CHKERRQ(ierr);
  extop->n++;
  ierr = BVSetActiveColumns(extop->X,0,extop->n);CHKERRQ(ierr);
  if (extop->n <= extop->szd) {
    /* update XpX */
    ierr = BVDotColumn(extop->X,extop->n-1,extop->XpX+(extop->n-1)*extop->szd);CHKERRQ(ierr);
    extop->XpX[(extop->n-1)*(1+extop->szd)] = 1.0;
    for (i=0;i<extop->n-1;i++) extop->XpX[i*extop->szd+extop->n-1] = PetscConj(extop->XpX[(extop->n-1)*extop->szd+i]);
    /* determine minimality index */
    extop->midx = PetscMin(extop->max_midx,extop->n);
    /* polynominal basis coeficients */
    for (i=0;i<extop->midx;i++) extop->bc[i] = extop->nep->target;
    /* evaluate the polynomial basis in H */
    ierr = NEPDeflationEvaluateBasisMat(extop,-extop->midx,PETSC_FALSE,NULL,extop->Hj,NULL);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationEvaluateHatFunction(NEP_EXT_OP extop, PetscInt idx,PetscScalar lambda,PetscScalar *y,PetscScalar *hfj,PetscScalar *hfjp,PetscInt ld)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,off,ini,fin,sz,ldh,n=extop->n;
  Mat            A,B;
  PetscScalar    *array;

  PetscFunctionBegin;
  if (idx<0) {ini = 0; fin = extop->nep->nt;}
  else {ini = idx; fin = idx+1;}
  if (y) sz = hfjp?n+2:n+1;
  else sz = hfjp?3*n:2*n;
  ldh = extop->szd+1;
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,sz,sz,NULL,&A);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,sz,sz,NULL,&B);CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&array);CHKERRQ(ierr);
  for (j=0;j<n;j++)
    for (i=0;i<n;i++) array[j*sz+i] = extop->H[j*ldh+i];
  ierr = MatDenseRestoreArray(A,&array);CHKERRQ(ierr);
  if (y) {
    ierr = MatDenseGetArray(A,&array);CHKERRQ(ierr);
    array[extop->n*(sz+1)] = lambda;
    if (hfjp) { array[(n+1)*sz+n] = 1.0; array[(n+1)*sz+n+1] = lambda;}
    for (i=0;i<n;i++) array[n*sz+i] = y[i];
    ierr = MatDenseRestoreArray(A,&array);CHKERRQ(ierr);
    for (j=ini;j<fin;j++) {
      ierr = FNEvaluateFunctionMat(extop->nep->f[j],A,B);CHKERRQ(ierr);
      ierr = MatDenseGetArray(B,&array);CHKERRQ(ierr);
      for (i=0;i<n;i++) hfj[j*ld+i] = array[n*sz+i];
      if (hfjp) for (i=0;i<n;i++) hfjp[j*ld+i] = array[(n+1)*sz+i];
      ierr = MatDenseRestoreArray(B,&array);CHKERRQ(ierr);
    }
  } else {
    off = idx<0?ld*n:0;
    ierr = MatDenseGetArray(A,&array);CHKERRQ(ierr);
    for (k=0;k<n;k++) {
      array[(n+k)*sz+k] = 1.0;
      array[(n+k)*sz+n+k] = lambda;
    }
    if (hfjp) for (k=0;k<n;k++) {
      array[(2*n+k)*sz+n+k] = 1.0;
      array[(2*n+k)*sz+2*n+k] = lambda;
    }
    ierr = MatDenseRestoreArray(A,&array);CHKERRQ(ierr);
    for (j=ini;j<fin;j++) {
      ierr = FNEvaluateFunctionMat(extop->nep->f[j],A,B);CHKERRQ(ierr);
      ierr = MatDenseGetArray(B,&array);CHKERRQ(ierr);
      for (i=0;i<n;i++) for (k=0;k<n;k++) hfj[j*off+i*ld+k] = array[n*sz+i*sz+k];
      if (hfjp) for (k=0;k<n;k++) for (i=0;i<n;i++) hfjp[j*off+i*ld+k] = array[2*n*sz+i*sz+k];
      ierr = MatDenseRestoreArray(B,&array);CHKERRQ(ierr);
    }
  }
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_NEPDeflation(Mat M,Vec x,Vec y)
{
  NEP_DEF_MATSHELL  *matctx;
  PetscErrorCode    ierr;
  NEP_EXT_OP        extop;
  Vec               x1,y1;
  PetscScalar       *yy,sone=1.0,zero=0.0;
  const PetscScalar *xx;
  PetscInt          nloc,i;
  PetscMPIInt       np;
  PetscBLASInt      n_,one=1,szd_;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)M),&np);CHKERRQ(ierr);
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  extop = matctx->extop;
  if (extop->ref) {
    ierr = VecZeroEntries(y);CHKERRQ(ierr);
  }
  if (extop->szd) {
    x1 = matctx->w[0]; y1 = matctx->w[1];
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecPlaceArray(x1,xx);CHKERRQ(ierr);
    ierr = VecGetArray(y,&yy);CHKERRQ(ierr);
    ierr = VecPlaceArray(y1,yy);CHKERRQ(ierr);
    ierr = MatMult(matctx->T,x1,y1);CHKERRQ(ierr);
    if (!extop->ref && extop->n) {
      ierr = VecGetLocalSize(x1,&nloc);CHKERRQ(ierr);
      /* copy for avoiding warning of constant array xx */
      for (i=0;i<extop->n;i++) matctx->work[i] = xx[nloc+i]*PetscSqrtReal(np);
      ierr = BVMultVec(matctx->U,1.0,1.0,y1,matctx->work);CHKERRQ(ierr);
      ierr = BVDotVec(extop->X,x1,matctx->work);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(extop->n,&n_);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(extop->szd,&szd_);CHKERRQ(ierr);
      PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,matctx->A,&szd_,matctx->work,&one,&zero,yy+nloc,&one));
      PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,matctx->B,&szd_,xx+nloc,&one,&sone,yy+nloc,&one));
      for (i=0;i<extop->n;i++) yy[nloc+i] /= PetscSqrtReal(np);
    }
    ierr = VecResetArray(x1);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecResetArray(y1);CHKERRQ(ierr);
    ierr = VecRestoreArray(y,&yy);CHKERRQ(ierr);
  } else {
    ierr = MatMult(matctx->T,x,y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatCreateVecs_NEPDeflation(Mat M,Vec *right,Vec *left)
{
  PetscErrorCode   ierr;
  NEP_DEF_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  if (right) {
    ierr = VecDuplicate(matctx->w[0],right);CHKERRQ(ierr);
  }
  if (left) {
    ierr = VecDuplicate(matctx->w[0],left);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDestroy_NEPDeflation(Mat M)
{
  PetscErrorCode   ierr;
  NEP_DEF_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  if (matctx->extop->szd) {
    ierr = BVDestroy(&matctx->U);CHKERRQ(ierr);
    ierr = PetscFree2(matctx->hfj,matctx->work);CHKERRQ(ierr);
    ierr = PetscFree2(matctx->A,matctx->B);CHKERRQ(ierr);
    ierr = VecDestroy(&matctx->w[0]);CHKERRQ(ierr);
    ierr = VecDestroy(&matctx->w[1]);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&matctx->T);CHKERRQ(ierr);
  ierr = PetscFree(matctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationEvaluateBasis(NEP_EXT_OP extop,PetscScalar lambda,PetscInt n,PetscScalar *val,PetscBool jacobian)
{
  PetscScalar p;
  PetscInt    i;

  PetscFunctionBegin;
  if (!jacobian) {
    val[0] = 1.0;
    for (i=1;i<extop->n;i++) val[i] = val[i-1]*(lambda-extop->bc[i-1]);
  } else {
    val[0] = 0.0;
    p = 1.0;
    for (i=1;i<extop->n;i++) {
      val[i] = val[i-1]*(lambda-extop->bc[i-1])+p;
      p *= (lambda-extop->bc[i-1]);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPDeflationComputeShellMat(NEP_EXT_OP extop,PetscScalar lambda,PetscBool jacobian,Mat *M)
{
  PetscErrorCode   ierr;
  NEP_DEF_MATSHELL *matctx,*matctxC;
  PetscInt         nloc,mloc,n=extop->n,j,i,szd=extop->szd,ldh=szd+1,k;
  Mat              F;
  Mat              Mshell,Mcomp;
  PetscBool        ini=PETSC_FALSE;
  PetscScalar      *hf,*hfj,*hfjp,sone=1.0,*hH,*hHprev,*pts,*B,*A,*Hj=extop->Hj,*basisv,zero=0.0;
  PetscBLASInt     n_,info,szd_;

  PetscFunctionBegin;
  if (!M) {
    Mshell = jacobian?extop->MJ:extop->MF;
  } else Mshell = *M;
  Mcomp  = jacobian?extop->MF:extop->MJ;
  if (!Mshell) {
    ini = PETSC_TRUE;
    ierr = PetscNew(&matctx);CHKERRQ(ierr);
    ierr = MatGetLocalSize(extop->nep->function,&mloc,&nloc);CHKERRQ(ierr);
    nloc += szd; mloc += szd;
    ierr = MatCreateShell(PetscObjectComm((PetscObject)extop->nep),nloc,mloc,PETSC_DETERMINE,PETSC_DETERMINE,matctx,&Mshell);CHKERRQ(ierr);
    ierr = MatShellSetOperation(Mshell,MATOP_MULT,(void(*)(void))MatMult_NEPDeflation);CHKERRQ(ierr);
    ierr = MatShellSetOperation(Mshell,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_NEPDeflation);CHKERRQ(ierr);
    ierr = MatShellSetOperation(Mshell,MATOP_DESTROY,(void(*)(void))MatDestroy_NEPDeflation);CHKERRQ(ierr);
    matctx->nep = extop->nep;
    matctx->extop = extop;
    if (!M) {
      if (jacobian) { matctx->jacob = PETSC_TRUE; matctx->T = extop->nep->jacobian; extop->MJ = Mshell; }
      else { matctx->jacob = PETSC_FALSE; matctx->T = extop->nep->function; extop->MF = Mshell; }
      ierr = PetscObjectReference((PetscObject)matctx->T);CHKERRQ(ierr);
    } else {
      matctx->jacob = jacobian;
      ierr = MatDuplicate(jacobian?extop->nep->jacobian:extop->nep->function,MAT_DO_NOT_COPY_VALUES, &matctx->T);CHKERRQ(ierr);
      *M = Mshell;
    }
    if (szd) {
      ierr = BVCreateVec(extop->nep->V,matctx->w);CHKERRQ(ierr);
      ierr = VecDuplicate(matctx->w[0],matctx->w+1);CHKERRQ(ierr);
      ierr = BVDuplicateResize(extop->nep->V,szd,&matctx->U);CHKERRQ(ierr);
      ierr = PetscMalloc2(extop->simpU?2*(szd)*(szd):2*(szd)*(szd)*extop->nep->nt,&matctx->hfj,szd,&matctx->work);CHKERRQ(ierr);
      ierr = PetscMalloc2(szd*szd,&matctx->A,szd*szd,&matctx->B);CHKERRQ(ierr);
    }
  } else {
    ierr = MatShellGetContext(Mshell,(void**)&matctx);CHKERRQ(ierr);
  }
  if (ini || matctx->theta != lambda || matctx->n != extop->n) {
    if (ini || matctx->theta != lambda) {
      if (jacobian) {
        ierr = NEPComputeJacobian(extop->nep,lambda,matctx->T);CHKERRQ(ierr);
      } else {
        ierr = NEPComputeFunction(extop->nep,lambda,matctx->T,matctx->T);CHKERRQ(ierr);
      }
    }
    if (n) {
      matctx->hfjset = PETSC_FALSE;
      if (!extop->simpU) {
        /* likely hfjp has been already computed */
        if (Mcomp) {
          ierr = MatShellGetContext(Mcomp,(void**)&matctxC);CHKERRQ(ierr);
          if (matctxC->hfjset && matctxC->theta == lambda && matctxC->n == extop->n) {
            ierr = PetscArraycpy(matctx->hfj,matctxC->hfj,2*extop->szd*extop->szd*extop->nep->nt);CHKERRQ(ierr);
            matctx->hfjset = PETSC_TRUE;
          }
        }
        hfj = matctx->hfj; hfjp = matctx->hfj+extop->szd*extop->szd*extop->nep->nt;
        if (!matctx->hfjset) {
          ierr = NEPDeflationEvaluateHatFunction(extop,-1,lambda,NULL,hfj,hfjp,n);CHKERRQ(ierr);
          matctx->hfjset = PETSC_TRUE;
        }
        ierr = BVSetActiveColumns(matctx->U,0,n);CHKERRQ(ierr);
        hf = jacobian?hfjp:hfj;
        ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,hf,&F);CHKERRQ(ierr);
        ierr = BVMatMult(extop->X,extop->nep->A[0],matctx->U);CHKERRQ(ierr);
        ierr = BVMultInPlace(matctx->U,F,0,n);CHKERRQ(ierr);
        ierr = BVSetActiveColumns(extop->W,0,extop->n);CHKERRQ(ierr);
        for (j=1;j<extop->nep->nt;j++) {
          ierr = BVMatMult(extop->X,extop->nep->A[j],extop->W);CHKERRQ(ierr);
          ierr = MatDensePlaceArray(F,hf+j*n*n);CHKERRQ(ierr);
          ierr = BVMult(matctx->U,1.0,1.0,extop->W,F);CHKERRQ(ierr);
          ierr = MatDenseResetArray(F);CHKERRQ(ierr);
        }
        ierr = MatDestroy(&F);CHKERRQ(ierr);
      } else {
        hfj = matctx->hfj;
        ierr = BVSetActiveColumns(matctx->U,0,n);CHKERRQ(ierr);
        ierr = BVMatMult(extop->X,matctx->T,matctx->U);CHKERRQ(ierr);
        for (j=0;j<n;j++) {
          for (i=0;i<n;i++) hfj[j*n+i] = -extop->H[j*ldh+i];
          hfj[j*(n+1)] += lambda;
        }
        ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
        PetscStackCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&n_,hfj,&n_,&info));
        SlepcCheckLapackInfo("trtri",info);
        ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,hfj,&F);CHKERRQ(ierr);
        ierr = BVMultInPlace(matctx->U,F,0,n);CHKERRQ(ierr);
        if (jacobian) {
          ierr = NEPDeflationComputeFunction(extop,lambda,NULL);CHKERRQ(ierr);
          ierr = MatShellGetContext(extop->MF,(void**)&matctxC);CHKERRQ(ierr);
          ierr = BVMult(matctx->U,-1.0,1.0,matctxC->U,F);CHKERRQ(ierr);
        }
        ierr = MatDestroy(&F);CHKERRQ(ierr);
      }
      ierr = PetscCalloc3(n,&basisv,szd*szd,&hH,szd*szd,&hHprev);CHKERRQ(ierr);
      ierr = NEPDeflationEvaluateBasis(extop,lambda,n,basisv,jacobian);CHKERRQ(ierr);
      A = matctx->A;
      ierr = PetscArrayzero(A,szd*szd);CHKERRQ(ierr);
      if (!jacobian) for (i=0;i<n;i++) A[i*(szd+1)] = 1.0;
      for (j=0;j<n;j++)
        for (i=0;i<n;i++)
          for (k=1;k<extop->midx;k++) A[j*szd+i] += basisv[k]*PetscConj(Hj[k*szd*szd+i*szd+j]);
      ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(szd,&szd_);CHKERRQ(ierr);
      B = matctx->B;
      ierr = PetscArrayzero(B,szd*szd);CHKERRQ(ierr);
      for (i=1;i<extop->midx;i++) {
        ierr = NEPDeflationEvaluateBasisMat(extop,i,PETSC_TRUE,basisv,hH,hHprev);CHKERRQ(ierr);
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,extop->XpX,&szd_,hH,&szd_,&zero,hHprev,&szd_));
        PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,extop->Hj+szd*szd*i,&szd_,hHprev,&szd_,&sone,B,&szd_));
        pts = hHprev; hHprev = hH; hH = pts;
      }
      ierr = PetscFree3(basisv,hH,hHprev);CHKERRQ(ierr);
    }
    matctx->theta = lambda;
    matctx->n = extop->n;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationComputeFunction(NEP_EXT_OP extop,PetscScalar lambda,Mat *F)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = NEPDeflationComputeShellMat(extop,lambda,PETSC_FALSE,NULL);CHKERRQ(ierr);
  if (F) *F = extop->MF;
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationComputeJacobian(NEP_EXT_OP extop,PetscScalar lambda,Mat *J)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = NEPDeflationComputeShellMat(extop,lambda,PETSC_TRUE,NULL);CHKERRQ(ierr);
  if (J) *J = extop->MJ;
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationSolveSetUp(NEP_EXT_OP extop,PetscScalar lambda)
{
  PetscErrorCode    ierr;
  NEP_DEF_MATSHELL  *matctx;
  NEP_DEF_FUN_SOLVE solve;
  PetscInt          i,j,n=extop->n;
  Vec               u,tu;
  Mat               F;
  PetscScalar       snone=-1.0,sone=1.0;
  PetscBLASInt      n_,szd_,ldh_,*p,info;
  Mat               Mshell;

  PetscFunctionBegin;
  solve = extop->solve;
  if (lambda!=solve->theta || n!=solve->n) {
    ierr = NEPDeflationComputeShellMat(extop,lambda,PETSC_FALSE,solve->sincf?NULL:&solve->T);CHKERRQ(ierr);
    Mshell = (solve->sincf)?extop->MF:solve->T;
    ierr = MatShellGetContext(Mshell,(void**)&matctx);CHKERRQ(ierr);
    ierr = KSPSetOperators(solve->ksp,matctx->T,matctx->T);CHKERRQ(ierr);
    if (!extop->ref && n) {
      ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(extop->szd,&szd_);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(extop->szd+1,&ldh_);CHKERRQ(ierr);
      if (!extop->simpU) {
        ierr = BVSetActiveColumns(solve->T_1U,0,n);CHKERRQ(ierr);
        for (i=0;i<n;i++) {
          ierr = BVGetColumn(matctx->U,i,&u);CHKERRQ(ierr);
          ierr = BVGetColumn(solve->T_1U,i,&tu);CHKERRQ(ierr);
          ierr = KSPSolve(solve->ksp,u,tu);CHKERRQ(ierr);
          ierr = BVRestoreColumn(solve->T_1U,i,&tu);CHKERRQ(ierr);
          ierr = BVRestoreColumn(matctx->U,i,&u);CHKERRQ(ierr);
        }
        ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,solve->work,&F);CHKERRQ(ierr);
        ierr = BVDot(solve->T_1U,extop->X,F);CHKERRQ(ierr);
        ierr = MatDestroy(&F);CHKERRQ(ierr);
      } else {
        for (j=0;j<n;j++)
          for (i=0;i<n;i++) solve->work[j*n+i] = extop->XpX[j*extop->szd+i];
        for (i=0;i<n;i++) extop->H[i*ldh_+i] -= lambda;
        PetscStackCallBLAS("BLAStrsm",BLAStrsm_("R","U","N","N",&n_,&n_,&snone,extop->H,&ldh_,solve->work,&n_));
        for (i=0;i<n;i++) extop->H[i*ldh_+i] += lambda;
      }
      ierr = PetscArraycpy(solve->M,matctx->B,extop->szd*extop->szd);CHKERRQ(ierr);
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&snone,matctx->A,&szd_,solve->work,&n_,&sone,solve->M,&szd_));
      ierr = PetscMalloc1(n,&p);CHKERRQ(ierr);
      PetscStackCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n_,&n_,solve->M,&szd_,p,&info));
      SlepcCheckLapackInfo("getrf",info);
      PetscStackCallBLAS("LAPACKgetri",LAPACKgetri_(&n_,solve->M,&szd_,p,solve->work,&n_,&info));
      SlepcCheckLapackInfo("getri",info);
      ierr = PetscFree(p);CHKERRQ(ierr);
    }
    solve->theta = lambda;
    solve->n = n;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationFunctionSolve(NEP_EXT_OP extop,Vec b,Vec x)
{
  PetscErrorCode    ierr;
  Vec               b1,x1;
  PetscScalar       *xx,*bb,*x2,*b2,*w,*w2,snone=-1.0,sone=1.0,zero=0.0;
  NEP_DEF_MATSHELL  *matctx;
  NEP_DEF_FUN_SOLVE solve=extop->solve;
  PetscBLASInt      one=1,szd_,n_,ldh_;
  PetscInt          nloc,i;
  PetscMPIInt       np,count;

  PetscFunctionBegin;
  if (extop->ref) {
    ierr = VecZeroEntries(x);CHKERRQ(ierr);
  }
  if (extop->szd) {
    x1 = solve->w[0]; b1 = solve->w[1];
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    ierr = VecPlaceArray(x1,xx);CHKERRQ(ierr);
    ierr = VecGetArray(b,&bb);CHKERRQ(ierr);
    ierr = VecPlaceArray(b1,bb);CHKERRQ(ierr);
  } else {
    b1 = b; x1 = x;
  }
  ierr = KSPSolve(extop->solve->ksp,b1,x1);CHKERRQ(ierr);
  if (!extop->ref && extop->n) {
    ierr = PetscBLASIntCast(extop->szd,&szd_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(extop->n,&n_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(extop->szd+1,&ldh_);CHKERRQ(ierr);
    ierr = BVGetSizes(extop->nep->V,&nloc,NULL,NULL);CHKERRQ(ierr);
    ierr = PetscMalloc2(extop->n,&b2,extop->n,&x2);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)b),&np);CHKERRQ(ierr);
    for (i=0;i<extop->n;i++) b2[i] = bb[nloc+i]*PetscSqrtReal(np);
    w = solve->work; w2 = solve->work+extop->n;
    ierr = MatShellGetContext(solve->sincf?extop->MF:solve->T,(void**)&matctx);CHKERRQ(ierr);
    ierr = PetscArraycpy(w2,b2,extop->n);CHKERRQ(ierr);
    ierr = BVDotVec(extop->X,x1,w);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&snone,matctx->A,&szd_,w,&one,&sone,w2,&one));
    PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,solve->M,&szd_,w2,&one,&zero,x2,&one));
    if (extop->simpU) {
      for (i=0;i<extop->n;i++) extop->H[i+i*(extop->szd+1)] -= solve->theta;
      for (i=0;i<extop->n;i++) w[i] = x2[i];
      PetscStackCallBLAS("BLAStrsm",BLAStrsm_("L","U","N","N",&n_,&one,&snone,extop->H,&ldh_,w,&n_));
      for (i=0;i<extop->n;i++) extop->H[i+i*(extop->szd+1)] += solve->theta;
      ierr = BVMultVec(extop->X,-1.0,1.0,x1,w);CHKERRQ(ierr);
    } else {
      ierr = BVMultVec(solve->T_1U,-1.0,1.0,x1,x2);CHKERRQ(ierr);
    }
    for (i=0;i<extop->n;i++) xx[i+nloc] = x2[i]/PetscSqrtReal(np);
    ierr = PetscMPIIntCast(extop->n,&count);CHKERRQ(ierr);
    ierr = MPI_Bcast(xx+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)b));CHKERRQ(ierr);
  }
  if (extop->szd) {
    ierr = VecResetArray(x1);CHKERRQ(ierr);
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
    ierr = VecResetArray(b1);CHKERRQ(ierr);
    ierr = VecRestoreArray(b,&bb);CHKERRQ(ierr);
    if (!extop->ref && extop->n) { ierr = PetscFree2(b2,x2);CHKERRQ(ierr);}
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationSetRefine(NEP_EXT_OP extop,PetscBool ref)
{
  PetscFunctionBegin;
  extop->ref = ref;
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationReset(NEP_EXT_OP extop)
{
  PetscErrorCode    ierr;
  PetscInt          j;
  NEP_DEF_FUN_SOLVE solve;

  PetscFunctionBegin;
  if (!extop) PetscFunctionReturn(0);
  ierr = PetscFree(extop->H);CHKERRQ(ierr);
  ierr = BVDestroy(&extop->X);CHKERRQ(ierr);
  if (extop->szd) {
    ierr = VecDestroy(&extop->w);CHKERRQ(ierr);
    ierr = PetscFree3(extop->Hj,extop->XpX,extop->bc);CHKERRQ(ierr);
    ierr = BVDestroy(&extop->W);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&extop->MF);CHKERRQ(ierr);
  ierr = MatDestroy(&extop->MJ);CHKERRQ(ierr);
  if (extop->solve) {
    solve = extop->solve;
    if (extop->szd) {
      if (!extop->simpU) {ierr = BVDestroy(&solve->T_1U);CHKERRQ(ierr);}
      ierr = PetscFree2(solve->M,solve->work);CHKERRQ(ierr);
      ierr = VecDestroy(&solve->w[0]);CHKERRQ(ierr);
      ierr = VecDestroy(&solve->w[1]);CHKERRQ(ierr);
    }
    if (!solve->sincf) {
      ierr = MatDestroy(&solve->T);CHKERRQ(ierr);
    }
    ierr = PetscFree(extop->solve);CHKERRQ(ierr);
  }
  if (extop->proj) {
    if (extop->szd) {
      for (j=0;j<extop->nep->nt;j++) {ierr = MatDestroy(&extop->proj->V1pApX[j]);CHKERRQ(ierr);}
      ierr = MatDestroy(&extop->proj->XpV1);CHKERRQ(ierr);
      ierr = PetscFree3(extop->proj->V2,extop->proj->V1pApX,extop->proj->work);CHKERRQ(ierr);
      ierr = VecDestroy(&extop->proj->w);CHKERRQ(ierr);
      ierr = BVDestroy(&extop->proj->V1);CHKERRQ(ierr);
    }
    ierr = PetscFree(extop->proj);CHKERRQ(ierr);
  }
  ierr = PetscFree(extop);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationInitialize(NEP nep,BV X,KSP ksp,PetscBool sincfun,PetscInt sz,NEP_EXT_OP *extop)
{
  PetscErrorCode    ierr;
  NEP_EXT_OP        op;
  NEP_DEF_FUN_SOLVE solve;
  PetscInt          szd;
  Vec               x;

  PetscFunctionBegin;
  ierr = NEPDeflationReset(*extop);CHKERRQ(ierr);
  ierr = PetscNew(&op);CHKERRQ(ierr);
  *extop  = op;
  op->nep = nep;
  op->n   = 0;
  op->szd = szd = sz-1;
  op->max_midx = PetscMin(MAX_MINIDX,szd);
  op->X = X;
  if (!X) { ierr = BVDuplicateResize(nep->V,sz,&op->X);CHKERRQ(ierr); }
  else { ierr = PetscObjectReference((PetscObject)X);CHKERRQ(ierr); }
  ierr = PetscCalloc1(sz*sz,&(op)->H);CHKERRQ(ierr);
  if (op->szd) {
    ierr = BVGetColumn(op->X,0,&x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&op->w);CHKERRQ(ierr);
    ierr = BVRestoreColumn(op->X,0,&x);CHKERRQ(ierr);
    op->simpU = PETSC_FALSE;
    if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
      /* undocumented option to use the simple expression for U = T*X*inv(lambda-H) */
      ierr = PetscOptionsGetBool(NULL,NULL,"-nep_deflation_simpleu",&op->simpU,NULL);CHKERRQ(ierr);
    } else {
      op->simpU = PETSC_TRUE;
    }
    ierr = PetscCalloc3(szd*szd*op->max_midx,&(op)->Hj,szd*szd,&(op)->XpX,szd,&op->bc);CHKERRQ(ierr);
    ierr = BVDuplicateResize(op->X,op->szd,&op->W);CHKERRQ(ierr);
  }
  if (ksp) {
    ierr = PetscNew(&solve);CHKERRQ(ierr);
    op->solve    = solve;
    solve->ksp   = ksp;
    solve->sincf = sincfun;
    solve->n     = -1;
    if (op->szd) {
      if (!op->simpU) {
        ierr = BVDuplicateResize(nep->V,szd,&solve->T_1U);CHKERRQ(ierr);
      }
      ierr = PetscMalloc2(szd*szd,&solve->M,2*szd*szd,&solve->work);CHKERRQ(ierr);
      ierr = BVCreateVec(nep->V,&solve->w[0]);CHKERRQ(ierr);
      ierr = VecDuplicate(solve->w[0],&solve->w[1]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationDSNEPComputeMatrix(DS ds,PetscScalar lambda,PetscBool deriv,DSMatType mat,void *ctx)
{
  PetscScalar     *T,*E,*w1,*w2,*w=NULL,*ww,*hH,*hHprev,*pts;
  PetscScalar     alpha,alpha2,*AB,sone=1.0,zero=0.0,*basisv,s;
  PetscInt        i,ldds,nwork=0,szd,nv,j,k,n;
  PetscBLASInt    inc=1,nv_,ldds_,dim_,dim2,szdk,szd_,n_,ldh_;
  PetscMPIInt     np;
  NEP_DEF_PROJECT proj=(NEP_DEF_PROJECT)ctx;
  NEP_EXT_OP      extop=proj->extop;
  NEP             nep=extop->nep;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = DSGetDimensions(ds,&nv,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(ds,&ldds);CHKERRQ(ierr);
  ierr = DSGetArray(ds,mat,&T);CHKERRQ(ierr);
  ierr = PetscArrayzero(T,ldds*nv);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldds*nv,&dim2);CHKERRQ(ierr);
  /* mat = V1^*T(lambda)V1 */
  for (i=0;i<nep->nt;i++) {
    if (deriv) {
      ierr = FNEvaluateDerivative(nep->f[i],lambda,&alpha);CHKERRQ(ierr);
    } else {
      ierr = FNEvaluateFunction(nep->f[i],lambda,&alpha);CHKERRQ(ierr);
    }
    ierr = DSGetArray(ds,DSMatExtra[i],&E);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&dim2,&alpha,E,&inc,T,&inc));
    ierr = DSRestoreArray(ds,DSMatExtra[i],&E);CHKERRQ(ierr);
  }
  if (!extop->ref && extop->n) {
    n = extop->n;
    szd = extop->szd;
    ierr = PetscArrayzero(proj->work,proj->lwork);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(nv,&nv_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(ldds,&ldds_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(szd,&szd_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(proj->dim,&dim_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(extop->szd+1,&ldh_);CHKERRQ(ierr);
    w1 = proj->work; w2 = proj->work+proj->dim*proj->dim;
    nwork += 2*proj->dim*proj->dim;

    /* mat = mat + V1^*U(lambda)V2 */
    for (i=0;i<nep->nt;i++) {
      ierr = MatDenseGetArray(proj->V1pApX[i],&E);CHKERRQ(ierr);
      if (extop->simpU) {
        if (deriv) {
          ierr = FNEvaluateDerivative(nep->f[i],lambda,&alpha);CHKERRQ(ierr);
        } else {
          ierr = FNEvaluateFunction(nep->f[i],lambda,&alpha);CHKERRQ(ierr);
        }
        ww = w1; w = w2;
        ierr = PetscArraycpy(ww,proj->V2,szd*nv);CHKERRQ(ierr);
        ierr = MPI_Comm_size(PetscObjectComm((PetscObject)ds),&np);CHKERRQ(ierr);
        for (j=0;j<szd*nv;j++) ww[j] *= PetscSqrtReal(np);
        for (j=0;j<n;j++) extop->H[j*ldh_+j] -= lambda;
        alpha = -alpha;
        PetscStackCallBLAS("BLAStrsm",BLAStrsm_("L","U","N","N",&n_,&nv_,&alpha,extop->H,&ldh_,ww,&szd_));
        if (deriv) {
          ierr = PetscBLASIntCast(szd*nv,&szdk);CHKERRQ(ierr);
          ierr = FNEvaluateFunction(nep->f[i],lambda,&alpha2);CHKERRQ(ierr);
          ierr = PetscArraycpy(w,proj->V2,szd*nv);CHKERRQ(ierr);
          for (j=0;j<szd*nv;j++) w[j] *= PetscSqrtReal(np);
          alpha2 = -alpha2;
          PetscStackCallBLAS("BLAStrsm",BLAStrsm_("L","U","N","N",&n_,&nv_,&alpha2,extop->H,&ldh_,w,&szd_));
          alpha2 = 1.0;
          PetscStackCallBLAS("BLAStrsm",BLAStrsm_("L","U","N","N",&n_,&nv_,&alpha2,extop->H,&ldh_,w,&szd_));
          PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&szdk,&sone,w,&inc,ww,&inc));
        }
        for (j=0;j<n;j++) extop->H[j*ldh_+j] += lambda;
      } else {
        ierr = NEPDeflationEvaluateHatFunction(extop,i,lambda,NULL,w1,w2,szd);CHKERRQ(ierr);
        w = deriv?w2:w1; ww = deriv?w1:w2;
        ierr = MPI_Comm_size(PetscObjectComm((PetscObject)ds),&np);CHKERRQ(ierr);
        s = PetscSqrtReal(np);
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&nv_,&n_,&s,w,&szd_,proj->V2,&szd_,&zero,ww,&szd_));
      }
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&nv_,&nv_,&n_,&sone,E,&dim_,ww,&szd_,&sone,T,&ldds_));
      ierr = MatDenseRestoreArray(proj->V1pApX[i],&E);CHKERRQ(ierr);
    }

    /* mat = mat + V2^*A(lambda)V1 */
    basisv = proj->work+nwork; nwork += szd;
    hH     = proj->work+nwork; nwork += szd*szd;
    hHprev = proj->work+nwork; nwork += szd*szd;
    AB     = proj->work+nwork;
    ierr = NEPDeflationEvaluateBasis(extop,lambda,n,basisv,deriv);CHKERRQ(ierr);
    if (!deriv) for (i=0;i<n;i++) AB[i*(szd+1)] = 1.0;
    for (j=0;j<n;j++)
      for (i=0;i<n;i++)
        for (k=1;k<extop->midx;k++) AB[j*szd+i] += basisv[k]*PetscConj(extop->Hj[k*szd*szd+i*szd+j]);
    ierr = MatDenseGetArray(proj->XpV1,&E);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&nv_,&n_,&sone,AB,&szd_,E,&szd_,&zero,w,&szd_));
    PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&nv_,&nv_,&n_,&sone,proj->V2,&szd_,w,&szd_,&sone,T,&ldds_));
    ierr = MatDenseRestoreArray(proj->XpV1,&E);CHKERRQ(ierr);

    /* mat = mat + V2^*B(lambda)V2 */
    ierr = PetscArrayzero(AB,szd*szd);CHKERRQ(ierr);
    for (i=1;i<extop->midx;i++) {
      ierr = NEPDeflationEvaluateBasisMat(extop,i,PETSC_TRUE,basisv,hH,hHprev);CHKERRQ(ierr);
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,extop->XpX,&szd_,hH,&szd_,&zero,hHprev,&szd_));
      PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,extop->Hj+szd*szd*i,&szd_,hHprev,&szd_,&sone,AB,&szd_));
      pts = hHprev; hHprev = hH; hH = pts;
    }
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&nv_,&n_,&sone,AB,&szd_,proj->V2,&szd_,&zero,w,&szd_));
    PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&nv_,&nv_,&n_,&sone,proj->V2,&szd_,w,&szd_,&sone,T,&ldds_));
  }
  ierr = DSRestoreArray(ds,mat,&T);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDeflationProjectOperator(NEP_EXT_OP extop,BV Vext,DS ds,PetscInt j0,PetscInt j1)
{
  PetscErrorCode  ierr;
  PetscInt        k,j,n=extop->n,dim;
  Vec             v,ve;
  BV              V1;
  Mat             G;
  NEP             nep=extop->nep;
  NEP_DEF_PROJECT proj;

  PetscFunctionBegin;
  NEPCheckSplit(extop->nep,1);
  proj = extop->proj;
  if (!proj) {
    /* Initialize the projection data structure */
    ierr = PetscNew(&proj);CHKERRQ(ierr);
    extop->proj = proj;
    proj->extop = extop;
    ierr = BVGetSizes(Vext,NULL,NULL,&dim);CHKERRQ(ierr);
    proj->dim = dim;
    if (extop->szd) {
      proj->lwork = 3*dim*dim+2*extop->szd*extop->szd+extop->szd;
      ierr = PetscMalloc3(dim*extop->szd,&proj->V2,nep->nt,&proj->V1pApX,proj->lwork,&proj->work);CHKERRQ(ierr);
      for (j=0;j<nep->nt;j++) {
        ierr =  MatCreateSeqDense(PETSC_COMM_SELF,proj->dim,extop->szd,NULL,&proj->V1pApX[j]);CHKERRQ(ierr);
      }
      ierr =  MatCreateSeqDense(PETSC_COMM_SELF,extop->szd,proj->dim,NULL,&proj->XpV1);CHKERRQ(ierr);
      ierr = BVCreateVec(extop->X,&proj->w);CHKERRQ(ierr);
      ierr = BVDuplicateResize(extop->X,proj->dim,&proj->V1);CHKERRQ(ierr);
    }
    ierr = DSNEPSetComputeMatrixFunction(ds,NEPDeflationDSNEPComputeMatrix,(void*)proj);CHKERRQ(ierr);
  }

  /* Split Vext in V1 and V2 */
  if (extop->szd) {
    for (j=j0;j<j1;j++) {
      ierr = BVGetColumn(Vext,j,&ve);CHKERRQ(ierr);
      ierr = BVGetColumn(proj->V1,j,&v);CHKERRQ(ierr);
      ierr = NEPDeflationCopyToExtendedVec(extop,v,proj->V2+j*extop->szd,ve,PETSC_TRUE);CHKERRQ(ierr);
      ierr = BVRestoreColumn(proj->V1,j,&v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(Vext,j,&ve);CHKERRQ(ierr);
    }
    V1 = proj->V1;
  } else V1 = Vext;

  /* Compute matrices V1^* A_i V1 */
  ierr = BVSetActiveColumns(V1,j0,j1);CHKERRQ(ierr);
  for (k=0;k<nep->nt;k++) {
    ierr = DSGetMat(ds,DSMatExtra[k],&G);CHKERRQ(ierr);
    ierr = BVMatProject(V1,nep->A[k],V1,G);CHKERRQ(ierr);
    ierr = DSRestoreMat(ds,DSMatExtra[k],&G);CHKERRQ(ierr);
  }

  if (extop->n) {
    if (extop->szd) {
      /* Compute matrices V1^* A_i X  and V1^* X */
      ierr = BVSetActiveColumns(extop->W,0,n);CHKERRQ(ierr);
      for (k=0;k<nep->nt;k++) {
        ierr = BVMatMult(extop->X,nep->A[k],extop->W);CHKERRQ(ierr);
        ierr = BVDot(extop->W,V1,proj->V1pApX[k]);CHKERRQ(ierr);
      }
      ierr = BVDot(V1,extop->X,proj->XpV1);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

