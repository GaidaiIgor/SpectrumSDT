/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Common subroutines for all Krylov-type PEP solvers
*/

#include <slepc/private/pepimpl.h>    /*I "slepcpep.h" I*/
#include <slepcblaslapack.h>
#include "pepkrylov.h"

PetscErrorCode PEPExtractVectors_TOAR(PEP pep)
{
  PetscErrorCode ierr;
  PetscInt       i,j,nq,deg=pep->nmat-1,lds,idxcpy=0,ldds,k,ld;
  PetscScalar    *X,*er,*ei,*SS,*vals,*ivals,sone=1.0,szero=0.0,*yi,*yr,*tr,*ti,alpha,t,*S,*pS0;
  PetscBLASInt   k_,nq_,lds_,one=1,ldds_;
  PetscBool      flg;
  PetscReal      norm,max,factor=1.0;
  Vec            xr,xi,w[4];
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;
  Mat            S0,MS;

  PetscFunctionBegin;
  ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
  ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
  ierr = BVGetSizes(pep->V,NULL,NULL,&ld);CHKERRQ(ierr);
  ierr = BVGetActiveColumns(pep->V,NULL,&nq);CHKERRQ(ierr);
  k = pep->nconv;
  if (k==0) PetscFunctionReturn(0);
  lds = deg*ld;
  ierr = DSGetLeadingDimension(pep->ds,&ldds);CHKERRQ(ierr);
  ierr = PetscCalloc5(k,&er,k,&ei,nq*k,&SS,pep->nmat,&vals,pep->nmat,&ivals);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (flg) factor = pep->sfactor;
  for (i=0;i<k;i++) {
    er[i] = factor*pep->eigr[i];
    ei[i] = factor*pep->eigi[i];
  }
  ierr = STBackTransform(pep->st,k,er,ei);CHKERRQ(ierr);

  ierr = DSVectors(pep->ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);
  ierr = DSGetArray(pep->ds,DS_MAT_X,&X);CHKERRQ(ierr);

  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(nq,&nq_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldds,&ldds_);CHKERRQ(ierr);

  if (pep->extract==PEP_EXTRACT_NONE || pep->refine==PEP_REFINE_MULTIPLE) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&nq_,&k_,&k_,&sone,S,&lds_,X,&ldds_,&szero,SS,&nq_));
  } else {
    switch (pep->extract) {
    case PEP_EXTRACT_NONE:
      break;
    case PEP_EXTRACT_NORM:
      for (i=0;i<k;i++) {
        ierr = PEPEvaluateBasis(pep,er[i],ei[i],vals,ivals);CHKERRQ(ierr);
        max = 1.0;
        for (j=1;j<deg;j++) {
          norm = SlepcAbsEigenvalue(vals[j],ivals[j]);
          if (max < norm) { max = norm; idxcpy = j; }
        }
        PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+idxcpy*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
#if !defined(PETSC_USE_COMPLEX)
        if (PetscRealPart(ei[i])!=0.0) {
          i++;
          PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+idxcpy*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
        }
#endif
      }
      break;
    case PEP_EXTRACT_RESIDUAL:
      ierr = VecDuplicate(pep->work[0],&xr);CHKERRQ(ierr);
      ierr = VecDuplicate(pep->work[0],&w[0]);CHKERRQ(ierr);
      ierr = VecDuplicate(pep->work[0],&w[1]);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
      ierr = VecDuplicate(pep->work[0],&w[2]);CHKERRQ(ierr);
      ierr = VecDuplicate(pep->work[0],&w[3]);CHKERRQ(ierr);
      ierr = VecDuplicate(pep->work[0],&xi);CHKERRQ(ierr);
#else
      xi = NULL;
#endif
      for (i=0;i<k;i++) {
        max = 0.0;
        for (j=0;j<deg;j++) {
          PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+j*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
          ierr = BVMultVec(pep->V,1.0,0.0,xr,SS+i*nq);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
          PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+j*ld,&lds_,X+(i+1)*ldds,&one,&szero,SS+i*nq,&one));
          ierr = BVMultVec(pep->V,1.0,0.0,xi,SS+i*nq);CHKERRQ(ierr);
#endif
          ierr = PEPComputeResidualNorm_Private(pep,er[i],ei[i],xr,xi,w,&norm);CHKERRQ(ierr);
          if (norm>max) { max = norm; idxcpy=j; }
        }
        PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+idxcpy*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
#if !defined(PETSC_USE_COMPLEX)
        if (PetscRealPart(ei[i])!=0.0) {
          i++;
          PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&sone,S+idxcpy*ld,&lds_,X+i*ldds,&one,&szero,SS+i*nq,&one));
        }
#endif
      }
      ierr = VecDestroy(&xr);CHKERRQ(ierr);
      ierr = VecDestroy(&w[0]);CHKERRQ(ierr);
      ierr = VecDestroy(&w[1]);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
      ierr = VecDestroy(&w[2]);CHKERRQ(ierr);
      ierr = VecDestroy(&w[3]);CHKERRQ(ierr);
      ierr = VecDestroy(&xi);CHKERRQ(ierr);
#endif
      break;
    case PEP_EXTRACT_STRUCTURED:
      ierr = PetscMalloc2(k,&tr,k,&ti);CHKERRQ(ierr);
      for (i=0;i<k;i++) {
        t = 0.0;
        ierr = PEPEvaluateBasis(pep,er[i],ei[i],vals,ivals);CHKERRQ(ierr);
        yr = X+i*ldds; yi = NULL;
#if !defined(PETSC_USE_COMPLEX)
        if (ei[i]!=0.0) { yr = tr; yi = ti; }
#endif
        for (j=0;j<deg;j++) {
          alpha = PetscConj(vals[j]);
#if !defined(PETSC_USE_COMPLEX)
          if (ei[i]!=0.0) {
            ierr = PetscArrayzero(tr,k);CHKERRQ(ierr);
            PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&k_,&vals[j],X+i*ldds,&one,tr,&one));
            PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&k_,&ivals[j],X+(i+1)*ldds,&one,tr,&one));
            ierr = PetscArrayzero(ti,k);CHKERRQ(ierr);
            PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&k_,&vals[j],X+(i+1)*ldds,&one,ti,&one));
            alpha = -ivals[j];
            PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&k_,&alpha,X+i*ldds,&one,ti,&one));
            alpha = 1.0;
          }
#endif
          PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&alpha,S+j*ld,&lds_,yr,&one,&sone,SS+i*nq,&one));
          t += SlepcAbsEigenvalue(vals[j],ivals[j])*SlepcAbsEigenvalue(vals[j],ivals[j]);
          if (yi) {
            PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&nq_,&k_,&alpha,S+j*ld,&lds_,yi,&one,&sone,SS+(i+1)*nq,&one));
          }
        }
        t = 1.0/t;
        PetscStackCallBLAS("BLASscal",BLASscal_(&nq_,&t,SS+i*nq,&one));
        if (yi) {
          PetscStackCallBLAS("BLASscal",BLASscal_(&nq_,&t,SS+(i+1)*nq,&one));
          i++;
        }
      }
      ierr = PetscFree2(tr,ti);CHKERRQ(ierr);
      break;
    }
  }

  /* update vectors V = V*S */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,nq,k,NULL,&S0);CHKERRQ(ierr);
  ierr = MatDenseGetArray(S0,&pS0);CHKERRQ(ierr);
  for (i=0;i<k;i++) {
    ierr = PetscArraycpy(pS0+i*nq,SS+i*nq,nq);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(S0,&pS0);CHKERRQ(ierr);
  ierr = BVMultInPlace(pep->V,S0,0,k);CHKERRQ(ierr);
  ierr = MatDestroy(&S0);CHKERRQ(ierr);
  ierr = PetscFree5(er,ei,SS,vals,ivals);CHKERRQ(ierr);
  if (ctx->V) {
    ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
    ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   PEPKrylovConvergence - This is the analogue to EPSKrylovConvergence, but
   for polynomial Krylov methods.

   Differences:
   - Always non-symmetric
   - Does not check for STSHIFT
   - No correction factor
   - No support for true residual
*/
PetscErrorCode PEPKrylovConvergence(PEP pep,PetscBool getall,PetscInt kini,PetscInt nits,PetscReal beta,PetscInt *kout)
{
  PetscErrorCode ierr;
  PetscInt       k,newk,marker,inside;
  PetscScalar    re,im;
  PetscReal      resnorm;
  PetscBool      istrivial;

  PetscFunctionBegin;
  ierr = RGIsTrivial(pep->rg,&istrivial);CHKERRQ(ierr);
  marker = -1;
  if (pep->trackall) getall = PETSC_TRUE;
  for (k=kini;k<kini+nits;k++) {
    /* eigenvalue */
    re = pep->eigr[k];
    im = pep->eigi[k];
    if (!istrivial) {
      ierr = STBackTransform(pep->st,1,&re,&im);CHKERRQ(ierr);
      ierr = RGCheckInside(pep->rg,1,&re,&im,&inside);CHKERRQ(ierr);
      if (marker==-1 && inside<0) marker = k;
      re = pep->eigr[k];
      im = pep->eigi[k];
    }
    newk = k;
    ierr = DSVectors(pep->ds,DS_MAT_X,&newk,&resnorm);CHKERRQ(ierr);
    resnorm *= beta;
    /* error estimate */
    ierr = (*pep->converged)(pep,re,im,resnorm,&pep->errest[k],pep->convergedctx);CHKERRQ(ierr);
    if (marker==-1 && pep->errest[k] >= pep->tol) marker = k;
    if (newk==k+1) {
      pep->errest[k+1] = pep->errest[k];
      k++;
    }
    if (marker!=-1 && !getall) break;
  }
  if (marker!=-1) k = marker;
  *kout = k;
  PetscFunctionReturn(0);
}

