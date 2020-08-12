/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "krylovschur"

   Method: Two-sided Arnoldi with Krylov-Schur restart (for left eigenvectors)

   References:

       [1] I.N. Zwaan and M.E. Hochstenbach, "Krylov-Schur-type restarts
           for the two-sided Arnoldi method", SIAM J. Matrix Anal. Appl.
           38(2):297-321, 2017.

*/

#include <slepc/private/epsimpl.h>
#include "krylovschur.h"
#include <slepcblaslapack.h>

static PetscErrorCode EPSTwoSidedRQUpdate1(EPS eps,Mat M,PetscInt nv)
{
  PetscErrorCode ierr;
  PetscScalar    *T,*S,*A,*w,*pM,beta;
  Vec            u;
  PetscInt       ld,ncv=eps->ncv,i,l,nnv;
  PetscBLASInt   info,n_,ncv_,*p,one=1;

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  ierr = PetscMalloc3(nv,&p,ncv*ncv,&A,ncv,&w);CHKERRQ(ierr);
  ierr = BVGetActiveColumns(eps->V,&l,&nnv);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,0,nv);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->W,0,nv);CHKERRQ(ierr);
  ierr = BVGetColumn(eps->V,nv,&u);CHKERRQ(ierr);
  ierr = BVDotVec(eps->W,u,w);CHKERRQ(ierr);
  ierr = BVRestoreColumn(eps->V,nv,&u);CHKERRQ(ierr);
  ierr = MatDenseGetArray(M,&pM);CHKERRQ(ierr);
  ierr = PetscArraycpy(A,pM,ncv*ncv);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(nv,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ncv,&ncv_);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n_,&n_,A,&ncv_,p,&info));
  SlepcCheckLapackInfo("getrf",info);
  ierr = PetscLogFlops(2.0*n_*n_*n_/3.0);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&n_,&one,A,&ncv_,p,w,&ncv_,&info));
  SlepcCheckLapackInfo("getrs",info);
  ierr = PetscLogFlops(2.0*n_*n_-n_);CHKERRQ(ierr);
  ierr = BVMultColumn(eps->V,-1.0,1.0,nv,w);CHKERRQ(ierr);
  ierr = DSGetArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
  beta = S[(nv-1)*ld+nv];
  for (i=0;i<nv;i++) S[(nv-1)*ld+i] += beta*w[i];
  ierr = DSRestoreArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
  ierr = BVGetColumn(eps->W,nv,&u);CHKERRQ(ierr);
  ierr = BVDotVec(eps->V,u,w);CHKERRQ(ierr);
  ierr = BVRestoreColumn(eps->W,nv,&u);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgetrs",LAPACKgetrs_("C",&n_,&one,A,&ncv_,p,w,&ncv_,&info));
  ierr = BVMultColumn(eps->W,-1.0,1.0,nv,w);CHKERRQ(ierr);
  ierr = DSGetArray(eps->dsts,DS_MAT_A,&T);CHKERRQ(ierr);
  beta = T[(nv-1)*ld+nv];
  for (i=0;i<nv;i++) T[(nv-1)*ld+i] += beta*w[i];
  ierr = DSRestoreArray(eps->dsts,DS_MAT_A,&T);CHKERRQ(ierr);
  ierr = PetscFree3(p,A,w);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,l,nnv);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->W,l,nnv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSTwoSidedRQUpdate2(EPS eps,Mat M,PetscInt k)
{
  PetscErrorCode ierr;
  PetscScalar    *Q,*pM,*w,zero=0.0,sone=1.0,*c,*A;
  PetscBLASInt   n_,ncv_,ld_;
  PetscReal      norm;
  PetscInt       l,nv,ncv=eps->ncv,ld,i,j;

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  ierr = BVGetActiveColumns(eps->V,&l,&nv);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,0,nv);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->W,0,nv);CHKERRQ(ierr);
  ierr = PetscMalloc2(ncv*ncv,&w,ncv,&c);CHKERRQ(ierr);
  /* u = u - V*V'*u */
  ierr = BVOrthogonalizeColumn(eps->V,k,c,&norm,NULL);CHKERRQ(ierr);
  ierr = BVScaleColumn(eps->V,k,1.0/norm);CHKERRQ(ierr);
  ierr = DSGetArray(eps->ds,DS_MAT_A,&A);CHKERRQ(ierr);
  /* H = H + V'*u*b' */
  for (j=l;j<k;j++) {
    for (i=0;i<k;i++) A[i+j*ld] += c[i]*A[k+j*ld];
    A[k+j*ld] *= norm;
  }
  ierr = DSRestoreArray(eps->ds,DS_MAT_A,&A);CHKERRQ(ierr);
  ierr = BVOrthogonalizeColumn(eps->W,k,c,&norm,NULL);CHKERRQ(ierr);
  ierr = BVScaleColumn(eps->W,k,1.0/norm);CHKERRQ(ierr);
  ierr = DSGetArray(eps->dsts,DS_MAT_A,&A);CHKERRQ(ierr);
  /* H = H + V'*u*b' */
  for (j=l;j<k;j++) {
    for (i=0;i<k;i++) A[i+j*ld] += c[i]*A[k+j*ld];
    A[k+j*ld] *= norm;
  }
  ierr = DSRestoreArray(eps->dsts,DS_MAT_A,&A);CHKERRQ(ierr);

  /* M = Q'*M*Q */
  ierr = MatDenseGetArray(M,&pM);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ncv,&ncv_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(nv,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ld,&ld_);CHKERRQ(ierr);
  ierr = DSGetArray(eps->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,pM,&ncv_,Q,&ld_,&zero,w,&ncv_));
  ierr = DSRestoreArray(eps->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
  ierr = DSGetArray(eps->dsts,DS_MAT_Q,&Q);CHKERRQ(ierr);
  PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,Q,&ld_,w,&ncv_,&zero,pM,&ncv_));
  ierr = DSRestoreArray(eps->dsts,DS_MAT_Q,&Q);CHKERRQ(ierr);
  ierr = PetscFree2(w,c);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,l,nv);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->W,l,nv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_KrylovSchur_TwoSided(EPS eps)
{
  PetscErrorCode  ierr;
  EPS_KRYLOVSCHUR *ctx = (EPS_KRYLOVSCHUR*)eps->data;
  Mat             M,U,Op,OpHT;
  PetscReal       norm,norm2,beta,betat,s,t;
  PetscScalar     *pM,*S,*T,*eigr,*eigi,*Q;
  PetscInt        ld,l,nv,nvt,ncv=eps->ncv,i,j,k,nconv,*p,cont,*idx,*idx2,id=0;
  PetscBool       breakdownt,breakdown,breakdownl;
#if defined(PETSC_USE_COMPLEX)
  Mat             A;
#endif

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  ierr = EPSGetStartVector(eps,0,NULL);CHKERRQ(ierr);
  ierr = EPSGetLeftStartVector(eps,0,NULL);CHKERRQ(ierr);
  l = 0;
  ierr = PetscMalloc6(ncv*ncv,&pM,ncv,&eigr,ncv,&eigi,ncv,&idx,ncv,&idx2,ncv,&p);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,eps->ncv,eps->ncv,pM,&M);CHKERRQ(ierr);

  ierr = STGetOperator(eps->st,&Op);CHKERRQ(ierr);
  ierr = MatCreateHermitianTranspose(Op,&OpHT);CHKERRQ(ierr);

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an nv-step Arnoldi factorization for Op */
    nv = PetscMin(eps->nconv+eps->mpd,eps->ncv);
    ierr = DSGetArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
    ierr = BVMatArnoldi(eps->V,Op,S,ld,eps->nconv+l,&nv,&beta,&breakdown);CHKERRQ(ierr);
    ierr = DSRestoreArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);

    /* Compute an nv-step Arnoldi factorization for Op' */
    nvt = nv;
    ierr = DSGetArray(eps->dsts,DS_MAT_A,&T);CHKERRQ(ierr);
    ierr = BVMatArnoldi(eps->W,OpHT,T,ld,eps->nconv+l,&nvt,&betat,&breakdownt);CHKERRQ(ierr);
    ierr = DSRestoreArray(eps->dsts,DS_MAT_A,&T);CHKERRQ(ierr);

    /* Make sure both factorizations have the same length */
    nv = PetscMin(nv,nvt);
    ierr = DSSetDimensions(eps->ds,nv,0,eps->nconv,eps->nconv+l);CHKERRQ(ierr);
    ierr = DSSetDimensions(eps->dsts,nv,0,eps->nconv,eps->nconv+l);CHKERRQ(ierr);
    if (l==0) {
      ierr = DSSetState(eps->ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
      ierr = DSSetState(eps->dsts,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
    } else {
      ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
      ierr = DSSetState(eps->dsts,DS_STATE_RAW);CHKERRQ(ierr);
    }
    breakdown = (breakdown || breakdownt)? PETSC_TRUE: PETSC_FALSE;

    /* Update M, modify Rayleigh quotients S and T */
    ierr = BVSetActiveColumns(eps->V,eps->nconv+l,nv);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(eps->W,eps->nconv+l,nv);CHKERRQ(ierr);
    ierr = BVMatProject(eps->V,NULL,eps->W,M);CHKERRQ(ierr);

    ierr = EPSTwoSidedRQUpdate1(eps,M,nv);CHKERRQ(ierr);

    /* Solve projected problem */
    ierr = DSSolve(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);
    ierr = DSSort(eps->ds,eps->eigr,eps->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DSSynchronize(eps->ds,eps->eigr,eps->eigi);CHKERRQ(ierr);
    ierr = DSSolve(eps->dsts,eigr,eigi);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = DSGetMat(eps->dsts,DS_MAT_A,&A);CHKERRQ(ierr);
    ierr = MatConjugate(A);CHKERRQ(ierr);
    ierr = DSRestoreMat(eps->dsts,DS_MAT_A,&A);CHKERRQ(ierr);
    ierr = DSGetMat(eps->dsts,DS_MAT_Q,&U);CHKERRQ(ierr);
    ierr = MatConjugate(U);CHKERRQ(ierr);
    ierr = DSRestoreMat(eps->dsts,DS_MAT_Q,&U);CHKERRQ(ierr);
    for (i=0;i<nv;i++) eigr[i] = PetscConj(eigr[i]);
#endif
    ierr = DSSort(eps->dsts,eigr,eigi,NULL,NULL,NULL);CHKERRQ(ierr);
    /* check correct eigenvalue correspondence */
    cont = 0;
    for (i=0;i<nv;i++) {
      if (SlepcAbsEigenvalue(eigr[i]-eps->eigr[i],eigi[i]-eps->eigi[i])>PETSC_SQRT_MACHINE_EPSILON) {idx2[cont] = i; idx[cont++] = i;}
      p[i] = -1;
    }
    if (cont) {
      for (i=0;i<cont;i++) {
        t = PETSC_MAX_REAL;
        for (j=0;j<cont;j++) if (idx2[j]!=-1 && (s=SlepcAbsEigenvalue(eigr[idx[j]]-eps->eigr[idx[i]],eigi[idx[j]]-eps->eigi[idx[i]]))<t) { id = j; t = s; }
        p[idx[i]] = idx[id];
        idx2[id] = -1;
      }
      for (i=0;i<nv;i++) if (p[i]==-1) p[i] = i;
      ierr = DSSortWithPermutation(eps->dsts,p,eigr,eigi);CHKERRQ(ierr);
    }
#if defined(PETSC_USE_COMPLEX)
    ierr = DSGetMat(eps->dsts,DS_MAT_A,&A);CHKERRQ(ierr);
    ierr = MatConjugate(A);CHKERRQ(ierr);
    ierr = DSRestoreMat(eps->dsts,DS_MAT_A,&A);CHKERRQ(ierr);
    ierr = DSGetMat(eps->dsts,DS_MAT_Q,&U);CHKERRQ(ierr);
    ierr = MatConjugate(U);CHKERRQ(ierr);
    ierr = DSRestoreMat(eps->dsts,DS_MAT_Q,&U);CHKERRQ(ierr);
#endif
    ierr = DSSynchronize(eps->dsts,eigr,eigi);CHKERRQ(ierr);

    /* Check convergence */
    ierr = BVNormColumn(eps->V,nv,NORM_2,&norm);CHKERRQ(ierr);
    ierr = BVNormColumn(eps->W,nv,NORM_2,&norm2);CHKERRQ(ierr);
    ierr = EPSKrylovConvergence(eps,PETSC_FALSE,eps->nconv,nv-eps->nconv,beta*norm,betat*norm2,1.0,&k);CHKERRQ(ierr);
    ierr = (*eps->stopping)(eps,eps->its,eps->max_it,k,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);
    nconv = k;

    /* Update l */
    if (eps->reason != EPS_CONVERGED_ITERATING || breakdown || k==nv) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
#if !defined(PETSC_USE_COMPLEX)
      ierr = DSGetArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
      if (S[k+l+(k+l-1)*ld] != 0.0) {
        if (k+l<nv-1) l = l+1;
        else l = l-1;
      }
      ierr = DSRestoreArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
#endif
    }
    if (!ctx->lock && l>0) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */

    /* Update the corresponding vectors V(:,idx) = V*Q(:,idx) */
    ierr = BVSetActiveColumns(eps->V,eps->nconv,nv);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(eps->W,eps->nconv,nv);CHKERRQ(ierr);
    ierr = DSGetMat(eps->ds,DS_MAT_Q,&U);CHKERRQ(ierr);
    ierr = BVMultInPlace(eps->V,U,eps->nconv,k+l);CHKERRQ(ierr);
    ierr = MatDestroy(&U);CHKERRQ(ierr);
    ierr = DSGetMat(eps->dsts,DS_MAT_Q,&U);CHKERRQ(ierr);
    ierr = BVMultInPlace(eps->W,U,eps->nconv,k+l);CHKERRQ(ierr);
    ierr = MatDestroy(&U);CHKERRQ(ierr);
    if (eps->reason == EPS_CONVERGED_ITERATING && !breakdown) {
      ierr = BVCopyColumn(eps->V,nv,k+l);CHKERRQ(ierr);
      ierr = BVCopyColumn(eps->W,nv,k+l);CHKERRQ(ierr);
    }

    if (eps->reason == EPS_CONVERGED_ITERATING) {
      if (breakdown || k==nv) {
        /* Start a new Arnoldi factorization */
        ierr = PetscInfo2(eps,"Breakdown in Krylov-Schur method (it=%D norm=%g)\n",eps->its,(double)beta);CHKERRQ(ierr);
        if (k<eps->nev) {
          ierr = EPSGetStartVector(eps,k,&breakdown);CHKERRQ(ierr);
          ierr = EPSGetLeftStartVector(eps,k,&breakdownl);CHKERRQ(ierr);
          if (breakdown || breakdownl) {
            eps->reason = EPS_DIVERGED_BREAKDOWN;
            ierr = PetscInfo(eps,"Unable to generate more start vectors\n");CHKERRQ(ierr);
          }
        }
      } else {
        /* Prepare the Rayleigh quotient for restart */
        ierr = DSGetArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
        ierr = DSGetArray(eps->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
        for (i=k;i<k+l;i++) S[k+l+i*ld] = Q[nv-1+i*ld]*beta;
        ierr = DSRestoreArray(eps->ds,DS_MAT_A,&S);CHKERRQ(ierr);
        ierr = DSRestoreArray(eps->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
        ierr = DSGetArray(eps->dsts,DS_MAT_A,&S);CHKERRQ(ierr);
        ierr = DSGetArray(eps->dsts,DS_MAT_Q,&Q);CHKERRQ(ierr);
        for (i=k;i<k+l;i++) S[k+l+i*ld] = Q[nv-1+i*ld]*betat;
        ierr = DSRestoreArray(eps->dsts,DS_MAT_A,&S);CHKERRQ(ierr);
        ierr = DSRestoreArray(eps->dsts,DS_MAT_Q,&Q);CHKERRQ(ierr);
      }
      ierr = EPSTwoSidedRQUpdate2(eps,M,k+l);CHKERRQ(ierr);
    }
    eps->nconv = k;
    ierr = EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,nv);CHKERRQ(ierr);
  }

  ierr = STRestoreOperator(eps->st,&Op);CHKERRQ(ierr);
  ierr = MatDestroy(&OpHT);CHKERRQ(ierr);

  /* truncate Schur decomposition and change the state to raw so that
     DSVectors() computes eigenvectors from scratch */
  ierr = DSSetDimensions(eps->ds,eps->nconv,0,0,0);CHKERRQ(ierr);
  ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
  ierr = DSSetDimensions(eps->dsts,eps->nconv,0,0,0);CHKERRQ(ierr);
  ierr = DSSetState(eps->dsts,DS_STATE_RAW);CHKERRQ(ierr);
  ierr = PetscFree6(pM,eigr,eigi,idx,idx2,p);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

