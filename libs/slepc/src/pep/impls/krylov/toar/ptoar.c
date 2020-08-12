/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc polynomial eigensolver: "toar"

   Method: TOAR

   Algorithm:

       Two-Level Orthogonal Arnoldi.

   References:

       [1] Y. Su, J. Zhang and Z. Bai, "A compact Arnoldi algorithm for
           polynomial eigenvalue problems", talk presented at RANMEP 2008.

       [2] C. Campos and J.E. Roman, "Parallel Krylov solvers for the
           polynomial eigenvalue problem in SLEPc", SIAM J. Sci. Comput.
           38(5):S385-S411, 2016.

       [3] D. Lu, Y. Su and Z. Bai, "Stability analysis of the two-level
           orthogonal Arnoldi procedure", SIAM J. Matrix Anal. App.
           37(1):195-214, 2016.
*/

#include <slepc/private/pepimpl.h>    /*I "slepcpep.h" I*/
#include "../src/pep/impls/krylov/pepkrylov.h"
#include <slepcblaslapack.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-pep,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"Parallel {Krylov} solvers for the polynomial eigenvalue problem in {SLEPc}\",\n"
  "   journal = \"{SIAM} J. Sci. Comput.\",\n"
  "   volume = \"38\",\n"
  "   number = \"5\",\n"
  "   pages = \"S385--S411\",\n"
  "   year = \"2016,\"\n"
  "   doi = \"https://doi.org/10.1137/15M1022458\"\n"
  "}\n";

PetscErrorCode PEPSetUp_TOAR(PEP pep)
{
  PetscErrorCode ierr;
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;
  PetscBool      shift,sinv,flg;
  PetscInt       i;

  PetscFunctionBegin;
  pep->lineariz = PETSC_TRUE;
  ierr = PEPSetDimensions_Default(pep,pep->nev,&pep->ncv,&pep->mpd);CHKERRQ(ierr);
  if (!ctx->lock && pep->mpd<pep->ncv) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Should not use mpd parameter in non-locking variant");
  if (pep->max_it==PETSC_DEFAULT) pep->max_it = PetscMax(100,2*(pep->nmat-1)*pep->n/pep->ncv);

  ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSHIFT,&shift);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv);CHKERRQ(ierr);
  if (!shift && !sinv) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Only STSHIFT and STSINVERT spectral transformations can be used");
  if (!pep->which) {
    if (sinv) pep->which = PEP_TARGET_MAGNITUDE;
    else pep->which = PEP_LARGEST_MAGNITUDE;
  }
  if (pep->which==PEP_ALL) SETERRQ(PetscObjectComm((PetscObject)pep),1,"Wrong value of pep->which");
  if (pep->problem_type!=PEP_GENERAL) {
    ierr = PetscInfo(pep,"Problem type ignored, performing a non-symmetric linearization\n");CHKERRQ(ierr);
  }

  if (!ctx->keep) ctx->keep = 0.5;

  ierr = PEPAllocateSolution(pep,pep->nmat-1);CHKERRQ(ierr);
  ierr = PEPSetWorkVecs(pep,3);CHKERRQ(ierr);
  ierr = DSSetType(pep->ds,DSNHEP);CHKERRQ(ierr);
  ierr = DSSetExtraRow(pep->ds,PETSC_TRUE);CHKERRQ(ierr);
  ierr = DSAllocate(pep->ds,pep->ncv+1);CHKERRQ(ierr);

  ierr = PEPBasisCoefficients(pep,pep->pbc);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (!flg) {
    ierr = PetscMalloc1(pep->nmat,&pep->solvematcoeffs);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)pep,pep->nmat*sizeof(PetscScalar));CHKERRQ(ierr);
    if (sinv) {
      ierr = PEPEvaluateBasis(pep,pep->target,0,pep->solvematcoeffs,NULL);CHKERRQ(ierr);
    } else {
      for (i=0;i<pep->nmat-1;i++) pep->solvematcoeffs[i] = 0.0;
      pep->solvematcoeffs[pep->nmat-1] = 1.0;
    }
  }
  ierr = BVDestroy(&ctx->V);CHKERRQ(ierr);
  ierr = BVCreateTensor(pep->V,pep->nmat-1,&ctx->V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Extend the TOAR basis by applying the the matrix operator
  over a vector which is decomposed in the TOAR way
  Input:
    - pbc: array containing the polynomial basis coefficients
    - S,V: define the latest Arnoldi vector (nv vectors in V)
  Output:
    - t: new vector extending the TOAR basis
    - r: temporary coefficients to compute the TOAR coefficients
         for the new Arnoldi vector
  Workspace: t_ (two vectors)
*/
static PetscErrorCode PEPTOARExtendBasis(PEP pep,PetscBool sinvert,PetscScalar sigma,PetscScalar *S,PetscInt ls,PetscInt nv,BV V,Vec t,PetscScalar *r,PetscInt lr,Vec *t_)
{
  PetscErrorCode ierr;
  PetscInt       nmat=pep->nmat,deg=nmat-1,k,j,off=0,lss;
  Vec            v=t_[0],ve=t_[1],q=t_[2];
  PetscScalar    alpha=1.0,*ss,a;
  PetscReal      *ca=pep->pbc,*cb=pep->pbc+nmat,*cg=pep->pbc+2*nmat;
  PetscBool      flg;

  PetscFunctionBegin;
  ierr = BVSetActiveColumns(pep->V,0,nv);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (sinvert) {
    for (j=0;j<nv;j++) {
      if (deg>1) r[lr+j] = S[j]/ca[0];
      if (deg>2) r[2*lr+j] = (S[ls+j]+(sigma-cb[1])*r[lr+j])/ca[1];
    }
    for (k=2;k<deg-1;k++) {
      for (j=0;j<nv;j++) r[(k+1)*lr+j] = (S[k*ls+j]+(sigma-cb[k])*r[k*lr+j]-cg[k]*r[(k-1)*lr+j])/ca[k];
    }
    k = deg-1;
    for (j=0;j<nv;j++) r[j] = (S[k*ls+j]+(sigma-cb[k])*r[k*lr+j]-cg[k]*r[(k-1)*lr+j])/ca[k];
    ss = r; lss = lr; off = 1; alpha = -1.0; a = pep->sfactor;
  } else {
    ss = S; lss = ls; off = 0; alpha = -ca[deg-1]; a = 1.0;
  }
  ierr = BVMultVec(V,1.0,0.0,v,ss+off*lss);CHKERRQ(ierr);
  if (pep->Dr) { /* balancing */
    ierr = VecPointwiseMult(v,v,pep->Dr);CHKERRQ(ierr);
  }
  ierr = STMatMult(pep->st,off,v,q);CHKERRQ(ierr);
  ierr = VecScale(q,a);CHKERRQ(ierr);
  for (j=1+off;j<deg+off-1;j++) {
    ierr = BVMultVec(V,1.0,0.0,v,ss+j*lss);CHKERRQ(ierr);
    if (pep->Dr) {
      ierr = VecPointwiseMult(v,v,pep->Dr);CHKERRQ(ierr);
    }
    ierr = STMatMult(pep->st,j,v,t);CHKERRQ(ierr);
    a *= pep->sfactor;
    ierr = VecAXPY(q,a,t);CHKERRQ(ierr);
  }
  if (sinvert) {
    ierr = BVMultVec(V,1.0,0.0,v,ss);CHKERRQ(ierr);
    if (pep->Dr) {
      ierr = VecPointwiseMult(v,v,pep->Dr);CHKERRQ(ierr);
    }
    ierr = STMatMult(pep->st,deg,v,t);CHKERRQ(ierr);
    a *= pep->sfactor;
    ierr = VecAXPY(q,a,t);CHKERRQ(ierr);
  } else {
    ierr = BVMultVec(V,1.0,0.0,ve,ss+(deg-1)*lss);CHKERRQ(ierr);
    if (pep->Dr) {
      ierr = VecPointwiseMult(ve,ve,pep->Dr);CHKERRQ(ierr);
    }
    a *= pep->sfactor;
    ierr = STMatMult(pep->st,deg-1,ve,t);CHKERRQ(ierr);
    ierr = VecAXPY(q,a,t);CHKERRQ(ierr);
    a *= pep->sfactor;
  }
  if (flg || !sinvert) alpha /= a;
  ierr = STMatSolve(pep->st,q,t);CHKERRQ(ierr);
  ierr = VecScale(t,alpha);CHKERRQ(ierr);
  if (!sinvert) {
    if (cg[deg-1]!=0) { ierr = VecAXPY(t,cg[deg-1],v);CHKERRQ(ierr); }
    if (cb[deg-1]!=0) { ierr = VecAXPY(t,cb[deg-1],ve);CHKERRQ(ierr); }
  }
  if (pep->Dr) {
    ierr = VecPointwiseDivide(t,t,pep->Dr);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  Compute TOAR coefficients of the blocks of the new Arnoldi vector computed
*/
static PetscErrorCode PEPTOARCoefficients(PEP pep,PetscBool sinvert,PetscScalar sigma,PetscInt nv,PetscScalar *S,PetscInt ls,PetscScalar *r,PetscInt lr,PetscScalar *x)
{
  PetscInt    k,j,nmat=pep->nmat,d=nmat-1;
  PetscReal   *ca=pep->pbc,*cb=pep->pbc+nmat,*cg=pep->pbc+2*nmat;
  PetscScalar t=1.0,tp=0.0,tt;

  PetscFunctionBegin;
  if (sinvert) {
    for (k=1;k<d;k++) {
      tt = t;
      t = ((sigma-cb[k-1])*t-cg[k-1]*tp)/ca[k-1]; /* k-th basis polynomial */
      tp = tt;
      for (j=0;j<=nv;j++) r[k*lr+j] += t*x[j];
    }
  } else {
    for (j=0;j<=nv;j++) r[j] = (cb[0]-sigma)*S[j]+ca[0]*S[ls+j];
    for (k=1;k<d-1;k++) {
      for (j=0;j<=nv;j++) r[k*lr+j] = (cb[k]-sigma)*S[k*ls+j]+ca[k]*S[(k+1)*ls+j]+cg[k]*S[(k-1)*ls+j];
    }
    if (sigma!=0.0) for (j=0;j<=nv;j++) r[(d-1)*lr+j] -= sigma*S[(d-1)*ls+j];
  }
  PetscFunctionReturn(0);
}

/*
  Compute a run of Arnoldi iterations dim(work)=ld
*/
static PetscErrorCode PEPTOARrun(PEP pep,PetscScalar sigma,PetscScalar *H,PetscInt ldh,PetscInt k,PetscInt *M,PetscBool *breakdown,Vec *t_)
{
  PetscErrorCode ierr;
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;
  PetscInt       j,m=*M,deg=pep->nmat-1,ld;
  PetscInt       lds,nqt,l;
  Vec            t;
  PetscReal      norm;
  PetscBool      flg,sinvert=PETSC_FALSE,lindep;
  PetscScalar    *x,*S;
  Mat            MS;

  PetscFunctionBegin;
  ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
  ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
  ierr = BVGetSizes(pep->V,NULL,NULL,&ld);CHKERRQ(ierr);
  lds = ld*deg;
  ierr = BVGetActiveColumns(pep->V,&l,&nqt);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (!flg) {
    /* spectral transformation handled by the solver */
    ierr = PetscObjectTypeCompareAny((PetscObject)pep->st,&flg,STSINVERT,STSHIFT,"");CHKERRQ(ierr);
    if (!flg) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"ST type not supported for TOAR without transforming matrices");
    ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinvert);CHKERRQ(ierr);
  }
  ierr = BVSetActiveColumns(ctx->V,0,m);CHKERRQ(ierr);
  for (j=k;j<m;j++) {
    /* apply operator */
    ierr = BVGetColumn(pep->V,nqt,&t);CHKERRQ(ierr);
    ierr = PEPTOARExtendBasis(pep,sinvert,sigma,S+j*lds,ld,nqt,pep->V,t,S+(j+1)*lds,ld,t_);CHKERRQ(ierr);
    ierr = BVRestoreColumn(pep->V,nqt,&t);CHKERRQ(ierr);

    /* orthogonalize */
    if (sinvert) x = S+(j+1)*lds;
    else x = S+(deg-1)*ld+(j+1)*lds;
    ierr = BVOrthogonalizeColumn(pep->V,nqt,x,&norm,&lindep);CHKERRQ(ierr);
    if (!lindep) {
      x[nqt] = norm;
      ierr = BVScaleColumn(pep->V,nqt,1.0/norm);CHKERRQ(ierr);
      nqt++;
    }

    ierr = PEPTOARCoefficients(pep,sinvert,sigma,nqt-1,S+j*lds,ld,S+(j+1)*lds,ld,x);CHKERRQ(ierr);

    /* level-2 orthogonalization */
    ierr = BVOrthogonalizeColumn(ctx->V,j+1,H+j*ldh,&norm,breakdown);CHKERRQ(ierr);
    H[j+1+ldh*j] = norm;
    if (*breakdown) {
      *M = j+1;
      break;
    }
    ierr = BVScaleColumn(ctx->V,j+1,1.0/norm);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(pep->V,l,nqt);CHKERRQ(ierr);
  }
  ierr = BVSetActiveColumns(ctx->V,0,*M);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
  ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Computes T_j = phi_idx(T). In T_j and T_p are phi_{idx-1}(T)
   and phi_{idx-2}(T) respectively or null if idx=0,1.
   Tp and Tj are input/output arguments
*/
static PetscErrorCode PEPEvaluateBasisM(PEP pep,PetscInt k,PetscScalar *T,PetscInt ldt,PetscInt idx,PetscScalar **Tp,PetscScalar **Tj)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      *ca,*cb,*cg;
  PetscScalar    *pt,g,a;
  PetscBLASInt   k_,ldt_;

  PetscFunctionBegin;
  if (idx==0) {
    ierr = PetscArrayzero(*Tj,k*k);CHKERRQ(ierr);
    ierr = PetscArrayzero(*Tp,k*k);CHKERRQ(ierr);
    for (i=0;i<k;i++) (*Tj)[i+i*k] = 1.0;
  } else {
    ierr = PetscBLASIntCast(ldt,&ldt_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
    ca = pep->pbc; cb = pep->pbc+pep->nmat; cg = pep->pbc+2*pep->nmat;
    for (i=0;i<k;i++) T[i*ldt+i] -= cb[idx-1];
    a = 1/ca[idx-1];
    g = (idx==1)?0.0:-cg[idx-1]/ca[idx-1];
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&a,T,&ldt_,*Tj,&k_,&g,*Tp,&k_));
    pt = *Tj; *Tj = *Tp; *Tp = pt;
    for (i=0;i<k;i++) T[i*ldt+i] += cb[idx-1];
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPExtractInvariantPair(PEP pep,PetscScalar sigma,PetscInt sr,PetscInt k,PetscScalar *S,PetscInt ld,PetscInt deg,PetscScalar *H,PetscInt ldh)
{
  PetscErrorCode ierr;
  PetscInt       i,j,jj,lds,ldt,d=pep->nmat-1,idxcpy=0;
  PetscScalar    *At,*Bt,*Hj,*Hp,*T,sone=1.0,g,a,*pM,*work;
  PetscBLASInt   k_,sr_,lds_,ldh_,info,*p,lwork,ldt_;
  PetscBool      transf=PETSC_FALSE,flg;
  PetscReal      nrm,norm,maxnrm,*rwork;
  BV             *R,Y;
  Mat            M,*A;
  Vec            v;

  PetscFunctionBegin;
  if (k==0) PetscFunctionReturn(0);
  lds = deg*ld;
  ierr = PetscCalloc6(k,&p,sr*k,&At,k*k,&Bt,k*k,&Hj,k*k,&Hp,sr*k,&work);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(sr,&sr_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldh,&ldh_);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (!flg) {
    ierr =  PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&flg);CHKERRQ(ierr);
    if (flg || sigma!=0.0) transf=PETSC_TRUE;
  }
  if (transf) {
    ierr = PetscMalloc1(k*k,&T);CHKERRQ(ierr);
    ldt = k;
    for (i=0;i<k;i++) {
      ierr = PetscArraycpy(T+k*i,H+i*ldh,k);CHKERRQ(ierr);
    }
    if (flg) {
      PetscStackCallBLAS("LAPACKgetrf",LAPACKgetrf_(&k_,&k_,T,&k_,p,&info));
      SlepcCheckLapackInfo("getrf",info);
      ierr = PetscBLASIntCast(sr*k,&lwork);CHKERRQ(ierr);
      PetscStackCallBLAS("LAPACKgetri",LAPACKgetri_(&k_,T,&k_,p,work,&lwork,&info));
      SlepcCheckLapackInfo("getri",info);
    }
    if (sigma!=0.0) for (i=0;i<k;i++) T[i+k*i] += sigma;
  } else {
    T = H; ldt = ldh;
  }
  ierr = PetscBLASIntCast(ldt,&ldt_);CHKERRQ(ierr);
  switch (pep->extract) {
  case PEP_EXTRACT_NONE:
    break;
  case PEP_EXTRACT_NORM:
    if (pep->basis == PEP_BASIS_MONOMIAL) {
      ierr = PetscBLASIntCast(ldt,&ldt_);CHKERRQ(ierr);
      ierr = PetscMalloc1(k,&rwork);CHKERRQ(ierr);
      norm = LAPACKlange_("F",&k_,&k_,T,&ldt_,rwork);
      PetscFree(rwork);CHKERRQ(ierr);
      if (norm>1.0) idxcpy = d-1;
    } else {
      ierr = PetscBLASIntCast(ldt,&ldt_);CHKERRQ(ierr);
      ierr = PetscMalloc1(k,&rwork);CHKERRQ(ierr);
      maxnrm = 0.0;
      for (i=0;i<pep->nmat-1;i++) {
        ierr = PEPEvaluateBasisM(pep,k,T,ldt,i,&Hp,&Hj);CHKERRQ(ierr);
        norm = LAPACKlange_("F",&k_,&k_,Hj,&k_,rwork);
        if (norm > maxnrm) {
          idxcpy = i;
          maxnrm = norm;
        }
      }
      PetscFree(rwork);CHKERRQ(ierr);
    }
    if (idxcpy>0) {
      /* copy block idxcpy of S to the first one */
      for (j=0;j<k;j++) {
        ierr = PetscArraycpy(S+j*lds,S+idxcpy*ld+j*lds,sr);CHKERRQ(ierr);
      }
    }
    break;
  case PEP_EXTRACT_RESIDUAL:
    ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = PetscMalloc1(pep->nmat,&A);CHKERRQ(ierr);
      for (i=0;i<pep->nmat;i++) {
        ierr = STGetMatrixTransformed(pep->st,i,A+i);CHKERRQ(ierr);
      }
    } else A = pep->A;
    ierr = PetscMalloc1(pep->nmat-1,&R);CHKERRQ(ierr);
    for (i=0;i<pep->nmat-1;i++) {
      ierr = BVDuplicateResize(pep->V,k,R+i);CHKERRQ(ierr);
    }
    ierr = BVDuplicateResize(pep->V,sr,&Y);CHKERRQ(ierr);
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,sr,k,NULL,&M);CHKERRQ(ierr);
    g = 0.0; a = 1.0;
    ierr = BVSetActiveColumns(pep->V,0,sr);CHKERRQ(ierr);
    for (j=0;j<pep->nmat;j++) {
      ierr = BVMatMult(pep->V,A[j],Y);CHKERRQ(ierr);
      ierr = PEPEvaluateBasisM(pep,k,T,ldt,i,&Hp,&Hj);CHKERRQ(ierr);
      for (i=0;i<pep->nmat-1;i++) {
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&sr_,&k_,&k_,&a,S+i*ld,&lds_,Hj,&k_,&g,At,&sr_));
        ierr = MatDenseGetArray(M,&pM);CHKERRQ(ierr);
        for (jj=0;jj<k;jj++) {
          ierr = PetscArraycpy(pM+jj*sr,At+jj*sr,sr);CHKERRQ(ierr);
        }
        ierr = MatDenseRestoreArray(M,&pM);CHKERRQ(ierr);
        ierr = BVMult(R[i],1.0,(i==0)?0.0:1.0,Y,M);CHKERRQ(ierr);
      }
    }

    /* frobenius norm */
    maxnrm = 0.0;
    for (i=0;i<pep->nmat-1;i++) {
      norm = 0.0;
      for (j=0;j<k;j++) {
        ierr = BVGetColumn(R[i],j,&v);CHKERRQ(ierr);
        ierr = VecNorm(v,NORM_2,&nrm);CHKERRQ(ierr);
        ierr = BVRestoreColumn(R[i],j,&v);CHKERRQ(ierr);
        norm += nrm*nrm;
      }
      norm = PetscSqrtReal(norm);
      if (maxnrm > norm) {
        maxnrm = norm;
        idxcpy = i;
      }
    }
    if (idxcpy>0) {
      /* copy block idxcpy of S to the first one */
      for (j=0;j<k;j++) {
        ierr = PetscArraycpy(S+j*lds,S+idxcpy*ld+j*lds,sr);CHKERRQ(ierr);
      }
    }
    if (flg) PetscFree(A);CHKERRQ(ierr);
    for (i=0;i<pep->nmat-1;i++) {
      ierr = BVDestroy(&R[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree(R);CHKERRQ(ierr);
    ierr = BVDestroy(&Y);CHKERRQ(ierr);
    ierr = MatDestroy(&M);CHKERRQ(ierr);
    break;
  case PEP_EXTRACT_STRUCTURED:
    for (j=0;j<k;j++) Bt[j+j*k] = 1.0;
    for (j=0;j<sr;j++) {
      for (i=0;i<k;i++) At[j*k+i] = PetscConj(S[i*lds+j]);
    }
    ierr = PEPEvaluateBasisM(pep,k,T,ldt,0,&Hp,&Hj);CHKERRQ(ierr);
    for (i=1;i<deg;i++) {
      ierr = PEPEvaluateBasisM(pep,k,T,ldt,i,&Hp,&Hj);CHKERRQ(ierr);
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&k_,&sr_,&k_,&sone,Hj,&k_,S+i*ld,&lds_,&sone,At,&k_));
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","C",&k_,&k_,&k_,&sone,Hj,&k_,Hj,&k_,&sone,Bt,&k_));
    }
    PetscStackCallBLAS("LAPACKgesv",LAPACKgesv_(&k_,&sr_,Bt,&k_,p,At,&k_,&info));
    SlepcCheckLapackInfo("gesv",info);
    for (j=0;j<sr;j++) {
      for (i=0;i<k;i++) S[i*lds+j] = PetscConj(At[j*k+i]);
    }
    break;
  }
  if (transf) { ierr = PetscFree(T);CHKERRQ(ierr); }
  ierr = PetscFree6(p,At,Bt,Hj,Hp,work);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSolve_TOAR(PEP pep)
{
  PetscErrorCode ierr;
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;
  PetscInt       i,j,k,l,nv=0,ld,lds,ldds,newn,nq=0,nconv=0;
  PetscInt       nmat=pep->nmat,deg=nmat-1;
  PetscScalar    *S,*H,sigma;
  PetscReal      beta;
  PetscBool      breakdown=PETSC_FALSE,flg,falselock=PETSC_FALSE,sinv=PETSC_FALSE;
  Mat            MS,MQ;

  PetscFunctionBegin;
  ierr = PetscCitationsRegister(citation,&cited);CHKERRQ(ierr);
  if (ctx->lock) {
    /* undocumented option to use a cheaper locking instead of the true locking */
    ierr = PetscOptionsGetBool(NULL,NULL,"-pep_toar_falselocking",&falselock,NULL);CHKERRQ(ierr);
  }
  ierr = DSGetLeadingDimension(pep->ds,&ldds);CHKERRQ(ierr);
  ierr = STGetShift(pep->st,&sigma);CHKERRQ(ierr);

  /* update polynomial basis coefficients */
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (pep->sfactor!=1.0) {
    for (i=0;i<nmat;i++) {
      pep->pbc[nmat+i] /= pep->sfactor;
      pep->pbc[2*nmat+i] /= pep->sfactor*pep->sfactor;
    }
    if (!flg) {
      pep->target /= pep->sfactor;
      ierr = RGPushScale(pep->rg,1.0/pep->sfactor);CHKERRQ(ierr);
      ierr = STScaleShift(pep->st,1.0/pep->sfactor);CHKERRQ(ierr);
      sigma /= pep->sfactor;
    } else {
      ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv);CHKERRQ(ierr);
      pep->target = sinv?pep->target*pep->sfactor:pep->target/pep->sfactor;
      ierr = RGPushScale(pep->rg,sinv?pep->sfactor:1.0/pep->sfactor);CHKERRQ(ierr);
      ierr = STScaleShift(pep->st,sinv?pep->sfactor:1.0/pep->sfactor);CHKERRQ(ierr);
    }
  }

  if (flg) sigma = 0.0;

  /* clean projected matrix (including the extra-arrow) */
  ierr = DSGetArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
  ierr = PetscArrayzero(H,ldds*ldds);CHKERRQ(ierr);
  ierr = DSRestoreArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);

  /* Get the starting Arnoldi vector */
  ierr = BVTensorBuildFirstColumn(ctx->V,pep->nini);CHKERRQ(ierr);

  /* restart loop */
  l = 0;
  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;

    /* compute an nv-step Lanczos factorization */
    nv = PetscMax(PetscMin(nconv+pep->mpd,pep->ncv),nv);
    ierr = DSGetArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
    ierr = PEPTOARrun(pep,sigma,H,ldds,pep->nconv+l,&nv,&breakdown,pep->work);CHKERRQ(ierr);
    beta = PetscAbsScalar(H[(nv-1)*ldds+nv]);
    ierr = DSRestoreArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
    ierr = DSSetDimensions(pep->ds,nv,0,pep->nconv,pep->nconv+l);CHKERRQ(ierr);
    if (l==0) {
      ierr = DSSetState(pep->ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
    } else {
      ierr = DSSetState(pep->ds,DS_STATE_RAW);CHKERRQ(ierr);
    }
    ierr = BVSetActiveColumns(ctx->V,pep->nconv,nv);CHKERRQ(ierr);

    /* solve projected problem */
    ierr = DSSolve(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);
    ierr = DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DSUpdateExtraRow(pep->ds);CHKERRQ(ierr);
    ierr = DSSynchronize(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);

    /* check convergence */
    ierr = PEPKrylovConvergence(pep,PETSC_FALSE,pep->nconv,nv-pep->nconv,beta,&k);CHKERRQ(ierr);
    ierr = (*pep->stopping)(pep,pep->its,pep->max_it,k,pep->nev,&pep->reason,pep->stoppingctx);CHKERRQ(ierr);

    /* update l */
    if (pep->reason != PEP_CONVERGED_ITERATING || breakdown) l = 0;
    else {
      l = (nv==k)?0:PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
      if (!breakdown) {
        /* prepare the Rayleigh quotient for restart */
        ierr = DSTruncate(pep->ds,k+l);CHKERRQ(ierr);
        ierr = DSGetDimensions(pep->ds,&newn,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
        l = newn-k;
      }
    }
    nconv = k;
    if (!ctx->lock && pep->reason == PEP_CONVERGED_ITERATING && !breakdown) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */

    /* update S */
    ierr = DSGetMat(pep->ds,DS_MAT_Q,&MQ);CHKERRQ(ierr);
    ierr = BVMultInPlace(ctx->V,MQ,pep->nconv,k+l);CHKERRQ(ierr);
    ierr = MatDestroy(&MQ);CHKERRQ(ierr);

    /* copy last column of S */
    ierr = BVCopyColumn(ctx->V,nv,k+l);CHKERRQ(ierr);

    if (breakdown && pep->reason == PEP_CONVERGED_ITERATING) {
      /* stop if breakdown */
      ierr = PetscInfo2(pep,"Breakdown TOAR method (it=%D norm=%g)\n",pep->its,(double)beta);CHKERRQ(ierr);
      pep->reason = PEP_DIVERGED_BREAKDOWN;
    }
    if (pep->reason != PEP_CONVERGED_ITERATING) l--;
    /* truncate S */
    ierr = BVGetActiveColumns(pep->V,NULL,&nq);CHKERRQ(ierr);
    if (k+l+deg<=nq) {
      ierr = BVSetActiveColumns(ctx->V,pep->nconv,k+l+1);CHKERRQ(ierr);
      if (!falselock && ctx->lock) {
        ierr = BVTensorCompress(ctx->V,k-pep->nconv);CHKERRQ(ierr);
      } else {
        ierr = BVTensorCompress(ctx->V,0);CHKERRQ(ierr);
      }
    }
    pep->nconv = k;
    ierr = PEPMonitor(pep,pep->its,nconv,pep->eigr,pep->eigi,pep->errest,nv);CHKERRQ(ierr);
  }
  if (pep->nconv>0) {
    /* {V*S_nconv^i}_{i=0}^{d-1} has rank nconv instead of nconv+d-1. Force zeros in each S_nconv^i block */
    ierr = BVSetActiveColumns(ctx->V,0,pep->nconv);CHKERRQ(ierr);
    ierr = BVGetActiveColumns(pep->V,NULL,&nq);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(pep->V,0,nq);CHKERRQ(ierr);
    if (nq>pep->nconv) {
      ierr = BVTensorCompress(ctx->V,pep->nconv);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(pep->V,0,pep->nconv);CHKERRQ(ierr);
      nq = pep->nconv;
    }

    /* perform Newton refinement if required */
    if (pep->refine==PEP_REFINE_MULTIPLE && pep->rits>0) {
      /* extract invariant pair */
      ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
      ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
      ierr = DSGetArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
      ierr = BVGetSizes(pep->V,NULL,NULL,&ld);CHKERRQ(ierr);
      lds = deg*ld;
      ierr = PEPExtractInvariantPair(pep,sigma,nq,pep->nconv,S,ld,deg,H,ldds);CHKERRQ(ierr);
      ierr = DSRestoreArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
      ierr = DSSetDimensions(pep->ds,pep->nconv,0,0,0);CHKERRQ(ierr);
      ierr = DSSetState(pep->ds,DS_STATE_RAW);CHKERRQ(ierr);
      ierr = PEPNewtonRefinement_TOAR(pep,sigma,&pep->rits,NULL,pep->nconv,S,lds);CHKERRQ(ierr);
      ierr = DSSolve(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);
      ierr = DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
      ierr = DSSynchronize(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);
      ierr = DSGetMat(pep->ds,DS_MAT_Q,&MQ);CHKERRQ(ierr);
      ierr = BVMultInPlace(ctx->V,MQ,0,pep->nconv);CHKERRQ(ierr);
      ierr = MatDestroy(&MQ);CHKERRQ(ierr);
      ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
      ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    }
  }
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (pep->refine!=PEP_REFINE_MULTIPLE || pep->rits==0) {
    if (!flg && pep->ops->backtransform) {
        ierr = (*pep->ops->backtransform)(pep);CHKERRQ(ierr);
    }
    if (pep->sfactor!=1.0) {
      for (j=0;j<pep->nconv;j++) {
        pep->eigr[j] *= pep->sfactor;
        pep->eigi[j] *= pep->sfactor;
      }
      /* restore original values */
      for (i=0;i<pep->nmat;i++){
        pep->pbc[pep->nmat+i] *= pep->sfactor;
        pep->pbc[2*pep->nmat+i] *= pep->sfactor*pep->sfactor;
      }
    }
  }
  /* restore original values */
  if (!flg) {
    pep->target *= pep->sfactor;
    ierr = STScaleShift(pep->st,pep->sfactor);CHKERRQ(ierr);
  } else {
    ierr = STScaleShift(pep->st,sinv?1.0/pep->sfactor:pep->sfactor);CHKERRQ(ierr);
    pep->target = (sinv)?pep->target/pep->sfactor:pep->target*pep->sfactor;
  }
  if (pep->sfactor!=1.0) { ierr = RGPopScale(pep->rg);CHKERRQ(ierr); }

  /* change the state to raw so that DSVectors() computes eigenvectors from scratch */
  ierr = DSSetDimensions(pep->ds,pep->nconv,0,0,0);CHKERRQ(ierr);
  ierr = DSSetState(pep->ds,DS_STATE_RAW);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPTOARSetRestart_TOAR(PEP pep,PetscReal keep)
{
  PEP_TOAR *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  if (keep==PETSC_DEFAULT) ctx->keep = 0.5;
  else {
    if (keep<0.1 || keep>0.9) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The keep argument must be in the range [0.1,0.9]");
    ctx->keep = keep;
  }
  PetscFunctionReturn(0);
}

/*@
   PEPTOARSetRestart - Sets the restart parameter for the TOAR
   method, in particular the proportion of basis vectors that must be kept
   after restart.

   Logically Collective on pep

   Input Parameters:
+  pep  - the eigenproblem solver context
-  keep - the number of vectors to be kept at restart

   Options Database Key:
.  -pep_toar_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,0.9]. The default is 0.5.

   Level: advanced

.seealso: PEPTOARGetRestart()
@*/
PetscErrorCode PEPTOARSetRestart(PEP pep,PetscReal keep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,keep,2);
  ierr = PetscTryMethod(pep,"PEPTOARSetRestart_C",(PEP,PetscReal),(pep,keep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPTOARGetRestart_TOAR(PEP pep,PetscReal *keep)
{
  PEP_TOAR *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  *keep = ctx->keep;
  PetscFunctionReturn(0);
}

/*@
   PEPTOARGetRestart - Gets the restart parameter used in the TOAR method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  keep - the restart parameter

   Level: advanced

.seealso: PEPTOARSetRestart()
@*/
PetscErrorCode PEPTOARGetRestart(PEP pep,PetscReal *keep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidRealPointer(keep,2);
  ierr = PetscUseMethod(pep,"PEPTOARGetRestart_C",(PEP,PetscReal*),(pep,keep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPTOARSetLocking_TOAR(PEP pep,PetscBool lock)
{
  PEP_TOAR *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(0);
}

/*@
   PEPTOARSetLocking - Choose between locking and non-locking variants of
   the TOAR method.

   Logically Collective on pep

   Input Parameters:
+  pep  - the eigenproblem solver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -pep_toar_locking - Sets the locking flag

   Notes:
   The default is to lock converged eigenpairs when the method restarts.
   This behaviour can be changed so that all directions are kept in the
   working subspace even if already converged to working accuracy (the
   non-locking variant).

   Level: advanced

.seealso: PEPTOARGetLocking()
@*/
PetscErrorCode PEPTOARSetLocking(PEP pep,PetscBool lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,lock,2);
  ierr = PetscTryMethod(pep,"PEPTOARSetLocking_C",(PEP,PetscBool),(pep,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPTOARGetLocking_TOAR(PEP pep,PetscBool *lock)
{
  PEP_TOAR *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(0);
}

/*@
   PEPTOARGetLocking - Gets the locking flag used in the TOAR method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: PEPTOARSetLocking()
@*/
PetscErrorCode PEPTOARGetLocking(PEP pep,PetscBool *lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidBoolPointer(lock,2);
  ierr = PetscUseMethod(pep,"PEPTOARGetLocking_C",(PEP,PetscBool*),(pep,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSetFromOptions_TOAR(PetscOptionItems *PetscOptionsObject,PEP pep)
{
  PetscErrorCode ierr;
  PetscBool      flg,lock;
  PetscReal      keep;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"PEP TOAR Options");CHKERRQ(ierr);

    ierr = PetscOptionsReal("-pep_toar_restart","Proportion of vectors kept after restart","PEPTOARSetRestart",0.5,&keep,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPTOARSetRestart(pep,keep);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-pep_toar_locking","Choose between locking and non-locking variants","PEPTOARSetLocking",PETSC_FALSE,&lock,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPTOARSetLocking(pep,lock);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPView_TOAR(PEP pep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  %d%% of basis vectors kept after restart\n",(int)(100*ctx->keep));CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  using the %slocking variant\n",ctx->lock?"":"non-");CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPDestroy_TOAR(PEP pep)
{
  PetscErrorCode ierr;
  PEP_TOAR       *ctx = (PEP_TOAR*)pep->data;

  PetscFunctionBegin;
  ierr = BVDestroy(&ctx->V);CHKERRQ(ierr);
  ierr = PetscFree(pep->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPTOARSetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPTOARGetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPTOARSetLocking_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPTOARGetLocking_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode PEPCreate_TOAR(PEP pep)
{
  PEP_TOAR       *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(pep,&ctx);CHKERRQ(ierr);
  pep->data = (void*)ctx;
  ctx->lock = PETSC_TRUE;

  pep->ops->solve          = PEPSolve_TOAR;
  pep->ops->setup          = PEPSetUp_TOAR;
  pep->ops->setfromoptions = PEPSetFromOptions_TOAR;
  pep->ops->destroy        = PEPDestroy_TOAR;
  pep->ops->view           = PEPView_TOAR;
  pep->ops->backtransform  = PEPBackTransform_Default;
  pep->ops->computevectors = PEPComputeVectors_Default;
  pep->ops->extractvectors = PEPExtractVectors_TOAR;

  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPTOARSetRestart_C",PEPTOARSetRestart_TOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPTOARGetRestart_C",PEPTOARGetRestart_TOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPTOARSetLocking_C",PEPTOARSetLocking_TOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPTOARGetLocking_C",PEPTOARGetLocking_TOAR);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

