/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Newton refinement for polynomial eigenproblems.

   References:

       [1] T. Betcke and D. Kressner, "Perturbation, extraction and refinement
           of invariant pairs for matrix polynomials", Linear Algebra Appl.
           435(3):514-536, 2011.

       [2] C. Campos and J.E. Roman, "Parallel iterative refinement in
           polynomial eigenvalue problems", Numer. Linear Algebra Appl. 23(4):
           730-745, 2016.
*/

#include <slepc/private/pepimpl.h>
#include <slepcblaslapack.h>

typedef struct {
  Mat          *A,M1;
  BV           V,M2,M3,W;
  PetscInt     k,nmat;
  PetscScalar  *fih,*work,*M4;
  PetscBLASInt *pM4;
  PetscBool    compM1;
  Vec          t;
} PEP_REFINE_MATSHELL;

typedef struct {
  Mat          E[2],M1;
  Vec          tN,ttN,t1,vseq;
  VecScatter   scatterctx;
  PetscBool    compM1;
  PetscInt     *map0,*map1,*idxg,*idxp;
  PetscSubcomm subc;
  VecScatter   scatter_sub;
  VecScatter   *scatter_id,*scatterp_id;
  Mat          *A;
  BV           V,W,M2,M3,Wt;
  PetscScalar  *M4,*w,*wt,*d,*dt;
  Vec          t,tg,Rv,Vi,tp,tpg;
  PetscInt     idx,*cols;
} PEP_REFINE_EXPLICIT;

static PetscErrorCode MatMult_FS(Mat M ,Vec x,Vec y)
{
  PetscErrorCode      ierr;
  PEP_REFINE_MATSHELL *ctx;
  PetscInt            k,i;
  PetscScalar         *c;
  PetscBLASInt        k_,one=1,info;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,&ctx);CHKERRQ(ierr);
  ierr = VecCopy(x,ctx->t);CHKERRQ(ierr);
  k    = ctx->k;
  c    = ctx->work;
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = MatMult(ctx->M1,x,y);CHKERRQ(ierr);
  ierr = VecConjugate(ctx->t);CHKERRQ(ierr);
  ierr = BVDotVec(ctx->M3,ctx->t,c);CHKERRQ(ierr);
  for (i=0;i<k;i++) c[i] = PetscConj(c[i]);
  PetscStackCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&k_,&one,ctx->M4,&k_,ctx->pM4,c,&k_,&info));
  SlepcCheckLapackInfo("getrs",info);
  ierr = BVMultVec(ctx->M2,-1.0,1.0,y,c);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Evaluates the first d elements of the polynomial basis
  on a given matrix H which is considered to be triangular
*/
static PetscErrorCode PEPEvaluateBasisforMatrix(PEP pep,PetscInt nm,PetscInt k,PetscScalar *H,PetscInt ldh,PetscScalar *fH)
{
  PetscErrorCode ierr;
  PetscInt       i,j,ldfh=nm*k,off,nmat=pep->nmat;
  PetscReal      *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat,t;
  PetscScalar    corr=0.0,alpha,beta;
  PetscBLASInt   k_,ldh_,ldfh_;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ldh,&ldh_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldfh,&ldfh_);CHKERRQ(ierr);
  ierr = PetscArrayzero(fH,nm*k*k);CHKERRQ(ierr);
  if (nm>0) for (j=0;j<k;j++) fH[j+j*ldfh] = 1.0;
  if (nm>1) {
    t = b[0]/a[0];
    off = k;
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) fH[off+i+j*ldfh] = H[i+j*ldh]/a[0];
      fH[j+j*ldfh] -= t;
    }
  }
  for (i=2;i<nm;i++) {
    off = i*k;
    if (i==2) {
      for (j=0;j<k;j++) {
        fH[off+j+j*ldfh] = 1.0;
        H[j+j*ldh] -= b[1];
      }
    } else {
      for (j=0;j<k;j++) {
        ierr = PetscArraycpy(fH+off+j*ldfh,fH+(i-2)*k+j*ldfh,k);CHKERRQ(ierr);
        H[j+j*ldh] += corr-b[i-1];
      }
    }
    corr  = b[i-1];
    beta  = -g[i-1]/a[i-1];
    alpha = 1/a[i-1];
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&alpha,H,&ldh_,fH+(i-1)*k,&ldfh_,&beta,fH+off,&ldfh_));
  }
  for (j=0;j<k;j++) H[j+j*ldh] += corr;
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefSysSetup_shell(PEP pep,PetscInt k,PetscScalar *fH,PetscScalar *S,PetscInt lds,PetscScalar *fh,PetscScalar h,PEP_REFINE_MATSHELL *ctx)
{
  PetscErrorCode    ierr;
  PetscScalar       *DHii,*T12,*Tr,*Ts,*array,s,ss,sone=1.0,zero=0.0,*M4=ctx->M4,t,*v,*T;
  const PetscScalar *m3,*m2;
  PetscInt          i,d,j,nmat=pep->nmat,lda=nmat*k,deg=nmat-1,nloc;
  PetscReal         *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat;
  PetscBLASInt      k_,lda_,lds_,nloc_,one=1,info;
  Mat               *A=ctx->A,Mk,M1=ctx->M1,P;
  BV                V=ctx->V,M2=ctx->M2,M3=ctx->M3,W=ctx->W;
  MatStructure      str;
  Vec               vc;

  PetscFunctionBegin;
  ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
  ierr = PetscMalloc3(nmat*k*k,&T12,k*k,&Tr,PetscMax(k*k,nmat),&Ts);CHKERRQ(ierr);
  DHii = T12;
  ierr = PetscArrayzero(DHii,k*k*nmat);CHKERRQ(ierr);
  for (i=0;i<k;i++) DHii[k+i+i*lda] = 1.0/a[0];
  for (d=2;d<nmat;d++) {
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) {
        DHii[d*k+i+j*lda] = ((h-b[d-1])*DHii[(d-1)*k+i+j*lda]+fH[(d-1)*k+i+j*lda]-g[d-1]*DHii[(d-2)*k+i+j*lda])/(a[d-1]);
      }
    }
  }
  /* T11 */
  if (!ctx->compM1) {
    ierr = MatCopy(A[0],M1,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = PEPEvaluateBasis(pep,h,0,Ts,NULL);CHKERRQ(ierr);
    for (j=1;j<nmat;j++) {
      ierr = MatAXPY(M1,Ts[j],A[j],str);CHKERRQ(ierr);
    }
  }

  /* T22 */
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lda,&lda_);CHKERRQ(ierr);
  PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,S,&lds_,S,&lds_,&zero,Tr,&k_));
  for (i=1;i<deg;i++) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,Tr,&k_,DHii+i*k,&lda_,&zero,Ts,&k_));
    s = (i==1)?0.0:1.0;
    PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,fH+i*k,&lda_,Ts,&k_,&s,M4,&k_));
  }
  for (i=0;i<k;i++) for (j=0;j<i;j++) { t=M4[i+j*k];M4[i+j*k]=M4[j+i*k];M4[j+i*k]=t; }

  /* T12 */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&Mk);CHKERRQ(ierr);
  for (i=1;i<nmat;i++) {
    ierr = MatDenseGetArray(Mk,&array);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,S,&lds_,DHii+i*k,&lda_,&zero,array,&k_));
    ierr = MatDenseRestoreArray(Mk,&array);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(W,0,k);CHKERRQ(ierr);
    ierr = BVMult(W,1.0,0.0,V,Mk);CHKERRQ(ierr);
    if (i==1) {
      ierr = BVMatMult(W,A[i],M2);CHKERRQ(ierr);
    } else {
      ierr = BVMatMult(W,A[i],M3);CHKERRQ(ierr); /* using M3 as work space */
      ierr = BVMult(M2,1.0,1.0,M3,NULL);CHKERRQ(ierr);
    }
  }

  /* T21 */
  ierr = MatDenseGetArray(Mk,&array);CHKERRQ(ierr);
  for (i=1;i<deg;i++) {
    s = (i==1)?0.0:1.0;
    ss = PetscConj(fh[i]);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&ss,S,&lds_,fH+i*k,&lda_,&s,array,&k_));
  }
  ierr = MatDenseRestoreArray(Mk,&array);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(M3,0,k);CHKERRQ(ierr);
  ierr = BVMult(M3,1.0,0.0,V,Mk);CHKERRQ(ierr);
  for (i=0;i<k;i++) {
    ierr = BVGetColumn(M3,i,&vc);CHKERRQ(ierr);
    ierr = VecConjugate(vc);CHKERRQ(ierr);
    ierr = BVRestoreColumn(M3,i,&vc);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&Mk);CHKERRQ(ierr);
  ierr = PetscFree3(T12,Tr,Ts);CHKERRQ(ierr);

  ierr = VecGetLocalSize(ctx->t,&nloc);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(nloc,&nloc_);CHKERRQ(ierr);
  ierr = PetscMalloc1(nloc*k,&T);CHKERRQ(ierr);
  ierr = KSPGetOperators(pep->refineksp,NULL,&P);CHKERRQ(ierr);
  if (!ctx->compM1) { ierr = MatCopy(ctx->M1,P,SAME_NONZERO_PATTERN);CHKERRQ(ierr); }
  ierr = BVGetArrayRead(ctx->M2,&m2);CHKERRQ(ierr);
  ierr = BVGetArrayRead(ctx->M3,&m3);CHKERRQ(ierr);
  ierr = VecGetArray(ctx->t,&v);CHKERRQ(ierr);
  for (i=0;i<nloc;i++) for (j=0;j<k;j++) T[j+i*k] = m3[i+j*nloc];
  PetscStackCallBLAS("LAPACKgesv",LAPACKgesv_(&k_,&nloc_,ctx->M4,&k_,ctx->pM4,T,&k_,&info));
  SlepcCheckLapackInfo("gesv",info);
  for (i=0;i<nloc;i++) v[i] = BLASdot_(&k_,m2+i,&nloc_,T+i*k,&one);
  ierr = VecRestoreArray(ctx->t,&v);CHKERRQ(ierr);
  ierr = BVRestoreArrayRead(ctx->M2,&m2);CHKERRQ(ierr);
  ierr = BVRestoreArrayRead(ctx->M3,&m3);CHKERRQ(ierr);
  ierr = MatDiagonalSet(P,ctx->t,ADD_VALUES);CHKERRQ(ierr);
  ierr = PetscFree(T);CHKERRQ(ierr);
  ierr = KSPSetUp(pep->refineksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefSysSolve_shell(KSP ksp,PetscInt nmat,Vec Rv,PetscScalar *Rh,PetscInt k,Vec dVi,PetscScalar *dHi)
{
  PetscErrorCode      ierr;
  PetscScalar         *t0;
  PetscBLASInt        k_,one=1,info,lda_;
  PetscInt            i,lda=nmat*k;
  Mat                 M;
  PEP_REFINE_MATSHELL *ctx;

  PetscFunctionBegin;
  ierr = KSPGetOperators(ksp,&M,NULL);CHKERRQ(ierr);
  ierr = MatShellGetContext(M,&ctx);CHKERRQ(ierr);
  ierr = PetscCalloc1(k,&t0);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lda,&lda_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  for (i=0;i<k;i++) t0[i] = Rh[i];
  PetscStackCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&k_,&one,ctx->M4,&k_,ctx->pM4,t0,&k_,&info));
  SlepcCheckLapackInfo("getrs",info);
  ierr = BVMultVec(ctx->M2,-1.0,1.0,Rv,t0);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,Rv,dVi);CHKERRQ(ierr);
  ierr = VecConjugate(dVi);CHKERRQ(ierr);
  ierr = BVDotVec(ctx->M3,dVi,dHi);CHKERRQ(ierr);
  ierr = VecConjugate(dVi);CHKERRQ(ierr);
  for (i=0;i<k;i++) dHi[i] = Rh[i]-PetscConj(dHi[i]);
  PetscStackCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&k_,&one,ctx->M4,&k_,ctx->pM4,dHi,&k_,&info));
  SlepcCheckLapackInfo("getrs",info);
  ierr = PetscFree(t0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Computes the residual P(H,V*S)*e_j for the polynomial
*/
static PetscErrorCode NRefRightSide(PetscInt nmat,PetscReal *pcf,Mat *A,PetscInt k,BV V,PetscScalar *S,PetscInt lds,PetscInt j,PetscScalar *H,PetscInt ldh,PetscScalar *fH,PetscScalar *DfH,PetscScalar *dH,BV dV,PetscScalar *dVS,PetscInt rds,Vec Rv,PetscScalar *Rh,BV W,Vec t)
{
  PetscErrorCode ierr;
  PetscScalar    *DS0,*DS1,*F,beta=0.0,sone=1.0,none=-1.0,tt=0.0,*h,zero=0.0,*Z,*c0;
  PetscReal      *a=pcf,*b=pcf+nmat,*g=b+nmat;
  PetscInt       i,ii,jj,lda;
  PetscBLASInt   lda_,k_,ldh_,lds_,nmat_,k2_,krds_,j_,one=1;
  Mat            M0;
  Vec            w;

  PetscFunctionBegin;
  ierr = PetscMalloc4(k*nmat,&h,k*k,&DS0,k*k,&DS1,k*k,&Z);CHKERRQ(ierr);
  lda = k*nmat;
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lda,&lda_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(nmat,&nmat_);CHKERRQ(ierr);
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&nmat_,&k_,&sone,S,&lds_,fH+j*lda,&k_,&zero,h,&k_));
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,nmat,h,&M0);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(W,0,nmat);CHKERRQ(ierr);
  ierr = BVMult(W,1.0,0.0,V,M0);CHKERRQ(ierr);
  ierr = MatDestroy(&M0);CHKERRQ(ierr);

  ierr = BVGetColumn(W,0,&w);CHKERRQ(ierr);
  ierr = MatMult(A[0],w,Rv);CHKERRQ(ierr);
  ierr = BVRestoreColumn(W,0,&w);CHKERRQ(ierr);
  for (i=1;i<nmat;i++) {
    ierr = BVGetColumn(W,i,&w);CHKERRQ(ierr);
    ierr = MatMult(A[i],w,t);CHKERRQ(ierr);
    ierr = BVRestoreColumn(W,i,&w);CHKERRQ(ierr);
    ierr = VecAXPY(Rv,1.0,t);CHKERRQ(ierr);
  }
  /* Update right-hand side */
  if (j) {
    ierr = PetscBLASIntCast(ldh,&ldh_);CHKERRQ(ierr);
    ierr = PetscArrayzero(Z,k*k);CHKERRQ(ierr);
    ierr = PetscArrayzero(DS0,k*k);CHKERRQ(ierr);
    ierr = PetscArraycpy(Z+(j-1)*k,dH+(j-1)*k,k);CHKERRQ(ierr);
    /* Update DfH */
    for (i=1;i<nmat;i++) {
      if (i>1) {
        beta = -g[i-1];
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,fH+(i-1)*k,&lda_,Z,&k_,&beta,DS0,&k_));
        tt += -b[i-1];
        for (ii=0;ii<k;ii++) H[ii+ii*ldh] += tt;
        tt = b[i-1];
        beta = 1.0/a[i-1];
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&beta,DS1,&k_,H,&ldh_,&beta,DS0,&k_));
        F = DS0; DS0 = DS1; DS1 = F;
      } else {
        ierr = PetscArrayzero(DS1,k*k);CHKERRQ(ierr);
        for (ii=0;ii<k;ii++) DS1[ii+(j-1)*k] = Z[ii+(j-1)*k]/a[0];
      }
      for (jj=j;jj<k;jj++) {
        for (ii=0;ii<k;ii++) DfH[k*i+ii+jj*lda] += DS1[ii+jj*k];
      }
    }
    for (ii=0;ii<k;ii++) H[ii+ii*ldh] += tt;
    /* Update right-hand side */
    ierr = PetscBLASIntCast(2*k,&k2_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(j,&j_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(k+rds,&krds_);CHKERRQ(ierr);
    c0 = DS0;
    ierr = PetscArrayzero(Rh,k);CHKERRQ(ierr);
    for (i=0;i<nmat;i++) {
      PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&krds_,&j_,&sone,dVS,&k2_,fH+j*lda+i*k,&one,&zero,h,&one));
      PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&k_,&k_,&sone,S,&lds_,DfH+i*k+j*lda,&one,&sone,h,&one));
      ierr = BVMultVec(V,1.0,0.0,t,h);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(dV,0,rds);CHKERRQ(ierr);
      ierr = BVMultVec(dV,1.0,1.0,t,h+k);CHKERRQ(ierr);
      ierr = BVGetColumn(W,i,&w);CHKERRQ(ierr);
      ierr = MatMult(A[i],t,w);CHKERRQ(ierr);
      ierr = BVRestoreColumn(W,i,&w);CHKERRQ(ierr);
      if (i>0 && i<nmat-1) {
        PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&k_,&k_,&sone,S,&lds_,h,&one,&zero,c0,&one));
        PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&k_,&k_,&none,fH+i*k,&lda_,c0,&one,&sone,Rh,&one));
      }
    }

    for (i=0;i<nmat;i++) h[i] = -1.0;
    ierr = BVMultVec(W,1.0,1.0,Rv,h);CHKERRQ(ierr);
  }
  ierr = PetscFree4(h,DS0,DS1,Z);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefSysSolve_mbe(PetscInt k,PetscInt sz,BV W,PetscScalar *w,BV Wt,PetscScalar *wt,PetscScalar *d,PetscScalar *dt,KSP ksp,BV T2,BV T3 ,PetscScalar *T4,PetscBool trans,Vec x1,PetscScalar *x2,Vec sol1,PetscScalar *sol2,Vec vw)
{
  PetscErrorCode ierr;
  PetscInt       i,j,incf,incc;
  PetscScalar    *y,*g,*xx2,*ww,y2,*dd;
  Vec            v,t,xx1;
  BV             WW,T;

  PetscFunctionBegin;
  ierr = PetscMalloc3(sz,&y,sz,&g,k,&xx2);CHKERRQ(ierr);
  if (trans) {
    WW = W; ww = w; dd = d; T = T3; incf = 0; incc = 1;
  } else {
    WW = Wt; ww = wt; dd = dt; T = T2; incf = 1; incc = 0;
  }
  xx1 = vw;
  ierr = VecCopy(x1,xx1);CHKERRQ(ierr);
  ierr = PetscArraycpy(xx2,x2,sz);CHKERRQ(ierr);
  ierr = PetscArrayzero(sol2,k);CHKERRQ(ierr);
  for (i=sz-1;i>=0;i--) {
    ierr = BVGetColumn(WW,i,&v);CHKERRQ(ierr);
    ierr = VecConjugate(v);CHKERRQ(ierr);
    ierr = VecDot(xx1,v,y+i);CHKERRQ(ierr);
    ierr = VecConjugate(v);CHKERRQ(ierr);
    ierr = BVRestoreColumn(WW,i,&v);CHKERRQ(ierr);
    for (j=0;j<i;j++) y[i] += ww[j+i*k]*xx2[j];
    y[i] = -(y[i]-xx2[i])/dd[i];
    ierr = BVGetColumn(T,i,&t);CHKERRQ(ierr);
    ierr = VecAXPY(xx1,-y[i],t);CHKERRQ(ierr);
    ierr = BVRestoreColumn(T,i,&t);CHKERRQ(ierr);
    for (j=0;j<=i;j++) xx2[j] -= y[i]*T4[j*incf+incc*i+(i*incf+incc*j)*k];
    g[i] = xx2[i];
  }
  if (trans) {
    ierr = KSPSolveTranspose(ksp,xx1,sol1);CHKERRQ(ierr);
  } else {
    ierr = KSPSolve(ksp,xx1,sol1);CHKERRQ(ierr);
  }
  if (trans) {
    WW = Wt; ww = wt; dd = dt; T = T2; incf = 1; incc = 0;
  } else {
    WW = W; ww = w; dd = d; T = T3; incf = 0; incc = 1;
  }
  for (i=0;i<sz;i++) {
    ierr = BVGetColumn(T,i,&t);CHKERRQ(ierr);
    ierr = VecConjugate(t);CHKERRQ(ierr);
    ierr = VecDot(sol1,t,&y2);CHKERRQ(ierr);
    ierr = VecConjugate(t);CHKERRQ(ierr);
    ierr = BVRestoreColumn(T,i,&t);CHKERRQ(ierr);
    for (j=0;j<i;j++) y2 += sol2[j]*T4[j*incf+incc*i+(i*incf+incc*j)*k];
    y2 = (g[i]-y2)/dd[i];
    ierr = BVGetColumn(WW,i,&v);CHKERRQ(ierr);
    ierr = VecAXPY(sol1,-y2,v);CHKERRQ(ierr);
    for (j=0;j<i;j++) sol2[j] -= ww[j+i*k]*y2;
    sol2[i] = y[i]+y2;
    ierr = BVRestoreColumn(WW,i,&v);CHKERRQ(ierr);
  }
  ierr = PetscFree3(y,g,xx2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefSysSetup_mbe(PEP pep,PetscInt k,KSP ksp,PetscScalar *fH,PetscScalar *S,PetscInt lds,PetscScalar *fh,PetscScalar h,BV V,PEP_REFINE_EXPLICIT *matctx)
{
  PetscErrorCode ierr;
  PetscInt       i,j,l,nmat=pep->nmat,lda=nmat*k,deg=nmat-1;
  Mat            M1=matctx->M1,*A,*At,Mk;
  PetscReal      *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat;
  PetscScalar    s,ss,*DHii,*T12,*array,*Ts,*Tr,*M4=matctx->M4,sone=1.0,zero=0.0;
  PetscScalar    *w=matctx->w,*wt=matctx->wt,*d=matctx->d,*dt=matctx->dt;
  PetscBLASInt   lds_,lda_,k_;
  MatStructure   str;
  PetscBool      flg;
  BV             M2=matctx->M2,M3=matctx->M3,W=matctx->W,Wt=matctx->Wt;
  Vec            vc,vc2;

  PetscFunctionBegin;
  ierr = PetscMalloc3(nmat*k*k,&T12,k*k,&Tr,PetscMax(k*k,nmat),&Ts);CHKERRQ(ierr);
  ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscMalloc1(pep->nmat,&At);CHKERRQ(ierr);
    for (i=0;i<pep->nmat;i++) {
      ierr = STGetMatrixTransformed(pep->st,i,&At[i]);CHKERRQ(ierr);
    }
  } else At = pep->A;
  if (matctx->subc) A = matctx->A;
  else A = At;
  /* Form the explicit system matrix */
  DHii = T12;
  ierr = PetscArrayzero(DHii,k*k*nmat);CHKERRQ(ierr);
  for (i=0;i<k;i++) DHii[k+i+i*lda] = 1.0/a[0];
  for (l=2;l<nmat;l++) {
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) {
        DHii[l*k+i+j*lda] = ((h-b[l-1])*DHii[(l-1)*k+i+j*lda]+fH[(l-1)*k+i+j*lda]-g[l-1]*DHii[(l-2)*k+i+j*lda])/a[l-1];
      }
    }
  }

  /* T11 */
  if (!matctx->compM1) {
    ierr = MatCopy(A[0],M1,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = PEPEvaluateBasis(pep,h,0,Ts,NULL);CHKERRQ(ierr);
    for (j=1;j<nmat;j++) {
      ierr = MatAXPY(M1,Ts[j],A[j],str);CHKERRQ(ierr);
    }
  }

  /* T22 */
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lda,&lda_);CHKERRQ(ierr);
  PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,S,&lds_,S,&lds_,&zero,Tr,&k_));
  for (i=1;i<deg;i++) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,Tr,&k_,DHii+i*k,&lda_,&zero,Ts,&k_));
    s = (i==1)?0.0:1.0;
    PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,fH+i*k,&lda_,Ts,&k_,&s,M4,&k_));
  }

  /* T12 */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&Mk);CHKERRQ(ierr);
  for (i=1;i<nmat;i++) {
    ierr = MatDenseGetArray(Mk,&array);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,S,&lds_,DHii+i*k,&lda_,&zero,array,&k_));
    ierr = MatDenseRestoreArray(Mk,&array);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(W,0,k);CHKERRQ(ierr);
    ierr = BVMult(W,1.0,0.0,V,Mk);CHKERRQ(ierr);
    if (i==1) {
      ierr = BVMatMult(W,A[i],M2);CHKERRQ(ierr);
    } else {
      ierr = BVMatMult(W,A[i],M3);CHKERRQ(ierr); /* using M3 as work space */
      ierr = BVMult(M2,1.0,1.0,M3,NULL);CHKERRQ(ierr);
    }
  }

  /* T21 */
  ierr = MatDenseGetArray(Mk,&array);CHKERRQ(ierr);
  for (i=1;i<deg;i++) {
    s = (i==1)?0.0:1.0;
    ss = PetscConj(fh[i]);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&ss,S,&lds_,fH+i*k,&lda_,&s,array,&k_));
  }
  ierr = MatDenseRestoreArray(Mk,&array);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(M3,0,k);CHKERRQ(ierr);
  ierr = BVMult(M3,1.0,0.0,V,Mk);CHKERRQ(ierr);
  for (i=0;i<k;i++) {
    ierr = BVGetColumn(M3,i,&vc);CHKERRQ(ierr);
    ierr = VecConjugate(vc);CHKERRQ(ierr);
    ierr = BVRestoreColumn(M3,i,&vc);CHKERRQ(ierr);
  }

  ierr = KSPSetOperators(ksp,M1,M1);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&Mk);CHKERRQ(ierr);

  /* Set up for BEMW */
  for (i=0;i<k;i++) {
    ierr = BVGetColumn(M2,i,&vc);CHKERRQ(ierr);
    ierr = BVGetColumn(W,i,&vc2);CHKERRQ(ierr);
    ierr = NRefSysSolve_mbe(k,i,W,w,Wt,wt,d,dt,ksp,M2,M3,M4,PETSC_FALSE,vc,M4+i*k,vc2,w+i*k,matctx->t);CHKERRQ(ierr);
    ierr = BVRestoreColumn(M2,i,&vc);CHKERRQ(ierr);
    ierr = BVGetColumn(M3,i,&vc);CHKERRQ(ierr);
    ierr = VecConjugate(vc);CHKERRQ(ierr);
    ierr = VecDot(vc2,vc,&d[i]);CHKERRQ(ierr);
    ierr = VecConjugate(vc);CHKERRQ(ierr);
    ierr = BVRestoreColumn(M3,i,&vc);CHKERRQ(ierr);
    for (j=0;j<i;j++) d[i] += M4[i+j*k]*w[j+i*k];
    d[i] = M4[i+i*k]-d[i];
    ierr = BVRestoreColumn(W,i,&vc2);CHKERRQ(ierr);

    ierr = BVGetColumn(M3,i,&vc);CHKERRQ(ierr);
    ierr = BVGetColumn(Wt,i,&vc2);CHKERRQ(ierr);
    for (j=0;j<=i;j++) Ts[j] = M4[i+j*k];
    ierr = NRefSysSolve_mbe(k,i,W,w,Wt,wt,d,dt,ksp,M2,M3,M4,PETSC_TRUE,vc,Ts,vc2,wt+i*k,matctx->t);CHKERRQ(ierr);
    ierr = BVRestoreColumn(M3,i,&vc);CHKERRQ(ierr);
    ierr = BVGetColumn(M2,i,&vc);CHKERRQ(ierr);
    ierr = VecConjugate(vc2);CHKERRQ(ierr);
    ierr = VecDot(vc,vc2,&dt[i]);CHKERRQ(ierr);
    ierr = VecConjugate(vc2);CHKERRQ(ierr);
    ierr = BVRestoreColumn(M2,i,&vc);CHKERRQ(ierr);
    for (j=0;j<i;j++) dt[i] += M4[j+i*k]*wt[j+i*k];
    dt[i] = M4[i+i*k]-dt[i];
    ierr = BVRestoreColumn(Wt,i,&vc2);CHKERRQ(ierr);
  }

  if (flg) {
    ierr = PetscFree(At);CHKERRQ(ierr);
  }
  ierr = PetscFree3(T12,Tr,Ts);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefSysSetup_explicit(PEP pep,PetscInt k,KSP ksp,PetscScalar *fH,PetscScalar *S,PetscInt lds,PetscScalar *fh,PetscScalar h,BV V,PEP_REFINE_EXPLICIT *matctx,BV W)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,d,n,n0,m0,n1,m1,nmat=pep->nmat,lda=nmat*k,deg=nmat-1;
  PetscInt          *idxg=matctx->idxg,*idxp=matctx->idxp,idx,ncols;
  Mat               M,*E=matctx->E,*A,*At,Mk,Md;
  PetscReal         *a=pep->pbc,*b=pep->pbc+nmat,*g=pep->pbc+2*nmat;
  PetscScalar       s,ss,*DHii,*T22,*T21,*T12,*Ts,*Tr,*array,*ts,sone=1.0,zero=0.0;
  PetscBLASInt      lds_,lda_,k_;
  const PetscInt    *idxmc;
  const PetscScalar *valsc,*carray;
  MatStructure      str;
  Vec               vc,vc0;
  PetscBool         flg;

  PetscFunctionBegin;
  ierr = PetscMalloc5(k*k,&T22,k*k,&T21,nmat*k*k,&T12,k*k,&Tr,k*k,&Ts);CHKERRQ(ierr);
  ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp,&M,NULL);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(E[1],&n1,&m1);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(E[0],&n0,&m0);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(M,&n,NULL);CHKERRQ(ierr);
  ierr = PetscMalloc1(nmat,&ts);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscMalloc1(pep->nmat,&At);CHKERRQ(ierr);
    for (i=0;i<pep->nmat;i++) {
      ierr = STGetMatrixTransformed(pep->st,i,&At[i]);CHKERRQ(ierr);
    }
  } else At = pep->A;
  if (matctx->subc) A = matctx->A;
  else A = At;
  /* Form the explicit system matrix */
  DHii = T12;
  ierr = PetscArrayzero(DHii,k*k*nmat);CHKERRQ(ierr);
  for (i=0;i<k;i++) DHii[k+i+i*lda] = 1.0/a[0];
  for (d=2;d<nmat;d++) {
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) {
        DHii[d*k+i+j*lda] = ((h-b[d-1])*DHii[(d-1)*k+i+j*lda]+fH[(d-1)*k+i+j*lda]-g[d-1]*DHii[(d-2)*k+i+j*lda])/a[d-1];
      }
    }
  }

  /* T11 */
  if (!matctx->compM1) {
    ierr = MatCopy(A[0],E[0],DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = PEPEvaluateBasis(pep,h,0,Ts,NULL);CHKERRQ(ierr);
    for (j=1;j<nmat;j++) {
      ierr = MatAXPY(E[0],Ts[j],A[j],str);CHKERRQ(ierr);
    }
  }
  for (i=n0;i<m0;i++) {
    ierr = MatGetRow(E[0],i,&ncols,&idxmc,&valsc);CHKERRQ(ierr);
    idx = n+i-n0;
    for (j=0;j<ncols;j++) {
      idxg[j] = matctx->map0[idxmc[j]];
    }
    ierr = MatSetValues(M,1,&idx,ncols,idxg,valsc,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatRestoreRow(E[0],i,&ncols,&idxmc,&valsc);CHKERRQ(ierr);
  }

  /* T22 */
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lda,&lda_);CHKERRQ(ierr);
  PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,S,&lds_,S,&lds_,&zero,Tr,&k_));
  for (i=1;i<deg;i++) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,Tr,&k_,DHii+i*k,&lda_,&zero,Ts,&k_));
    s = (i==1)?0.0:1.0;
    PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&k_,&k_,&k_,&sone,fH+i*k,&lda_,Ts,&k_,&s,T22,&k_));
  }
  for (j=0;j<k;j++) idxp[j] = matctx->map1[j];
  for (i=0;i<m1-n1;i++) {
    idx = n+m0-n0+i;
    for (j=0;j<k;j++) {
      Tr[j] = T22[n1+i+j*k];
    }
    ierr = MatSetValues(M,1,&idx,k,idxp,Tr,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* T21 */
  for (i=1;i<deg;i++) {
    s = (i==1)?0.0:1.0;
    ss = PetscConj(fh[i]);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&ss,S,&lds_,fH+i*k,&lda_,&s,T21,&k_));
  }
  ierr = BVSetActiveColumns(W,0,k);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,T21,&Mk);CHKERRQ(ierr);
  ierr = BVMult(W,1.0,0.0,V,Mk);CHKERRQ(ierr);
  for (i=0;i<k;i++) {
    ierr = BVGetColumn(W,i,&vc);CHKERRQ(ierr);
    ierr = VecConjugate(vc);CHKERRQ(ierr);
    ierr = VecGetArrayRead(vc,&carray);CHKERRQ(ierr);
    idx = matctx->map1[i];
    ierr = MatSetValues(M,1,&idx,m0-n0,matctx->map0+n0,carray,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(vc,&carray);CHKERRQ(ierr);
    ierr = BVRestoreColumn(W,i,&vc);CHKERRQ(ierr);
  }

  /* T12 */
  for (i=1;i<nmat;i++) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,S,&lds_,DHii+i*k,&lda_,&zero,Ts,&k_));
    for (j=0;j<k;j++) {
      ierr = PetscArraycpy(T12+i*k+j*lda,Ts+j*k,k);CHKERRQ(ierr);
    }
  }
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,nmat-1,NULL,&Md);CHKERRQ(ierr);
  for (i=0;i<nmat;i++) ts[i] = 1.0;
  for (j=0;j<k;j++) {
    ierr = MatDenseGetArray(Md,&array);CHKERRQ(ierr);
    ierr = PetscArraycpy(array,T12+k+j*lda,(nmat-1)*k);CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(Md,&array);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(W,0,nmat-1);CHKERRQ(ierr);
    ierr = BVMult(W,1.0,0.0,V,Md);CHKERRQ(ierr);
    for (i=nmat-1;i>0;i--) {
      ierr = BVGetColumn(W,i-1,&vc0);CHKERRQ(ierr);
      ierr = BVGetColumn(W,i,&vc);CHKERRQ(ierr);
      ierr = MatMult(A[i],vc0,vc);CHKERRQ(ierr);
      ierr = BVRestoreColumn(W,i-1,&vc0);CHKERRQ(ierr);
      ierr = BVRestoreColumn(W,i,&vc);CHKERRQ(ierr);
    }
    ierr = BVSetActiveColumns(W,1,nmat);CHKERRQ(ierr);
    ierr = BVGetColumn(W,0,&vc0);CHKERRQ(ierr);
    ierr = BVMultVec(W,1.0,0.0,vc0,ts);CHKERRQ(ierr);
    ierr = VecGetArrayRead(vc0,&carray);CHKERRQ(ierr);
    idx = matctx->map1[j];
    ierr = MatSetValues(M,m0-n0,matctx->map0+n0,1,&idx,carray,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(vc0,&carray);CHKERRQ(ierr);
    ierr = BVRestoreColumn(W,0,&vc0);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,M,M);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  ierr = PetscFree(ts);CHKERRQ(ierr);
  ierr = MatDestroy(&Mk);CHKERRQ(ierr);
  ierr = MatDestroy(&Md);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscFree(At);CHKERRQ(ierr);
  }
  ierr = PetscFree5(T22,T21,T12,Tr,Ts);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefSysSolve_explicit(PetscInt k,KSP ksp,Vec Rv,PetscScalar *Rh,Vec dVi,PetscScalar *dHi,PEP_REFINE_EXPLICIT *matctx)
{
  PetscErrorCode    ierr;
  PetscInt          n0,m0,n1,m1,i;
  PetscScalar       *arrayV;
  const PetscScalar *array;

  PetscFunctionBegin;
  ierr = MatGetOwnershipRange(matctx->E[1],&n1,&m1);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(matctx->E[0],&n0,&m0);CHKERRQ(ierr);

  /* Right side */
  ierr = VecGetArrayRead(Rv,&array);CHKERRQ(ierr);
  ierr = VecSetValues(matctx->tN,m0-n0,matctx->map0+n0,array,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Rv,&array);CHKERRQ(ierr);
  ierr = VecSetValues(matctx->tN,m1-n1,matctx->map1+n1,Rh+n1,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(matctx->tN);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(matctx->tN);CHKERRQ(ierr);

  /* Solve */
  ierr = KSPSolve(ksp,matctx->tN,matctx->ttN);CHKERRQ(ierr);

  /* Retrieve solution */
  ierr = VecGetArray(dVi,&arrayV);CHKERRQ(ierr);
  ierr = VecGetArrayRead(matctx->ttN,&array);CHKERRQ(ierr);
  ierr = PetscArraycpy(arrayV,array,m0-n0);CHKERRQ(ierr);
  ierr = VecRestoreArray(dVi,&arrayV);CHKERRQ(ierr);
  if (!matctx->subc) {
    ierr = VecGetArray(matctx->t1,&arrayV);CHKERRQ(ierr);
    for (i=0;i<m1-n1;i++) arrayV[i] =  array[m0-n0+i];
    ierr = VecRestoreArray(matctx->t1,&arrayV);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(matctx->ttN,&array);CHKERRQ(ierr);
    ierr = VecScatterBegin(matctx->scatterctx,matctx->t1,matctx->vseq,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(matctx->scatterctx,matctx->t1,matctx->vseq,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecGetArrayRead(matctx->vseq,&array);CHKERRQ(ierr);
    for (i=0;i<k;i++) dHi[i] = array[i];
    ierr = VecRestoreArrayRead(matctx->vseq,&array);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefSysIter(PetscInt i,PEP pep,PetscInt k,KSP ksp,PetscScalar *fH,PetscScalar *S,PetscInt lds,PetscScalar *fh,PetscScalar *H,PetscInt ldh,Vec Rv,PetscScalar *Rh,BV V,Vec dVi,PetscScalar *dHi,PEP_REFINE_EXPLICIT *matctx,BV W)
{
  PetscErrorCode      ierr;
  PetscInt            j,m,lda=pep->nmat*k,n0,m0,idx;
  PetscMPIInt         root,len;
  PetscScalar         *array2,h;
  const PetscScalar   *array;
  Vec                 R,Vi;
  PEP_REFINE_MATSHELL *ctx;
  Mat                 M;

  PetscFunctionBegin;
  if (!matctx || !matctx->subc) {
    for (j=0;j<pep->nmat;j++) fh[j] = fH[j*k+i+i*lda];
    h   = H[i+i*ldh];
    idx = i;
    R   = Rv;
    Vi  = dVi;
    switch (pep->scheme) {
    case PEP_REFINE_SCHEME_EXPLICIT:
      ierr = NRefSysSetup_explicit(pep,k,ksp,fH,S,lds,fh,h,V,matctx,W);CHKERRQ(ierr);
      matctx->compM1 = PETSC_FALSE;
      break;
    case PEP_REFINE_SCHEME_MBE:
      ierr = NRefSysSetup_mbe(pep,k,ksp,fH,S,lds,fh,h,V,matctx);CHKERRQ(ierr);
      matctx->compM1 = PETSC_FALSE;
      break;
    case PEP_REFINE_SCHEME_SCHUR:
      ierr = KSPGetOperators(ksp,&M,NULL);CHKERRQ(ierr);
      ierr = MatShellGetContext(M,&ctx);CHKERRQ(ierr);
      ierr = NRefSysSetup_shell(pep,k,fH,S,lds,fh,h,ctx);CHKERRQ(ierr);
      ctx->compM1 = PETSC_FALSE;
      break;
    }
  } else {
    if (i%matctx->subc->n==0 && (idx=i+matctx->subc->color)<k) {
      for (j=0;j<pep->nmat;j++) fh[j] = fH[j*k+idx+idx*lda];
      h = H[idx+idx*ldh];
      matctx->idx = idx;
      switch (pep->scheme) {
      case PEP_REFINE_SCHEME_EXPLICIT:
        ierr = NRefSysSetup_explicit(pep,k,ksp,fH,S,lds,fh,h,matctx->V,matctx,matctx->W);CHKERRQ(ierr);
        matctx->compM1 = PETSC_FALSE;
        break;
      case PEP_REFINE_SCHEME_MBE:
        ierr = NRefSysSetup_mbe(pep,k,ksp,fH,S,lds,fh,h,matctx->V,matctx);CHKERRQ(ierr);
        matctx->compM1 = PETSC_FALSE;
        break;
      case PEP_REFINE_SCHEME_SCHUR:
        break;
      }
    } else idx = matctx->idx;
    ierr = VecScatterBegin(matctx->scatter_id[i%matctx->subc->n],Rv,matctx->tg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(matctx->scatter_id[i%matctx->subc->n],Rv,matctx->tg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecGetArrayRead(matctx->tg,&array);CHKERRQ(ierr);
    ierr = VecPlaceArray(matctx->t,array);CHKERRQ(ierr);
    ierr = VecCopy(matctx->t,matctx->Rv);CHKERRQ(ierr);
    ierr = VecResetArray(matctx->t);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(matctx->tg,&array);CHKERRQ(ierr);
    R  = matctx->Rv;
    Vi = matctx->Vi;
  }
  if (idx==i && idx<k) {
    switch (pep->scheme) {
      case PEP_REFINE_SCHEME_EXPLICIT:
        ierr = NRefSysSolve_explicit(k,ksp,R,Rh,Vi,dHi,matctx);CHKERRQ(ierr);
        break;
      case PEP_REFINE_SCHEME_MBE:
        ierr = NRefSysSolve_mbe(k,k,matctx->W,matctx->w,matctx->Wt,matctx->wt,matctx->d,matctx->dt,ksp,matctx->M2,matctx->M3 ,matctx->M4,PETSC_FALSE,R,Rh,Vi,dHi,matctx->t);CHKERRQ(ierr);
        break;
      case PEP_REFINE_SCHEME_SCHUR:
        ierr = NRefSysSolve_shell(ksp,pep->nmat,R,Rh,k,Vi,dHi);CHKERRQ(ierr);
        break;
    }
  }
  if (matctx && matctx->subc) {
    ierr = VecGetLocalSize(Vi,&m);CHKERRQ(ierr);
    ierr = VecGetArrayRead(Vi,&array);CHKERRQ(ierr);
    ierr = VecGetArray(matctx->tg,&array2);CHKERRQ(ierr);
    ierr = PetscArraycpy(array2,array,m);CHKERRQ(ierr);
    ierr = VecRestoreArray(matctx->tg,&array2);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(Vi,&array);CHKERRQ(ierr);
    ierr = VecScatterBegin(matctx->scatter_id[i%matctx->subc->n],matctx->tg,dVi,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
    ierr = VecScatterEnd(matctx->scatter_id[i%matctx->subc->n],matctx->tg,dVi,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
    switch (pep->scheme) {
    case PEP_REFINE_SCHEME_EXPLICIT:
      ierr = MatGetOwnershipRange(matctx->E[0],&n0,&m0);CHKERRQ(ierr);
      ierr = VecGetArrayRead(matctx->ttN,&array);CHKERRQ(ierr);
      ierr = VecPlaceArray(matctx->tp,array+m0-n0);CHKERRQ(ierr);
      ierr = VecScatterBegin(matctx->scatterp_id[i%matctx->subc->n],matctx->tp,matctx->tpg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = VecScatterEnd(matctx->scatterp_id[i%matctx->subc->n],matctx->tp,matctx->tpg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = VecResetArray(matctx->tp);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(matctx->ttN,&array);CHKERRQ(ierr);
      ierr = VecGetArrayRead(matctx->tpg,&array);CHKERRQ(ierr);
      for (j=0;j<k;j++) dHi[j] = array[j];
      ierr = VecRestoreArrayRead(matctx->tpg,&array);CHKERRQ(ierr);
      break;
     case PEP_REFINE_SCHEME_MBE:
      root = 0;
      for (j=0;j<i%matctx->subc->n;j++) root += matctx->subc->subsize[j];
      ierr = PetscMPIIntCast(k,&len);CHKERRQ(ierr);
      ierr = MPI_Bcast(dHi,len,MPIU_SCALAR,root,matctx->subc->dupparent);CHKERRQ(ierr);
      break;
    case PEP_REFINE_SCHEME_SCHUR:
      break;
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPNRefForwardSubstitution(PEP pep,PetscInt k,PetscScalar *S,PetscInt lds,PetscScalar *H,PetscInt ldh,PetscScalar *fH,BV dV,PetscScalar *dVS,PetscInt *rds,PetscScalar *dH,PetscInt lddh,KSP ksp,PEP_REFINE_EXPLICIT *matctx)
{
  PetscErrorCode      ierr;
  PetscInt            i,nmat=pep->nmat,lda=nmat*k;
  PetscScalar         *fh,*Rh,*DfH;
  PetscReal           norm;
  BV                  W;
  Vec                 Rv,t,dvi;
  PEP_REFINE_MATSHELL *ctx;
  Mat                 M,*At;
  PetscBool           flg,lindep;

  PetscFunctionBegin;
  ierr = PetscMalloc2(nmat*k*k,&DfH,k,&Rh);CHKERRQ(ierr);
  *rds = 0;
  ierr = BVCreateVec(pep->V,&Rv);CHKERRQ(ierr);
  switch (pep->scheme) {
  case PEP_REFINE_SCHEME_EXPLICIT:
    ierr = BVCreateVec(pep->V,&t);CHKERRQ(ierr);
    ierr = BVDuplicateResize(pep->V,PetscMax(k,nmat),&W);CHKERRQ(ierr);
    ierr = PetscMalloc1(nmat,&fh);CHKERRQ(ierr);
    break;
  case PEP_REFINE_SCHEME_MBE:
    if (matctx->subc) {
      ierr = BVCreateVec(pep->V,&t);CHKERRQ(ierr);
      ierr = BVDuplicateResize(pep->V,PetscMax(k,nmat),&W);CHKERRQ(ierr);
    } else {
      W = matctx->W;
      ierr = PetscObjectReference((PetscObject)W);CHKERRQ(ierr);
      t = matctx->t;
      ierr = PetscObjectReference((PetscObject)t);CHKERRQ(ierr);
    }
    ierr = BVScale(matctx->W,0.0);CHKERRQ(ierr);
    ierr = BVScale(matctx->Wt,0.0);CHKERRQ(ierr);
    ierr = BVScale(matctx->M2,0.0);CHKERRQ(ierr);
    ierr = BVScale(matctx->M3,0.0);CHKERRQ(ierr);
    ierr = PetscMalloc1(nmat,&fh);CHKERRQ(ierr);
    break;
  case PEP_REFINE_SCHEME_SCHUR:
    ierr = KSPGetOperators(ksp,&M,NULL);CHKERRQ(ierr);
    ierr = MatShellGetContext(M,&ctx);CHKERRQ(ierr);
    ierr = BVCreateVec(pep->V,&t);CHKERRQ(ierr);
    ierr = BVDuplicateResize(pep->V,PetscMax(k,nmat),&W);CHKERRQ(ierr);
    fh = ctx->fih;
    break;
  }
  ierr = PetscArrayzero(dVS,2*k*k);CHKERRQ(ierr);
  ierr = PetscArrayzero(DfH,lda*k);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscMalloc1(pep->nmat,&At);CHKERRQ(ierr);
    for (i=0;i<pep->nmat;i++) {
      ierr = STGetMatrixTransformed(pep->st,i,&At[i]);CHKERRQ(ierr);
    }
  } else At = pep->A;

  /* Main loop for computing the ith columns of dX and dS */
  for (i=0;i<k;i++) {
    /* Compute and update i-th column of the right hand side */
    ierr = PetscArrayzero(Rh,k);CHKERRQ(ierr);
    ierr = NRefRightSide(nmat,pep->pbc,At,k,pep->V,S,lds,i,H,ldh,fH,DfH,dH,dV,dVS,*rds,Rv,Rh,W,t);CHKERRQ(ierr);

    /* Update and solve system */
    ierr = BVGetColumn(dV,i,&dvi);CHKERRQ(ierr);
    ierr = NRefSysIter(i,pep,k,ksp,fH,S,lds,fh,H,ldh,Rv,Rh,pep->V,dvi,dH+i*k,matctx,W);CHKERRQ(ierr);
    /* Orthogonalize computed solution */
    ierr = BVOrthogonalizeVec(pep->V,dvi,dVS+i*2*k,&norm,&lindep);CHKERRQ(ierr);
    ierr = BVRestoreColumn(dV,i,&dvi);CHKERRQ(ierr);
    if (!lindep) {
      ierr = BVOrthogonalizeColumn(dV,i,dVS+k+i*2*k,&norm,&lindep);CHKERRQ(ierr);
      if (!lindep) {
        dVS[k+i+i*2*k] = norm;
        ierr = BVScaleColumn(dV,i,1.0/norm);CHKERRQ(ierr);
        (*rds)++;
      }
    }
  }
  ierr = BVSetActiveColumns(dV,0,*rds);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = VecDestroy(&Rv);CHKERRQ(ierr);
  ierr = BVDestroy(&W);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscFree(At);CHKERRQ(ierr);
  }
  ierr = PetscFree2(DfH,Rh);CHKERRQ(ierr);
  if (pep->scheme!=PEP_REFINE_SCHEME_SCHUR) { ierr = PetscFree(fh);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefOrthogStep(PEP pep,PetscInt k,PetscScalar *H,PetscInt ldh,PetscScalar *fH,PetscScalar *S,PetscInt lds)
{
  PetscErrorCode ierr;
  PetscInt       j,nmat=pep->nmat,deg=nmat-1,lda=nmat*k,ldg;
  PetscScalar    *G,*tau,sone=1.0,zero=0.0,*work;
  PetscBLASInt   lds_,k_,ldh_,info,ldg_,lda_;

  PetscFunctionBegin;
  ierr = PetscMalloc3(k,&tau,k,&work,deg*k*k,&G);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lda,&lda_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);

  /* Form auxiliary matrix for the orthogonalization step */
  ldg = deg*k;
  ierr = PEPEvaluateBasisforMatrix(pep,nmat,k,H,ldh,fH);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldg,&ldg_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldh,&ldh_);CHKERRQ(ierr);
  for (j=0;j<deg;j++) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&k_,&k_,&k_,&sone,S,&lds_,fH+j*k,&lda_,&zero,G+j*k,&ldg_));
  }
  /* Orthogonalize and update S */
  PetscStackCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&ldg_,&k_,G,&ldg_,tau,work,&k_,&info));
  SlepcCheckLapackInfo("geqrf",info);
  PetscStackCallBLAS("BLAStrsm",BLAStrsm_("R","U","N","N",&k_,&k_,&sone,G,&ldg_,S,&lds_));

  /* Update H */
  PetscStackCallBLAS("BLAStrmm",BLAStrmm_("L","U","N","N",&k_,&k_,&sone,G,&ldg_,H,&ldh_));
  PetscStackCallBLAS("BLAStrsm",BLAStrsm_("R","U","N","N",&k_,&k_,&sone,G,&ldg_,H,&ldh_));
  ierr = PetscFree3(tau,work,G);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPNRefUpdateInvPair(PEP pep,PetscInt k,PetscScalar *H,PetscInt ldh,PetscScalar *fH,PetscScalar *dH,PetscScalar *S,PetscInt lds,BV dV,PetscScalar *dVS,PetscInt rds)
{
  PetscErrorCode ierr;
  PetscInt       i,j,nmat=pep->nmat,lda=nmat*k;
  PetscScalar    *tau,*array,*work;
  PetscBLASInt   lds_,k_,lda_,ldh_,kdrs_,info,k2_;
  Mat            M0;

  PetscFunctionBegin;
  ierr = PetscMalloc2(k,&tau,k,&work);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(lda,&lda_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldh,&ldh_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(2*k,&k2_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast((k+rds),&kdrs_);CHKERRQ(ierr);
  /* Update H */
  for (j=0;j<k;j++) {
    for (i=0;i<k;i++) H[i+j*ldh] -= dH[i+j*k];
  }
  /* Update V */
  for (j=0;j<k;j++) {
    for (i=0;i<k;i++) dVS[i+j*2*k] = -dVS[i+j*2*k]+S[i+j*lds];
    for (i=k;i<2*k;i++) dVS[i+j*2*k] = -dVS[i+j*2*k];
  }
  PetscStackCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&kdrs_,&k_,dVS,&k2_,tau,work,&k_,&info));
  SlepcCheckLapackInfo("geqrf",info);
  /* Copy triangular matrix in S */
  for (j=0;j<k;j++) {
    for (i=0;i<=j;i++) S[i+j*lds] = dVS[i+j*2*k];
    for (i=j+1;i<k;i++) S[i+j*lds] = 0.0;
  }
  PetscStackCallBLAS("LAPACKorgqr",LAPACKorgqr_(&k2_,&k_,&k_,dVS,&k2_,tau,work,&k_,&info));
  SlepcCheckLapackInfo("orgqr",info);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M0);CHKERRQ(ierr);
  ierr = MatDenseGetArray(M0,&array);CHKERRQ(ierr);
  for (j=0;j<k;j++) {
    ierr = PetscArraycpy(array+j*k,dVS+j*2*k,k);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(M0,&array);CHKERRQ(ierr);
  ierr = BVMultInPlace(pep->V,M0,0,k);CHKERRQ(ierr);
  if (rds) {
    ierr = MatDenseGetArray(M0,&array);CHKERRQ(ierr);
    for (j=0;j<k;j++) {
      ierr = PetscArraycpy(array+j*k,dVS+k+j*2*k,rds);CHKERRQ(ierr);
    }
    ierr = MatDenseRestoreArray(M0,&array);CHKERRQ(ierr);
    ierr = BVMultInPlace(dV,M0,0,k);CHKERRQ(ierr);
    ierr = BVMult(pep->V,1.0,1.0,dV,NULL);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&M0);CHKERRQ(ierr);
  ierr = NRefOrthogStep(pep,k,H,ldh,fH,S,lds);CHKERRQ(ierr);
  ierr = PetscFree2(tau,work);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPNRefSetUp(PEP pep,PetscInt k,PetscScalar *H,PetscInt ldh,PEP_REFINE_EXPLICIT *matctx,PetscBool ini)
{
  PetscErrorCode      ierr;
  PEP_REFINE_MATSHELL *ctx;
  PetscScalar         t,*coef;
  const PetscScalar   *array;
  MatStructure        str;
  PetscInt            j,nmat=pep->nmat,n0,m0,n1,m1,n0_,m0_,n1_,m1_,N0,N1,p,*idx1,*idx2,count,si,i,l0;
  MPI_Comm            comm;
  PetscMPIInt         np;
  const PetscInt      *rgs0,*rgs1;
  Mat                 B,C,*E,*A,*At;
  IS                  is1,is2;
  Vec                 v;
  PetscBool           flg;
  Mat                 M,P;

  PetscFunctionBegin;
  ierr = PetscMalloc1(nmat,&coef);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscMalloc1(pep->nmat,&At);CHKERRQ(ierr);
    for (i=0;i<pep->nmat;i++) {
      ierr = STGetMatrixTransformed(pep->st,i,&At[i]);CHKERRQ(ierr);
    }
  } else At = pep->A;
  switch (pep->scheme) {
  case PEP_REFINE_SCHEME_EXPLICIT:
    if (ini) {
      if (matctx->subc) {
        A = matctx->A;
        comm = PetscSubcommChild(matctx->subc);
      } else {
        A = At;
        ierr = PetscObjectGetComm((PetscObject)pep,&comm);CHKERRQ(ierr);
      }
      E = matctx->E;
      ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
      ierr = MatDuplicate(A[0],MAT_COPY_VALUES,&E[0]);CHKERRQ(ierr);
      j = (matctx->subc)?matctx->subc->color:0;
      ierr = PEPEvaluateBasis(pep,H[j+j*ldh],0,coef,NULL);CHKERRQ(ierr);
      for (j=1;j<nmat;j++) {
        ierr = MatAXPY(E[0],coef[j],A[j],str);CHKERRQ(ierr);
      }
      ierr = MatCreateDense(comm,PETSC_DECIDE,PETSC_DECIDE,k,k,NULL,&E[1]);CHKERRQ(ierr);
      ierr = MatAssemblyBegin(E[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(E[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatGetOwnershipRange(E[0],&n0,&m0);CHKERRQ(ierr);
      ierr = MatGetOwnershipRange(E[1],&n1,&m1);CHKERRQ(ierr);
      ierr = MatGetOwnershipRangeColumn(E[0],&n0_,&m0_);CHKERRQ(ierr);
      ierr = MatGetOwnershipRangeColumn(E[1],&n1_,&m1_);CHKERRQ(ierr);
      /* T12 and T21 are computed from V and V*, so,
         they must have the same column and row ranges */
      if (m0_-n0_ != m0-n0) SETERRQ(PETSC_COMM_SELF,1,"Inconsistent dimensions");
      ierr = MatCreateDense(comm,m0-n0,m1_-n1_,PETSC_DECIDE,PETSC_DECIDE,NULL,&B);CHKERRQ(ierr);
      ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatCreateDense(comm,m1-n1,m0_-n0_,PETSC_DECIDE,PETSC_DECIDE,NULL,&C);CHKERRQ(ierr);
      ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatCreateTile(1.0,E[0],1.0,B,1.0,C,1.0,E[1],&M);CHKERRQ(ierr);
      ierr = MatDestroy(&B);CHKERRQ(ierr);
      ierr = MatDestroy(&C);CHKERRQ(ierr);
      matctx->compM1 = PETSC_TRUE;
      ierr = MatGetSize(E[0],NULL,&N0);CHKERRQ(ierr);
      ierr = MatGetSize(E[1],NULL,&N1);CHKERRQ(ierr);
      ierr = MPI_Comm_size(PetscObjectComm((PetscObject)M),&np);CHKERRQ(ierr);
      ierr = MatGetOwnershipRanges(E[0],&rgs0);CHKERRQ(ierr);
      ierr = MatGetOwnershipRanges(E[1],&rgs1);CHKERRQ(ierr);
      ierr = PetscMalloc4(PetscMax(k,N1),&matctx->idxp,N0,&matctx->idxg,N0,&matctx->map0,N1,&matctx->map1);CHKERRQ(ierr);
      /* Create column (and row) mapping */
      for (p=0;p<np;p++) {
        for (j=rgs0[p];j<rgs0[p+1];j++) matctx->map0[j] = j+rgs1[p];
        for (j=rgs1[p];j<rgs1[p+1];j++) matctx->map1[j] = j+rgs0[p+1];
      }
      ierr = MatCreateVecs(M,NULL,&matctx->tN);CHKERRQ(ierr);
      ierr = MatCreateVecs(matctx->E[1],NULL,&matctx->t1);CHKERRQ(ierr);
      ierr = VecDuplicate(matctx->tN,&matctx->ttN);CHKERRQ(ierr);
      if (matctx->subc) {
        ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np);CHKERRQ(ierr);
        count = np*k;
        ierr = PetscMalloc2(count,&idx1,count,&idx2);CHKERRQ(ierr);
        ierr = VecCreateMPI(PetscObjectComm((PetscObject)pep),m1-n1,PETSC_DECIDE,&matctx->tp);CHKERRQ(ierr);
        ierr = VecGetOwnershipRange(matctx->tp,&l0,NULL);CHKERRQ(ierr);
        ierr = VecCreateMPI(PetscObjectComm((PetscObject)pep),k,PETSC_DECIDE,&matctx->tpg);CHKERRQ(ierr);
        for (si=0;si<matctx->subc->n;si++) {
          if (matctx->subc->color==si) {
            j=0;
            if (matctx->subc->color==si) {
              for (p=0;p<np;p++) {
                for (i=n1;i<m1;i++) {
                  idx1[j] = l0+i-n1;
                  idx2[j++] =p*k+i;
                }
              }
            }
            count = np*(m1-n1);
          } else count =0;
          ierr = ISCreateGeneral(PetscObjectComm((PetscObject)pep),count,idx1,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
          ierr = ISCreateGeneral(PetscObjectComm((PetscObject)pep),count,idx2,PETSC_COPY_VALUES,&is2);CHKERRQ(ierr);
          ierr = VecScatterCreate(matctx->tp,is1,matctx->tpg,is2,&matctx->scatterp_id[si]);CHKERRQ(ierr);
          ierr = ISDestroy(&is1);CHKERRQ(ierr);
          ierr = ISDestroy(&is2);CHKERRQ(ierr);
        }
        ierr = PetscFree2(idx1,idx2);CHKERRQ(ierr);
      } else {
        ierr = VecScatterCreateToAll(matctx->t1,&matctx->scatterctx,&matctx->vseq);CHKERRQ(ierr);
      }
      P = M;
    } else {
      if (matctx->subc) {
        /* Scatter vectors pep->V */
        for (i=0;i<k;i++) {
          ierr = BVGetColumn(pep->V,i,&v);CHKERRQ(ierr);
          ierr = VecScatterBegin(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
          ierr = VecScatterEnd(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
          ierr = BVRestoreColumn(pep->V,i,&v);CHKERRQ(ierr);
          ierr = VecGetArrayRead(matctx->tg,&array);CHKERRQ(ierr);
          ierr = VecPlaceArray(matctx->t,(const PetscScalar*)array);CHKERRQ(ierr);
          ierr = BVInsertVec(matctx->V,i,matctx->t);CHKERRQ(ierr);
          ierr = VecResetArray(matctx->t);CHKERRQ(ierr);
          ierr = VecRestoreArrayRead(matctx->tg,&array);CHKERRQ(ierr);
        }
      }
    }
    break;
  case PEP_REFINE_SCHEME_MBE:
    if (ini) {
      if (matctx->subc) {
        A = matctx->A;
        comm = PetscSubcommChild(matctx->subc);
      } else {
        matctx->V = pep->V;
        A = At;
        ierr = PetscObjectGetComm((PetscObject)pep,&comm);CHKERRQ(ierr);
        ierr = MatCreateVecs(pep->A[0],&matctx->t,NULL);CHKERRQ(ierr);
      }
      ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
      ierr = MatDuplicate(A[0],MAT_COPY_VALUES,&matctx->M1);CHKERRQ(ierr);
      j = (matctx->subc)?matctx->subc->color:0;
      ierr = PEPEvaluateBasis(pep,H[j+j*ldh],0,coef,NULL);CHKERRQ(ierr);
      for (j=1;j<nmat;j++) {
        ierr = MatAXPY(matctx->M1,coef[j],A[j],str);CHKERRQ(ierr);
      }
      ierr = BVDuplicateResize(matctx->V,PetscMax(k,pep->nmat),&matctx->W);CHKERRQ(ierr);
      ierr = BVDuplicateResize(matctx->V,k,&matctx->M2);CHKERRQ(ierr);
      ierr = BVDuplicate(matctx->M2,&matctx->M3);CHKERRQ(ierr);
      ierr = BVDuplicate(matctx->M2,&matctx->Wt);CHKERRQ(ierr);
      ierr = PetscMalloc5(k*k,&matctx->M4,k*k,&matctx->w,k*k,&matctx->wt,k,&matctx->d,k,&matctx->dt);CHKERRQ(ierr);
      matctx->compM1 = PETSC_TRUE;
      M = matctx->M1;
      P = M;
    }
    break;
  case PEP_REFINE_SCHEME_SCHUR:
    if (ini) {
      ierr = PetscObjectGetComm((PetscObject)pep,&comm);CHKERRQ(ierr);
      ierr = MatGetSize(At[0],&m0,&n0);CHKERRQ(ierr);
      ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);
      ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
      /* Create a shell matrix to solve the linear system */
      ctx->V = pep->V;
      ctx->k = k; ctx->nmat = nmat;
      ierr = PetscMalloc5(nmat,&ctx->A,k*k,&ctx->M4,k,&ctx->pM4,2*k*k,&ctx->work,nmat,&ctx->fih);CHKERRQ(ierr);
      for (i=0;i<nmat;i++) ctx->A[i] = At[i];
      ierr = PetscArrayzero(ctx->M4,k*k);CHKERRQ(ierr);
      ierr = MatCreateShell(comm,PETSC_DECIDE,PETSC_DECIDE,m0,n0,ctx,&M);CHKERRQ(ierr);
      ierr = MatShellSetOperation(M,MATOP_MULT,(void(*)(void))MatMult_FS);CHKERRQ(ierr);
      ierr = BVDuplicateResize(ctx->V,PetscMax(k,pep->nmat),&ctx->W);CHKERRQ(ierr);
      ierr = BVDuplicateResize(ctx->V,k,&ctx->M2);CHKERRQ(ierr);
      ierr = BVDuplicate(ctx->M2,&ctx->M3);CHKERRQ(ierr);
      ierr = BVCreateVec(pep->V,&ctx->t);CHKERRQ(ierr);
      ierr = MatDuplicate(At[0],MAT_COPY_VALUES,&ctx->M1);CHKERRQ(ierr);
      ierr = PEPEvaluateBasis(pep,H[0],0,coef,NULL);CHKERRQ(ierr);
      for (j=1;j<nmat;j++) {
        ierr = MatAXPY(ctx->M1,coef[j],At[j],str);CHKERRQ(ierr);
      }
      ierr = MatDuplicate(At[0],MAT_COPY_VALUES,&P);CHKERRQ(ierr);
      /* Compute a precond matrix for the system */
      t = H[0];
      ierr = PEPEvaluateBasis(pep,t,0,coef,NULL);CHKERRQ(ierr);
      for (j=1;j<nmat;j++) {
        ierr = MatAXPY(P,coef[j],At[j],str);CHKERRQ(ierr);
      }
      ctx->compM1 = PETSC_TRUE;
    }
    break;
  }
  if (ini) {
    ierr = PEPRefineGetKSP(pep,&pep->refineksp);CHKERRQ(ierr);
    ierr = KSPSetErrorIfNotConverged(pep->refineksp,PETSC_TRUE);CHKERRQ(ierr);
    ierr = KSPSetOperators(pep->refineksp,M,P);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(pep->refineksp);CHKERRQ(ierr);
  }

  if (!ini && matctx && matctx->subc) {
     /* Scatter vectors pep->V */
    for (i=0;i<k;i++) {
      ierr = BVGetColumn(pep->V,i,&v);CHKERRQ(ierr);
      ierr = VecScatterBegin(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = VecScatterEnd(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pep->V,i,&v);CHKERRQ(ierr);
      ierr = VecGetArrayRead(matctx->tg,&array);CHKERRQ(ierr);
      ierr = VecPlaceArray(matctx->t,(const PetscScalar*)array);CHKERRQ(ierr);
      ierr = BVInsertVec(matctx->V,i,matctx->t);CHKERRQ(ierr);
      ierr = VecResetArray(matctx->t);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(matctx->tg,&array);CHKERRQ(ierr);
    }
   }
  ierr = PetscFree(coef);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscFree(At);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefSubcommSetup(PEP pep,PetscInt k,PEP_REFINE_EXPLICIT *matctx,PetscInt nsubc)
{
  PetscErrorCode    ierr;
  PetscInt          i,si,j,m0,n0,nloc0,nloc_sub,*idx1,*idx2;
  IS                is1,is2;
  BVType            type;
  Vec               v;
  const PetscScalar *array;
  Mat               *A;
  PetscBool         flg;

  PetscFunctionBegin;
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscMalloc1(pep->nmat,&A);CHKERRQ(ierr);
    for (i=0;i<pep->nmat;i++) {
      ierr = STGetMatrixTransformed(pep->st,i,&A[i]);CHKERRQ(ierr);
    }
  } else A = pep->A;

  /* Duplicate pep matrices */
  ierr = PetscMalloc3(pep->nmat,&matctx->A,nsubc,&matctx->scatter_id,nsubc,&matctx->scatterp_id);CHKERRQ(ierr);
  for (i=0;i<pep->nmat;i++) {
    ierr = MatCreateRedundantMatrix(A[i],0,PetscSubcommChild(matctx->subc),MAT_INITIAL_MATRIX,&matctx->A[i]);CHKERRQ(ierr);
  }

  /* Create Scatter */
  ierr = MatCreateVecs(matctx->A[0],&matctx->t,NULL);CHKERRQ(ierr);
  ierr = MatGetLocalSize(matctx->A[0],&nloc_sub,NULL);CHKERRQ(ierr);
  ierr = VecCreateMPI(PetscSubcommContiguousParent(matctx->subc),nloc_sub,PETSC_DECIDE,&matctx->tg);CHKERRQ(ierr);
  ierr = BVGetColumn(pep->V,0,&v);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(v,&n0,&m0);CHKERRQ(ierr);
  nloc0 = m0-n0;
  ierr = PetscMalloc2(matctx->subc->n*nloc0,&idx1,matctx->subc->n*nloc0,&idx2);CHKERRQ(ierr);
  j = 0;
  for (si=0;si<matctx->subc->n;si++) {
    for (i=n0;i<m0;i++) {
      idx1[j]   = i;
      idx2[j++] = i+pep->n*si;
    }
  }
  ierr = ISCreateGeneral(PetscObjectComm((PetscObject)pep),matctx->subc->n*nloc0,idx1,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
  ierr = ISCreateGeneral(PetscObjectComm((PetscObject)pep),matctx->subc->n*nloc0,idx2,PETSC_COPY_VALUES,&is2);CHKERRQ(ierr);
  ierr = VecScatterCreate(v,is1,matctx->tg,is2,&matctx->scatter_sub);CHKERRQ(ierr);
  ierr = ISDestroy(&is1);CHKERRQ(ierr);
  ierr = ISDestroy(&is2);CHKERRQ(ierr);
  for (si=0;si<matctx->subc->n;si++) {
    j=0;
    for (i=n0;i<m0;i++) {
      idx1[j] = i;
      idx2[j++] = i+pep->n*si;
    }
    ierr = ISCreateGeneral(PetscObjectComm((PetscObject)pep),nloc0,idx1,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
    ierr = ISCreateGeneral(PetscObjectComm((PetscObject)pep),nloc0,idx2,PETSC_COPY_VALUES,&is2);CHKERRQ(ierr);
    ierr = VecScatterCreate(v,is1,matctx->tg,is2,&matctx->scatter_id[si]);CHKERRQ(ierr);
    ierr = ISDestroy(&is1);CHKERRQ(ierr);
    ierr = ISDestroy(&is2);CHKERRQ(ierr);
  }
  ierr = BVRestoreColumn(pep->V,0,&v);CHKERRQ(ierr);
  ierr = PetscFree2(idx1,idx2);CHKERRQ(ierr);

  /* Duplicate pep->V vecs */
  ierr = BVGetType(pep->V,&type);CHKERRQ(ierr);
  ierr = BVCreate(PetscSubcommChild(matctx->subc),&matctx->V);CHKERRQ(ierr);
  ierr = BVSetType(matctx->V,type);CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(matctx->V,matctx->t,k);CHKERRQ(ierr);
  if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT) {
    ierr = BVDuplicateResize(matctx->V,PetscMax(k,pep->nmat),&matctx->W);CHKERRQ(ierr);
  }
  for (i=0;i<k;i++) {
    ierr = BVGetColumn(pep->V,i,&v);CHKERRQ(ierr);
    ierr = VecScatterBegin(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(matctx->scatter_sub,v,matctx->tg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = BVRestoreColumn(pep->V,i,&v);CHKERRQ(ierr);
    ierr = VecGetArrayRead(matctx->tg,&array);CHKERRQ(ierr);
    ierr = VecPlaceArray(matctx->t,(const PetscScalar*)array);CHKERRQ(ierr);
    ierr = BVInsertVec(matctx->V,i,matctx->t);CHKERRQ(ierr);
    ierr = VecResetArray(matctx->t);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(matctx->tg,&array);CHKERRQ(ierr);
  }

  ierr = VecDuplicate(matctx->t,&matctx->Rv);CHKERRQ(ierr);
  ierr = VecDuplicate(matctx->t,&matctx->Vi);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscFree(A);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NRefSubcommDestroy(PEP pep,PEP_REFINE_EXPLICIT *matctx)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = VecScatterDestroy(&matctx->scatter_sub);CHKERRQ(ierr);
  for (i=0;i<matctx->subc->n;i++) {
    ierr = VecScatterDestroy(&matctx->scatter_id[i]);CHKERRQ(ierr);
  }
  for (i=0;i<pep->nmat;i++) {
    ierr = MatDestroy(&matctx->A[i]);CHKERRQ(ierr);
  }
  if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT) {
    for (i=0;i<matctx->subc->n;i++) {
      ierr = VecScatterDestroy(&matctx->scatterp_id[i]);CHKERRQ(ierr);
    }
    ierr = VecDestroy(&matctx->tp);CHKERRQ(ierr);
    ierr = VecDestroy(&matctx->tpg);CHKERRQ(ierr);
    ierr = BVDestroy(&matctx->W);CHKERRQ(ierr);
  }
  ierr = PetscFree3(matctx->A,matctx->scatter_id,matctx->scatterp_id);CHKERRQ(ierr);
  ierr = BVDestroy(&matctx->V);CHKERRQ(ierr);
  ierr = VecDestroy(&matctx->t);CHKERRQ(ierr);
  ierr = VecDestroy(&matctx->tg);CHKERRQ(ierr);
  ierr = VecDestroy(&matctx->Rv);CHKERRQ(ierr);
  ierr = VecDestroy(&matctx->Vi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPNewtonRefinement_TOAR(PEP pep,PetscScalar sigma,PetscInt *maxits,PetscReal *tol,PetscInt k,PetscScalar *S,PetscInt lds)
{
  PetscErrorCode      ierr;
  PetscScalar         *H,*work,*dH,*fH,*dVS;
  PetscInt            ldh,i,j,its=1,nmat=pep->nmat,nsubc=pep->npart,rds;
  PetscBLASInt        k_,ld_,*p,info;
  BV                  dV;
  PetscBool           sinvert,flg;
  PEP_REFINE_EXPLICIT *matctx=NULL;
  Vec                 v;
  Mat                 M,P;
  PEP_REFINE_MATSHELL *ctx;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(PEP_Refine,pep,0,0,0);CHKERRQ(ierr);
  if (k > pep->n) SETERRQ1(PetscObjectComm((PetscObject)pep),1,"Multiple Refinement available only for invariant pairs of dimension smaller than n=%D",pep->n);
  /* the input tolerance is not being taken into account (by the moment) */
  its = *maxits;
  ierr = PetscMalloc3(k*k,&dH,nmat*k*k,&fH,k,&work);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(pep->ds,&ldh);CHKERRQ(ierr);
  ierr = DSGetArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
  ierr = DSRestoreArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
  ierr = PetscMalloc1(2*k*k,&dVS);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (!flg && pep->st && pep->ops->backtransform) { /* BackTransform */
    ierr = PetscBLASIntCast(k,&k_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(ldh,&ld_);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinvert);CHKERRQ(ierr);
    if (sinvert) {
      ierr = DSGetArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
      ierr = PetscMalloc1(k,&p);CHKERRQ(ierr);
      PetscStackCallBLAS("LAPACKgetrf",LAPACKgetrf_(&k_,&k_,H,&ld_,p,&info));
      SlepcCheckLapackInfo("getrf",info);
      PetscStackCallBLAS("LAPACKgetri",LAPACKgetri_(&k_,H,&ld_,p,work,&k_,&info));
      SlepcCheckLapackInfo("getri",info);
      ierr = DSRestoreArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
      pep->ops->backtransform = NULL;
    }
    if (sigma!=0.0) {
      ierr = DSGetArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
      for (i=0;i<k;i++) H[i+ldh*i] += sigma;
      ierr = DSRestoreArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
      pep->ops->backtransform = NULL;
    }
  }
  if ((pep->scale==PEP_SCALE_BOTH || pep->scale==PEP_SCALE_SCALAR) && pep->sfactor!=1.0) {
    ierr = DSGetArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
    for (j=0;j<k;j++) {
      for (i=0;i<k;i++) H[i+j*ldh] *= pep->sfactor;
    }
    ierr = DSRestoreArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
    if (!flg) {
      /* Restore original values */
      for (i=0;i<pep->nmat;i++){
        pep->pbc[pep->nmat+i] *= pep->sfactor;
        pep->pbc[2*pep->nmat+i] *= pep->sfactor*pep->sfactor;
      }
    }
  }
  if ((pep->scale==PEP_SCALE_DIAGONAL || pep->scale==PEP_SCALE_BOTH) && pep->Dr) {
    for (i=0;i<k;i++) {
      ierr = BVGetColumn(pep->V,i,&v);CHKERRQ(ierr);
      ierr = VecPointwiseMult(v,v,pep->Dr);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pep->V,i,&v);CHKERRQ(ierr);
    }
  }
  ierr = DSGetArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);

  ierr = NRefOrthogStep(pep,k,H,ldh,fH,S,lds);CHKERRQ(ierr);
  /* check if H is in Schur form */
  for (i=0;i<k-1;i++) {
    if (H[i+1+i*ldh]!=0.0) {
#if !defined(PETSC_USE_COMPLEX)
      SETERRQ(PetscObjectComm((PetscObject)pep),1,"Iterative Refinement requires the complex Schur form of the projected matrix");
#else
      SETERRQ(PetscObjectComm((PetscObject)pep),1,"Iterative Refinement requires an upper triangular projected matrix");
#endif
    }
  }
  if (nsubc>k) SETERRQ(PetscObjectComm((PetscObject)pep),1,"Amount of subcommunicators should not be larger than the invariant pair dimension");
  ierr = BVSetActiveColumns(pep->V,0,k);CHKERRQ(ierr);
  ierr = BVDuplicateResize(pep->V,k,&dV);CHKERRQ(ierr);
  ierr = PetscLogObjectParent((PetscObject)pep,(PetscObject)dV);CHKERRQ(ierr);
  if (pep->scheme!=PEP_REFINE_SCHEME_SCHUR) {
    ierr = PetscMalloc1(1,&matctx);CHKERRQ(ierr);
    if (nsubc>1) { /* spliting in subcommunicators */
      matctx->subc = pep->refinesubc;
      ierr = NRefSubcommSetup(pep,k,matctx,nsubc);CHKERRQ(ierr);
    } else matctx->subc=NULL;
  }

  /* Loop performing iterative refinements */
  for (i=0;i<its;i++) {
    /* Pre-compute the polynomial basis evaluated in H */
    ierr = PEPEvaluateBasisforMatrix(pep,nmat,k,H,ldh,fH);CHKERRQ(ierr);
    ierr = PEPNRefSetUp(pep,k,H,ldh,matctx,PetscNot(i));CHKERRQ(ierr);
    /* Solve the linear system */
    ierr = PEPNRefForwardSubstitution(pep,k,S,lds,H,ldh,fH,dV,dVS,&rds,dH,k,pep->refineksp,matctx);CHKERRQ(ierr);
    /* Update X (=V*S) and H, and orthogonalize [X;X*fH1;...;XfH(deg-1)] */
    ierr = PEPNRefUpdateInvPair(pep,k,H,ldh,fH,dH,S,lds,dV,dVS,rds);CHKERRQ(ierr);
  }
  ierr = DSRestoreArray(pep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
  if (!flg && sinvert) {
    ierr = PetscFree(p);CHKERRQ(ierr);
  }
  ierr = PetscFree3(dH,fH,work);CHKERRQ(ierr);
  ierr = PetscFree(dVS);CHKERRQ(ierr);
  ierr = BVDestroy(&dV);CHKERRQ(ierr);
  switch (pep->scheme) {
  case PEP_REFINE_SCHEME_EXPLICIT:
    for (i=0;i<2;i++) {
      ierr = MatDestroy(&matctx->E[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree4(matctx->idxp,matctx->idxg,matctx->map0,matctx->map1);CHKERRQ(ierr);
    ierr = VecDestroy(&matctx->tN);CHKERRQ(ierr);
    ierr = VecDestroy(&matctx->ttN);CHKERRQ(ierr);
    ierr = VecDestroy(&matctx->t1);CHKERRQ(ierr);
    if (nsubc>1) {
      ierr = NRefSubcommDestroy(pep,matctx);CHKERRQ(ierr);
    } else {
      ierr = VecDestroy(&matctx->vseq);CHKERRQ(ierr);
      ierr = VecScatterDestroy(&matctx->scatterctx);CHKERRQ(ierr);
    }
    ierr = PetscFree(matctx);CHKERRQ(ierr);
    ierr = KSPGetOperators(pep->refineksp,&M,NULL);CHKERRQ(ierr);
    ierr = MatDestroy(&M);CHKERRQ(ierr);
    break;
  case PEP_REFINE_SCHEME_MBE:
    ierr = BVDestroy(&matctx->W);CHKERRQ(ierr);
    ierr = BVDestroy(&matctx->Wt);CHKERRQ(ierr);
    ierr = BVDestroy(&matctx->M2);CHKERRQ(ierr);
    ierr = BVDestroy(&matctx->M3);CHKERRQ(ierr);
    ierr = MatDestroy(&matctx->M1);CHKERRQ(ierr);
    ierr = VecDestroy(&matctx->t);CHKERRQ(ierr);
    ierr = PetscFree5(matctx->M4,matctx->w,matctx->wt,matctx->d,matctx->dt);CHKERRQ(ierr);
    if (nsubc>1) {
      ierr = NRefSubcommDestroy(pep,matctx);CHKERRQ(ierr);
    }
    ierr = PetscFree(matctx);CHKERRQ(ierr);
    break;
  case PEP_REFINE_SCHEME_SCHUR:
    ierr = KSPGetOperators(pep->refineksp,&M,&P);CHKERRQ(ierr);
    ierr = MatShellGetContext(M,&ctx);CHKERRQ(ierr);
    ierr = PetscFree5(ctx->A,ctx->M4,ctx->pM4,ctx->work,ctx->fih);CHKERRQ(ierr);
    ierr = MatDestroy(&ctx->M1);CHKERRQ(ierr);
    ierr = BVDestroy(&ctx->M2);CHKERRQ(ierr);
    ierr = BVDestroy(&ctx->M3);CHKERRQ(ierr);
    ierr = BVDestroy(&ctx->W);CHKERRQ(ierr);
    ierr = VecDestroy(&ctx->t);CHKERRQ(ierr);
    ierr = PetscFree(ctx);CHKERRQ(ierr);
    ierr = MatDestroy(&M);CHKERRQ(ierr);
    ierr = MatDestroy(&P);CHKERRQ(ierr);
    break;
  }
  ierr = PetscLogEventEnd(PEP_Refine,pep,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

