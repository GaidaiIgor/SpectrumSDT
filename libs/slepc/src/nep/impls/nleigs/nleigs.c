/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "nleigs"

   Method: NLEIGS

   Algorithm:

       Fully rational Krylov method for nonlinear eigenvalue problems.

   References:

       [1] S. Guttel et al., "NLEIGS: A class of robust fully rational Krylov
           method for nonlinear eigenvalue problems", SIAM J. Sci. Comput.
           36(6):A2842-A2864, 2014.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <slepcblaslapack.h>
#include "nleigs.h"

PetscErrorCode NEPNLEIGSBackTransform(PetscObject ob,PetscInt n,PetscScalar *valr,PetscScalar *vali)
{
  NEP         nep;
  PetscInt    j;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar t;
#endif

  PetscFunctionBegin;
  nep = (NEP)ob;
#if !defined(PETSC_USE_COMPLEX)
  for (j=0;j<n;j++) {
    if (vali[j] == 0) valr[j] = 1.0 / valr[j] + nep->target;
    else {
      t = valr[j] * valr[j] + vali[j] * vali[j];
      valr[j] = valr[j] / t + nep->target;
      vali[j] = - vali[j] / t;
    }
  }
#else
  for (j=0;j<n;j++) {
    valr[j] = 1.0 / valr[j] + nep->target;
  }
#endif
  PetscFunctionReturn(0);
}

/* Computes the roots of a polynomial */
static PetscErrorCode NEPNLEIGSAuxiliarPRootFinder(PetscInt deg,PetscScalar *polcoeffs,PetscScalar *wr,PetscScalar *wi,PetscBool *avail)
{
  PetscErrorCode ierr;
  PetscScalar    *C;
  PetscBLASInt   n_,lwork;
  PetscInt       i;
#if defined(PETSC_USE_COMPLEX) || defined(PETSC_HAVE_ESSL)
  PetscReal      *rwork=NULL;
#endif
#if defined(PETSC_HAVE_ESSL)
  PetscScalar    sdummy,*wri;
  PetscBLASInt   idummy,io=0;
#else
  PetscScalar    *work;
  PetscBLASInt   info;
#endif

  PetscFunctionBegin;
  *avail = PETSC_TRUE;
  if (deg>0) {
    ierr = PetscCalloc1(deg*deg,&C);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(deg,&n_);CHKERRQ(ierr);
    for (i=0;i<deg-1;i++) {
      C[(deg+1)*i+1]   = 1.0;
      C[(deg-1)*deg+i] = -polcoeffs[deg-i]/polcoeffs[0];
    }
    C[deg*deg+-1] = -polcoeffs[1]/polcoeffs[0];
    ierr = PetscBLASIntCast(3*deg,&lwork);CHKERRQ(ierr);

    ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
#if !defined(PETSC_HAVE_ESSL)
#if !defined(PETSC_USE_COMPLEX)
    ierr = PetscMalloc1(lwork,&work);CHKERRQ(ierr);
    PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_("N","N",&n_,C,&n_,wr,wi,NULL,&n_,NULL,&n_,work,&lwork,&info));
    if (info) *avail = PETSC_FALSE;
    ierr = PetscFree(work);CHKERRQ(ierr);
#else
    ierr = PetscMalloc2(2*deg,&rwork,lwork,&work);CHKERRQ(ierr);
    PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_("N","N",&n_,C,&n_,wr,NULL,&n_,NULL,&n_,work,&lwork,rwork,&info));
    if (info) *avail = PETSC_FALSE;
    ierr = PetscFree2(rwork,work);CHKERRQ(ierr);
#endif
#else /* defined(PETSC_HAVE_ESSL) */
    ierr = PetscMalloc2(lwork,&rwork,2*deg,&wri);CHKERRQ(ierr);
    PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_(&io,C,&n_,wri,&sdummy,&idummy,&idummy,&n_,rwork,&lwork));
#if !defined(PETSC_USE_COMPLEX)
    for (i=0;i<deg;i++) {
      wr[i] = wri[2*i];
      wi[i] = wri[2*i+1];
    }
#else
    for (i=0;i<deg;i++) wr[i] = wri[i];
#endif
    ierr = PetscFree2(rwork,wri);CHKERRQ(ierr);
#endif
    ierr = PetscFPTrapPop();CHKERRQ(ierr);
    ierr = PetscFree(C);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSAuxiliarRmDuplicates(PetscInt nin,PetscScalar *pin,PetscInt *nout,PetscScalar *pout,PetscInt max)
{
  PetscInt i,j;

  PetscFunctionBegin;
  for (i=0;i<nin;i++) {
    if (max && *nout>=max) break;
    pout[(*nout)++] = pin[i];
    for (j=0;j<*nout-1;j++)
      if (PetscAbsScalar(pin[i]-pout[j])<PETSC_MACHINE_EPSILON*100) {
        (*nout)--;
        break;
      }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSFNSingularities(FN f,PetscInt *nisol,PetscScalar **isol,PetscBool *rational)
{
  PetscErrorCode ierr;
  FNCombineType  ctype;
  FN             f1,f2;
  PetscInt       i,nq,nisol1,nisol2;
  PetscScalar    *qcoeff,*wr,*wi,*isol1,*isol2;
  PetscBool      flg,avail,rat1,rat2;

  PetscFunctionBegin;
  *rational = PETSC_FALSE;
  ierr = PetscObjectTypeCompare((PetscObject)f,FNRATIONAL,&flg);CHKERRQ(ierr);
  if (flg) {
    *rational = PETSC_TRUE;
    ierr = FNRationalGetDenominator(f,&nq,&qcoeff);CHKERRQ(ierr);
    if (nq>1) {
      ierr = PetscMalloc2(nq-1,&wr,nq-1,&wi);CHKERRQ(ierr);
      ierr = NEPNLEIGSAuxiliarPRootFinder(nq-1,qcoeff,wr,wi,&avail);CHKERRQ(ierr);
      if (avail) {
        ierr = PetscCalloc1(nq-1,isol);CHKERRQ(ierr);
        *nisol = 0;
        for (i=0;i<nq-1;i++)
#if !defined(PETSC_USE_COMPLEX)
          if (wi[i]==0)
#endif
            (*isol)[(*nisol)++] = wr[i];
        nq = *nisol; *nisol = 0;
        for (i=0;i<nq;i++) wr[i] = (*isol)[i];
        ierr = NEPNLEIGSAuxiliarRmDuplicates(nq,wr,nisol,*isol,0);CHKERRQ(ierr);
        ierr = PetscFree2(wr,wi);CHKERRQ(ierr);
      } else { *nisol=0; *isol = NULL; }
    } else { *nisol = 0; *isol = NULL; }
    ierr = PetscFree(qcoeff);CHKERRQ(ierr);
  }
  ierr = PetscObjectTypeCompare((PetscObject)f,FNCOMBINE,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = FNCombineGetChildren(f,&ctype,&f1,&f2);CHKERRQ(ierr);
    if (ctype != FN_COMBINE_COMPOSE && ctype != FN_COMBINE_DIVIDE) {
      ierr = NEPNLEIGSFNSingularities(f1,&nisol1,&isol1,&rat1);CHKERRQ(ierr);
      ierr = NEPNLEIGSFNSingularities(f2,&nisol2,&isol2,&rat2);CHKERRQ(ierr);
      if (nisol1+nisol2>0) {
        ierr = PetscCalloc1(nisol1+nisol2,isol);CHKERRQ(ierr);
        *nisol = 0;
        ierr = NEPNLEIGSAuxiliarRmDuplicates(nisol1,isol1,nisol,*isol,0);CHKERRQ(ierr);
        ierr = NEPNLEIGSAuxiliarRmDuplicates(nisol2,isol2,nisol,*isol,0);CHKERRQ(ierr);
      }
      *rational = (rat1&&rat2)?PETSC_TRUE:PETSC_FALSE;
      ierr = PetscFree(isol1);CHKERRQ(ierr);
      ierr = PetscFree(isol2);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSRationalSingularities(NEP nep,PetscInt *ndptx,PetscScalar *dxi,PetscBool *rational)
{
  PetscErrorCode ierr;
  PetscInt       nt,i,nisol;
  FN             f;
  PetscScalar    *isol;
  PetscBool      rat;

  PetscFunctionBegin;
  *rational = PETSC_TRUE;
  *ndptx = 0;
  ierr = NEPGetSplitOperatorInfo(nep,&nt,NULL);CHKERRQ(ierr);
  for (i=0;i<nt;i++) {
    ierr = NEPGetSplitOperatorTerm(nep,i,NULL,&f);CHKERRQ(ierr);
    ierr = NEPNLEIGSFNSingularities(f,&nisol,&isol,&rat);CHKERRQ(ierr);
    if (nisol) {
      ierr = NEPNLEIGSAuxiliarRmDuplicates(nisol,isol,ndptx,dxi,0);CHKERRQ(ierr);
      ierr = PetscFree(isol);CHKERRQ(ierr);
    }
    *rational = ((*rational)&&rat)?PETSC_TRUE:PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}

/*  Adaptive Anderson-Antoulas algorithm */
static PetscErrorCode NEPNLEIGSAAAComputation(NEP nep,PetscInt ndpt,PetscScalar *ds,PetscScalar *F,PetscInt *ndptx,PetscScalar *dxi)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscScalar    mean=0.0,*z,*f,*C,*A,*VT,*work,*ww,szero=0.0,sone=1.0;
  PetscScalar    *N,*D;
  PetscReal      *S,norm,err,*R;
  PetscInt       i,k,j,idx=0,cont;
  PetscBLASInt   n_,m_,lda_,lwork,info,one=1;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      *rwork;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(8*ndpt,&lwork);CHKERRQ(ierr);
  ierr = PetscMalloc5(ndpt,&R,ndpt,&z,ndpt,&f,ndpt*ndpt,&C,ndpt,&ww);CHKERRQ(ierr);
  ierr = PetscMalloc6(ndpt*ndpt,&A,ndpt,&S,ndpt*ndpt,&VT,lwork,&work,ndpt,&D,ndpt,&N);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscMalloc1(8*ndpt,&rwork);CHKERRQ(ierr);
#endif
  norm = 0.0;
  for (i=0;i<ndpt;i++) {
    mean += F[i];
    norm = PetscMax(PetscAbsScalar(F[i]),norm);
  }
  mean /= ndpt;
  ierr = PetscBLASIntCast(ndpt,&lda_);CHKERRQ(ierr);
  for (i=0;i<ndpt;i++) R[i] = PetscAbsScalar(F[i]-mean);
  /* next support point */
  err = 0.0;
  for (i=0;i<ndpt;i++) if (R[i]>=err) {idx = i; err = R[i];}
  for (k=0;k<ndpt-1;k++) {
    z[k] = ds[idx]; f[k] = F[idx]; R[idx] = -1.0;
    /* next column of Cauchy matrix */
    for (i=0;i<ndpt;i++) {
      C[i+k*ndpt] = 1.0/(ds[i]-ds[idx]);
    }

    ierr = PetscArrayzero(A,ndpt*ndpt);CHKERRQ(ierr);
    cont = 0;
    for (i=0;i<ndpt;i++) {
      if (R[i]!=-1.0) {
        for (j=0;j<=k;j++)A[cont+j*ndpt] = C[i+j*ndpt]*F[i]-C[i+j*ndpt]*f[j];
        cont++;
      }
    }
    ierr = PetscBLASIntCast(cont,&m_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(k+1,&n_);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","A",&m_,&n_,A,&lda_,S,NULL,&lda_,VT,&lda_,work,&lwork,rwork,&info));
#else
    PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("N","A",&m_,&n_,A,&lda_,S,NULL,&lda_,VT,&lda_,work,&lwork,&info));
#endif
    SlepcCheckLapackInfo("gesvd",info);
    for (i=0;i<=k;i++) {
      ww[i] = PetscConj(VT[i*ndpt+k]);
      D[i] = ww[i]*f[i];
    }
    PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&lda_,&n_,&sone,C,&lda_,D,&one,&szero,N,&one));
    PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&lda_,&n_,&sone,C,&lda_,ww,&one,&szero,D,&one));
    for (i=0;i<ndpt;i++) if (R[i]>=0) R[i] = PetscAbsScalar(F[i]-N[i]/D[i]);
    /* next support point */
    err = 0.0;
    for (i=0;i<ndpt;i++) if (R[i]>=err) {idx = i; err = R[i];}
    if (err <= ctx->ddtol*norm) break;
  }

  if (k==ndpt-1) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Failed to determine singularities automatically in general problem");
  /* poles */
  ierr = PetscArrayzero(C,ndpt*ndpt);CHKERRQ(ierr);
  ierr = PetscArrayzero(A,ndpt*ndpt);CHKERRQ(ierr);
  for (i=0;i<=k;i++) {
    C[i+ndpt*i] = 1.0;
    A[(i+1)*ndpt] = ww[i];
    A[i+1] = 1.0;
    A[i+1+(i+1)*ndpt] = z[i];
  }
  C[0] = 0.0; C[k+1+(k+1)*ndpt] = 1.0;
  n_++;
#if defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKggev",LAPACKggev_("N","N",&n_,A,&lda_,C,&lda_,D,N,NULL,&lda_,NULL,&lda_,work,&lwork,rwork,&info));
#else
  PetscStackCallBLAS("LAPACKggev",LAPACKggev_("N","N",&n_,A,&lda_,C,&lda_,D,VT,N,NULL,&lda_,NULL,&lda_,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("ggev",info);
  cont = 0.0;
  for (i=0;i<n_;i++) if (N[i]!=0.0) {
    dxi[cont++] = D[i]/N[i];
  }
  *ndptx = cont;
  ierr = PetscFree5(R,z,f,C,ww);CHKERRQ(ierr);
  ierr = PetscFree6(A,S,VT,work,D,N);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscFree(rwork);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

/*  Singularities using Adaptive Anderson-Antoulas algorithm */
static PetscErrorCode NEPNLEIGSAAASingularities(NEP nep,PetscInt ndpt,PetscScalar *ds,PetscInt *ndptx,PetscScalar *dxi)
{
  PetscErrorCode ierr;
  Vec            u,v,w;
  PetscRandom    rand;
  PetscScalar    *F,*isol;
  PetscInt       i,k,nisol,nt;
  Mat            T;
  FN             f;

  PetscFunctionBegin;
  ierr = PetscMalloc1(ndpt,&F);CHKERRQ(ierr);
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    ierr = PetscMalloc1(ndpt,&isol);CHKERRQ(ierr);
    *ndptx = 0;
    ierr = NEPGetSplitOperatorInfo(nep,&nt,NULL);CHKERRQ(ierr);
    nisol = *ndptx;
    for (k=0;k<nt;k++) {
      ierr = NEPGetSplitOperatorTerm(nep,k,NULL,&f);CHKERRQ(ierr);
      for (i=0;i<ndpt;i++) {
        ierr = FNEvaluateFunction(f,ds[i],&F[i]);CHKERRQ(ierr);
      }
      ierr = NEPNLEIGSAAAComputation(nep,ndpt,ds,F,&nisol,isol);CHKERRQ(ierr);
      if (nisol) {
        ierr = NEPNLEIGSAuxiliarRmDuplicates(nisol,isol,ndptx,dxi,ndpt);CHKERRQ(ierr);
      }
    }
    ierr = PetscFree(isol);CHKERRQ(ierr);
  } else {
    ierr = MatCreateVecs(nep->function,&u,NULL);CHKERRQ(ierr);
    ierr = VecDuplicate(u,&v);CHKERRQ(ierr);
    ierr = VecDuplicate(u,&w);CHKERRQ(ierr);
    ierr = BVGetRandomContext(nep->V,&rand);CHKERRQ(ierr);
    ierr = VecSetRandom(u,rand);CHKERRQ(ierr);
    ierr = VecNormalize(u,NULL);CHKERRQ(ierr);
    ierr = VecSetRandom(v,rand);CHKERRQ(ierr);
    ierr = VecNormalize(v,NULL);CHKERRQ(ierr);
    T = nep->function;
    for (i=0;i<ndpt;i++) {
      ierr = NEPComputeFunction(nep,ds[i],T,T);CHKERRQ(ierr);
      ierr = MatMult(T,v,w);CHKERRQ(ierr);
      ierr = VecDot(w,u,&F[i]);CHKERRQ(ierr);
    }
    ierr = NEPNLEIGSAAAComputation(nep,ndpt,ds,F,ndptx,dxi);CHKERRQ(ierr);
    ierr = VecDestroy(&u);CHKERRQ(ierr);
    ierr = VecDestroy(&v);CHKERRQ(ierr);
    ierr = VecDestroy(&w);CHKERRQ(ierr);
  }
  ierr = PetscFree(F);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSLejaBagbyPoints(NEP nep)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       i,k,ndpt=NDPOINTS,ndptx=NDPOINTS;
  PetscScalar    *ds,*dsi,*dxi,*nrs,*nrxi,*s=ctx->s,*xi=ctx->xi,*beta=ctx->beta;
  PetscReal      maxnrs,minnrxi;
  PetscBool      rational;
#if !defined(PETSC_USE_COMPLEX)
  PetscReal      a,b,h;
#endif

  PetscFunctionBegin;
  if (!ctx->computesingularities && nep->problem_type!=NEP_RATIONAL) ndpt = ndptx = LBPOINTS;
  ierr = PetscMalloc5(ndpt+1,&ds,ndpt+1,&dsi,ndpt,&dxi,ndpt+1,&nrs,ndpt,&nrxi);CHKERRQ(ierr);

  /* Discretize the target region boundary */
  ierr = RGComputeContour(nep->rg,ndpt,ds,dsi);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  for (i=0;i<ndpt;i++) if (dsi[i]!=0.0) break;
  if (i<ndpt) {
    if (nep->problem_type==NEP_RATIONAL) {
      /* Select a segment in the real axis */
      ierr = RGComputeBoundingBox(nep->rg,&a,&b,NULL,NULL);CHKERRQ(ierr);
      if (a<=-PETSC_MAX_REAL || b>=PETSC_MAX_REAL) SETERRQ(PetscObjectComm((PetscObject)nep),1,"NLEIGS requires a bounded target set");
      h = (b-a)/ndpt;
      for (i=0;i<ndpt;i++) {ds[i] = a+h*i; dsi[i] = 0.0;}
    } else SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"NLEIGS with real arithmetic requires the target set to be included in the real axis");
  }
#endif
  /* Discretize the singularity region */
  if (ctx->computesingularities) {
    ierr = (ctx->computesingularities)(nep,&ndptx,dxi,ctx->singularitiesctx);CHKERRQ(ierr);
  } else {
    if (nep->problem_type==NEP_RATIONAL) {
      ierr = NEPNLEIGSRationalSingularities(nep,&ndptx,dxi,&rational);CHKERRQ(ierr);
      if (!rational) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Failed to determine singularities automatically in rational problem; consider solving the problem as general");
    } else {
      /* AAA algorithm */
      ierr = NEPNLEIGSAAASingularities(nep,ndpt,ds,&ndptx,dxi);CHKERRQ(ierr);
    }
  }
  /* Look for Leja-Bagby points in the discretization sets */
  s[0]    = ds[0];
  xi[0]   = (ndptx>0)?dxi[0]:PETSC_INFINITY;
  if (PetscAbsScalar(xi[0])<10*PETSC_MACHINE_EPSILON) SETERRQ2(PetscObjectComm((PetscObject)nep),1,"Singularity point %D is nearly zero: %g; consider removing the singularity or shifting the problem",0,(double)PetscAbsScalar(xi[0]));
  beta[0] = 1.0; /* scaling factors are also computed here */
  for (i=0;i<ndpt;i++) {
    nrs[i] = 1.0;
    nrxi[i] = 1.0;
  }
  for (k=1;k<ctx->ddmaxit;k++) {
    maxnrs = 0.0;
    minnrxi = PETSC_MAX_REAL;
    for (i=0;i<ndpt;i++) {
      nrs[i] *= ((ds[i]-s[k-1])/(1.0-ds[i]/xi[k-1]))/beta[k-1];
      if (PetscAbsScalar(nrs[i])>maxnrs) {maxnrs = PetscAbsScalar(nrs[i]); s[k] = ds[i];}
    }
    if (ndptx>k) {
      for (i=1;i<ndptx;i++) {
        nrxi[i] *= ((dxi[i]-s[k-1])/(1.0-dxi[i]/xi[k-1]))/beta[k-1];
        if (PetscAbsScalar(nrxi[i])<minnrxi) {minnrxi = PetscAbsScalar(nrxi[i]); xi[k] = dxi[i];}
      }
      if (PetscAbsScalar(xi[k])<10*PETSC_MACHINE_EPSILON) SETERRQ2(PetscObjectComm((PetscObject)nep),1,"Singularity point %D is nearly zero: %g; consider removing the singularity or shifting the problem",k,(double)PetscAbsScalar(xi[k]));
    } else xi[k] = PETSC_INFINITY;
    beta[k] = maxnrs;
  }
  ierr = PetscFree5(ds,dsi,dxi,nrs,nrxi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPNLEIGSEvalNRTFunct(NEP nep,PetscInt k,PetscScalar sigma,PetscScalar *b)
{
  NEP_NLEIGS  *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt    i;
  PetscScalar *beta=ctx->beta,*s=ctx->s,*xi=ctx->xi;

  PetscFunctionBegin;
  b[0] = 1.0/beta[0];
  for (i=0;i<k;i++) {
    b[i+1] = ((sigma-s[i])*b[i])/(beta[i+1]*(1.0-sigma/xi[i]));
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_Fun(Mat A,Vec x,Vec y)
{
  PetscErrorCode      ierr;
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            i;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = MatMult(ctx->A[0],x,y);CHKERRQ(ierr);
  if (ctx->coeff[0]!=1.0) { ierr = VecScale(y,ctx->coeff[0]);CHKERRQ(ierr); }
  for (i=1;i<ctx->nmat;i++) {
    ierr = MatMult(ctx->A[i],x,ctx->t);CHKERRQ(ierr);
    ierr = VecAXPY(y,ctx->coeff[i],ctx->t);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMultTranspose_Fun(Mat A,Vec x,Vec y)
{
  PetscErrorCode      ierr;
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            i;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = MatMultTranspose(ctx->A[0],x,y);CHKERRQ(ierr);
  if (ctx->coeff[0]!=1.0) { ierr = VecScale(y,ctx->coeff[0]);CHKERRQ(ierr); }
  for (i=1;i<ctx->nmat;i++) {
    ierr = MatMultTranspose(ctx->A[i],x,ctx->t);CHKERRQ(ierr);
    ierr = VecAXPY(y,ctx->coeff[i],ctx->t);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatGetDiagonal_Fun(Mat A,Vec diag)
{
  PetscErrorCode      ierr;
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            i;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = MatGetDiagonal(ctx->A[0],diag);CHKERRQ(ierr);
  if (ctx->coeff[0]!=1.0) { ierr = VecScale(diag,ctx->coeff[0]);CHKERRQ(ierr); }
  for (i=1;i<ctx->nmat;i++) {
    ierr = MatGetDiagonal(ctx->A[i],ctx->t);CHKERRQ(ierr);
    ierr = VecAXPY(diag,ctx->coeff[i],ctx->t);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDuplicate_Fun(Mat A,MatDuplicateOption op,Mat *B)
{
  PetscInt            n,i;
  NEP_NLEIGS_MATSHELL *ctxnew,*ctx;
  void                (*fun)(void);
  PetscErrorCode      ierr;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = PetscNew(&ctxnew);CHKERRQ(ierr);
  ctxnew->nmat = ctx->nmat;
  ctxnew->maxnmat = ctx->maxnmat;
  ierr = PetscMalloc2(ctxnew->maxnmat,&ctxnew->A,ctxnew->maxnmat,&ctxnew->coeff);CHKERRQ(ierr);
  for (i=0;i<ctx->nmat;i++) {
    ierr = PetscObjectReference((PetscObject)ctx->A[i]);CHKERRQ(ierr);
    ctxnew->A[i] = ctx->A[i];
    ctxnew->coeff[i] = ctx->coeff[i];
  }
  ierr = MatGetSize(ctx->A[0],&n,NULL);CHKERRQ(ierr);
  ierr = VecDuplicate(ctx->t,&ctxnew->t);CHKERRQ(ierr);
  ierr = MatCreateShell(PetscObjectComm((PetscObject)A),n,n,n,n,(void*)ctxnew,B);CHKERRQ(ierr);
  ierr = MatShellSetManageScalingShifts(*B);CHKERRQ(ierr);
  ierr = MatShellGetOperation(A,MATOP_MULT,&fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_MULT,fun);CHKERRQ(ierr);
  ierr = MatShellGetOperation(A,MATOP_MULT_TRANSPOSE,&fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_MULT_TRANSPOSE,fun);CHKERRQ(ierr);
  ierr = MatShellGetOperation(A,MATOP_GET_DIAGONAL,&fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_GET_DIAGONAL,fun);CHKERRQ(ierr);
  ierr = MatShellGetOperation(A,MATOP_DUPLICATE,&fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_DUPLICATE,fun);CHKERRQ(ierr);
  ierr = MatShellGetOperation(A,MATOP_DESTROY,&fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_DESTROY,fun);CHKERRQ(ierr);
  ierr = MatShellGetOperation(A,MATOP_AXPY,&fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_AXPY,fun);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDestroy_Fun(Mat A)
{
  NEP_NLEIGS_MATSHELL *ctx;
  PetscErrorCode      ierr;
  PetscInt            i;

  PetscFunctionBeginUser;
  if (A) {
    ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
    for (i=0;i<ctx->nmat;i++) {
      ierr = MatDestroy(&ctx->A[i]);CHKERRQ(ierr);
    }
    ierr = VecDestroy(&ctx->t);CHKERRQ(ierr);
    ierr = PetscFree2(ctx->A,ctx->coeff);CHKERRQ(ierr);
    ierr = PetscFree(ctx);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatAXPY_Fun(Mat Y,PetscScalar a,Mat X,MatStructure str)
{
  NEP_NLEIGS_MATSHELL *ctxY,*ctxX;
  PetscErrorCode      ierr;
  PetscInt            i,j;
  PetscBool           found;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(Y,(void**)&ctxY);CHKERRQ(ierr);
  ierr = MatShellGetContext(X,(void**)&ctxX);CHKERRQ(ierr);
  for (i=0;i<ctxX->nmat;i++) {
    found = PETSC_FALSE;
    for (j=0;!found&&j<ctxY->nmat;j++) {
      if (ctxX->A[i]==ctxY->A[j]) {
        found = PETSC_TRUE;
        ctxY->coeff[j] += a*ctxX->coeff[i];
      }
    }
    if (!found) {
      ctxY->coeff[ctxY->nmat] = a*ctxX->coeff[i];
      ctxY->A[ctxY->nmat++] = ctxX->A[i];
      ierr = PetscObjectReference((PetscObject)ctxX->A[i]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatScale_Fun(Mat M,PetscScalar a)
{
  NEP_NLEIGS_MATSHELL *ctx;
  PetscErrorCode      ierr;
  PetscInt            i;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(M,(void**)&ctx);CHKERRQ(ierr);
  for (i=0;i<ctx->nmat;i++) ctx->coeff[i] *= a;
  PetscFunctionReturn(0);
}

static PetscErrorCode NLEIGSMatToMatShellArray(Mat M,Mat *Ms,PetscInt maxnmat)
{
  PetscErrorCode      ierr;
  NEP_NLEIGS_MATSHELL *ctx;
  PetscInt            n;
  PetscBool           has;

  PetscFunctionBegin;
  ierr = MatHasOperation(M,MATOP_DUPLICATE,&has);CHKERRQ(ierr);
  if (!has) SETERRQ(PetscObjectComm((PetscObject)M),1,"MatDuplicate operation required");
  ierr = PetscNew(&ctx);CHKERRQ(ierr);
  ctx->maxnmat = maxnmat;
  ierr = PetscMalloc2(ctx->maxnmat,&ctx->A,ctx->maxnmat,&ctx->coeff);CHKERRQ(ierr);
  ierr = MatDuplicate(M,MAT_COPY_VALUES,&ctx->A[0]);CHKERRQ(ierr);
  ctx->nmat = 1;
  ctx->coeff[0] = 1.0;
  ierr = MatCreateVecs(M,&ctx->t,NULL);CHKERRQ(ierr);
  ierr = MatGetSize(M,&n,NULL);CHKERRQ(ierr);
  ierr = MatCreateShell(PetscObjectComm((PetscObject)M),n,n,n,n,(void*)ctx,Ms);CHKERRQ(ierr);
  ierr = MatShellSetManageScalingShifts(*Ms);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Ms,MATOP_MULT,(void(*)(void))MatMult_Fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Ms,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMultTranspose_Fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Ms,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Ms,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_Fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Ms,MATOP_DESTROY,(void(*)(void))MatDestroy_Fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Ms,MATOP_AXPY,(void(*)(void))MatAXPY_Fun);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*Ms,MATOP_SCALE,(void(*)(void))MatScale_Fun);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSNormEstimation(NEP nep,Mat M,PetscReal *norm,Vec *w)
{
  PetscScalar    *z,*x,*y;
  PetscReal      tr;
  Vec            X=w[0],Y=w[1];
  PetscInt       n,i;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscRandom    rand;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!ctx->vrn) {
    /* generate a random vector with normally distributed entries with the Box-Muller transform */
    ierr = BVGetRandomContext(nep->V,&rand);CHKERRQ(ierr);
    ierr = MatCreateVecs(M,&ctx->vrn,NULL);CHKERRQ(ierr);
    ierr = VecSetRandom(X,rand);CHKERRQ(ierr);
    ierr = VecSetRandom(Y,rand);CHKERRQ(ierr);
    ierr = VecGetLocalSize(ctx->vrn,&n);CHKERRQ(ierr);
    ierr = VecGetArray(ctx->vrn,&z);CHKERRQ(ierr);
    ierr = VecGetArray(X,&x);CHKERRQ(ierr);
    ierr = VecGetArray(Y,&y);CHKERRQ(ierr);
    for (i=0;i<n;i++) {
#if defined(PETSC_USE_COMPLEX)
      z[i] = PetscCMPLX(PetscSqrtReal(-2.0*PetscLogReal(PetscRealPart(x[i])))*PetscCosReal(2.0*PETSC_PI*PetscRealPart(y[i])),PetscSqrtReal(-2.0*PetscLogReal(PetscImaginaryPart(x[i])))*PetscCosReal(2.0*PETSC_PI*PetscImaginaryPart(y[i])));
#else
      z[i] = PetscSqrtReal(-2.0*PetscLogReal(x[i]))*PetscCosReal(2.0*PETSC_PI*y[i]);
#endif
    }
    ierr = VecRestoreArray(ctx->vrn,&z);CHKERRQ(ierr);
    ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
    ierr = VecRestoreArray(Y,&y);CHKERRQ(ierr);
    ierr = VecNorm(ctx->vrn,NORM_2,&tr);CHKERRQ(ierr);
    ierr = VecScale(ctx->vrn,1/tr);CHKERRQ(ierr);
  }
  /* matrix-free norm estimator of Ipsen https://ipsen.math.ncsu.edu/ps/slides_ima.pdf */
  ierr = MatGetSize(M,&n,NULL);CHKERRQ(ierr);
  ierr = MatMult(M,ctx->vrn,X);CHKERRQ(ierr);
  ierr = VecNorm(X,NORM_2,norm);CHKERRQ(ierr);
  *norm *= PetscSqrtReal((PetscReal)n);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSDividedDifferences_split(NEP nep)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       k,j,i,maxnmat,nmax;
  PetscReal      norm0,norm,*matnorm;
  PetscScalar    *s=ctx->s,*beta=ctx->beta,*xi=ctx->xi,*b,alpha,*coeffs,*pK,*pH,sone=1.0;
  Mat            T,Ts,K,H;
  PetscBool      shell,hasmnorm=PETSC_FALSE,matrix=PETSC_TRUE;
  PetscBLASInt   n_;

  PetscFunctionBegin;
  nmax = ctx->ddmaxit;
  ierr = PetscMalloc1(nep->nt*nmax,&ctx->coeffD);CHKERRQ(ierr);
  ierr = PetscMalloc3(nmax+1,&b,nmax+1,&coeffs,nep->nt,&matnorm);CHKERRQ(ierr);
  for (j=0;j<nep->nt;j++) {
    ierr = MatHasOperation(nep->A[j],MATOP_NORM,&hasmnorm);CHKERRQ(ierr);
    if (!hasmnorm) break;
    ierr = MatNorm(nep->A[j],NORM_INFINITY,matnorm+j);CHKERRQ(ierr);
  }
  /* Try matrix functions scheme */
  ierr = PetscCalloc2(nmax*nmax,&pK,nmax*nmax,&pH);CHKERRQ(ierr);
  for (i=0;i<nmax-1;i++) {
    pK[(nmax+1)*i]   = 1.0;
    pK[(nmax+1)*i+1] = beta[i+1]/xi[i];
    pH[(nmax+1)*i]   = s[i];
    pH[(nmax+1)*i+1] = beta[i+1];
  }
  pH[nmax*nmax-1] = s[nmax-1];
  pK[nmax*nmax-1] = 1.0;
  ierr = PetscBLASIntCast(nmax,&n_);CHKERRQ(ierr);
  PetscStackCallBLAS("BLAStrsm",BLAStrsm_("R","L","N","U",&n_,&n_,&sone,pK,&n_,pH,&n_));
  /* The matrix to be used is in H. K will be a work-space matrix */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,nmax,nmax,pH,&H);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,nmax,nmax,pK,&K);CHKERRQ(ierr);
  for (j=0;matrix&&j<nep->nt;j++) {
    PetscPushErrorHandler(PetscIgnoreErrorHandler,NULL);
    ierr = FNEvaluateFunctionMat(nep->f[j],H,K);
    PetscPopErrorHandler();
    if (!ierr) {
      for (i=0;i<nmax;i++) { ctx->coeffD[j+i*nep->nt] = pK[i]*beta[0]; }
    } else {
      matrix = PETSC_FALSE;
      ierr = PetscFPTrapPop();CHKERRQ(ierr);
    }
  }
  ierr = MatDestroy(&H);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  if (!matrix) {
    for (j=0;j<nep->nt;j++) {
      ierr = FNEvaluateFunction(nep->f[j],s[0],ctx->coeffD+j);CHKERRQ(ierr);
      ctx->coeffD[j] *= beta[0];
    }
  }
  if (hasmnorm) {
    norm0 = 0.0;
    for (j=0;j<nep->nt;j++) norm0 += matnorm[j]*PetscAbsScalar(ctx->coeffD[j]);
  } else {
    norm0 = 0.0;
    for (j=0;j<nep->nt;j++) norm0 = PetscMax(PetscAbsScalar(ctx->coeffD[j]),norm0);
  }
  ctx->nmat = ctx->ddmaxit;
  for (k=1;k<ctx->ddmaxit;k++) {
    if (!matrix) {
      ierr = NEPNLEIGSEvalNRTFunct(nep,k,s[k],b);CHKERRQ(ierr);
      for (i=0;i<nep->nt;i++) {
        ierr = FNEvaluateFunction(nep->f[i],s[k],ctx->coeffD+k*nep->nt+i);CHKERRQ(ierr);
        for (j=0;j<k;j++) {
          ctx->coeffD[k*nep->nt+i] -= b[j]*ctx->coeffD[i+nep->nt*j];
        }
        ctx->coeffD[k*nep->nt+i] /= b[k];
      }
    }
    if (hasmnorm) {
      norm = 0.0;
      for (j=0;j<nep->nt;j++) norm += matnorm[j]*PetscAbsScalar(ctx->coeffD[k*nep->nt+j]);
    } else {
      norm = 0.0;
      for (j=0;j<nep->nt;j++) norm = PetscMax(PetscAbsScalar(ctx->coeffD[k*nep->nt+j]),norm);
    }
    if (k>1 && norm/norm0 < ctx->ddtol) {
      ctx->nmat = k+1;
      break;
    }
  }
  if (!ctx->ksp) { ierr = NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp);CHKERRQ(ierr); }
  ierr = PetscObjectTypeCompare((PetscObject)nep->A[0],MATSHELL,&shell);CHKERRQ(ierr);
  maxnmat = PetscMax(ctx->ddmaxit,nep->nt);
  for (i=0;i<ctx->nshiftsw;i++) {
    ierr = NEPNLEIGSEvalNRTFunct(nep,ctx->nmat-1,ctx->shifts[i],coeffs);CHKERRQ(ierr);
    if (!shell) {
      ierr = MatDuplicate(nep->A[0],MAT_COPY_VALUES,&T);CHKERRQ(ierr);
    } else {
      ierr = NLEIGSMatToMatShellArray(nep->A[0],&T,maxnmat);CHKERRQ(ierr);
    }
    alpha = 0.0;
    for (j=0;j<ctx->nmat;j++) alpha += coeffs[j]*ctx->coeffD[j*nep->nt];
    ierr = MatScale(T,alpha);CHKERRQ(ierr);
    for (k=1;k<nep->nt;k++) {
      alpha = 0.0;
      for (j=0;j<ctx->nmat;j++) alpha += coeffs[j]*ctx->coeffD[j*nep->nt+k];
      if (shell) {
        ierr = NLEIGSMatToMatShellArray(nep->A[k],&Ts,maxnmat);CHKERRQ(ierr);
      }
      ierr = MatAXPY(T,alpha,shell?Ts:nep->A[k],nep->mstr);CHKERRQ(ierr);
      if (shell) {
        ierr = MatDestroy(&Ts);CHKERRQ(ierr);
      }
    }
    ierr = KSPSetOperators(ctx->ksp[i],T,T);CHKERRQ(ierr);
    ierr = KSPSetUp(ctx->ksp[i]);CHKERRQ(ierr);
    ierr = MatDestroy(&T);CHKERRQ(ierr);
  }
  ierr = PetscFree3(b,coeffs,matnorm);CHKERRQ(ierr);
  ierr = PetscFree2(pK,pH);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSDividedDifferences_callback(NEP nep)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       k,j,i,maxnmat;
  PetscReal      norm0,norm;
  PetscScalar    *s=ctx->s,*beta=ctx->beta,*b,*coeffs;
  Mat            *D=ctx->D,T;
  PetscBool      shell,has,vec=PETSC_FALSE;
  Vec            w[2];

  PetscFunctionBegin;
  ierr = PetscMalloc2(ctx->ddmaxit+1,&b,ctx->ddmaxit+1,&coeffs);CHKERRQ(ierr);
  T = nep->function;
  ierr = NEPComputeFunction(nep,s[0],T,T);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)T,MATSHELL,&shell);CHKERRQ(ierr);
  maxnmat = PetscMax(ctx->ddmaxit,nep->nt);
  if (!shell) {
    ierr = MatDuplicate(T,MAT_COPY_VALUES,&D[0]);CHKERRQ(ierr);
  } else {
    ierr = NLEIGSMatToMatShellArray(T,&D[0],maxnmat);CHKERRQ(ierr);
  }
  if (beta[0]!=1.0) {
    ierr = MatScale(D[0],1.0/beta[0]);CHKERRQ(ierr);
  }
  ierr = MatHasOperation(D[0],MATOP_NORM,&has);CHKERRQ(ierr);
  if (has) {
    ierr = MatNorm(D[0],NORM_FROBENIUS,&norm0);CHKERRQ(ierr);
  } else {
    ierr = MatCreateVecs(D[0],NULL,&w[0]);CHKERRQ(ierr);
    ierr = VecDuplicate(w[0],&w[1]);CHKERRQ(ierr);
    vec = PETSC_TRUE;
    ierr = NEPNLEIGSNormEstimation(nep,D[0],&norm0,w);CHKERRQ(ierr);
  }
  ctx->nmat = ctx->ddmaxit;
  for (k=1;k<ctx->ddmaxit;k++) {
    ierr = NEPNLEIGSEvalNRTFunct(nep,k,s[k],b);CHKERRQ(ierr);
    ierr = NEPComputeFunction(nep,s[k],T,T);CHKERRQ(ierr);
    if (!shell) {
      ierr = MatDuplicate(T,MAT_COPY_VALUES,&D[k]);CHKERRQ(ierr);
    } else {
      ierr = NLEIGSMatToMatShellArray(T,&D[k],maxnmat);CHKERRQ(ierr);
    }
    for (j=0;j<k;j++) {
      ierr = MatAXPY(D[k],-b[j],D[j],nep->mstr);CHKERRQ(ierr);
    }
    ierr = MatScale(D[k],1.0/b[k]);CHKERRQ(ierr);
    ierr = MatHasOperation(D[k],MATOP_NORM,&has);CHKERRQ(ierr);
    if (has) {
      ierr = MatNorm(D[k],NORM_FROBENIUS,&norm);CHKERRQ(ierr);
    } else {
      if (!vec) {
        ierr = MatCreateVecs(D[k],NULL,&w[0]);CHKERRQ(ierr);
        ierr = VecDuplicate(w[0],&w[1]);CHKERRQ(ierr);
        vec = PETSC_TRUE;
      }
      ierr = NEPNLEIGSNormEstimation(nep,D[k],&norm,w);CHKERRQ(ierr);
    }
    if (k>1 && norm/norm0 < ctx->ddtol && k>1) {
      ctx->nmat = k+1;
      break;
    }
  }
  if (!ctx->ksp) { ierr = NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp);CHKERRQ(ierr); }
  for (i=0;i<ctx->nshiftsw;i++) {
    ierr = NEPNLEIGSEvalNRTFunct(nep,ctx->nmat-1,ctx->shifts[i],coeffs);CHKERRQ(ierr);
    ierr = MatDuplicate(ctx->D[0],MAT_COPY_VALUES,&T);CHKERRQ(ierr);
    if (coeffs[0]!=1.0) { ierr = MatScale(T,coeffs[0]);CHKERRQ(ierr); }
    for (j=1;j<ctx->nmat;j++) {
      ierr = MatAXPY(T,coeffs[j],ctx->D[j],nep->mstr);CHKERRQ(ierr);
    }
    ierr = KSPSetOperators(ctx->ksp[i],T,T);CHKERRQ(ierr);
    ierr = KSPSetUp(ctx->ksp[i]);CHKERRQ(ierr);
    ierr = MatDestroy(&T);CHKERRQ(ierr);
  }
  ierr = PetscFree2(b,coeffs);CHKERRQ(ierr);
  if (vec) {
    ierr = VecDestroy(&w[0]);CHKERRQ(ierr);
    ierr = VecDestroy(&w[1]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   NEPKrylovConvergence - This is the analogue to EPSKrylovConvergence.
*/
PetscErrorCode NEPNLEIGSKrylovConvergence(NEP nep,PetscBool getall,PetscInt kini,PetscInt nits,PetscReal betah,PetscScalar betak,PetscInt *kout,Vec *w)
{
  PetscErrorCode ierr;
  PetscInt       k,newk,marker,inside;
  PetscScalar    re,im;
  PetscReal      resnorm,tt;
  PetscBool      istrivial;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  ierr = RGIsTrivial(nep->rg,&istrivial);CHKERRQ(ierr);
  marker = -1;
  if (nep->trackall) getall = PETSC_TRUE;
  for (k=kini;k<kini+nits;k++) {
    /* eigenvalue */
    re = nep->eigr[k];
    im = nep->eigi[k];
    if (!istrivial) {
      if (!ctx->nshifts) {
        ierr = NEPNLEIGSBackTransform((PetscObject)nep,1,&re,&im);CHKERRQ(ierr);
      }
      ierr = RGCheckInside(nep->rg,1,&re,&im,&inside);CHKERRQ(ierr);
      if (marker==-1 && inside<0) marker = k;
    }
    newk = k;
    ierr = DSVectors(nep->ds,DS_MAT_X,&newk,&resnorm);CHKERRQ(ierr);
    tt = (ctx->nshifts)?SlepcAbsEigenvalue(betak-nep->eigr[k]*betah,nep->eigi[k]*betah):betah;
    resnorm *=  PetscAbsReal(tt);
    /* error estimate */
    ierr = (*nep->converged)(nep,nep->eigr[k],nep->eigi[k],resnorm,&nep->errest[k],nep->convergedctx);CHKERRQ(ierr);
    if (marker==-1 && nep->errest[k] >= nep->tol) marker = k;
    if (newk==k+1) {
      nep->errest[k+1] = nep->errest[k];
      k++;
    }
    if (marker!=-1 && !getall) break;
  }
  if (marker!=-1) k = marker;
  *kout = k;
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSetUp_NLEIGS(NEP nep)
{
  PetscErrorCode ierr;
  PetscInt       k,in;
  PetscScalar    zero=0.0;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  SlepcSC        sc;
  PetscBool      istrivial;

  PetscFunctionBegin;
  ierr = NEPSetDimensions_Default(nep,nep->nev,&nep->ncv,&nep->mpd);CHKERRQ(ierr);
  if (nep->ncv>nep->nev+nep->mpd) SETERRQ(PetscObjectComm((PetscObject)nep),1,"The value of ncv must not be larger than nev+mpd");
  if (nep->max_it==PETSC_DEFAULT) nep->max_it = PetscMax(5000,2*nep->n/nep->ncv);
  if (!ctx->ddmaxit) ctx->ddmaxit = LBPOINTS;
  ierr = RGIsTrivial(nep->rg,&istrivial);CHKERRQ(ierr);
  if (istrivial) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"NEPNLEIGS requires a nontrivial region defining the target set");
  if (!nep->which) nep->which = NEP_TARGET_MAGNITUDE;
  if (nep->which!=NEP_TARGET_MAGNITUDE && nep->which!=NEP_TARGET_REAL && nep->which!=NEP_TARGET_IMAGINARY && nep->which!=NEP_WHICH_USER) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Should set a target selection in NEPSetWhichEigenvalues()");

  /* Initialize the NLEIGS context structure */
  k = ctx->ddmaxit;
  ierr = PetscMalloc4(k,&ctx->s,k,&ctx->xi,k,&ctx->beta,k,&ctx->D);CHKERRQ(ierr);
  nep->data = ctx;
  if (nep->tol==PETSC_DEFAULT) nep->tol = SLEPC_DEFAULT_TOL;
  if (ctx->ddtol==PETSC_DEFAULT) ctx->ddtol = nep->tol/10.0;
  if (!ctx->keep) ctx->keep = 0.5;

  /* Compute Leja-Bagby points and scaling values */
  ierr = NEPNLEIGSLejaBagbyPoints(nep);CHKERRQ(ierr);
  if (nep->problem_type!=NEP_RATIONAL) {
    ierr = RGCheckInside(nep->rg,1,&nep->target,&zero,&in);CHKERRQ(ierr);
    if (in<0) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The target is not inside the target set");
  }

  /* Compute the divided difference matrices */
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    ierr = NEPNLEIGSDividedDifferences_split(nep);CHKERRQ(ierr);
  } else {
    ierr = NEPNLEIGSDividedDifferences_callback(nep);CHKERRQ(ierr);
  }
  ierr = NEPAllocateSolution(nep,ctx->nmat-1);CHKERRQ(ierr);
  ierr = NEPSetWorkVecs(nep,4);CHKERRQ(ierr);
  if (!ctx->fullbasis) {
    if (nep->twosided) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Two-sided variant requires the full-basis option, rerun with -nep_nleigs_full_basis");
    /* set-up DS and transfer split operator functions */
    ierr = DSSetType(nep->ds,ctx->nshifts?DSGNHEP:DSNHEP);CHKERRQ(ierr);
    ierr = DSAllocate(nep->ds,nep->ncv+1);CHKERRQ(ierr);
    ierr = DSGetSlepcSC(nep->ds,&sc);CHKERRQ(ierr);
    if (!ctx->nshifts) {
      sc->map = NEPNLEIGSBackTransform;
      ierr = DSSetExtraRow(nep->ds,PETSC_TRUE);CHKERRQ(ierr);
    }
    sc->mapobj        = (PetscObject)nep;
    sc->rg            = nep->rg;
    sc->comparison    = nep->sc->comparison;
    sc->comparisonctx = nep->sc->comparisonctx;
    ierr = BVDestroy(&ctx->V);CHKERRQ(ierr);
    ierr = BVCreateTensor(nep->V,ctx->nmat-1,&ctx->V);CHKERRQ(ierr);
    nep->ops->solve          = NEPSolve_NLEIGS;
    nep->ops->computevectors = NEPComputeVectors_Schur;
  } else {
    ierr = NEPSetUp_NLEIGS_FullBasis(nep);CHKERRQ(ierr);
    nep->ops->solve          = NEPSolve_NLEIGS_FullBasis;
    nep->ops->computevectors = NULL;
  }
  PetscFunctionReturn(0);
}

/*
  Extend the TOAR basis by applying the the matrix operator
  over a vector which is decomposed on the TOAR way
  Input:
    - S,V: define the latest Arnoldi vector (nv vectors in V)
  Output:
    - t: new vector extending the TOAR basis
    - r: temporally coefficients to compute the TOAR coefficients
         for the new Arnoldi vector
  Workspace: t_ (two vectors)
*/
static PetscErrorCode NEPTOARExtendBasis(NEP nep,PetscInt idxrktg,PetscScalar *S,PetscInt ls,PetscInt nv,BV W,BV V,Vec t,PetscScalar *r,PetscInt lr,Vec *t_)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       deg=ctx->nmat-1,k,j;
  Vec            v=t_[0],q=t_[1],w;
  PetscScalar    *beta=ctx->beta,*s=ctx->s,*xi=ctx->xi,*coeffs,sigma;

  PetscFunctionBegin;
  if (!ctx->ksp) { ierr = NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp);CHKERRQ(ierr); }
  sigma = ctx->shifts[idxrktg];
  ierr = BVSetActiveColumns(nep->V,0,nv);CHKERRQ(ierr);
  ierr = PetscMalloc1(ctx->nmat,&coeffs);CHKERRQ(ierr);
  if (PetscAbsScalar(s[deg-2]-sigma)<100*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_SELF,1,"Breakdown in NLEIGS");
  /* i-part stored in (i-1) position */
  for (j=0;j<nv;j++) {
    r[(deg-2)*lr+j] = (S[(deg-2)*ls+j]+(beta[deg-1]/xi[deg-2])*S[(deg-1)*ls+j])/(s[deg-2]-sigma);
  }
  ierr = BVSetActiveColumns(W,0,deg);CHKERRQ(ierr);
  ierr = BVGetColumn(W,deg-1,&w);CHKERRQ(ierr);
  ierr = BVMultVec(V,1.0/beta[deg],0,w,S+(deg-1)*ls);CHKERRQ(ierr);
  ierr = BVRestoreColumn(W,deg-1,&w);CHKERRQ(ierr);
  ierr = BVGetColumn(W,deg-2,&w);CHKERRQ(ierr);
  ierr = BVMultVec(V,1.0,0.0,w,r+(deg-2)*lr);CHKERRQ(ierr);
  ierr = BVRestoreColumn(W,deg-2,&w);CHKERRQ(ierr);
  for (k=deg-2;k>0;k--) {
    if (PetscAbsScalar(s[k-1]-sigma)<100*PETSC_MACHINE_EPSILON) SETERRQ(PETSC_COMM_SELF,1,"Breakdown in NLEIGS");
    for (j=0;j<nv;j++) r[(k-1)*lr+j] = (S[(k-1)*ls+j]+(beta[k]/xi[k-1])*S[k*ls+j]-beta[k]*(1.0-sigma/xi[k-1])*r[(k)*lr+j])/(s[k-1]-sigma);
    ierr = BVGetColumn(W,k-1,&w);CHKERRQ(ierr);
    ierr = BVMultVec(V,1.0,0.0,w,r+(k-1)*lr);CHKERRQ(ierr);
    ierr = BVRestoreColumn(W,k-1,&w);CHKERRQ(ierr);
  }
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    for (j=0;j<ctx->nmat-2;j++) coeffs[j] = ctx->coeffD[nep->nt*j];
    coeffs[ctx->nmat-2] = ctx->coeffD[nep->nt*(ctx->nmat-1)];
    ierr = BVMultVec(W,1.0,0.0,v,coeffs);CHKERRQ(ierr);
    ierr = MatMult(nep->A[0],v,q);CHKERRQ(ierr);
    for (k=1;k<nep->nt;k++) {
      for (j=0;j<ctx->nmat-2;j++) coeffs[j] = ctx->coeffD[nep->nt*j+k];
      coeffs[ctx->nmat-2] = ctx->coeffD[nep->nt*(ctx->nmat-1)+k];
      ierr = BVMultVec(W,1.0,0,v,coeffs);CHKERRQ(ierr);
      ierr = MatMult(nep->A[k],v,t);CHKERRQ(ierr);
      ierr = VecAXPY(q,1.0,t);CHKERRQ(ierr);
    }
    ierr = KSPSolve(ctx->ksp[idxrktg],q,t);CHKERRQ(ierr);
    ierr = VecScale(t,-1.0);CHKERRQ(ierr);
  } else {
    for (k=0;k<deg-1;k++) {
      ierr = BVGetColumn(W,k,&w);CHKERRQ(ierr);
      ierr = MatMult(ctx->D[k],w,q);CHKERRQ(ierr);
      ierr = BVRestoreColumn(W,k,&w);CHKERRQ(ierr);
      ierr = BVInsertVec(W,k,q);CHKERRQ(ierr);
    }
    ierr = BVGetColumn(W,deg-1,&w);CHKERRQ(ierr);
    ierr = MatMult(ctx->D[deg],w,q);CHKERRQ(ierr);
    ierr = BVRestoreColumn(W,k,&w);CHKERRQ(ierr);
    ierr = BVInsertVec(W,k,q);CHKERRQ(ierr);
    for (j=0;j<ctx->nmat-1;j++) coeffs[j] = 1.0;
    ierr = BVMultVec(W,1.0,0.0,q,coeffs);CHKERRQ(ierr);
    ierr = KSPSolve(ctx->ksp[idxrktg],q,t);CHKERRQ(ierr);
    ierr = VecScale(t,-1.0);CHKERRQ(ierr);
  }
  ierr = PetscFree(coeffs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Compute TOAR coefficients of the blocks of the new Arnoldi vector computed
*/
static PetscErrorCode NEPTOARCoefficients(NEP nep,PetscScalar sigma,PetscInt nv,PetscScalar *S,PetscInt ls,PetscScalar *r,PetscInt lr,PetscScalar *x,PetscScalar *work)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       k,j,d=ctx->nmat-1;
  PetscScalar    *t=work;

  PetscFunctionBegin;
  ierr = NEPNLEIGSEvalNRTFunct(nep,d-1,sigma,t);CHKERRQ(ierr);
  for (k=0;k<d-1;k++) {
    for (j=0;j<=nv;j++) r[k*lr+j] += t[k]*x[j];
  }
  for (j=0;j<=nv;j++) r[(d-1)*lr+j] = t[d-1]*x[j];
  PetscFunctionReturn(0);
}

/*
  Compute continuation vector coefficients for the Rational-Krylov run.
  dim(work) >= (end-ini)*(end-ini+1) + end+1 + 2*(end-ini+1), dim(t) = end.
*/
static PetscErrorCode NEPNLEIGS_RKcontinuation(NEP nep,PetscInt ini,PetscInt end,PetscScalar *K,PetscScalar *H,PetscInt ld,PetscScalar sigma,PetscScalar *S,PetscInt lds,PetscScalar *cont,PetscScalar *t,PetscScalar *work)
{
  PetscErrorCode ierr;
  PetscScalar    *x,*W,*tau,sone=1.0,szero=0.0;
  PetscInt       i,j,n1,n,nwu=0;
  PetscBLASInt   info,n_,n1_,one=1,dim,lds_;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (!ctx->nshifts || !end) {
    t[0] = 1;
    ierr = PetscArraycpy(cont,S+end*lds,lds);CHKERRQ(ierr);
  } else {
    n   = end-ini;
    n1  = n+1;
    x   = work+nwu;
    nwu += end+1;
    tau = work+nwu;
    nwu += n;
    W   = work+nwu;
    nwu += n1*n;
    for (j=ini;j<end;j++) {
      for (i=ini;i<=end;i++) W[(j-ini)*n1+i-ini] = K[j*ld+i] -H[j*ld+i]*sigma;
    }
    ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(n1,&n1_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(end+1,&dim);CHKERRQ(ierr);
    PetscStackCallBLAS("LAPACKgeqrf",LAPACKgeqrf_(&n1_,&n_,W,&n1_,tau,work+nwu,&n1_,&info));
    SlepcCheckLapackInfo("geqrf",info);
    for (i=0;i<end;i++) t[i] = 0.0;
    t[end] = 1.0;
    for (j=n-1;j>=0;j--) {
      for (i=0;i<ini+j;i++) x[i] = 0.0;
      x[ini+j] = 1.0;
      for (i=j+1;i<n1;i++) x[i+ini] = W[i+n1*j];
      tau[j] = PetscConj(tau[j]);
      PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&dim,&one,x,&one,tau+j,t,&dim,work+nwu));
    }
    ierr = PetscBLASIntCast(lds,&lds_);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&lds_,&n1_,&sone,S,&lds_,t,&one,&szero,cont,&one));
  }
  PetscFunctionReturn(0);
}

/*
  Compute a run of Arnoldi iterations
*/
PetscErrorCode NEPNLEIGSTOARrun(NEP nep,PetscScalar *K,PetscScalar *H,PetscInt ldh,BV W,PetscInt k,PetscInt *M,PetscBool *breakdown,Vec *t_)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;
  PetscInt       i,j,m=*M,lwa,deg=ctx->nmat-1,lds,nqt,ld,l;
  Vec            t;
  PetscReal      norm;
  PetscScalar    *x,*work,*tt,sigma,*cont,*S;
  PetscBool      lindep;
  Mat            MS;

  PetscFunctionBegin;
  ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
  ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
  ierr = BVGetSizes(nep->V,NULL,NULL,&ld);CHKERRQ(ierr);
  lds = ld*deg;
  ierr = BVGetActiveColumns(nep->V,&l,&nqt);CHKERRQ(ierr);
  lwa = PetscMax(ld,deg)+(m+1)*(m+1)+4*(m+1);
  ierr = PetscMalloc4(ld,&x,lwa,&work,m+1,&tt,lds,&cont);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(ctx->V,0,m);CHKERRQ(ierr);
  for (j=k;j<m;j++) {
    sigma = ctx->shifts[(++(ctx->idxrk))%ctx->nshiftsw];

    /* Continuation vector */
    ierr = NEPNLEIGS_RKcontinuation(nep,0,j,K,H,ldh,sigma,S,lds,cont,tt,work);CHKERRQ(ierr);

    /* apply operator */
    ierr = BVGetColumn(nep->V,nqt,&t);CHKERRQ(ierr);
    ierr = NEPTOARExtendBasis(nep,(ctx->idxrk)%ctx->nshiftsw,cont,ld,nqt,W,nep->V,t,S+(j+1)*lds,ld,t_);CHKERRQ(ierr);
    ierr = BVRestoreColumn(nep->V,nqt,&t);CHKERRQ(ierr);

    /* orthogonalize */
    ierr = BVOrthogonalizeColumn(nep->V,nqt,x,&norm,&lindep);CHKERRQ(ierr);
    if (!lindep) {
      x[nqt] = norm;
      ierr = BVScaleColumn(nep->V,nqt,1.0/norm);CHKERRQ(ierr);
      nqt++;
    } else x[nqt] = 0.0;

    ierr = NEPTOARCoefficients(nep,sigma,nqt-1,cont,ld,S+(j+1)*lds,ld,x,work);CHKERRQ(ierr);

    /* Level-2 orthogonalization */
    ierr = BVOrthogonalizeColumn(ctx->V,j+1,H+j*ldh,&norm,breakdown);CHKERRQ(ierr);
    H[j+1+ldh*j] = norm;
    if (ctx->nshifts) {
      for (i=0;i<=j;i++) K[i+ldh*j] = sigma*H[i+ldh*j] + tt[i];
      K[j+1+ldh*j] = sigma*H[j+1+ldh*j];
    }
    if (*breakdown) {
      *M = j+1;
      break;
    }
    ierr = BVScaleColumn(ctx->V,j+1,1.0/norm);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(nep->V,l,nqt);CHKERRQ(ierr);
  }
  ierr = PetscFree4(x,work,tt,cont);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
  ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSolve_NLEIGS(NEP nep)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;
  PetscInt       i,k=0,l,nv=0,ld,lds,ldds,nq,newn;
  PetscInt       deg=ctx->nmat-1,nconv=0;
  PetscScalar    *S,*Q,*H,*pU,*K,betak=0,*eigr,*eigi;
  PetscReal      betah;
  PetscBool      falselock=PETSC_FALSE,breakdown=PETSC_FALSE;
  BV             W;
  Mat            MS,MQ,U;

  PetscFunctionBegin;
  if (ctx->lock) {
    /* undocumented option to use a cheaper locking instead of the true locking */
    ierr = PetscOptionsGetBool(NULL,NULL,"-nep_nleigs_falselocking",&falselock,NULL);CHKERRQ(ierr);
  }

  ierr = BVGetSizes(nep->V,NULL,NULL,&ld);CHKERRQ(ierr);
  lds = deg*ld;
  ierr = DSGetLeadingDimension(nep->ds,&ldds);CHKERRQ(ierr);
  if (!ctx->nshifts) {
    ierr = PetscMalloc2(nep->ncv,&eigr,nep->ncv,&eigi);CHKERRQ(ierr);
  } else { eigr = nep->eigr; eigi = nep->eigi; }
  ierr = BVDuplicateResize(nep->V,PetscMax(nep->nt-1,ctx->nmat-1),&W);CHKERRQ(ierr);


  /* clean projected matrix (including the extra-arrow) */
  ierr = DSGetArray(nep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
  ierr = PetscArrayzero(H,ldds*ldds);CHKERRQ(ierr);
  ierr = DSRestoreArray(nep->ds,DS_MAT_A,&H);CHKERRQ(ierr);
  if (ctx->nshifts) {
    ierr = DSGetArray(nep->ds,DS_MAT_B,&H);CHKERRQ(ierr);
    ierr = PetscArrayzero(H,ldds*ldds);CHKERRQ(ierr);
    ierr = DSRestoreArray(nep->ds,DS_MAT_B,&H);CHKERRQ(ierr);
  }

  /* Get the starting Arnoldi vector */
  ierr = BVTensorBuildFirstColumn(ctx->V,nep->nini);CHKERRQ(ierr);

  /* Restart loop */
  l = 0;
  while (nep->reason == NEP_CONVERGED_ITERATING) {
    nep->its++;

    /* Compute an nv-step Krylov relation */
    nv = PetscMin(nep->nconv+nep->mpd,nep->ncv);
    if (ctx->nshifts) { ierr = DSGetArray(nep->ds,DS_MAT_A,&K);CHKERRQ(ierr); }
    ierr = DSGetArray(nep->ds,ctx->nshifts?DS_MAT_B:DS_MAT_A,&H);CHKERRQ(ierr);
    ierr = NEPNLEIGSTOARrun(nep,K,H,ldds,W,nep->nconv+l,&nv,&breakdown,nep->work);CHKERRQ(ierr);
    betah = PetscAbsScalar(H[(nv-1)*ldds+nv]);
    ierr = DSRestoreArray(nep->ds,ctx->nshifts?DS_MAT_B:DS_MAT_A,&H);CHKERRQ(ierr);
    if (ctx->nshifts) {
      betak = K[(nv-1)*ldds+nv];
      ierr = DSRestoreArray(nep->ds,DS_MAT_A,&K);CHKERRQ(ierr);
    }
    ierr = DSSetDimensions(nep->ds,nv,0,nep->nconv,nep->nconv+l);CHKERRQ(ierr);
    if (l==0) {
      ierr = DSSetState(nep->ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
    } else {
      ierr = DSSetState(nep->ds,DS_STATE_RAW);CHKERRQ(ierr);
    }

    /* Solve projected problem */
    ierr = DSSolve(nep->ds,nep->eigr,nep->eigi);CHKERRQ(ierr);
    ierr = DSSort(nep->ds,nep->eigr,nep->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
    if (!ctx->nshifts) {
      ierr = DSUpdateExtraRow(nep->ds);CHKERRQ(ierr);
    }
    ierr = DSSynchronize(nep->ds,nep->eigr,nep->eigi);CHKERRQ(ierr);

    /* Check convergence */
    ierr = NEPNLEIGSKrylovConvergence(nep,PETSC_FALSE,nep->nconv,nv-nep->nconv,betah,betak,&k,nep->work);CHKERRQ(ierr);
    ierr = (*nep->stopping)(nep,nep->its,nep->max_it,k,nep->nev,&nep->reason,nep->stoppingctx);CHKERRQ(ierr);

    /* Update l */
    if (nep->reason != NEP_CONVERGED_ITERATING || breakdown) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)*ctx->keep));
      if (!breakdown) {
        /* Prepare the Rayleigh quotient for restart */
        if (!ctx->nshifts) {
          ierr = DSTruncate(nep->ds,k+l);CHKERRQ(ierr);
          ierr = DSGetDimensions(nep->ds,&newn,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
          l = newn-k;
        } else {
          ierr = DSGetArray(nep->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
          ierr = DSGetArray(nep->ds,DS_MAT_B,&H);CHKERRQ(ierr);
          ierr = DSGetArray(nep->ds,DS_MAT_A,&K);CHKERRQ(ierr);
          for (i=ctx->lock?k:0;i<k+l;i++) {
            H[k+l+i*ldds] = betah*Q[nv-1+i*ldds];
            K[k+l+i*ldds] = betak*Q[nv-1+i*ldds];
          }
          ierr = DSRestoreArray(nep->ds,DS_MAT_B,&H);CHKERRQ(ierr);
          ierr = DSRestoreArray(nep->ds,DS_MAT_A,&K);CHKERRQ(ierr);
          ierr = DSRestoreArray(nep->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
        }
      }
    }
    nconv = k;
    if (!ctx->lock && nep->reason == NEP_CONVERGED_ITERATING && !breakdown) { l += k; k = 0; }

    /* Update S */
    ierr = DSGetMat(nep->ds,ctx->nshifts?DS_MAT_Z:DS_MAT_Q,&MQ);CHKERRQ(ierr);
    ierr = BVMultInPlace(ctx->V,MQ,nep->nconv,k+l);CHKERRQ(ierr);
    ierr = MatDestroy(&MQ);CHKERRQ(ierr);

    /* Copy last column of S */
    ierr = BVCopyColumn(ctx->V,nv,k+l);CHKERRQ(ierr);

    if (breakdown && nep->reason == NEP_CONVERGED_ITERATING) {
      /* Stop if breakdown */
      ierr = PetscInfo2(nep,"Breakdown (it=%D norm=%g)\n",nep->its,(double)betah);CHKERRQ(ierr);
      nep->reason = NEP_DIVERGED_BREAKDOWN;
    }
    if (nep->reason != NEP_CONVERGED_ITERATING) l--;
    /* truncate S */
    ierr = BVGetActiveColumns(nep->V,NULL,&nq);CHKERRQ(ierr);
    if (k+l+deg<=nq) {
      ierr = BVSetActiveColumns(ctx->V,nep->nconv,k+l+1);CHKERRQ(ierr);
      if (!falselock && ctx->lock) {
        ierr = BVTensorCompress(ctx->V,k-nep->nconv);CHKERRQ(ierr);
      } else {
        ierr = BVTensorCompress(ctx->V,0);CHKERRQ(ierr);
      }
    }
    nep->nconv = k;
    if (!ctx->nshifts) {
      for (i=0;i<nv;i++) { eigr[i] = nep->eigr[i]; eigi[i] = nep->eigi[i]; }
      ierr = NEPNLEIGSBackTransform((PetscObject)nep,nv,eigr,eigi);CHKERRQ(ierr);
    }
    ierr = NEPMonitor(nep,nep->its,nconv,eigr,eigi,nep->errest,nv);CHKERRQ(ierr);
  }
  nep->nconv = nconv;
  if (nep->nconv>0) {
    ierr = BVSetActiveColumns(ctx->V,0,nep->nconv);CHKERRQ(ierr);
    ierr = BVGetActiveColumns(nep->V,NULL,&nq);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(nep->V,0,nq);CHKERRQ(ierr);
    if (nq>nep->nconv) {
      ierr = BVTensorCompress(ctx->V,nep->nconv);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(nep->V,0,nep->nconv);CHKERRQ(ierr);
      nq = nep->nconv;
    }
    if (ctx->nshifts) {
      ierr = DSGetMat(nep->ds,DS_MAT_B,&MQ);CHKERRQ(ierr);
      ierr = BVMultInPlace(ctx->V,MQ,0,nep->nconv);CHKERRQ(ierr);
      ierr = MatDestroy(&MQ);CHKERRQ(ierr);
    }
    ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
    ierr = PetscMalloc1(nq*nep->nconv,&pU);CHKERRQ(ierr);
    for (i=0;i<nep->nconv;i++) {
      ierr = PetscArraycpy(pU+i*nq,S+i*lds,nq);CHKERRQ(ierr);
    }
    ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
    ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,nq,nep->nconv,pU,&U);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(nep->V,0,nq);CHKERRQ(ierr);
    ierr = BVMultInPlace(nep->V,U,0,nep->nconv);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(nep->V,0,nep->nconv);CHKERRQ(ierr);
    ierr = MatDestroy(&U);CHKERRQ(ierr);
    ierr = PetscFree(pU);CHKERRQ(ierr);
    /* truncate Schur decomposition and change the state to raw so that
       DSVectors() computes eigenvectors from scratch */
    ierr = DSSetDimensions(nep->ds,nep->nconv,0,0,0);CHKERRQ(ierr);
    ierr = DSSetState(nep->ds,DS_STATE_RAW);CHKERRQ(ierr);
  }

  /* Map eigenvalues back to the original problem */
  if (!ctx->nshifts) {
    ierr = NEPNLEIGSBackTransform((PetscObject)nep,nep->nconv,nep->eigr,nep->eigi);CHKERRQ(ierr);
    ierr = PetscFree2(eigr,eigi);CHKERRQ(ierr);
  }
  ierr = BVDestroy(&W);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSSetSingularitiesFunction_NLEIGS(NEP nep,PetscErrorCode (*fun)(NEP,PetscInt*,PetscScalar*,void*),void *ctx)
{
  NEP_NLEIGS *nepctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (fun) nepctx->computesingularities = fun;
  if (ctx) nepctx->singularitiesctx     = ctx;
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@C
   NEPNLEIGSSetSingularitiesFunction - Sets a user function to compute a discretization
   of the singularity set (where T(.) is not analytic).

   Logically Collective on nep

   Input Parameters:
+  nep - the NEP context
.  fun - user function (if NULL then NEP retains any previously set value)
-  ctx - [optional] user-defined context for private data for the function
         (may be NULL, in which case NEP retains any previously set value)

   Calling Sequence of fun:
$   fun(NEP nep,PetscInt *maxnp,PetscScalar *xi,void *ctx)

+   nep   - the NEP context
.   maxnp - on input number of requested points in the discretization (can be set)
.   xi    - computed values of the discretization
-   ctx   - optional context, as set by NEPNLEIGSSetSingularitiesFunction()

   Notes:
   The user-defined function can set a smaller value of maxnp if necessary.
   It is wrong to return a larger value.

   If the problem type has been set to rational with NEPSetProblemType(),
   then it is not necessary to set the singularities explicitly since the
   solver will try to determine them automatically.

   Level: intermediate

.seealso: NEPNLEIGSGetSingularitiesFunction(), NEPSetProblemType()
@*/
PetscErrorCode NEPNLEIGSSetSingularitiesFunction(NEP nep,PetscErrorCode (*fun)(NEP,PetscInt*,PetscScalar*,void*),void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  ierr = PetscTryMethod(nep,"NEPNLEIGSSetSingularitiesFunction_C",(NEP,PetscErrorCode(*)(NEP,PetscInt*,PetscScalar*,void*),void*),(nep,fun,ctx));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSGetSingularitiesFunction_NLEIGS(NEP nep,PetscErrorCode (**fun)(NEP,PetscInt*,PetscScalar*,void*),void **ctx)
{
  NEP_NLEIGS *nepctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (fun) *fun = nepctx->computesingularities;
  if (ctx) *ctx = nepctx->singularitiesctx;
  PetscFunctionReturn(0);
}

/*@C
   NEPNLEIGSGetSingularitiesFunction - Returns the Function and optionally the user
   provided context for computing a discretization of the singularity set.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameters:
+  fun - location to put the function (or NULL)
-  ctx - location to stash the function context (or NULL)

   Level: advanced

.seealso: NEPNLEIGSSetSingularitiesFunction()
@*/
PetscErrorCode NEPNLEIGSGetSingularitiesFunction(NEP nep,PetscErrorCode (**fun)(NEP,PetscInt*,PetscScalar*,void*),void **ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  ierr = PetscUseMethod(nep,"NEPNLEIGSGetSingularitiesFunction_C",(NEP,PetscErrorCode(**)(NEP,PetscInt*,PetscScalar*,void*),void**),(nep,fun,ctx));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSSetRestart_NLEIGS(NEP nep,PetscReal keep)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (keep==PETSC_DEFAULT) ctx->keep = 0.5;
  else {
    if (keep<0.1 || keep>0.9) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"The keep argument must be in the range [0.1,0.9]");
    ctx->keep = keep;
  }
  PetscFunctionReturn(0);
}

/*@
   NEPNLEIGSSetRestart - Sets the restart parameter for the NLEIGS
   method, in particular the proportion of basis vectors that must be kept
   after restart.

   Logically Collective on nep

   Input Parameters:
+  nep  - the nonlinear eigensolver context
-  keep - the number of vectors to be kept at restart

   Options Database Key:
.  -nep_nleigs_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,0.9]. The default is 0.5.

   Level: advanced

.seealso: NEPNLEIGSGetRestart()
@*/
PetscErrorCode NEPNLEIGSSetRestart(NEP nep,PetscReal keep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,keep,2);
  ierr = PetscTryMethod(nep,"NEPNLEIGSSetRestart_C",(NEP,PetscReal),(nep,keep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSGetRestart_NLEIGS(NEP nep,PetscReal *keep)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  *keep = ctx->keep;
  PetscFunctionReturn(0);
}

/*@
   NEPNLEIGSGetRestart - Gets the restart parameter used in the NLEIGS method.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  keep - the restart parameter

   Level: advanced

.seealso: NEPNLEIGSSetRestart()
@*/
PetscErrorCode NEPNLEIGSGetRestart(NEP nep,PetscReal *keep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidRealPointer(keep,2);
  ierr = PetscUseMethod(nep,"NEPNLEIGSGetRestart_C",(NEP,PetscReal*),(nep,keep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSSetLocking_NLEIGS(NEP nep,PetscBool lock)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(0);
}

/*@
   NEPNLEIGSSetLocking - Choose between locking and non-locking variants of
   the NLEIGS method.

   Logically Collective on nep

   Input Parameters:
+  nep  - the nonlinear eigensolver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -nep_nleigs_locking - Sets the locking flag

   Notes:
   The default is to lock converged eigenpairs when the method restarts.
   This behaviour can be changed so that all directions are kept in the
   working subspace even if already converged to working accuracy (the
   non-locking variant).

   Level: advanced

.seealso: NEPNLEIGSGetLocking()
@*/
PetscErrorCode NEPNLEIGSSetLocking(NEP nep,PetscBool lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,lock,2);
  ierr = PetscTryMethod(nep,"NEPNLEIGSSetLocking_C",(NEP,PetscBool),(nep,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSGetLocking_NLEIGS(NEP nep,PetscBool *lock)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(0);
}

/*@
   NEPNLEIGSGetLocking - Gets the locking flag used in the NLEIGS method.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: NEPNLEIGSSetLocking()
@*/
PetscErrorCode NEPNLEIGSGetLocking(NEP nep,PetscBool *lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidBoolPointer(lock,2);
  ierr = PetscUseMethod(nep,"NEPNLEIGSGetLocking_C",(NEP,PetscBool*),(nep,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSSetInterpolation_NLEIGS(NEP nep,PetscReal tol,PetscInt degree)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (tol == PETSC_DEFAULT) {
    ctx->ddtol = PETSC_DEFAULT;
    nep->state = NEP_STATE_INITIAL;
  } else {
    if (tol <= 0.0) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    ctx->ddtol = tol;
  }
  if (degree == PETSC_DEFAULT || degree == PETSC_DECIDE) {
    ctx->ddmaxit = 0;
    if (nep->state) { ierr = NEPReset(nep);CHKERRQ(ierr); }
    nep->state = NEP_STATE_INITIAL;
  } else {
    if (degree <= 0) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of degree. Must be > 0");
    if (ctx->ddmaxit != degree) {
      ctx->ddmaxit = degree;
      if (nep->state) { ierr = NEPReset(nep);CHKERRQ(ierr); }
      nep->state = NEP_STATE_INITIAL;
    }
  }
  PetscFunctionReturn(0);
}

/*@
   NEPNLEIGSSetInterpolation - Sets the tolerance and maximum degree
   when building the interpolation via divided differences.

   Logically Collective on nep

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  tol    - tolerance to stop computing divided differences
-  degree - maximum degree of interpolation

   Options Database Key:
+  -nep_nleigs_interpolation_tol <tol> - Sets the tolerance to stop computing divided differences
-  -nep_nleigs_interpolation_degree <degree> - Sets the maximum degree of interpolation

   Notes:
   Use PETSC_DEFAULT for either argument to assign a reasonably good value.

   Level: advanced

.seealso: NEPNLEIGSGetInterpolation()
@*/
PetscErrorCode NEPNLEIGSSetInterpolation(NEP nep,PetscReal tol,PetscInt degree)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,tol,2);
  PetscValidLogicalCollectiveInt(nep,degree,3);
  ierr = PetscTryMethod(nep,"NEPNLEIGSSetInterpolation_C",(NEP,PetscReal,PetscInt),(nep,tol,degree));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSGetInterpolation_NLEIGS(NEP nep,PetscReal *tol,PetscInt *degree)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (tol)    *tol    = ctx->ddtol;
  if (degree) *degree = ctx->ddmaxit;
  PetscFunctionReturn(0);
}

/*@
   NEPNLEIGSGetInterpolation - Gets the tolerance and maximum degree
   when building the interpolation via divided differences.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
+  tol    - tolerance to stop computing divided differences
-  degree - maximum degree of interpolation

   Level: advanced

.seealso: NEPNLEIGSSetInterpolation()
@*/
PetscErrorCode NEPNLEIGSGetInterpolation(NEP nep,PetscReal *tol,PetscInt *degree)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  ierr = PetscTryMethod(nep,"NEPNLEIGSGetInterpolation_C",(NEP,PetscReal*,PetscInt*),(nep,tol,degree));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSSetRKShifts_NLEIGS(NEP nep,PetscInt ns,PetscScalar *shifts)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       i;

  PetscFunctionBegin;
  if (ns<0) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_WRONG,"Number of shifts must be non-negative");
  if (ctx->nshifts) { ierr = PetscFree(ctx->shifts);CHKERRQ(ierr); }
  for (i=0;i<ctx->nshiftsw;i++) { ierr = KSPDestroy(&ctx->ksp[i]);CHKERRQ(ierr); }
  ierr = PetscFree(ctx->ksp);CHKERRQ(ierr);
  ctx->ksp = NULL;
  if (ns) {
    ierr = PetscMalloc1(ns,&ctx->shifts);CHKERRQ(ierr);
    for (i=0;i<ns;i++) ctx->shifts[i] = shifts[i];
  }
  ctx->nshifts = ns;
  nep->state   = NEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   NEPNLEIGSSetRKShifts - Sets a list of shifts to be used in the Rational
   Krylov method.

   Logically Collective on nep

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  ns     - number of shifts
-  shifts - array of scalar values specifying the shifts

   Options Database Key:
.  -nep_nleigs_rk_shifts - Sets the list of shifts

   Notes:
   If only one shift is provided, the built subspace built is equivalent to
   shift-and-invert Krylov-Schur (provided that the absolute convergence
   criterion is used).

   In the case of real scalars, complex shifts are not allowed. In the
   command line, a comma-separated list of complex values can be provided with
   the format [+/-][realnumber][+/-]realnumberi with no spaces, e.g.
   -nep_nleigs_rk_shifts 1.0+2.0i,1.5+2.0i,1.0+1.5i

   Use ns=0 to remove previously set shifts.

   Level: advanced

.seealso: NEPNLEIGSGetRKShifts()
@*/
PetscErrorCode NEPNLEIGSSetRKShifts(NEP nep,PetscInt ns,PetscScalar shifts[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,ns,2);
  if (ns) PetscValidScalarPointer(nep,shifts);
  ierr = PetscTryMethod(nep,"NEPNLEIGSSetRKShifts_C",(NEP,PetscInt,PetscScalar*),(nep,ns,shifts));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSGetRKShifts_NLEIGS(NEP nep,PetscInt *ns,PetscScalar **shifts)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscInt       i;

  PetscFunctionBegin;
  *ns = ctx->nshifts;
  if (ctx->nshifts) {
    ierr = PetscMalloc1(ctx->nshifts,shifts);CHKERRQ(ierr);
    for (i=0;i<ctx->nshifts;i++) (*shifts)[i] = ctx->shifts[i];
  }
  PetscFunctionReturn(0);
}

/*@C
   NEPNLEIGSGetRKShifts - Gets the list of shifts used in the Rational
   Krylov method.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
+  ns     - number of shifts
-  shifts - array of shifts

   Note:
   The user is responsible for deallocating the returned array.

   Level: advanced

.seealso: NEPNLEIGSSetRKShifts()
@*/
PetscErrorCode NEPNLEIGSGetRKShifts(NEP nep,PetscInt *ns,PetscScalar *shifts[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidIntPointer(nep,ns);
  PetscValidPointer(nep,shifts);
  ierr = PetscTryMethod(nep,"NEPNLEIGSGetRKShifts_C",(NEP,PetscInt*,PetscScalar**),(nep,ns,shifts));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSGetKSPs_NLEIGS(NEP nep,PetscInt *nsolve,KSP **ksp)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;
  PetscInt       i;
  PC             pc;

  PetscFunctionBegin;
  if (!ctx->ksp) {
    ierr = NEPNLEIGSSetShifts(nep,&ctx->nshiftsw);CHKERRQ(ierr);
    ierr = PetscMalloc1(ctx->nshiftsw,&ctx->ksp);CHKERRQ(ierr);
    for (i=0;i<ctx->nshiftsw;i++) {
      ierr = KSPCreate(PetscObjectComm((PetscObject)nep),&ctx->ksp[i]);CHKERRQ(ierr);
      ierr = PetscObjectIncrementTabLevel((PetscObject)ctx->ksp[i],(PetscObject)nep,1);CHKERRQ(ierr);
      ierr = KSPSetOptionsPrefix(ctx->ksp[i],((PetscObject)nep)->prefix);CHKERRQ(ierr);
      ierr = KSPAppendOptionsPrefix(ctx->ksp[i],"nep_nleigs_");CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->ksp[i]);CHKERRQ(ierr);
      ierr = PetscObjectSetOptions((PetscObject)ctx->ksp[i],((PetscObject)nep)->options);CHKERRQ(ierr);
      ierr = KSPSetErrorIfNotConverged(ctx->ksp[i],PETSC_TRUE);CHKERRQ(ierr);
      ierr = KSPSetTolerances(ctx->ksp[i],SLEPC_DEFAULT_TOL,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
      ierr = KSPGetPC(ctx->ksp[i],&pc);CHKERRQ(ierr);
      ierr = KSPSetType(ctx->ksp[i],KSPPREONLY);CHKERRQ(ierr);
      ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
    }
  }
  if (nsolve) *nsolve = ctx->nshiftsw;
  if (ksp)    *ksp    = ctx->ksp;
  PetscFunctionReturn(0);
}

/*@C
   NEPNLEIGSGetKSPs - Retrieve the array of linear solver objects associated with
   the nonlinear eigenvalue solver.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameters:
+  nsolve - number of returned KSP objects
-  ksp - array of linear solver object

   Notes:
   The number of KSP objects is equal to the number of shifts provided by the user,
   or 1 if the user did not provide shifts.

   Level: advanced

.seealso: NEPNLEIGSSetRKShifts()
@*/
PetscErrorCode NEPNLEIGSGetKSPs(NEP nep,PetscInt *nsolve,KSP **ksp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  ierr = PetscUseMethod(nep,"NEPNLEIGSGetKSPs_C",(NEP,PetscInt*,KSP**),(nep,nsolve,ksp));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSSetFullBasis_NLEIGS(NEP nep,PetscBool fullbasis)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (fullbasis!=ctx->fullbasis) {
    ctx->fullbasis = fullbasis;
    nep->state     = NEP_STATE_INITIAL;
    nep->useds     = PetscNot(fullbasis);
  }
  PetscFunctionReturn(0);
}

/*@
   NEPNLEIGSSetFullBasis - Choose between TOAR-basis (default) and full-basis
   variants of the NLEIGS method.

   Logically Collective on nep

   Input Parameters:
+  nep  - the nonlinear eigensolver context
-  fullbasis - true if the full-basis variant must be selected

   Options Database Key:
.  -nep_nleigs_full_basis - Sets the full-basis flag

   Notes:
   The default is to use a compact representation of the Krylov basis, that is,
   V = (I otimes U) S, with a tensor BV. This behaviour can be changed so that
   the full basis V is explicitly stored and operated with. This variant is more
   expensive in terms of memory and computation, but is necessary in some cases,
   particularly for two-sided computations, see NEPSetTwoSided().

   In the full-basis variant, the NLEIGS solver uses an EPS object to explicitly
   solve the linearized eigenproblem, see NEPNLEIGSGetEPS().

   Level: advanced

.seealso: NEPNLEIGSGetFullBasis(), NEPNLEIGSGetEPS(), NEPSetTwoSided(), BVCreateTensor()
@*/
PetscErrorCode NEPNLEIGSSetFullBasis(NEP nep,PetscBool fullbasis)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(nep,fullbasis,2);
  ierr = PetscTryMethod(nep,"NEPNLEIGSSetFullBasis_C",(NEP,PetscBool),(nep,fullbasis));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPNLEIGSGetFullBasis_NLEIGS(NEP nep,PetscBool *fullbasis)
{
  NEP_NLEIGS *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  *fullbasis = ctx->fullbasis;
  PetscFunctionReturn(0);
}

/*@
   NEPNLEIGSGetFullBasis - Gets the flag that indicates if NLEIGS is using the
   full-basis variant.

   Not Collective

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Output Parameter:
.  fullbasis - the flag

   Level: advanced

.seealso: NEPNLEIGSSetFullBasis()
@*/
PetscErrorCode NEPNLEIGSGetFullBasis(NEP nep,PetscBool *fullbasis)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidBoolPointer(fullbasis,2);
  ierr = PetscUseMethod(nep,"NEPNLEIGSGetFullBasis_C",(NEP,PetscBool*),(nep,fullbasis));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define SHIFTMAX 30

PetscErrorCode NEPSetFromOptions_NLEIGS(PetscOptionItems *PetscOptionsObject,NEP nep)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;
  PetscInt       i=0,k;
  PetscBool      flg1,flg2,b;
  PetscReal      r;
  PetscScalar    array[SHIFTMAX];

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"NEP NLEIGS Options");CHKERRQ(ierr);

    ierr = PetscOptionsReal("-nep_nleigs_restart","Proportion of vectors kept after restart","NEPNLEIGSSetRestart",0.5,&r,&flg1);CHKERRQ(ierr);
    if (flg1) { ierr = NEPNLEIGSSetRestart(nep,r);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-nep_nleigs_locking","Choose between locking and non-locking variants","NEPNLEIGSSetLocking",PETSC_FALSE,&b,&flg1);CHKERRQ(ierr);
    if (flg1) { ierr = NEPNLEIGSSetLocking(nep,b);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-nep_nleigs_full_basis","Choose between TOAR and full-basis variants","NEPNLEIGSSetFullBasis",PETSC_FALSE,&b,&flg1);CHKERRQ(ierr);
    if (flg1) { ierr = NEPNLEIGSSetFullBasis(nep,b);CHKERRQ(ierr); }

    ierr = NEPNLEIGSGetInterpolation(nep,&r,&i);CHKERRQ(ierr);
    if (!i) i = PETSC_DEFAULT;
    ierr = PetscOptionsInt("-nep_nleigs_interpolation_degree","Maximum number of terms for interpolation via divided differences","NEPNLEIGSSetInterpolation",i,&i,&flg1);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-nep_nleigs_interpolation_tol","Tolerance for interpolation via divided differences","NEPNLEIGSSetInterpolation",r,&r,&flg2);CHKERRQ(ierr);
    if (flg1 || flg2) { ierr = NEPNLEIGSSetInterpolation(nep,r,i);CHKERRQ(ierr); }

    k = SHIFTMAX;
    for (i=0;i<k;i++) array[i] = 0;
    ierr = PetscOptionsScalarArray("-nep_nleigs_rk_shifts","Shifts for Rational Krylov","NEPNLEIGSSetRKShifts",array,&k,&flg1);CHKERRQ(ierr);
    if (flg1) { ierr = NEPNLEIGSSetRKShifts(nep,k,array);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);

  if (!ctx->ksp) { ierr = NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp);CHKERRQ(ierr); }
  for (i=0;i<ctx->nshiftsw;i++) {
    ierr = KSPSetFromOptions(ctx->ksp[i]);CHKERRQ(ierr);
  }

  if (ctx->fullbasis) {
    if (!ctx->eps) { ierr = NEPNLEIGSGetEPS(nep,&ctx->eps);CHKERRQ(ierr); }
    ierr = EPSSetFromOptions(ctx->eps);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPView_NLEIGS(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;
  PetscBool      isascii;
  PetscInt       i;
  char           str[50];

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  %d%% of basis vectors kept after restart\n",(int)(100*ctx->keep));CHKERRQ(ierr);
    if (ctx->fullbasis) {
      ierr = PetscViewerASCIIPrintf(viewer,"  using the full-basis variant\n");CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"  using the %slocking variant\n",ctx->lock?"":"non-");CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"  divided difference terms: used=%D, max=%D\n",ctx->nmat,ctx->ddmaxit);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  tolerance for divided difference convergence: %g\n",(double)ctx->ddtol);CHKERRQ(ierr);
    if (ctx->nshifts) {
      ierr = PetscViewerASCIIPrintf(viewer,"  RK shifts: ");CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
      for (i=0;i<ctx->nshifts;i++) {
        ierr = SlepcSNPrintfScalar(str,50,ctx->shifts[i],PETSC_FALSE);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"%s%s",str,(i<ctx->nshifts-1)?",":"");CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    }
    if (!ctx->ksp) { ierr = NEPNLEIGSGetKSPs(nep,&ctx->nshiftsw,&ctx->ksp);CHKERRQ(ierr); }
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = KSPView(ctx->ksp[0],viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    if (ctx->fullbasis) {
      if (!ctx->eps) { ierr = NEPNLEIGSGetEPS(nep,&ctx->eps);CHKERRQ(ierr); }
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = EPSView(ctx->eps,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPReset_NLEIGS(NEP nep)
{
  PetscErrorCode ierr;
  PetscInt       k;
  NEP_NLEIGS     *ctx=(NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (nep->fui==NEP_USER_INTERFACE_SPLIT) {
    ierr = PetscFree(ctx->coeffD);CHKERRQ(ierr);
  } else {
    for (k=0;k<ctx->nmat;k++) { ierr = MatDestroy(&ctx->D[k]);CHKERRQ(ierr); }
  }
  ierr = PetscFree4(ctx->s,ctx->xi,ctx->beta,ctx->D);CHKERRQ(ierr);
  for (k=0;k<ctx->nshiftsw;k++) { ierr = KSPReset(ctx->ksp[k]);CHKERRQ(ierr); }
  if (ctx->vrn) {
    ierr = VecDestroy(&ctx->vrn);CHKERRQ(ierr);
  }
  if (ctx->fullbasis) {
    ierr = MatDestroy(&ctx->A);CHKERRQ(ierr);
    ierr = EPSReset(ctx->eps);CHKERRQ(ierr);
    for (k=0;k<4;k++) { ierr = VecDestroy(&ctx->w[k]);CHKERRQ(ierr); }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDestroy_NLEIGS(NEP nep)
{
  PetscErrorCode ierr;
  PetscInt       k;
  NEP_NLEIGS     *ctx = (NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  ierr = BVDestroy(&ctx->V);CHKERRQ(ierr);
  for (k=0;k<ctx->nshiftsw;k++) { ierr = KSPDestroy(&ctx->ksp[k]);CHKERRQ(ierr); }
  ierr = PetscFree(ctx->ksp);CHKERRQ(ierr);
  if (ctx->nshifts) { ierr = PetscFree(ctx->shifts);CHKERRQ(ierr); }
  if (ctx->fullbasis) { ierr = EPSDestroy(&ctx->eps);CHKERRQ(ierr); }
  ierr = PetscFree(nep->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetSingularitiesFunction_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetSingularitiesFunction_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetLocking_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetLocking_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetInterpolation_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetInterpolation_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetRKShifts_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetRKShifts_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetKSPs_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetFullBasis_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetFullBasis_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetEPS_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetEPS_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_NLEIGS(NEP nep)
{
  PetscErrorCode ierr;
  NEP_NLEIGS     *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(nep,&ctx);CHKERRQ(ierr);
  nep->data  = (void*)ctx;
  ctx->lock  = PETSC_TRUE;
  ctx->ddtol = PETSC_DEFAULT;

  nep->useds = PETSC_TRUE;
  nep->hasts = PETSC_TRUE;

  nep->ops->setup          = NEPSetUp_NLEIGS;
  nep->ops->setfromoptions = NEPSetFromOptions_NLEIGS;
  nep->ops->view           = NEPView_NLEIGS;
  nep->ops->destroy        = NEPDestroy_NLEIGS;
  nep->ops->reset          = NEPReset_NLEIGS;

  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetSingularitiesFunction_C",NEPNLEIGSSetSingularitiesFunction_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetSingularitiesFunction_C",NEPNLEIGSGetSingularitiesFunction_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetRestart_C",NEPNLEIGSSetRestart_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetRestart_C",NEPNLEIGSGetRestart_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetLocking_C",NEPNLEIGSSetLocking_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetLocking_C",NEPNLEIGSGetLocking_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetInterpolation_C",NEPNLEIGSSetInterpolation_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetInterpolation_C",NEPNLEIGSGetInterpolation_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetRKShifts_C",NEPNLEIGSSetRKShifts_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetRKShifts_C",NEPNLEIGSGetRKShifts_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetKSPs_C",NEPNLEIGSGetKSPs_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetFullBasis_C",NEPNLEIGSSetFullBasis_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetFullBasis_C",NEPNLEIGSGetFullBasis_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSSetEPS_C",NEPNLEIGSSetEPS_NLEIGS);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPNLEIGSGetEPS_C",NEPNLEIGSGetEPS_NLEIGS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

