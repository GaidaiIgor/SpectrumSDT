/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc polynomial eigensolver: "stoar"

   Method: S-TOAR

   Algorithm:

       Symmetric Two-Level Orthogonal Arnoldi.

   References:

       [1] C. Campos and J.E. Roman, "Restarted Q-Arnoldi-type methods
           exploiting symmetry in quadratic eigenvalue problems", BIT
           Numer. Math. 56(4):1213-1236, 2016.
*/

#include <slepc/private/pepimpl.h>         /*I "slepcpep.h" I*/
#include "../src/pep/impls/krylov/pepkrylov.h"
#include <slepcblaslapack.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-stoar,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"Restarted {Q-Arnoldi-type} methods exploiting symmetry in quadratic eigenvalue problems\",\n"
  "   journal = \"{BIT} Numer. Math.\",\n"
  "   volume = \"56\",\n"
  "   number = \"4\",\n"
  "   pages = \"1213--1236\",\n"
  "   year = \"2016,\"\n"
  "   doi = \"https://doi.org/10.1007/s10543-016-0601-5\"\n"
  "}\n";

typedef struct {
  PetscReal   scal[2];
  Mat         A[2];
  Vec         t;
} PEP_STOAR_MATSHELL;

static PetscErrorCode MatMult_STOAR(Mat A,Vec x,Vec y)
{
  PetscErrorCode     ierr;
  PEP_STOAR_MATSHELL *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,&ctx);CHKERRQ(ierr);
  ierr = MatMult(ctx->A[0],x,y);CHKERRQ(ierr);
  ierr = VecScale(y,ctx->scal[0]);CHKERRQ(ierr);
  if (ctx->scal[1]) {
    ierr = MatMult(ctx->A[1],x,ctx->t);CHKERRQ(ierr);
    ierr = VecAXPY(y,ctx->scal[1],ctx->t);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDestroy_STOAR(Mat A)
{
  PEP_STOAR_MATSHELL *ctx;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->t);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSTOARSetUpInnerMatrix(PEP pep,Mat *B)
{
  Mat                pB[4],Bs[3],D[3];
  PetscInt           i,j,n,m;
  PEP_STOAR_MATSHELL *ctxMat[3];
  PEP_STOAR          *ctx=(PEP_STOAR*)pep->data;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  for (i=0;i<3;i++) {
    ierr = STGetMatrixTransformed(pep->st,i,&D[i]);CHKERRQ(ierr); /* D[2] = M */
  }
  ierr = MatGetLocalSize(D[2],&m,&n);CHKERRQ(ierr);

  for (j=0;j<3;j++) {
    ierr = PetscNew(ctxMat+j);CHKERRQ(ierr);
    ierr = MatCreateShell(PetscObjectComm((PetscObject)pep),m,n,PETSC_DETERMINE,PETSC_DETERMINE,ctxMat[j],&Bs[j]);CHKERRQ(ierr);
    ierr = MatShellSetOperation(Bs[j],MATOP_MULT,(void(*)(void))MatMult_STOAR);CHKERRQ(ierr);
    ierr = MatShellSetOperation(Bs[j],MATOP_DESTROY,(void(*)(void))MatDestroy_STOAR);CHKERRQ(ierr);
  }
  for (i=0;i<4;i++) pB[i] = NULL;
  if (ctx->alpha) {
    ctxMat[0]->A[0] = D[0]; ctxMat[0]->scal[0] = ctx->alpha; ctxMat[0]->scal[1] = 0.0;
    ctxMat[2]->A[0] = D[2]; ctxMat[2]->scal[0] = -ctx->alpha*pep->sfactor*pep->sfactor; ctxMat[2]->scal[1] = 0.0;
    pB[0] = Bs[0]; pB[3] = Bs[2];
  }
  if (ctx->beta) {
    i = (ctx->alpha)?1:0;
    ctxMat[0]->scal[1] = 0.0;
    ctxMat[0]->A[i] = D[1]; ctxMat[0]->scal[i] = -ctx->beta*pep->sfactor;
    ctxMat[1]->A[0] = D[2]; ctxMat[1]->scal[0] = -ctx->beta*pep->sfactor*pep->sfactor; ctxMat[1]->scal[1] = 0.0;
    pB[0] = Bs[0]; pB[1] = pB[2] = Bs[1];
  }
  ierr = BVCreateVec(pep->V,&ctxMat[0]->t);CHKERRQ(ierr);
  ierr = MatCreateNest(PetscObjectComm((PetscObject)pep),2,NULL,2,NULL,pB,B);CHKERRQ(ierr);
  for (j=0;j<3;j++) { ierr = MatDestroy(&Bs[j]);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSetUp_STOAR(PEP pep)
{
  PetscErrorCode    ierr;
  PetscBool         shift,sinv,flg;
  PEP_STOAR         *ctx = (PEP_STOAR*)pep->data;
  PetscInt          ld,i;
  PetscReal         eta;
  BVOrthogType      otype;
  BVOrthogBlockType obtype;

  PetscFunctionBegin;
  if (pep->problem_type!=PEP_HERMITIAN && pep->problem_type!=PEP_HYPERBOLIC) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Requested method is only available for Hermitian (or hyperbolic) problems");
  if (pep->nmat!=3) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Solver only available for quadratic problems");
  if (pep->basis!=PEP_BASIS_MONOMIAL) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Solver not implemented for non-monomial bases");
  /* spectrum slicing requires special treatment of default values */
  ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv);CHKERRQ(ierr);
  if (pep->which==PEP_ALL) {
    if (pep->stopping!=PEPStoppingBasic) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Spectrum slicing does not support user-defined stopping test");
    pep->ops->solve = PEPSolve_STOAR_QSlice;
    pep->ops->extractvectors = NULL;
    pep->ops->setdefaultst   = NULL;
    ierr = PEPSetUp_STOAR_QSlice(pep);CHKERRQ(ierr);
  } else {
    ierr = PEPSetDimensions_Default(pep,pep->nev,&pep->ncv,&pep->mpd);CHKERRQ(ierr);
    if (!ctx->lock && pep->mpd<pep->ncv) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Should not use mpd parameter in non-locking variant");
    if (pep->max_it==PETSC_DEFAULT) pep->max_it = PetscMax(100,2*(pep->nmat-1)*pep->n/pep->ncv);
    pep->ops->solve = PEPSolve_STOAR;
    ld   = pep->ncv+2;
    ierr = DSSetType(pep->ds,DSGHIEP);CHKERRQ(ierr);
    ierr = DSSetCompact(pep->ds,PETSC_TRUE);CHKERRQ(ierr);
    ierr = DSAllocate(pep->ds,ld);CHKERRQ(ierr);
    ierr = PEPBasisCoefficients(pep,pep->pbc);CHKERRQ(ierr);
    ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
    if (!flg) {
      ierr = PetscFree(pep->solvematcoeffs);CHKERRQ(ierr);
      ierr = PetscMalloc1(pep->nmat,&pep->solvematcoeffs);CHKERRQ(ierr);
      ierr = PetscLogObjectMemory((PetscObject)pep,pep->nmat*sizeof(PetscScalar));CHKERRQ(ierr);
      if (sinv) {
        ierr = PEPEvaluateBasis(pep,pep->target,0,pep->solvematcoeffs,NULL);CHKERRQ(ierr);
      } else {
        for (i=0;i<pep->nmat-1;i++) pep->solvematcoeffs[i] = 0.0;
        pep->solvematcoeffs[pep->nmat-1] = 1.0;
      }
    }
  }

  pep->lineariz = PETSC_TRUE;
  ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSHIFT,&shift);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv);CHKERRQ(ierr);
  if (!shift && !sinv) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Only STSHIFT and STSINVERT spectral transformations can be used");
  if (!pep->which) {
    if (sinv) pep->which = PEP_TARGET_MAGNITUDE;
    else pep->which = PEP_LARGEST_MAGNITUDE;
  }

  ierr = PEPAllocateSolution(pep,2);CHKERRQ(ierr);
  ierr = PEPSetWorkVecs(pep,4);CHKERRQ(ierr);
  ierr = BVDestroy(&ctx->V);CHKERRQ(ierr);
  ierr = BVCreateTensor(pep->V,pep->nmat-1,&ctx->V);CHKERRQ(ierr);
  ierr = BVGetOrthogonalization(pep->V,&otype,NULL,&eta,&obtype);CHKERRQ(ierr);
  ierr = BVSetOrthogonalization(ctx->V,otype,BV_ORTHOG_REFINE_ALWAYS,eta,obtype);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Compute a run of Lanczos iterations. dim(work)=(ctx->ld)*4
*/
static PetscErrorCode PEPSTOARrun(PEP pep,PetscReal *a,PetscReal *b,PetscReal *omega,PetscInt k,PetscInt *M,PetscBool *breakdown,PetscBool *symmlost,Vec *t_)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscInt       i,j,m=*M,l,lock;
  PetscInt       lds,d,ld,offq,nqt,ldds;
  Vec            v=t_[0],t=t_[1],q=t_[2];
  PetscReal      norm,sym=0.0,fro=0.0,*f;
  PetscScalar    *y,*S,*x,sigma;
  PetscBLASInt   j_,one=1;
  PetscBool      lindep,flg,sinvert=PETSC_FALSE;
  Mat            MS;

  PetscFunctionBegin;
  ierr = PetscMalloc1(*M,&y);CHKERRQ(ierr);
  ierr = BVGetSizes(pep->V,NULL,NULL,&ld);CHKERRQ(ierr);
  ierr = BVTensorGetDegree(ctx->V,&d);CHKERRQ(ierr);
  ierr = BVGetActiveColumns(pep->V,&lock,&nqt);CHKERRQ(ierr);
  lds = d*ld;
  offq = ld;
  ierr = DSGetLeadingDimension(pep->ds,&ldds);CHKERRQ(ierr);
  *breakdown = PETSC_FALSE; /* ----- */
  ierr = DSGetDimensions(pep->ds,NULL,NULL,&l,NULL,NULL);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(ctx->V,0,m);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(pep->V,0,nqt);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (!flg) {
    /* spectral transformation handled by the solver */
    ierr = PetscObjectTypeCompareAny((PetscObject)pep->st,&flg,STSINVERT,STSHIFT,"");CHKERRQ(ierr);
    if (!flg) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"ST type not supported for TOAR without transforming matrices");
    ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinvert);CHKERRQ(ierr);
    ierr = STGetShift(pep->st,&sigma);CHKERRQ(ierr);
  }
  for (j=k;j<m;j++) {
    /* apply operator */
    ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
    ierr = BVGetColumn(pep->V,nqt,&t);CHKERRQ(ierr);
    ierr = BVMultVec(pep->V,1.0,0.0,v,S+j*lds);CHKERRQ(ierr);
    if (!sinvert) {
      ierr = STMatMult(pep->st,0,v,q);CHKERRQ(ierr);
      ierr = BVMultVec(pep->V,1.0,0.0,v,S+offq+j*lds);CHKERRQ(ierr);
      ierr = STMatMult(pep->st,1,v,t);CHKERRQ(ierr);
      ierr = VecAXPY(q,pep->sfactor,t);CHKERRQ(ierr);
      if (ctx->beta && ctx->alpha) {
        ierr = STMatMult(pep->st,2,v,t);CHKERRQ(ierr);
        ierr = VecAXPY(q,-pep->sfactor*pep->sfactor*ctx->beta/ctx->alpha,t);CHKERRQ(ierr);
      }
      ierr = STMatSolve(pep->st,q,t);CHKERRQ(ierr);
      ierr = VecScale(t,-1.0/(pep->sfactor*pep->sfactor));CHKERRQ(ierr);
    } else {
      ierr = STMatMult(pep->st,1,v,q);CHKERRQ(ierr);
      ierr = STMatMult(pep->st,2,v,t);CHKERRQ(ierr);
      ierr = VecAXPY(q,sigma*pep->sfactor,t);CHKERRQ(ierr);
      ierr = VecScale(q,pep->sfactor);CHKERRQ(ierr);
      ierr = BVMultVec(pep->V,1.0,0.0,v,S+offq+j*lds);CHKERRQ(ierr);
      ierr = STMatMult(pep->st,2,v,t);CHKERRQ(ierr);
      ierr = VecAXPY(q,pep->sfactor*pep->sfactor,t);CHKERRQ(ierr);
      ierr = STMatSolve(pep->st,q,t);CHKERRQ(ierr);
      ierr = VecScale(t,-1.0);CHKERRQ(ierr);
    }
    ierr = BVRestoreColumn(pep->V,nqt,&t);CHKERRQ(ierr);

    /* orthogonalize */
    if (!sinvert) x = S+offq+(j+1)*lds;
    else x = S+(j+1)*lds;
    ierr = BVOrthogonalizeColumn(pep->V,nqt,x,&norm,&lindep);CHKERRQ(ierr);

    if (!lindep) {
      if (!sinvert) *(S+offq+(j+1)*lds+nqt) = norm;
      else *(S+(j+1)*lds+nqt) = norm;
      ierr = BVScaleColumn(pep->V,nqt,1.0/norm);CHKERRQ(ierr);
      nqt++;
    }
    if (!sinvert) {
      for (i=0;i<=nqt-1;i++) *(S+(j+1)*lds+i) = *(S+offq+j*lds+i);
      if (ctx->beta && ctx->alpha) {
        for (i=0;i<=nqt-1;i++) *(S+(j+1)*lds+offq+i) -= *(S+(j+1)*lds+i)*ctx->beta/ctx->alpha;
      }
    } else for (i=0;i<nqt;i++) *(S+(j+1)*lds+offq+i) = *(S+j*lds+i)+sigma*(*(S+(j+1)*lds+i));
    ierr = BVSetActiveColumns(pep->V,0,nqt);CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
    ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);

    /* level-2 orthogonalization */
    ierr = BVOrthogonalizeColumn(ctx->V,j+1,y,&norm,&lindep);CHKERRQ(ierr);
    a[j] = PetscRealPart(y[j]);
    omega[j+1] = (norm > 0)?1.0:-1.0;
    ierr = BVScaleColumn(ctx->V,j+1,1.0/norm);CHKERRQ(ierr);
    b[j] = PetscAbsReal(norm);

    /* check symmetry */
    ierr = DSGetArrayReal(pep->ds,DS_MAT_T,&f);CHKERRQ(ierr);
    if (j==k) {
      for (i=l;i<j-1;i++) y[i] = PetscAbsScalar(y[i])-PetscAbsReal(f[2*ldds+i]);
      for (i=0;i<l;i++) y[i] = 0.0;
    }
    ierr = DSRestoreArrayReal(pep->ds,DS_MAT_T,&f);CHKERRQ(ierr);
    if (j>0) y[j-1] = PetscAbsScalar(y[j-1])-PetscAbsReal(b[j-1]);
    ierr = PetscBLASIntCast(j,&j_);CHKERRQ(ierr);
    sym = SlepcAbs(BLASnrm2_(&j_,y,&one),sym);
    fro = SlepcAbs(fro,SlepcAbs(a[j],b[j]));
    if (j>0) fro = SlepcAbs(fro,b[j-1]);
    if (sym/fro>PetscMax(PETSC_SQRT_MACHINE_EPSILON,10*pep->tol)) {
      *symmlost = PETSC_TRUE;
      *M=j;
      break;
    }
  }
  ierr = BVSetActiveColumns(pep->V,lock,nqt);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(ctx->V,0,*M);CHKERRQ(ierr);
  ierr = PetscFree(y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if 0
static PetscErrorCode PEPSTOARpreKConvergence(PEP pep,PetscInt nv,PetscReal *norm,Vec *w)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscBLASInt   n_,one=1;
  PetscInt       lds=2*ctx->ld;
  PetscReal      t1,t2;
  PetscScalar    *S=ctx->S;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(nv+2,&n_);CHKERRQ(ierr);
  t1 = BLASnrm2_(&n_,S+nv*2*ctx->ld,&one);
  t2 = BLASnrm2_(&n_,S+(nv*2+1)*ctx->ld,&one);
  *norm = SlepcAbs(t1,t2);
  ierr = BVSetActiveColumns(pep->V,0,nv+2);CHKERRQ(ierr);
  ierr = BVMultVec(pep->V,1.0,0.0,w[1],S+nv*lds);CHKERRQ(ierr);
  ierr = STMatMult(pep->st,0,w[1],w[2]);CHKERRQ(ierr);
  ierr = VecNorm(w[2],NORM_2,&t1);CHKERRQ(ierr);
  ierr = BVMultVec(pep->V,1.0,0.0,w[1],S+ctx->ld+nv*lds);CHKERRQ(ierr);
  ierr = STMatMult(pep->st,2,w[1],w[2]);CHKERRQ(ierr);
  ierr = VecNorm(w[2],NORM_2,&t2);CHKERRQ(ierr);
  t2 *= pep->sfactor*pep->sfactor;
  *norm = PetscMax(*norm,SlepcAbs(t1,t2));
  PetscFunctionReturn(0);
}
#endif

PetscErrorCode PEPSolve_STOAR(PEP pep)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscInt       j,k,l,nv=0,ld,ldds,t,nq=0;
  PetscInt       nconv=0,deg=pep->nmat-1;
  PetscScalar    *Q,*om,sigma;
  PetscReal      beta,norm=1.0,*omega,*a,*b,*r;
  PetscBool      breakdown,symmlost=PETSC_FALSE,sinv=PETSC_FALSE,falselock=PETSC_TRUE,flg;
  Mat            MQ,A;
  Vec            vomega;

  PetscFunctionBegin;
  ierr = PetscCitationsRegister(citation,&cited);CHKERRQ(ierr);
  ierr = PEPSTOARSetUpInnerMatrix(pep,&A);CHKERRQ(ierr);
  ierr = BVSetMatrix(ctx->V,A,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  if (ctx->lock) {
    /* undocumented option to use a cheaper locking instead of the true locking */
    ierr = PetscOptionsGetBool(NULL,NULL,"-pep_stoar_falselocking",&falselock,NULL);CHKERRQ(ierr);
  }
  ierr = BVGetSizes(pep->V,NULL,NULL,&ld);CHKERRQ(ierr);
  ierr = STGetShift(pep->st,&sigma);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (pep->sfactor!=1.0) {
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

  /* Get the starting Arnoldi vector */
  ierr = BVTensorBuildFirstColumn(ctx->V,pep->nini);CHKERRQ(ierr);
  ierr = DSGetArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF,1,&vomega);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(ctx->V,0,1);CHKERRQ(ierr);
  ierr = BVGetSignature(ctx->V,vomega);CHKERRQ(ierr);
  ierr = VecGetArray(vomega,&om);CHKERRQ(ierr);
  omega[0] = PetscRealPart(om[0]);
  ierr = VecRestoreArray(vomega,&om);CHKERRQ(ierr);
  ierr = DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
  ierr = VecDestroy(&vomega);CHKERRQ(ierr);

  /* Restart loop */
  l = 0;
  ierr = DSGetLeadingDimension(pep->ds,&ldds);CHKERRQ(ierr);
  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;
    ierr = DSGetArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
    b = a+ldds;
    ierr = DSGetArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);

    /* Compute an nv-step Lanczos factorization */
    nv = PetscMin(pep->nconv+pep->mpd,pep->ncv);
    ierr = PEPSTOARrun(pep,a,b,omega,pep->nconv+l,&nv,&breakdown,&symmlost,pep->work);CHKERRQ(ierr);
    beta = b[nv-1];
    if (symmlost && nv==pep->nconv+l) {
      pep->reason = PEP_DIVERGED_SYMMETRY_LOST;
      pep->nconv = nconv;
      if (falselock || !ctx->lock) {
       ierr = BVSetActiveColumns(ctx->V,0,pep->nconv);CHKERRQ(ierr);
       ierr = BVTensorCompress(ctx->V,0);CHKERRQ(ierr);
      }
      break;
    }
    ierr = DSRestoreArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
    ierr = DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
    ierr = DSSetDimensions(pep->ds,nv,0,pep->nconv,pep->nconv+l);CHKERRQ(ierr);
    if (l==0) {
      ierr = DSSetState(pep->ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
    } else {
      ierr = DSSetState(pep->ds,DS_STATE_RAW);CHKERRQ(ierr);
    }

    /* Solve projected problem */
    ierr = DSSolve(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);
    ierr = DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DSSynchronize(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);

    /* Check convergence */
    /* ierr = PEPSTOARpreKConvergence(pep,nv,&norm,pep->work);CHKERRQ(ierr);*/
    norm = 1.0;
    ierr = DSGetDimensions(pep->ds,NULL,NULL,NULL,NULL,&t);CHKERRQ(ierr);
    ierr = PEPKrylovConvergence(pep,PETSC_FALSE,pep->nconv,t-pep->nconv,PetscAbsReal(beta)*norm,&k);CHKERRQ(ierr);
    ierr = (*pep->stopping)(pep,pep->its,pep->max_it,k,pep->nev,&pep->reason,pep->stoppingctx);CHKERRQ(ierr);

    /* Update l */
    if (pep->reason != PEP_CONVERGED_ITERATING || breakdown) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)/2));
      l = PetscMin(l,t);
      if (!breakdown) {
        ierr = DSGetArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
        if (*(a+ldds+k+l-1)!=0) {
          if (k+l<nv-1) l = l+1;
          else l = l-1;
        }
        /* Prepare the Rayleigh quotient for restart */
        ierr = DSGetArray(pep->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
        ierr = DSGetArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
        r = a + 2*ldds;
        for (j=k;j<k+l;j++) {
          r[j] = PetscRealPart(Q[nv-1+j*ldds]*beta);
        }
        b = a+ldds;
        b[k+l-1] = r[k+l-1];
        omega[k+l] = omega[nv];
        ierr = DSRestoreArray(pep->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
        ierr = DSRestoreArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
        ierr = DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
      }
    }
    nconv = k;
    if (!ctx->lock && pep->reason == PEP_CONVERGED_ITERATING && !breakdown) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */

    /* Update S */
    ierr = DSGetMat(pep->ds,DS_MAT_Q,&MQ);CHKERRQ(ierr);
    ierr = BVMultInPlace(ctx->V,MQ,pep->nconv,k+l);CHKERRQ(ierr);
    ierr = MatDestroy(&MQ);CHKERRQ(ierr);

    /* Copy last column of S */
    ierr = BVCopyColumn(ctx->V,nv,k+l);CHKERRQ(ierr);
    ierr = DSGetArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,k+l,&vomega);CHKERRQ(ierr);
    ierr = VecGetArray(vomega,&om);CHKERRQ(ierr);
    for (j=0;j<k+l;j++) om[j] = omega[j];
    ierr = VecRestoreArray(vomega,&om);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(ctx->V,0,k+l);CHKERRQ(ierr);
    ierr = BVSetSignature(ctx->V,vomega);CHKERRQ(ierr);
    ierr = VecDestroy(&vomega);CHKERRQ(ierr);
    ierr = DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);

    if (breakdown && pep->reason == PEP_CONVERGED_ITERATING) {
      /* stop if breakdown */
      ierr = PetscInfo2(pep,"Breakdown TOAR method (it=%D norm=%g)\n",pep->its,(double)beta);CHKERRQ(ierr);
      pep->reason = PEP_DIVERGED_BREAKDOWN;
    }
    if (pep->reason != PEP_CONVERGED_ITERATING) l--;
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
    ierr = BVSetActiveColumns(ctx->V,0,pep->nconv);CHKERRQ(ierr);
    ierr = BVGetActiveColumns(pep->V,NULL,&nq);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(pep->V,0,nq);CHKERRQ(ierr);
    if (nq>pep->nconv) {
      ierr = BVTensorCompress(ctx->V,pep->nconv);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(pep->V,0,pep->nconv);CHKERRQ(ierr);
    }
  }
  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (!flg && pep->ops->backtransform) {
      ierr = (*pep->ops->backtransform)(pep);CHKERRQ(ierr);
  }
  if (pep->sfactor!=1.0) {
    for (j=0;j<pep->nconv;j++) {
      pep->eigr[j] *= pep->sfactor;
      pep->eigi[j] *= pep->sfactor;
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

  /* truncate Schur decomposition and change the state to raw so that
     DSVectors() computes eigenvectors from scratch */
  ierr = DSSetDimensions(pep->ds,pep->nconv,0,0,0);CHKERRQ(ierr);
  ierr = DSSetState(pep->ds,DS_STATE_RAW);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSetFromOptions_STOAR(PetscOptionItems *PetscOptionsObject,PEP pep)
{
  PetscErrorCode ierr;
  PetscBool      flg,lock,b,f1,f2,f3;
  PetscInt       i,j,k;
  PetscReal      array[2]={0,0};
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"PEP STOAR Options");CHKERRQ(ierr);

    ierr = PetscOptionsBool("-pep_stoar_locking","Choose between locking and non-locking variants","PEPSTOARSetLocking",PETSC_FALSE,&lock,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSTOARSetLocking(pep,lock);CHKERRQ(ierr); }

    b = ctx->detect;
    ierr = PetscOptionsBool("-pep_stoar_detect_zeros","Check zeros during factorizations at interval boundaries","PEPSTOARSetDetectZeros",ctx->detect,&b,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSTOARSetDetectZeros(pep,b);CHKERRQ(ierr); }

    i = 1;
    j = k = PETSC_DECIDE;
    ierr = PetscOptionsInt("-pep_stoar_nev","Number of eigenvalues to compute in each subsolve (only for spectrum slicing)","PEPSTOARSetDimensions",20,&i,&f1);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-pep_stoar_ncv","Number of basis vectors in each subsolve (only for spectrum slicing)","PEPSTOARSetDimensions",40,&j,&f2);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-pep_stoar_mpd","Maximum dimension of projected problem in each subsolve (only for spectrum slicing)","PEPSTOARSetDimensions",40,&k,&f3);CHKERRQ(ierr);
    if (f1 || f2 || f3) { ierr = PEPSTOARSetDimensions(pep,i,j,k);CHKERRQ(ierr); }

    k = 2;
    ierr = PetscOptionsRealArray("-pep_stoar_linearization","Parameters of the linearization","PEPSTOARSetLinearization",array,&k,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = PEPSTOARSetLinearization(pep,array[0],array[1]);CHKERRQ(ierr);
    }

    b = ctx->checket;
    ierr = PetscOptionsBool("-pep_stoar_check_eigenvalue_type","Check eigenvalue type during spectrum slicing","PEPSTOARSetCheckEigenvalueType",ctx->checket,&b,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSTOARSetCheckEigenvalueType(pep,b);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARSetLocking_STOAR(PEP pep,PetscBool lock)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARSetLocking - Choose between locking and non-locking variants of
   the STOAR method.

   Logically Collective on pep

   Input Parameters:
+  pep  - the eigenproblem solver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -pep_stoar_locking - Sets the locking flag

   Notes:
   The default is to lock converged eigenpairs when the method restarts.
   This behaviour can be changed so that all directions are kept in the
   working subspace even if already converged to working accuracy (the
   non-locking variant).

   Level: advanced

.seealso: PEPSTOARGetLocking()
@*/
PetscErrorCode PEPSTOARSetLocking(PEP pep,PetscBool lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,lock,2);
  ierr = PetscTryMethod(pep,"PEPSTOARSetLocking_C",(PEP,PetscBool),(pep,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARGetLocking_STOAR(PEP pep,PetscBool *lock)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARGetLocking - Gets the locking flag used in the STOAR method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: PEPSTOARSetLocking()
@*/
PetscErrorCode PEPSTOARGetLocking(PEP pep,PetscBool *lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidBoolPointer(lock,2);
  ierr = PetscUseMethod(pep,"PEPSTOARGetLocking_C",(PEP,PetscBool*),(pep,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARGetInertias_STOAR(PEP pep,PetscInt *n,PetscReal **shifts,PetscInt **inertias)
{
  PetscErrorCode ierr;
  PetscInt       i,numsh;
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PEP_SR         sr = ctx->sr;

  PetscFunctionBegin;
  if (!pep->state) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONGSTATE,"Must call PEPSetUp() first");
  if (!ctx->sr) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see PEPSetInterval()");
  switch (pep->state) {
  case PEP_STATE_INITIAL:
    break;
  case PEP_STATE_SETUP:
    if (n) *n = 2;
    if (shifts) {
      ierr = PetscMalloc1(2,shifts);CHKERRQ(ierr);
      (*shifts)[0] = pep->inta;
      (*shifts)[1] = pep->intb;
    }
    if (inertias) {
      ierr = PetscMalloc1(2,inertias);CHKERRQ(ierr);
      (*inertias)[0] = (sr->dir==1)?sr->inertia0:sr->inertia1;
      (*inertias)[1] = (sr->dir==1)?sr->inertia1:sr->inertia0;
    }
    break;
  case PEP_STATE_SOLVED:
  case PEP_STATE_EIGENVECTORS:
    numsh = ctx->nshifts;
    if (n) *n = numsh;
    if (shifts) {
      ierr = PetscMalloc1(numsh,shifts);CHKERRQ(ierr);
      for (i=0;i<numsh;i++) (*shifts)[i] = ctx->shifts[i];
    }
    if (inertias) {
      ierr = PetscMalloc1(numsh,inertias);CHKERRQ(ierr);
      for (i=0;i<numsh;i++) (*inertias)[i] = ctx->inertias[i];
    }
    break;
  }
  PetscFunctionReturn(0);
}

/*@C
   PEPSTOARGetInertias - Gets the values of the shifts and their
   corresponding inertias in case of doing spectrum slicing for a
   computational interval.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameters:
+  n        - number of shifts, including the endpoints of the interval
.  shifts   - the values of the shifts used internally in the solver
-  inertias - the values of the inertia in each shift

   Notes:
   If called after PEPSolve(), all shifts used internally by the solver are
   returned (including both endpoints and any intermediate ones). If called
   before PEPSolve() and after PEPSetUp() then only the information of the
   endpoints of subintervals is available.

   This function is only available for spectrum slicing runs.

   The returned arrays should be freed by the user. Can pass NULL in any of
   the two arrays if not required.

   Fortran Notes:
   The calling sequence from Fortran is
.vb
   PEPSTOARGetInertias(pep,n,shifts,inertias,ierr)
   integer n
   double precision shifts(*)
   integer inertias(*)
.ve
   The arrays should be at least of length n. The value of n can be determined
   by an initial call
.vb
   PEPSTOARGetInertias(pep,n,PETSC_NULL_REAL,PETSC_NULL_INTEGER,ierr)
.ve

   Level: advanced

.seealso: PEPSetInterval()
@*/
PetscErrorCode PEPSTOARGetInertias(PEP pep,PetscInt *n,PetscReal **shifts,PetscInt **inertias)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidIntPointer(n,2);
  ierr = PetscUseMethod(pep,"PEPSTOARGetInertias_C",(PEP,PetscInt*,PetscReal**,PetscInt**),(pep,n,shifts,inertias));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARSetDetectZeros_STOAR(PEP pep,PetscBool detect)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  ctx->detect = detect;
  pep->state  = PEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARSetDetectZeros - Sets a flag to enforce detection of
   zeros during the factorizations throughout the spectrum slicing computation.

   Logically Collective on pep

   Input Parameters:
+  pep    - the eigenproblem solver context
-  detect - check for zeros

   Options Database Key:
.  -pep_stoar_detect_zeros - Check for zeros; this takes an optional
   bool value (0/1/no/yes/true/false)

   Notes:
   A zero in the factorization indicates that a shift coincides with an eigenvalue.

   This flag is turned off by default, and may be necessary in some cases.
   This feature currently requires an external package for factorizations
   with support for zero detection, e.g. MUMPS.

   Level: advanced

.seealso: PEPSetInterval()
@*/
PetscErrorCode PEPSTOARSetDetectZeros(PEP pep,PetscBool detect)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,detect,2);
  ierr = PetscTryMethod(pep,"PEPSTOARSetDetectZeros_C",(PEP,PetscBool),(pep,detect));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARGetDetectZeros_STOAR(PEP pep,PetscBool *detect)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  *detect = ctx->detect;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARGetDetectZeros - Gets the flag that enforces zero detection
   in spectrum slicing.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  detect - whether zeros detection is enforced during factorizations

   Level: advanced

.seealso: PEPSTOARSetDetectZeros()
@*/
PetscErrorCode PEPSTOARGetDetectZeros(PEP pep,PetscBool *detect)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidBoolPointer(detect,2);
  ierr = PetscUseMethod(pep,"PEPSTOARGetDetectZeros_C",(PEP,PetscBool*),(pep,detect));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARSetLinearization_STOAR(PEP pep,PetscReal alpha,PetscReal beta)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  if (beta==0.0 && alpha==0.0) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"Parameters alpha and beta cannot be zero simultaneously");
  ctx->alpha = alpha;
  ctx->beta  = beta;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARSetLinearization - Set the coefficients that define
   the linearization of a quadratic eigenproblem.

   Logically Collective on pep

   Input Parameters:
+  pep   - polynomial eigenvalue solver
.  alpha - first parameter of the linearization
-  beta  - second parameter of the linearization

   Options Database Key:
.  -pep_stoar_linearization <alpha,beta> - Sets the coefficients

   Notes:
   Cannot pass zero for both alpha and beta. The default values are
   alpha=1 and beta=0.

   Level: advanced

.seealso: PEPSTOARGetLinearization()
@*/
PetscErrorCode PEPSTOARSetLinearization(PEP pep,PetscReal alpha,PetscReal beta)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,alpha,2);
  PetscValidLogicalCollectiveReal(pep,beta,3);
  ierr = PetscTryMethod(pep,"PEPSTOARSetLinearization_C",(PEP,PetscReal,PetscReal),(pep,alpha,beta));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARGetLinearization_STOAR(PEP pep,PetscReal *alpha,PetscReal *beta)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  if (alpha) *alpha = ctx->alpha;
  if (beta)  *beta  = ctx->beta;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARGetLinearization - Returns the coefficients that define
   the linearization of a quadratic eigenproblem.

   Not Collective

   Input Parameter:
.  pep  - polynomial eigenvalue solver

   Output Parameters:
+  alpha - the first parameter of the linearization
-  beta  - the second parameter of the linearization

   Level: advanced

.seealso: PEPSTOARSetLinearization()
@*/
PetscErrorCode PEPSTOARGetLinearization(PEP pep,PetscReal *alpha,PetscReal *beta)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  ierr = PetscUseMethod(pep,"PEPSTOARGetLinearization_C",(PEP,PetscReal*,PetscReal*),(pep,alpha,beta));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARSetDimensions_STOAR(PEP pep,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  if (nev<1) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nev. Must be > 0");
  ctx->nev = nev;
  if (ncv == PETSC_DECIDE || ncv == PETSC_DEFAULT) {
    ctx->ncv = PETSC_DEFAULT;
  } else {
    if (ncv<1) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    ctx->ncv = ncv;
  }
  if (mpd == PETSC_DECIDE || mpd == PETSC_DEFAULT) {
    ctx->mpd = PETSC_DEFAULT;
  } else {
    if (mpd<1) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    ctx->mpd = mpd;
  }
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARSetDimensions - Sets the dimensions used for each subsolve
   step in case of doing spectrum slicing for a computational interval.
   The meaning of the parameters is the same as in PEPSetDimensions().

   Logically Collective on pep

   Input Parameters:
+  pep - the eigenproblem solver context
.  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the subsolve
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Key:
+  -eps_stoar_nev <nev> - Sets the number of eigenvalues
.  -eps_stoar_ncv <ncv> - Sets the dimension of the subspace
-  -eps_stoar_mpd <mpd> - Sets the maximum projected dimension

   Level: advanced

.seealso: PEPSTOARGetDimensions(), PEPSetDimensions(), PEPSetInterval()
@*/
PetscErrorCode PEPSTOARSetDimensions(PEP pep,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,nev,2);
  PetscValidLogicalCollectiveInt(pep,ncv,3);
  PetscValidLogicalCollectiveInt(pep,mpd,4);
  ierr = PetscTryMethod(pep,"PEPSTOARSetDimensions_C",(PEP,PetscInt,PetscInt,PetscInt),(pep,nev,ncv,mpd));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARGetDimensions_STOAR(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  if (nev) *nev = ctx->nev;
  if (ncv) *ncv = ctx->ncv;
  if (mpd) *mpd = ctx->mpd;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARGetDimensions - Gets the dimensions used for each subsolve
   step in case of doing spectrum slicing for a computational interval.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameters:
+  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the subsolve
-  mpd - the maximum dimension allowed for the projected problem

   Level: advanced

.seealso: PEPSTOARSetDimensions()
@*/
PetscErrorCode PEPSTOARGetDimensions(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  ierr = PetscUseMethod(pep,"PEPSTOARGetDimensions_C",(PEP,PetscInt*,PetscInt*,PetscInt*),(pep,nev,ncv,mpd));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARSetCheckEigenvalueType_STOAR(PEP pep,PetscBool checket)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  ctx->checket = checket;
  pep->state   = PEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARSetCheckEigenvalueType - Sets a flag to check that all the eigenvalues
   obtained throughout the spectrum slicing computation have the same definite type.

   Logically Collective on pep

   Input Parameters:
+  pep     - the eigenproblem solver context
-  checket - check eigenvalue type

   Options Database Key:
.  -pep_stoar_check_eigenvalue_type - Check eigenvalue type; this takes an optional
   bool value (0/1/no/yes/true/false)

   Notes:
   This option is relevant only for spectrum slicing computations, but it is
   ignored if the problem type is PEP_HYPERBOLIC.

   This flag is turned on by default, to guarantee that the computed eigenvalues
   have the same type (otherwise the computed solution might be wrong). But since
   the check is computationally quite expensive, the check may be turned off if
   the user knows for sure that all eigenvalues in the requested interval have
   the same type.

   Level: advanced

.seealso: PEPSetProblemType(), PEPSetInterval()
@*/
PetscErrorCode PEPSTOARSetCheckEigenvalueType(PEP pep,PetscBool checket)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,checket,2);
  ierr = PetscTryMethod(pep,"PEPSTOARSetCheckEigenvalueType_C",(PEP,PetscBool),(pep,checket));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOARGetCheckEigenvalueType_STOAR(PEP pep,PetscBool *checket)
{
  PEP_STOAR *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  *checket = ctx->checket;
  PetscFunctionReturn(0);
}

/*@
   PEPSTOARGetCheckEigenvalueType - Gets the flag for the eigenvalue type
   check in spectrum slicing.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  checket - whether eigenvalue type must be checked during spectrum slcing

   Level: advanced

.seealso: PEPSTOARSetCheckEigenvalueType()
@*/
PetscErrorCode PEPSTOARGetCheckEigenvalueType(PEP pep,PetscBool *checket)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidBoolPointer(checket,2);
  ierr = PetscUseMethod(pep,"PEPSTOARGetCheckEigenvalueType_C",(PEP,PetscBool*),(pep,checket));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPView_STOAR(PEP pep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  using the %slocking variant\n",ctx->lock?"":"non-");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  linearization parameters: alpha=%g beta=%g\n",(double)ctx->alpha,(double)ctx->beta);CHKERRQ(ierr);
    if (pep->which==PEP_ALL && !ctx->hyperbolic) {
      ierr = PetscViewerASCIIPrintf(viewer,"  checking eigenvalue type: %s\n",ctx->checket?"enabled":"disabled");CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPReset_STOAR(PEP pep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (pep->which==PEP_ALL) {
    ierr = PEPReset_STOAR_QSlice(pep);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPDestroy_STOAR(PEP pep)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  ierr = BVDestroy(&ctx->V);CHKERRQ(ierr);
  ierr = PetscFree(pep->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetLocking_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetLocking_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetDetectZeros_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetDetectZeros_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetInertias_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetDimensions_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetDimensions_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetLinearization_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetLinearization_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetCheckEigenvalueType_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetCheckEigenvalueType_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode PEPCreate_STOAR(PEP pep)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(pep,&ctx);CHKERRQ(ierr);
  pep->data = (void*)ctx;
  ctx->lock    = PETSC_TRUE;
  ctx->nev     = 1;
  ctx->ncv     = PETSC_DEFAULT;
  ctx->mpd     = PETSC_DEFAULT;
  ctx->alpha   = 1.0;
  ctx->beta    = 0.0;
  ctx->checket = PETSC_TRUE;

  pep->ops->setup          = PEPSetUp_STOAR;
  pep->ops->setfromoptions = PEPSetFromOptions_STOAR;
  pep->ops->destroy        = PEPDestroy_STOAR;
  pep->ops->view           = PEPView_STOAR;
  pep->ops->backtransform  = PEPBackTransform_Default;
  pep->ops->computevectors = PEPComputeVectors_Default;
  pep->ops->extractvectors = PEPExtractVectors_TOAR;
  pep->ops->reset          = PEPReset_STOAR;

  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetLocking_C",PEPSTOARSetLocking_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetLocking_C",PEPSTOARGetLocking_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetDetectZeros_C",PEPSTOARSetDetectZeros_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetDetectZeros_C",PEPSTOARGetDetectZeros_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetInertias_C",PEPSTOARGetInertias_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetDimensions_C",PEPSTOARGetDimensions_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetDimensions_C",PEPSTOARSetDimensions_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetLinearization_C",PEPSTOARSetLinearization_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetLinearization_C",PEPSTOARGetLinearization_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARSetCheckEigenvalueType_C",PEPSTOARSetCheckEigenvalueType_STOAR);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPSTOARGetCheckEigenvalueType_C",PEPSTOARGetCheckEigenvalueType_STOAR);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

