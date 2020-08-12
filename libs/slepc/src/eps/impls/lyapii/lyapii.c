/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "lyapii"

   Method: Lyapunov inverse iteration

   Algorithm:

       Lyapunov inverse iteration using LME solvers

   References:

       [1] H.C. Elman and M. Wu, "Lyapunov inverse iteration for computing a
           few rightmost eigenvalues of large generalized eigenvalue problems",
           SIAM J. Matrix Anal. Appl. 34(4):1685-1707, 2013.

       [2] K. Meerbergen and A. Spence, "Inverse iteration for purely imaginary
           eigenvalues with application to the detection of Hopf bifurcations in
           large-scale problems", SIAM J. Matrix Anal. Appl. 31:1982-1999, 2010.
*/

#include <slepc/private/epsimpl.h>          /*I "slepceps.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  LME      lme;      /* Lyapunov solver */
  DS       ds;       /* used to compute the SVD for compression */
  PetscInt rkl;      /* prescribed rank for the Lyapunov solver */
  PetscInt rkc;      /* the compressed rank, cannot be larger than rkl */
} EPS_LYAPII;

typedef struct {
  Mat      S;        /* the operator matrix, S=A^{-1}*B */
  BV       Q;        /* orthogonal basis of converged eigenvectors */
} EPS_LYAPII_MATSHELL;

typedef struct {
  Mat      S;        /* the matrix from which the implicit operator is built */
  PetscInt n;        /* the size of matrix S, the operator is nxn */
  LME      lme;      /* dummy LME object */
#if defined(PETSC_USE_COMPLEX)
  Mat      A,B,F;
  Vec      w;
#endif
} EPS_EIG_MATSHELL;

PetscErrorCode EPSSetUp_LyapII(EPS eps)
{
  PetscErrorCode ierr;
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;
  PetscBool      istrivial,issinv;

  PetscFunctionBegin;
  if (eps->ncv!=PETSC_DEFAULT) {
    if (eps->ncv<eps->nev+1) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv must be at least nev+1");
  } else eps->ncv = eps->nev+1;
  if (eps->mpd!=PETSC_DEFAULT) { ierr = PetscInfo(eps,"Warning: parameter mpd ignored\n");CHKERRQ(ierr); }
  if (!ctx->rkc) ctx->rkc = 10;
  if (!ctx->rkl) ctx->rkl = 3*ctx->rkc;
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(1000*eps->nev,100*eps->n);
  if (!eps->which) eps->which=EPS_LARGEST_REAL;
  if (eps->which!=EPS_LARGEST_REAL) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Wrong value of eps->which");
  if (eps->extraction) { ierr = PetscInfo(eps,"Warning: extraction type ignored\n");CHKERRQ(ierr); }
  if (eps->balance!=EPS_BALANCE_NONE) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Balancing not supported in this solver");
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");
  ierr = RGIsTrivial(eps->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support region filtering");

  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&issinv);CHKERRQ(ierr);
  if (!issinv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Must use STSINVERT spectral transformation");

  if (!ctx->lme) { ierr = EPSLyapIIGetLME(eps,&ctx->lme);CHKERRQ(ierr); }
  ierr = LMESetProblemType(ctx->lme,LME_LYAPUNOV);CHKERRQ(ierr);
  ierr = LMESetErrorIfNotConverged(ctx->lme,PETSC_TRUE);CHKERRQ(ierr);

  if (!ctx->ds) {
    ierr = DSCreate(PetscObjectComm((PetscObject)eps),&ctx->ds);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)ctx->ds);CHKERRQ(ierr);
    ierr = DSSetType(ctx->ds,DSSVD);CHKERRQ(ierr);
  }
  ierr = DSAllocate(ctx->ds,ctx->rkl);CHKERRQ(ierr);

  ierr = DSSetType(eps->ds,DSNHEP);CHKERRQ(ierr);
  ierr = DSAllocate(eps->ds,eps->ncv);CHKERRQ(ierr);

  ierr = EPSAllocateSolution(eps,0);CHKERRQ(ierr);
  ierr = EPSSetWorkVecs(eps,3);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_EPSLyapIIOperator(Mat M,Vec x,Vec r)
{
  PetscErrorCode      ierr;
  EPS_LYAPII_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  ierr = MatMult(matctx->S,x,r);CHKERRQ(ierr);
  ierr = BVOrthogonalizeVec(matctx->Q,r,NULL,NULL,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDestroy_EPSLyapIIOperator(Mat M)
{
  PetscErrorCode      ierr;
  EPS_LYAPII_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
  ierr = MatDestroy(&matctx->S);CHKERRQ(ierr);
  ierr = PetscFree(matctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_EigOperator(Mat M,Vec x,Vec y)
{
  PetscErrorCode    ierr;
  EPS_EIG_MATSHELL  *matctx;
#if !defined(PETSC_USE_COMPLEX)
  PetscInt          n;
  PetscScalar       *S,*Y,*C,zero=0.0,done=1.0,dtwo=2.0;
  const PetscScalar *X;
  PetscBLASInt      n_;
#endif

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
  ierr = MatMult(matctx->B,x,matctx->w);CHKERRQ(ierr);
  ierr = MatSolve(matctx->F,matctx->w,y);CHKERRQ(ierr);
#else
  ierr = VecGetArrayRead(x,&X);CHKERRQ(ierr);
  ierr = VecGetArray(y,&Y);CHKERRQ(ierr);
  ierr = MatDenseGetArray(matctx->S,&S);CHKERRQ(ierr);

  n = matctx->n;
  ierr = PetscCalloc1(n*n,&C);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);

  /* C = 2*S*X*S.' */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&dtwo,S,&n_,X,&n_,&zero,Y,&n_));
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","T",&n_,&n_,&n_,&done,Y,&n_,S,&n_,&zero,C,&n_));

  /* Solve S*Y + Y*S' = -C */
  ierr = LMEDenseLyapunov(matctx->lme,n,S,n,C,n,Y,n);CHKERRQ(ierr);

  ierr = PetscFree(C);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(x,&X);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&Y);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(matctx->S,&S);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDestroy_EigOperator(Mat M)
{
  PetscErrorCode   ierr;
  EPS_EIG_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,(void**)&matctx);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = MatDestroy(&matctx->A);CHKERRQ(ierr);
  ierr = MatDestroy(&matctx->B);CHKERRQ(ierr);
  ierr = MatDestroy(&matctx->F);CHKERRQ(ierr);
  ierr = VecDestroy(&matctx->w);CHKERRQ(ierr);
#endif
  ierr = PetscFree(matctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   EV2x2: solve the eigenproblem for a 2x2 matrix M
 */
static PetscErrorCode EV2x2(PetscScalar *M,PetscInt ld,PetscScalar *wr,PetscScalar *wi,PetscScalar *vec)
{
  PetscErrorCode ierr;
  PetscScalar    work[10];
  PetscBLASInt   lwork=10,ld_;
#if !defined(PETSC_HAVE_ESSL)
  PetscBLASInt   two=2,info;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      rwork[4];
#endif
#else
  PetscInt       i;
  PetscBLASInt   idummy,io=1;
  PetscScalar    wri[4];
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(ld,&ld_);CHKERRQ(ierr);
  ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
#if !defined(PETSC_HAVE_ESSL)
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_("N","V",&two,M,&ld_,wr,wi,NULL,&ld_,vec,&ld_,work,&lwork,&info));
#else
  PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_("N","V",&two,M,&ld_,wr,NULL,&ld_,vec,&ld_,work,&lwork,rwork,&info));
#endif
  SlepcCheckLapackInfo("geev",info);
#else /* defined(PETSC_HAVE_ESSL) */
  PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_(&io,M,&ld_,wri,vec,ld_,&idummy,&ld_,work,&lwork));
#if !defined(PETSC_USE_COMPLEX)
  for (i=0;i<2;i++) {
    wr[i] = wri[2*i];
    wi[i] = wri[2*i+1];
  }
#else
  for (i=0;i<2;i++) wr[i] = wri[i];
#endif
#endif
  ierr = PetscFPTrapPop();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   LyapIIBuildRHS: prepare the right-hand side of the Lyapunov equation SY + YS' = -2*S*Z*S'
   in factored form:
      if (V)  U=sqrt(2)*S*V    (uses 1 work vector)
      else    U=sqrt(2)*S*U    (uses 2 work vectors)
   where U,V are assumed to have rk columns.
 */
static PetscErrorCode LyapIIBuildRHS(Mat S,PetscInt rk,Mat U,BV V,Vec *work)
{
  PetscErrorCode ierr;
  PetscScalar    *array,*uu;
  PetscInt       i,nloc;
  Vec            v,u=work[0];

  PetscFunctionBegin;
  ierr = MatGetLocalSize(U,&nloc,NULL);CHKERRQ(ierr);
  for (i=0;i<rk;i++) {
    ierr = MatDenseGetColumn(U,i,&array);CHKERRQ(ierr);
    if (V) {
      ierr = BVGetColumn(V,i,&v);CHKERRQ(ierr);
    } else {
      v = work[1];
      ierr = VecPlaceArray(v,array);CHKERRQ(ierr);
    }
    ierr = MatMult(S,v,u);CHKERRQ(ierr);
    if (V) {
      ierr = BVRestoreColumn(V,i,&v);CHKERRQ(ierr);
    } else {
      ierr = VecResetArray(v);CHKERRQ(ierr);
    }
    ierr = VecScale(u,PetscSqrtReal(2.0));CHKERRQ(ierr);
    ierr = VecGetArray(u,&uu);CHKERRQ(ierr);
    ierr = PetscMemcpy(array,uu,nloc*sizeof(PetscScalar));CHKERRQ(ierr);
    ierr = VecRestoreArray(u,&uu);CHKERRQ(ierr);
    ierr = MatDenseRestoreColumn(U,&array);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   LyapIIBuildEigenMat: create shell matrix Op=A\B with A = kron(I,S)+kron(S,I), B = -2*kron(S,S)
   where S is a sequential square dense matrix of order n.
   v0 is the initial vector, should have the form v0 = w*w' (for instance 1*1')
 */
static PetscErrorCode LyapIIBuildEigenMat(LME lme,Mat S,Mat *Op,Vec *v0)
{
  PetscErrorCode   ierr;
  PetscInt         n,m;
  PetscBool        create=PETSC_FALSE;
  EPS_EIG_MATSHELL *matctx;
#if defined(PETSC_USE_COMPLEX)
  PetscScalar      theta,*aa,*bb,*ss;
  PetscInt         i,j,f,c,off,ld;
  IS               perm;
#endif

  PetscFunctionBegin;
  ierr = MatGetSize(S,&n,NULL);CHKERRQ(ierr);
  if (!*Op) create=PETSC_TRUE;
  else {
    ierr = MatGetSize(*Op,&m,NULL);CHKERRQ(ierr);
    if (m!=n*n) create=PETSC_TRUE;
  }
  if (create) {
    ierr = MatDestroy(Op);CHKERRQ(ierr);
    ierr = VecDestroy(v0);CHKERRQ(ierr);
    ierr = PetscNew(&matctx);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,n*n,n*n,NULL,&matctx->A);CHKERRQ(ierr);
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,n*n,n*n,NULL,&matctx->B);CHKERRQ(ierr);
    ierr = MatCreateVecs(matctx->A,NULL,&matctx->w);CHKERRQ(ierr);
#endif
    ierr = MatCreateShell(PETSC_COMM_SELF,n*n,n*n,PETSC_DETERMINE,PETSC_DETERMINE,matctx,Op);CHKERRQ(ierr);
    ierr = MatShellSetOperation(*Op,MATOP_MULT,(void(*)(void))MatMult_EigOperator);CHKERRQ(ierr);
    ierr = MatShellSetOperation(*Op,MATOP_DESTROY,(void(*)(void))MatDestroy_EigOperator);CHKERRQ(ierr);
    ierr = MatCreateVecs(*Op,NULL,v0);CHKERRQ(ierr);
  } else {
    ierr = MatShellGetContext(*Op,(void**)&matctx);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = MatZeroEntries(matctx->A);CHKERRQ(ierr);
#endif
  }
#if defined(PETSC_USE_COMPLEX)
  ierr = MatDenseGetArray(matctx->A,&aa);CHKERRQ(ierr);
  ierr = MatDenseGetArray(matctx->B,&bb);CHKERRQ(ierr);
  ierr = MatDenseGetArray(S,&ss);CHKERRQ(ierr);
  ld = n*n;
  for (f=0;f<n;f++) {
    off = f*n+f*n*ld;
    for (i=0;i<n;i++) for (j=0;j<n;j++) aa[off+i+j*ld] = ss[i+j*n];
    for (c=0;c<n;c++) {
      off = f*n+c*n*ld;
      theta = ss[f+c*n];
      for (i=0;i<n;i++) aa[off+i+i*ld] += theta;
      for (i=0;i<n;i++) for (j=0;j<n;j++) bb[off+i+j*ld] = -2*theta*ss[i+j*n];
    }
  }
  ierr = MatDenseRestoreArray(matctx->A,&aa);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(matctx->B,&bb);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(S,&ss);CHKERRQ(ierr);
  ierr = ISCreateStride(PETSC_COMM_SELF,n*n,0,1,&perm);CHKERRQ(ierr);
  ierr = MatDestroy(&matctx->F);CHKERRQ(ierr);
  ierr = MatDuplicate(matctx->A,MAT_COPY_VALUES,&matctx->F);CHKERRQ(ierr);
  ierr = MatLUFactor(matctx->F,perm,perm,0);CHKERRQ(ierr);
  ierr = ISDestroy(&perm);CHKERRQ(ierr);
#endif
  matctx->lme = lme;
  matctx->S = S;
  matctx->n = n;
  ierr = VecSet(*v0,1.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_LyapII(EPS eps)
{
  PetscErrorCode      ierr;
  EPS_LYAPII          *ctx = (EPS_LYAPII*)eps->data;
  PetscInt            i,ldds,rk,nloc,mloc,nv,idx,k;
  Vec                 v,w,z=eps->work[0],v0=NULL;
  Mat                 S,C,Ux[2],Y,Y1,R,U,W,X,Op=NULL;
  BV                  V;
  BVOrthogType        type;
  BVOrthogRefineType  refine;
  PetscScalar         eigr[2],eigi[2],*array,er,ei,*uu,*s,*xx,*aa,pM[4],vec[4];
  PetscReal           eta;
  EPS                 epsrr;
  PetscReal           norm;
  EPS_LYAPII_MATSHELL *matctx;

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(ctx->ds,&ldds);CHKERRQ(ierr);

  /* Operator for the Lyapunov equation */
  ierr = PetscNew(&matctx);CHKERRQ(ierr);
  ierr = STGetOperator(eps->st,&matctx->S);CHKERRQ(ierr);
  ierr = MatGetLocalSize(matctx->S,&mloc,&nloc);CHKERRQ(ierr);
  ierr = MatCreateShell(PetscObjectComm((PetscObject)eps),mloc,nloc,PETSC_DETERMINE,PETSC_DETERMINE,matctx,&S);CHKERRQ(ierr);
  matctx->Q = eps->V;
  ierr = MatShellSetOperation(S,MATOP_MULT,(void(*)(void))MatMult_EPSLyapIIOperator);CHKERRQ(ierr);
  ierr = MatShellSetOperation(S,MATOP_DESTROY,(void(*)(void))MatDestroy_EPSLyapIIOperator);CHKERRQ(ierr);
  ierr = LMESetCoefficients(ctx->lme,S,NULL,NULL,NULL);CHKERRQ(ierr);

  /* Right-hand side */
  ierr = BVDuplicateResize(eps->V,ctx->rkl,&V);CHKERRQ(ierr);
  ierr = BVGetOrthogonalization(V,&type,&refine,&eta,NULL);CHKERRQ(ierr);
  ierr = BVSetOrthogonalization(V,type,refine,eta,BV_ORTHOG_BLOCK_TSQR);CHKERRQ(ierr);
  ierr = MatCreateDense(PetscObjectComm((PetscObject)eps),eps->nloc,PETSC_DECIDE,PETSC_DECIDE,1,NULL,&Ux[0]);CHKERRQ(ierr);
  ierr = MatCreateDense(PetscObjectComm((PetscObject)eps),eps->nloc,PETSC_DECIDE,PETSC_DECIDE,2,NULL,&Ux[1]);CHKERRQ(ierr);
  nv = ctx->rkl;
  ierr = PetscMalloc1(nv,&s);CHKERRQ(ierr);

  /* Initialize first column */
  ierr = EPSGetStartVector(eps,0,NULL);CHKERRQ(ierr);
  ierr = BVGetColumn(eps->V,0,&v);CHKERRQ(ierr);
  ierr = BVInsertVec(V,0,v);CHKERRQ(ierr);
  ierr = BVRestoreColumn(eps->V,0,&v);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,0,0);CHKERRQ(ierr);  /* no deflation at the beginning */
  ierr = LyapIIBuildRHS(S,1,Ux[0],V,eps->work);CHKERRQ(ierr);
  idx = 0;

  /* EPS for rank reduction */
  ierr = EPSCreate(PETSC_COMM_SELF,&epsrr);CHKERRQ(ierr);
  ierr = EPSSetOptionsPrefix(epsrr,((PetscObject)eps)->prefix);CHKERRQ(ierr);
  ierr = EPSAppendOptionsPrefix(epsrr,"eps_lyapii_");CHKERRQ(ierr);
  ierr = EPSSetDimensions(epsrr,1,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = EPSSetTolerances(epsrr,PETSC_MACHINE_EPSILON*100,PETSC_DEFAULT);CHKERRQ(ierr);

  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Matrix for placing the solution of the Lyapunov equation (an alias of V) */
    ierr = BVSetActiveColumns(V,0,nv);CHKERRQ(ierr);
    ierr = BVGetMat(V,&Y1);CHKERRQ(ierr);
    ierr = MatZeroEntries(Y1);CHKERRQ(ierr);
    ierr = MatCreateLRC(NULL,Y1,NULL,NULL,&Y);CHKERRQ(ierr);
    ierr = LMESetSolution(ctx->lme,Y);CHKERRQ(ierr);

    /* Solve the Lyapunov equation SY + YS' = -2*S*Z*S' */
    ierr = MatCreateLRC(NULL,Ux[idx],NULL,NULL,&C);CHKERRQ(ierr);
    ierr = LMESetRHS(ctx->lme,C);CHKERRQ(ierr);
    ierr = MatDestroy(&C);CHKERRQ(ierr);
    ierr = LMESolve(ctx->lme);CHKERRQ(ierr);
    ierr = BVRestoreMat(V,&Y1);CHKERRQ(ierr);
    ierr = MatDestroy(&Y);CHKERRQ(ierr);

    /* SVD of the solution: [Q,R]=qr(V); [U,Sigma,~]=svd(R) */
    ierr = DSSetDimensions(ctx->ds,nv,nv,0,0);CHKERRQ(ierr);
    ierr = DSGetMat(ctx->ds,DS_MAT_A,&R);CHKERRQ(ierr);
    ierr = BVOrthogonalize(V,R);CHKERRQ(ierr);
    ierr = DSRestoreMat(ctx->ds,DS_MAT_A,&R);CHKERRQ(ierr);
    ierr = DSSetState(ctx->ds,DS_STATE_RAW);CHKERRQ(ierr);
    ierr = DSSolve(ctx->ds,s,NULL);CHKERRQ(ierr);

    /* Determine rank */
    rk = nv;
    for (i=1;i<nv;i++) if (PetscAbsScalar(s[i]/s[0])<PETSC_SQRT_MACHINE_EPSILON) {rk=i; break;}
    ierr = PetscInfo1(eps,"The computed solution of the Lyapunov equation has rank %D\n",rk);CHKERRQ(ierr);
    rk = PetscMin(rk,ctx->rkc);
    ierr = DSGetMat(ctx->ds,DS_MAT_U,&U);CHKERRQ(ierr);
    ierr = BVMultInPlace(V,U,0,rk);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(V,0,rk);CHKERRQ(ierr);
    ierr = MatDestroy(&U);CHKERRQ(ierr);

    /* Rank reduction */
    ierr = DSSetDimensions(ctx->ds,rk,rk,0,0);CHKERRQ(ierr);
    ierr = DSGetMat(ctx->ds,DS_MAT_A,&W);CHKERRQ(ierr);
    ierr = BVMatProject(V,S,V,W);CHKERRQ(ierr);
    ierr = LyapIIBuildEigenMat(ctx->lme,W,&Op,&v0);CHKERRQ(ierr); /* Op=A\B, A=kron(I,S)+kron(S,I), B=-2*kron(S,S) */
    ierr = EPSSetOperators(epsrr,Op,NULL);CHKERRQ(ierr);
    ierr = EPSSetInitialSpace(epsrr,1,&v0);CHKERRQ(ierr);
    ierr = EPSSolve(epsrr);CHKERRQ(ierr);
    ierr = MatDestroy(&W);CHKERRQ(ierr);
    ierr = EPSComputeVectors(epsrr);CHKERRQ(ierr);
    /* Copy first eigenvector, vec(A)=x */
    ierr = BVGetArray(epsrr->V,&xx);CHKERRQ(ierr);
    ierr = DSGetArray(ctx->ds,DS_MAT_A,&aa);CHKERRQ(ierr);
    for (i=0;i<rk;i++) {
      ierr = PetscMemcpy(aa+i*ldds,xx+i*rk,rk*sizeof(PetscScalar));CHKERRQ(ierr);
    }
    ierr = DSRestoreArray(ctx->ds,DS_MAT_A,&aa);CHKERRQ(ierr);
    ierr = BVRestoreArray(epsrr->V,&xx);CHKERRQ(ierr);
    ierr = DSSetState(ctx->ds,DS_STATE_RAW);CHKERRQ(ierr);
    /* Compute [U,Sigma,~] = svd(A), its rank should be 1 or 2 */
    ierr = DSSolve(ctx->ds,s,NULL);CHKERRQ(ierr);
    if (PetscAbsScalar(s[1]/s[0])<PETSC_SQRT_MACHINE_EPSILON) rk=1;
    else rk = 2;
    ierr = PetscInfo1(eps,"The eigenvector has rank %D\n",rk);CHKERRQ(ierr);
    ierr = DSGetMat(ctx->ds,DS_MAT_U,&U);CHKERRQ(ierr);
    ierr = BVMultInPlace(V,U,0,rk);CHKERRQ(ierr);
    ierr = MatDestroy(&U);CHKERRQ(ierr);

    /* Save V in Ux */
    idx = (rk==2)?1:0;
    for (i=0;i<rk;i++) {
      ierr = BVGetColumn(V,i,&v);CHKERRQ(ierr);
      ierr = VecGetArray(v,&uu);CHKERRQ(ierr);
      ierr = MatDenseGetColumn(Ux[idx],i,&array);CHKERRQ(ierr);
      ierr = PetscMemcpy(array,uu,eps->nloc*sizeof(PetscScalar));CHKERRQ(ierr);
      ierr = MatDenseRestoreColumn(Ux[idx],&array);CHKERRQ(ierr);
      ierr = VecRestoreArray(v,&uu);CHKERRQ(ierr);
      ierr = BVRestoreColumn(V,i,&v);CHKERRQ(ierr);
    }

    /* Eigenpair approximation */
    ierr = BVGetColumn(V,0,&v);CHKERRQ(ierr);
    ierr = MatMult(S,v,z);CHKERRQ(ierr);
    ierr = VecDot(z,v,pM);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,0,&v);CHKERRQ(ierr);
    if (rk>1) {
      ierr = BVGetColumn(V,1,&w);CHKERRQ(ierr);
      ierr = VecDot(z,w,pM+1);CHKERRQ(ierr);
      ierr = MatMult(S,w,z);CHKERRQ(ierr);
      ierr = VecDot(z,w,pM+3);CHKERRQ(ierr);
      ierr = BVGetColumn(V,0,&v);CHKERRQ(ierr);
      ierr = VecDot(z,v,pM+2);CHKERRQ(ierr);
      ierr = BVRestoreColumn(V,0,&v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(V,1,&w);CHKERRQ(ierr);
      ierr = EV2x2(pM,2,eigr,eigi,vec);CHKERRQ(ierr);
      ierr = MatCreateSeqDense(PETSC_COMM_SELF,2,2,vec,&X);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(V,0,rk);CHKERRQ(ierr);
      ierr = BVMultInPlace(V,X,0,rk);CHKERRQ(ierr);
      ierr = MatDestroy(&X);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
      norm = eigr[0]*eigr[0]+eigi[0]*eigi[0];
      er = eigr[0]/norm; ei = -eigi[0]/norm;
#else
      er =1.0/eigr[0]; ei = 0.0;
#endif
    } else {
      eigr[0] = pM[0]; eigi[0] = 0.0;
      er = 1.0/eigr[0]; ei = 0.0;
    }
    ierr = BVGetColumn(V,0,&v);CHKERRQ(ierr);
    if (eigi[0]!=0.0) {
      ierr = BVGetColumn(V,1,&w);CHKERRQ(ierr);
    } else w = NULL;
    eps->eigr[eps->nconv] = eigr[0]; eps->eigi[eps->nconv] = eigi[0];
    ierr = EPSComputeResidualNorm_Private(eps,PETSC_FALSE,er,ei,v,w,eps->work,&norm);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,0,&v);CHKERRQ(ierr);
    if (w) {
      ierr = BVRestoreColumn(V,1,&w);CHKERRQ(ierr);
    }
    ierr = (*eps->converged)(eps,er,ei,norm,&eps->errest[eps->nconv],eps->convergedctx);CHKERRQ(ierr);
    k = 0;
    if (eps->errest[eps->nconv]<eps->tol) {
      k++;
      if (rk==2) {
#if !defined (PETSC_USE_COMPLEX)
        eps->eigr[eps->nconv+k] = eigr[0]; eps->eigi[eps->nconv+k] = -eigi[0];
#else
        eps->eigr[eps->nconv+k] = PetscConj(eps->eigr[eps->nconv]);
#endif
        k++;
      }
      /* Store converged eigenpairs and vectors for deflation */
      for (i=0;i<k;i++) {
        ierr = BVGetColumn(V,i,&v);CHKERRQ(ierr);
        ierr = BVInsertVec(eps->V,eps->nconv+i,v);CHKERRQ(ierr);
        ierr = BVRestoreColumn(V,i,&v);CHKERRQ(ierr);
      }
      eps->nconv += k;
      ierr = BVSetActiveColumns(eps->V,eps->nconv-rk,eps->nconv);CHKERRQ(ierr);
      ierr = BVOrthogonalize(eps->V,NULL);CHKERRQ(ierr);
      ierr = DSSetDimensions(eps->ds,eps->nconv,0,0,0);CHKERRQ(ierr);
      ierr = DSGetMat(eps->ds,DS_MAT_A,&W);CHKERRQ(ierr);
      ierr = BVMatProject(eps->V,matctx->S,eps->V,W);CHKERRQ(ierr);
      ierr = DSRestoreMat(eps->ds,DS_MAT_A,&W);CHKERRQ(ierr);
      if (eps->nconv<eps->nev) {
        idx = 0;
        ierr = BVSetRandomColumn(V,0);CHKERRQ(ierr);
        ierr = BVNormColumn(V,0,NORM_2,&norm);CHKERRQ(ierr);
        ierr = BVScaleColumn(V,0,1.0/norm);CHKERRQ(ierr);
        ierr = LyapIIBuildRHS(S,1,Ux[idx],V,eps->work);CHKERRQ(ierr);
      }
    } else {
      /* Prepare right-hand side */
      ierr = LyapIIBuildRHS(S,rk,Ux[idx],NULL,eps->work);CHKERRQ(ierr);
    }
    ierr = (*eps->stopping)(eps,eps->its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);
    ierr = EPSMonitor(eps,eps->its,eps->nconv,eps->eigr,eps->eigi,eps->errest,eps->nconv+1);CHKERRQ(ierr);
  }
  ierr = STRestoreOperator(eps->st,&matctx->S);CHKERRQ(ierr);
  ierr = MatDestroy(&S);CHKERRQ(ierr);
  ierr = MatDestroy(&Ux[0]);CHKERRQ(ierr);
  ierr = MatDestroy(&Ux[1]);CHKERRQ(ierr);
  ierr = MatDestroy(&Op);CHKERRQ(ierr);
  ierr = VecDestroy(&v0);CHKERRQ(ierr);
  ierr = BVDestroy(&V);CHKERRQ(ierr);
  ierr = EPSDestroy(&epsrr);CHKERRQ(ierr);
  ierr = PetscFree(s);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_LyapII(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode ierr;
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;
  PetscInt       k,array[2]={PETSC_DEFAULT,PETSC_DEFAULT};
  PetscBool      flg;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS Lyapunov Inverse Iteration Options");CHKERRQ(ierr);

    k = 2;
    ierr = PetscOptionsIntArray("-eps_lyapii_ranks","Ranks for Lyapunov equation (one or two comma-separated integers)","EPSLyapIISetRanks",array,&k,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = EPSLyapIISetRanks(eps,array[0],array[1]);CHKERRQ(ierr);
    }

  ierr = PetscOptionsTail();CHKERRQ(ierr);

  if (!ctx->lme) { ierr = EPSLyapIIGetLME(eps,&ctx->lme);CHKERRQ(ierr); }
  ierr = LMESetFromOptions(ctx->lme);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLyapIISetRanks_LyapII(EPS eps,PetscInt rkc,PetscInt rkl)
{
  EPS_LYAPII *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  if (rkc==PETSC_DEFAULT) rkc = 10;
  if (rkc<2) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The compressed rank %D must be larger than 1",rkc);
  if (rkl==PETSC_DEFAULT) rkl = 3*rkc;
  if (rkl<rkc) SETERRQ2(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The Lyapunov rank %D cannot be smaller than the compressed rank %D",rkl,rkc);
  if (rkc != ctx->rkc) {
    ctx->rkc   = rkc;
    eps->state = EPS_STATE_INITIAL;
  }
  if (rkl != ctx->rkl) {
    ctx->rkl   = rkl;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSLyapIISetRanks - Set the ranks used in the solution of the Lyapunov equation.

   Collective on EPS

   Input Parameters:
+  eps - the eigenproblem solver context
.  rkc - the compressed rank
-  rkl - the Lyapunov rank

   Options Database Key:
.  -eps_lyapii_ranks <rkc,rkl> - Sets the rank parameters

   Notes:
   Lyapunov inverse iteration needs to solve a large-scale Lyapunov equation
   at each iteration of the eigensolver. For this, an iterative solver (LME)
   is used, which requires to prescribe the rank of the solution matrix X. This
   is the meaning of parameter rkl. Later, this matrix is compressed into
   another matrix of rank rkc. If not provided, rkl is a small multiple of rkc.

   Level: intermediate

.seealso: EPSLyapIIGetRanks()
@*/
PetscErrorCode EPSLyapIISetRanks(EPS eps,PetscInt rkc,PetscInt rkl)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,rkc,2);
  PetscValidLogicalCollectiveInt(eps,rkl,3);
  ierr = PetscTryMethod(eps,"EPSLyapIISetRanks_C",(EPS,PetscInt,PetscInt),(eps,rkc,rkl));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLyapIIGetRanks_LyapII(EPS eps,PetscInt *rkc,PetscInt *rkl)
{
  EPS_LYAPII *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  if (rkc) *rkc = ctx->rkc;
  if (rkl) *rkl = ctx->rkl;
  PetscFunctionReturn(0);
}

/*@
   EPSLyapIIGetRanks - Return the rank values used for the Lyapunov step.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
+  rkc - the compressed rank
-  rkl - the Lyapunov rank

   Level: intermediate

.seealso: EPSLyapIISetRanks()
@*/
PetscErrorCode EPSLyapIIGetRanks(EPS eps,PetscInt *rkc,PetscInt *rkl)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = PetscUseMethod(eps,"EPSLyapIIGetRanks_C",(EPS,PetscInt*,PetscInt*),(eps,rkc,rkl));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLyapIISetLME_LyapII(EPS eps,LME lme)
{
  PetscErrorCode ierr;
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)lme);CHKERRQ(ierr);
  ierr = LMEDestroy(&ctx->lme);CHKERRQ(ierr);
  ctx->lme = lme;
  ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)ctx->lme);CHKERRQ(ierr);
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSLyapIISetLME - Associate a linear matrix equation solver object (LME) to the
   eigenvalue solver.

   Collective on EPS

   Input Parameters:
+  eps - the eigenproblem solver context
-  lme - the linear matrix equation solver object

   Level: advanced

.seealso: EPSLyapIIGetLME()
@*/
PetscErrorCode EPSLyapIISetLME(EPS eps,LME lme)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(lme,LME_CLASSID,2);
  PetscCheckSameComm(eps,1,lme,2);
  ierr = PetscTryMethod(eps,"EPSLyapIISetLME_C",(EPS,LME),(eps,lme));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLyapIIGetLME_LyapII(EPS eps,LME *lme)
{
  PetscErrorCode ierr;
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  if (!ctx->lme) {
    ierr = LMECreate(PetscObjectComm((PetscObject)eps),&ctx->lme);CHKERRQ(ierr);
    ierr = LMESetOptionsPrefix(ctx->lme,((PetscObject)eps)->prefix);CHKERRQ(ierr);
    ierr = LMEAppendOptionsPrefix(ctx->lme,"eps_lyapii_");CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)ctx->lme,(PetscObject)eps,1);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)ctx->lme);CHKERRQ(ierr);
  }
  *lme = ctx->lme;
  PetscFunctionReturn(0);
}

/*@
   EPSLyapIIGetLME - Retrieve the linear matrix equation solver object (LME)
   associated with the eigenvalue solver.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  lme - the linear matrix equation solver object

   Level: advanced

.seealso: EPSLyapIISetLME()
@*/
PetscErrorCode EPSLyapIIGetLME(EPS eps,LME *lme)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(lme,2);
  ierr = PetscUseMethod(eps,"EPSLyapIIGetLME_C",(EPS,LME*),(eps,lme));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_LyapII(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  ranks: for Lyapunov solver=%D, after compression=%D\n",ctx->rkl,ctx->rkc);CHKERRQ(ierr);
    if (!ctx->lme) { ierr = EPSLyapIIGetLME(eps,&ctx->lme);CHKERRQ(ierr); }
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = LMEView(ctx->lme,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_LyapII(EPS eps)
{
  PetscErrorCode ierr;
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  if (!ctx->lme) { ierr = LMEReset(ctx->lme);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_LyapII(EPS eps)
{
  PetscErrorCode ierr;
  EPS_LYAPII     *ctx = (EPS_LYAPII*)eps->data;

  PetscFunctionBegin;
  ierr = LMEDestroy(&ctx->lme);CHKERRQ(ierr);
  ierr = DSDestroy(&ctx->ds);CHKERRQ(ierr);
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIISetLME_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIIGetLME_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIISetRanks_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIIGetRanks_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetDefaultST_LyapII(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) {
    ierr = STSetType(eps->st,STSINVERT);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_LyapII(EPS eps)
{
  EPS_LYAPII     *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&ctx);CHKERRQ(ierr);
  eps->data = (void*)ctx;

  eps->useds = PETSC_TRUE;

  eps->ops->solve          = EPSSolve_LyapII;
  eps->ops->setup          = EPSSetUp_LyapII;
  eps->ops->setfromoptions = EPSSetFromOptions_LyapII;
  eps->ops->reset          = EPSReset_LyapII;
  eps->ops->destroy        = EPSDestroy_LyapII;
  eps->ops->view           = EPSView_LyapII;
  eps->ops->setdefaultst   = EPSSetDefaultST_LyapII;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_Schur;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIISetLME_C",EPSLyapIISetLME_LyapII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIIGetLME_C",EPSLyapIIGetLME_LyapII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIISetRanks_C",EPSLyapIISetRanks_LyapII);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLyapIIGetRanks_C",EPSLyapIIGetRanks_LyapII);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

