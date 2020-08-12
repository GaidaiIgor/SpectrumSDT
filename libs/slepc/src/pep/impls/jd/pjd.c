/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc polynomial eigensolver: "jd"

   Method: Jacobi-Davidson

   Algorithm:

       Jacobi-Davidson for polynomial eigenvalue problems.

   References:

       [1] C. Campos and J.E. Roman, "A polynomial Jacobi-Davidson solver
           with support for non-monomial bases and deflation", BIT Numer.
           Math. (in press), 2019.

       [2] G.L.G. Sleijpen et al., "Jacobi-Davidson type methods for
           generalized eigenproblems and polynomial eigenproblems", BIT
           36(3):595-633, 1996.

       [3] Feng-Nan Hwang, Zih-Hao Wei, Tsung-Ming Huang, Weichung Wang,
           "A Parallel Additive Schwarz Preconditioned Jacobi-Davidson
           Algorithm for Polynomial Eigenvalue Problems in Quantum Dot
           Simulation", J. Comput. Phys. 229(8):2932-2947, 2010.
*/

#include <slepc/private/pepimpl.h>    /*I "slepcpep.h" I*/
#include <slepcblaslapack.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-slice-qep,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"A polynomial {Jacobi-Davidson} solver with support for non-monomial bases and deflation\",\n"
  "   journal = \"{BIT} Numer. Math.\",\n"
  "   volume = \"IP\",\n"
  "   number = \"-\",\n"
  "   pages = \"1--24\",\n"
  "   year = \"2019,\"\n"
  "   doi = \"https://doi.org/10.1007/s10543-019-00778-z\"\n"
  "}\n";

typedef struct {
  PetscReal   keep;          /* restart parameter */
  PetscReal   fix;           /* fix parameter */
  PetscBool   reusepc;       /* flag indicating whether pc is rebuilt or not */
  BV          V;             /* work basis vectors to store the search space */
  BV          W;             /* work basis vectors to store the test space */
  BV          *TV;           /* work basis vectors to store T*V (each TV[i] is the coefficient for \lambda^i of T*V for the extended T) */
  BV          *AX;           /* work basis vectors to store A_i*X for locked eigenvectors */
  BV          N[2];          /* auxiliary work BVs */
  BV          X;             /* locked eigenvectors */
  PetscScalar *T;            /* matrix of the invariant pair */
  PetscScalar *Tj;           /* matrix containing the powers of the invariant pair matrix */
  PetscScalar *XpX;          /* X^H*X */
  PetscInt    ld;            /* leading dimension for Tj and XpX */
  PC          pcshell;       /* preconditioner including basic precond+projector */
  Mat         Pshell;        /* auxiliary shell matrix */
  PetscInt    nlock;         /* number of locked vectors in the invariant pair */
  Vec         vtempl;        /* reference nested vector */
  PetscInt    midx;          /* minimality index */
  PetscInt    mmidx;         /* maximum allowed minimality index */
  PEPJDProjection proj;      /* projection type (orthogonal, harmonic) */
} PEP_JD;

typedef struct {
  PEP         pep;
  PC          pc;            /* basic preconditioner */
  Vec         Bp[2];         /* preconditioned residual of derivative polynomial, B\p */
  Vec         u[2];          /* Ritz vector */
  PetscScalar gamma[2];      /* precomputed scalar u'*B\p */
  PetscScalar theta;
  PetscScalar *M;
  PetscScalar *ps;
  PetscInt    ld;
  Vec         *work;
  Mat         PPr;
  BV          X;
  PetscInt    n;
} PEP_JD_PCSHELL;

typedef struct {
  Mat         Pr,Pi;         /* matrix polynomial evaluated at theta */
  PEP         pep;
  Vec         *work;
  PetscScalar theta[2];
} PEP_JD_MATSHELL;

/*
   Duplicate and resize auxiliary basis
*/
static PetscErrorCode PEPJDDuplicateBasis(PEP pep,BV *basis)
{
  PetscErrorCode     ierr;
  PEP_JD             *pjd = (PEP_JD*)pep->data;
  PetscInt           nloc,m;
  BVType             type;
  BVOrthogType       otype;
  BVOrthogRefineType oref;
  PetscReal          oeta;
  BVOrthogBlockType  oblock;

  PetscFunctionBegin;
  if (pjd->ld>1) {
    ierr = BVCreate(PetscObjectComm((PetscObject)pep),basis);CHKERRQ(ierr);
    ierr = BVGetSizes(pep->V,&nloc,NULL,&m);CHKERRQ(ierr);
    nloc += pjd->ld-1;
    ierr = BVSetSizes(*basis,nloc,PETSC_DECIDE,m);CHKERRQ(ierr);
    ierr = BVGetType(pep->V,&type);CHKERRQ(ierr);
    ierr = BVSetType(*basis,type);CHKERRQ(ierr);
    ierr = BVGetOrthogonalization(pep->V,&otype,&oref,&oeta,&oblock);CHKERRQ(ierr);
    ierr = BVSetOrthogonalization(*basis,otype,oref,oeta,oblock);CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease((PetscObject)*basis);CHKERRQ(ierr);
  } else {
    ierr = BVDuplicate(pep->V,basis);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSetUp_JD(PEP pep)
{
  PetscErrorCode ierr;
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscBool      isprecond,flg;
  PetscInt       i;

  PetscFunctionBegin;
  pep->lineariz = PETSC_FALSE;
  ierr = PEPSetDimensions_Default(pep,pep->nev,&pep->ncv,&pep->mpd);CHKERRQ(ierr);
  if (pep->max_it==PETSC_DEFAULT) pep->max_it = PetscMax(100,2*pep->n/pep->ncv);
  if (!pep->which) pep->which = PEP_TARGET_MAGNITUDE;
  if (pep->which!=PEP_TARGET_MAGNITUDE && pep->which!=PEP_TARGET_REAL && pep->which!=PEP_TARGET_IMAGINARY) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Wrong value of pep->which");

  ierr = PetscObjectTypeCompare((PetscObject)pep->st,STPRECOND,&isprecond);CHKERRQ(ierr);
  if (!isprecond) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"JD only works with PRECOND spectral transformation");

  ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
  if (flg) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Solver requires the ST transformation flag unset, see STSetTransform()");

  if (!pjd->mmidx) pjd->mmidx = pep->nmat-1;
  pjd->mmidx = PetscMin(pjd->mmidx,pep->nmat-1);
  if (!pjd->keep) pjd->keep = 0.5;
  ierr = PEPBasisCoefficients(pep,pep->pbc);CHKERRQ(ierr);
  ierr = PEPAllocateSolution(pep,0);CHKERRQ(ierr);
  ierr = PEPSetWorkVecs(pep,5);CHKERRQ(ierr);
  pjd->ld = pep->nev;
#if !defined (PETSC_USE_COMPLEX)
  pjd->ld++;
#endif
  ierr = PetscMalloc2(pep->nmat,&pjd->TV,pep->nmat,&pjd->AX);CHKERRQ(ierr);
  for (i=0;i<pep->nmat;i++) {
    ierr = PEPJDDuplicateBasis(pep,pjd->TV+i);CHKERRQ(ierr);
  }
  if (pjd->ld>1) {
    ierr = PEPJDDuplicateBasis(pep,&pjd->V);CHKERRQ(ierr);
    ierr = BVSetFromOptions(pjd->V);CHKERRQ(ierr);
    for (i=0;i<pep->nmat;i++) {
      ierr = BVDuplicateResize(pep->V,pjd->ld-1,pjd->AX+i);CHKERRQ(ierr);
    }
    ierr = BVDuplicateResize(pep->V,pjd->ld-1,pjd->N);CHKERRQ(ierr);
    ierr = BVDuplicateResize(pep->V,pjd->ld-1,pjd->N+1);CHKERRQ(ierr);
    pjd->X = pep->V;
    ierr = PetscCalloc3((pjd->ld)*(pjd->ld),&pjd->XpX,pep->ncv*pep->ncv,&pjd->T,pjd->ld*pjd->ld*pep->nmat,&pjd->Tj);CHKERRQ(ierr);
  } else pjd->V = pep->V;
  if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) { ierr = PEPJDDuplicateBasis(pep,&pjd->W);CHKERRQ(ierr); }
  else pjd->W = pjd->V;
  ierr = DSSetType(pep->ds,DSPEP);CHKERRQ(ierr);
  ierr = DSPEPSetDegree(pep->ds,pep->nmat-1);CHKERRQ(ierr);
  if (pep->basis!=PEP_BASIS_MONOMIAL) {
    ierr = DSPEPSetCoefficients(pep->ds,pep->pbc);CHKERRQ(ierr);
  }
  ierr = DSAllocate(pep->ds,pep->ncv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Updates columns (low to (high-1)) of TV[i]
*/
static PetscErrorCode PEPJDUpdateTV(PEP pep,PetscInt low,PetscInt high,Vec *w)
{
  PetscErrorCode ierr;
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscInt       pp,col,i,nloc,nconv;
  Vec            v1,v2,t1,t2;
  PetscScalar    *array1,*array2,*x2,*xx,*N,*Np,*y2=NULL,zero=0.0,sone=1.0,*pT,fact,*psc;
  PetscReal      *cg,*ca,*cb;
  PetscMPIInt    rk,np;
  PetscBLASInt   n_,ld_,one=1;
  Mat            T;
  BV             pbv;

  PetscFunctionBegin;
  ca = pep->pbc; cb = ca+pep->nmat; cg = cb + pep->nmat;
  nconv = pjd->nlock;
  ierr = PetscMalloc5(nconv,&x2,nconv,&xx,nconv*nconv,&pT,nconv*nconv,&N,nconv*nconv,&Np);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)pep),&rk);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np);CHKERRQ(ierr);
  ierr = BVGetSizes(pep->V,&nloc,NULL,NULL);CHKERRQ(ierr);
  t1 = w[0];
  t2 = w[1];
  ierr = PetscBLASIntCast(pjd->nlock,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(pjd->ld,&ld_);CHKERRQ(ierr);
  if (nconv){
    for (i=0;i<nconv;i++) {
      ierr = PetscArraycpy(pT+i*nconv,pjd->T+i*pep->ncv,nconv);CHKERRQ(ierr);
    }
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,nconv,nconv,pT,&T);CHKERRQ(ierr);
  }
  for (col=low;col<high;col++) {
    ierr = BVGetColumn(pjd->V,col,&v1);CHKERRQ(ierr);
    ierr = VecGetArray(v1,&array1);CHKERRQ(ierr);
    if (nconv>0) {
      for (i=0;i<nconv;i++) x2[i] = array1[nloc+i]* PetscSqrtReal(np);
    }
    ierr = VecPlaceArray(t1,array1);CHKERRQ(ierr);
    if (nconv) {
      ierr = BVSetActiveColumns(pjd->N[0],0,nconv);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(pjd->N[1],0,nconv);CHKERRQ(ierr);
      ierr = BVDotVec(pjd->X,t1,xx);CHKERRQ(ierr);
    }
    for (pp=pep->nmat-1;pp>=0;pp--) {
      ierr = BVGetColumn(pjd->TV[pp],col,&v2);CHKERRQ(ierr);
      ierr = VecGetArray(v2,&array2);CHKERRQ(ierr);
      ierr = VecPlaceArray(t2,array2);CHKERRQ(ierr);
      ierr = MatMult(pep->A[pp],t1,t2);CHKERRQ(ierr);
      if (nconv) {
        if (pp<pep->nmat-3) {
          ierr = BVMult(pjd->N[0],1.0,-cg[pp+2],pjd->AX[pp+1],NULL);CHKERRQ(ierr);
          ierr = MatShift(T,-cb[pp+1]);CHKERRQ(ierr);
          ierr = BVMult(pjd->N[0],1.0/ca[pp],1.0/ca[pp],pjd->N[1],T);CHKERRQ(ierr);
          pbv = pjd->N[0]; pjd->N[0] = pjd->N[1]; pjd->N[1] = pbv;
          ierr = BVMultVec(pjd->N[1],1.0,1.0,t2,x2);CHKERRQ(ierr);
          ierr = MatShift(T,cb[pp+1]);CHKERRQ(ierr);
        } else if (pp==pep->nmat-3) {
          ierr = BVCopy(pjd->AX[pp+2],pjd->N[0]);CHKERRQ(ierr);
          ierr = BVScale(pjd->N[0],1/ca[pp+1]);CHKERRQ(ierr);
          ierr = BVCopy(pjd->AX[pp+1],pjd->N[1]);CHKERRQ(ierr);
          ierr = MatShift(T,-cb[pp+1]);CHKERRQ(ierr);
          ierr = BVMult(pjd->N[1],1.0/ca[pp],1.0/ca[pp],pjd->N[0],T);CHKERRQ(ierr);
          ierr = BVMultVec(pjd->N[1],1.0,1.0,t2,x2);CHKERRQ(ierr);
          ierr = MatShift(T,cb[pp+1]);CHKERRQ(ierr);
        } else if (pp==pep->nmat-2) {
          ierr = BVMultVec(pjd->AX[pp+1],1.0/ca[pp],1.0,t2,x2);CHKERRQ(ierr);
        }
        if (pp<pjd->midx) {
          y2 = array2+nloc;
          PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&n_,&n_,&sone,pjd->Tj+pjd->ld*pjd->ld*pp,&ld_,xx,&one,&zero,y2,&one));
          if (pp<pjd->midx-2) {
            fact = -cg[pp+2];
            PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,pjd->Tj+(pp+1)*pjd->ld*pjd->ld,&ld_,pjd->XpX,&ld_,&fact,Np,&n_));
            fact = 1/ca[pp];
            ierr = MatShift(T,-cb[pp+1]);CHKERRQ(ierr);
            PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&fact,N,&n_,pT,&n_,&fact,Np,&n_));
            ierr = MatShift(T,cb[pp+1]);CHKERRQ(ierr);
            psc = Np; Np = N; N = psc;
            PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,N,&n_,x2,&one,&sone,y2,&one));
          } else if (pp==pjd->midx-2) {
            fact = 1/ca[pp];
            PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&fact,pjd->Tj+(pp+1)*pjd->ld*pjd->ld,&ld_,pjd->XpX,&ld_,&zero,N,&n_));
            PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,N,&n_,x2,&one,&sone,y2,&one));
          } else if (pp==pjd->midx-1) {
            ierr = PetscArrayzero(Np,nconv*nconv);CHKERRQ(ierr);
          }
        }
        for (i=0;i<nconv;i++) array2[nloc+i] /= PetscSqrtReal(np);
      }
      ierr = VecResetArray(t2);CHKERRQ(ierr);
      ierr = VecRestoreArray(v2,&array2);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pjd->TV[pp],col,&v2);CHKERRQ(ierr);
    }
    ierr = VecResetArray(t1);CHKERRQ(ierr);
    ierr = VecRestoreArray(v1,&array1);CHKERRQ(ierr);
    ierr = BVRestoreColumn(pjd->V,col,&v1);CHKERRQ(ierr);
  }
  if (nconv) {ierr = MatDestroy(&T);CHKERRQ(ierr);}
  ierr = PetscFree5(x2,xx,pT,N,Np);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   RRQR of X. Xin*P=Xou*R. Rank of R is rk
*/
static PetscErrorCode PEPJDOrthogonalize(PetscInt row,PetscInt col,PetscScalar *X,PetscInt ldx,PetscInt *rk,PetscInt *P,PetscScalar *R,PetscInt ldr)
{
  PetscErrorCode ierr;
  PetscInt       i,j,n,r;
  PetscBLASInt   row_,col_,ldx_,*p,lwork,info,n_;
  PetscScalar    *tau,*work;
  PetscReal      tol,*rwork;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(row,&row_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(col,&col_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldx,&ldx_);CHKERRQ(ierr);
  n = PetscMin(row,col);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  lwork = 3*col_+1;
  ierr = PetscMalloc4(col,&p,n,&tau,lwork,&work,2*col,&rwork);CHKERRQ(ierr);
  for (i=1;i<col;i++) p[i] = 0;
  p[0] = 1;

  /* rank revealing QR */
#if defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgeqp3",LAPACKgeqp3_(&row_,&col_,X,&ldx_,p,tau,work,&lwork,rwork,&info));
#else
  PetscStackCallBLAS("LAPACKgeqp3",LAPACKgeqp3_(&row_,&col_,X,&ldx_,p,tau,work,&lwork,&info));
#endif
  SlepcCheckLapackInfo("geqp3",info);
  if (P) for (i=0;i<col;i++) P[i] = p[i]-1;

  /* rank computation */
  tol = PetscMax(row,col)*PETSC_MACHINE_EPSILON*PetscAbsScalar(X[0]);
  r = 1;
  for (i=1;i<n;i++) {
    if (PetscAbsScalar(X[i+ldx*i])>tol) r++;
    else break;
  }
  if (rk) *rk=r;

  /* copy upper triangular matrix if requested */
  if (R) {
     for (i=0;i<r;i++) {
       ierr = PetscArrayzero(R+i*ldr,r);CHKERRQ(ierr);
       for (j=0;j<=i;j++) R[i*ldr+j] = X[i*ldx+j];
     }
  }
  PetscStackCallBLAS("LAPACKorgqr",LAPACKorgqr_(&row_,&n_,&n_,X,&ldx_,tau,work,&lwork,&info));
  SlepcCheckLapackInfo("orgqr",info);
  ierr = PetscFree4(p,tau,work,rwork);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Application of extended preconditioner
*/
static PetscErrorCode PEPJDExtendedPCApply(PC pc,Vec x,Vec y)
{
  PetscInt          i,j,nloc,n,ld=0;
  PetscMPIInt       np;
  Vec               tx,ty;
  PEP_JD_PCSHELL    *ctx;
  PetscErrorCode    ierr;
  const PetscScalar *array1;
  PetscScalar       *x2=NULL,*t=NULL,*ps=NULL,*array2,zero=0.0,sone=1.0;
  PetscBLASInt      one=1,ld_,n_,ncv_;
  PEP_JD            *pjd=NULL;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pc),&np);CHKERRQ(ierr);
  ierr = PCShellGetContext(pc,(void**)&ctx);CHKERRQ(ierr);
  n  = ctx->n;
  if (n) {
    pjd = (PEP_JD*)ctx->pep->data;
    ps = ctx->ps;
    ld = pjd->ld;
    ierr = PetscMalloc2(n,&x2,n,&t);CHKERRQ(ierr);
    ierr = VecGetLocalSize(ctx->work[0],&nloc);CHKERRQ(ierr);
    ierr = VecGetArrayRead(x,&array1);CHKERRQ(ierr);
    for (i=0;i<n;i++) x2[i] = array1[nloc+i]* PetscSqrtReal(np);
    ierr = VecRestoreArrayRead(x,&array1);CHKERRQ(ierr);
  }

  /* y = B\x apply PC */
  tx = ctx->work[0];
  ty = ctx->work[1];
  ierr = VecGetArrayRead(x,&array1);CHKERRQ(ierr);
  ierr = VecPlaceArray(tx,array1);CHKERRQ(ierr);
  ierr = VecGetArray(y,&array2);CHKERRQ(ierr);
  ierr = VecPlaceArray(ty,array2);CHKERRQ(ierr);
  ierr = PCApply(ctx->pc,tx,ty);CHKERRQ(ierr);
  if (n) {
    ierr = PetscBLASIntCast(ld,&ld_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
    for (i=0;i<n;i++) {
      t[i] = 0.0;
      for (j=0;j<n;j++) t[i] += ctx->M[i+j*ld]*x2[j];
    }
    if (pjd->midx==1) {
      ierr = PetscBLASIntCast(ctx->pep->ncv,&ncv_);CHKERRQ(ierr);
      for (i=0;i<n;i++) pjd->T[i*(1+ctx->pep->ncv)] -= ctx->theta;
      PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,pjd->T,&ncv_,t,&one,&zero,x2,&one));
      for (i=0;i<n;i++) pjd->T[i*(1+ctx->pep->ncv)] += ctx->theta;
      for (i=0;i<n;i++) array2[nloc+i] = x2[i];
      for (i=0;i<n;i++) x2[i] = -t[i];
    } else {
      for (i=0;i<n;i++) array2[nloc+i] = t[i];
      PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,ps,&ld_,t,&one,&zero,x2,&one));
    }
    for (i=0;i<n;i++) array2[nloc+i] /= PetscSqrtReal(np);
    ierr = BVSetActiveColumns(pjd->X,0,n);CHKERRQ(ierr);
    ierr = BVMultVec(pjd->X,-1.0,1.0,ty,x2);CHKERRQ(ierr);
    ierr = PetscFree2(x2,t);CHKERRQ(ierr);
  }
  ierr = VecResetArray(tx);CHKERRQ(ierr);
  ierr = VecResetArray(ty);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(x,&array1);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&array2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Application of shell preconditioner:
      y = B\x - eta*B\p,  with eta = (u'*B\x)/(u'*B\p)
*/
static PetscErrorCode PCShellApply_PEPJD(PC pc,Vec x,Vec y)
{
  PetscErrorCode ierr;
  PetscScalar    rr,eta;
  PEP_JD_PCSHELL *ctx;
  PetscInt       sz;
  const Vec      *xs,*ys;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    rx,xr,xx;
#endif

  PetscFunctionBegin;
  ierr = PCShellGetContext(pc,(void**)&ctx);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(x,&sz,&xs);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(y,NULL,&ys);CHKERRQ(ierr);
  /* y = B\x apply extended PC */
  ierr = PEPJDExtendedPCApply(pc,xs[0],ys[0]);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    ierr = PEPJDExtendedPCApply(pc,xs[1],ys[1]);CHKERRQ(ierr);
  }
#endif

  /* Compute eta = u'*y / u'*Bp */
  ierr = VecDot(ys[0],ctx->u[0],&rr);CHKERRQ(ierr);
  eta  = -rr*ctx->gamma[0];
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    ierr = VecDot(ys[0],ctx->u[1],&xr);CHKERRQ(ierr);
    ierr = VecDot(ys[1],ctx->u[0],&rx);CHKERRQ(ierr);
    ierr = VecDot(ys[1],ctx->u[1],&xx);CHKERRQ(ierr);
    eta += -ctx->gamma[0]*xx-ctx->gamma[1]*(-xr+rx);
  }
#endif
  eta /= ctx->gamma[0]*ctx->gamma[0]+ctx->gamma[1]*ctx->gamma[1];

  /* y = y - eta*Bp */
  ierr = VecAXPY(ys[0],eta,ctx->Bp[0]);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    ierr = VecAXPY(ys[1],eta,ctx->Bp[1]);CHKERRQ(ierr);
    eta = -ctx->gamma[1]*(rr+xx)+ctx->gamma[0]*(-xr+rx);
    eta /= ctx->gamma[0]*ctx->gamma[0]+ctx->gamma[1]*ctx->gamma[1];
    ierr = VecAXPY(ys[0],eta,ctx->Bp[1]);CHKERRQ(ierr);
    ierr = VecAXPY(ys[1],-eta,ctx->Bp[0]);CHKERRQ(ierr);
  }
#endif
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPJDCopyToExtendedVec(PEP pep,Vec v,PetscScalar *a,PetscInt na,PetscInt off,Vec vex,PetscBool back)
{
  PetscErrorCode ierr;
  PetscMPIInt    np,rk,count;
  PetscScalar    *array1,*array2;
  PetscInt       nloc;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)pep),&rk);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np);CHKERRQ(ierr);
  ierr = BVGetSizes(pep->V,&nloc,NULL,NULL);CHKERRQ(ierr);
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
      ierr = PetscArraycpy(a,array2+nloc+off,na);CHKERRQ(ierr);
      ierr = PetscMPIIntCast(na,&count);CHKERRQ(ierr);
      ierr = MPI_Bcast(a,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep));CHKERRQ(ierr);
    } else {
      ierr = PetscArraycpy(array2+nloc+off,a,na);CHKERRQ(ierr);
      ierr = PetscMPIIntCast(na,&count);CHKERRQ(ierr);
      ierr = MPI_Bcast(array2+nloc+off,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep));CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(vex,&array2);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* Computes Phi^hat(lambda) times a vector or its derivative (depends on beval)
     if no vector is provided returns a matrix
 */
static PetscErrorCode PEPJDEvaluateHatBasis(PEP pep,PetscInt n,PetscScalar *H,PetscInt ldh,PetscScalar *beval,PetscScalar *t,PetscInt idx,PetscScalar *qpp,PetscScalar *qp,PetscScalar *q)
{
  PetscErrorCode ierr;
  PetscInt       j,i;
  PetscBLASInt   n_,ldh_,one=1;
  PetscReal      *a,*b,*g;
  PetscScalar    sone=1.0,zero=0.0;

  PetscFunctionBegin;
  a = pep->pbc; b=a+pep->nmat; g=b+pep->nmat;
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldh,&ldh_);CHKERRQ(ierr);
  if (idx<1) {
    ierr = PetscArrayzero(q,t?n:n*n);CHKERRQ(ierr);
  } else if (idx==1) {
    if (t) {for (j=0;j<n;j++) q[j] = t[j]*beval[idx-1]/a[0];}
    else {
      ierr = PetscArrayzero(q,n*n);CHKERRQ(ierr);
      for (j=0;j<n;j++) q[(j+1)*n] = beval[idx-1]/a[0];
    }
  } else {
    if (t) {
      PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n_,&n_,&sone,H,&ldh_,qp,&one,&zero,q,&one));
      for (j=0;j<n;j++) {
        q[j] += beval[idx-1]*t[j]-b[idx-1]*qp[j]-g[idx-1]*qpp[j];
        q[j] /= a[idx-1];
      }
    } else {
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,H,&ldh_,qp,&n_,&zero,q,&n_));
      for (j=0;j<n;j++) {
        q[j+n*j] += beval[idx-1];
        for (i=0;i<n;i++) {
          q[i+n*j] += -b[idx-1]*qp[j*n+i]-g[idx-1]*qpp[j*n+i];
          q[i+n*j] /= a[idx-1];
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPJDComputeResidual(PEP pep,PetscBool derivative,PetscInt sz,Vec *u,PetscScalar *theta,Vec *p,Vec *work)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscErrorCode ierr;
  PetscMPIInt    rk,np,count;
  Vec            tu,tp,w;
  PetscScalar    *dval,*dvali,*array1,*array2,*x2=NULL,*y2,*qj=NULL,*tt=NULL,*xx=NULL,*xxi=NULL,sone=1.0;
  PetscInt       i,j,nconv,nloc;
  PetscBLASInt   n,ld,one=1;
#if !defined(PETSC_USE_COMPLEX)
  Vec            tui=NULL,tpi=NULL;
  PetscScalar    *x2i=NULL,*qji=NULL,*qq,*y2i,*arrayi1,*arrayi2;
#endif

  PetscFunctionBegin;
  nconv = pjd->nlock;
  if (!nconv) {
    ierr = PetscMalloc1(2*sz*pep->nmat,&dval);CHKERRQ(ierr);
  } else {
    ierr = PetscMalloc5(2*pep->nmat,&dval,2*nconv,&xx,nconv,&tt,sz*nconv,&x2,(sz==2?3:1)*nconv*pep->nmat,&qj);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)pep),&rk);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np);CHKERRQ(ierr);
    ierr = BVGetSizes(pep->V,&nloc,NULL,NULL);CHKERRQ(ierr);
    ierr = VecGetArray(u[0],&array1);CHKERRQ(ierr);
    for (i=0;i<nconv;i++) x2[i] = array1[nloc+i]*PetscSqrtReal(np);
    ierr = VecRestoreArray(u[0],&array1);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      x2i = x2+nconv;
      ierr = VecGetArray(u[1],&arrayi1);CHKERRQ(ierr);
      for (i=0;i<nconv;i++) x2i[i] = arrayi1[nloc+i]*PetscSqrtReal(np);
      ierr = VecRestoreArray(u[1],&arrayi1);CHKERRQ(ierr);
    }
#endif
  }
  dvali = dval+pep->nmat;
  tu = work[0];
  tp = work[1];
  w  = work[2];
  ierr = VecGetArray(u[0],&array1);CHKERRQ(ierr);
  ierr = VecPlaceArray(tu,array1);CHKERRQ(ierr);
  ierr = VecGetArray(p[0],&array2);CHKERRQ(ierr);
  ierr = VecPlaceArray(tp,array2);CHKERRQ(ierr);
  ierr = VecSet(tp,0.0);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    tui = work[3];
    tpi = work[4];
    ierr = VecGetArray(u[1],&arrayi1);CHKERRQ(ierr);
    ierr = VecPlaceArray(tui,arrayi1);CHKERRQ(ierr);
    ierr = VecGetArray(p[1],&arrayi2);CHKERRQ(ierr);
    ierr = VecPlaceArray(tpi,arrayi2);CHKERRQ(ierr);
    ierr = VecSet(tpi,0.0);CHKERRQ(ierr);
  }
#endif
  if (derivative) {
    ierr = PEPEvaluateBasisDerivative(pep,theta[0],theta[1],dval,dvali);CHKERRQ(ierr);
  } else {
    ierr = PEPEvaluateBasis(pep,theta[0],theta[1],dval,dvali);CHKERRQ(ierr);
  }
  for (i=derivative?1:0;i<pep->nmat;i++) {
    ierr = MatMult(pep->A[i],tu,w);CHKERRQ(ierr);
    ierr = VecAXPY(tp,dval[i],w);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      ierr = VecAXPY(tpi,dvali[i],w);CHKERRQ(ierr);
      ierr = MatMult(pep->A[i],tui,w);CHKERRQ(ierr);
      ierr = VecAXPY(tpi,dval[i],w);CHKERRQ(ierr);
      ierr = VecAXPY(tp,-dvali[i],w);CHKERRQ(ierr);
    }
#endif
  }
  if (nconv) {
    for (i=0;i<pep->nmat;i++) {
      ierr = PEPJDEvaluateHatBasis(pep,nconv,pjd->T,pep->ncv,dval,x2,i,i>1?qj+(i-2)*nconv:NULL,i>0?qj+(i-1)*nconv:NULL,qj+i*nconv);CHKERRQ(ierr);
    }
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      qji = qj+nconv*pep->nmat;
      qq = qji+nconv*pep->nmat;
      for (i=0;i<pep->nmat;i++) {
        ierr = PEPJDEvaluateHatBasis(pep,nconv,pjd->T,pep->ncv,dvali,x2i,i,i>1?qji+(i-2)*nconv:NULL,i>0?qji+(i-1)*nconv:NULL,qji+i*nconv);CHKERRQ(ierr);
      }
      for (i=0;i<nconv*pep->nmat;i++) qj[i] -= qji[i];
      for (i=0;i<pep->nmat;i++) {
        ierr = PEPJDEvaluateHatBasis(pep,nconv,pjd->T,pep->ncv,dval,x2i,i,i>1?qji+(i-2)*nconv:NULL,i>0?qji+(i-1)*nconv:NULL,qji+i*nconv);CHKERRQ(ierr);
        ierr = PEPJDEvaluateHatBasis(pep,nconv,pjd->T,pep->ncv,dvali,x2,i,i>1?qq+(i-2)*nconv:NULL,i>0?qq+(i-1)*nconv:NULL,qq+i*nconv);CHKERRQ(ierr);
      }
      for (i=0;i<nconv*pep->nmat;i++) qji[i] += qq[i];
      for (i=derivative?2:1;i<pep->nmat;i++) {
        ierr = BVMultVec(pjd->AX[i],1.0,1.0,tpi,qji+i*nconv);CHKERRQ(ierr);
      }
    }
#endif
    for (i=derivative?2:1;i<pep->nmat;i++) {
      ierr = BVMultVec(pjd->AX[i],1.0,1.0,tp,qj+i*nconv);CHKERRQ(ierr);
    }

    /* extended vector part */
    ierr = BVSetActiveColumns(pjd->X,0,nconv);CHKERRQ(ierr);
    ierr = BVDotVec(pjd->X,tu,xx);CHKERRQ(ierr);
    xxi = xx+nconv;
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) { ierr = BVDotVec(pjd->X,tui,xxi);CHKERRQ(ierr); }
#endif
    if (sz==1) { ierr = PetscArrayzero(xxi,nconv);CHKERRQ(ierr); }
    if (rk==np-1) {
      ierr = PetscBLASIntCast(nconv,&n);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(pjd->ld,&ld);CHKERRQ(ierr);
      y2  = array2+nloc;
      ierr = PetscArrayzero(y2,nconv);CHKERRQ(ierr);
      for (j=derivative?1:0;j<pjd->midx;j++) {
        for (i=0;i<nconv;i++) tt[i] = dval[j]*xx[i]-dvali[j]*xxi[i];
        PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,pjd->XpX,&ld,qj+j*nconv,&one,&sone,tt,&one));
        PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&sone,pjd->Tj+j*ld*ld,&ld,tt,&one,&sone,y2,&one));
      }
      for (i=0;i<nconv;i++) array2[nloc+i] /= PetscSqrtReal(np);
#if !defined(PETSC_USE_COMPLEX)
      if (sz==2) {
        y2i = arrayi2+nloc;
        ierr = PetscArrayzero(y2i,nconv);CHKERRQ(ierr);
        for (j=derivative?1:0;j<pjd->midx;j++) {
          for (i=0;i<nconv;i++) tt[i] = dval[j]*xxi[i]+dvali[j]*xx[i];
          PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,pjd->XpX,&ld,qji+j*nconv,&one,&sone,tt,&one));
          PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&sone,pjd->Tj+j*ld*ld,&ld,tt,&one,&sone,y2i,&one));
        }
        for (i=0;i<nconv;i++) arrayi2[nloc+i] /= PetscSqrtReal(np);
      }
#endif
    }
    ierr = PetscMPIIntCast(nconv,&count);CHKERRQ(ierr);
    ierr = MPI_Bcast(array2+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep));CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      ierr = MPI_Bcast(arrayi2+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep));CHKERRQ(ierr);
    }
#endif
  }
  if (nconv) {
    ierr = PetscFree5(dval,xx,tt,x2,qj);CHKERRQ(ierr);
  } else {
    ierr = PetscFree(dval);CHKERRQ(ierr);
  }
  ierr = VecResetArray(tu);CHKERRQ(ierr);
  ierr = VecRestoreArray(u[0],&array1);CHKERRQ(ierr);
  ierr = VecResetArray(tp);CHKERRQ(ierr);
  ierr = VecRestoreArray(p[0],&array2);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    ierr = VecResetArray(tui);CHKERRQ(ierr);
    ierr = VecRestoreArray(u[1],&arrayi1);CHKERRQ(ierr);
    ierr = VecResetArray(tpi);CHKERRQ(ierr);
    ierr = VecRestoreArray(p[1],&arrayi2);CHKERRQ(ierr);
  }
#endif
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPJDProcessInitialSpace(PEP pep,Vec *w)
{
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscErrorCode ierr;
  PetscScalar    *tt,target[2];
  Vec            vg,wg;
  PetscInt       i;
  PetscReal      norm;

  PetscFunctionBegin;
  ierr = PetscMalloc1(pjd->ld-1,&tt);CHKERRQ(ierr);
  if (pep->nini==0) {
    ierr = BVSetRandomColumn(pjd->V,0);CHKERRQ(ierr);
    for (i=0;i<pjd->ld-1;i++) tt[i] = 0.0;
    ierr = BVGetColumn(pjd->V,0,&vg);CHKERRQ(ierr);
    ierr = PEPJDCopyToExtendedVec(pep,NULL,tt,pjd->ld-1,0,vg,PETSC_FALSE);CHKERRQ(ierr);
    ierr = BVRestoreColumn(pjd->V,0,&vg);CHKERRQ(ierr);
    ierr = BVNormColumn(pjd->V,0,NORM_2,&norm);CHKERRQ(ierr);
    ierr = BVScaleColumn(pjd->V,0,1.0/norm);CHKERRQ(ierr);
    if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) {
      ierr = BVGetColumn(pjd->V,0,&vg);CHKERRQ(ierr);
      ierr = BVGetColumn(pjd->W,0,&wg);CHKERRQ(ierr);
      ierr = VecSet(wg,0.0);CHKERRQ(ierr);
      target[0] = pep->target; target[1] = 0.0;
      ierr = PEPJDComputeResidual(pep,PETSC_TRUE,1,&vg,target,&wg,w);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pjd->W,0,&wg);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pjd->V,0,&vg);CHKERRQ(ierr);
      ierr = BVNormColumn(pjd->W,0,NORM_2,&norm);CHKERRQ(ierr);
      ierr = BVScaleColumn(pjd->W,0,1.0/norm);CHKERRQ(ierr);
    }
  } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Support for initial vectors not implemented yet");
  ierr = PetscFree(tt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_PEPJD(Mat P,Vec x,Vec y)
{
  PetscErrorCode  ierr;
  PEP_JD_MATSHELL *matctx;
  PEP_JD          *pjd;
  PetscInt        i,j,nconv,nloc,nmat,ldt,ncv,sz;
  Vec             tx,ty;
  const Vec       *xs,*ys;
  PetscScalar     *array1,*array2,*x2=NULL,*y2,*tt=NULL,*xx=NULL,*xxi,theta[2],sone=1.0,*qj,*val,*vali=NULL;
  PetscBLASInt    n,ld,one=1;
  PetscMPIInt     np;
#if !defined(PETSC_USE_COMPLEX)
  Vec             txi=NULL,tyi=NULL;
  PetscScalar     *x2i=NULL,*qji=NULL,*qq,*y2i,*arrayi1,*arrayi2;
#endif

  PetscFunctionBegin;
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)P),&np);CHKERRQ(ierr);
  ierr  = MatShellGetContext(P,(void**)&matctx);CHKERRQ(ierr);
  pjd   = (PEP_JD*)(matctx->pep->data);
  nconv = pjd->nlock;
  nmat  = matctx->pep->nmat;
  ncv   = matctx->pep->ncv;
  ldt   = pjd->ld;
  ierr = VecCompGetSubVecs(x,&sz,&xs);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(y,NULL,&ys);CHKERRQ(ierr);
  theta[0] = matctx->theta[0];
  theta[1] = (sz==2)?matctx->theta[1]:0.0;
  if (nconv>0) {
    ierr = PetscMalloc5(nconv,&tt,sz*nconv,&x2,(sz==2?3:1)*nconv*nmat,&qj,2*nconv,&xx,2*nmat,&val);CHKERRQ(ierr);
    ierr = BVGetSizes(matctx->pep->V,&nloc,NULL,NULL);CHKERRQ(ierr);
    ierr = VecGetArray(xs[0],&array1);CHKERRQ(ierr);
    for (i=0;i<nconv;i++) x2[i] = array1[nloc+i]* PetscSqrtReal(np);
    ierr = VecRestoreArray(xs[0],&array1);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      x2i = x2+nconv;
      ierr = VecGetArray(xs[1],&arrayi1);CHKERRQ(ierr);
      for (i=0;i<nconv;i++) x2i[i] = arrayi1[nloc+i]* PetscSqrtReal(np);
      ierr = VecRestoreArray(xs[1],&arrayi1);CHKERRQ(ierr);
    }
#endif
    vali = val+nmat;
  }
  tx = matctx->work[0];
  ty = matctx->work[1];
  ierr = VecGetArray(xs[0],&array1);CHKERRQ(ierr);
  ierr = VecPlaceArray(tx,array1);CHKERRQ(ierr);
  ierr = VecGetArray(ys[0],&array2);CHKERRQ(ierr);
  ierr = VecPlaceArray(ty,array2);CHKERRQ(ierr);
  ierr = MatMult(matctx->Pr,tx,ty);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    txi = matctx->work[2];
    tyi = matctx->work[3];
    ierr = VecGetArray(xs[1],&arrayi1);CHKERRQ(ierr);
    ierr = VecPlaceArray(txi,arrayi1);CHKERRQ(ierr);
    ierr = VecGetArray(ys[1],&arrayi2);CHKERRQ(ierr);
    ierr = VecPlaceArray(tyi,arrayi2);CHKERRQ(ierr);
    ierr = MatMult(matctx->Pr,txi,tyi);CHKERRQ(ierr);
    if (theta[1]!=0.0) {
      ierr = MatMult(matctx->Pi,txi,matctx->work[4]);CHKERRQ(ierr);
      ierr = VecAXPY(ty,-1.0,matctx->work[4]);CHKERRQ(ierr);
      ierr = MatMult(matctx->Pi,tx,matctx->work[4]);CHKERRQ(ierr);
      ierr = VecAXPY(tyi,1.0,matctx->work[4]);CHKERRQ(ierr);
    }
  }
#endif
  if (nconv>0) {
    ierr = PEPEvaluateBasis(matctx->pep,theta[0],theta[1],val,vali);CHKERRQ(ierr);
    for (i=0;i<nmat;i++) {
      ierr = PEPJDEvaluateHatBasis(matctx->pep,nconv,pjd->T,ncv,val,x2,i,i>1?qj+(i-2)*nconv:NULL,i>0?qj+(i-1)*nconv:NULL,qj+i*nconv);CHKERRQ(ierr);
    }
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      qji = qj+nconv*nmat;
      qq = qji+nconv*nmat;
      for (i=0;i<nmat;i++) {
        ierr = PEPJDEvaluateHatBasis(matctx->pep,nconv,pjd->T,matctx->pep->ncv,vali,x2i,i,i>1?qji+(i-2)*nconv:NULL,i>0?qji+(i-1)*nconv:NULL,qji+i*nconv);CHKERRQ(ierr);
      }
      for (i=0;i<nconv*nmat;i++) qj[i] -= qji[i];
      for (i=0;i<nmat;i++) {
        ierr = PEPJDEvaluateHatBasis(matctx->pep,nconv,pjd->T,matctx->pep->ncv,val,x2i,i,i>1?qji+(i-2)*nconv:NULL,i>0?qji+(i-1)*nconv:NULL,qji+i*nconv);CHKERRQ(ierr);
        ierr = PEPJDEvaluateHatBasis(matctx->pep,nconv,pjd->T,matctx->pep->ncv,vali,x2,i,i>1?qq+(i-2)*nconv:NULL,i>0?qq+(i-1)*nconv:NULL,qq+i*nconv);CHKERRQ(ierr);
      }
      for (i=0;i<nconv*nmat;i++) qji[i] += qq[i];
      for (i=1;i<matctx->pep->nmat;i++) {
        ierr = BVMultVec(pjd->AX[i],1.0,1.0,tyi,qji+i*nconv);CHKERRQ(ierr);
      }
    }
#endif
    for (i=1;i<nmat;i++) {
      ierr = BVMultVec(pjd->AX[i],1.0,1.0,ty,qj+i*nconv);CHKERRQ(ierr);
    }

    /* extended vector part */
    ierr = BVSetActiveColumns(pjd->X,0,nconv);CHKERRQ(ierr);
    ierr = BVDotVec(pjd->X,tx,xx);CHKERRQ(ierr);
    xxi = xx+nconv;
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) { ierr = BVDotVec(pjd->X,txi,xxi);CHKERRQ(ierr); }
#endif
    if (sz==1) { ierr = PetscArrayzero(xxi,nconv);CHKERRQ(ierr); }
    ierr = PetscBLASIntCast(pjd->nlock,&n);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(ldt,&ld);CHKERRQ(ierr);
    y2 = array2+nloc;
    ierr = PetscArrayzero(y2,nconv);CHKERRQ(ierr);
    for (j=0;j<pjd->midx;j++) {
      for (i=0;i<nconv;i++) tt[i] = val[j]*xx[i]-vali[j]*xxi[i];
      PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,pjd->XpX,&ld,qj+j*nconv,&one,&sone,tt,&one));
      PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&sone,pjd->Tj+j*ld*ld,&ld,tt,&one,&sone,y2,&one));
    }
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      y2i = arrayi2+nloc;
      ierr = PetscArrayzero(y2i,nconv);CHKERRQ(ierr);
      for (j=0;j<pjd->midx;j++) {
        for (i=0;i<nconv;i++) tt[i] = val[j]*xxi[i]+vali[j]*xx[i];
        PetscStackCallBLAS("BLASgemv",BLASgemv_("N",&n,&n,&sone,pjd->XpX,&ld,qji+j*nconv,&one,&sone,tt,&one));
        PetscStackCallBLAS("BLASgemv",BLASgemv_("C",&n,&n,&sone,pjd->Tj+j*ld*ld,&ld,tt,&one,&sone,y2i,&one));
      }
      for (i=0;i<nconv;i++) arrayi2[nloc+i] /= PetscSqrtReal(np);
    }
#endif
    for (i=0;i<nconv;i++) array2[nloc+i] /= PetscSqrtReal(np);
    ierr = PetscFree5(tt,x2,qj,xx,val);CHKERRQ(ierr);
  }
  ierr = VecResetArray(tx);CHKERRQ(ierr);
  ierr = VecRestoreArray(xs[0],&array1);CHKERRQ(ierr);
  ierr = VecResetArray(ty);CHKERRQ(ierr);
  ierr = VecRestoreArray(ys[0],&array2);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    ierr = VecResetArray(txi);CHKERRQ(ierr);
    ierr = VecRestoreArray(xs[1],&arrayi1);CHKERRQ(ierr);
    ierr = VecResetArray(tyi);CHKERRQ(ierr);
    ierr = VecRestoreArray(ys[1],&arrayi2);CHKERRQ(ierr);
  }
#endif
  PetscFunctionReturn(0);
}

static PetscErrorCode MatCreateVecs_PEPJD(Mat A,Vec *right,Vec *left)
{
  PetscErrorCode  ierr;
  PEP_JD_MATSHELL *matctx;
  PEP_JD          *pjd;
  PetscInt        kspsf=1,i;
  Vec             v[2];

  PetscFunctionBegin;
  ierr  = MatShellGetContext(A,(void**)&matctx);CHKERRQ(ierr);
  pjd   = (PEP_JD*)(matctx->pep->data);
#if !defined (PETSC_USE_COMPLEX)
  kspsf = 2;
#endif
  for (i=0;i<kspsf;i++){
    ierr = BVCreateVec(pjd->V,v+i);CHKERRQ(ierr);
  }
  if (right) {
    ierr = VecCreateCompWithVecs(v,kspsf,pjd->vtempl,right);CHKERRQ(ierr);
  }
  if (left) {
    ierr = VecCreateCompWithVecs(v,kspsf,pjd->vtempl,left);CHKERRQ(ierr);
  }
  for (i=0;i<kspsf;i++) {
    ierr = VecDestroy(&v[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPJDUpdateExtendedPC(PEP pep,PetscScalar theta)
{
  PetscErrorCode ierr;
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PEP_JD_PCSHELL *pcctx;
  PetscInt       i,j,k,n=pjd->nlock,ld=pjd->ld,deg=pep->nmat-1;
  PetscScalar    *M,*ps,*work,*U,*V,*S,*Sp,*Spp,snone=-1.0,sone=1.0,zero=0.0,*val;
  PetscReal      tol,maxeig=0.0,*sg,*rwork;
  PetscBLASInt   n_,info,ld_,*p,lw_,rk=0;

  PetscFunctionBegin;
  if (n) {
    ierr = PCShellGetContext(pjd->pcshell,(void**)&pcctx);CHKERRQ(ierr);
    pcctx->theta = theta;
    pcctx->n = n;
    M  = pcctx->M;
    ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(ld,&ld_);CHKERRQ(ierr);
    if (pjd->midx==1) {
      ierr = PetscArraycpy(M,pjd->XpX,ld*ld);CHKERRQ(ierr);
      ierr = PetscCalloc2(10*n,&work,n,&p);CHKERRQ(ierr);
    } else {
      ps = pcctx->ps;
      ierr = PetscCalloc7(2*n*n,&U,3*n*n,&S,n,&sg,10*n,&work,5*n,&rwork,n,&p,deg+1,&val);CHKERRQ(ierr);
      V = U+n*n;
      /* pseudo-inverse */
      for (j=0;j<n;j++) {
        for (i=0;i<n;i++) S[n*j+i] = -pjd->T[pep->ncv*j+i];
        S[n*j+j] += theta;
      }
      lw_ = 10*n_;
#if !defined (PETSC_USE_COMPLEX)
      PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&n_,&n_,S,&n_,sg,U,&n_,V,&n_,work,&lw_,&info));
#else
      PetscStackCallBLAS("LAPACKgesvd",LAPACKgesvd_("S","S",&n_,&n_,S,&n_,sg,U,&n_,V,&n_,work,&lw_,rwork,&info));
#endif
      SlepcCheckLapackInfo("gesvd",info);
      for (i=0;i<n;i++) maxeig = PetscMax(maxeig,sg[i]);
      tol = 10*PETSC_MACHINE_EPSILON*n*maxeig;
      for (j=0;j<n;j++) {
        if (sg[j]>tol) {
          for (i=0;i<n;i++) U[j*n+i] /= sg[j];
          rk++;
        } else break;
      }
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&rk,&sone,U,&n_,V,&n_,&zero,ps,&ld_));

      /* compute M */
      ierr = PEPEvaluateBasis(pep,theta,0.0,val,NULL);CHKERRQ(ierr);
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&snone,pjd->XpX,&ld_,ps,&ld_,&zero,M,&ld_));
      ierr = PetscArrayzero(S,2*n*n);CHKERRQ(ierr);
      Sp = S+n*n;
      for (j=0;j<n;j++) S[j*(n+1)] = 1.0;
      for (k=1;k<pjd->midx;k++) {
        for (j=0;j<n;j++) for (i=0;i<n;i++) V[j*n+i] = S[j*n+i] - ps[j*ld+i]*val[k];
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,pjd->XpX,&ld_,V,&n_,&zero,U,&n_));
        PetscStackCallBLAS("BLASgemm",BLASgemm_("C","N",&n_,&n_,&n_,&sone,pjd->Tj+k*ld*ld,&ld_,U,&n_,&sone,M,&ld_));
        Spp = Sp; Sp = S;
        ierr = PEPJDEvaluateHatBasis(pep,n,pjd->T,pep->ncv,val,NULL,k+1,Spp,Sp,S);CHKERRQ(ierr);
      }
    }
    /* inverse */
    PetscStackCallBLAS("LAPACKgetrf",LAPACKgetrf_(&n_,&n_,M,&ld_,p,&info));
    SlepcCheckLapackInfo("getrf",info);
    PetscStackCallBLAS("LAPACKgetri",LAPACKgetri_(&n_,M,&ld_,p,work,&n_,&info));
    SlepcCheckLapackInfo("getri",info);
    if (pjd->midx==1) {
      ierr = PetscFree2(work,p);CHKERRQ(ierr);
    } else {
      ierr = PetscFree7(U,S,sg,work,rwork,p,val);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPJDMatSetUp(PEP pep,PetscInt sz,PetscScalar *theta)
{
  PetscErrorCode  ierr;
  PEP_JD          *pjd = (PEP_JD*)pep->data;
  PEP_JD_MATSHELL *matctx;
  PEP_JD_PCSHELL  *pcctx;
  MatStructure    str;
  PetscScalar     *vals,*valsi;
  PetscBool       skipmat=PETSC_FALSE;
  PetscInt        i;
  Mat             Pr=NULL;

  PetscFunctionBegin;
  if (sz==2 && theta[1]==0.0) sz = 1;
  ierr = MatShellGetContext(pjd->Pshell,(void**)&matctx);CHKERRQ(ierr);
  ierr = PCShellGetContext(pjd->pcshell,(void**)&pcctx);CHKERRQ(ierr);
  if (matctx->Pr && matctx->theta[0]==theta[0] && ((!matctx->Pi && sz==1) || (sz==2 && matctx->theta[1]==theta[1]))) {
    if (pcctx->n == pjd->nlock) PetscFunctionReturn(0);
    skipmat = PETSC_TRUE;
  }
  if (!skipmat) {
    ierr = PetscMalloc2(pep->nmat,&vals,pep->nmat,&valsi);CHKERRQ(ierr);
    ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
    ierr = PEPEvaluateBasis(pep,theta[0],theta[1],vals,valsi);CHKERRQ(ierr);
    if (!matctx->Pr) {
      ierr = MatDuplicate(pep->A[0],MAT_COPY_VALUES,&matctx->Pr);CHKERRQ(ierr);
    } else {
      ierr = MatCopy(pep->A[0],matctx->Pr,str);CHKERRQ(ierr);
    }
    for (i=1;i<pep->nmat;i++) {
      ierr = MatAXPY(matctx->Pr,vals[i],pep->A[i],str);CHKERRQ(ierr);
    }
    if (!pjd->reusepc) {
      if (pcctx->PPr && sz==2) {
        ierr = MatCopy(matctx->Pr,pcctx->PPr,str);CHKERRQ(ierr);
        Pr = pcctx->PPr;
      } else Pr = matctx->Pr;
    }
    matctx->theta[0] = theta[0];
#if !defined(PETSC_USE_COMPLEX)
    if (sz==2) {
      if (!matctx->Pi) {
        ierr = MatDuplicate(pep->A[0],MAT_COPY_VALUES,&matctx->Pi);CHKERRQ(ierr);
      } else {
        ierr = MatCopy(pep->A[1],matctx->Pi,str);CHKERRQ(ierr);
      }
      ierr = MatScale(matctx->Pi,valsi[1]);CHKERRQ(ierr);
      for (i=2;i<pep->nmat;i++) {
        ierr = MatAXPY(matctx->Pi,valsi[i],pep->A[i],str);CHKERRQ(ierr);
      }
      matctx->theta[1] = theta[1];
    }
#endif
    ierr = PetscFree2(vals,valsi);CHKERRQ(ierr);
  }
  if (!pjd->reusepc) {
    if (!skipmat) {
      ierr = PCSetOperators(pcctx->pc,Pr,Pr);CHKERRQ(ierr);
      ierr = PCSetUp(pcctx->pc);CHKERRQ(ierr);
    }
    ierr = PEPJDUpdateExtendedPC(pep,theta[0]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPJDCreateShellPC(PEP pep,Vec *ww)
{
  PEP_JD          *pjd = (PEP_JD*)pep->data;
  PEP_JD_PCSHELL  *pcctx;
  PEP_JD_MATSHELL *matctx;
  KSP             ksp;
  PetscInt        nloc,mloc,kspsf=1;
  Vec             v[2];
  PetscScalar     target[2];
  PetscErrorCode  ierr;
  Mat             Pr;

  PetscFunctionBegin;
  /* Create the reference vector */
  ierr = BVGetColumn(pjd->V,0,&v[0]);CHKERRQ(ierr);
  v[1] = v[0];
#if !defined (PETSC_USE_COMPLEX)
  kspsf = 2;
#endif
  ierr = VecCreateCompWithVecs(v,kspsf,NULL,&pjd->vtempl);CHKERRQ(ierr);
  ierr = BVRestoreColumn(pjd->V,0,&v[0]);CHKERRQ(ierr);
  ierr = PetscLogObjectParent((PetscObject)pep,(PetscObject)pjd->vtempl);CHKERRQ(ierr);

  /* Replace preconditioner with one containing projectors */
  ierr = PCCreate(PetscObjectComm((PetscObject)pep),&pjd->pcshell);CHKERRQ(ierr);
  ierr = PCSetType(pjd->pcshell,PCSHELL);CHKERRQ(ierr);
  ierr = PCShellSetName(pjd->pcshell,"PCPEPJD");CHKERRQ(ierr);
  ierr = PCShellSetApply(pjd->pcshell,PCShellApply_PEPJD);CHKERRQ(ierr);
  ierr = PetscNew(&pcctx);CHKERRQ(ierr);
  ierr = PCShellSetContext(pjd->pcshell,pcctx);CHKERRQ(ierr);
  ierr = STGetKSP(pep->st,&ksp);CHKERRQ(ierr);
  ierr = BVCreateVec(pjd->V,&pcctx->Bp[0]);CHKERRQ(ierr);
  ierr = VecDuplicate(pcctx->Bp[0],&pcctx->Bp[1]);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pcctx->pc);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)pcctx->pc);CHKERRQ(ierr);
  ierr = MatGetLocalSize(pep->A[0],&mloc,&nloc);CHKERRQ(ierr);
  if (pjd->ld>1) {
    nloc += pjd->ld-1; mloc += pjd->ld-1;
  }
  ierr = PetscNew(&matctx);CHKERRQ(ierr);
  ierr = MatCreateShell(PetscObjectComm((PetscObject)pep),kspsf*nloc,kspsf*mloc,PETSC_DETERMINE,PETSC_DETERMINE,matctx,&pjd->Pshell);CHKERRQ(ierr);
  ierr = MatShellSetOperation(pjd->Pshell,MATOP_MULT,(void(*)(void))MatMult_PEPJD);CHKERRQ(ierr);
  ierr = MatShellSetOperation(pjd->Pshell,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_PEPJD);CHKERRQ(ierr);
  matctx->pep = pep;
  target[0] = pep->target; target[1] = 0.0;
  ierr = PEPJDMatSetUp(pep,1,target);CHKERRQ(ierr);
  Pr = matctx->Pr;
  pcctx->PPr = NULL;
#if !defined(PETSC_USE_COMPLEX)
  if (!pjd->reusepc) {
    ierr = MatDuplicate(matctx->Pr,MAT_COPY_VALUES,&pcctx->PPr);CHKERRQ(ierr);
    Pr = pcctx->PPr;
  }
#endif
  ierr = PCSetOperators(pcctx->pc,Pr,Pr);CHKERRQ(ierr);
  ierr = PCSetErrorIfFailure(pcctx->pc,PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPSetPC(ksp,pjd->pcshell);CHKERRQ(ierr);
  if (pjd->reusepc) {
    ierr = PCSetReusePreconditioner(pcctx->pc,PETSC_TRUE);CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(ksp,PETSC_TRUE);CHKERRQ(ierr);
  }
  ierr = KSPSetOperators(ksp,pjd->Pshell,pjd->Pshell);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  if (pjd->ld>1) {
    ierr = PetscMalloc2(pjd->ld*pjd->ld,&pcctx->M,pjd->ld*pjd->ld,&pcctx->ps);CHKERRQ(ierr);
    pcctx->pep = pep;
  }
  matctx->work = ww;
  pcctx->work  = ww;
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPJDEigenvectors(PEP pep)
{
  PetscErrorCode ierr;
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscBLASInt   ld,nconv,info,nc;
  PetscScalar    *Z,*w;
  PetscReal      *wr,norm;
  PetscInt       i;
  Mat            U;
#if !defined(PETSC_USE_COMPLEX)
  Vec            v,v1;
#endif

  PetscFunctionBegin;
  ierr = PetscMalloc3(pep->nconv*pep->nconv,&Z,3*pep->ncv,&wr,2*pep->ncv,&w);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(pep->ncv,&ld);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(pep->nconv,&nconv);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKtrevc",LAPACKtrevc_("R","A",NULL,&nconv,pjd->T,&ld,NULL,&nconv,Z,&nconv,&nconv,&nc,wr,&info));
#else
  PetscStackCallBLAS("LAPACKtrevc",LAPACKtrevc_("R","A",NULL,&nconv,pjd->T,&ld,NULL,&nconv,Z,&nconv,&nconv,&nc,w,wr,&info));
#endif
  SlepcCheckLapackInfo("trevc",info);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,nconv,nconv,Z,&U);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(pjd->X,0,pep->nconv);CHKERRQ(ierr);
  ierr = BVMultInPlace(pjd->X,U,0,pep->nconv);CHKERRQ(ierr);
  for (i=0;i<pep->nconv;i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (pep->eigi[i]!=0.0) {   /* first eigenvalue of a complex conjugate pair */
      ierr = BVGetColumn(pjd->X,i,&v);CHKERRQ(ierr);
      ierr = BVGetColumn(pjd->X,i+1,&v1);CHKERRQ(ierr);
      ierr = VecNormalizeComplex(v,v1,PETSC_TRUE,NULL);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pjd->X,i,&v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pjd->X,i+1,&v1);CHKERRQ(ierr);
      i++;
    } else   /* real eigenvalue */
#endif
    {
      ierr = BVNormColumn(pjd->X,i,NORM_2,&norm);CHKERRQ(ierr);
      ierr = BVScaleColumn(pjd->X,i,1.0/norm);CHKERRQ(ierr);
    }
  }
  ierr = MatDestroy(&U);CHKERRQ(ierr);
  ierr = PetscFree3(Z,wr,w);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPJDLockConverged(PEP pep,PetscInt *nv,PetscInt sz)
{
  PetscErrorCode ierr;
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscInt       j,i,*P,ldds,rk=0,nvv=*nv;
  Vec            v,x,w;
  PetscScalar    *R,*r,*pX,target[2];
  Mat            X;
  PetscBLASInt   sz_,rk_,nv_,info;
  PetscMPIInt    np;

  PetscFunctionBegin;
  /* update AX and XpX */
  for (i=sz;i>0;i--) {
    ierr = BVGetColumn(pjd->X,pjd->nlock-i,&x);CHKERRQ(ierr);
    for (j=0;j<pep->nmat;j++) {
      ierr = BVGetColumn(pjd->AX[j],pjd->nlock-i,&v);CHKERRQ(ierr);
      ierr = MatMult(pep->A[j],x,v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pjd->AX[j],pjd->nlock-i,&v);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(pjd->AX[j],0,pjd->nlock-i+1);CHKERRQ(ierr);
    }
    ierr = BVRestoreColumn(pjd->X,pjd->nlock-i,&x);CHKERRQ(ierr);
    ierr = BVDotColumn(pjd->X,(pjd->nlock-i),pjd->XpX+(pjd->nlock-i)*(pjd->ld));CHKERRQ(ierr);
    pjd->XpX[(pjd->nlock-i)*(1+pjd->ld)] = 1.0;
    for (j=0;j<pjd->nlock-i;j++) pjd->XpX[j*(pjd->ld)+pjd->nlock-i] = PetscConj(pjd->XpX[(pjd->nlock-i)*(pjd->ld)+j]);
  }

  /* minimality index */
  pjd->midx = PetscMin(pjd->mmidx,pjd->nlock);

  /* evaluate the polynomial basis in T */
  ierr = PetscArrayzero(pjd->Tj,pjd->ld*pjd->ld*pep->nmat);CHKERRQ(ierr);
  for (j=0;j<pep->nmat;j++) {
    ierr = PEPEvaluateBasisMat(pep,pjd->nlock,pjd->T,pep->ncv,j,(j>1)?pjd->Tj+(j-2)*pjd->ld*pjd->ld:NULL,pjd->ld,j?pjd->Tj+(j-1)*pjd->ld*pjd->ld:NULL,pjd->ld,pjd->Tj+j*pjd->ld*pjd->ld,pjd->ld);CHKERRQ(ierr);
  }

  /* Extend search space */
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np);CHKERRQ(ierr);
  ierr = PetscCalloc3(nvv,&P,nvv*nvv,&R,nvv*sz,&r);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(pep->ds,&ldds);CHKERRQ(ierr);
  ierr = DSGetArray(pep->ds,DS_MAT_X,&pX);CHKERRQ(ierr);
  ierr = PEPJDOrthogonalize(nvv,nvv,pX,ldds,&rk,P,R,nvv);CHKERRQ(ierr);
  for (j=0;j<sz;j++) {
    for (i=0;i<rk;i++) r[i*sz+j] = PetscConj(R[nvv*i+j]*pep->eigr[P[i]]); /* first row scaled with permuted diagonal */
  }
  ierr = PetscBLASIntCast(rk,&rk_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(sz,&sz_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(nvv,&nv_);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&rk_,R,&nv_,&info));
  SlepcCheckLapackInfo("trtri",info);
  for (i=0;i<sz;i++) PetscStackCallBLAS("BLAStrmv",BLAStrmv_("U","C","N",&rk_,R,&nv_,r+i,&sz_));
  for (i=0;i<sz*rk;i++) r[i] = PetscConj(r[i])/PetscSqrtReal(np); /* revert */
  ierr = BVSetActiveColumns(pjd->V,0,nvv);CHKERRQ(ierr);
  rk -= sz;
  for (j=0;j<rk;j++) {
    ierr = PetscArraycpy(R+j*nvv,pX+(j+sz)*ldds,nvv);CHKERRQ(ierr);
  }
  ierr = DSRestoreArray(pep->ds,DS_MAT_X,&pX);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,nvv,rk,R,&X);CHKERRQ(ierr);
  ierr = BVMultInPlace(pjd->V,X,0,rk);CHKERRQ(ierr);
  ierr = MatDestroy(&X);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(pjd->V,0,rk);CHKERRQ(ierr);
  for (j=0;j<rk;j++) {
    ierr = BVGetColumn(pjd->V,j,&v);CHKERRQ(ierr);
    ierr = PEPJDCopyToExtendedVec(pep,NULL,r+sz*(j+sz),sz,pjd->nlock-sz,v,PETSC_FALSE);CHKERRQ(ierr);
    ierr = BVRestoreColumn(pjd->V,j,&v);CHKERRQ(ierr);
  }
  ierr = BVOrthogonalize(pjd->V,NULL);CHKERRQ(ierr);

  if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) {
    for (j=0;j<rk;j++) {
      /* W = P(target)*V */
      ierr = BVGetColumn(pjd->W,j,&w);CHKERRQ(ierr);
      ierr = BVGetColumn(pjd->V,j,&v);CHKERRQ(ierr);
      target[0] = pep->target; target[1] = 0.0;
      ierr = PEPJDComputeResidual(pep,PETSC_FALSE,1,&v,target,&w,pep->work);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pjd->V,j,&v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pjd->W,j,&w);CHKERRQ(ierr);
    }
    ierr = BVSetActiveColumns(pjd->W,0,rk);CHKERRQ(ierr);
    ierr = BVOrthogonalize(pjd->W,NULL);CHKERRQ(ierr);
  }
  *nv = rk;
  ierr = PetscFree3(P,R,r);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDSystemSetUp(PEP pep,PetscInt sz,PetscScalar *theta,Vec *u,Vec *p,Vec *ww)
{
  PetscErrorCode ierr;
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PEP_JD_PCSHELL *pcctx;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    s[2];
#endif

  PetscFunctionBegin;
  ierr = PCShellGetContext(pjd->pcshell,(void**)&pcctx);CHKERRQ(ierr);
  ierr = PEPJDMatSetUp(pep,sz,theta);CHKERRQ(ierr);
  pcctx->u[0] = u[0]; pcctx->u[1] = u[1];
  /* Compute r'. p is a work space vector */
  ierr = PEPJDComputeResidual(pep,PETSC_TRUE,sz,u,theta,p,ww);CHKERRQ(ierr);
  ierr = PEPJDExtendedPCApply(pjd->pcshell,p[0],pcctx->Bp[0]);CHKERRQ(ierr);
  ierr = VecDot(pcctx->Bp[0],u[0],pcctx->gamma);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  if (sz==2) {
    ierr = PEPJDExtendedPCApply(pjd->pcshell,p[1],pcctx->Bp[1]);CHKERRQ(ierr);
    ierr = VecDot(pcctx->Bp[0],u[1],pcctx->gamma+1);CHKERRQ(ierr);
    ierr = VecMDot(pcctx->Bp[1],2,u,s);CHKERRQ(ierr);
    pcctx->gamma[0] += s[1];
    pcctx->gamma[1] = -pcctx->gamma[1]+s[0];
  }
#endif
  if (sz==1) {
    ierr = VecZeroEntries(pcctx->Bp[1]);CHKERRQ(ierr);
    pcctx->gamma[1] = 0.0;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSolve_JD(PEP pep)
{
  PetscErrorCode  ierr;
  PEP_JD          *pjd = (PEP_JD*)pep->data;
  PetscInt        k,nv,nvc,ld,minv,dim,bupdated=0,sz=1,kspsf=1,idx,off,maxits,nloc;
  PetscMPIInt     np,count;
  PetscScalar     theta[2]={0.0,0.0},ritz[2]={0.0,0.0},*pX,*eig,*eigi,*array;
  PetscReal       norm,*res,tol=0.0,rtol,abstol, dtol;
  PetscBool       lindep,ini=PETSC_TRUE;
  Vec             tc,t[2]={NULL,NULL},u[2]={NULL,NULL},p[2]={NULL,NULL};
  Vec             rc,rr[2],r[2]={NULL,NULL},*ww=pep->work,v[2];
  Mat             G,X,Y;
  KSP             ksp;
  PEP_JD_PCSHELL  *pcctx;
  PEP_JD_MATSHELL *matctx;
#if !defined(PETSC_USE_COMPLEX)
  PetscReal       norm1;
#endif

  PetscFunctionBegin;
  ierr = PetscCitationsRegister(citation,&cited);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pep),&np);CHKERRQ(ierr);
  ierr = BVGetSizes(pep->V,&nloc,NULL,NULL);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(pep->ds,&ld);CHKERRQ(ierr);
  ierr = PetscCalloc3(pep->ncv+pep->nev,&eig,pep->ncv+pep->nev,&eigi,pep->ncv+pep->nev,&res);CHKERRQ(ierr);
  pjd->nlock = 0;
  ierr = STGetKSP(pep->st,&ksp);CHKERRQ(ierr);
  ierr = KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&maxits);CHKERRQ(ierr);
#if !defined (PETSC_USE_COMPLEX)
  kspsf = 2;
#endif
  ierr = PEPJDProcessInitialSpace(pep,ww);CHKERRQ(ierr);
  nv = (pep->nini)?pep->nini:1;

  /* Replace preconditioner with one containing projectors */
  ierr = PEPJDCreateShellPC(pep,ww);CHKERRQ(ierr);
  ierr = PCShellGetContext(pjd->pcshell,(void**)&pcctx);CHKERRQ(ierr);

  /* Create auxiliar vectors */
  ierr = BVCreateVec(pjd->V,&u[0]);CHKERRQ(ierr);
  ierr = VecDuplicate(u[0],&p[0]);CHKERRQ(ierr);
  ierr = VecDuplicate(u[0],&r[0]);CHKERRQ(ierr);
#if !defined (PETSC_USE_COMPLEX)
  ierr = VecDuplicate(u[0],&u[1]);CHKERRQ(ierr);
  ierr = VecDuplicate(u[0],&p[1]);CHKERRQ(ierr);
  ierr = VecDuplicate(u[0],&r[1]);CHKERRQ(ierr);
#endif

  /* Restart loop */
  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;
    ierr = DSSetDimensions(pep->ds,nv,0,0,0);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(pjd->V,bupdated,nv);CHKERRQ(ierr);
    ierr = PEPJDUpdateTV(pep,bupdated,nv,ww);CHKERRQ(ierr);
    if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) { ierr = BVSetActiveColumns(pjd->W,bupdated,nv);CHKERRQ(ierr); }
    for (k=0;k<pep->nmat;k++) {
      ierr = BVSetActiveColumns(pjd->TV[k],bupdated,nv);CHKERRQ(ierr);
      ierr = DSGetMat(pep->ds,DSMatExtra[k],&G);CHKERRQ(ierr);
      ierr = BVMatProject(pjd->TV[k],NULL,pjd->W,G);CHKERRQ(ierr);
      ierr = DSRestoreMat(pep->ds,DSMatExtra[k],&G);CHKERRQ(ierr);
    }
    ierr = BVSetActiveColumns(pjd->V,0,nv);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(pjd->W,0,nv);CHKERRQ(ierr);

    /* Solve projected problem */
    ierr = DSSetState(pep->ds,DS_STATE_RAW);CHKERRQ(ierr);
    ierr = DSSolve(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);
    ierr = DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DSSynchronize(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);
    idx = 0;
    do {
      ritz[0] = pep->eigr[idx];
#if !defined(PETSC_USE_COMPLEX)
      ritz[1] = pep->eigi[idx];
      sz = (ritz[1]==0.0)?1:2;
#endif
      /* Compute Ritz vector u=V*X(:,1) */
      ierr = DSGetArray(pep->ds,DS_MAT_X,&pX);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(pjd->V,0,nv);CHKERRQ(ierr);
      ierr = BVMultVec(pjd->V,1.0,0.0,u[0],pX+idx*ld);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
      if (sz==2) {
        ierr = BVMultVec(pjd->V,1.0,0.0,u[1],pX+(idx+1)*ld);CHKERRQ(ierr);
      }
#endif
      ierr = DSRestoreArray(pep->ds,DS_MAT_X,&pX);CHKERRQ(ierr);
      ierr = PEPJDComputeResidual(pep,PETSC_FALSE,sz,u,ritz,r,ww);CHKERRQ(ierr);
      /* Check convergence */
      ierr = VecNorm(r[0],NORM_2,&norm);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
      if (sz==2) {
        ierr = VecNorm(r[1],NORM_2,&norm1);CHKERRQ(ierr);
        norm = SlepcAbs(norm,norm1);
      }
#endif
      ierr = (*pep->converged)(pep,ritz[0],ritz[1],norm,&pep->errest[pep->nconv],pep->convergedctx);CHKERRQ(ierr);
      if (sz==2) pep->errest[pep->nconv+1] = pep->errest[pep->nconv];
      if (ini) {
        tol = PetscMin(.1,pep->errest[pep->nconv]); ini = PETSC_FALSE;
      } else tol = PetscMin(pep->errest[pep->nconv],tol/2);
      ierr = (*pep->stopping)(pep,pep->its,pep->max_it,(pep->errest[pep->nconv]<pep->tol)?pep->nconv+sz:pep->nconv,pep->nev,&pep->reason,pep->stoppingctx);CHKERRQ(ierr);
      if (pep->errest[pep->nconv]<pep->tol) {
        /* Ritz pair converged */
        ini = PETSC_TRUE;
        minv = PetscMin(nv,(PetscInt)(pjd->keep*pep->ncv));
        if (pjd->ld>1) {
          ierr = BVGetColumn(pjd->X,pep->nconv,&v[0]);CHKERRQ(ierr);
          ierr = PEPJDCopyToExtendedVec(pep,v[0],pjd->T+pep->ncv*pep->nconv,pjd->ld-1,0,u[0],PETSC_TRUE);CHKERRQ(ierr);
          ierr = BVRestoreColumn(pjd->X,pep->nconv,&v[0]);CHKERRQ(ierr);
          ierr = BVSetActiveColumns(pjd->X,0,pep->nconv+1);CHKERRQ(ierr);
          ierr = BVNormColumn(pjd->X,pep->nconv,NORM_2,&norm);CHKERRQ(ierr);
          ierr = BVScaleColumn(pjd->X,pep->nconv,1.0/norm);CHKERRQ(ierr);
          for (k=0;k<pep->nconv;k++) pjd->T[pep->ncv*pep->nconv+k] *= PetscSqrtReal(np)/norm;
          pjd->T[(pep->ncv+1)*pep->nconv] = ritz[0];
          eig[pep->nconv] = ritz[0];
          idx++;
#if !defined(PETSC_USE_COMPLEX)
          if (sz==2) {
            ierr = BVGetColumn(pjd->X,pep->nconv+1,&v[0]);CHKERRQ(ierr);
            ierr = PEPJDCopyToExtendedVec(pep,v[0],pjd->T+pep->ncv*(pep->nconv+1),pjd->ld-1,0,u[1],PETSC_TRUE);CHKERRQ(ierr);
            ierr = BVRestoreColumn(pjd->X,pep->nconv+1,&v[0]);CHKERRQ(ierr);
            ierr = BVSetActiveColumns(pjd->X,0,pep->nconv+2);CHKERRQ(ierr);
            ierr = BVNormColumn(pjd->X,pep->nconv+1,NORM_2,&norm1);CHKERRQ(ierr);
            ierr = BVScaleColumn(pjd->X,pep->nconv+1,1.0/norm1);CHKERRQ(ierr);
            for (k=0;k<pep->nconv;k++) pjd->T[pep->ncv*(pep->nconv+1)+k] *= PetscSqrtReal(np)/norm1;
            pjd->T[(pep->ncv+1)*(pep->nconv+1)] = ritz[0];
            pjd->T[(pep->ncv+1)*pep->nconv+1] = -ritz[1]*norm1/norm;
            pjd->T[(pep->ncv+1)*(pep->nconv+1)-1] = ritz[1]*norm/norm1;
            eig[pep->nconv+1] = ritz[0];
            eigi[pep->nconv] = ritz[1]; eigi[pep->nconv+1] = -ritz[1];
            idx++;
          }
#endif
        } else {
          ierr = BVInsertVec(pep->V,pep->nconv,u[0]);CHKERRQ(ierr);
        }
        pep->nconv += sz;
      }
    } while (pep->errest[pep->nconv]<pep->tol && pep->nconv<nv);

    if (pep->reason==PEP_CONVERGED_ITERATING) {
      nvc = nv;
      if (idx) {
        pjd->nlock +=idx;
        ierr = PEPJDLockConverged(pep,&nv,idx);CHKERRQ(ierr);
      }
      if (nv+sz>=pep->ncv-1) {
        /* Basis full, force restart */
        minv = PetscMin(nv,(PetscInt)(pjd->keep*pep->ncv));
        ierr = DSGetDimensions(pep->ds,&dim,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
        ierr = DSGetArray(pep->ds,DS_MAT_X,&pX);CHKERRQ(ierr);
        ierr = PEPJDOrthogonalize(dim,minv,pX,ld,&minv,NULL,NULL,ld);CHKERRQ(ierr);
        ierr = DSRestoreArray(pep->ds,DS_MAT_X,&pX);CHKERRQ(ierr);
        ierr = DSGetArray(pep->ds,DS_MAT_Y,&pX);CHKERRQ(ierr);
        ierr = PEPJDOrthogonalize(dim,minv,pX,ld,&minv,NULL,NULL,ld);CHKERRQ(ierr);
        ierr = DSRestoreArray(pep->ds,DS_MAT_Y,&pX);CHKERRQ(ierr);
        ierr = DSGetMat(pep->ds,DS_MAT_X,&X);CHKERRQ(ierr);
        ierr = BVMultInPlace(pjd->V,X,0,minv);CHKERRQ(ierr);
        if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) {
         ierr = DSGetMat(pep->ds,DS_MAT_Y,&Y);CHKERRQ(ierr);
         ierr = BVMultInPlace(pjd->W,Y,pep->nconv,minv);CHKERRQ(ierr);
         ierr = DSRestoreMat(pep->ds,DS_MAT_Y,&Y);CHKERRQ(ierr);
        }
        ierr = MatDestroy(&X);CHKERRQ(ierr);
        nv = minv;
        bupdated = 0;
      } else {
        if (!idx && pep->errest[pep->nconv]<pjd->fix) {theta[0] = ritz[0]; theta[1] = ritz[1];}
        else {theta[0] = pep->target; theta[1] = 0.0;}
        /* Update system mat */
        ierr = PEPJDSystemSetUp(pep,sz,theta,u,p,ww);CHKERRQ(ierr);
        /* Solve correction equation to expand basis */
        ierr = BVGetColumn(pjd->V,nv,&t[0]);CHKERRQ(ierr);
        rr[0] = r[0];
        if (sz==2) {
          ierr = BVGetColumn(pjd->V,nv+1,&t[1]);CHKERRQ(ierr);
          rr[1] = r[1];
        } else {
          t[1] = NULL;
          rr[1] = NULL;
        }
        ierr = VecCreateCompWithVecs(t,kspsf,pjd->vtempl,&tc);CHKERRQ(ierr);
        ierr = VecCreateCompWithVecs(rr,kspsf,pjd->vtempl,&rc);CHKERRQ(ierr);
        ierr = VecCompSetSubVecs(pjd->vtempl,sz,NULL);CHKERRQ(ierr);
        tol  = PetscMax(rtol,tol/2);
        ierr = KSPSetTolerances(ksp,tol,abstol,dtol,maxits);CHKERRQ(ierr);
        ierr = KSPSolve(ksp,rc,tc);CHKERRQ(ierr);
        ierr = VecDestroy(&tc);CHKERRQ(ierr);
        ierr = VecDestroy(&rc);CHKERRQ(ierr);
        ierr = VecGetArray(t[0],&array);CHKERRQ(ierr);
        ierr = PetscMPIIntCast(pep->nconv,&count);CHKERRQ(ierr);
        ierr = MPI_Bcast(array+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep));CHKERRQ(ierr);
        ierr = VecRestoreArray(t[0],&array);CHKERRQ(ierr);
        ierr = BVRestoreColumn(pjd->V,nv,&t[0]);CHKERRQ(ierr);
        ierr = BVOrthogonalizeColumn(pjd->V,nv,NULL,&norm,&lindep);CHKERRQ(ierr);
        if (lindep || norm==0.0) {
          if (sz==1) SETERRQ(PETSC_COMM_SELF,1,"Linearly dependent continuation vector");
          else off = 1;
        } else {
          off = 0;
          ierr = BVScaleColumn(pjd->V,nv,1.0/norm);CHKERRQ(ierr);
        }
#if !defined(PETSC_USE_COMPLEX)
        if (sz==2) {
          ierr = VecGetArray(t[1],&array);CHKERRQ(ierr);
          ierr = MPI_Bcast(array+nloc,count,MPIU_SCALAR,np-1,PetscObjectComm((PetscObject)pep));CHKERRQ(ierr);
          ierr = VecRestoreArray(t[1],&array);CHKERRQ(ierr);
          ierr = BVRestoreColumn(pjd->V,nv+1,&t[1]);CHKERRQ(ierr);
          if (off) {
            ierr = BVCopyColumn(pjd->V,nv+1,nv);CHKERRQ(ierr);
          }
          ierr = BVOrthogonalizeColumn(pjd->V,nv+1-off,NULL,&norm,&lindep);CHKERRQ(ierr);
          if (lindep || norm==0.0) {
            if (off) SETERRQ(PETSC_COMM_SELF,1,"Linearly dependent continuation vector");
            else off = 1;
          } else {
            ierr = BVScaleColumn(pjd->V,nv+1-off,1.0/norm);CHKERRQ(ierr);
          }
        }
#endif
        if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) {
          ierr = BVInsertVec(pjd->W,nv,r[0]);CHKERRQ(ierr);
          if (sz==2 && !off) {
            ierr = BVInsertVec(pjd->W,nv+1,r[1]);CHKERRQ(ierr);
          }
          ierr = BVOrthogonalizeColumn(pjd->W,nv,NULL,&norm,&lindep);CHKERRQ(ierr);
          if (lindep || norm==0.0) SETERRQ(PETSC_COMM_SELF,1,"Linearly dependent continuation vector");
          ierr = BVScaleColumn(pjd->W,nv,1.0/norm);CHKERRQ(ierr);
          if (sz==2 && !off) {
            ierr = BVOrthogonalizeColumn(pjd->W,nv+1,NULL,&norm,&lindep);CHKERRQ(ierr);
            if (lindep || norm==0.0) SETERRQ(PETSC_COMM_SELF,1,"Linearly dependent continuation vector");
            ierr = BVScaleColumn(pjd->W,nv+1,1.0/norm);CHKERRQ(ierr);
          }
        }
        bupdated = idx?0:nv;
        nv += sz-off;
      }
      for (k=0;k<nvc;k++) {
        eig[pep->nconv-idx+k] = pep->eigr[k];
#if !defined(PETSC_USE_COMPLEX)
        eigi[pep->nconv-idx+k] = pep->eigi[k];
#endif
      }
      ierr = PEPMonitor(pep,pep->its,pep->nconv,eig,eigi,pep->errest,pep->nconv+1);CHKERRQ(ierr);
    }
  }
  if (pjd->ld>1) {
    for (k=0;k<pep->nconv;k++) {
      pep->eigr[k] = eig[k];
      pep->eigi[k] = eigi[k];
    }
    if (pep->nconv>0) { ierr = PEPJDEigenvectors(pep);CHKERRQ(ierr); }
    ierr = PetscFree2(pcctx->M,pcctx->ps);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&u[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&r[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&p[0]);CHKERRQ(ierr);
#if !defined (PETSC_USE_COMPLEX)
  ierr = VecDestroy(&u[1]);CHKERRQ(ierr);
  ierr = VecDestroy(&r[1]);CHKERRQ(ierr);
  ierr = VecDestroy(&p[1]);CHKERRQ(ierr);
#endif
  ierr = KSPSetTolerances(ksp,rtol,abstol,dtol,maxits);CHKERRQ(ierr);
  ierr = KSPSetPC(ksp,pcctx->pc);CHKERRQ(ierr);
  ierr = VecDestroy(&pcctx->Bp[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&pcctx->Bp[1]);CHKERRQ(ierr);
  ierr = MatShellGetContext(pjd->Pshell,(void**)&matctx);CHKERRQ(ierr);
  ierr = MatDestroy(&matctx->Pr);CHKERRQ(ierr);
  ierr = MatDestroy(&matctx->Pi);CHKERRQ(ierr);
  ierr = MatDestroy(&pjd->Pshell);CHKERRQ(ierr);
  ierr = MatDestroy(&pcctx->PPr);CHKERRQ(ierr);
  ierr = PCDestroy(&pcctx->pc);CHKERRQ(ierr);
  ierr = PetscFree(pcctx);CHKERRQ(ierr);
  ierr = PetscFree(matctx);CHKERRQ(ierr);
  ierr = PCDestroy(&pjd->pcshell);CHKERRQ(ierr);
  ierr = PetscFree3(eig,eigi,res);CHKERRQ(ierr);
  ierr = VecDestroy(&pjd->vtempl);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDSetRestart_JD(PEP pep,PetscReal keep)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  if (keep==PETSC_DEFAULT) pjd->keep = 0.5;
  else {
    if (keep<0.1 || keep>0.9) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"The keep argument must be in the range [0.1,0.9]");
    pjd->keep = keep;
  }
  PetscFunctionReturn(0);
}

/*@
   PEPJDSetRestart - Sets the restart parameter for the Jacobi-Davidson
   method, in particular the proportion of basis vectors that must be kept
   after restart.

   Logically Collective on pep

   Input Parameters:
+  pep  - the eigenproblem solver context
-  keep - the number of vectors to be kept at restart

   Options Database Key:
.  -pep_jd_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,0.9]. The default is 0.5.

   Level: advanced

.seealso: PEPJDGetRestart()
@*/
PetscErrorCode PEPJDSetRestart(PEP pep,PetscReal keep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,keep,2);
  ierr = PetscTryMethod(pep,"PEPJDSetRestart_C",(PEP,PetscReal),(pep,keep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDGetRestart_JD(PEP pep,PetscReal *keep)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *keep = pjd->keep;
  PetscFunctionReturn(0);
}

/*@
   PEPJDGetRestart - Gets the restart parameter used in the Jacobi-Davidson method.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  keep - the restart parameter

   Level: advanced

.seealso: PEPJDSetRestart()
@*/
PetscErrorCode PEPJDGetRestart(PEP pep,PetscReal *keep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidRealPointer(keep,2);
  ierr = PetscUseMethod(pep,"PEPJDGetRestart_C",(PEP,PetscReal*),(pep,keep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDSetFix_JD(PEP pep,PetscReal fix)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  if (fix == PETSC_DEFAULT || fix == PETSC_DECIDE) pjd->fix = 0.01;
  else {
    if (fix < 0.0) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid fix value");
    pjd->fix = fix;
  }
  PetscFunctionReturn(0);
}

/*@
   PEPJDSetFix - Sets the threshold for changing the target in the correction
   equation.

   Logically Collective on pep

   Input Parameters:
+  pep - the eigenproblem solver context
-  fix - threshold for changing the target

   Options Database Key:
.  -pep_jd_fix - the fix value

   Note:
   The target in the correction equation is fixed at the first iterations.
   When the norm of the residual vector is lower than the fix value,
   the target is set to the corresponding eigenvalue.

   Level: advanced

.seealso: PEPJDGetFix()
@*/
PetscErrorCode PEPJDSetFix(PEP pep,PetscReal fix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,fix,2);
  ierr = PetscTryMethod(pep,"PEPJDSetFix_C",(PEP,PetscReal),(pep,fix));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDGetFix_JD(PEP pep,PetscReal *fix)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *fix = pjd->fix;
  PetscFunctionReturn(0);
}

/*@
   PEPJDGetFix - Returns the threshold for changing the target in the correction
   equation.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  fix - threshold for changing the target

   Note:
   The target in the correction equation is fixed at the first iterations.
   When the norm of the residual vector is lower than the fix value,
   the target is set to the corresponding eigenvalue.

   Level: advanced

.seealso: PEPJDSetFix()
@*/
PetscErrorCode PEPJDGetFix(PEP pep,PetscReal *fix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidRealPointer(fix,2);
  ierr = PetscUseMethod(pep,"PEPJDGetFix_C",(PEP,PetscReal*),(pep,fix));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDSetReusePreconditioner_JD(PEP pep,PetscBool reusepc)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  pjd->reusepc = reusepc;
  PetscFunctionReturn(0);
}

/*@
   PEPJDSetReusePreconditioner - Sets a flag indicating whether the preconditioner
   must be reused or not.

   Logically Collective on pep

   Input Parameters:
+  pep     - the eigenproblem solver context
-  reusepc - the reuse flag

   Options Database Key:
.  -pep_jd_reuse_preconditioner - the reuse flag

   Note:
   The default value is False. If set to True, the preconditioner is built
   only at the beginning, using the target value. Otherwise, it may be rebuilt
   (depending on the fix parameter) at each iteration from the Ritz value.

   Level: advanced

.seealso: PEPJDGetReusePreconditioner(), PEPJDSetFix()
@*/
PetscErrorCode PEPJDSetReusePreconditioner(PEP pep,PetscBool reusepc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,reusepc,2);
  ierr = PetscTryMethod(pep,"PEPJDSetReusePreconditioner_C",(PEP,PetscBool),(pep,reusepc));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDGetReusePreconditioner_JD(PEP pep,PetscBool *reusepc)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *reusepc = pjd->reusepc;
  PetscFunctionReturn(0);
}

/*@
   PEPJDGetReusePreconditioner - Returns the flag for reusing the preconditioner.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  reusepc - the reuse flag

   Level: advanced

.seealso: PEPJDSetReusePreconditioner()
@*/
PetscErrorCode PEPJDGetReusePreconditioner(PEP pep,PetscBool *reusepc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidBoolPointer(reusepc,2);
  ierr = PetscUseMethod(pep,"PEPJDGetReusePreconditioner_C",(PEP,PetscBool*),(pep,reusepc));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDSetMinimalityIndex_JD(PEP pep,PetscInt mmidx)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  if (mmidx == PETSC_DEFAULT || mmidx == PETSC_DECIDE) pjd->mmidx = 1;
  else {
    if (mmidx < 1) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmidx value");
    pjd->mmidx = mmidx;
    pep->state = PEP_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   PEPJDSetMinimalityIndex - Sets the maximum allowed value for the minimality index.

   Logically Collective on pep

   Input Parameters:
+  pep   - the eigenproblem solver context
-  mmidx - maximum minimality index

   Options Database Key:
.  -pep_jd_minimality_index - the minimality index value

   Note:
   The default value is equal to the degree of the polynomial. A smaller value
   can be used if the wanted eigenvectors are known to be linearly independent.

   Level: advanced

.seealso: PEPJDGetMinimalityIndex()
@*/
PetscErrorCode PEPJDSetMinimalityIndex(PEP pep,PetscInt mmidx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,mmidx,2);
  ierr = PetscTryMethod(pep,"PEPJDSetMinimalityIndex_C",(PEP,PetscInt),(pep,mmidx));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDGetMinimalityIndex_JD(PEP pep,PetscInt *mmidx)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *mmidx = pjd->mmidx;
  PetscFunctionReturn(0);
}

/*@
   PEPJDGetMinimalityIndex - Returns the maximum allowed value of the minimality
   index.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  mmidx - minimality index

   Level: advanced

.seealso: PEPJDSetMinimalityIndex()
@*/
PetscErrorCode PEPJDGetMinimalityIndex(PEP pep,PetscInt *mmidx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidIntPointer(mmidx,2);
  ierr = PetscUseMethod(pep,"PEPJDGetMinimalityIndex_C",(PEP,PetscInt*),(pep,mmidx));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDSetProjection_JD(PEP pep,PEPJDProjection proj)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  switch (proj) {
    case PEP_JD_PROJECTION_HARMONIC:
    case PEP_JD_PROJECTION_ORTHOGONAL:
      if (pjd->proj != proj) {
        pep->state = PEP_STATE_INITIAL;
        pjd->proj = proj;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'proj' value");
  }
  PetscFunctionReturn(0);
}

/*@
   PEPJDSetProjection - Sets the type of projection to be used in the Jacobi-Davidson solver.

   Logically Collective on pep

   Input Parameters:
+  pep  - the eigenproblem solver context
-  proj - the type of projection

   Options Database Key:
.  -pep_jd_projection - the projection type, either orthogonal or harmonic

   Level: advanced

.seealso: PEPJDGetProjection()
@*/
PetscErrorCode PEPJDSetProjection(PEP pep,PEPJDProjection proj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,proj,2);
  ierr = PetscTryMethod(pep,"PEPJDSetProjection_C",(PEP,PEPJDProjection),(pep,proj));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPJDGetProjection_JD(PEP pep,PEPJDProjection *proj)
{
  PEP_JD *pjd = (PEP_JD*)pep->data;

  PetscFunctionBegin;
  *proj = pjd->proj;
  PetscFunctionReturn(0);
}

/*@
   PEPJDGetProjection - Returns the type of projection used by the Jacobi-Davidson solver.

   Not Collective

   Input Parameter:
.  pep - the eigenproblem solver context

   Output Parameter:
.  proj - the type of projection

   Level: advanced

.seealso: PEPJDSetProjection()
@*/
PetscErrorCode PEPJDGetProjection(PEP pep,PEPJDProjection *proj)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(proj,2);
  ierr = PetscUseMethod(pep,"PEPJDGetProjection_C",(PEP,PEPJDProjection*),(pep,proj));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSetFromOptions_JD(PetscOptionItems *PetscOptionsObject,PEP pep)
{
  PetscErrorCode  ierr;
  PetscBool       flg,b1;
  PetscReal       r1;
  PetscInt        i1;
  PEPJDProjection proj;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"PEP JD Options");CHKERRQ(ierr);

    ierr = PetscOptionsReal("-pep_jd_restart","Proportion of vectors kept after restart","PEPJDSetRestart",0.5,&r1,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPJDSetRestart(pep,r1);CHKERRQ(ierr); }

    ierr = PetscOptionsReal("-pep_jd_fix","Tolerance for changing the target in the correction equation","PEPJDSetFix",0.01,&r1,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPJDSetFix(pep,r1);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-pep_jd_reuse_preconditioner","Whether to reuse the preconditioner","PEPJDSetReusePreconditoiner",PETSC_FALSE,&b1,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPJDSetReusePreconditioner(pep,b1);CHKERRQ(ierr); }

    ierr = PetscOptionsInt("-pep_jd_minimality_index","Maximum allowed minimality index","PEPJDSetMinimalityIndex",1,&i1,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPJDSetMinimalityIndex(pep,i1);CHKERRQ(ierr); }

    ierr = PetscOptionsEnum("-pep_jd_projection","Type of projection","PEPJDSetProjection",PEPJDProjectionTypes,(PetscEnum)PEP_JD_PROJECTION_HARMONIC,(PetscEnum*)&proj,&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPJDSetProjection(pep,proj);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPView_JD(PEP pep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  %d%% of basis vectors kept after restart\n",(int)(100*pjd->keep));CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  threshold for changing the target in the correction equation (fix): %g\n",(double)pjd->fix);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  projection type: %s\n",PEPJDProjectionTypes[pjd->proj]);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum allowed minimality index: %d\n",pjd->mmidx);CHKERRQ(ierr);
    if (pjd->reusepc) { ierr = PetscViewerASCIIPrintf(viewer,"  reusing the preconditioner\n");CHKERRQ(ierr); }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSetDefaultST_JD(PEP pep)
{
  PetscErrorCode ierr;
  KSP            ksp;

  PetscFunctionBegin;
  if (!((PetscObject)pep->st)->type_name) {
    ierr = STSetType(pep->st,STPRECOND);CHKERRQ(ierr);
    ierr = STPrecondSetKSPHasMat(pep->st,PETSC_TRUE);CHKERRQ(ierr);
  }
  ierr = STSetTransform(pep->st,PETSC_FALSE);CHKERRQ(ierr);
  ierr = STGetKSP(pep->st,&ksp);CHKERRQ(ierr);
  if (!((PetscObject)ksp)->type_name) {
    ierr = KSPSetType(ksp,KSPBCGSL);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1e-5,PETSC_DEFAULT,PETSC_DEFAULT,100);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPReset_JD(PEP pep)
{
  PetscErrorCode ierr;
  PEP_JD         *pjd = (PEP_JD*)pep->data;
  PetscInt       i;

  PetscFunctionBegin;
  for (i=0;i<pep->nmat;i++) {
    ierr = BVDestroy(pjd->TV+i);CHKERRQ(ierr);
  }
  if (pjd->proj==PEP_JD_PROJECTION_HARMONIC) { ierr = BVDestroy(&pjd->W);CHKERRQ(ierr); }
  if (pjd->ld>1) {
    ierr = BVDestroy(&pjd->V);CHKERRQ(ierr);
    for (i=0;i<pep->nmat;i++) {
      ierr = BVDestroy(pjd->AX+i);CHKERRQ(ierr);
    }
    ierr = BVDestroy(&pjd->N[0]);CHKERRQ(ierr);
    ierr = BVDestroy(&pjd->N[1]);CHKERRQ(ierr);
    ierr = PetscFree3(pjd->XpX,pjd->T,pjd->Tj);CHKERRQ(ierr);
  }
  ierr = PetscFree2(pjd->TV,pjd->AX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPDestroy_JD(PEP pep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(pep->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetFix_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetFix_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetReusePreconditioner_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetReusePreconditioner_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetMinimalityIndex_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetMinimalityIndex_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetProjection_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetProjection_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode PEPCreate_JD(PEP pep)
{
  PEP_JD         *pjd;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(pep,&pjd);CHKERRQ(ierr);
  pep->data = (void*)pjd;

  pjd->fix   = 0.01;
  pjd->mmidx = 0;

  pep->ops->solve          = PEPSolve_JD;
  pep->ops->setup          = PEPSetUp_JD;
  pep->ops->setfromoptions = PEPSetFromOptions_JD;
  pep->ops->destroy        = PEPDestroy_JD;
  pep->ops->reset          = PEPReset_JD;
  pep->ops->view           = PEPView_JD;
  pep->ops->setdefaultst   = PEPSetDefaultST_JD;

  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetRestart_C",PEPJDSetRestart_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetRestart_C",PEPJDGetRestart_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetFix_C",PEPJDSetFix_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetFix_C",PEPJDGetFix_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetReusePreconditioner_C",PEPJDSetReusePreconditioner_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetReusePreconditioner_C",PEPJDGetReusePreconditioner_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetMinimalityIndex_C",PEPJDSetMinimalityIndex_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetMinimalityIndex_C",PEPJDGetMinimalityIndex_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDSetProjection_C",PEPJDSetProjection_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)pep,"PEPJDGetProjection_C",PEPJDGetProjection_JD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

