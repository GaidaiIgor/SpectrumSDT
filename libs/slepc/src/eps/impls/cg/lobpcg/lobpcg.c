/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "lobpcg"

   Method: Locally Optimal Block Preconditioned Conjugate Gradient

   Algorithm:

       LOBPCG with soft and hard locking. Follows the implementation
       in BLOPEX [2].

   References:

       [1] A. V. Knyazev, "Toward the optimal preconditioned eigensolver:
           locally optimal block preconditioned conjugate gradient method",
           SIAM J. Sci. Comput. 23(2):517-541, 2001.

       [2] A. V. Knyazev et al., "Block Locally Optimal Preconditioned
           Eigenvalue Xolvers (BLOPEX) in Hypre and PETSc", SIAM J. Sci.
           Comput. 29(5):2224-2239, 2007.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/

typedef struct {
  PetscInt  bs;        /* block size */
  PetscBool lock;      /* soft locking active/inactive */
  PetscReal restart;   /* restart parameter */
  PetscInt  guard;     /* number of guard vectors */
} EPS_LOBPCG;

PetscErrorCode EPSSetDimensions_LOBPCG(EPS eps,PetscInt nev,PetscInt *ncv,PetscInt *mpd)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;
  PetscInt   k;

  PetscFunctionBegin;
  k = PetscMax(3*ctx->bs,((eps->nev-1)/ctx->bs+3)*ctx->bs);
  if (*ncv!=PETSC_DEFAULT) { /* ncv set */
    if (*ncv<k) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv is not sufficiently large");
  } else *ncv = k;
  if (*mpd==PETSC_DEFAULT) *mpd = 3*ctx->bs;
  else if (*mpd!=3*ctx->bs) SETERRQ(PetscObjectComm((PetscObject)eps),1,"This solver does not allow a value of mpd different from 3*blocksize");
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetUp_LOBPCG(EPS eps)
{
  PetscErrorCode ierr;
  EPS_LOBPCG     *ctx = (EPS_LOBPCG*)eps->data;
  PetscBool      istrivial;

  PetscFunctionBegin;
  if (!eps->ishermitian || (eps->isgeneralized && !eps->ispositive)) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"LOBPCG only works for Hermitian problems");
  if (!ctx->bs) ctx->bs = PetscMin(16,eps->nev);
  if (eps->n-eps->nds<5*ctx->bs) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The problem size is too small relative to the block size");
  ierr = EPSSetDimensions_LOBPCG(eps,eps->nev,&eps->ncv,&eps->mpd);CHKERRQ(ierr);
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) eps->which = EPS_SMALLEST_REAL;
  if (eps->which!=EPS_SMALLEST_REAL && eps->which!=EPS_LARGEST_REAL) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Wrong value of eps->which");
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");
  if (eps->extraction) { ierr = PetscInfo(eps,"Warning: extraction type ignored\n");CHKERRQ(ierr); }
  ierr = RGIsTrivial(eps->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support region filtering");

  if (!ctx->restart) ctx->restart = 0.9;

  /* number of guard vectors */
  if (ctx->bs==1) ctx->guard = 0;
  else ctx->guard = PetscMin((PetscInt)((1.0-ctx->restart)*ctx->bs+0.45),ctx->bs-1);

  ierr = EPSAllocateSolution(eps,0);CHKERRQ(ierr);
  ierr = EPS_SetInnerProduct(eps);CHKERRQ(ierr);
  ierr = DSSetType(eps->ds,DSGHEP);CHKERRQ(ierr);
  ierr = DSAllocate(eps->ds,eps->mpd);CHKERRQ(ierr);
  ierr = EPSSetWorkVecs(eps,1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_LOBPCG(EPS eps)
{
  PetscErrorCode ierr;
  EPS_LOBPCG     *ctx = (EPS_LOBPCG*)eps->data;
  PetscInt       i,j,k,ld,nv,ini,nmat,nc,nconv,locked,its;
  PetscReal      norm;
  PetscScalar    *eigr,dot;
  PetscBool      breakdown,countc,flip=PETSC_FALSE,checkprecond=PETSC_FALSE;
  Mat            A,B,M;
  Vec            v,z,w=eps->work[0];
  BV             X,Y=NULL,Z,R,P,AX,BX;
  SlepcSC        sc;

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  ierr = STGetNumMatrices(eps->st,&nmat);CHKERRQ(ierr);
  ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
  if (nmat>1) { ierr = STGetMatrix(eps->st,1,&B);CHKERRQ(ierr); }
  else B = NULL;

  if (eps->which==EPS_LARGEST_REAL) {  /* flip spectrum */
    flip = PETSC_TRUE;
    ierr = DSGetSlepcSC(eps->ds,&sc);CHKERRQ(ierr);
    sc->comparison = SlepcCompareSmallestReal;
  }

  /* undocumented option to check for a positive-definite preconditioner (turn-off by default) */
  ierr = PetscOptionsGetBool(NULL,NULL,"-eps_lobpcg_checkprecond",&checkprecond,NULL);CHKERRQ(ierr);

  /* 1. Allocate memory */
  ierr = PetscCalloc1(3*ctx->bs,&eigr);CHKERRQ(ierr);
  ierr = BVDuplicateResize(eps->V,3*ctx->bs,&Z);CHKERRQ(ierr);
  ierr = BVDuplicateResize(eps->V,ctx->bs,&X);CHKERRQ(ierr);
  ierr = BVDuplicateResize(eps->V,ctx->bs,&R);CHKERRQ(ierr);
  ierr = BVDuplicateResize(eps->V,ctx->bs,&P);CHKERRQ(ierr);
  ierr = BVDuplicateResize(eps->V,ctx->bs,&AX);CHKERRQ(ierr);
  if (B) {
    ierr = BVDuplicateResize(eps->V,ctx->bs,&BX);CHKERRQ(ierr);
  }
  nc = eps->nds;
  if (nc>0 || eps->nev>ctx->bs-ctx->guard) {
    ierr = BVDuplicateResize(eps->V,nc+eps->nev,&Y);CHKERRQ(ierr);
  }
  if (nc>0) {
    for (j=0;j<nc;j++) {
      ierr = BVGetColumn(eps->V,-nc+j,&v);CHKERRQ(ierr);
      ierr = BVInsertVec(Y,j,v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(eps->V,-nc+j,&v);CHKERRQ(ierr);
    }
    ierr = BVSetActiveColumns(Y,0,nc);CHKERRQ(ierr);
  }

  /* 2. Apply the constraints to the initial vectors */
  /* 3. B-orthogonalize initial vectors */
  for (k=eps->nini;k<eps->ncv-ctx->bs;k++) { /* Generate more initial vectors if necessary */
    ierr = BVSetRandomColumn(eps->V,k);CHKERRQ(ierr);
    ierr = BVOrthonormalizeColumn(eps->V,k,PETSC_TRUE,NULL,NULL);CHKERRQ(ierr);
  }
  nv = ctx->bs;
  ierr = BVSetActiveColumns(eps->V,0,nv);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(Z,0,nv);CHKERRQ(ierr);
  ierr = BVCopy(eps->V,Z);CHKERRQ(ierr);
  ierr = BVCopy(Z,X);CHKERRQ(ierr);

  /* 4. Compute initial Ritz vectors */
  ierr = BVMatMult(X,A,AX);CHKERRQ(ierr);
  ierr = DSSetDimensions(eps->ds,nv,0,0,0);CHKERRQ(ierr);
  ierr = DSGetMat(eps->ds,DS_MAT_A,&M);CHKERRQ(ierr);
  ierr = BVMatProject(AX,NULL,X,M);CHKERRQ(ierr);
  if (flip) { ierr = MatScale(M,-1.0);CHKERRQ(ierr); }
  ierr = DSRestoreMat(eps->ds,DS_MAT_A,&M);CHKERRQ(ierr);
  ierr = DSSetIdentity(eps->ds,DS_MAT_B);CHKERRQ(ierr);
  ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
  ierr = DSSolve(eps->ds,eigr,NULL);CHKERRQ(ierr);
  ierr = DSSort(eps->ds,eigr,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DSSynchronize(eps->ds,eigr,NULL);CHKERRQ(ierr);
  for (j=0;j<nv;j++) eps->eigr[j] = flip? -eigr[j]: eigr[j];
  ierr = DSVectors(eps->ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);
  ierr = DSGetMat(eps->ds,DS_MAT_X,&M);CHKERRQ(ierr);
  ierr = BVMultInPlace(X,M,0,nv);CHKERRQ(ierr);
  ierr = BVMultInPlace(AX,M,0,nv);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);

  /* 5. Initialize range of active iterates */
  locked = 0;  /* hard-locked vectors, the leading locked columns of V are eigenvectors */
  nconv  = 0;  /* number of converged eigenvalues in the current block */
  its    = 0;  /* iterations for the current block */

  /* 6. Main loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {

    if (ctx->lock) {
      ierr = BVSetActiveColumns(R,nconv,ctx->bs);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(AX,nconv,ctx->bs);CHKERRQ(ierr);
      if (B) {
        ierr = BVSetActiveColumns(BX,nconv,ctx->bs);CHKERRQ(ierr);
      }
    }

    /* 7. Compute residuals */
    ini = (ctx->lock)? nconv: 0;
    ierr = BVCopy(AX,R);CHKERRQ(ierr);
    if (B) { ierr = BVMatMult(X,B,BX);CHKERRQ(ierr); }
    for (j=ini;j<ctx->bs;j++) {
      ierr = BVGetColumn(R,j,&v);CHKERRQ(ierr);
      ierr = BVGetColumn(B?BX:X,j,&z);CHKERRQ(ierr);
      ierr = VecAXPY(v,-eps->eigr[locked+j],z);CHKERRQ(ierr);
      ierr = BVRestoreColumn(R,j,&v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(B?BX:X,j,&z);CHKERRQ(ierr);
    }

    /* 8. Compute residual norms and update index set of active iterates */
    k = ini;
    countc = PETSC_TRUE;
    for (j=ini;j<ctx->bs;j++) {
      i = locked+j;
      ierr = BVGetColumn(R,j,&v);CHKERRQ(ierr);
      ierr = VecNorm(v,NORM_2,&norm);CHKERRQ(ierr);
      ierr = BVRestoreColumn(R,j,&v);CHKERRQ(ierr);
      ierr = (*eps->converged)(eps,eps->eigr[i],eps->eigi[i],norm,&eps->errest[i],eps->convergedctx);CHKERRQ(ierr);
      if (countc) {
        if (eps->errest[i] < eps->tol) k++;
        else countc = PETSC_FALSE;
      }
      if (!countc && !eps->trackall) break;
    }
    nconv = k;
    eps->nconv = locked + nconv;
    if (its) {
      ierr = EPSMonitor(eps,eps->its+its,eps->nconv,eps->eigr,eps->eigi,eps->errest,locked+ctx->bs);CHKERRQ(ierr);
    }
    ierr = (*eps->stopping)(eps,eps->its+its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);
    if (eps->reason != EPS_CONVERGED_ITERATING || nconv >= ctx->bs-ctx->guard) {
      ierr = BVSetActiveColumns(eps->V,locked,eps->nconv);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(X,0,nconv);CHKERRQ(ierr);
      ierr = BVCopy(X,eps->V);CHKERRQ(ierr);
    }
    if (eps->reason != EPS_CONVERGED_ITERATING) {
      break;
    } else if (nconv >= ctx->bs-ctx->guard) {
      eps->its += its-1;
      its = 0;
    } else its++;

    if (nconv >= ctx->bs-ctx->guard) {  /* force hard locking of vectors and compute new R */

      /* extend constraints */
      ierr = BVSetActiveColumns(Y,nc+locked,nc+locked+nconv);CHKERRQ(ierr);
      ierr = BVCopy(X,Y);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(Y,0,nc+locked+nconv);CHKERRQ(ierr);

      /* shift work BV's */
      for (j=nconv;j<ctx->bs;j++) {
        ierr = BVCopyColumn(X,j,j-nconv);CHKERRQ(ierr);
        ierr = BVCopyColumn(R,j,j-nconv);CHKERRQ(ierr);
        ierr = BVCopyColumn(P,j,j-nconv);CHKERRQ(ierr);
        ierr = BVCopyColumn(AX,j,j-nconv);CHKERRQ(ierr);
        if (B) {
          ierr = BVCopyColumn(BX,j,j-nconv);CHKERRQ(ierr);
        }
      }

      /* set new initial vectors */
      ierr = BVSetActiveColumns(eps->V,locked+ctx->bs,locked+ctx->bs+nconv);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(X,ctx->bs-nconv,ctx->bs);CHKERRQ(ierr);
      ierr = BVCopy(eps->V,X);CHKERRQ(ierr);
      for (j=ctx->bs-nconv;j<ctx->bs;j++) {
        ierr = BVGetColumn(X,j,&v);CHKERRQ(ierr);
        ierr = BVOrthogonalizeVec(Y,v,NULL,&norm,&breakdown);CHKERRQ(ierr);
        if (norm>0.0 && !breakdown) {
          ierr = VecScale(v,1.0/norm);CHKERRQ(ierr);
        } else {
          ierr = PetscInfo(eps,"Orthogonalization of initial vector failed\n");CHKERRQ(ierr);
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          goto diverged;
        }
        ierr = BVRestoreColumn(X,j,&v);CHKERRQ(ierr);
      }
      locked += nconv;
      nconv = 0;
      ierr = BVSetActiveColumns(X,nconv,ctx->bs);CHKERRQ(ierr);

      /* B-orthogonalize initial vectors */
      ierr = BVOrthogonalize(X,NULL);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(Z,nconv,ctx->bs);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(AX,nconv,ctx->bs);CHKERRQ(ierr);
      ierr = BVCopy(X,Z);CHKERRQ(ierr);

      /* compute initial Ritz vectors */
      nv = ctx->bs;
      ierr = BVMatMult(X,A,AX);CHKERRQ(ierr);
      ierr = DSSetDimensions(eps->ds,nv,0,0,0);CHKERRQ(ierr);
      ierr = DSGetMat(eps->ds,DS_MAT_A,&M);CHKERRQ(ierr);
      ierr = BVMatProject(AX,NULL,X,M);CHKERRQ(ierr);
      if (flip) { ierr = MatScale(M,-1.0);CHKERRQ(ierr); }
      ierr = DSRestoreMat(eps->ds,DS_MAT_A,&M);CHKERRQ(ierr);
      ierr = DSSetIdentity(eps->ds,DS_MAT_B);CHKERRQ(ierr);
      ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
      ierr = DSSolve(eps->ds,eigr,NULL);CHKERRQ(ierr);
      ierr = DSSort(eps->ds,eigr,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
      ierr = DSSynchronize(eps->ds,eigr,NULL);CHKERRQ(ierr);
      for (j=0;j<nv;j++) if (locked+j<eps->ncv) eps->eigr[locked+j] = flip? -eigr[j]: eigr[j];
      ierr = DSVectors(eps->ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);
      ierr = DSGetMat(eps->ds,DS_MAT_X,&M);CHKERRQ(ierr);
      ierr = BVMultInPlace(X,M,0,nv);CHKERRQ(ierr);
      ierr = BVMultInPlace(AX,M,0,nv);CHKERRQ(ierr);
      ierr = MatDestroy(&M);CHKERRQ(ierr);

      continue;   /* skip the rest of the iteration */
    }

    ini = (ctx->lock)? nconv: 0;
    if (ctx->lock) {
      ierr = BVSetActiveColumns(R,nconv,ctx->bs);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(P,nconv,ctx->bs);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(AX,nconv,ctx->bs);CHKERRQ(ierr);
      if (B) {
        ierr = BVSetActiveColumns(BX,nconv,ctx->bs);CHKERRQ(ierr);
      }
    }

    /* 9. Apply preconditioner to the residuals */
    for (j=ini;j<ctx->bs;j++) {
      ierr = BVGetColumn(R,j,&v);CHKERRQ(ierr);
      ierr = STApply(eps->st,v,w);CHKERRQ(ierr);
      if (checkprecond) {
        ierr = VecDot(v,w,&dot);CHKERRQ(ierr);
        if (PetscRealPart(dot)<0.0) {
          ierr = PetscInfo(eps,"The preconditioner is not positive-definite\n");CHKERRQ(ierr);
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          goto diverged;
        }
      }
      if (nc+locked>0) {
        ierr = BVOrthogonalizeVec(Y,w,NULL,&norm,&breakdown);CHKERRQ(ierr);
        if (norm>0.0 && !breakdown) {
          ierr = VecScale(w,1.0/norm);CHKERRQ(ierr);
        } else {
          ierr = PetscInfo(eps,"Orthogonalization of preconditioned residual failed\n");CHKERRQ(ierr);
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          goto diverged;
        }
      }
      ierr = VecCopy(w,v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(R,j,&v);CHKERRQ(ierr);
    }

    /* 11. B-orthonormalize preconditioned residuals */
    ierr = BVOrthogonalize(R,NULL);CHKERRQ(ierr);

    /* 13-16. B-orthonormalize conjugate directions */
    if (its>1) {
      ierr = BVOrthogonalize(P,NULL);CHKERRQ(ierr);
    }

    /* 17-23. Compute symmetric Gram matrices */
    ierr = BVSetActiveColumns(Z,0,ctx->bs);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(X,0,ctx->bs);CHKERRQ(ierr);
    ierr = BVCopy(X,Z);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(Z,ctx->bs,2*ctx->bs-ini);CHKERRQ(ierr);
    ierr = BVCopy(R,Z);CHKERRQ(ierr);
    if (its>1) {
      ierr = BVSetActiveColumns(Z,2*ctx->bs-ini,3*ctx->bs-2*ini);CHKERRQ(ierr);
      ierr = BVCopy(P,Z);CHKERRQ(ierr);
    }

    if (its>1) nv = 3*ctx->bs-2*ini;
    else nv = 2*ctx->bs-ini;

    ierr = BVSetActiveColumns(Z,0,nv);CHKERRQ(ierr);
    ierr = DSSetDimensions(eps->ds,nv,0,0,0);CHKERRQ(ierr);
    ierr = DSGetMat(eps->ds,DS_MAT_A,&M);CHKERRQ(ierr);
    ierr = BVMatProject(Z,A,Z,M);CHKERRQ(ierr);
    if (flip) { ierr = MatScale(M,-1.0);CHKERRQ(ierr); }
    ierr = DSRestoreMat(eps->ds,DS_MAT_A,&M);CHKERRQ(ierr);
    ierr = DSGetMat(eps->ds,DS_MAT_B,&M);CHKERRQ(ierr);
    ierr = BVMatProject(Z,B,Z,M);CHKERRQ(ierr); /* covers also the case B=NULL */
    ierr = DSRestoreMat(eps->ds,DS_MAT_B,&M);CHKERRQ(ierr);

    /* 24. Solve the generalized eigenvalue problem */
    ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
    ierr = DSSolve(eps->ds,eigr,NULL);CHKERRQ(ierr);
    ierr = DSSort(eps->ds,eigr,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DSSynchronize(eps->ds,eigr,NULL);CHKERRQ(ierr);
    for (j=0;j<nv;j++) if (locked+j<eps->ncv) eps->eigr[locked+j] = flip? -eigr[j]: eigr[j];
    ierr = DSVectors(eps->ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);

    /* 25-33. Compute Ritz vectors */
    ierr = DSGetMat(eps->ds,DS_MAT_X,&M);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(Z,ctx->bs,nv);CHKERRQ(ierr);
    if (ctx->lock) {
      ierr = BVSetActiveColumns(P,0,ctx->bs);CHKERRQ(ierr);
    }
    ierr = BVMult(P,1.0,0.0,Z,M);CHKERRQ(ierr);
    ierr = BVCopy(P,X);CHKERRQ(ierr);
    if (ctx->lock) {
      ierr = BVSetActiveColumns(P,nconv,ctx->bs);CHKERRQ(ierr);
    }
    ierr = BVSetActiveColumns(Z,0,ctx->bs);CHKERRQ(ierr);
    ierr = BVMult(X,1.0,1.0,Z,M);CHKERRQ(ierr);
    if (ctx->lock) {
      ierr = BVSetActiveColumns(X,nconv,ctx->bs);CHKERRQ(ierr);
    }
    ierr = BVMatMult(X,A,AX);CHKERRQ(ierr);
    ierr = MatDestroy(&M);CHKERRQ(ierr);
  }

diverged:
  eps->its += its;

  if (flip) sc->comparison = SlepcCompareLargestReal;
  ierr = PetscFree(eigr);CHKERRQ(ierr);
  ierr = BVDestroy(&Z);CHKERRQ(ierr);
  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = BVDestroy(&R);CHKERRQ(ierr);
  ierr = BVDestroy(&P);CHKERRQ(ierr);
  ierr = BVDestroy(&AX);CHKERRQ(ierr);
  if (B) {
    ierr = BVDestroy(&BX);CHKERRQ(ierr);
  }
  if (nc>0 || eps->nev>ctx->bs-ctx->guard) {
    ierr = BVDestroy(&Y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLOBPCGSetBlockSize_LOBPCG(EPS eps,PetscInt bs)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  if (bs == PETSC_DEFAULT || bs == PETSC_DECIDE) bs = 1;
  if (bs <= 0) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid block size %D",bs);
  if (ctx->bs != bs) {
    ctx->bs = bs;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSLOBPCGSetBlockSize - Sets the block size of the LOBPCG method.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  bs  - the block size

   Options Database Key:
.  -eps_lobpcg_blocksize - Sets the block size

   Level: advanced

.seealso: EPSLOBPCGGetBlockSize()
@*/
PetscErrorCode EPSLOBPCGSetBlockSize(EPS eps,PetscInt bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,bs,2);
  ierr = PetscTryMethod(eps,"EPSLOBPCGSetBlockSize_C",(EPS,PetscInt),(eps,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLOBPCGGetBlockSize_LOBPCG(EPS eps,PetscInt *bs)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  *bs = ctx->bs;
  PetscFunctionReturn(0);
}

/*@
   EPSLOBPCGGetBlockSize - Gets the block size used in the LOBPCG method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  bs - the block size

   Level: advanced

.seealso: EPSLOBPCGSetBlockSize()
@*/
PetscErrorCode EPSLOBPCGGetBlockSize(EPS eps,PetscInt *bs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(bs,2);
  ierr = PetscUseMethod(eps,"EPSLOBPCGGetBlockSize_C",(EPS,PetscInt*),(eps,bs));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLOBPCGSetRestart_LOBPCG(EPS eps,PetscReal restart)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  if (restart==PETSC_DEFAULT) restart = 0.9;
  if (restart<0.1 || restart>1.0) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The restart argument %g must be in the range [0.1,1.0]",(double)restart);
  if (restart != ctx->restart) {
    ctx->restart = restart;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSLOBPCGSetRestart - Sets the restart parameter for the LOBPCG method.
   The meaning of this parameter is the proportion of vectors within the
   current block iterate that must have converged in order to force a
   restart with hard locking.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  restart - the percentage of the block of vectors to force a restart

   Options Database Key:
.  -eps_lobpcg_restart - Sets the restart parameter

   Notes:
   Allowed values are in the range [0.1,1.0]. The default is 0.6.

   Level: advanced

.seealso: EPSLOBPCGGetRestart()
@*/
PetscErrorCode EPSLOBPCGSetRestart(EPS eps,PetscReal restart)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,restart,2);
  ierr = PetscTryMethod(eps,"EPSLOBPCGSetRestart_C",(EPS,PetscReal),(eps,restart));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLOBPCGGetRestart_LOBPCG(EPS eps,PetscReal *restart)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  *restart = ctx->restart;
  PetscFunctionReturn(0);
}

/*@
   EPSLOBPCGGetRestart - Gets the restart parameter used in the LOBPCG method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  restart - the restart parameter

   Level: advanced

.seealso: EPSLOBPCGSetRestart()
@*/
PetscErrorCode EPSLOBPCGGetRestart(EPS eps,PetscReal *restart)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidRealPointer(restart,2);
  ierr = PetscUseMethod(eps,"EPSLOBPCGGetRestart_C",(EPS,PetscReal*),(eps,restart));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLOBPCGSetLocking_LOBPCG(EPS eps,PetscBool lock)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  ctx->lock = lock;
  PetscFunctionReturn(0);
}

/*@
   EPSLOBPCGSetLocking - Choose between locking and non-locking variants of
   the LOBPCG method.

   Logically Collective on eps

   Input Parameters:
+  eps  - the eigenproblem solver context
-  lock - true if the locking variant must be selected

   Options Database Key:
.  -eps_lobpcg_locking - Sets the locking flag

   Notes:
   This flag refers to soft locking (converged vectors within the current
   block iterate), since hard locking is always used (when nev is larger
   than the block size).

   Level: advanced

.seealso: EPSLOBPCGGetLocking()
@*/
PetscErrorCode EPSLOBPCGSetLocking(EPS eps,PetscBool lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,lock,2);
  ierr = PetscTryMethod(eps,"EPSLOBPCGSetLocking_C",(EPS,PetscBool),(eps,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLOBPCGGetLocking_LOBPCG(EPS eps,PetscBool *lock)
{
  EPS_LOBPCG *ctx = (EPS_LOBPCG*)eps->data;

  PetscFunctionBegin;
  *lock = ctx->lock;
  PetscFunctionReturn(0);
}

/*@
   EPSLOBPCGGetLocking - Gets the locking flag used in the LOBPCG method.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  lock - the locking flag

   Level: advanced

.seealso: EPSLOBPCGSetLocking()
@*/
PetscErrorCode EPSLOBPCGGetLocking(EPS eps,PetscBool *lock)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(lock,2);
  ierr = PetscUseMethod(eps,"EPSLOBPCGGetLocking_C",(EPS,PetscBool*),(eps,lock));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_LOBPCG(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  EPS_LOBPCG     *ctx = (EPS_LOBPCG*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  block size %D\n",ctx->bs);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  restart parameter=%g (using %D guard vectors)\n",(double)ctx->restart,ctx->guard);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  soft locking %sactivated\n",ctx->lock?"":"de");CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_LOBPCG(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      lock,flg;
  PetscInt       bs;
  PetscReal      restart;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS LOBPCG Options");CHKERRQ(ierr);

    ierr = PetscOptionsInt("-eps_lobpcg_blocksize","Block size","EPSLOBPCGSetBlockSize",20,&bs,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSLOBPCGSetBlockSize(eps,bs);CHKERRQ(ierr); }

    ierr = PetscOptionsReal("-eps_lobpcg_restart","Percentage of the block of vectors to force a restart","EPSLOBPCGSetRestart",0.5,&restart,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSLOBPCGSetRestart(eps,restart);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-eps_lobpcg_locking","Choose between locking and non-locking variants","EPSLOBPCGSetLocking",PETSC_TRUE,&lock,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSLOBPCGSetLocking(eps,lock);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_LOBPCG(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetLocking_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetLocking_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_LOBPCG(EPS eps)
{
  EPS_LOBPCG     *lobpcg;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&lobpcg);CHKERRQ(ierr);
  eps->data = (void*)lobpcg;
  lobpcg->lock = PETSC_TRUE;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_LOBPCG;
  eps->ops->setup          = EPSSetUp_LOBPCG;
  eps->ops->setfromoptions = EPSSetFromOptions_LOBPCG;
  eps->ops->destroy        = EPSDestroy_LOBPCG;
  eps->ops->view           = EPSView_LOBPCG;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->setdefaultst   = EPSSetDefaultST_GMRES;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetBlockSize_C",EPSLOBPCGSetBlockSize_LOBPCG);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetBlockSize_C",EPSLOBPCGGetBlockSize_LOBPCG);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetRestart_C",EPSLOBPCGSetRestart_LOBPCG);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetRestart_C",EPSLOBPCGGetRestart_LOBPCG);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGSetLocking_C",EPSLOBPCGSetLocking_LOBPCG);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLOBPCGGetLocking_C",EPSLOBPCGGetLocking_LOBPCG);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

