/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Newton refinement for NEP, simple version
*/

#include <slepc/private/nepimpl.h>
#include <slepcblaslapack.h>

#define NREF_MAXIT 10

typedef struct {
  VecScatter    *scatter_id,nst;
  Mat           *A;
  Vec           nv,vg,v,w;
  FN            *fn;
} NEPSimpNRefctx;

typedef struct {
  Mat          M1;
  Vec          M2,M3;
  PetscScalar  M4,m3;
} NEP_REFINE_MATSHELL;

static PetscErrorCode MatMult_FS(Mat M ,Vec x,Vec y)
{
  PetscErrorCode      ierr;
  NEP_REFINE_MATSHELL *ctx;
  PetscScalar         t;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,&ctx);CHKERRQ(ierr);
  ierr = VecDot(x,ctx->M3,&t);CHKERRQ(ierr);
  t *= ctx->m3/ctx->M4;
  ierr = MatMult(ctx->M1,x,y);CHKERRQ(ierr);
  ierr = VecAXPY(y,-t,ctx->M2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSimpleNRefSetUp(NEP nep,NEPSimpNRefctx **ctx_)
{
  PetscErrorCode ierr;
  PetscInt       i,si,j,n0,m0,nloc,*idx1,*idx2,ne;
  IS             is1,is2;
  NEPSimpNRefctx *ctx;
  Vec            v;
  PetscMPIInt    rank,size;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,ctx_);CHKERRQ(ierr);
  ctx = *ctx_;
  if (nep->npart==1) {
    ctx->A  = nep->A;
    ctx->fn = nep->f;
    nep->refinesubc = NULL;
    ctx->scatter_id = NULL;
  } else {
    ierr = PetscMalloc2(nep->nt,&ctx->A,nep->npart,&ctx->scatter_id);CHKERRQ(ierr);

    /* Duplicate matrices */
    for (i=0;i<nep->nt;i++) {
      ierr = MatCreateRedundantMatrix(nep->A[i],0,PetscSubcommChild(nep->refinesubc),MAT_INITIAL_MATRIX,&ctx->A[i]);CHKERRQ(ierr);
    }
    ierr = MatCreateVecs(ctx->A[0],&ctx->v,NULL);CHKERRQ(ierr);

    /* Duplicate FNs */
    ierr = PetscMalloc1(nep->nt,&ctx->fn);CHKERRQ(ierr);
    for (i=0;i<nep->nt;i++) {
      ierr = FNDuplicate(nep->f[i],PetscSubcommChild(nep->refinesubc),&ctx->fn[i]);CHKERRQ(ierr);
    }

    /* Create scatters for sending vectors to each subcommucator */
    ierr = BVGetColumn(nep->V,0,&v);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(v,&n0,&m0);CHKERRQ(ierr);
    ierr = BVRestoreColumn(nep->V,0,&v);CHKERRQ(ierr);
    ierr = VecGetLocalSize(ctx->v,&nloc);CHKERRQ(ierr);
    ierr = PetscMalloc2(m0-n0,&idx1,m0-n0,&idx2);CHKERRQ(ierr);
    ierr = VecCreateMPI(PetscObjectComm((PetscObject)nep),nloc,PETSC_DECIDE,&ctx->vg);CHKERRQ(ierr);
    for (si=0;si<nep->npart;si++) {
      j = 0;
      for (i=n0;i<m0;i++) {
        idx1[j]   = i;
        idx2[j++] = i+nep->n*si;
      }
      ierr = ISCreateGeneral(PetscObjectComm((PetscObject)nep),(m0-n0),idx1,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
      ierr = ISCreateGeneral(PetscObjectComm((PetscObject)nep),(m0-n0),idx2,PETSC_COPY_VALUES,&is2);CHKERRQ(ierr);
      ierr = BVGetColumn(nep->V,0,&v);CHKERRQ(ierr);
      ierr = VecScatterCreate(v,is1,ctx->vg,is2,&ctx->scatter_id[si]);CHKERRQ(ierr);
      ierr = BVRestoreColumn(nep->V,0,&v);CHKERRQ(ierr);
      ierr = ISDestroy(&is1);CHKERRQ(ierr);
      ierr = ISDestroy(&is2);CHKERRQ(ierr);
    }
    ierr = PetscFree2(idx1,idx2);CHKERRQ(ierr);
  }
  if (nep->scheme==NEP_REFINE_SCHEME_EXPLICIT) {
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)ctx->A[0]),&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)ctx->A[0]),&size);CHKERRQ(ierr);
    if (size>1) {
      if (nep->npart==1) {
        ierr = BVGetColumn(nep->V,0,&v);CHKERRQ(ierr);
      } else v = ctx->v;
      ierr = VecGetOwnershipRange(v,&n0,&m0);CHKERRQ(ierr);
      ne = (rank == size-1)?nep->n:0;
      ierr = VecCreateMPI(PetscObjectComm((PetscObject)ctx->A[0]),ne,PETSC_DECIDE,&ctx->nv);CHKERRQ(ierr);
      ierr = PetscMalloc1(m0-n0,&idx1);CHKERRQ(ierr);
      for (i=n0;i<m0;i++) {
        idx1[i-n0] = i;
      }
      ierr = ISCreateGeneral(PetscObjectComm((PetscObject)ctx->A[0]),(m0-n0),idx1,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
      ierr = VecScatterCreate(v,is1,ctx->nv,is1,&ctx->nst);CHKERRQ(ierr);
      if (nep->npart==1) {
        ierr = BVRestoreColumn(nep->V,0,&v);CHKERRQ(ierr);
      }
      ierr = PetscFree(idx1);CHKERRQ(ierr);
      ierr = ISDestroy(&is1);CHKERRQ(ierr);
    }
  }  PetscFunctionReturn(0);
}

/*
  Gather Eigenpair idx from subcommunicator with color sc
*/
static PetscErrorCode NEPSimpleNRefGatherEigenpair(NEP nep,NEPSimpNRefctx *ctx,PetscInt sc,PetscInt idx,PetscInt *fail)
{
  PetscErrorCode ierr;
  PetscMPIInt    nproc,p;
  MPI_Comm       comm=((PetscObject)nep)->comm;
  Vec            v;
  PetscScalar    *array;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);
  p = (nproc/nep->npart)*(sc+1)+PetscMin(sc+1,nproc%nep->npart)-1;
  if (nep->npart>1) {
    /* Communicate convergence successful */
    ierr = MPI_Bcast(fail,1,MPIU_INT,p,comm);CHKERRQ(ierr);
    if (!(*fail)) {
      /* Process 0 of subcommunicator sc broadcasts the eigenvalue */
      ierr = MPI_Bcast(&nep->eigr[idx],1,MPIU_SCALAR,p,comm);CHKERRQ(ierr);
      /* Gather nep->V[idx] from the subcommuniator sc */
      ierr = BVGetColumn(nep->V,idx,&v);CHKERRQ(ierr);
      if (nep->refinesubc->color==sc) {
        ierr = VecGetArray(ctx->v,&array);CHKERRQ(ierr);
        ierr = VecPlaceArray(ctx->vg,array);CHKERRQ(ierr);
      }
      ierr = VecScatterBegin(ctx->scatter_id[sc],ctx->vg,v,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx->scatter_id[sc],ctx->vg,v,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
      if (nep->refinesubc->color==sc) {
        ierr = VecResetArray(ctx->vg);CHKERRQ(ierr);
        ierr = VecRestoreArray(ctx->v,&array);CHKERRQ(ierr);
      }
      ierr = BVRestoreColumn(nep->V,idx,&v);CHKERRQ(ierr);
    }
  } else {
    if (nep->scheme==NEP_REFINE_SCHEME_EXPLICIT && !(*fail)) {
      ierr = MPI_Bcast(&nep->eigr[idx],1,MPIU_SCALAR,p,comm);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSimpleNRefScatterEigenvector(NEP nep,NEPSimpNRefctx *ctx,PetscInt sc,PetscInt idx)
{
  PetscErrorCode ierr;
  Vec            v;
  PetscScalar    *array;

  PetscFunctionBegin;
  if (nep->npart>1) {
    ierr = BVGetColumn(nep->V,idx,&v);CHKERRQ(ierr);
    if (nep->refinesubc->color==sc) {
      ierr = VecGetArray(ctx->v,&array);CHKERRQ(ierr);
      ierr = VecPlaceArray(ctx->vg,array);CHKERRQ(ierr);
    }
    ierr = VecScatterBegin(ctx->scatter_id[sc],v,ctx->vg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatter_id[sc],v,ctx->vg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    if (nep->refinesubc->color==sc) {
      ierr = VecResetArray(ctx->vg);CHKERRQ(ierr);
      ierr = VecRestoreArray(ctx->v,&array);CHKERRQ(ierr);
    }
    ierr = BVRestoreColumn(nep->V,idx,&v);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPSimpleNRefSetUpSystem(NEP nep,NEPSimpNRefctx *ctx,Mat *A,PetscInt idx,Mat *Mt,Mat *T,Mat *P,PetscBool ini,Vec t,Vec v)
{
  PetscErrorCode      ierr;
  PetscInt            i,st,ml,m0,n0,m1,mg;
  PetscInt            *dnz,*onz,ncols,*cols2=NULL,*nnz,nt=nep->nt;
  PetscScalar         zero=0.0,*coeffs,*coeffs2;
  PetscMPIInt         rank,size;
  MPI_Comm            comm;
  const PetscInt      *cols;
  const PetscScalar   *vals,*array;
  NEP_REFINE_MATSHELL *fctx;
  Vec                 w=ctx->w;
  Mat                 M;

  PetscFunctionBegin;
  ierr = PetscMalloc2(nt,&coeffs,nt,&coeffs2);CHKERRQ(ierr);
  switch (nep->scheme) {
  case NEP_REFINE_SCHEME_SCHUR:
    if (ini) {
      ierr = PetscCalloc1(1,&fctx);CHKERRQ(ierr);
      ierr = MatGetSize(A[0],&m0,&n0);CHKERRQ(ierr);
      ierr = MatCreateShell(PetscObjectComm((PetscObject)A[0]),PETSC_DECIDE,PETSC_DECIDE,m0,n0,fctx,T);CHKERRQ(ierr);
      ierr = MatShellSetOperation(*T,MATOP_MULT,(void(*)(void))MatMult_FS);CHKERRQ(ierr);
    } else {
      ierr = MatShellGetContext(*T,&fctx);CHKERRQ(ierr);
    }
    M=fctx->M1;
    break;
  case NEP_REFINE_SCHEME_MBE:
    M=*T;
    break;
  case NEP_REFINE_SCHEME_EXPLICIT:
    M=*Mt;
    break;
  }
  if (ini) {
    ierr = MatDuplicate(A[0],MAT_COPY_VALUES,&M);CHKERRQ(ierr);
  } else {
    ierr = MatCopy(A[0],M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  }
  for (i=0;i<nt;i++) {
    ierr = FNEvaluateFunction(ctx->fn[i],nep->eigr[idx],coeffs+i);CHKERRQ(ierr);
  }
  if (coeffs[0]!=1.0) {
    ierr = MatScale(M,coeffs[0]);CHKERRQ(ierr);
  }
  for (i=1;i<nt;i++) {
    ierr = MatAXPY(M,coeffs[i],A[i],(ini)?nep->mstr:SUBSET_NONZERO_PATTERN);CHKERRQ(ierr);
  }
  for (i=0;i<nt;i++) {
    ierr = FNEvaluateDerivative(ctx->fn[i],nep->eigr[idx],coeffs2+i);CHKERRQ(ierr);
  }
  st = 0;
  for (i=0;i<nt && PetscAbsScalar(coeffs2[i])==0.0;i++) st++;
  ierr = MatMult(A[st],v,w);CHKERRQ(ierr);
  if (coeffs2[st]!=1.0) {
    ierr = VecScale(w,coeffs2[st]);CHKERRQ(ierr);
  }
  for (i=st+1;i<nt;i++) {
    ierr = MatMult(A[i],v,t);CHKERRQ(ierr);
    ierr = VecAXPY(w,coeffs2[i],t);CHKERRQ(ierr);
  }

  switch (nep->scheme) {
  case NEP_REFINE_SCHEME_EXPLICIT:
    comm = PetscObjectComm((PetscObject)A[0]);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
    ierr = MatGetSize(M,&mg,NULL);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(M,&m0,&m1);CHKERRQ(ierr);
    if (ini) {
      ierr = MatCreate(comm,T);CHKERRQ(ierr);
      ierr = MatGetLocalSize(M,&ml,NULL);CHKERRQ(ierr);
      if (rank==size-1) ml++;
      ierr = MatSetSizes(*T,ml,ml,mg+1,mg+1);CHKERRQ(ierr);
      ierr = MatSetFromOptions(*T);CHKERRQ(ierr);
      ierr = MatSetUp(*T);CHKERRQ(ierr);
      /* Preallocate M */
      if (size>1) {
        ierr = MatPreallocateInitialize(comm,ml,ml,dnz,onz);CHKERRQ(ierr);
        for (i=m0;i<m1;i++) {
          ierr = MatGetRow(M,i,&ncols,&cols,NULL);CHKERRQ(ierr);
          ierr = MatPreallocateSet(i,ncols,cols,dnz,onz);CHKERRQ(ierr);
          ierr = MatPreallocateSet(i,1,&mg,dnz,onz);CHKERRQ(ierr);
          ierr = MatRestoreRow(M,i,&ncols,&cols,NULL);CHKERRQ(ierr);
        }
        if (rank==size-1) {
          ierr = PetscCalloc1(mg+1,&cols2);CHKERRQ(ierr);
          for (i=0;i<mg+1;i++) cols2[i]=i;
          ierr = MatPreallocateSet(m1,mg+1,cols2,dnz,onz);CHKERRQ(ierr);
          ierr = PetscFree(cols2);CHKERRQ(ierr);
        }
        ierr = MatMPIAIJSetPreallocation(*T,0,dnz,0,onz);CHKERRQ(ierr);
        ierr = MatPreallocateFinalize(dnz,onz);CHKERRQ(ierr);
      } else {
        ierr = PetscCalloc1(mg+1,&nnz);CHKERRQ(ierr);
        for (i=0;i<mg;i++) {
          ierr = MatGetRow(M,i,&ncols,NULL,NULL);CHKERRQ(ierr);
          nnz[i] = ncols+1;
          ierr = MatRestoreRow(M,i,&ncols,NULL,NULL);CHKERRQ(ierr);
        }
        nnz[mg] = mg+1;
        ierr = MatSeqAIJSetPreallocation(*T,0,nnz);CHKERRQ(ierr);
        ierr = PetscFree(nnz);CHKERRQ(ierr);
      }
      *Mt = M;
      *P  = *T;
    }

    /* Set values */
    ierr = VecGetArrayRead(w,&array);CHKERRQ(ierr);
    for (i=m0;i<m1;i++) {
      ierr = MatGetRow(M,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      ierr = MatSetValues(*T,1,&i,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatRestoreRow(M,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      ierr = MatSetValues(*T,1,&i,1,&mg,array+i-m0,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(w,&array);CHKERRQ(ierr);
    ierr = VecConjugate(v);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)A[0]),&size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)A[0]),&rank);CHKERRQ(ierr);
    if (size>1) {
      if (rank==size-1) {
        ierr = PetscMalloc1(nep->n,&cols2);CHKERRQ(ierr);
        for (i=0;i<nep->n;i++) cols2[i]=i;
      }
      ierr = VecScatterBegin(ctx->nst,v,ctx->nv,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx->nst,v,ctx->nv,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = VecGetArrayRead(ctx->nv,&array);CHKERRQ(ierr);
      if (rank==size-1) {
        ierr = MatSetValues(*T,1,&mg,nep->n,cols2,array,INSERT_VALUES);CHKERRQ(ierr);
        ierr = MatSetValues(*T,1,&mg,1,&mg,&zero,INSERT_VALUES);CHKERRQ(ierr);
      }
      ierr = VecRestoreArrayRead(ctx->nv,&array);CHKERRQ(ierr);
    } else {
      ierr = PetscMalloc1(m1-m0,&cols2);CHKERRQ(ierr);
      for (i=0;i<m1-m0;i++) cols2[i]=m0+i;
      ierr = VecGetArrayRead(v,&array);CHKERRQ(ierr);
      ierr = MatSetValues(*T,1,&mg,m1-m0,cols2,array,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatSetValues(*T,1,&mg,1,&mg,&zero,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(v,&array);CHKERRQ(ierr);
    }
    ierr = VecConjugate(v);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(*T,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*T,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = PetscFree(cols2);CHKERRQ(ierr);
    break;
  case NEP_REFINE_SCHEME_SCHUR:
    fctx->M2 = ctx->w;
    fctx->M3 = v;
    fctx->m3 = 1.0+PetscConj(nep->eigr[idx])*nep->eigr[idx];
    fctx->M4 = PetscConj(nep->eigr[idx]);
    fctx->M1 = M;
    if (ini) {
      ierr = MatDuplicate(M,MAT_COPY_VALUES,P);CHKERRQ(ierr);
    } else {
      ierr = MatCopy(M,*P,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
    }
    if (fctx->M4!=0.0) {
      ierr = VecConjugate(v);CHKERRQ(ierr);
      ierr = VecPointwiseMult(t,v,w);CHKERRQ(ierr);
      ierr = VecConjugate(v);CHKERRQ(ierr);
      ierr = VecScale(t,-fctx->m3/fctx->M4);CHKERRQ(ierr);
      ierr = MatDiagonalSet(*P,t,ADD_VALUES);CHKERRQ(ierr);
    }
    break;
  case NEP_REFINE_SCHEME_MBE:
    *T = M;
    *P = M;
    break;
  }
  ierr = PetscFree2(coeffs,coeffs2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPNewtonRefinementSimple(NEP nep,PetscInt *maxits,PetscReal tol,PetscInt k)
{
  PetscErrorCode      ierr;
  PetscInt            i,n,its,idx=0,*idx_sc,*its_sc,color,*fail_sc;
  PetscMPIInt         rank,size;
  Mat                 Mt=NULL,T=NULL,P=NULL;
  MPI_Comm            comm;
  Vec                 r,v,dv,rr=NULL,dvv=NULL,t[2];
  const PetscScalar   *array;
  PetscScalar         *array2,deig=0.0,tt[2],ttt;
  PetscReal           norm,error;
  PetscBool           ini=PETSC_TRUE,sc_pend,solved=PETSC_FALSE;
  NEPSimpNRefctx      *ctx;
  NEP_REFINE_MATSHELL *fctx=NULL;
  KSPConvergedReason  reason;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(NEP_Refine,nep,0,0,0);CHKERRQ(ierr);
  ierr = NEPSimpleNRefSetUp(nep,&ctx);CHKERRQ(ierr);
  its = (maxits)?*maxits:NREF_MAXIT;
  comm = (nep->npart==1)?PetscObjectComm((PetscObject)nep):PetscSubcommChild(nep->refinesubc);
  if (!nep->refineksp) { ierr = NEPRefineGetKSP(nep,&nep->refineksp);CHKERRQ(ierr); }
  if (nep->npart==1) {
    ierr = BVGetColumn(nep->V,0,&v);CHKERRQ(ierr);
  } else v = ctx->v;
  ierr = VecDuplicate(v,&ctx->w);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&dv);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&t[0]);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&t[1]);CHKERRQ(ierr);
  if (nep->npart==1) { ierr = BVRestoreColumn(nep->V,0,&v);CHKERRQ(ierr); }
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = VecGetLocalSize(r,&n);CHKERRQ(ierr);
  ierr = PetscMalloc3(nep->npart,&idx_sc,nep->npart,&its_sc,nep->npart,&fail_sc);CHKERRQ(ierr);
  for (i=0;i<nep->npart;i++) fail_sc[i] = 0;
  for (i=0;i<nep->npart;i++) its_sc[i] = 0;
  color = (nep->npart==1)?0:nep->refinesubc->color;

  /* Loop performing iterative refinements */
  while (!solved) {
    for (i=0;i<nep->npart;i++) {
      sc_pend = PETSC_TRUE;
      if (its_sc[i]==0) {
        idx_sc[i] = idx++;
        if (idx_sc[i]>=k) {
          sc_pend = PETSC_FALSE;
        } else {
          ierr = NEPSimpleNRefScatterEigenvector(nep,ctx,i,idx_sc[i]);CHKERRQ(ierr);
        }
      }  else { /* Gather Eigenpair from subcommunicator i */
        ierr = NEPSimpleNRefGatherEigenpair(nep,ctx,i,idx_sc[i],&fail_sc[i]);CHKERRQ(ierr);
      }
      while (sc_pend) {
        if (!fail_sc[i]) {
          ierr = NEPComputeError(nep,idx_sc[i],NEP_ERROR_RELATIVE,&error);CHKERRQ(ierr);
        }
        if (error<=tol || its_sc[i]>=its || fail_sc[i]) {
          idx_sc[i] = idx++;
          its_sc[i] = 0;
          fail_sc[i] = 0;
          if (idx_sc[i]<k) { ierr = NEPSimpleNRefScatterEigenvector(nep,ctx,i,idx_sc[i]);CHKERRQ(ierr); }
        } else {
          sc_pend = PETSC_FALSE;
          its_sc[i]++;
        }
        if (idx_sc[i]>=k) sc_pend = PETSC_FALSE;
      }
    }
    solved = PETSC_TRUE;
    for (i=0;i<nep->npart&&solved;i++) solved = PetscNot(idx_sc[i]<k);
    if (idx_sc[color]<k) {
#if !defined(PETSC_USE_COMPLEX)
      if (nep->eigi[idx_sc[color]]!=0.0) SETERRQ(PetscObjectComm((PetscObject)nep),1,"Simple Refinement not implemented in real scalar for complex eigenvalues");
#endif
      if (nep->npart==1) {
        ierr = BVGetColumn(nep->V,idx_sc[color],&v);CHKERRQ(ierr);
      } else v = ctx->v;
      ierr = NEPSimpleNRefSetUpSystem(nep,ctx,ctx->A,idx_sc[color],&Mt,&T,&P,ini,t[0],v);CHKERRQ(ierr);
      ierr = KSPSetOperators(nep->refineksp,T,P);CHKERRQ(ierr);
      if (ini) {
        ierr = KSPSetFromOptions(nep->refineksp);CHKERRQ(ierr);
        if (nep->scheme==NEP_REFINE_SCHEME_EXPLICIT) {
          ierr = MatCreateVecs(T,&dvv,NULL);CHKERRQ(ierr);
          ierr = VecDuplicate(dvv,&rr);CHKERRQ(ierr);
        }
        ini = PETSC_FALSE;
      }
      switch (nep->scheme) {
      case NEP_REFINE_SCHEME_EXPLICIT:
        ierr = MatMult(Mt,v,r);CHKERRQ(ierr);
        ierr = VecGetArrayRead(r,&array);CHKERRQ(ierr);
        if (rank==size-1) {
          ierr = VecGetArray(rr,&array2);CHKERRQ(ierr);
          ierr = PetscArraycpy(array2,array,n);CHKERRQ(ierr);
          array2[n] = 0.0;
          ierr = VecRestoreArray(rr,&array2);CHKERRQ(ierr);
        } else {
          ierr = VecPlaceArray(rr,array);CHKERRQ(ierr);
        }
        ierr = KSPSolve(nep->refineksp,rr,dvv);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(nep->refineksp,&reason);CHKERRQ(ierr);
        if (reason>0) {
          if (rank != size-1) {
            ierr = VecResetArray(rr);CHKERRQ(ierr);
          }
          ierr = VecRestoreArrayRead(r,&array);CHKERRQ(ierr);
          ierr = VecGetArrayRead(dvv,&array);CHKERRQ(ierr);
          ierr = VecPlaceArray(dv,array);CHKERRQ(ierr);
          ierr = VecAXPY(v,-1.0,dv);CHKERRQ(ierr);
          ierr = VecNorm(v,NORM_2,&norm);CHKERRQ(ierr);
          ierr = VecScale(v,1.0/norm);CHKERRQ(ierr);
          ierr = VecResetArray(dv);CHKERRQ(ierr);
          if (rank==size-1) nep->eigr[idx_sc[color]] -= array[n];
          ierr = VecRestoreArrayRead(dvv,&array);CHKERRQ(ierr);
        } else fail_sc[color] = 1;
        break;
      case NEP_REFINE_SCHEME_MBE:
        ierr = MatMult(T,v,r);CHKERRQ(ierr);
        /* Mixed block elimination */
        ierr = VecConjugate(v);CHKERRQ(ierr);
        ierr = KSPSolveTranspose(nep->refineksp,v,t[0]);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(nep->refineksp,&reason);CHKERRQ(ierr);
        if (reason>0) {
          ierr = VecConjugate(t[0]);CHKERRQ(ierr);
          ierr = VecDot(ctx->w,t[0],&tt[0]);CHKERRQ(ierr);
          ierr = KSPSolve(nep->refineksp,ctx->w,t[1]);CHKERRQ(ierr);
          ierr = KSPGetConvergedReason(nep->refineksp,&reason);CHKERRQ(ierr);
          if (reason>0) {
            ierr = VecDot(t[1],v,&tt[1]);CHKERRQ(ierr);
            ierr = VecDot(r,t[0],&ttt);CHKERRQ(ierr);
            tt[0] = ttt/tt[0];
            ierr = VecAXPY(r,-tt[0],ctx->w);CHKERRQ(ierr);
            ierr = KSPSolve(nep->refineksp,r,dv);CHKERRQ(ierr);
            ierr = KSPGetConvergedReason(nep->refineksp,&reason);CHKERRQ(ierr);
            if (reason>0) {
              ierr = VecDot(dv,v,&ttt);CHKERRQ(ierr);
              tt[1] = ttt/tt[1];
              ierr = VecAXPY(dv,-tt[1],t[1]);CHKERRQ(ierr);
              deig = tt[0]+tt[1];
            }
          }
          ierr = VecConjugate(v);CHKERRQ(ierr);
          ierr = VecAXPY(v,-1.0,dv);CHKERRQ(ierr);
          ierr = VecNorm(v,NORM_2,&norm);CHKERRQ(ierr);
          ierr = VecScale(v,1.0/norm);CHKERRQ(ierr);
          nep->eigr[idx_sc[color]] -= deig;
          fail_sc[color] = 0;
        } else {
          ierr = VecConjugate(v);CHKERRQ(ierr);
          fail_sc[color] = 1;
        }
        break;
      case NEP_REFINE_SCHEME_SCHUR:
        fail_sc[color] = 1;
        ierr = MatShellGetContext(T,&fctx);CHKERRQ(ierr);
        if (fctx->M4!=0.0) {
          ierr = MatMult(fctx->M1,v,r);CHKERRQ(ierr);
          ierr = KSPSolve(nep->refineksp,r,dv);CHKERRQ(ierr);
          ierr = KSPGetConvergedReason(nep->refineksp,&reason);CHKERRQ(ierr);
          if (reason>0) {
            ierr = VecDot(dv,v,&deig);CHKERRQ(ierr);
            deig *= -fctx->m3/fctx->M4;
            ierr = VecAXPY(v,-1.0,dv);CHKERRQ(ierr);
            ierr = VecNorm(v,NORM_2,&norm);CHKERRQ(ierr);
            ierr = VecScale(v,1.0/norm);CHKERRQ(ierr);
            nep->eigr[idx_sc[color]] -= deig;
            fail_sc[color] = 0;
          }
        }
        break;
      }
      if (nep->npart==1) { ierr = BVRestoreColumn(nep->V,idx_sc[color],&v);CHKERRQ(ierr); }
    }
  }
  ierr = VecDestroy(&t[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&t[1]);CHKERRQ(ierr);
  ierr = VecDestroy(&dv);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->w);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = PetscFree3(idx_sc,its_sc,fail_sc);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&ctx->nst);CHKERRQ(ierr);
  if (nep->npart>1) {
    ierr = VecDestroy(&ctx->vg);CHKERRQ(ierr);
    ierr = VecDestroy(&ctx->v);CHKERRQ(ierr);
    for (i=0;i<nep->nt;i++) {
      ierr = MatDestroy(&ctx->A[i]);CHKERRQ(ierr);
    }
    for (i=0;i<nep->npart;i++) {
      ierr = VecScatterDestroy(&ctx->scatter_id[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree2(ctx->A,ctx->scatter_id);CHKERRQ(ierr);
  }
  if (fctx && nep->scheme==NEP_REFINE_SCHEME_SCHUR) {
    ierr = MatDestroy(&P);CHKERRQ(ierr);
    ierr = MatDestroy(&fctx->M1);CHKERRQ(ierr);
    ierr = PetscFree(fctx);CHKERRQ(ierr);
  }
  if (nep->scheme==NEP_REFINE_SCHEME_EXPLICIT) {
    ierr = MatDestroy(&Mt);CHKERRQ(ierr);
    ierr = VecDestroy(&dvv);CHKERRQ(ierr);
    ierr = VecDestroy(&rr);CHKERRQ(ierr);
    ierr = VecDestroy(&ctx->nv);CHKERRQ(ierr);
    if (nep->npart>1) {
      for (i=0;i<nep->nt;i++) { ierr = FNDestroy(&ctx->fn[i]);CHKERRQ(ierr); }
      ierr = PetscFree(ctx->fn);CHKERRQ(ierr);
    }
  }
  if (nep->scheme==NEP_REFINE_SCHEME_MBE) {
    if (nep->npart>1) {
      for (i=0;i<nep->nt;i++) { ierr = FNDestroy(&ctx->fn[i]);CHKERRQ(ierr); }
      ierr = PetscFree(ctx->fn);CHKERRQ(ierr);
    }
  }
  ierr = MatDestroy(&T);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(NEP_Refine,nep,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

