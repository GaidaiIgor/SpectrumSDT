/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Newton refinement for PEP, simple version
*/

#include <slepc/private/pepimpl.h>
#include <slepcblaslapack.h>

#define NREF_MAXIT 10

typedef struct {
  VecScatter *scatter_id,nst;
  Mat        *A;
  Vec        nv,vg,v,w;
} PEPSimpNRefctx;

typedef struct {
  Mat          M1;
  Vec          M2,M3;
  PetscScalar  M4,m3;
} PEP_REFINES_MATSHELL;

static PetscErrorCode MatMult_FS(Mat M ,Vec x,Vec y)
{
  PetscErrorCode       ierr;
  PEP_REFINES_MATSHELL *ctx;
  PetscScalar          t;

  PetscFunctionBegin;
  ierr = MatShellGetContext(M,&ctx);CHKERRQ(ierr);
  ierr = VecDot(x,ctx->M3,&t);CHKERRQ(ierr);
  t *= ctx->m3/ctx->M4;
  ierr = MatMult(ctx->M1,x,y);CHKERRQ(ierr);
  ierr = VecAXPY(y,-t,ctx->M2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSimpleNRefSetUp(PEP pep,PEPSimpNRefctx **ctx_)
{
  PetscErrorCode ierr;
  PetscInt       i,si,j,n0,m0,nloc,*idx1,*idx2,ne;
  IS             is1,is2;
  PEPSimpNRefctx *ctx;
  Vec            v;
  PetscMPIInt    rank,size;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,ctx_);CHKERRQ(ierr);
  ctx = *ctx_;
  if (pep->npart==1) {
    pep->refinesubc = NULL;
    ctx->scatter_id = NULL;
    ctx->A = pep->A;
  } else {
    ierr = PetscMalloc2(pep->nmat,&ctx->A,pep->npart,&ctx->scatter_id);CHKERRQ(ierr);

    /* Duplicate matrices */
    for (i=0;i<pep->nmat;i++) {
      ierr = MatCreateRedundantMatrix(pep->A[i],0,PetscSubcommChild(pep->refinesubc),MAT_INITIAL_MATRIX,&ctx->A[i]);CHKERRQ(ierr);
    }
    ierr = MatCreateVecs(ctx->A[0],&ctx->v,NULL);CHKERRQ(ierr);

    /* Create scatters for sending vectors to each subcommucator */
    ierr = BVGetColumn(pep->V,0,&v);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(v,&n0,&m0);CHKERRQ(ierr);
    ierr = BVRestoreColumn(pep->V,0,&v);CHKERRQ(ierr);
    ierr = VecGetLocalSize(ctx->v,&nloc);CHKERRQ(ierr);
    ierr = PetscMalloc2(m0-n0,&idx1,m0-n0,&idx2);CHKERRQ(ierr);
    ierr = VecCreateMPI(PetscObjectComm((PetscObject)pep),nloc,PETSC_DECIDE,&ctx->vg);CHKERRQ(ierr);
    for (si=0;si<pep->npart;si++) {
      j = 0;
      for (i=n0;i<m0;i++) {
        idx1[j]   = i;
        idx2[j++] = i+pep->n*si;
      }
      ierr = ISCreateGeneral(PetscObjectComm((PetscObject)pep),(m0-n0),idx1,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
      ierr = ISCreateGeneral(PetscObjectComm((PetscObject)pep),(m0-n0),idx2,PETSC_COPY_VALUES,&is2);CHKERRQ(ierr);
      ierr = BVGetColumn(pep->V,0,&v);CHKERRQ(ierr);
      ierr = VecScatterCreate(v,is1,ctx->vg,is2,&ctx->scatter_id[si]);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pep->V,0,&v);CHKERRQ(ierr);
      ierr = ISDestroy(&is1);CHKERRQ(ierr);
      ierr = ISDestroy(&is2);CHKERRQ(ierr);
    }
    ierr = PetscFree2(idx1,idx2);CHKERRQ(ierr);
  }
  if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT){
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)ctx->A[0]),&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)ctx->A[0]),&size);CHKERRQ(ierr);
    if (size>1) {
      if (pep->npart==1) {
        ierr = BVGetColumn(pep->V,0,&v);CHKERRQ(ierr);
      } else v = ctx->v;
      ierr = VecGetOwnershipRange(v,&n0,&m0);CHKERRQ(ierr);
      ne = (rank == size-1)?pep->n:0;
      ierr = VecCreateMPI(PetscObjectComm((PetscObject)ctx->A[0]),ne,PETSC_DECIDE,&ctx->nv);CHKERRQ(ierr);
      ierr = PetscMalloc1(m0-n0,&idx1);CHKERRQ(ierr);
      for (i=n0;i<m0;i++) idx1[i-n0] = i;
      ierr = ISCreateGeneral(PetscObjectComm((PetscObject)ctx->A[0]),(m0-n0),idx1,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
      ierr = VecScatterCreate(v,is1,ctx->nv,is1,&ctx->nst);CHKERRQ(ierr);
      if (pep->npart==1) {
        ierr = BVRestoreColumn(pep->V,0,&v);CHKERRQ(ierr);
      }
      ierr = PetscFree(idx1);CHKERRQ(ierr);
      ierr = ISDestroy(&is1);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*
  Gather Eigenpair idx from subcommunicator with color sc
*/
static PetscErrorCode PEPSimpleNRefGatherEigenpair(PEP pep,PEPSimpNRefctx *ctx,PetscInt sc,PetscInt idx,PetscInt *fail)
{
  PetscErrorCode    ierr;
  PetscMPIInt       nproc,p;
  MPI_Comm          comm=((PetscObject)pep)->comm;
  Vec               v;
  const PetscScalar *array;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(comm,&nproc);CHKERRQ(ierr);
  p = (nproc/pep->npart)*(sc+1)+PetscMin(nproc%pep->npart,sc+1)-1;
  if (pep->npart>1) {
    /* Communicate convergence successful */
    ierr = MPI_Bcast(fail,1,MPIU_INT,p,comm);CHKERRQ(ierr);
    if (!(*fail)) {
      /* Process 0 of subcommunicator sc broadcasts the eigenvalue */
      ierr = MPI_Bcast(&pep->eigr[idx],1,MPIU_SCALAR,p,comm);CHKERRQ(ierr);
      /* Gather pep->V[idx] from the subcommuniator sc */
      ierr = BVGetColumn(pep->V,idx,&v);CHKERRQ(ierr);
      if (pep->refinesubc->color==sc) {
        ierr = VecGetArrayRead(ctx->v,&array);CHKERRQ(ierr);
        ierr = VecPlaceArray(ctx->vg,array);CHKERRQ(ierr);
      }
      ierr = VecScatterBegin(ctx->scatter_id[sc],ctx->vg,v,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx->scatter_id[sc],ctx->vg,v,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
      if (pep->refinesubc->color==sc) {
        ierr = VecResetArray(ctx->vg);CHKERRQ(ierr);
        ierr = VecRestoreArrayRead(ctx->v,&array);CHKERRQ(ierr);
      }
      ierr = BVRestoreColumn(pep->V,idx,&v);CHKERRQ(ierr);
    }
  } else {
    if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT && !(*fail)) {
      ierr = MPI_Bcast(&pep->eigr[idx],1,MPIU_SCALAR,p,comm);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSimpleNRefScatterEigenvector(PEP pep,PEPSimpNRefctx *ctx,PetscInt sc,PetscInt idx)
{
  PetscErrorCode    ierr;
  Vec               v;
  const PetscScalar *array;

  PetscFunctionBegin;
  if (pep->npart>1) {
    ierr = BVGetColumn(pep->V,idx,&v);CHKERRQ(ierr);
    if (pep->refinesubc->color==sc) {
      ierr = VecGetArrayRead(ctx->v,&array);CHKERRQ(ierr);
      ierr = VecPlaceArray(ctx->vg,array);CHKERRQ(ierr);
    }
    ierr = VecScatterBegin(ctx->scatter_id[sc],v,ctx->vg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatter_id[sc],v,ctx->vg,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    if (pep->refinesubc->color==sc) {
      ierr = VecResetArray(ctx->vg);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(ctx->v,&array);CHKERRQ(ierr);
    }
    ierr = BVRestoreColumn(pep->V,idx,&v);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPEvaluateFunctionDerivatives(PEP pep,PetscScalar alpha,PetscScalar *vals)
{
  PetscInt    i,nmat=pep->nmat;
  PetscScalar a0,a1,a2;
  PetscReal   *a=pep->pbc,*b=a+nmat,*g=b+nmat;

  PetscFunctionBegin;
  a0 = 0.0;
  a1 = 1.0;
  vals[0] = 0.0;
  if (nmat>1) vals[1] = 1/a[0];
  for (i=2;i<nmat;i++) {
    a2 = ((alpha-b[i-2])*a1-g[i-2]*a0)/a[i-2];
    vals[i] = (a2+(alpha-b[i-1])*vals[i-1]-g[i-1]*vals[i-2])/a[i-1];
    a0 = a1; a1 = a2;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSimpleNRefSetUpSystem(PEP pep,Mat *A,PEPSimpNRefctx *ctx,PetscInt idx,Mat *Mt,Mat *T,Mat *P,PetscBool ini,Vec t,Vec v)
{
  PetscErrorCode       ierr;
  PetscInt             i,nmat=pep->nmat,ml,m0,n0,m1,mg;
  PetscInt             *dnz,*onz,ncols,*cols2=NULL,*nnz;
  PetscScalar          zero=0.0,*coeffs,*coeffs2;
  PetscMPIInt          rank,size;
  MPI_Comm             comm;
  const PetscInt       *cols;
  const PetscScalar    *vals,*array;
  MatStructure         str;
  PEP_REFINES_MATSHELL *fctx;
  Vec                  w=ctx->w;
  Mat                  M;

  PetscFunctionBegin;
  ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
  ierr = PetscMalloc2(nmat,&coeffs,nmat,&coeffs2);CHKERRQ(ierr);
  switch (pep->scheme) {
  case PEP_REFINE_SCHEME_SCHUR:
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
  case PEP_REFINE_SCHEME_MBE:
    M=*T;
    break;
  case PEP_REFINE_SCHEME_EXPLICIT:
    M=*Mt;
    break;
  }
  if (ini) {
    ierr = MatDuplicate(A[0],MAT_COPY_VALUES,&M);CHKERRQ(ierr);
  } else {
    ierr = MatCopy(A[0],M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  }
  ierr = PEPEvaluateBasis(pep,pep->eigr[idx],0,coeffs,NULL);CHKERRQ(ierr);
  ierr = MatScale(M,coeffs[0]);CHKERRQ(ierr);
  for (i=1;i<nmat;i++) {
    ierr = MatAXPY(M,coeffs[i],A[i],(ini)?str:SUBSET_NONZERO_PATTERN);CHKERRQ(ierr);
  }
  ierr = PEPEvaluateFunctionDerivatives(pep,pep->eigr[idx],coeffs2);CHKERRQ(ierr);
  for (i=0;i<nmat && PetscAbsScalar(coeffs2[i])==0.0;i++);
  ierr = MatMult(A[i],v,w);CHKERRQ(ierr);
  if (coeffs2[i]!=1.0) {
    ierr = VecScale(w,coeffs2[i]);CHKERRQ(ierr);
  }
  for (i++;i<nmat;i++) {
    ierr = MatMult(A[i],v,t);CHKERRQ(ierr);
    ierr = VecAXPY(w,coeffs2[i],t);CHKERRQ(ierr);
  }
  switch (pep->scheme) {
  case PEP_REFINE_SCHEME_EXPLICIT:
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
        ierr = PetscMalloc1(pep->n,&cols2);CHKERRQ(ierr);
        for (i=0;i<pep->n;i++) cols2[i]=i;
      }
      ierr = VecScatterBegin(ctx->nst,v,ctx->nv,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx->nst,v,ctx->nv,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = VecGetArrayRead(ctx->nv,&array);CHKERRQ(ierr);
      if (rank==size-1) {
        ierr = MatSetValues(*T,1,&mg,pep->n,cols2,array,INSERT_VALUES);CHKERRQ(ierr);
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
  case PEP_REFINE_SCHEME_SCHUR:
    fctx->M2 = ctx->w;
    fctx->M3 = v;
    fctx->m3 = 0.0;
    for (i=1;i<nmat-1;i++) fctx->m3 += PetscConj(coeffs[i])*coeffs[i];
    fctx->M4 = 0.0;
    for (i=1;i<nmat-1;i++) fctx->M4 += PetscConj(coeffs[i])*coeffs2[i];
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
  case PEP_REFINE_SCHEME_MBE:
    *T = M;
    *P = M;
    break;
  }
  ierr = PetscFree2(coeffs,coeffs2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPNewtonRefinementSimple(PEP pep,PetscInt *maxits,PetscReal tol,PetscInt k)
{
  PetscErrorCode       ierr;
  PetscInt             i,n,its,idx=0,*idx_sc,*its_sc,color,*fail_sc;
  PetscMPIInt          rank,size;
  Mat                  Mt=NULL,T=NULL,P=NULL;
  MPI_Comm             comm;
  Vec                  r,v,dv,rr=NULL,dvv=NULL,t[2];
  PetscScalar          *array2,deig=0.0,tt[2],ttt;
  const PetscScalar    *array;
  PetscReal            norm,error;
  PetscBool            ini=PETSC_TRUE,sc_pend,solved=PETSC_FALSE;
  PEPSimpNRefctx       *ctx;
  PEP_REFINES_MATSHELL *fctx=NULL;
  KSPConvergedReason   reason;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(PEP_Refine,pep,0,0,0);CHKERRQ(ierr);
  ierr = PEPSimpleNRefSetUp(pep,&ctx);CHKERRQ(ierr);
  its = (maxits)?*maxits:NREF_MAXIT;
  if (!pep->refineksp) { ierr = PEPRefineGetKSP(pep,&pep->refineksp);CHKERRQ(ierr); }
  comm = (pep->npart==1)?PetscObjectComm((PetscObject)pep):PetscSubcommChild(pep->refinesubc);
  if (pep->npart==1) {
    ierr = BVGetColumn(pep->V,0,&v);CHKERRQ(ierr);
  } else v = ctx->v;
  ierr = VecDuplicate(v,&ctx->w);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&dv);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&t[0]);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&t[1]);CHKERRQ(ierr);
  if (pep->npart==1) { ierr = BVRestoreColumn(pep->V,0,&v);CHKERRQ(ierr); }
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = VecGetLocalSize(r,&n);CHKERRQ(ierr);
  ierr = PetscMalloc3(pep->npart,&idx_sc,pep->npart,&its_sc,pep->npart,&fail_sc);CHKERRQ(ierr);
  for (i=0;i<pep->npart;i++) fail_sc[i] = 0;
  for (i=0;i<pep->npart;i++) its_sc[i] = 0;
  color = (pep->npart==1)?0:pep->refinesubc->color;

  /* Loop performing iterative refinements */
  while (!solved) {
    for (i=0;i<pep->npart;i++) {
      sc_pend = PETSC_TRUE;
      if (its_sc[i]==0) {
        idx_sc[i] = idx++;
        if (idx_sc[i]>=k) {
          sc_pend = PETSC_FALSE;
        } else {
          ierr = PEPSimpleNRefScatterEigenvector(pep,ctx,i,idx_sc[i]);CHKERRQ(ierr);
        }
      }  else { /* Gather Eigenpair from subcommunicator i */
        ierr = PEPSimpleNRefGatherEigenpair(pep,ctx,i,idx_sc[i],&fail_sc[i]);CHKERRQ(ierr);
      }
      while (sc_pend) {
        if (!fail_sc[i]) {
          ierr = PEPComputeError(pep,idx_sc[i],PEP_ERROR_BACKWARD,&error);CHKERRQ(ierr);
        }
        if (error<=tol || its_sc[i]>=its || fail_sc[i]) {
          idx_sc[i] = idx++;
          its_sc[i] = 0;
          fail_sc[i] = 0;
          if (idx_sc[i]<k) { ierr = PEPSimpleNRefScatterEigenvector(pep,ctx,i,idx_sc[i]);CHKERRQ(ierr); }
        } else {
          sc_pend = PETSC_FALSE;
          its_sc[i]++;
        }
        if (idx_sc[i]>=k) sc_pend = PETSC_FALSE;
      }
    }
    solved = PETSC_TRUE;
    for (i=0;i<pep->npart&&solved;i++) solved = PetscNot(idx_sc[i]<k);
    if (idx_sc[color]<k) {
#if !defined(PETSC_USE_COMPLEX)
      if (pep->eigi[idx_sc[color]]!=0.0) SETERRQ(PetscObjectComm((PetscObject)pep),1,"Simple Refinement not implemented in real scalars for complex eigenvalues");
#endif
      if (pep->npart==1) {
        ierr = BVGetColumn(pep->V,idx_sc[color],&v);CHKERRQ(ierr);
      } else v = ctx->v;
      ierr = PEPSimpleNRefSetUpSystem(pep,ctx->A,ctx,idx_sc[color],&Mt,&T,&P,ini,t[0],v);CHKERRQ(ierr);
      ierr = KSPSetOperators(pep->refineksp,T,P);CHKERRQ(ierr);
      if (ini) {
        ierr = KSPSetFromOptions(pep->refineksp);CHKERRQ(ierr);
        if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT) {
          ierr = MatCreateVecs(T,&dvv,NULL);CHKERRQ(ierr);
          ierr = VecDuplicate(dvv,&rr);CHKERRQ(ierr);
        }
        ini = PETSC_FALSE;
      }

      switch (pep->scheme) {
      case PEP_REFINE_SCHEME_EXPLICIT:
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
        ierr = KSPSolve(pep->refineksp,rr,dvv);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(pep->refineksp,&reason);CHKERRQ(ierr);
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
          if (rank==size-1) pep->eigr[idx_sc[color]] -= array[n];
          ierr = VecRestoreArrayRead(dvv,&array);CHKERRQ(ierr);
        } else fail_sc[color] = 1;
        break;
      case PEP_REFINE_SCHEME_MBE:
        ierr = MatMult(T,v,r);CHKERRQ(ierr);
        /* Mixed block elimination */
        ierr = VecConjugate(v);CHKERRQ(ierr);
        ierr = KSPSolveTranspose(pep->refineksp,v,t[0]);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(pep->refineksp,&reason);CHKERRQ(ierr);
        if (reason>0) {
          ierr = VecConjugate(t[0]);CHKERRQ(ierr);
          ierr = VecDot(ctx->w,t[0],&tt[0]);CHKERRQ(ierr);
          ierr = KSPSolve(pep->refineksp,ctx->w,t[1]);CHKERRQ(ierr);
          ierr = KSPGetConvergedReason(pep->refineksp,&reason);CHKERRQ(ierr);
          if (reason>0) {
            ierr = VecDot(t[1],v,&tt[1]);CHKERRQ(ierr);
            ierr = VecDot(r,t[0],&ttt);CHKERRQ(ierr);
            tt[0] = ttt/tt[0];
            ierr = VecAXPY(r,-tt[0],ctx->w);CHKERRQ(ierr);
            ierr = KSPSolve(pep->refineksp,r,dv);CHKERRQ(ierr);
            ierr = KSPGetConvergedReason(pep->refineksp,&reason);CHKERRQ(ierr);
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
          pep->eigr[idx_sc[color]] -= deig;
          fail_sc[color] = 0;
        } else {
          ierr = VecConjugate(v);CHKERRQ(ierr);
          fail_sc[color] = 1;
        }
        break;
      case PEP_REFINE_SCHEME_SCHUR:
        fail_sc[color] = 1;
        ierr = MatShellGetContext(T,&fctx);CHKERRQ(ierr);
        if (fctx->M4!=0.0) {
          ierr = MatMult(fctx->M1,v,r);CHKERRQ(ierr);
          ierr = KSPSolve(pep->refineksp,r,dv);CHKERRQ(ierr);
          ierr = KSPGetConvergedReason(pep->refineksp,&reason);CHKERRQ(ierr);
          if (reason>0) {
            ierr = VecDot(dv,v,&deig);CHKERRQ(ierr);
            deig *= -fctx->m3/fctx->M4;
            ierr = VecAXPY(v,-1.0,dv);CHKERRQ(ierr);
            ierr = VecNorm(v,NORM_2,&norm);CHKERRQ(ierr);
            ierr = VecScale(v,1.0/norm);CHKERRQ(ierr);
            pep->eigr[idx_sc[color]] -= deig;
            fail_sc[color] = 0;
          }
        }
        break;
      }
      if (pep->npart==1) { ierr = BVRestoreColumn(pep->V,idx_sc[color],&v);CHKERRQ(ierr); }
    }
  }
  ierr = VecDestroy(&t[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&t[1]);CHKERRQ(ierr);
  ierr = VecDestroy(&dv);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->w);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = PetscFree3(idx_sc,its_sc,fail_sc);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&ctx->nst);CHKERRQ(ierr);
  if (pep->npart>1) {
    ierr = VecDestroy(&ctx->vg);CHKERRQ(ierr);
    ierr = VecDestroy(&ctx->v);CHKERRQ(ierr);
    for (i=0;i<pep->nmat;i++) {
      ierr = MatDestroy(&ctx->A[i]);CHKERRQ(ierr);
    }
    for (i=0;i<pep->npart;i++) {
      ierr = VecScatterDestroy(&ctx->scatter_id[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree2(ctx->A,ctx->scatter_id);CHKERRQ(ierr);
  }
  if (fctx && pep->scheme==PEP_REFINE_SCHEME_SCHUR) {
    ierr = MatDestroy(&P);CHKERRQ(ierr);
    ierr = MatDestroy(&fctx->M1);CHKERRQ(ierr);
    ierr = PetscFree(fctx);CHKERRQ(ierr);
  }
  if (pep->scheme==PEP_REFINE_SCHEME_EXPLICIT) {
    ierr = MatDestroy(&Mt);CHKERRQ(ierr);
    ierr = VecDestroy(&dvv);CHKERRQ(ierr);
    ierr = VecDestroy(&rr);CHKERRQ(ierr);
    ierr = VecDestroy(&ctx->nv);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&T);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(PEP_Refine,pep,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

