/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "davidson"

   Step: improve the eigenvectors X with GD2
*/

#include "davidson.h"

typedef struct {
  PetscInt size_X;
} dvdImprovex_gd2;

static PetscErrorCode dvd_improvex_gd2_d(dvdDashboard *d)
{
  PetscErrorCode  ierr;
  dvdImprovex_gd2 *data = (dvdImprovex_gd2*)d->improveX_data;

  PetscFunctionBegin;
  /* Free local data and objects */
  ierr = PetscFree(data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_improvex_gd2_gen(dvdDashboard *d,PetscInt r_s,PetscInt r_e,PetscInt *size_D)
{
  dvdImprovex_gd2 *data = (dvdImprovex_gd2*)d->improveX_data;
  PetscErrorCode  ierr;
  PetscInt        i,j,n,s,ld,lv,kv,max_size_D;
  PetscInt        oldnpreconv = d->npreconv;
  PetscScalar     *pX,*b;
  Vec             *Ax,*Bx,v,*x;
  Mat             M;
  BV              X;

  PetscFunctionBegin;
  /* Compute the number of pairs to improve */
  ierr = BVGetActiveColumns(d->eps->V,&lv,&kv);CHKERRQ(ierr);
  max_size_D = d->eps->ncv-kv;
  n = PetscMin(PetscMin(data->size_X*2,max_size_D),(r_e-r_s)*2)/2;

  /* Quick exit */
  if (max_size_D == 0 || r_e-r_s <= 0 || n == 0) {
   *size_D = 0;
    PetscFunctionReturn(0);
  }

  ierr = BVDuplicateResize(d->eps->V,4,&X);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,4,2,NULL,&M);CHKERRQ(ierr);

  /* Compute the eigenvectors of the selected pairs */
  for (i=r_s;i<r_s+n; i++) {
    ierr = DSVectors(d->eps->ds,DS_MAT_X,&i,NULL);CHKERRQ(ierr);
  }
  ierr = DSGetArray(d->eps->ds,DS_MAT_X,&pX);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(d->eps->ds,&ld);CHKERRQ(ierr);

  ierr = SlepcVecPoolGetVecs(d->auxV,n,&Ax);CHKERRQ(ierr);
  ierr = SlepcVecPoolGetVecs(d->auxV,n,&Bx);CHKERRQ(ierr);

  /* Bx <- B*X(i) */
  if (d->BX) {
    /* Compute the norms of the eigenvectors */
    if (d->correctXnorm) {
      ierr = dvd_improvex_compute_X(d,r_s,r_s+n,Bx,pX,ld);CHKERRQ(ierr);
    } else {
      for (i=0;i<n;i++) d->nX[r_s+i] = 1.0;
    }
    for (i=0;i<n;i++) {
      ierr = BVMultVec(d->BX,1.0,0.0,Bx[i],&pX[ld*(r_s+i)]);CHKERRQ(ierr);
    }
  } else if (d->B) {
    ierr = SlepcVecPoolGetVecs(d->auxV,1,&x);CHKERRQ(ierr);
    for (i=0;i<n;i++) {
      /* auxV(0) <- X(i) */
      ierr = dvd_improvex_compute_X(d,r_s+i,r_s+i+1,x,pX,ld);CHKERRQ(ierr);
      /* Bx(i) <- B*auxV(0) */
      ierr = MatMult(d->B,x[0],Bx[i]);CHKERRQ(ierr);
    }
    ierr = SlepcVecPoolRestoreVecs(d->auxV,1,&x);CHKERRQ(ierr);
  } else {
    /* Bx <- X */
    ierr = dvd_improvex_compute_X(d,r_s,r_s+n,Bx,pX,ld);CHKERRQ(ierr);
  }

  /* Ax <- A*X(i) */
  for (i=0;i<n;i++) {
    ierr = BVMultVec(d->AX,1.0,0.0,Ax[i],&pX[ld*(i+r_s)]);CHKERRQ(ierr);
  }

  ierr = DSRestoreArray(d->eps->ds,DS_MAT_X,&pX);CHKERRQ(ierr);

  for (i=0,s=0;i<n;i+=s) {
#if !defined(PETSC_USE_COMPLEX)
    if (d->eigi[r_s+i] != 0.0 && i+2<=n) {
       /* [Ax_i Ax_i+1 Bx_i Bx_i+1]*= [   1        0
                                          0        1
                                       -eigr_i -eigi_i
                                        eigi_i -eigr_i] */
      ierr = MatDenseGetArray(M,&b);CHKERRQ(ierr);
      b[0] = b[5] = 1.0/d->nX[r_s+i];
      b[2] = b[7] = -d->eigr[r_s+i]/d->nX[r_s+i];
      b[6] = -(b[3] = d->eigi[r_s+i]/d->nX[r_s+i]);
      b[1] = b[4] = 0.0;
      ierr = MatDenseRestoreArray(M,&b);CHKERRQ(ierr);
      ierr = BVInsertVec(X,0,Ax[i]);CHKERRQ(ierr);
      ierr = BVInsertVec(X,1,Ax[i+1]);CHKERRQ(ierr);
      ierr = BVInsertVec(X,2,Bx[i]);CHKERRQ(ierr);
      ierr = BVInsertVec(X,3,Bx[i+1]);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(X,0,4);CHKERRQ(ierr);
      ierr = BVMultInPlace(X,M,0,2);CHKERRQ(ierr);
      ierr = BVCopyVec(X,0,Ax[i]);CHKERRQ(ierr);
      ierr = BVCopyVec(X,1,Ax[i+1]);CHKERRQ(ierr);
      s = 2;
    } else
#endif
    {
      /* [Ax_i Bx_i]*= [ 1/nX_i    conj(eig_i/nX_i)
                       -eig_i/nX_i     1/nX_i       ] */
      ierr = MatDenseGetArray(M,&b);CHKERRQ(ierr);
      b[0] = 1.0/d->nX[r_s+i];
      b[1] = -d->eigr[r_s+i]/d->nX[r_s+i];
      b[4] = PetscConj(d->eigr[r_s+i]/d->nX[r_s+i]);
      b[5] = 1.0/d->nX[r_s+i];
      ierr = MatDenseRestoreArray(M,&b);CHKERRQ(ierr);
      ierr = BVInsertVec(X,0,Ax[i]);CHKERRQ(ierr);
      ierr = BVInsertVec(X,1,Bx[i]);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(X,0,2);CHKERRQ(ierr);
      ierr = BVMultInPlace(X,M,0,2);CHKERRQ(ierr);
      ierr = BVCopyVec(X,0,Ax[i]);CHKERRQ(ierr);
      ierr = BVCopyVec(X,1,Bx[i]);CHKERRQ(ierr);
      s = 1;
    }
    for (j=0;j<s;j++) d->nX[r_s+i+j] = 1.0;

    /* Ax = R <- P*(Ax - eig_i*Bx) */
    ierr = d->calcpairs_proj_res(d,r_s+i,r_s+i+s,&Ax[i]);CHKERRQ(ierr);

    /* Check if the first eigenpairs are converged */
    if (i == 0) {
      ierr = d->preTestConv(d,0,r_s+s,r_s+s,&d->npreconv);CHKERRQ(ierr);
      if (d->npreconv > oldnpreconv) break;
    }
  }

  /* D <- K*[Ax Bx] */
  if (d->npreconv <= oldnpreconv) {
    for (i=0;i<n;i++) {
      ierr = BVGetColumn(d->eps->V,kv+i,&v);CHKERRQ(ierr);
      ierr = d->improvex_precond(d,r_s+i,Ax[i],v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(d->eps->V,kv+i,&v);CHKERRQ(ierr);
    }
    for (i=n;i<n*2;i++) {
      ierr = BVGetColumn(d->eps->V,kv+i,&v);CHKERRQ(ierr);
      ierr = d->improvex_precond(d,r_s+i-n,Bx[i-n],v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(d->eps->V,kv+i,&v);CHKERRQ(ierr);
    }
    *size_D = 2*n;
#if !defined(PETSC_USE_COMPLEX)
    if (d->eigi[r_s] != 0.0) {
      s = 4;
    } else
#endif
    {
      s = 2;
    }
    /* Prevent that short vectors are discarded in the orthogonalization */
    for (i=0; i<s && i<*size_D; i++) {
      if (d->eps->errest[d->nconv+r_s+i] > PETSC_MACHINE_EPSILON && d->eps->errest[d->nconv+r_s+i] < PETSC_MAX_REAL) {
        ierr = BVScaleColumn(d->eps->V,i+kv,1.0/d->eps->errest[d->nconv+r_s+i]);CHKERRQ(ierr);
      }
    }
  } else *size_D = 0;

  ierr = SlepcVecPoolRestoreVecs(d->auxV,n,&Bx);CHKERRQ(ierr);
  ierr = SlepcVecPoolRestoreVecs(d->auxV,n,&Ax);CHKERRQ(ierr);
  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_improvex_gd2(dvdDashboard *d,dvdBlackboard *b,KSP ksp,PetscInt max_bs)
{
  PetscErrorCode  ierr;
  dvdImprovex_gd2 *data;
  PC              pc;

  PetscFunctionBegin;
  /* Setting configuration constrains */
  /* If the arithmetic is real and the problem is not Hermitian, then
     the block size is incremented in one */
#if !defined(PETSC_USE_COMPLEX)
  if (!DVD_IS(d->sEP, DVD_EP_HERMITIAN)) {
    max_bs++;
    b->max_size_P = PetscMax(b->max_size_P,2);
  } else
#endif
  {
    b->max_size_P = PetscMax(b->max_size_P,1);
  }
  b->max_size_X = PetscMax(b->max_size_X,max_bs);

  /* Setup the preconditioner */
  if (ksp) {
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = dvd_static_precond_PC(d,b,pc);CHKERRQ(ierr);
  } else {
    ierr = dvd_static_precond_PC(d,b,0);CHKERRQ(ierr);
  }

  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    ierr = PetscNewLog(d->eps,&data);CHKERRQ(ierr);
    d->improveX_data = data;
    data->size_X = b->max_size_X;
    d->improveX = dvd_improvex_gd2_gen;

    ierr = EPSDavidsonFLAdd(&d->destroyList,dvd_improvex_gd2_d);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

