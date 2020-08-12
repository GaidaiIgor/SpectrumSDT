/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the LAPACK SVD subroutines
*/

#include <slepc/private/svdimpl.h>

PetscErrorCode SVDSetUp_LAPACK(SVD svd)
{
  PetscErrorCode ierr;
  PetscInt       M,N;

  PetscFunctionBegin;
  ierr = SVDMatGetSize(svd,&M,&N);CHKERRQ(ierr);
  svd->ncv = N;
  if (svd->mpd!=PETSC_DEFAULT) { ierr = PetscInfo(svd,"Warning: parameter mpd ignored\n");CHKERRQ(ierr); }
  if (svd->stop!=SVD_STOP_BASIC) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"User-defined stopping test not supported in this solver");
  if (svd->max_it==PETSC_DEFAULT) svd->max_it = 1;
  svd->leftbasis = PETSC_TRUE;
  ierr = SVDAllocateSolution(svd,0);CHKERRQ(ierr);
  ierr = DSSetType(svd->ds,DSSVD);CHKERRQ(ierr);
  ierr = DSAllocate(svd->ds,PetscMax(M,N));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSolve_LAPACK(SVD svd)
{
  PetscErrorCode ierr;
  PetscInt       M,N,n,i,j,k,ld,lowu,lowv,highu,highv;
  Mat            Ar,mat;
  Vec            u,v;
  PetscScalar    *pU,*pVT,*pmat,*pu,*pv,*A,*w;

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(svd->ds,&ld);CHKERRQ(ierr);
  ierr = MatCreateRedundantMatrix(svd->OP,0,PETSC_COMM_SELF,MAT_INITIAL_MATRIX,&Ar);CHKERRQ(ierr);
  ierr = MatConvert(Ar,MATSEQDENSE,MAT_INITIAL_MATRIX,&mat);CHKERRQ(ierr);
  ierr = MatDestroy(&Ar);CHKERRQ(ierr);
  ierr = MatGetSize(mat,&M,&N);CHKERRQ(ierr);
  ierr = DSSetDimensions(svd->ds,M,N,0,0);CHKERRQ(ierr);
  ierr = MatDenseGetArray(mat,&pmat);CHKERRQ(ierr);
  ierr = DSGetArray(svd->ds,DS_MAT_A,&A);CHKERRQ(ierr);
  for (i=0;i<M;i++)
    for (j=0;j<N;j++)
      A[i+j*ld] = pmat[i+j*M];
  ierr = DSRestoreArray(svd->ds,DS_MAT_A,&A);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(mat,&pmat);CHKERRQ(ierr);
  ierr = DSSetState(svd->ds,DS_STATE_RAW);CHKERRQ(ierr);

  n = PetscMin(M,N);
  ierr = PetscMalloc1(n,&w);CHKERRQ(ierr);
  ierr = DSSolve(svd->ds,w,NULL);CHKERRQ(ierr);
  ierr = DSSort(svd->ds,w,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DSSynchronize(svd->ds,w,NULL);CHKERRQ(ierr);

  /* copy singular vectors */
  ierr = DSGetArray(svd->ds,DS_MAT_U,&pU);CHKERRQ(ierr);
  ierr = DSGetArray(svd->ds,DS_MAT_VT,&pVT);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
    if (svd->which == SVD_SMALLEST) k = n - i - 1;
    else k = i;
    svd->sigma[k] = PetscRealPart(w[i]);
    ierr = BVGetColumn(svd->U,k,&u);CHKERRQ(ierr);
    ierr = BVGetColumn(svd->V,k,&v);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(u,&lowu,&highu);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(v,&lowv,&highv);CHKERRQ(ierr);
    ierr = VecGetArray(u,&pu);CHKERRQ(ierr);
    ierr = VecGetArray(v,&pv);CHKERRQ(ierr);
    if (M>=N) {
      for (j=lowu;j<highu;j++) pu[j-lowu] = pU[i*ld+j];
      for (j=lowv;j<highv;j++) pv[j-lowv] = PetscConj(pVT[j*ld+i]);
    } else {
      for (j=lowu;j<highu;j++) pu[j-lowu] = PetscConj(pVT[j*ld+i]);
      for (j=lowv;j<highv;j++) pv[j-lowv] = pU[i*ld+j];
    }
    ierr = VecRestoreArray(u,&pu);CHKERRQ(ierr);
    ierr = VecRestoreArray(v,&pv);CHKERRQ(ierr);
    ierr = BVRestoreColumn(svd->U,k,&u);CHKERRQ(ierr);
    ierr = BVRestoreColumn(svd->V,k,&v);CHKERRQ(ierr);
  }
  ierr = DSRestoreArray(svd->ds,DS_MAT_U,&pU);CHKERRQ(ierr);
  ierr = DSRestoreArray(svd->ds,DS_MAT_VT,&pVT);CHKERRQ(ierr);

  svd->nconv = n;
  svd->reason = SVD_CONVERGED_TOL;

  ierr = MatDestroy(&mat);CHKERRQ(ierr);
  ierr = PetscFree(w);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_LAPACK(SVD svd)
{
  PetscFunctionBegin;
  svd->ops->setup   = SVDSetUp_LAPACK;
  svd->ops->solve   = SVDSolve_LAPACK;
  PetscFunctionReturn(0);
}

