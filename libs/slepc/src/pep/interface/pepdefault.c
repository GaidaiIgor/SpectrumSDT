/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Simple default routines for common PEP operations
*/

#include <slepc/private/pepimpl.h>     /*I "slepcpep.h" I*/

/*@
   PEPSetWorkVecs - Sets a number of work vectors into a PEP object.

   Collective on pep

   Input Parameters:
+  pep - polynomial eigensolver context
-  nw  - number of work vectors to allocate

   Developers Note:
   This is SLEPC_EXTERN because it may be required by user plugin PEP
   implementations.

   Level: developer
@*/
PetscErrorCode PEPSetWorkVecs(PEP pep,PetscInt nw)
{
  PetscErrorCode ierr;
  Vec            t;

  PetscFunctionBegin;
  if (pep->nwork < nw) {
    ierr = VecDestroyVecs(pep->nwork,&pep->work);CHKERRQ(ierr);
    pep->nwork = nw;
    ierr = BVGetColumn(pep->V,0,&t);CHKERRQ(ierr);
    ierr = VecDuplicateVecs(t,nw,&pep->work);CHKERRQ(ierr);
    ierr = BVRestoreColumn(pep->V,0,&t);CHKERRQ(ierr);
    ierr = PetscLogObjectParents(pep,nw,pep->work);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  PEPConvergedRelative - Checks convergence relative to the eigenvalue.
*/
PetscErrorCode PEPConvergedRelative(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal w;

  PetscFunctionBegin;
  w = SlepcAbsEigenvalue(eigr,eigi);
  *errest = res/w;
  PetscFunctionReturn(0);
}

/*
  PEPConvergedNorm - Checks convergence relative to the matrix norms.
*/
PetscErrorCode PEPConvergedNorm(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscReal      w=0.0,t;
  PetscInt       j;
  PetscBool      flg;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* initialization of matrix norms */
  if (!pep->nrma[pep->nmat-1]) {
    for (j=0;j<pep->nmat;j++) {
      ierr = MatHasOperation(pep->A[j],MATOP_NORM,&flg);CHKERRQ(ierr);
      if (!flg) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"The convergence test related to the matrix norms requires a matrix norm operation");
      ierr = MatNorm(pep->A[j],NORM_INFINITY,&pep->nrma[j]);CHKERRQ(ierr);
    }
  }
  t = SlepcAbsEigenvalue(eigr,eigi);
  for (j=pep->nmat-1;j>=0;j--) {
    w = w*t+pep->nrma[j];
  }
  *errest = res/w;
  PetscFunctionReturn(0);
}

/*
  PEPConvergedAbsolute - Checks convergence absolutely.
*/
PetscErrorCode PEPConvergedAbsolute(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = res;
  PetscFunctionReturn(0);
}

/*@C
   PEPStoppingBasic - Default routine to determine whether the outer eigensolver
   iteration must be stopped.

   Collective on pep

   Input Parameters:
+  pep    - eigensolver context obtained from PEPCreate()
.  its    - current number of iterations
.  max_it - maximum number of iterations
.  nconv  - number of currently converged eigenpairs
.  nev    - number of requested eigenpairs
-  ctx    - context (not used here)

   Output Parameter:
.  reason - result of the stopping test

   Notes:
   A positive value of reason indicates that the iteration has finished successfully
   (converged), and a negative value indicates an error condition (diverged). If
   the iteration needs to be continued, reason must be set to PEP_CONVERGED_ITERATING
   (zero).

   PEPStoppingBasic() will stop if all requested eigenvalues are converged, or if
   the maximum number of iterations has been reached.

   Use PEPSetStoppingTest() to provide your own test instead of using this one.

   Level: advanced

.seealso: PEPSetStoppingTest(), PEPConvergedReason, PEPGetConvergedReason()
@*/
PetscErrorCode PEPStoppingBasic(PEP pep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,PEPConvergedReason *reason,void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *reason = PEP_CONVERGED_ITERATING;
  if (nconv >= nev) {
    ierr = PetscInfo2(pep,"Polynomial eigensolver finished successfully: %D eigenpairs converged at iteration %D\n",nconv,its);CHKERRQ(ierr);
    *reason = PEP_CONVERGED_TOL;
  } else if (its >= max_it) {
    *reason = PEP_DIVERGED_ITS;
    ierr = PetscInfo1(pep,"Polynomial eigensolver iteration reached maximum number of iterations (%D)\n",its);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPBackTransform_Default(PEP pep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = STBackTransform(pep->st,pep->nconv,pep->eigr,pep->eigi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PEPComputeVectors_Default(PEP pep)
{
  PetscErrorCode ierr;
  PetscInt       i;
  Vec            v;
#if !defined(PETSC_USE_COMPLEX)
  Vec            v1;
#endif

  PetscFunctionBegin;
  ierr = PEPExtractVectors(pep);CHKERRQ(ierr);

  /* Fix eigenvectors if balancing was used */
  if ((pep->scale==PEP_SCALE_DIAGONAL || pep->scale==PEP_SCALE_BOTH) && pep->Dr && (pep->refine!=PEP_REFINE_MULTIPLE)) {
    for (i=0;i<pep->nconv;i++) {
      ierr = BVGetColumn(pep->V,i,&v);CHKERRQ(ierr);
      ierr = VecPointwiseMult(v,v,pep->Dr);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pep->V,i,&v);CHKERRQ(ierr);
    }
  }

  /* normalization */
  for (i=0;i<pep->nconv;i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (pep->eigi[i]!=0.0) {   /* first eigenvalue of a complex conjugate pair */
      ierr = BVGetColumn(pep->V,i,&v);CHKERRQ(ierr);
      ierr = BVGetColumn(pep->V,i+1,&v1);CHKERRQ(ierr);
      ierr = VecNormalizeComplex(v,v1,PETSC_TRUE,NULL);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pep->V,i,&v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pep->V,i+1,&v1);CHKERRQ(ierr);
      i++;
    } else   /* real eigenvalue */
#endif
    {
      ierr = BVGetColumn(pep->V,i,&v);CHKERRQ(ierr);
      ierr = VecNormalizeComplex(v,NULL,PETSC_FALSE,NULL);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pep->V,i,&v);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*
  PEPBuildDiagonalScaling - compute two diagonal matrices to be applied for balancing
  in polynomial eigenproblems.
*/
PetscErrorCode PEPBuildDiagonalScaling(PEP pep)
{
  PetscErrorCode ierr;
  PetscInt       it,i,j,k,nmat,nr,e,nz,lst,lend,nc=0,*cols,emax,emin,emaxl,eminl;
  const PetscInt *cidx,*ridx;
  Mat            M,*T,A;
  PetscMPIInt    n;
  PetscBool      cont=PETSC_TRUE,flg=PETSC_FALSE;
  PetscScalar    *array,*Dr,*Dl,t;
  PetscReal      l2,d,*rsum,*aux,*csum,w=1.0;
  MatStructure   str;
  MatInfo        info;

  PetscFunctionBegin;
  l2 = 2*PetscLogReal(2.0);
  nmat = pep->nmat;
  ierr = PetscMPIIntCast(pep->n,&n);CHKERRQ(ierr);
  ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
  ierr = PetscMalloc1(nmat,&T);CHKERRQ(ierr);
  for (k=0;k<nmat;k++) {
    ierr = STGetMatrixTransformed(pep->st,k,&T[k]);CHKERRQ(ierr);
  }
  /* Form local auxiliar matrix M */
  ierr = PetscObjectTypeCompareAny((PetscObject)T[0],&cont,MATMPIAIJ,MATSEQAIJ,"");CHKERRQ(ierr);
  if (!cont) SETERRQ(PetscObjectComm((PetscObject)T[0]),PETSC_ERR_SUP,"Only for MPIAIJ or SEQAIJ matrix types");
  ierr = PetscObjectTypeCompare((PetscObject)T[0],MATMPIAIJ,&cont);CHKERRQ(ierr);
  if (cont) {
    ierr = MatMPIAIJGetLocalMat(T[0],MAT_INITIAL_MATRIX,&M);CHKERRQ(ierr);
    flg = PETSC_TRUE;
  } else {
    ierr = MatDuplicate(T[0],MAT_COPY_VALUES,&M);CHKERRQ(ierr);
  }
  ierr = MatGetInfo(M,MAT_LOCAL,&info);CHKERRQ(ierr);
  nz = (PetscInt)info.nz_used;
  ierr = MatSeqAIJGetArray(M,&array);CHKERRQ(ierr);
  for (i=0;i<nz;i++) {
    t = PetscAbsScalar(array[i]);
    array[i] = t*t;
  }
  ierr = MatSeqAIJRestoreArray(M,&array);CHKERRQ(ierr);
  for (k=1;k<nmat;k++) {
    if (flg) {
      ierr = MatMPIAIJGetLocalMat(T[k],MAT_INITIAL_MATRIX,&A);CHKERRQ(ierr);
    } else {
      if (str==SAME_NONZERO_PATTERN) {
        ierr = MatCopy(T[k],A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
      } else {
        ierr = MatDuplicate(T[k],MAT_COPY_VALUES,&A);CHKERRQ(ierr);
      }
    }
    ierr = MatGetInfo(A,MAT_LOCAL,&info);CHKERRQ(ierr);
    nz = (PetscInt)info.nz_used;
    ierr = MatSeqAIJGetArray(A,&array);CHKERRQ(ierr);
    for (i=0;i<nz;i++) {
      t = PetscAbsScalar(array[i]);
      array[i] = t*t;
    }
    ierr = MatSeqAIJRestoreArray(A,&array);CHKERRQ(ierr);
    w *= pep->slambda*pep->slambda*pep->sfactor;
    ierr = MatAXPY(M,w,A,str);CHKERRQ(ierr);
    if (flg || str!=SAME_NONZERO_PATTERN || k==nmat-2) {
      ierr = MatDestroy(&A);CHKERRQ(ierr);
    }
  }
  ierr = MatGetRowIJ(M,0,PETSC_FALSE,PETSC_FALSE,&nr,&ridx,&cidx,&cont);CHKERRQ(ierr);
  if (!cont) SETERRQ(PetscObjectComm((PetscObject)T[0]),PETSC_ERR_SUP,"It is not possible to compute scaling diagonals for these PEP matrices");
  ierr = MatGetInfo(M,MAT_LOCAL,&info);CHKERRQ(ierr);
  nz = (PetscInt)info.nz_used;
  ierr = VecGetOwnershipRange(pep->Dl,&lst,&lend);CHKERRQ(ierr);
  ierr = PetscMalloc4(nr,&rsum,pep->n,&csum,pep->n,&aux,PetscMin(pep->n-lend+lst,nz),&cols);CHKERRQ(ierr);
  ierr = VecSet(pep->Dr,1.0);CHKERRQ(ierr);
  ierr = VecSet(pep->Dl,1.0);CHKERRQ(ierr);
  ierr = VecGetArray(pep->Dl,&Dl);CHKERRQ(ierr);
  ierr = VecGetArray(pep->Dr,&Dr);CHKERRQ(ierr);
  ierr = MatSeqAIJGetArray(M,&array);CHKERRQ(ierr);
  ierr = PetscArrayzero(aux,pep->n);CHKERRQ(ierr);
  for (j=0;j<nz;j++) {
    /* Search non-zero columns outsize lst-lend */
    if (aux[cidx[j]]==0 && (cidx[j]<lst || lend<=cidx[j])) cols[nc++] = cidx[j];
    /* Local column sums */
    aux[cidx[j]] += PetscAbsScalar(array[j]);
  }
  for (it=0;it<pep->sits && cont;it++) {
    emaxl = 0; eminl = 0;
    /* Column sum  */
    if (it>0) { /* it=0 has been already done*/
      ierr = MatSeqAIJGetArray(M,&array);CHKERRQ(ierr);
      ierr = PetscArrayzero(aux,pep->n);CHKERRQ(ierr);
      for (j=0;j<nz;j++) aux[cidx[j]] += PetscAbsScalar(array[j]);
      ierr = MatSeqAIJRestoreArray(M,&array);CHKERRQ(ierr);
    }
    ierr = MPI_Allreduce(aux,csum,n,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)pep->Dr));CHKERRQ(ierr);
    /* Update Dr */
    for (j=lst;j<lend;j++) {
      d = PetscLogReal(csum[j])/l2;
      e = -(PetscInt)((d < 0)?(d-0.5):(d+0.5));
      d = PetscPowReal(2.0,e);
      Dr[j-lst] *= d;
      aux[j] = d*d;
      emaxl = PetscMax(emaxl,e);
      eminl = PetscMin(eminl,e);
    }
    for (j=0;j<nc;j++) {
      d = PetscLogReal(csum[cols[j]])/l2;
      e = -(PetscInt)((d < 0)?(d-0.5):(d+0.5));
      d = PetscPowReal(2.0,e);
      aux[cols[j]] = d*d;
      emaxl = PetscMax(emaxl,e);
      eminl = PetscMin(eminl,e);
    }
    /* Scale M */
    ierr = MatSeqAIJGetArray(M,&array);CHKERRQ(ierr);
    for (j=0;j<nz;j++) {
      array[j] *= aux[cidx[j]];
    }
    ierr = MatSeqAIJRestoreArray(M,&array);CHKERRQ(ierr);
    /* Row sum */
    ierr = PetscArrayzero(rsum,nr);CHKERRQ(ierr);
    ierr = MatSeqAIJGetArray(M,&array);CHKERRQ(ierr);
    for (i=0;i<nr;i++) {
      for (j=ridx[i];j<ridx[i+1];j++) rsum[i] += PetscAbsScalar(array[j]);
      /* Update Dl */
      d = PetscLogReal(rsum[i])/l2;
      e = -(PetscInt)((d < 0)?(d-0.5):(d+0.5));
      d = PetscPowReal(2.0,e);
      Dl[i] *= d;
      /* Scale M */
      for (j=ridx[i];j<ridx[i+1];j++) array[j] *= d*d;
      emaxl = PetscMax(emaxl,e);
      eminl = PetscMin(eminl,e);
    }
    ierr = MatSeqAIJRestoreArray(M,&array);CHKERRQ(ierr);
    /* Compute global max and min */
    ierr = MPI_Allreduce(&emaxl,&emax,1,MPIU_INT,MPI_MAX,PetscObjectComm((PetscObject)pep->Dl));CHKERRQ(ierr);
    ierr = MPI_Allreduce(&eminl,&emin,1,MPIU_INT,MPI_MIN,PetscObjectComm((PetscObject)pep->Dl));CHKERRQ(ierr);
    if (emax<=emin+2) cont = PETSC_FALSE;
  }
  ierr = VecRestoreArray(pep->Dr,&Dr);CHKERRQ(ierr);
  ierr = VecRestoreArray(pep->Dl,&Dl);CHKERRQ(ierr);
  /* Free memory*/
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = PetscFree4(rsum,csum,aux,cols);CHKERRQ(ierr);
  ierr = PetscFree(T);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   PEPComputeScaleFactor - compute sfactor as described in [Betcke 2008].
*/
PetscErrorCode PEPComputeScaleFactor(PEP pep)
{
  PetscErrorCode ierr;
  PetscBool      has0,has1,flg;
  PetscReal      norm0,norm1;
  Mat            T[2];
  PEPBasis       basis;
  PetscInt       i;

  PetscFunctionBegin;
  if (pep->scale==PEP_SCALE_NONE || pep->scale==PEP_SCALE_DIAGONAL) {  /* no scalar scaling */
    pep->sfactor = 1.0;
    pep->dsfactor = 1.0;
    PetscFunctionReturn(0);
  }
  if (pep->sfactor_set) PetscFunctionReturn(0);  /* user provided value */
  pep->sfactor = 1.0;
  pep->dsfactor = 1.0;
  ierr = PEPGetBasis(pep,&basis);CHKERRQ(ierr);
  if (basis==PEP_BASIS_MONOMIAL) {
    ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = STGetMatrixTransformed(pep->st,0,&T[0]);CHKERRQ(ierr);
      ierr = STGetMatrixTransformed(pep->st,pep->nmat-1,&T[1]);CHKERRQ(ierr);
    } else {
      T[0] = pep->A[0];
      T[1] = pep->A[pep->nmat-1];
    }
    if (pep->nmat>2) {
      ierr = MatHasOperation(T[0],MATOP_NORM,&has0);CHKERRQ(ierr);
      ierr = MatHasOperation(T[1],MATOP_NORM,&has1);CHKERRQ(ierr);
      if (has0 && has1) {
        ierr = MatNorm(T[0],NORM_INFINITY,&norm0);CHKERRQ(ierr);
        ierr = MatNorm(T[1],NORM_INFINITY,&norm1);CHKERRQ(ierr);
        pep->sfactor = PetscPowReal(norm0/norm1,1.0/(pep->nmat-1));
        pep->dsfactor = norm1;
        for (i=pep->nmat-2;i>0;i--) {
          ierr = STGetMatrixTransformed(pep->st,i,&T[1]);CHKERRQ(ierr);
          ierr = MatHasOperation(T[1],MATOP_NORM,&has1);CHKERRQ(ierr);
          if (has1) {
            ierr = MatNorm(T[1],NORM_INFINITY,&norm1);CHKERRQ(ierr);
            pep->dsfactor = pep->dsfactor*pep->sfactor+norm1;
          } else break;
        }
        if (has1) {
          pep->dsfactor = pep->dsfactor*pep->sfactor+norm0;
          pep->dsfactor = pep->nmat/pep->dsfactor;
        } else pep->dsfactor = 1.0;
      }
    }
  }
  PetscFunctionReturn(0);
}

/*
   PEPBasisCoefficients - compute polynomial basis coefficients
*/
PetscErrorCode PEPBasisCoefficients(PEP pep,PetscReal *pbc)
{
  PetscReal *ca,*cb,*cg;
  PetscInt  k,nmat=pep->nmat;

  PetscFunctionBegin;
  ca = pbc;
  cb = pbc+nmat;
  cg = pbc+2*nmat;
  switch (pep->basis) {
  case PEP_BASIS_MONOMIAL:
    for (k=0;k<nmat;k++) {
      ca[k] = 1.0; cb[k] = 0.0; cg[k] = 0.0;
    }
    break;
  case PEP_BASIS_CHEBYSHEV1:
    ca[0] = 1.0; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = .5; cb[k] = 0.0; cg[k] = .5;
    }
    break;
  case PEP_BASIS_CHEBYSHEV2:
    ca[0] = .5; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = .5; cb[k] = 0.0; cg[k] = .5;
    }
    break;
  case PEP_BASIS_LEGENDRE:
    ca[0] = 1.0; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = k+1; cb[k] = -2*k; cg[k] = k;
    }
    break;
  case PEP_BASIS_LAGUERRE:
    ca[0] = -1.0; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = -(k+1); cb[k] = 2*k+1; cg[k] = -k;
    }
    break;
  case PEP_BASIS_HERMITE:
    ca[0] = .5; cb[0] = 0.0; cg[0] = 0.0;
    for (k=1;k<nmat;k++) {
      ca[k] = .5; cb[k] = 0.0; cg[k] = -k;
    }
    break;
  }
  PetscFunctionReturn(0);
}

