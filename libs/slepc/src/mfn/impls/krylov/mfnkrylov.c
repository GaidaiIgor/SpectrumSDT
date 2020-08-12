/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc matrix function solver: "krylov"

   Method: Arnoldi with Eiermann-Ernst restart

   Algorithm:

       Build Arnoldi approximations using f(H) for the Hessenberg matrix H,
       restart by discarding the Krylov basis but keeping H.

   References:

       [1] M. Eiermann and O. Ernst, "A restarted Krylov subspace method
           for the evaluation of matrix functions", SIAM J. Numer. Anal.
           44(6):2481-2504, 2006.
*/

#include <slepc/private/mfnimpl.h>
#include <slepcblaslapack.h>

PetscErrorCode MFNSetUp_Krylov(MFN mfn)
{
  PetscErrorCode ierr;
  PetscInt       N;

  PetscFunctionBegin;
  ierr = MatGetSize(mfn->A,&N,NULL);CHKERRQ(ierr);
  if (mfn->ncv==PETSC_DEFAULT) mfn->ncv = PetscMin(30,N);
  if (mfn->max_it==PETSC_DEFAULT) mfn->max_it = 100;
  ierr = MFNAllocateSolution(mfn,1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MFNSolve_Krylov(MFN mfn,Vec b,Vec x)
{
  PetscErrorCode ierr;
  PetscInt       n=0,m,ld,ldh,j;
  PetscBLASInt   m_,inc=1;
  Mat            G=NULL,H=NULL;
  Vec            F=NULL;
  PetscScalar    *array,*farray,*garray,*harray;
  PetscReal      beta,betaold=0.0,nrm=1.0;
  PetscBool      breakdown,set,flg,symm=PETSC_FALSE;

  PetscFunctionBegin;
  m  = mfn->ncv;
  ld = m+1;
  ierr = PetscCalloc1(ld*ld,&array);CHKERRQ(ierr);

  /* set initial vector to b/||b|| */
  ierr = BVInsertVec(mfn->V,0,b);CHKERRQ(ierr);
  ierr = BVScaleColumn(mfn->V,0,1.0/mfn->bnorm);CHKERRQ(ierr);
  ierr = VecSet(x,0.0);CHKERRQ(ierr);

  /* Restart loop */
  while (mfn->reason == MFN_CONVERGED_ITERATING) {
    mfn->its++;

    /* compute Arnoldi factorization */
    ierr = BVMatArnoldi(mfn->V,mfn->transpose_solve?mfn->AT:mfn->A,array,ld,0,&m,&beta,&breakdown);CHKERRQ(ierr);

    /* save previous Hessenberg matrix in G; allocate new storage for H and f(H) */
    if (mfn->its>1) { G = H; H = NULL; }
    ldh = n+m;
    ierr = MFN_CreateVec(ldh,&F);CHKERRQ(ierr);
    ierr = MFN_CreateDenseMat(ldh,&H);CHKERRQ(ierr);

    /* glue together the previous H and the new H obtained with Arnoldi */
    ierr = MatDenseGetArray(H,&harray);CHKERRQ(ierr);
    for (j=0;j<m;j++) {
      ierr = PetscArraycpy(harray+n+(j+n)*ldh,array+j*ld,m);CHKERRQ(ierr);
    }
    if (mfn->its>1) {
      ierr = MatDenseGetArray(G,&garray);CHKERRQ(ierr);
      for (j=0;j<n;j++) {
        ierr = PetscArraycpy(harray+j*ldh,garray+j*n,n);CHKERRQ(ierr);
      }
      ierr = MatDenseRestoreArray(G,&garray);CHKERRQ(ierr);
      ierr = MatDestroy(&G);CHKERRQ(ierr);
      harray[n+(n-1)*ldh] = betaold;
    }
    ierr = MatDenseRestoreArray(H,&harray);CHKERRQ(ierr);

    if (mfn->its==1) {
      /* set symmetry flag of H from A */
      ierr = MatIsHermitianKnown(mfn->A,&set,&flg);CHKERRQ(ierr);
      symm = set? flg: PETSC_FALSE;
      if (symm) {
        ierr = MatSetOption(H,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
      }
    }

    /* evaluate f(H) */
    ierr = FNEvaluateFunctionMatVec(mfn->fn,H,F);CHKERRQ(ierr);

    /* x += ||b||*V*f(H)*e_1 */
    ierr = VecGetArray(F,&farray);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(m,&m_);CHKERRQ(ierr);
    nrm = BLASnrm2_(&m_,farray+n,&inc);   /* relative norm of the update ||u||/||b|| */
    ierr = MFNMonitor(mfn,mfn->its,nrm);CHKERRQ(ierr);
    for (j=0;j<m;j++) farray[j+n] *= mfn->bnorm;
    ierr = BVSetActiveColumns(mfn->V,0,m);CHKERRQ(ierr);
    ierr = BVMultVec(mfn->V,1.0,1.0,x,farray+n);CHKERRQ(ierr);
    ierr = VecRestoreArray(F,&farray);CHKERRQ(ierr);

    /* check convergence */
    if (mfn->its >= mfn->max_it) mfn->reason = MFN_DIVERGED_ITS;
    if (mfn->its>1) {
      if (m<mfn->ncv || breakdown || beta==0.0 || nrm<mfn->tol) mfn->reason = MFN_CONVERGED_TOL;
    }

    /* restart with vector v_{m+1} */
    if (mfn->reason == MFN_CONVERGED_ITERATING) {
      ierr = BVCopyColumn(mfn->V,m,0);CHKERRQ(ierr);
      n += m;
      betaold = beta;
    }
  }

  ierr = MatDestroy(&H);CHKERRQ(ierr);
  ierr = MatDestroy(&G);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = PetscFree(array);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode MFNCreate_Krylov(MFN mfn)
{
  PetscFunctionBegin;
  mfn->ops->solve          = MFNSolve_Krylov;
  mfn->ops->setup          = MFNSetUp_Krylov;
  PetscFunctionReturn(0);
}
