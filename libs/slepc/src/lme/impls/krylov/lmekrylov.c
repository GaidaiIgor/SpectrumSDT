/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc matrix equation solver: "krylov"

   Method: Arnoldi with Eiermann-Ernst restart

   Algorithm:

       Project the equation onto the Arnoldi basis and solve the compressed
       equation the Hessenberg matrix H, restart by discarding the Krylov
       basis but keeping H.

   References:

       [1] Y. Saad, "Numerical solution of large Lyapunov equations", in
           Signal processing, scattering and operator theory, and numerical
           methods, vol. 5 of Progr. Systems Control Theory, pages 503-511,
           1990.

       [2] D. Kressner, "Memory-efficient Krylov subspace techniques for
           solving large-scale Lyapunov equations", in 2008 IEEE Int. Conf.
           Computer-Aided Control Systems, pages 613-618, 2008.
*/

#include <slepc/private/lmeimpl.h>
#include <slepcblaslapack.h>

PetscErrorCode LMESetUp_Krylov(LME lme)
{
  PetscErrorCode ierr;
  PetscInt       N;

  PetscFunctionBegin;
  ierr = MatGetSize(lme->A,&N,NULL);CHKERRQ(ierr);
  if (lme->ncv==PETSC_DEFAULT) lme->ncv = PetscMin(30,N);
  if (lme->max_it==PETSC_DEFAULT) lme->max_it = 100;
  ierr = LMEAllocateSolution(lme,1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode LMESolve_Krylov_Lyapunov_Vec(LME lme,Vec b,PetscBool fixed,PetscInt rrank,BV C1,BV *X1,PetscInt *col,PetscBool *fail,PetscInt *totalits)
{
  PetscErrorCode ierr;
  PetscInt       n=0,m,ldh,ldg=0,i,j,rank=0,lrank,pass,nouter=0,its;
  PetscReal      bnorm,beta,errest;
  PetscBool      breakdown;
  PetscScalar    *H,*G=NULL,*Gnew=NULL,*L,*U,*r,*Qarray,sone=1.0,zero=0.0;
  PetscBLASInt   n_,m_,rk_;
  Mat            Q;

  PetscFunctionBegin;
  *fail = PETSC_FALSE;
  its = 0;
  m  = lme->ncv;
  ldh = m+1;
  ierr = PetscCalloc1(ldh*m,&H);CHKERRQ(ierr);

  ierr = VecNorm(b,NORM_2,&bnorm);CHKERRQ(ierr);
  if (!bnorm) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONG,"Cannot process a zero vector in the right-hand side");

  for (pass=0;pass<2;pass++) {

    /* set initial vector to b/||b|| */
    ierr = BVInsertVec(lme->V,0,b);CHKERRQ(ierr);
    ierr = BVScaleColumn(lme->V,0,1.0/bnorm);CHKERRQ(ierr);

    /* Restart loop */
    while ((pass==0 && !*fail) || (pass==1 && its+1<nouter)) {
      its++;

      /* compute Arnoldi factorization */
      ierr = BVMatArnoldi(lme->V,lme->A,H,ldh,0,&m,&beta,&breakdown);CHKERRQ(ierr);

      if (pass==0) {
        /* glue together the previous H and the new H obtained with Arnoldi */
        ldg = n+m+1;
        ierr = PetscCalloc1(ldg*(n+m),&Gnew);CHKERRQ(ierr);
        for (j=0;j<m;j++) {
          ierr = PetscArraycpy(Gnew+n+(j+n)*ldg,H+j*ldh,m);CHKERRQ(ierr);
        }
        Gnew[n+m+(n+m-1)*ldg] = beta;
        if (G) {
          for (j=0;j<n;j++) {
            ierr = PetscArraycpy(Gnew+j*ldg,G+j*(n+1),n+1);CHKERRQ(ierr);
          }
          ierr = PetscFree(G);CHKERRQ(ierr);
        }
        G = Gnew;
        n += m;
      } else {
        /* update Z = Z + V(:,1:m)*Q    with   Q=U(blk,:)*P(1:nrk,:)'  */
        ierr = MatCreateSeqDense(PETSC_COMM_SELF,m,*col+rank,NULL,&Q);CHKERRQ(ierr);
        ierr = MatDenseGetArray(Q,&Qarray);CHKERRQ(ierr);
        ierr = PetscBLASIntCast(m,&m_);CHKERRQ(ierr);
        ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
        ierr = PetscBLASIntCast(rank,&rk_);CHKERRQ(ierr);
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&m_,&rk_,&rk_,&sone,U+its*m,&n_,L,&n_,&zero,Qarray+(*col)*m,&m_));
        ierr = MatDenseRestoreArray(Q,&Qarray);CHKERRQ(ierr);
        ierr = BVSetActiveColumns(*X1,*col,*col+rank);CHKERRQ(ierr);
        ierr = BVMult(*X1,1.0,1.0,lme->V,Q);CHKERRQ(ierr);
        ierr = MatDestroy(&Q);CHKERRQ(ierr);
      }

      if (pass==0) {
        /* solve compressed Lyapunov equation */
        ierr = PetscCalloc1(n,&r);CHKERRQ(ierr);
        ierr = PetscCalloc1(n*n,&L);CHKERRQ(ierr);
        r[0] = bnorm;
        errest = PetscAbsScalar(G[n+(n-1)*ldg]);
        ierr = LMEDenseHessLyapunovChol(lme,n,G,ldg,1,r,n,L,n,&errest);CHKERRQ(ierr);
        ierr = LMEMonitor(lme,*totalits+its,errest);CHKERRQ(ierr);
        ierr = PetscFree(r);CHKERRQ(ierr);

        /* check convergence */
        if (errest<lme->tol) {
          lme->errest += errest;
          ierr = PetscMalloc1(n*n,&U);CHKERRQ(ierr);
          /* transpose L */
          for (j=0;j<n;j++) {
            for (i=j+1;i<n;i++) {
              L[i+j*n] = PetscConj(L[j+i*n]);
              L[j+i*n] = 0.0;
            }
          }
          ierr = LMEDenseRankSVD(lme,n,L,n,U,n,&lrank);CHKERRQ(ierr);
          ierr = PetscInfo1(lme,"Rank of the Cholesky factor = %D\n",lrank);CHKERRQ(ierr);
          nouter = its;
          its = -1;
          if (!fixed) {  /* X1 was not set by user, allocate it with rank columns */
            rank = lrank;
            if (*col) {
              ierr = BVResize(*X1,*col+rank,PETSC_TRUE);CHKERRQ(ierr);
            } else {
              ierr = BVDuplicateResize(C1,rank,X1);CHKERRQ(ierr);
            }
          } else rank = PetscMin(lrank,rrank);
          ierr = PetscFree(G);CHKERRQ(ierr);
          break;
        } else {
          ierr = PetscFree(L);CHKERRQ(ierr);
          if (*totalits+its>=lme->max_it) *fail = PETSC_TRUE;
        }
      }

      /* restart with vector v_{m+1} */
      if (!*fail) {
        ierr = BVCopyColumn(lme->V,m,0);CHKERRQ(ierr);
      }
    }
  }

  *col += rank;
  *totalits += its+1;
  ierr = PetscFree(H);CHKERRQ(ierr);
  if (L) { ierr = PetscFree(L);CHKERRQ(ierr); }
  if (U) { ierr = PetscFree(U);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode LMESolve_Krylov_Lyapunov(LME lme)
{
  PetscErrorCode ierr;
  PetscBool      fail,fixed = lme->X? PETSC_TRUE: PETSC_FALSE;
  PetscInt       i,k,rank=0,col=0;
  Vec            b;
  BV             X1=NULL,C1;
  Mat            X1m,X1t,C1m;

  PetscFunctionBegin;
  ierr = MatLRCGetMats(lme->C,NULL,&C1m,NULL,NULL);CHKERRQ(ierr);
  ierr = BVCreateFromMat(C1m,&C1);CHKERRQ(ierr);
  ierr = BVSetFromOptions(C1);CHKERRQ(ierr);
  ierr = BVGetActiveColumns(C1,NULL,&k);CHKERRQ(ierr);
  if (fixed) {
    ierr = MatLRCGetMats(lme->X,NULL,&X1m,NULL,NULL);CHKERRQ(ierr);
    ierr = BVCreateFromMat(X1m,&X1);CHKERRQ(ierr);
    ierr = BVSetFromOptions(X1);CHKERRQ(ierr);
    ierr = BVGetActiveColumns(X1,NULL,&rank);CHKERRQ(ierr);
    rank = rank/k;
  }
  for (i=0;i<k;i++) {
    ierr = BVGetColumn(C1,i,&b);CHKERRQ(ierr);
    ierr = LMESolve_Krylov_Lyapunov_Vec(lme,b,fixed,rank,C1,&X1,&col,&fail,&lme->its);CHKERRQ(ierr);
    ierr = BVRestoreColumn(C1,i,&b);CHKERRQ(ierr);
    if (fail) {
      lme->reason = LME_DIVERGED_ITS;
      break;
    }
  }
  if (lme->reason==LME_CONVERGED_ITERATING) lme->reason = LME_CONVERGED_TOL;
  ierr = BVCreateMat(X1,&X1t);CHKERRQ(ierr);
  if (fixed) {
    ierr = MatCopy(X1t,X1m,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  } else {
    ierr = MatCreateLRC(NULL,X1t,NULL,NULL,&lme->X);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&X1t);CHKERRQ(ierr);
  ierr = BVDestroy(&C1);CHKERRQ(ierr);
  ierr = BVDestroy(&X1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode LMECreate_Krylov(LME lme)
{
  PetscFunctionBegin;
  lme->ops->solve[LME_LYAPUNOV]      = LMESolve_Krylov_Lyapunov;
  lme->ops->setup                    = LMESetUp_Krylov;
  PetscFunctionReturn(0);
}
