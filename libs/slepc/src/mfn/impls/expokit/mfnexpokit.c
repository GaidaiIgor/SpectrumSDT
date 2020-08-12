/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc matrix function solver: "expokit"

   Method: Arnoldi method tailored for the matrix exponential

   Algorithm:

       Uses Arnoldi relations to compute exp(t_step*A)*v_last for
       several time steps.

   References:

       [1] R. Sidje, "Expokit: a software package for computing matrix
           exponentials", ACM Trans. Math. Softw. 24(1):130-156, 1998.
*/

#include <slepc/private/mfnimpl.h>

PetscErrorCode MFNSetUp_Expokit(MFN mfn)
{
  PetscErrorCode ierr;
  PetscInt       N;
  PetscBool      isexp;

  PetscFunctionBegin;
  ierr = MatGetSize(mfn->A,&N,NULL);CHKERRQ(ierr);
  if (mfn->ncv==PETSC_DEFAULT) mfn->ncv = PetscMin(30,N);
  if (mfn->max_it==PETSC_DEFAULT) mfn->max_it = 100;
  ierr = MFNAllocateSolution(mfn,2);CHKERRQ(ierr);

  ierr = PetscObjectTypeCompare((PetscObject)mfn->fn,FNEXP,&isexp);CHKERRQ(ierr);
  if (!isexp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"This solver only supports the exponential function");
  PetscFunctionReturn(0);
}

PetscErrorCode MFNSolve_Expokit(MFN mfn,Vec b,Vec x)
{
  PetscErrorCode ierr;
  PetscInt       mxstep,mxrej,m,mb,ld,i,j,ireject,mx,k1;
  Vec            v,r;
  Mat            M=NULL,K=NULL;
  FN             fn;
  PetscScalar    *H,*B,*F,*betaF,t,sgn,sfactor;
  PetscReal      anorm,avnorm,tol,err_loc,rndoff;
  PetscReal      t_out,t_new,t_now,t_step;
  PetscReal      xm,fact,s,p1,p2;
  PetscReal      beta,beta2,gamma,delta;
  PetscBool      breakdown;

  PetscFunctionBegin;
  m   = mfn->ncv;
  tol = mfn->tol;
  mxstep = mfn->max_it;
  mxrej = 10;
  gamma = 0.9;
  delta = 1.2;
  mb    = m;
  ierr = FNGetScale(mfn->fn,&t,&sfactor);CHKERRQ(ierr);
  ierr = FNDuplicate(mfn->fn,PetscObjectComm((PetscObject)mfn->fn),&fn);CHKERRQ(ierr);
  ierr = FNSetScale(fn,1.0,1.0);CHKERRQ(ierr);
  t_out = PetscAbsScalar(t);
  t_now = 0.0;
  ierr = MatNorm(mfn->A,NORM_INFINITY,&anorm);CHKERRQ(ierr);
  rndoff = anorm*PETSC_MACHINE_EPSILON;

  k1 = 2;
  xm = 1.0/(PetscReal)m;
  beta = mfn->bnorm;
  fact = PetscPowRealInt((m+1)/2.72,m+1)*PetscSqrtReal(2*PETSC_PI*(m+1));
  t_new = (1.0/anorm)*PetscPowReal((fact*tol)/(4.0*beta*anorm),xm);
  s = PetscPowReal(10.0,PetscFloorReal(PetscLog10Real(t_new))-1);
  t_new = PetscCeilReal(t_new/s)*s;
  sgn = t/PetscAbsScalar(t);

  ierr = VecCopy(b,x);CHKERRQ(ierr);
  ld = m+2;
  ierr = PetscCalloc3(m+1,&betaF,ld*ld,&H,ld*ld,&B);CHKERRQ(ierr);

  while (mfn->reason == MFN_CONVERGED_ITERATING) {
    mfn->its++;
    if (PetscIsInfOrNanReal(t_new)) t_new = PETSC_MAX_REAL;
    t_step = PetscMin(t_out-t_now,t_new);
    ierr = BVInsertVec(mfn->V,0,x);CHKERRQ(ierr);
    ierr = BVScaleColumn(mfn->V,0,1.0/beta);CHKERRQ(ierr);
    ierr = BVMatArnoldi(mfn->V,mfn->transpose_solve?mfn->AT:mfn->A,H,ld,0,&mb,&beta2,&breakdown);CHKERRQ(ierr);
    if (breakdown) {
      k1 = 0;
      t_step = t_out-t_now;
    }
    if (k1!=0) {
      H[m+1+ld*m] = 1.0;
      ierr = BVGetColumn(mfn->V,m,&v);CHKERRQ(ierr);
      ierr = BVGetColumn(mfn->V,m+1,&r);CHKERRQ(ierr);
      ierr = MatMult(mfn->transpose_solve?mfn->AT:mfn->A,v,r);CHKERRQ(ierr);
      ierr = BVRestoreColumn(mfn->V,m,&v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(mfn->V,m+1,&r);CHKERRQ(ierr);
      ierr = BVNormColumn(mfn->V,m+1,NORM_2,&avnorm);CHKERRQ(ierr);
    }
    ierr = PetscArraycpy(B,H,ld*ld);CHKERRQ(ierr);

    ireject = 0;
    while (ireject <= mxrej) {
      mx = mb + k1;
      for (i=0;i<mx;i++) {
        for (j=0;j<mx;j++) {
          H[i+j*ld] = sgn*B[i+j*ld]*t_step;
        }
      }
      ierr = MFN_CreateDenseMat(mx,&M);CHKERRQ(ierr);
      ierr = MFN_CreateDenseMat(mx,&K);CHKERRQ(ierr);
      ierr = MatDenseGetArray(M,&F);CHKERRQ(ierr);
      for (j=0;j<mx;j++) {
        ierr = PetscArraycpy(F+j*mx,H+j*ld,mx);CHKERRQ(ierr);
      }
      ierr = MatDenseRestoreArray(M,&F);CHKERRQ(ierr);
      ierr = FNEvaluateFunctionMat(fn,M,K);CHKERRQ(ierr);

      if (k1==0) {
        err_loc = tol;
        break;
      } else {
        ierr = MatDenseGetArray(K,&F);CHKERRQ(ierr);
        p1 = PetscAbsScalar(beta*F[m]);
        p2 = PetscAbsScalar(beta*F[m+1]*avnorm);
        ierr = MatDenseRestoreArray(K,&F);CHKERRQ(ierr);
        if (p1 > 10*p2) {
          err_loc = p2;
          xm = 1.0/(PetscReal)m;
        } else if (p1 > p2) {
          err_loc = (p1*p2)/(p1-p2);
          xm = 1.0/(PetscReal)m;
        } else {
          err_loc = p1;
          xm = 1.0/(PetscReal)(m-1);
        }
      }
      if (err_loc <= delta*t_step*tol) break;
      else {
        t_step = gamma*t_step*PetscPowReal(t_step*tol/err_loc,xm);
        s = PetscPowReal(10.0,PetscFloorReal(PetscLog10Real(t_step))-1);
        t_step = PetscCeilReal(t_step/s)*s;
        ireject = ireject+1;
      }
    }

    mx = mb + PetscMax(0,k1-1);
    ierr = MatDenseGetArray(K,&F);CHKERRQ(ierr);
    for (j=0;j<mx;j++) betaF[j] = beta*F[j];
    ierr = MatDenseRestoreArray(K,&F);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(mfn->V,0,mx);CHKERRQ(ierr);
    ierr = BVMultVec(mfn->V,1.0,0.0,x,betaF);CHKERRQ(ierr);
    ierr = VecNorm(x,NORM_2,&beta);CHKERRQ(ierr);

    t_now = t_now+t_step;
    if (t_now>=t_out) mfn->reason = MFN_CONVERGED_TOL;
    else {
      t_new = gamma*t_step*PetscPowReal((t_step*tol)/err_loc,xm);
      s = PetscPowReal(10.0,PetscFloorReal(PetscLog10Real(t_new))-1);
      t_new = PetscCeilReal(t_new/s)*s;
    }
    err_loc = PetscMax(err_loc,rndoff);
    if (mfn->its==mxstep) mfn->reason = MFN_DIVERGED_ITS;
    ierr = MFNMonitor(mfn,mfn->its,err_loc);CHKERRQ(ierr);
  }
  ierr = VecScale(x,sfactor);CHKERRQ(ierr);

  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  ierr = FNDestroy(&fn);CHKERRQ(ierr);
  ierr = PetscFree3(betaF,H,B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode MFNCreate_Expokit(MFN mfn)
{
  PetscFunctionBegin;
  mfn->ops->solve          = MFNSolve_Expokit;
  mfn->ops->setup          = MFNSetUp_Expokit;
  PetscFunctionReturn(0);
}
