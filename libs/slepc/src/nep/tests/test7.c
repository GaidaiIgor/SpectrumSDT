/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test the NLEIGS solver with shell matrices.\n\n"
  "This is based on ex27.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = matrix dimension.\n"
  "  -split <0/1>, to select the split form in the problem definition (enabled by default)\n";

/*
   Solve T(lambda)x=0 using NLEIGS solver
      with T(lambda) = -D+sqrt(lambda)*I
      where D is the Laplacian operator in 1 dimension
      and with the interpolation interval [.01,16]
*/

#include <slepcnep.h>

/* User-defined routines */
PetscErrorCode FormFunction(NEP,PetscScalar,Mat,Mat,void*);
PetscErrorCode ComputeSingularities(NEP,PetscInt*,PetscScalar*,void*);
PetscErrorCode MatMult_A0(Mat,Vec,Vec);
PetscErrorCode MatGetDiagonal_A0(Mat,Vec);
PetscErrorCode MatDuplicate_A0(Mat,MatDuplicateOption,Mat*);
PetscErrorCode MatMult_A1(Mat,Vec,Vec);
PetscErrorCode MatGetDiagonal_A1(Mat,Vec);
PetscErrorCode MatDuplicate_A1(Mat,MatDuplicateOption,Mat*);
PetscErrorCode MatMult_F(Mat,Vec,Vec);
PetscErrorCode MatGetDiagonal_F(Mat,Vec);
PetscErrorCode MatDuplicate_F(Mat,MatDuplicateOption,Mat*);
PetscErrorCode MatDestroy_F(Mat);

typedef struct {
  PetscScalar t;  /* square root of lambda */
} MatCtx;

int main(int argc,char **argv)
{
  NEP            nep;
  KSP            *ksp;
  PC             pc;
  Mat            F,A[2];
  NEPType        type;
  PetscInt       i,n=100,nev,its,nsolve;
  PetscReal      keep,tol=PETSC_SQRT_MACHINE_EPSILON/10;
  PetscErrorCode ierr;
  RG             rg;
  FN             f[2];
  PetscMPIInt    size;
  PetscBool      terse,flg,lock,split=PETSC_TRUE;
  PetscScalar    coeffs;
  MatCtx         *ctx;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-split",&split,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSquare root eigenproblem, n=%D%s\n\n",n,split?" (in split form)":"");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create NEP context, configure NLEIGS with appropriate options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = NEPCreate(PETSC_COMM_WORLD,&nep);CHKERRQ(ierr);
  ierr = NEPSetType(nep,NEPNLEIGS);CHKERRQ(ierr);
  ierr = NEPNLEIGSSetSingularitiesFunction(nep,ComputeSingularities,NULL);CHKERRQ(ierr);
  ierr = NEPGetRG(nep,&rg);CHKERRQ(ierr);
  ierr = RGSetType(rg,RGINTERVAL);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = RGIntervalSetEndpoints(rg,0.01,16.0,-0.001,0.001);CHKERRQ(ierr);
#else
  ierr = RGIntervalSetEndpoints(rg,0.01,16.0,0,0);CHKERRQ(ierr);
#endif
  ierr = NEPSetTarget(nep,1.1);CHKERRQ(ierr);
  ierr = NEPNLEIGSGetKSPs(nep,&nsolve,&ksp);CHKERRQ(ierr);
  for (i=0;i<nsolve;i++) {
   ierr = KSPSetType(ksp[i],KSPBICG);CHKERRQ(ierr);
   ierr = KSPGetPC(ksp[i],&pc);CHKERRQ(ierr);
   ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
   ierr = KSPSetTolerances(ksp[i],tol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Define the nonlinear problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  if (split) {
    /* Create matrix A0 (tridiagonal) */
    ierr = MatCreateShell(PETSC_COMM_WORLD,n,n,n,n,NULL,&A[0]);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A[0],MATOP_MULT,(void(*)(void))MatMult_A0);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A[0],MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_A0);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A[0],MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_A0);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A[0],MATOP_DUPLICATE,(void(*)(void))MatDuplicate_A0);CHKERRQ(ierr);

    /* Create matrix A0 (identity) */
    ierr = MatCreateShell(PETSC_COMM_WORLD,n,n,n,n,NULL,&A[1]);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A[1],MATOP_MULT,(void(*)(void))MatMult_A1);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A[1],MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_A1);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A[1],MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_A1);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A[1],MATOP_DUPLICATE,(void(*)(void))MatDuplicate_A1);CHKERRQ(ierr);

    /* Define funcions for the split form */
    ierr = FNCreate(PETSC_COMM_WORLD,&f[0]);CHKERRQ(ierr);
    ierr = FNSetType(f[0],FNRATIONAL);CHKERRQ(ierr);
    coeffs = 1.0;
    ierr = FNRationalSetNumerator(f[0],1,&coeffs);CHKERRQ(ierr);
    ierr = FNCreate(PETSC_COMM_WORLD,&f[1]);CHKERRQ(ierr);
    ierr = FNSetType(f[1],FNSQRT);CHKERRQ(ierr);
    ierr = NEPSetSplitOperator(nep,2,A,f,SUBSET_NONZERO_PATTERN);CHKERRQ(ierr);
  } else {
    /* Callback form: create shell matrix for F=A0+sqrt(lambda)*A1  */
    ierr = PetscNew(&ctx);CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,n,n,n,n,(void*)ctx,&F);CHKERRQ(ierr);
    ierr = MatShellSetOperation(F,MATOP_MULT,(void(*)(void))MatMult_F);CHKERRQ(ierr);
    ierr = MatShellSetOperation(F,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_F);CHKERRQ(ierr);
    ierr = MatShellSetOperation(F,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_F);CHKERRQ(ierr);
    ierr = MatShellSetOperation(F,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_F);CHKERRQ(ierr);
    ierr = MatShellSetOperation(F,MATOP_DESTROY,(void(*)(void))MatDestroy_F);CHKERRQ(ierr);
    /* Set Function evaluation routine */
    ierr = NEPSetFunction(nep,F,F,FormFunction,NULL);CHKERRQ(ierr);
  }

  /* Set solver parameters at runtime */
  ierr = NEPSetFromOptions(nep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = NEPSolve(nep);CHKERRQ(ierr);
  ierr = NEPGetType(nep,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type);CHKERRQ(ierr);
  ierr = NEPGetDimensions(nep,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)nep,NEPNLEIGS,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = NEPNLEIGSGetRestart(nep,&keep);CHKERRQ(ierr);
    ierr = NEPNLEIGSGetLocking(nep,&lock);CHKERRQ(ierr);
    ierr = NEPNLEIGSGetInterpolation(nep,&tol,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Restart factor is %3.2f",(double)keep);CHKERRQ(ierr);
    if (lock) { ierr = PetscPrintf(PETSC_COMM_WORLD," (locking activated)");CHKERRQ(ierr); }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Divided diferences with tol=%6.2g maxit=%D\n",(double)tol,its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = NEPErrorView(nep,NEP_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = NEPReasonView(nep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = NEPErrorView(nep,NEP_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = NEPDestroy(&nep);CHKERRQ(ierr);
  if (split) {
    ierr = MatDestroy(&A[0]);CHKERRQ(ierr);
    ierr = MatDestroy(&A[1]);CHKERRQ(ierr);
    ierr = FNDestroy(&f[0]);CHKERRQ(ierr);
    ierr = FNDestroy(&f[1]);CHKERRQ(ierr);
  } else {
    ierr = MatDestroy(&F);CHKERRQ(ierr);
  }
  ierr = SlepcFinalize();
  return ierr;
}

/*
   FormFunction - Computes Function matrix  T(lambda)
*/
PetscErrorCode FormFunction(NEP nep,PetscScalar lambda,Mat fun,Mat B,void *ctx)
{
  PetscErrorCode ierr;
  MatCtx         *ctxF;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(fun,(void**)&ctxF);CHKERRQ(ierr);
  ctxF->t = PetscSqrtScalar(lambda);
  PetscFunctionReturn(0);
}

/*
   ComputeSingularities - Computes maxnp points (at most) in the complex plane where
   the function T(.) is not analytic.

   In this case, we discretize the singularity region (-inf,0)~(-10e+6,-10e-6)
*/
PetscErrorCode ComputeSingularities(NEP nep,PetscInt *maxnp,PetscScalar *xi,void *pt)
{
  PetscReal h;
  PetscInt  i;

  PetscFunctionBeginUser;
  h = 12.0/(*maxnp-1);
  xi[0] = -1e-6; xi[*maxnp-1] = -1e+6;
  for (i=1;i<*maxnp-1;i++) xi[i] = -PetscPowReal(10,-6+h*i);
  PetscFunctionReturn(0);
}

/* -------------------------------- A0 ----------------------------------- */

PetscErrorCode MatMult_A0(Mat A,Vec x,Vec y)
{
  PetscErrorCode    ierr;
  PetscInt          i,n;
  const PetscScalar *px;
  PetscScalar       *py;

  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecGetArray(y,&py);CHKERRQ(ierr);
  ierr = VecGetSize(x,&n);CHKERRQ(ierr);
  py[0] = -2.0*px[0]+px[1];
  for (i=1;i<n-1;i++) py[i] = px[i-1]-2.0*px[i]+px[i+1];
  py[n-1] = px[n-2]-2.0*px[n-1];
  ierr = VecRestoreArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatGetDiagonal_A0(Mat A,Vec diag)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecSet(diag,-2.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatDuplicate_A0(Mat A,MatDuplicateOption op,Mat *B)
{
  PetscInt       n;
  MPI_Comm       comm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatGetSize(A,&n,NULL);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
  ierr = MatCreateShell(comm,n,n,n,n,NULL,B);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_MULT,(void(*)(void))MatMult_A0);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_A0);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_A0);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_A0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* -------------------------------- A1 ----------------------------------- */

PetscErrorCode MatMult_A1(Mat A,Vec x,Vec y)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecCopy(x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatGetDiagonal_A1(Mat A,Vec diag)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecSet(diag,1.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatDuplicate_A1(Mat A,MatDuplicateOption op,Mat *B)
{
  PetscInt       n;
  MPI_Comm       comm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatGetSize(A,&n,NULL);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
  ierr = MatCreateShell(comm,n,n,n,n,NULL,B);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_MULT,(void(*)(void))MatMult_A1);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_A1);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_A1);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_A1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* -------------------------------- F ----------------------------------- */

PetscErrorCode MatMult_F(Mat A,Vec x,Vec y)
{
  PetscErrorCode    ierr;
  PetscInt          i,n;
  const PetscScalar *px;
  PetscScalar       *py,d;
  MatCtx            *ctx;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = VecGetArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecGetArray(y,&py);CHKERRQ(ierr);
  ierr = VecGetSize(x,&n);CHKERRQ(ierr);
  d = -2.0+ctx->t;
  py[0] = d*px[0]+px[1];
  for (i=1;i<n-1;i++) py[i] = px[i-1]+d*px[i]+px[i+1];
  py[n-1] = px[n-2]+d*px[n-1];
  ierr = VecRestoreArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatGetDiagonal_F(Mat A,Vec diag)
{
  PetscErrorCode ierr;
  MatCtx         *ctx;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = VecSet(diag,-2.0+ctx->t);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatDuplicate_F(Mat A,MatDuplicateOption op,Mat *B)
{
  MatCtx         *actx,*bctx;
  PetscInt       n;
  MPI_Comm       comm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,(void**)&actx);CHKERRQ(ierr);
  ierr = MatGetSize(A,&n,NULL);CHKERRQ(ierr);
  ierr = PetscNew(&bctx);CHKERRQ(ierr);
  bctx->t = actx->t;
  ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
  ierr = MatCreateShell(comm,n,n,n,n,(void*)bctx,B);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_MULT,(void(*)(void))MatMult_F);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_F);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_F);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_F);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*B,MATOP_DESTROY,(void(*)(void))MatDestroy_F);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatDestroy_F(Mat A)
{
  MatCtx         *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      args: -nep_nev 3 -nep_tol 1e-8 -nep_nleigs_locking 0 -nep_nleigs_interpolation_degree 90 -nep_nleigs_interpolation_tol 1e-8 -nep_nleigs_restart 0.4 -terse
      requires: !single

   test:
      suffix: 2
      args: -split 0 -nep_nev 3 -nep_tol 1e-8 -terse
      requires: !single

TEST*/
