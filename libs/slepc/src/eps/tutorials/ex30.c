/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates the use of a region for filtering; the number of wanted eigenvalues in not known a priori.\n\n"
  "The problem is the Brusselator wave model as in ex9.c.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = block dimension of the 2x2 block matrix.\n"
  "  -L <L>, where <L> = bifurcation parameter.\n"
  "  -alpha <alpha>, -beta <beta>, -delta1 <delta1>,  -delta2 <delta2>,\n"
  "       where <alpha> <beta> <delta1> <delta2> = model parameters.\n\n";

#include <slepceps.h>

/*
   This example tries to compute all eigenvalues lying outside the real axis.
   This could be achieved by computing LARGEST_IMAGINARY eigenvalues, but
   here we take a different route: define a region of the complex plane where
   eigenvalues must be emphasized (eigenvalues outside the region are filtered
   out). In this case, we select the region as the complement of a thin stripe
   around the real axis.
 */

PetscErrorCode MatMult_Brussel(Mat,Vec,Vec);
PetscErrorCode MatGetDiagonal_Brussel(Mat,Vec);
PetscErrorCode MyStoppingTest(EPS,PetscInt,PetscInt,PetscInt,PetscInt,EPSConvergedReason*,void*);

typedef struct {
  Mat         T;
  Vec         x1,x2,y1,y2;
  PetscScalar alpha,beta,tau1,tau2,sigma;
  PetscInt    lastnconv;      /* last value of nconv; used in stopping test */
  PetscInt    nreps;          /* number of repetitions of nconv; used in stopping test */
} CTX_BRUSSEL;

int main(int argc,char **argv)
{
  Mat            A;               /* eigenvalue problem matrix */
  EPS            eps;             /* eigenproblem solver context */
  RG             rg;              /* region object */
  PetscScalar    delta1,delta2,L,h;
  PetscInt       N=30,n,i,Istart,Iend,mpd;
  CTX_BRUSSEL    *ctx;
  PetscBool      terse;
  PetscViewer    viewer;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nBrusselator wave model, n=%D\n\n",N);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create shell matrix context and set default parameters
  */
  ierr = PetscNew(&ctx);CHKERRQ(ierr);
  ctx->alpha = 2.0;
  ctx->beta  = 5.45;
  delta1     = 0.008;
  delta2     = 0.004;
  L          = 0.51302;

  /*
     Look the command line for user-provided parameters
  */
  ierr = PetscOptionsGetScalar(NULL,NULL,"-L",&L,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-alpha",&ctx->alpha,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-beta",&ctx->beta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-delta1",&delta1,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-delta2",&delta2,NULL);CHKERRQ(ierr);

  /*
     Create matrix T
  */
  ierr = MatCreate(PETSC_COMM_WORLD,&ctx->T);CHKERRQ(ierr);
  ierr = MatSetSizes(ctx->T,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(ctx->T);CHKERRQ(ierr);
  ierr = MatSetUp(ctx->T);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(ctx->T,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(ctx->T,i,i-1,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<N-1) { ierr = MatSetValue(ctx->T,i,i+1,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(ctx->T,i,i,-2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(ctx->T,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(ctx->T,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatGetLocalSize(ctx->T,&n,NULL);CHKERRQ(ierr);

  /*
     Fill the remaining information in the shell matrix context
     and create auxiliary vectors
  */
  h = 1.0 / (PetscReal)(N+1);
  ctx->tau1 = delta1 / ((h*L)*(h*L));
  ctx->tau2 = delta2 / ((h*L)*(h*L));
  ctx->sigma = 0.0;
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->x1);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->x2);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->y1);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->y2);CHKERRQ(ierr);

  /*
     Create the shell matrix
  */
  ierr = MatCreateShell(PETSC_COMM_WORLD,2*n,2*n,2*N,2*N,(void*)ctx,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatMult_Brussel);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Brussel);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and configure the region
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr);

  /*
     Define the region containing the eigenvalues of interest
  */
  ierr = EPSGetRG(eps,&rg);CHKERRQ(ierr);
  ierr = RGSetType(rg,RGINTERVAL);CHKERRQ(ierr);
  ierr = RGIntervalSetEndpoints(rg,-PETSC_INFINITY,PETSC_INFINITY,-0.01,0.01);CHKERRQ(ierr);
  ierr = RGSetComplement(rg,PETSC_TRUE);CHKERRQ(ierr);
  /* sort eigenvalue approximations wrt a target, otherwise convergence will be erratic */
  ierr = EPSSetTarget(eps,0.0);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);CHKERRQ(ierr);

  /*
     Set solver options. In particular, we must allocate sufficient
     storage for all eigenpairs that may converge (ncv). This is
     application-dependent.
  */
  mpd = 40;
  ierr = EPSSetDimensions(eps,2*mpd,3*mpd,mpd);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,1e-7,2000);CHKERRQ(ierr);
  ctx->lastnconv = 0;
  ctx->nreps     = 0;
  ierr = EPSSetStoppingTestFunction(eps,MyStoppingTest,(void*)ctx,NULL);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Solve the eigensystem and display solution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /* show detailed info unless -terse option is given by user */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
  ierr = EPSReasonView(eps,viewer);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (!terse) {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->T);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->x1);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->x2);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->y1);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->y2);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

PetscErrorCode MatMult_Brussel(Mat A,Vec x,Vec y)
{
  PetscInt          n;
  const PetscScalar *px;
  PetscScalar       *py;
  CTX_BRUSSEL       *ctx;
  PetscErrorCode    ierr;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = MatGetLocalSize(ctx->T,&n,NULL);CHKERRQ(ierr);
  ierr = VecGetArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecGetArray(y,&py);CHKERRQ(ierr);
  ierr = VecPlaceArray(ctx->x1,px);CHKERRQ(ierr);
  ierr = VecPlaceArray(ctx->x2,px+n);CHKERRQ(ierr);
  ierr = VecPlaceArray(ctx->y1,py);CHKERRQ(ierr);
  ierr = VecPlaceArray(ctx->y2,py+n);CHKERRQ(ierr);

  ierr = MatMult(ctx->T,ctx->x1,ctx->y1);CHKERRQ(ierr);
  ierr = VecScale(ctx->y1,ctx->tau1);CHKERRQ(ierr);
  ierr = VecAXPY(ctx->y1,ctx->beta - 1.0 + ctx->sigma,ctx->x1);CHKERRQ(ierr);
  ierr = VecAXPY(ctx->y1,ctx->alpha * ctx->alpha,ctx->x2);CHKERRQ(ierr);

  ierr = MatMult(ctx->T,ctx->x2,ctx->y2);CHKERRQ(ierr);
  ierr = VecScale(ctx->y2,ctx->tau2);CHKERRQ(ierr);
  ierr = VecAXPY(ctx->y2,-ctx->beta,ctx->x1);CHKERRQ(ierr);
  ierr = VecAXPY(ctx->y2,-ctx->alpha * ctx->alpha + ctx->sigma,ctx->x2);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py);CHKERRQ(ierr);
  ierr = VecResetArray(ctx->x1);CHKERRQ(ierr);
  ierr = VecResetArray(ctx->x2);CHKERRQ(ierr);
  ierr = VecResetArray(ctx->y1);CHKERRQ(ierr);
  ierr = VecResetArray(ctx->y2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatGetDiagonal_Brussel(Mat A,Vec diag)
{
  Vec            d1,d2;
  PetscInt       n;
  PetscScalar    *pd;
  MPI_Comm       comm;
  CTX_BRUSSEL    *ctx;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,(void**)&ctx);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
  ierr = MatGetLocalSize(ctx->T,&n,NULL);CHKERRQ(ierr);
  ierr = VecGetArray(diag,&pd);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(comm,1,n,PETSC_DECIDE,pd,&d1);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(comm,1,n,PETSC_DECIDE,pd+n,&d2);CHKERRQ(ierr);

  ierr = VecSet(d1,-2.0*ctx->tau1 + ctx->beta - 1.0 + ctx->sigma);CHKERRQ(ierr);
  ierr = VecSet(d2,-2.0*ctx->tau2 - ctx->alpha*ctx->alpha + ctx->sigma);CHKERRQ(ierr);

  ierr = VecDestroy(&d1);CHKERRQ(ierr);
  ierr = VecDestroy(&d2);CHKERRQ(ierr);
  ierr = VecRestoreArray(diag,&pd);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    Function for user-defined stopping test.

    Ignores the value of nev. It only takes into account the number of
    eigenpairs that have converged in recent outer iterations (restarts);
    if no new eigenvalus have converged in the last few restarts,
    we stop the iteration, assuming that no more eigenvalues are present
    inside the region.
*/
PetscErrorCode MyStoppingTest(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ptr)
{
  PetscErrorCode ierr;
  CTX_BRUSSEL    *ctx = (CTX_BRUSSEL*)ptr;

  PetscFunctionBeginUser;
  /* check usual termination conditions, but ignoring the case nconv>=nev */
  ierr = EPSStoppingBasic(eps,its,max_it,nconv,PETSC_MAX_INT,reason,NULL);CHKERRQ(ierr);
  if (*reason==EPS_CONVERGED_ITERATING) {
    /* check if nconv is the same as before */
    if (nconv==ctx->lastnconv) ctx->nreps++;
    else {
      ctx->lastnconv = nconv;
      ctx->nreps     = 0;
    }
    /* check if no eigenvalues converged in last 10 restarts */
    if (nconv && ctx->nreps>10) *reason = EPS_CONVERGED_USER;
  }
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      args: -n 100 -terse
      requires: !single

TEST*/
