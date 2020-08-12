/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates how to obtain invariant subspaces. "
  "Based on ex9.\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = block dimension of the 2x2 block matrix.\n"
  "  -L <L>, where <L> = bifurcation parameter.\n"
  "  -alpha <alpha>, -beta <beta>, -delta1 <delta1>,  -delta2 <delta2>,\n"
  "       where <alpha> <beta> <delta1> <delta2> = model parameters.\n\n";

#include <slepceps.h>

/*
   This example computes the eigenvalues with largest real part of the
   following matrix

        A = [ tau1*T+(beta-1)*I     alpha^2*I
                  -beta*I        tau2*T-alpha^2*I ],

   where

        T = tridiag{1,-2,1}
        h = 1/(n+1)
        tau1 = delta1/(h*L)^2
        tau2 = delta2/(h*L)^2
 */

/* Matrix operations */
PetscErrorCode MatMult_Brussel(Mat,Vec,Vec);
PetscErrorCode MatGetDiagonal_Brussel(Mat,Vec);

typedef struct {
  Mat         T;
  Vec         x1,x2,y1,y2;
  PetscScalar alpha,beta,tau1,tau2,sigma;
} CTX_BRUSSEL;

int main(int argc,char **argv)
{
  EPS            eps;
  Mat            A;
  Vec            *Q,v;
  PetscScalar    delta1,delta2,L,h,kr,ki;
  PetscReal      errest,tol,re,im,lev;
  PetscInt       N=30,n,i,j,Istart,Iend,nev,nconv;
  CTX_BRUSSEL    *ctx;
  PetscBool      errok,trueres;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nBrusselator wave model, n=%D\n\n",N);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PetscNew(&ctx);CHKERRQ(ierr);
  ctx->alpha = 2.0;
  ctx->beta  = 5.45;
  delta1     = 0.008;
  delta2     = 0.004;
  L          = 0.51302;

  ierr = PetscOptionsGetScalar(NULL,NULL,"-L",&L,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-alpha",&ctx->alpha,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-beta",&ctx->beta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-delta1",&delta1,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-delta2",&delta2,NULL);CHKERRQ(ierr);

  /* Create matrix T */
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

  /* Fill the remaining information in the shell matrix context
     and create auxiliary vectors */
  h = 1.0 / (PetscReal)(N+1);
  ctx->tau1 = delta1 / ((h*L)*(h*L));
  ctx->tau2 = delta2 / ((h*L)*(h*L));
  ctx->sigma = 0.0;
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->x1);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->x2);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->y1);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,n,PETSC_DECIDE,NULL,&ctx->y2);CHKERRQ(ierr);

  /* Create the shell matrix */
  ierr = MatCreateShell(PETSC_COMM_WORLD,2*n,2*n,2*N,2*N,(void*)ctx,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatMult_Brussel);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Brussel);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL);CHKERRQ(ierr);
  ierr = EPSSetTrueResidual(eps,PETSC_FALSE);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);

  ierr = EPSGetTrueResidual(eps,&trueres);CHKERRQ(ierr);
  /*if (trueres) { ierr = PetscPrintf(PETSC_COMM_WORLD," Computing true residuals explicitly\n\n");CHKERRQ(ierr); }*/

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,NULL);CHKERRQ(ierr);
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  if (nconv<nev) {
    ierr = PetscPrintf(PETSC_COMM_WORLD," Problem: less than %D eigenvalues converged\n\n",nev);CHKERRQ(ierr);
  } else {
    /* Check that all converged eigenpairs satisfy the requested tolerance
       (in this example we use the solver's error estimate instead of computing
       the residual norm explicitly) */
    errok = PETSC_TRUE;
    for (i=0;i<nev;i++) {
      ierr = EPSGetErrorEstimate(eps,i,&errest);CHKERRQ(ierr);
      errok = (errok && errest<5.0*tol)? PETSC_TRUE: PETSC_FALSE;
    }
    if (!errok) {
      ierr = PetscPrintf(PETSC_COMM_WORLD," Problem: some of the first %D relative errors are higher than the tolerance\n\n",nev);CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD," All requested eigenvalues computed up to the required tolerance:");CHKERRQ(ierr);
      for (i=0;i<=(nev-1)/8;i++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n     ");CHKERRQ(ierr);
        for (j=0;j<PetscMin(8,nev-8*i);j++) {
          ierr = EPSGetEigenpair(eps,8*i+j,&kr,&ki,NULL,NULL);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
          re = PetscRealPart(kr);
          im = PetscImaginaryPart(kr);
#else
          re = kr;
          im = ki;
#endif
          if (PetscAbs(re)/PetscAbs(im)<PETSC_SMALL) re = 0.0;
          if (PetscAbs(im)/PetscAbs(re)<PETSC_SMALL) im = 0.0;
          if (im!=0.0) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"%.5f%+.5fi",(double)re,(double)im);CHKERRQ(ierr);
          } else {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"%.5f",(double)re);CHKERRQ(ierr);
          }
          if (8*i+j+1<nev) { ierr = PetscPrintf(PETSC_COMM_WORLD,", ");CHKERRQ(ierr); }
        }
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n");CHKERRQ(ierr);
    }
  }

  /* Get an orthogonal basis of the invariant subspace and check it is indeed
     orthogonal (note that eigenvectors are not orthogonal in this case) */
  if (nconv>1) {
    ierr = MatCreateVecs(A,&v,NULL);CHKERRQ(ierr);
    ierr = VecDuplicateVecs(v,nconv,&Q);CHKERRQ(ierr);
    ierr = EPSGetInvariantSubspace(eps,Q);CHKERRQ(ierr);
    ierr = VecCheckOrthogonality(Q,nconv,NULL,nconv,NULL,NULL,&lev);CHKERRQ(ierr);
    if (lev<10*tol) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality below the tolerance\n");CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Level of orthogonality: %g\n",(double)lev);CHKERRQ(ierr);
    }
    ierr = VecDestroyVecs(nconv,&Q);CHKERRQ(ierr);
    ierr = VecDestroy(&v);CHKERRQ(ierr);
  }

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

/*TEST

   test:
      suffix: 1
      args: -eps_nev 4 -eps_true_residual {{0 1}}
      requires: !single
      output_file: output/test22_1.out

   test:
      suffix: 2
      args: -eps_nev 4 -eps_true_residual -eps_balance oneside -eps_tol 1e-7
      requires: !single

TEST*/
