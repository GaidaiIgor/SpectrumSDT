/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Spectrum folding for a standard symmetric eigenproblem.\n\n"
  "The problem matrix is the 2-D Laplacian.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n";

#include <slepceps.h>

/*
   User context for spectrum folding
*/
typedef struct {
  Mat       A;
  Vec       w;
  PetscReal target;
} CTX_FOLD;

/*
   Auxiliary routines
*/
PetscErrorCode MatMult_Fold(Mat,Vec,Vec);
PetscErrorCode RayleighQuotient(Mat,Vec,PetscScalar*);
PetscErrorCode ComputeResidualNorm(Mat,PetscScalar,Vec,PetscReal*);

int main(int argc,char **argv)
{
  Mat            A,M,P;       /* problem matrix, shell matrix and preconditioner */
  Vec            x;           /* eigenvector */
  EPS            eps;         /* eigenproblem solver context */
  ST             st;          /* spectral transformation */
  KSP            ksp;
  PC             pc;
  EPSType        type;
  CTX_FOLD       *ctx;
  PetscInt       nconv,N,n=10,m,nloc,mloc,Istart,Iend,II,i,j;
  PetscReal      *error,*evals,target=0.0,tol;
  PetscScalar    lambda;
  PetscBool      flag,terse,errok,hasmat;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag);CHKERRQ(ierr);
  if (!flag) m=n;
  ierr = PetscOptionsGetReal(NULL,NULL,"-target",&target,NULL);CHKERRQ(ierr);
  N = n*m;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSpectrum Folding, N=%D (%Dx%D grid) target=%f\n\n",N,n,m,(double)target);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the 5-point stencil Laplacian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) { ierr = MatSetValue(A,II,II-n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<m-1) { ierr = MatSetValue(A,II,II+n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j>0) { ierr = MatSetValue(A,II,II-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(A,II,II+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,II,II,4.0,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,&x,NULL);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&nloc,&mloc);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create shell matrix to perform spectrum folding
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscNew(&ctx);CHKERRQ(ierr);
  ctx->A = A;
  ctx->target = target;
  ierr = VecDuplicate(x,&ctx->w);CHKERRQ(ierr);

  ierr = MatCreateShell(PETSC_COMM_WORLD,nloc,mloc,N,N,ctx,&M);CHKERRQ(ierr);
  ierr = MatShellSetOperation(M,MATOP_MULT,(void(*)(void))MatMult_Fold);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,M,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  ierr = PetscObjectTypeCompareAny((PetscObject)eps,&flag,EPSGD,EPSJD,EPSBLOPEX,EPSLOBPCG,EPSRQCG,"");CHKERRQ(ierr);
  if (flag) {
    /*
       Build preconditioner specific for this application (diagonal of A^2)
    */
    ierr = MatGetRowSum(A,x);CHKERRQ(ierr);
    ierr = VecScale(x,-1.0);CHKERRQ(ierr);
    ierr = VecShift(x,20.0);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&P);CHKERRQ(ierr);
    ierr = MatSetSizes(P,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
    ierr = MatSetFromOptions(P);CHKERRQ(ierr);
    ierr = MatSetUp(P);CHKERRQ(ierr);
    ierr = MatDiagonalSet(P,x,INSERT_VALUES);CHKERRQ(ierr);
    /*
       Set diagonal preconditioner
    */
    ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
    ierr = STSetType(st,STPRECOND);CHKERRQ(ierr);
    ierr = STPrecondSetMatForPC(st,P);CHKERRQ(ierr);
    ierr = MatDestroy(&P);CHKERRQ(ierr);
    ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    ierr = STPrecondGetKSPHasMat(st,&hasmat);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Preconditioned solver, hasmat=%s\n",hasmat?"true":"false");CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);
  if (nconv>0) {
    ierr = PetscMalloc2(nconv,&evals,nconv,&error);CHKERRQ(ierr);
    for (i=0;i<nconv;i++) {
      /*  Get i-th eigenvector, compute eigenvalue approximation from
          Rayleigh quotient and compute residual norm */
      ierr = EPSGetEigenpair(eps,i,NULL,NULL,x,NULL);CHKERRQ(ierr);
      ierr = RayleighQuotient(A,x,&lambda);CHKERRQ(ierr);
      ierr = ComputeResidualNorm(A,lambda,x,&error[i]);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
      evals[i] = PetscRealPart(lambda);
#else
      evals[i] = lambda;
#endif
    }
    ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
    if (!terse) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
           "           k              ||Ax-kx||\n"
           "   ----------------- ------------------\n");CHKERRQ(ierr);
      for (i=0;i<nconv;i++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12.2g\n",(double)evals[i],(double)error[i]);CHKERRQ(ierr);
      }
    } else {
      errok = PETSC_TRUE;
      for (i=0;i<nconv;i++) errok = (errok && error[i]<5.0*tol)? PETSC_TRUE: PETSC_FALSE;
      if (!errok) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," Problem: some of the first %D relative errors are higher than the tolerance\n\n",nconv);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD," nconv=%D eigenvalues computed up to the required tolerance:",nconv);CHKERRQ(ierr);
        for (i=0;i<nconv;i++) {
          ierr = PetscPrintf(PETSC_COMM_WORLD," %.5f",(double)evals[i]);CHKERRQ(ierr);
        }
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    ierr = PetscFree2(evals,error);CHKERRQ(ierr);
  }

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->w);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*
    Matrix-vector product subroutine for the spectrum folding.
       y <-- (A-t*I)^2*x
 */
PetscErrorCode MatMult_Fold(Mat M,Vec x,Vec y)
{
  CTX_FOLD       *ctx;
  PetscScalar    sigma;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(M,&ctx);CHKERRQ(ierr);
  sigma = -ctx->target;
  ierr = MatMult(ctx->A,x,ctx->w);CHKERRQ(ierr);
  ierr = VecAXPY(ctx->w,sigma,x);CHKERRQ(ierr);
  ierr = MatMult(ctx->A,ctx->w,y);CHKERRQ(ierr);
  ierr = VecAXPY(y,sigma,ctx->w);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    Computes the Rayleigh quotient of a vector x
       r <-- x^T*A*x       (assumes x has unit norm)
 */
PetscErrorCode RayleighQuotient(Mat A,Vec x,PetscScalar *r)
{
  Vec            Ax;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecDuplicate(x,&Ax);CHKERRQ(ierr);
  ierr = MatMult(A,x,Ax);CHKERRQ(ierr);
  ierr = VecDot(Ax,x,r);CHKERRQ(ierr);
  ierr = VecDestroy(&Ax);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    Computes the residual norm of an approximate eigenvector x, |A*x-lambda*x|
 */
PetscErrorCode ComputeResidualNorm(Mat A,PetscScalar lambda,Vec x,PetscReal *r)
{
  Vec            Ax;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecDuplicate(x,&Ax);CHKERRQ(ierr);
  ierr = MatMult(A,x,Ax);CHKERRQ(ierr);
  ierr = VecAXPY(Ax,-lambda,x);CHKERRQ(ierr);
  ierr = VecNorm(Ax,NORM_2,r);CHKERRQ(ierr);
  ierr = VecDestroy(&Ax);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   testset:
      args: -n 15 -eps_nev 1 -eps_ncv 12 -eps_max_it 1000 -eps_tol 1e-5 -terse
      filter: grep -v Solution
      test:
         suffix: 1
      test:
         suffix: 1_lobpcg
         args: -eps_type lobpcg
         requires: !single
      test:
         suffix: 1_gd
         args: -eps_type gd
         requires: !single

TEST*/
