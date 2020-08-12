/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Standard symmetric eigenproblem for the 3-D Laplacian built with the DM interface.\n\n"
"Use -seed <k> to modify the random initial vector.\n"
"Use -da_grid_x <nx> etc. to change the problem size.\n\n";

#include <slepceps.h>
#include <petscdmda.h>
#include <petsctime.h>

PetscErrorCode GetExactEigenvalues(PetscInt M,PetscInt N,PetscInt P,PetscInt nconv,PetscReal *exact)
{
  PetscInt       n,i,j,k,l;
  PetscReal      *evals,ax,ay,az,sx,sy,sz;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ax = PETSC_PI/2/(M+1);
  ay = PETSC_PI/2/(N+1);
  az = PETSC_PI/2/(P+1);
  n = PetscCeilReal(PetscPowReal((PetscReal)nconv,0.33333)+1);
  ierr = PetscMalloc1(n*n*n,&evals);CHKERRQ(ierr);
  l = 0;
  for (i=1;i<=n;i++) {
    sx = PetscSinReal(ax*i);
    for (j=1;j<=n;j++) {
      sy = PetscSinReal(ay*j);
      for (k=1;k<=n;k++) {
        sz = PetscSinReal(az*k);
        evals[l++] = 4.0*(sx*sx+sy*sy+sz*sz);
      }
    }
  }
  ierr = PetscSortReal(n*n*n,evals);CHKERRQ(ierr);
  for (i=0;i<nconv;i++) exact[i] = evals[i];
  ierr = PetscFree(evals);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FillMatrix(DM da,Mat A)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,idx;
  PetscScalar    v[7];
  MatStencil     row,col[7];

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  for (k=zs;k<zs+zm;k++) {
    for (j=ys;j<ys+ym;j++) {
      for (i=xs;i<xs+xm;i++) {
        row.i=i; row.j=j; row.k=k;
        col[0].i=row.i; col[0].j=row.j; col[0].k=row.k;
        v[0]=6.0;
        idx=1;
        if (k>0) { v[idx]=-1.0; col[idx].i=i; col[idx].j=j; col[idx].k=k-1; idx++; }
        if (j>0) { v[idx]=-1.0; col[idx].i=i; col[idx].j=j-1; col[idx].k=k; idx++; }
        if (i>0) { v[idx]=-1.0; col[idx].i=i-1; col[idx].j=j; col[idx].k=k; idx++; }
        if (i<mx-1) { v[idx]=-1.0; col[idx].i=i+1; col[idx].j=j; col[idx].k=k; idx++; }
        if (j<my-1) { v[idx]=-1.0; col[idx].i=i; col[idx].j=j+1; col[idx].k=k; idx++; }
        if (k<mz-1) { v[idx]=-1.0; col[idx].i=i; col[idx].j=j; col[idx].k=k+1; idx++; }
        ierr = MatSetValuesStencil(A,1,&row,idx,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  DM             da;
  Vec            v0;
  PetscReal      error,tol,re,im,*exact;
  PetscScalar    kr,ki;
  PetscInt       M,N,P,m,n,p,nev,maxit,i,its,nconv,seed;
  PetscLogDouble t1,t2,t3;
  PetscBool      flg,terse;
  PetscRandom    rctx;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n3-D Laplacian Eigenproblem\n\n");CHKERRQ(ierr);

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
                      DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,10,10,10,
                      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      1,1,NULL,NULL,NULL,&da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);

  /* print DM information */
  ierr = DMDAGetInfo(da,NULL,&M,&N,&P,&m,&n,&p,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Grid partitioning: %D %D %D\n",m,n,p);CHKERRQ(ierr);

  /* create and fill the matrix */
  ierr = DMCreateMatrix(da,&A);CHKERRQ(ierr);
  ierr = FillMatrix(da,A);CHKERRQ(ierr);

  /* create random initial vector */
  seed = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-seed",&seed,NULL);CHKERRQ(ierr);
  if (seed<0) SETERRQ(PETSC_COMM_WORLD,1,"Seed must be >=0");
  ierr = MatCreateVecs(A,&v0,NULL);CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
  for (i=0;i<seed;i++) {   /* simulate different seeds in the random generator */
    ierr = VecSetRandom(v0,rctx);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

  /*
     Set specific solver options
  */
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,1e-8,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = EPSSetInitialSpace(eps,1,&v0);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PetscTime(&t1);CHKERRQ(ierr);
  ierr = EPSSetUp(eps);CHKERRQ(ierr);
  ierr = PetscTime(&t2);CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscTime(&t3);CHKERRQ(ierr);
  if (!terse) {
    ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);

    /*
       Optional: Get some information from the solver and display it
    */
    ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
    ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  if (terse) {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    /*
       Get number of converged approximate eigenpairs
    */
    ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);

    if (nconv>0) {
      ierr = PetscMalloc1(nconv,&exact);CHKERRQ(ierr);
      ierr = GetExactEigenvalues(M,N,P,nconv,exact);CHKERRQ(ierr);
      /*
         Display eigenvalues and relative errors
      */
      ierr = PetscPrintf(PETSC_COMM_WORLD,
           "           k          ||Ax-kx||/||kx||   Eigenvalue Error \n"
           "   ----------------- ------------------ ------------------\n");CHKERRQ(ierr);

      for (i=0;i<nconv;i++) {
        /*
          Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
          ki (imaginary part)
        */
        ierr = EPSGetEigenpair(eps,i,&kr,&ki,NULL,NULL);CHKERRQ(ierr);
        /*
           Compute the relative error associated to each eigenpair
        */
        ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
        re = PetscRealPart(kr);
        im = PetscImaginaryPart(kr);
#else
        re = kr;
        im = ki;
#endif
        if (im!=0.0) SETERRQ(PETSC_COMM_WORLD,1,"Eigenvalue should be real");
        else {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12g       %12g        %12g\n",(double)re,(double)error,(double)PetscAbsReal(re-exact[i]));CHKERRQ(ierr);
        }
      }
      ierr = PetscFree(exact);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    }
  }

  /*
     Show computing times
  */
  ierr = PetscOptionsHasName(NULL,NULL,"-showtimes",&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD," Elapsed time: %g (setup), %g (solve)\n",(double)(t2-t1),(double)(t3-t2));CHKERRQ(ierr);
  }

  /*
     Free work space
  */
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&v0);CHKERRQ(ierr);
  ierr = PetscRandomDestroy(&rctx);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      args: -eps_nev 8 -terse
      requires: double
      output_file: output/ex19_1.out
      test:
         suffix: 1_krylovschur
         args: -eps_type krylovschur -eps_ncv 64
      test:
         suffix: 1_lobpcg
         args: -eps_type lobpcg -eps_tol 1e-7
      test:
         suffix: 1_blopex
         args: -eps_type blopex -eps_tol 1e-7 -eps_blopex_blocksize 4
         requires: blopex

TEST*/
