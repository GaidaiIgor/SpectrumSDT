/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Computes exp(t*A)*v for an advection diffusion operator with Peclet number.\n\n"
  "The command line options are:\n"
  "  -n <idim>, where <idim> = number of subdivisions of the mesh in each spatial direction.\n"
  "  -t <sval>, where <sval> = scalar value that multiplies the argument.\n"
  "  -peclet <sval>, where <sval> = Peclet value.\n"
  "  -steps <ival>, where <ival> = number of time steps.\n\n";

#include <slepcmfn.h>

int main(int argc,char **argv)
{
  Mat                A;           /* problem matrix */
  MFN                mfn;
  FN                 f;
  PetscInt           i,j,Istart,Iend,II,m,n=10,N,steps=5,its,totits=0,ncv,maxit;
  PetscReal          tol,norm,h,h2,peclet=0.5,epsilon=1.0,c,i1h,j1h;
  PetscScalar        t=1e-4,sone=1.0,value,upper,diag,lower;
  Vec                v;
  PetscErrorCode     ierr;
  MFNConvergedReason reason;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetScalar(NULL,NULL,"-t",&t,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-peclet",&peclet,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-steps",&steps,NULL);CHKERRQ(ierr);
  m = n;
  N = m*n;
  /* interval [0,1], homogeneous Dirichlet boundary conditions */
  h = 1.0/(n+1.0);
  h2 = h*h;
  c = 2.0*epsilon*peclet/h;
  upper = epsilon/h2+c/(2.0*h);
  diag = 2.0*(-2.0*epsilon/h2);
  lower = epsilon/h2-c/(2.0*h);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nAdvection diffusion via y=exp(%g*A), n=%D, steps=%D, Peclet=%g\n\n",(double)PetscRealPart(t),n,steps,(double)peclet);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Generate matrix A
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) { ierr = MatSetValue(A,II,II-n,lower,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<m-1) { ierr = MatSetValue(A,II,II+n,upper,INSERT_VALUES);CHKERRQ(ierr); }
    if (j>0) { ierr = MatSetValue(A,II,II-1,lower,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(A,II,II+1,upper,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,II,II,diag,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&v);CHKERRQ(ierr);

  /*
     Set initial condition v = 256*i^2*(1-i)^2*j^2*(1-j)^2
  */
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    i1h = (i+1)*h; j1h = (j+1)*h;
    value = 256.0*i1h*i1h*(1.0-i1h)*(1.0-i1h)*(j1h*j1h)*(1.0-j1h)*(1.0-j1h);
    ierr = VecSetValue(v,i+j*n,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(v);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MFNCreate(PETSC_COMM_WORLD,&mfn);CHKERRQ(ierr);
  ierr = MFNSetOperator(mfn,A);CHKERRQ(ierr);
  ierr = MFNGetFN(mfn,&f);CHKERRQ(ierr);
  ierr = FNSetType(f,FNEXP);CHKERRQ(ierr);
  ierr = FNSetScale(f,t,sone);CHKERRQ(ierr);
  ierr = MFNSetFromOptions(mfn);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the problem, y=exp(t*A)*v
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for (i=0;i<steps;i++) {
    ierr = MFNSolve(mfn,v,v);CHKERRQ(ierr);
    ierr = MFNGetConvergedReason(mfn,&reason);CHKERRQ(ierr);
    if (reason<0) SETERRQ(PETSC_COMM_WORLD,1,"Solver did not converge");
    ierr = MFNGetIterationNumber(mfn,&its);CHKERRQ(ierr);
    totits += its;
  }

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",totits);CHKERRQ(ierr);
  ierr = MFNGetDimensions(mfn,&ncv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Subspace dimension: %D\n",ncv);CHKERRQ(ierr);
  ierr = MFNGetTolerances(mfn,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);
  ierr = VecNorm(v,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Computed vector at time t=%.4g has norm %g\n",(double)PetscRealPart(t)*steps,(double)norm);CHKERRQ(ierr);

  /*
     Free work space
  */
  ierr = MFNDestroy(&mfn);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -mfn_tol 1e-6

TEST*/
