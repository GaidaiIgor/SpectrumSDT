/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates the computation of left eigenvectors.\n\n"
  "The problem is the Markov model as in ex5.c.\n"
  "The command line options are:\n"
  "  -m <m>, where <m> = number of grid subdivisions in each dimension.\n\n";

#include <slepceps.h>

/*
   User-defined routines
*/
PetscErrorCode MatMarkovModel(PetscInt,Mat);
PetscErrorCode ComputeResidualNorm(Mat,PetscBool,PetscScalar,PetscScalar,Vec,Vec,Vec,PetscReal*);

int main(int argc,char **argv)
{
  Vec            v0,w0;           /* initial vectors */
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscInt       i,N,m=15,nconv;
  PetscBool      twosided;
  PetscReal      nrmr,nrml=0.0,re,im,lev;
  PetscScalar    *kr,*ki;
  Vec            t,*xr,*xi,*yr,*yi;
  PetscMPIInt    rank;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  N = m*(m+1)/2;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nMarkov Model, N=%D (m=%D)\n\n",N,m);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = MatMarkovModel(m,A);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&t);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr);

  /* use a two-sided algorithm to compute left eigenvectors as well */
  ierr = EPSSetTwoSided(eps,PETSC_TRUE);CHKERRQ(ierr);

  /* allow user to change settings at run time */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  ierr = EPSGetTwoSided(eps,&twosided);CHKERRQ(ierr);

  /*
     Set the initial vectors. This is optional, if not done the initial
     vectors are set to random values
  */
  ierr = MatCreateVecs(A,&v0,&w0);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  if (!rank) {
    ierr = VecSetValue(v0,0,1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(v0,1,1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(v0,2,1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(w0,0,2.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(w0,2,0.5,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(v0);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(w0);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(v0);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(w0);CHKERRQ(ierr);
  ierr = EPSSetInitialSpace(eps,1,&v0);CHKERRQ(ierr);
  ierr = EPSSetLeftInitialSpace(eps,1,&w0);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Get number of converged approximate eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);
  ierr = PetscMalloc2(nconv,&kr,nconv,&ki);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(t,nconv,&xr);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(t,nconv,&xi);CHKERRQ(ierr);
  if (twosided) {
    ierr = VecDuplicateVecs(t,nconv,&yr);CHKERRQ(ierr);
    ierr = VecDuplicateVecs(t,nconv,&yi);CHKERRQ(ierr);
  }

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k            ||Ax-kx||         ||y'A-y'k||\n"
         "   ---------------- ------------------ ------------------\n");CHKERRQ(ierr);

    for (i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,i,&kr[i],&ki[i],xr[i],xi[i]);CHKERRQ(ierr);
      if (twosided) {
        ierr = EPSGetLeftEigenvector(eps,i,yr[i],yi[i]);CHKERRQ(ierr);
      }
      /*
         Compute the residual norms associated to each eigenpair
      */
      ierr = ComputeResidualNorm(A,PETSC_FALSE,kr[i],ki[i],xr[i],xi[i],t,&nrmr);CHKERRQ(ierr);
      if (twosided) {
        ierr = ComputeResidualNorm(A,PETSC_TRUE,kr[i],ki[i],yr[i],yi[i],t,&nrml);CHKERRQ(ierr);
      }

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr[i]);
      im = PetscImaginaryPart(kr[i]);
#else
      re = kr[i];
      im = ki[i];
#endif
      if (im!=0.0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," %8f%+8fi %12g %12g\n",(double)re,(double)im,(double)nrmr,(double)nrml);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g       %12g\n",(double)re,(double)nrmr,(double)nrml);CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    /*
       Check bi-orthogonality of eigenvectors
    */
    if (twosided) {
      ierr = VecCheckOrthogonality(xr,nconv,yr,nconv,NULL,NULL,&lev);CHKERRQ(ierr);
      if (lev<100*PETSC_MACHINE_EPSILON) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of bi-orthogonality of eigenvectors < 100*eps\n\n");CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of bi-orthogonality of eigenvectors: %g\n\n",(double)lev);CHKERRQ(ierr);
      }
    }
  }

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&v0);CHKERRQ(ierr);
  ierr = VecDestroy(&w0);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = PetscFree2(kr,ki);CHKERRQ(ierr);
  ierr = VecDestroyVecs(nconv,&xr);CHKERRQ(ierr);
  ierr = VecDestroyVecs(nconv,&xi);CHKERRQ(ierr);
  if (twosided) {
    ierr = VecDestroyVecs(nconv,&yr);CHKERRQ(ierr);
    ierr = VecDestroyVecs(nconv,&yi);CHKERRQ(ierr);
  }
  ierr = SlepcFinalize();
  return ierr;
}

/*
    Matrix generator for a Markov model of a random walk on a triangular grid.

    This subroutine generates a test matrix that models a random walk on a
    triangular grid. This test example was used by G. W. Stewart ["{SRRIT} - a
    FORTRAN subroutine to calculate the dominant invariant subspaces of a real
    matrix", Tech. report. TR-514, University of Maryland (1978).] and in a few
    papers on eigenvalue problems by Y. Saad [see e.g. LAA, vol. 34, pp. 269-295
    (1980) ]. These matrices provide reasonably easy test problems for eigenvalue
    algorithms. The transpose of the matrix  is stochastic and so it is known
    that one is an exact eigenvalue. One seeks the eigenvector of the transpose
    associated with the eigenvalue unity. The problem is to calculate the steady
    state probability distribution of the system, which is the eigevector
    associated with the eigenvalue one and scaled in such a way that the sum all
    the components is equal to one.

    Note: the code will actually compute the transpose of the stochastic matrix
    that contains the transition probabilities.
*/
PetscErrorCode MatMarkovModel(PetscInt m,Mat A)
{
  const PetscReal cst = 0.5/(PetscReal)(m-1);
  PetscReal       pd,pu;
  PetscInt        Istart,Iend,i,j,jmax,ix=0;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=1;i<=m;i++) {
    jmax = m-i+1;
    for (j=1;j<=jmax;j++) {
      ix = ix + 1;
      if (ix-1<Istart || ix>Iend) continue;  /* compute only owned rows */
      if (j!=jmax) {
        pd = cst*(PetscReal)(i+j-1);
        /* north */
        if (i==1) {
          ierr = MatSetValue(A,ix-1,ix,2*pd,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          ierr = MatSetValue(A,ix-1,ix,pd,INSERT_VALUES);CHKERRQ(ierr);
        }
        /* east */
        if (j==1) {
          ierr = MatSetValue(A,ix-1,ix+jmax-1,2*pd,INSERT_VALUES);CHKERRQ(ierr);
        } else {
          ierr = MatSetValue(A,ix-1,ix+jmax-1,pd,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
      /* south */
      pu = 0.5 - cst*(PetscReal)(i+j-3);
      if (j>1) {
        ierr = MatSetValue(A,ix-1,ix-2,pu,INSERT_VALUES);CHKERRQ(ierr);
      }
      /* west */
      if (i>1) {
        ierr = MatSetValue(A,ix-1,ix-jmax-2,pu,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   ComputeResidualNorm - Computes the norm of the residual vector
   associated with an eigenpair.

   Input Parameters:
     trans - whether A' must be used instead of A
     kr,ki - eigenvalue
     xr,xi - eigenvector
     u     - work vector
*/
PetscErrorCode ComputeResidualNorm(Mat A,PetscBool trans,PetscScalar kr,PetscScalar ki,Vec xr,Vec xi,Vec u,PetscReal *norm)
{
  PetscErrorCode ierr;
#if !defined(PETSC_USE_COMPLEX)
  PetscReal      ni,nr;
#endif
  PetscErrorCode (*matmult)(Mat,Vec,Vec) = trans? MatMultTranspose: MatMult;

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  if (ki == 0 || PetscAbsScalar(ki) < PetscAbsScalar(kr*PETSC_MACHINE_EPSILON)) {
#endif
    ierr = (*matmult)(A,xr,u);CHKERRQ(ierr);
    if (PetscAbsScalar(kr) > PETSC_MACHINE_EPSILON) {
      ierr = VecAXPY(u,-kr,xr);CHKERRQ(ierr);
    }
    ierr = VecNorm(u,NORM_2,norm);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  } else {
    ierr = (*matmult)(A,xr,u);CHKERRQ(ierr);
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      ierr = VecAXPY(u,-kr,xr);CHKERRQ(ierr);
      ierr = VecAXPY(u,ki,xi);CHKERRQ(ierr);
    }
    ierr = VecNorm(u,NORM_2,&nr);CHKERRQ(ierr);
    ierr = (*matmult)(A,xi,u);CHKERRQ(ierr);
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      ierr = VecAXPY(u,-kr,xi);CHKERRQ(ierr);
      ierr = VecAXPY(u,-ki,xr);CHKERRQ(ierr);
    }
    ierr = VecNorm(u,NORM_2,&ni);CHKERRQ(ierr);
    *norm = SlepcAbsEigenvalue(nr,ni);
  }
#endif
  PetscFunctionReturn(0);
}

/*TEST

   testset:
      args: -st_type sinvert -eps_target 1.1 -eps_nev 4
      filter: grep -v method | sed -e "s/[+-]0.0*i//" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: !single
      output_file: output/ex41_1.out
      test:
         suffix: 1
         args: -eps_type {{power krylovschur}}
      test:
         suffix: 1_balance
         args: -eps_balance {{oneside twoside}} -eps_ncv 18 -eps_krylovschur_locking 0

TEST*/
