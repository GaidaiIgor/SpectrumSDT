/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Tests solving an eigenproblem defined with MatNest. "
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

int main(int argc,char **argv)
{
  EPS            eps;
  Mat            A,T1,T2,D1,D2,mats[4];
  PetscScalar    alpha,beta,tau1,tau2,delta1,delta2,L,h;
  PetscInt       N=30,i,Istart,Iend;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nBrusselator wave model, n=%D\n\n",N);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Generate the matrix
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  alpha  = 2.0;
  beta   = 5.45;
  delta1 = 0.008;
  delta2 = 0.004;
  L      = 0.51302;

  ierr = PetscOptionsGetScalar(NULL,NULL,"-L",&L,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-alpha",&alpha,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-beta",&beta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-delta1",&delta1,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-delta2",&delta2,NULL);CHKERRQ(ierr);

  h    = 1.0 / (PetscReal)(N+1);
  tau1 = delta1 / ((h*L)*(h*L));
  tau2 = delta2 / ((h*L)*(h*L));

  /* Create matrices T1, T2 */
  ierr = MatCreate(PETSC_COMM_WORLD,&T1);CHKERRQ(ierr);
  ierr = MatSetSizes(T1,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(T1);CHKERRQ(ierr);
  ierr = MatSetUp(T1);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(T1,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(T1,i,i-1,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<N-1) { ierr = MatSetValue(T1,i,i+1,1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(T1,i,i,-2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(T1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(T1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatDuplicate(T1,MAT_COPY_VALUES,&T2);CHKERRQ(ierr);
  ierr = MatScale(T1,tau1);CHKERRQ(ierr);
  ierr = MatShift(T1,beta-1.0);CHKERRQ(ierr);
  ierr = MatScale(T2,tau2);CHKERRQ(ierr);
  ierr = MatShift(T2,-alpha*alpha);CHKERRQ(ierr);

  /* Create matrices D1, D2 */
  ierr = MatCreate(PETSC_COMM_WORLD,&D1);CHKERRQ(ierr);
  ierr = MatSetSizes(D1,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(D1);CHKERRQ(ierr);
  ierr = MatSetUp(D1);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(D1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(D1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatDuplicate(D1,MAT_COPY_VALUES,&D2);CHKERRQ(ierr);
  ierr = MatShift(D1,alpha*alpha);CHKERRQ(ierr);
  ierr = MatShift(D2,-beta);CHKERRQ(ierr);

  /* Create the nest matrix */
  mats[0] = T1;
  mats[1] = D1;
  mats[2] = D2;
  mats[3] = T2;
  ierr = MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,mats,&A);CHKERRQ(ierr);

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

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&T1);CHKERRQ(ierr);
  ierr = MatDestroy(&T2);CHKERRQ(ierr);
  ierr = MatDestroy(&D1);CHKERRQ(ierr);
  ierr = MatDestroy(&D2);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -eps_nev 4
      requires: !single
      filter: sed -e "s/0.00019-2.13938i, 0.00019+2.13938i/0.00019+2.13938i, 0.00019-2.13938i/" | sed -e "s/-0.67192-2.52712i, -0.67192+2.52712i/-0.67192+2.52712i, -0.67192-2.52712i/"

TEST*/
