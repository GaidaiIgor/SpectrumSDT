/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This example implements one of the problems found at
       NLEVP: A Collection of Nonlinear Eigenvalue Problems,
       The University of Manchester.
   The details of the collection can be found at:
       [1] T. Betcke et al., "NLEVP: A Collection of Nonlinear Eigenvalue
           Problems", ACM Trans. Math. Software 39(2), Article 7, 2013.

   The planar_waveguide problem is a quartic PEP with symmetric matrices,
   arising from a finite element solution of the propagation constants in a
   six-layer planar waveguide.
*/

static char help[] = "FEM solution of the propagation constants in a six-layer planar waveguide.\n\n"
  "The command line options are:\n"
  "  -n <n>, the dimension of the matrices.\n\n";

#include <slepcpep.h>

#define NMAT 5
#define NL   6

int main(int argc,char **argv)
{
  Mat            A[NMAT];         /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PetscInt       n=128,nlocal,k,Istart,Iend,i,j,start_ct,end_ct;
  PetscReal      w=9.92918,a=0.0,b=2.0,h,deltasq;
  PetscReal      nref[NL],K2[NL],q[NL],*md,*supd,*subd;
  PetscScalar    v,alpha;
  PetscBool      terse;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  n = (n/4)*4;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPlanar waveguide, n=%D\n\n",n+1);CHKERRQ(ierr);
  h = (b-a)/n;
  nlocal = (n/4)-1;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Set waveguide parameters used in construction of matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* refractive indices in each layer */
  nref[0] = 1.5;
  nref[1] = 1.66;
  nref[2] = 1.6;
  nref[3] = 1.53;
  nref[4] = 1.66;
  nref[5] = 1.0;

  for (i=0;i<NL;i++) K2[i] = w*w*nref[i]*nref[i];
  deltasq = K2[0] - K2[NL-1];
  for (i=0;i<NL;i++) q[i] = K2[i] - (K2[0] + K2[NL-1])/2;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                     Compute the polynomial matrices
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* initialize matrices */
  for (i=0;i<NMAT;i++) {
    ierr = MatCreate(PETSC_COMM_WORLD,&A[i]);CHKERRQ(ierr);
    ierr = MatSetSizes(A[i],PETSC_DECIDE,PETSC_DECIDE,n+1,n+1);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A[i]);CHKERRQ(ierr);
    ierr = MatSetUp(A[i]);CHKERRQ(ierr);
  }
  ierr = MatGetOwnershipRange(A[0],&Istart,&Iend);CHKERRQ(ierr);

  /* A0 */
  alpha = (h/6)*(deltasq*deltasq/16);
  for (i=Istart;i<Iend;i++) {
    v = 4.0;
    if (i==0 || i==n) v = 2.0;
    ierr = MatSetValue(A[0],i,i,v*alpha,INSERT_VALUES);CHKERRQ(ierr);
    if (i>0) { ierr = MatSetValue(A[0],i,i-1,alpha,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n) { ierr = MatSetValue(A[0],i,i+1,alpha,INSERT_VALUES);CHKERRQ(ierr); }
  }

  /* A1 */
  if (Istart==0) { ierr = MatSetValue(A[1],0,0,-deltasq/4,INSERT_VALUES);CHKERRQ(ierr); }
  if (Iend==n+1) { ierr = MatSetValue(A[1],n,n,deltasq/4,INSERT_VALUES);CHKERRQ(ierr); }

  /* A2 */
  alpha = 1.0/h;
  for (i=Istart;i<Iend;i++) {
    v = 2.0;
    if (i==0 || i==n) v = 1.0;
    ierr = MatSetValue(A[2],i,i,v*alpha,ADD_VALUES);CHKERRQ(ierr);
    if (i>0) { ierr = MatSetValue(A[2],i,i-1,-alpha,ADD_VALUES);CHKERRQ(ierr); }
    if (i<n) { ierr = MatSetValue(A[2],i,i+1,-alpha,ADD_VALUES);CHKERRQ(ierr); }
  }
  ierr = PetscMalloc3(n+1,&md,n+1,&supd,n+1,&subd);CHKERRQ(ierr);

  md[0]   = 2.0*q[1];
  supd[1] = q[1];
  subd[0] = q[1];

  for (k=1;k<=NL-2;k++) {

    end_ct = k*(nlocal+1);
    start_ct = end_ct-nlocal;

    for (j=start_ct;j<end_ct;j++) {
      md[j] = 4*q[k];
      supd[j+1] = q[k];
      subd[j] = q[k];
    }

    if (k < 4) {  /* interface points */
      md[end_ct] = 4*(q[k] + q[k+1])/2.0;
      supd[end_ct+1] = q[k+1];
      subd[end_ct] = q[k+1];
    }

  }

  md[n] = 2*q[NL-2];
  supd[n] = q[NL-2];
  subd[n] = q[NL-2];

  alpha = -h/6.0;
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(A[2],i,i,md[i]*alpha,ADD_VALUES);CHKERRQ(ierr);
    if (i>0) { ierr = MatSetValue(A[2],i,i-1,subd[i-1]*alpha,ADD_VALUES);CHKERRQ(ierr); }
    if (i<n) { ierr = MatSetValue(A[2],i,i+1,supd[i+1]*alpha,ADD_VALUES);CHKERRQ(ierr); }
  }
  ierr = PetscFree3(md,supd,subd);CHKERRQ(ierr);

  /* A3 */
  if (Istart==0) { ierr = MatSetValue(A[3],0,0,1.0,INSERT_VALUES);CHKERRQ(ierr); }
  if (Iend==n+1) { ierr = MatSetValue(A[3],n,n,1.0,INSERT_VALUES);CHKERRQ(ierr); }

  /* A4 */
  alpha = (h/6);
  for (i=Istart;i<Iend;i++) {
    v = 4.0;
    if (i==0 || i==n) v = 2.0;
    ierr = MatSetValue(A[4],i,i,v*alpha,INSERT_VALUES);CHKERRQ(ierr);
    if (i>0) { ierr = MatSetValue(A[4],i,i-1,alpha,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n) { ierr = MatSetValue(A[4],i,i+1,alpha,INSERT_VALUES);CHKERRQ(ierr); }
  }

  /* assemble matrices */
  for (i=0;i<NMAT;i++) {
    ierr = MatAssemblyBegin(A[i],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  for (i=0;i<NMAT;i++) {
    ierr = MatAssemblyEnd(A[i],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  ierr = PEPSetOperators(pep,NMAT,A);CHKERRQ(ierr);
  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);
  ierr = PEPSolve(pep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = PEPReasonView(pep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  for (i=0;i<NMAT;i++) {
    ierr = MatDestroy(&A[i]);CHKERRQ(ierr);
  }
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -pep_type {{toar linear}} -pep_nev 4 -st_type sinvert -terse
      requires: !single

TEST*/
