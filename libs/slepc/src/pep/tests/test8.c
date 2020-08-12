/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test interface functions of polynomial JD.\n\n"
  "This is based on ex16.c. The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n"
  "  -m <m>, where <m> = number of grid subdivisions in y dimension.\n\n";

#include <slepcpep.h>

int main(int argc,char **argv)
{
  Mat             M,C,K,A[3];      /* problem matrices */
  PEP             pep;             /* polynomial eigenproblem solver context */
  PetscInt        N,n=10,m,Istart,Iend,II,i,j,midx;
  PetscReal       restart,fix;
  PetscBool       flag,reuse;
  PEPJDProjection proj;
  PetscErrorCode  ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,&flag);CHKERRQ(ierr);
  if (!flag) m=n;
  N = n*m;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nQuadratic Eigenproblem, N=%D (%Dx%D grid)\n\n",N,n,m);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is the 2-D Laplacian */
  ierr = MatCreate(PETSC_COMM_WORLD,&K);CHKERRQ(ierr);
  ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(K);CHKERRQ(ierr);
  ierr = MatSetUp(K);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(K,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) { ierr = MatSetValue(K,II,II-n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<m-1) { ierr = MatSetValue(K,II,II+n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j>0) { ierr = MatSetValue(K,II,II-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(K,II,II+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(K,II,II,4.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* C is the 1-D Laplacian on horizontal lines */
  ierr = MatCreate(PETSC_COMM_WORLD,&C);CHKERRQ(ierr);
  ierr = MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(C);CHKERRQ(ierr);
  ierr = MatSetUp(C);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(C,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (j>0) { ierr = MatSetValue(C,II,II-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(C,II,II+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(C,II,II,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* M is a diagonal matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&M);CHKERRQ(ierr);
  ierr = MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(M);CHKERRQ(ierr);
  ierr = MatSetUp(M);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(M,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    ierr = MatSetValue(M,II,II,(PetscReal)(II+1),INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  A[0] = K; A[1] = C; A[2] = M;
  ierr = PEPSetOperators(pep,3,A);CHKERRQ(ierr);
  ierr = PEPSetType(pep,PEPJD);CHKERRQ(ierr);

  /*
     Test interface functions of STOAR solver
  */
  ierr = PEPJDGetRestart(pep,&restart);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Restart parameter before changing = %g",(double)restart);CHKERRQ(ierr);
  ierr = PEPJDSetRestart(pep,0.3);CHKERRQ(ierr);
  ierr = PEPJDGetRestart(pep,&restart);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %g\n",(double)restart);CHKERRQ(ierr);

  ierr = PEPJDGetFix(pep,&fix);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Fix parameter before changing = %g",(double)fix);CHKERRQ(ierr);
  ierr = PEPJDSetFix(pep,0.001);CHKERRQ(ierr);
  ierr = PEPJDGetFix(pep,&fix);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %g\n",(double)fix);CHKERRQ(ierr);

  ierr = PEPJDGetReusePreconditioner(pep,&reuse);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Reuse preconditioner flag before changing = %d",(int)reuse);CHKERRQ(ierr);
  ierr = PEPJDSetReusePreconditioner(pep,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PEPJDGetReusePreconditioner(pep,&reuse);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)reuse);CHKERRQ(ierr);

  ierr = PEPJDGetProjection(pep,&proj);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Projection type before changing = %d",(int)proj);CHKERRQ(ierr);
  ierr = PEPJDSetProjection(pep,PEP_JD_PROJECTION_ORTHOGONAL);CHKERRQ(ierr);
  ierr = PEPJDGetProjection(pep,&proj);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %d\n",(int)proj);CHKERRQ(ierr);

  ierr = PEPJDGetMinimalityIndex(pep,&midx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Minimality index before changing = %D",midx);CHKERRQ(ierr);
  ierr = PEPJDSetMinimalityIndex(pep,2);CHKERRQ(ierr);
  ierr = PEPJDGetMinimalityIndex(pep,&midx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," ... changed to %D\n",midx);CHKERRQ(ierr);

  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPSolve(pep);CHKERRQ(ierr);
  ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL);CHKERRQ(ierr);
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      args: -n 12 -pep_nev 2 -pep_ncv 21 -pep_conv_abs

TEST*/
