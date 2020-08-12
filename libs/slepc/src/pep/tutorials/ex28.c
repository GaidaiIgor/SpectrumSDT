/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "A quadratic eigenproblem defined using shell matrices.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x and y dimensions.\n\n";

#include <slepcpep.h>

/*
   User-defined routines
*/
PetscErrorCode MatMult_Laplacian2D(Mat A,Vec x,Vec y);
PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A,Vec diag);
PetscErrorCode MatMult_Zero(Mat A,Vec x,Vec y);
PetscErrorCode MatGetDiagonal_Zero(Mat A,Vec diag);
PetscErrorCode MatMult_Identity(Mat A,Vec x,Vec y);
PetscErrorCode MatGetDiagonal_Identity(Mat A,Vec diag);

int main(int argc,char **argv)
{
  Mat            M,C,K,A[3];      /* problem matrices */
  PEP            pep;             /* polynomial eigenproblem solver context */
  PEPType        type;
  PetscInt       N,n=10,nev;
  PetscMPIInt    size;
  PetscBool      terse;
  PetscErrorCode ierr;
  ST             st;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only");

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  N = n*n;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nQuadratic Eigenproblem with shell matrices, N=%D (%Dx%D grid)\n\n",N,n,n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is the 2-D Laplacian */
  ierr = MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,&n,&K);CHKERRQ(ierr);
  ierr = MatShellSetOperation(K,MATOP_MULT,(void(*)(void))MatMult_Laplacian2D);CHKERRQ(ierr);
  ierr = MatShellSetOperation(K,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_Laplacian2D);CHKERRQ(ierr);
  ierr = MatShellSetOperation(K,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Laplacian2D);CHKERRQ(ierr);

  /* C is the zero matrix */
  ierr = MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,NULL,&C);CHKERRQ(ierr);
  ierr = MatShellSetOperation(C,MATOP_MULT,(void(*)(void))MatMult_Zero);CHKERRQ(ierr);
  ierr = MatShellSetOperation(C,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_Zero);CHKERRQ(ierr);
  ierr = MatShellSetOperation(C,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Zero);CHKERRQ(ierr);

  /* M is the identity matrix */
  ierr = MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,NULL,&M);CHKERRQ(ierr);
  ierr = MatShellSetOperation(M,MATOP_MULT,(void(*)(void))MatMult_Identity);CHKERRQ(ierr);
  ierr = MatShellSetOperation(M,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_Identity);CHKERRQ(ierr);
  ierr = MatShellSetOperation(M,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Identity);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);

  /*
     Set matrices and problem type
  */
  A[0] = K; A[1] = C; A[2] = M;
  ierr = PEPSetOperators(pep,3,A);CHKERRQ(ierr);
  ierr = PEPGetST(pep,&st);CHKERRQ(ierr);
  ierr = STSetMatMode(st,ST_MATMODE_SHELL);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PEPSolve(pep);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = PEPGetType(pep,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = PEPGetDimensions(pep,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = PEPErrorView(pep,PEP_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = PEPReasonView(pep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PEPErrorView(pep,PEP_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = MatDestroy(&K);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*
    Compute the matrix vector multiplication y<---T*x where T is a nx by nx
    tridiagonal matrix with DD on the diagonal, DL on the subdiagonal, and
    DU on the superdiagonal.
 */
static void tv(int nx,const PetscScalar *x,PetscScalar *y)
{
  PetscScalar dd,dl,du;
  int         j;

  dd  = 4.0;
  dl  = -1.0;
  du  = -1.0;

  y[0] =  dd*x[0] + du*x[1];
  for (j=1;j<nx-1;j++)
    y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
  y[nx-1] = dl*x[nx-2] + dd*x[nx-1];
}

/*
    Matrix-vector product subroutine for the 2D Laplacian.

    The matrix used is the 2 dimensional discrete Laplacian on unit square with
    zero Dirichlet boundary condition.

    Computes y <-- A*x, where A is the block tridiagonal matrix

                 | T -I          |
                 |-I  T -I       |
             A = |   -I  T       |
                 |        ...  -I|
                 |           -I T|

    The subroutine TV is called to compute y<--T*x.
 */
PetscErrorCode MatMult_Laplacian2D(Mat A,Vec x,Vec y)
{
  void              *ctx;
  int               nx,lo,i,j;
  const PetscScalar *px;
  PetscScalar       *py;
  PetscErrorCode    ierr;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,&ctx);CHKERRQ(ierr);
  nx = *(int*)ctx;
  ierr = VecGetArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecGetArray(y,&py);CHKERRQ(ierr);

  tv(nx,&px[0],&py[0]);
  for (i=0;i<nx;i++) py[i] -= px[nx+i];

  for (j=2;j<nx;j++) {
    lo = (j-1)*nx;
    tv(nx,&px[lo],&py[lo]);
    for (i=0;i<nx;i++) py[lo+i] -= px[lo-nx+i] + px[lo+nx+i];
  }

  lo = (nx-1)*nx;
  tv(nx,&px[lo],&py[lo]);
  for (i=0;i<nx;i++) py[lo+i] -= px[lo-nx+i];

  ierr = VecRestoreArrayRead(x,&px);CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A,Vec diag)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecSet(diag,4.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    Matrix-vector product subroutine for the Null matrix.
 */
PetscErrorCode MatMult_Zero(Mat A,Vec x,Vec y)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecSet(y,0.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatGetDiagonal_Zero(Mat A,Vec diag)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecSet(diag,0.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    Matrix-vector product subroutine for the Identity matrix.
 */
PetscErrorCode MatMult_Identity(Mat A,Vec x,Vec y)
{
  PetscErrorCode    ierr;

  PetscFunctionBeginUser;
  ierr = VecCopy(x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatGetDiagonal_Identity(Mat A,Vec diag)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecSet(diag,1.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      args: -pep_type {{toar qarnoldi linear}} -pep_nev 4 -terse
      filter: grep -v Solution
      requires: !single
      output_file: output/ex28_1.out

TEST*/
