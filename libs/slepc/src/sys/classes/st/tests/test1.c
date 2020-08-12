/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test ST with shell matrices.\n\n";

#include <slepcst.h>

static PetscErrorCode MatGetDiagonal_Shell(Mat S,Vec diag);
static PetscErrorCode MatMultTranspose_Shell(Mat S,Vec x,Vec y);
static PetscErrorCode MatMult_Shell(Mat S,Vec x,Vec y);
static PetscErrorCode MatDuplicate_Shell(Mat S,MatDuplicateOption op,Mat *M);

static PetscErrorCode MyShellMatCreate(Mat *A,Mat *M)
{
  PetscErrorCode ierr;
  MPI_Comm       comm;
  PetscInt       n;

  PetscFunctionBeginUser;
  ierr = MatGetSize(*A,&n,NULL);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)*A,&comm);CHKERRQ(ierr);
  ierr = MatCreateShell(comm,PETSC_DECIDE,PETSC_DECIDE,n,n,A,M);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*M,MATOP_MULT,(void(*)(void))MatMult_Shell);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*M,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMultTranspose_Shell);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*M,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Shell);CHKERRQ(ierr);
  ierr = MatShellSetOperation(*M,MATOP_DUPLICATE,(void(*)(void))MatDuplicate_Shell);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  Mat            A,S,mat[1];
  ST             st;
  Vec            v,w;
  STType         type;
  KSP            ksp;
  PC             pc;
  PetscScalar    sigma;
  PetscInt       n=10,i,Istart,Iend;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian with shell matrices, n=%D\n\n",n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix for the 1-D Laplacian
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) { ierr = MatSetValue(A,i,i-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A,i,i+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,i,i,2.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* create the shell version of A */
  ierr = MyShellMatCreate(&A,&S);CHKERRQ(ierr);

  /* work vectors */
  ierr = MatCreateVecs(A,&v,&w);CHKERRQ(ierr);
  ierr = VecSet(v,1.0);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the spectral transformation object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = STCreate(PETSC_COMM_WORLD,&st);CHKERRQ(ierr);
  mat[0] = S;
  ierr = STSetMatrices(st,1,mat);CHKERRQ(ierr);
  ierr = STSetTransform(st,PETSC_TRUE);CHKERRQ(ierr);
  ierr = STSetFromOptions(st);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              Apply the transformed operator for several ST's
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* shift, sigma=0.0 */
  ierr = STSetUp(st);CHKERRQ(ierr);
  ierr = STGetType(st,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"ST type %s\n",type);CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);
  ierr = STApplyTranspose(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* shift, sigma=0.1 */
  sigma = 0.1;
  ierr = STSetShift(st,sigma);CHKERRQ(ierr);
  ierr = STGetShift(st,&sigma);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma));CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* sinvert, sigma=0.1 */
  ierr = STPostSolve(st);CHKERRQ(ierr);   /* undo changes if inplace */
  ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = STGetType(st,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"ST type %s\n",type);CHKERRQ(ierr);
  ierr = STGetShift(st,&sigma);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma));CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* sinvert, sigma=-0.5 */
  sigma = -0.5;
  ierr = STSetShift(st,sigma);CHKERRQ(ierr);
  ierr = STGetShift(st,&sigma);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma));CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  ierr = STDestroy(&st);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&S);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

static PetscErrorCode MatMult_Shell(Mat S,Vec x,Vec y)
{
  PetscErrorCode    ierr;
  Mat               *A;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(S,&A);CHKERRQ(ierr);
  ierr = MatMult(*A,x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMultTranspose_Shell(Mat S,Vec x,Vec y)
{
  PetscErrorCode    ierr;
  Mat               *A;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(S,&A);CHKERRQ(ierr);
  ierr = MatMultTranspose(*A,x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatGetDiagonal_Shell(Mat S,Vec diag)
{
  PetscErrorCode    ierr;
  Mat               *A;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(S,&A);CHKERRQ(ierr);
  ierr = MatGetDiagonal(*A,diag);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatDuplicate_Shell(Mat S,MatDuplicateOption op,Mat *M)
{
  PetscErrorCode ierr;
  Mat            *A;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(S,&A);CHKERRQ(ierr);
  ierr = MyShellMatCreate(A,M);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      args: -st_matmode {{inplace shell}}
      output_file: output/test1_1.out
      requires: !single

TEST*/
