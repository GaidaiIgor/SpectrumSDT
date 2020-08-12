/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test ST with one matrix.\n\n";

#include <slepcst.h>

int main(int argc,char **argv)
{
  Mat            A,B,mat[1];
  ST             st;
  Vec            v,w;
  STType         type;
  PetscScalar    sigma,tau;
  PetscInt       n=10,i,Istart,Iend;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian, n=%D\n\n",n);CHKERRQ(ierr);

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
  ierr = MatCreateVecs(A,&v,&w);CHKERRQ(ierr);
  ierr = VecSet(v,1.0);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the spectral transformation object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = STCreate(PETSC_COMM_WORLD,&st);CHKERRQ(ierr);
  mat[0] = A;
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

  /* shift, sigma=0.1 */
  sigma = 0.1;
  ierr = STSetShift(st,sigma);CHKERRQ(ierr);
  ierr = STGetShift(st,&sigma);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma));CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);
  ierr = STApplyTranspose(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* sinvert, sigma=0.1 */
  ierr = STPostSolve(st);CHKERRQ(ierr);   /* undo changes if inplace */
  ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
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

  /* cayley, sigma=-0.5, tau=-0.5 (equal to sigma by default) */
  ierr = STPostSolve(st);CHKERRQ(ierr);   /* undo changes if inplace */
  ierr = STSetType(st,STCAYLEY);CHKERRQ(ierr);
  ierr = STSetUp(st);CHKERRQ(ierr);
  ierr = STGetType(st,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"ST type %s\n",type);CHKERRQ(ierr);
  ierr = STGetShift(st,&sigma);CHKERRQ(ierr);
  ierr = STCayleyGetAntishift(st,&tau);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With shift=%g, antishift=%g\n",(double)PetscRealPart(sigma),(double)PetscRealPart(tau));CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* cayley, sigma=1.1, tau=1.1 (still equal to sigma) */
  sigma = 1.1;
  ierr = STSetShift(st,sigma);CHKERRQ(ierr);
  ierr = STGetShift(st,&sigma);CHKERRQ(ierr);
  ierr = STCayleyGetAntishift(st,&tau);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With shift=%g, antishift=%g\n",(double)PetscRealPart(sigma),(double)PetscRealPart(tau));CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* cayley, sigma=1.1, tau=-1.0 */
  tau = -1.0;
  ierr = STCayleySetAntishift(st,tau);CHKERRQ(ierr);
  ierr = STGetShift(st,&sigma);CHKERRQ(ierr);
  ierr = STCayleyGetAntishift(st,&tau);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With shift=%g, antishift=%g\n",(double)PetscRealPart(sigma),(double)PetscRealPart(tau));CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* check inner product matrix in Cayley */
  ierr = STGetBilinearForm(st,&B);CHKERRQ(ierr);
  ierr = MatMult(B,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  /* check again sinvert, sigma=0.1 */
  ierr = STPostSolve(st);CHKERRQ(ierr);   /* undo changes if inplace */
  ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
  ierr = STGetType(st,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"ST type %s\n",type);CHKERRQ(ierr);
  sigma = 0.1;
  ierr = STSetShift(st,sigma);CHKERRQ(ierr);
  ierr = STGetShift(st,&sigma);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"With shift=%g\n",(double)PetscRealPart(sigma));CHKERRQ(ierr);
  ierr = STApply(st,v,w);CHKERRQ(ierr);
  ierr = VecView(w,NULL);CHKERRQ(ierr);

  ierr = STDestroy(&st);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -st_matmode {{copy inplace shell}}
      output_file: output/test2_1.out
      requires: !single

TEST*/
