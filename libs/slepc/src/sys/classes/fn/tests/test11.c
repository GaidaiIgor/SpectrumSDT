/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Define the function

        f(x) = (exp(x)-1)/x    (the phi_1 function)

   with the following tree:

            f(x)                  f(x)              (combined by division)
           /    \                 p(x) = x          (polynomial)
        a(x)    p(x)              a(x)              (combined by addition)
       /    \                     e(x) = exp(x)     (exponential)
     e(x)   c(x)                  c(x) = -1         (constant)
*/

static char help[] = "Another test of a combined function.\n\n";

#include <slepcfn.h>

/*
   Compute matrix function B = A\(exp(A)-I)
 */
PetscErrorCode TestMatCombine(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace)
{
  PetscErrorCode ierr;
  PetscBool      set,flg;
  PetscInt       n;
  Mat            F;
  Vec            v,f0;
  PetscReal      nrm;

  PetscFunctionBeginUser;
  ierr = MatGetSize(A,&n,NULL);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&F);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)F,"F");CHKERRQ(ierr);
  /* compute matrix function */
  if (inplace) {
    ierr = MatCopy(A,F,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = MatIsHermitianKnown(A,&set,&flg);CHKERRQ(ierr);
    if (set && flg) { ierr = MatSetOption(F,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr); }
    ierr = FNEvaluateFunctionMat(fn,F,NULL);CHKERRQ(ierr);
  } else {
    ierr = FNEvaluateFunctionMat(fn,A,F);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(A,viewer);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed f(A) - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(F,viewer);CHKERRQ(ierr);
  }
  /* print matrix norm for checking */
  ierr = MatNorm(F,NORM_1,&nrm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The 1-norm of f(A) is %6.3f\n",(double)nrm);CHKERRQ(ierr);
  /* check FNEvaluateFunctionMatVec() */
  ierr = MatCreateVecs(A,&v,&f0);CHKERRQ(ierr);
  ierr = MatGetColumnVector(F,f0,0);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMatVec(fn,A,v);CHKERRQ(ierr);
  ierr = VecAXPY(v,-1.0,f0);CHKERRQ(ierr);
  ierr = VecNorm(v,NORM_2,&nrm);CHKERRQ(ierr);
  if (nrm>100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-v is %g\n",(double)nrm);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&F);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&f0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  FN             f,p,a,e,c,f1,f2;
  FNCombineType  ctype;
  Mat            A;
  PetscInt       i,j,n=10,np;
  PetscScalar    x,y,yp,*As,coeffs[10];
  char           strx[50],str[50];
  PetscViewer    viewer;
  PetscBool      verbose,inplace;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-inplace",&inplace);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Phi1 via a combined function, n=%D.\n",n);CHKERRQ(ierr);

  /* Create function */

  /* e(x) = exp(x) */
  ierr = FNCreate(PETSC_COMM_WORLD,&e);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)e,"e");CHKERRQ(ierr);
  ierr = FNSetType(e,FNEXP);CHKERRQ(ierr);
  ierr = FNSetFromOptions(e);CHKERRQ(ierr);
  /* c(x) = -1 */
  ierr = FNCreate(PETSC_COMM_WORLD,&c);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)c,"c");CHKERRQ(ierr);
  ierr = FNSetType(c,FNRATIONAL);CHKERRQ(ierr);
  ierr = FNSetFromOptions(c);CHKERRQ(ierr);
  np = 1;
  coeffs[0] = -1.0;
  ierr = FNRationalSetNumerator(c,np,coeffs);CHKERRQ(ierr);
  /* a(x) */
  ierr = FNCreate(PETSC_COMM_WORLD,&a);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)a,"a");CHKERRQ(ierr);
  ierr = FNSetType(a,FNCOMBINE);CHKERRQ(ierr);
  ierr = FNSetFromOptions(a);CHKERRQ(ierr);
  ierr = FNCombineSetChildren(a,FN_COMBINE_ADD,e,c);CHKERRQ(ierr);
  /* p(x) = x */
  ierr = FNCreate(PETSC_COMM_WORLD,&p);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)p,"p");CHKERRQ(ierr);
  ierr = FNSetType(p,FNRATIONAL);CHKERRQ(ierr);
  ierr = FNSetFromOptions(p);CHKERRQ(ierr);
  np = 2;
  coeffs[0] = 1.0; coeffs[1] = 0.0;
  ierr = FNRationalSetNumerator(p,np,coeffs);CHKERRQ(ierr);
  /* f(x) */
  ierr = FNCreate(PETSC_COMM_WORLD,&f);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)f,"f");CHKERRQ(ierr);
  ierr = FNSetType(f,FNCOMBINE);CHKERRQ(ierr);
  ierr = FNSetFromOptions(f);CHKERRQ(ierr);
  ierr = FNCombineSetChildren(f,FN_COMBINE_DIVIDE,a,p);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = FNCombineGetChildren(f,&ctype,&f1,&f2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Two functions combined with division:\n");CHKERRQ(ierr);
  ierr = FNView(f1,viewer);CHKERRQ(ierr);
  ierr = FNView(f2,viewer);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Scalar evaluation */
  x = 2.2;
  ierr = SlepcSNPrintfScalar(strx,50,x,PETSC_FALSE);CHKERRQ(ierr);
  ierr = FNEvaluateFunction(f,x,&y);CHKERRQ(ierr);
  ierr = FNEvaluateDerivative(f,x,&yp);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,y,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f(%s)=%s\n",strx,str);CHKERRQ(ierr);
  ierr = SlepcSNPrintfScalar(str,50,yp,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  f'(%s)=%s\n",strx,str);CHKERRQ(ierr);

  /* Create matrices */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A,"A");CHKERRQ(ierr);

  /* Fill A with 1-D Laplacian matrix */
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  for (i=0;i<n;i++) As[i+i*n]=2.0;
  j=1;
  for (i=0;i<n-j;i++) { As[i+(i+j)*n]=-1.0; As[(i+j)+i*n]=-1.0; }
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = TestMatCombine(f,A,viewer,verbose,inplace);CHKERRQ(ierr);

  /* Repeat with same matrix as non-symmetric */
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE);CHKERRQ(ierr);
  ierr = TestMatCombine(f,A,viewer,verbose,inplace);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = FNDestroy(&f);CHKERRQ(ierr);
  ierr = FNDestroy(&p);CHKERRQ(ierr);
  ierr = FNDestroy(&a);CHKERRQ(ierr);
  ierr = FNDestroy(&e);CHKERRQ(ierr);
  ierr = FNDestroy(&c);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1

   test:
      suffix: 2
      nsize: 1
      args: -inplace
      output_file: output/test11_1.out

TEST*/
