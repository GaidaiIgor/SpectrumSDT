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

        f(x) = (1-x^2) exp(-x/(1+x^2))

   with the following tree:

            f(x)                  f(x)              (combined by product)
           /    \                 g(x) = 1-x^2      (polynomial)
        g(x)    h(x)              h(x)              (combined by composition)
               /    \             r(x) = -x/(1+x^2) (rational)
             r(x)   e(x)          e(x) = exp(x)     (exponential)
*/

static char help[] = "Test combined function.\n\n";

#include <slepcfn.h>

/*
   Compute matrix function B = (I-A^2) exp(-(I+A^2)\A)
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
  FN             f,g,h,e,r,fcopy;
  Mat            A;
  PetscInt       i,j,n=10,np,nq;
  PetscScalar    x,y,yp,*As,p[10],q[10];
  char           strx[50],str[50];
  PetscViewer    viewer;
  PetscBool      verbose,inplace;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-inplace",&inplace);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Combined function, n=%D.\n",n);CHKERRQ(ierr);

  /* Create function */

  /* e(x) = exp(x) */
  ierr = FNCreate(PETSC_COMM_WORLD,&e);CHKERRQ(ierr);
  ierr = FNSetType(e,FNEXP);CHKERRQ(ierr);
  ierr = FNSetFromOptions(e);CHKERRQ(ierr);
  /* r(x) = x/(1+x^2) */
  ierr = FNCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
  ierr = FNSetType(r,FNRATIONAL);CHKERRQ(ierr);
  ierr = FNSetFromOptions(r);CHKERRQ(ierr);
  np = 2; nq = 3;
  p[0] = -1.0; p[1] = 0.0;
  q[0] = 1.0; q[1] = 0.0; q[2] = 1.0;
  ierr = FNRationalSetNumerator(r,np,p);CHKERRQ(ierr);
  ierr = FNRationalSetDenominator(r,nq,q);CHKERRQ(ierr);
  /* h(x) */
  ierr = FNCreate(PETSC_COMM_WORLD,&h);CHKERRQ(ierr);
  ierr = FNSetType(h,FNCOMBINE);CHKERRQ(ierr);
  ierr = FNSetFromOptions(h);CHKERRQ(ierr);
  ierr = FNCombineSetChildren(h,FN_COMBINE_COMPOSE,r,e);CHKERRQ(ierr);
  /* g(x) = 1-x^2 */
  ierr = FNCreate(PETSC_COMM_WORLD,&g);CHKERRQ(ierr);
  ierr = FNSetType(g,FNRATIONAL);CHKERRQ(ierr);
  ierr = FNSetFromOptions(g);CHKERRQ(ierr);
  np = 3;
  p[0] = -1.0; p[1] = 0.0; p[2] = 1.0;
  ierr = FNRationalSetNumerator(g,np,p);CHKERRQ(ierr);
  /* f(x) */
  ierr = FNCreate(PETSC_COMM_WORLD,&f);CHKERRQ(ierr);
  ierr = FNSetType(f,FNCOMBINE);CHKERRQ(ierr);
  ierr = FNSetFromOptions(f);CHKERRQ(ierr);
  ierr = FNCombineSetChildren(f,FN_COMBINE_MULTIPLY,g,h);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = FNView(f,viewer);CHKERRQ(ierr);
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

  /* Test duplication */
  ierr = FNDuplicate(f,PetscObjectComm((PetscObject)f),&fcopy);CHKERRQ(ierr);

  /* Create matrices */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A,"A");CHKERRQ(ierr);

  /* Fill A with a symmetric Toeplitz matrix */
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  for (i=0;i<n;i++) As[i+i*n]=2.0;
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) { As[i+(i+j)*n]=1.0; As[(i+j)+i*n]=1.0; }
  }
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = TestMatCombine(fcopy,A,viewer,verbose,inplace);CHKERRQ(ierr);

  /* Repeat with same matrix as non-symmetric */
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE);CHKERRQ(ierr);
  ierr = TestMatCombine(fcopy,A,viewer,verbose,inplace);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = FNDestroy(&f);CHKERRQ(ierr);
  ierr = FNDestroy(&fcopy);CHKERRQ(ierr);
  ierr = FNDestroy(&g);CHKERRQ(ierr);
  ierr = FNDestroy(&h);CHKERRQ(ierr);
  ierr = FNDestroy(&e);CHKERRQ(ierr);
  ierr = FNDestroy(&r);CHKERRQ(ierr);
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
      output_file: output/test6_1.out

TEST*/
