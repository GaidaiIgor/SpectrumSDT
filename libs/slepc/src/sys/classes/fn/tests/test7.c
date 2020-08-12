/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test matrix square root.\n\n";

#include <slepcfn.h>

/*
   Compute matrix square root B = sqrtm(A)
   Check result as norm(B*B-A)
 */
PetscErrorCode TestMatSqrt(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace)
{
  PetscErrorCode ierr;
  PetscScalar    tau,eta;
  PetscReal      nrm;
  PetscBool      set,flg;
  PetscInt       n;
  Mat            S,R;
  Vec            v,f0;

  PetscFunctionBeginUser;
  ierr = MatGetSize(A,&n,NULL);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&S);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)S,"S");CHKERRQ(ierr);
  ierr = FNGetScale(fn,&tau,&eta);CHKERRQ(ierr);
  /* compute square root */
  if (inplace) {
    ierr = MatCopy(A,S,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = MatIsHermitianKnown(A,&set,&flg);CHKERRQ(ierr);
    if (set && flg) { ierr = MatSetOption(S,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr); }
    ierr = FNEvaluateFunctionMat(fn,S,NULL);CHKERRQ(ierr);
  } else {
    ierr = FNEvaluateFunctionMat(fn,A,S);CHKERRQ(ierr);
  }
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix A - - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(A,viewer);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed sqrtm(A) - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(S,viewer);CHKERRQ(ierr);
  }
  /* check error ||S*S-A||_F */
  ierr = MatMatMult(S,S,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&R);CHKERRQ(ierr);
  if (eta!=1.0) {
    ierr = MatScale(R,1.0/(eta*eta));CHKERRQ(ierr);
  }
  ierr = MatAXPY(R,-tau,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatNorm(R,NORM_FROBENIUS,&nrm);CHKERRQ(ierr);
  if (nrm<100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||S*S-A||_F < 100*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||S*S-A||_F = %g\n",(double)nrm);CHKERRQ(ierr);
  }
  /* check FNEvaluateFunctionMatVec() */
  ierr = MatCreateVecs(A,&v,&f0);CHKERRQ(ierr);
  ierr = MatGetColumnVector(S,f0,0);CHKERRQ(ierr);
  ierr = FNEvaluateFunctionMatVec(fn,A,v);CHKERRQ(ierr);
  ierr = VecAXPY(v,-1.0,f0);CHKERRQ(ierr);
  ierr = VecNorm(v,NORM_2,&nrm);CHKERRQ(ierr);
  if (nrm>100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: the norm of f(A)*e_1-v is %g\n",(double)nrm);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&S);CHKERRQ(ierr);
  ierr = MatDestroy(&R);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&f0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  FN             fn;
  Mat            A;
  PetscInt       i,j,n=10;
  PetscScalar    *As;
  PetscViewer    viewer;
  PetscBool      verbose,inplace;
  PetscRandom    myrand;
  PetscReal      v;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-inplace",&inplace);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix square root, n=%D.\n",n);CHKERRQ(ierr);

  /* Create function object */
  ierr = FNCreate(PETSC_COMM_WORLD,&fn);CHKERRQ(ierr);
  ierr = FNSetType(fn,FNSQRT);CHKERRQ(ierr);
  ierr = FNSetFromOptions(fn);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = FNView(fn,viewer);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Create matrix */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A,"A");CHKERRQ(ierr);

  /* Compute square root of a symmetric matrix A */
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  for (i=0;i<n;i++) As[i+i*n]=2.5;
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) { As[i+(i+j)*n]=1.0; As[(i+j)+i*n]=1.0; }
  }
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = TestMatSqrt(fn,A,viewer,verbose,inplace);CHKERRQ(ierr);

  /* Repeat with upper triangular A */
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) As[(i+j)+i*n]=0.0;
  }
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE);CHKERRQ(ierr);
  ierr = TestMatSqrt(fn,A,viewer,verbose,inplace);CHKERRQ(ierr);

  /* Repeat with non-symmetic A */
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&myrand);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(myrand);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(myrand,0.0,1.0);CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
  for (j=1;j<3;j++) {
    for (i=0;i<n-j;i++) {
      ierr = PetscRandomGetValueReal(myrand,&v);CHKERRQ(ierr);
      As[(i+j)+i*n]=v;
    }
  }
  ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  ierr = PetscRandomDestroy(&myrand);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_FALSE);CHKERRQ(ierr);
  ierr = TestMatSqrt(fn,A,viewer,verbose,inplace);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = FNDestroy(&fn);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -fn_scale .05,2 -n 100 -fn_method {{0 1 2}shared output}
      filter: grep -v "computing matrix functions"
      output_file: output/test7_1.out
      timeoutfactor: 2

   test:
      suffix: 1_sadeghi
      nsize: 1
      args: -fn_scale .05,2 -n 100 -fn_method 3
      requires: !single
      filter: grep -v "computing matrix functions"
      output_file: output/test7_1.out

   test:
      suffix: 2
      nsize: 1
      args: -fn_scale .05,2 -n 100 -inplace -fn_method {{0 1 2}shared output}
      filter: grep -v "computing matrix functions"
      output_file: output/test7_1.out
      timeoutfactor: 2

   test:
      suffix: 2_sadeghi
      nsize: 1
      args: -fn_scale .05,2 -n 100 -inplace -fn_method 3
      requires: !single
      filter: grep -v "computing matrix functions"
      output_file: output/test7_1.out

   test:
      suffix: 3
      nsize: 3
      args: -fn_scale .05,2 -n 100 -fn_parallel synchronized
      filter: grep -v "computing matrix functions" | grep -v "SYNCHRONIZED" | sed -e "s/3 MPI/1 MPI/g"
      output_file: output/test7_1.out

   test:
      suffix: 4
      nsize: 3
      args: -fn_scale .05,2 -n 100 -inplace -fn_parallel synchronized
      filter: grep -v "computing matrix functions" | grep -v "SYNCHRONIZED" | sed -e "s/3 MPI/1 MPI/g"
      output_file: output/test7_1.out

TEST*/
