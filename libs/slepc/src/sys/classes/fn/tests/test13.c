/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test matrix logarithm.\n\n";

#include <slepcfn.h>

/*
   Compute matrix logarithm B = logm(A)
 */
PetscErrorCode TestMatLog(FN fn,Mat A,PetscViewer viewer,PetscBool verbose,PetscBool inplace)
{
  PetscErrorCode ierr;
  PetscBool      set,flg;
  PetscScalar    tau,eta;
  PetscInt       n;
  Mat            F,R;
  Vec            v,f0;
  FN             fnexp;
  PetscReal      nrm;

  PetscFunctionBeginUser;
  ierr = MatGetSize(A,&n,NULL);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&F);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)F,"F");CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&R);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)R,"R");CHKERRQ(ierr);
  ierr = FNGetScale(fn,&tau,&eta);CHKERRQ(ierr);
  /* compute matrix logarithm */
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
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed logm(A) - - - - - - -\n");CHKERRQ(ierr);
    ierr = MatView(F,viewer);CHKERRQ(ierr);
  }
  /* check error ||expm(F)-A||_F */
  ierr = FNCreate(PETSC_COMM_WORLD,&fnexp);CHKERRQ(ierr);
  ierr = FNSetType(fnexp,FNEXP);CHKERRQ(ierr);
  ierr = MatCopy(F,R,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  if (eta!=1.0) {
    ierr = MatScale(R,1.0/eta);CHKERRQ(ierr);
  }
  ierr = FNEvaluateFunctionMat(fnexp,R,NULL);CHKERRQ(ierr);
  ierr = FNDestroy(&fnexp);CHKERRQ(ierr);
  ierr = MatAXPY(R,-tau,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatNorm(R,NORM_FROBENIUS,&nrm);CHKERRQ(ierr);
  if (nrm<100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||expm(F)-A||_F < 100*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"||expm(F)-A||_F = %g\n",(double)nrm);CHKERRQ(ierr);
  }
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
  PetscBool      verbose,inplace,random,triang;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-inplace",&inplace);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-random",&random);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-triang",&triang);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix logarithm, n=%D.\n",n);CHKERRQ(ierr);

  /* Create logarithm function object */
  ierr = FNCreate(PETSC_COMM_WORLD,&fn);CHKERRQ(ierr);
  ierr = FNSetType(fn,FNLOG);CHKERRQ(ierr);
  ierr = FNSetFromOptions(fn);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = FNView(fn,viewer);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Create matrices */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,NULL,&A);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)A,"A");CHKERRQ(ierr);

  if (random) {
    ierr = MatSetRandom(A,NULL);CHKERRQ(ierr);
  } else {
    /* Fill A with a non-symmetric Toeplitz matrix */
    ierr = MatDenseGetArray(A,&As);CHKERRQ(ierr);
    for (i=0;i<n;i++) As[i+i*n]=2.0;
    for (j=1;j<3;j++) {
      for (i=0;i<n-j;i++) {
        As[i+(i+j)*n]=1.0;
        if (!triang) As[(i+j)+i*n]=-1.0;
      }
    }
    As[(n-1)*n] = -5.0;
    As[0] = 2.01;
    ierr = MatDenseRestoreArray(A,&As);CHKERRQ(ierr);
  }
  ierr = TestMatLog(fn,A,viewer,verbose,inplace);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = FNDestroy(&fn);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   testset:
      filter: grep -v "computing matrix functions"
      output_file: output/test13_1.out
      test:
         suffix: 1
         args: -fn_scale .04,2 -n 75
         requires: c99_complex !__float128
      test:
         suffix: 1_triang
         args: -fn_scale .04,2 -n 75 -triang
         requires: c99_complex !__float128
      test:
         suffix: 1_random
         args: -fn_scale .04,2 -n 75 -random
         requires: complex

TEST*/
