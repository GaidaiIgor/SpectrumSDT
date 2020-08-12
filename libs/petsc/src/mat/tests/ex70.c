#include <petscmat.h>

static char help[] = "Tests MatMatMult with MAT_REUSE_MATRIX and already allocated dense result.\n\n";

static PetscErrorCode CheckLocal(Mat A, Mat B, PetscScalar *a, PetscScalar *b)
{
  PetscErrorCode ierr;
  PetscBool      wA = PETSC_FALSE, wB = PETSC_FALSE;

  PetscFunctionBegin;
  if (a) {
    const PetscScalar *Aa;
    ierr = MatDenseGetArrayRead(A,&Aa);CHKERRQ(ierr);
    wA   = (PetscBool)(a != Aa);
    ierr = MatDenseRestoreArrayRead(A,&Aa);CHKERRQ(ierr);
  }
  if (b) {
    const PetscScalar *Bb;
    ierr = MatDenseGetArrayRead(B,&Bb);CHKERRQ(ierr);
    wB   = (PetscBool)(b != Bb);
    ierr = MatDenseRestoreArrayRead(B,&Bb);CHKERRQ(ierr);
  }
  if (wA || wB) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Wrong array in A Mat %d, Wrong array in B Mat %d",wA,wB);
  PetscFunctionReturn(0);
}

int main(int argc,char **args)
{
  Mat            X,B,A,T;
  PetscInt       m,n,M = 10,N = 10,K = 5, ldx = 0, ldb = 0;
  const char     *deft = MATAIJ;
  char           mattype[256];
  PetscBool      flg,symm = PETSC_FALSE,testtt = PETSC_TRUE, testnest = PETSC_TRUE, testtranspose = PETSC_TRUE, testcircular = PETSC_FALSE, local = PETSC_TRUE;
  PetscScalar    *dataX = NULL,*dataB = NULL;
  PetscScalar    *aX,*aB;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc,&args,NULL,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-K",&K,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-symm",&symm,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-local",&local,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ldx",&ldx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ldb",&ldb,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-testtranspose",&testtranspose,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-testnest",&testnest,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-testtt",&testtt,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-testcircular",&testcircular,NULL);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  testtranspose = PETSC_FALSE;
  testtt = PETSC_FALSE;
#endif
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetType(A,MATAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = MatSetRandom(A,NULL);CHKERRQ(ierr);
  if (M==N && symm) {
    Mat AT;

    ierr = MatHermitianTranspose(A,MAT_INITIAL_MATRIX,&AT);CHKERRQ(ierr);
    ierr = MatAXPY(A,1.0,AT,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = MatDestroy(&AT);CHKERRQ(ierr);
    ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
  }
  ierr = MatViewFromOptions(A,NULL,"-A_init_view");CHKERRQ(ierr);
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","","");CHKERRQ(ierr);
  ierr = PetscOptionsFList("-mat_type","Matrix type","MatSetType",MatList,deft,mattype,256,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (flg) {
    Mat A2;

    ierr = MatConvert(A,mattype,MAT_INITIAL_MATRIX,&A2);CHKERRQ(ierr);
    ierr = MatMultEqual(A,A2,10,&flg);CHKERRQ(ierr);
    if (!flg) {
      Mat AE,A2E;

      ierr = PetscPrintf(PETSC_COMM_WORLD,"Error with convert\n");CHKERRQ(ierr);
      ierr = MatComputeOperator(A,MATDENSE,&AE);CHKERRQ(ierr);
      ierr = MatComputeOperator(A2,MATDENSE,&A2E);CHKERRQ(ierr);
      ierr = MatView(AE,NULL);CHKERRQ(ierr);
      ierr = MatView(A2E,NULL);CHKERRQ(ierr);
      ierr = MatAXPY(A2E,-1.0,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
      ierr = MatView(A2E,NULL);CHKERRQ(ierr);
      ierr = MatDestroy(&A2E);CHKERRQ(ierr);
      ierr = MatDestroy(&AE);CHKERRQ(ierr);
    }
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    A = A2;
  }
  ierr = MatViewFromOptions(A,NULL,"-A_view");CHKERRQ(ierr);

  ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
  if (local) {
    ierr = PetscMalloc1(PetscMax((m+ldx)*K,1),&dataX);CHKERRQ(ierr);
    ierr = PetscMalloc1(PetscMax((n+ldb)*K,1),&dataB);CHKERRQ(ierr);
  }
  ierr = MatCreateDense(PETSC_COMM_WORLD,n,PETSC_DECIDE,N,K,dataB,&B);CHKERRQ(ierr);
  ierr = MatCreateDense(PETSC_COMM_WORLD,m,PETSC_DECIDE,M,K,dataX,&X);CHKERRQ(ierr);

  /* store pointer to dense data for testing */
  ierr = MatDenseGetArrayRead(B,(const PetscScalar**)&dataB);CHKERRQ(ierr);
  ierr = MatDenseGetArrayRead(X,(const PetscScalar**)&dataX);CHKERRQ(ierr);
  aX   = dataX;
  aB   = dataB;
  ierr = MatDenseRestoreArrayRead(B,(const PetscScalar**)&dataB);CHKERRQ(ierr);
  ierr = MatDenseRestoreArrayRead(X,(const PetscScalar**)&dataX);CHKERRQ(ierr);
  if (local) {
    dataX = aX;
    dataB = aB;
  }
  ierr = MatSetRandom(B,NULL);CHKERRQ(ierr);

  /* Test reusing a previously allocated dense buffer */
  ierr = MatMatMult(A,B,MAT_REUSE_MATRIX,PETSC_DEFAULT,&X);CHKERRQ(ierr);
  ierr = CheckLocal(B,X,aB,aX);CHKERRQ(ierr);
  ierr = MatMatMultEqual(A,B,X,10,&flg);CHKERRQ(ierr);
  if (!flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Error with reusage\n");CHKERRQ(ierr);
    ierr = MatMatMult(A,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&T);CHKERRQ(ierr);
    ierr = MatAXPY(T,-1.0,X,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = MatView(T,NULL);CHKERRQ(ierr);
    ierr = MatDestroy(&T);CHKERRQ(ierr);
  }

  if (testnest) { /* test with MatNest */
    Mat NA;

    ierr = MatCreateNest(PETSC_COMM_WORLD,1,NULL,1,NULL,&A,&NA);CHKERRQ(ierr);
    ierr = MatViewFromOptions(NA,NULL,"-NA_view");CHKERRQ(ierr);
    ierr = MatMatMult(NA,B,MAT_REUSE_MATRIX,PETSC_DEFAULT,&X);CHKERRQ(ierr);
    ierr = CheckLocal(B,X,aB,aX);CHKERRQ(ierr);
    ierr = MatMatMultEqual(NA,B,X,10,&flg);CHKERRQ(ierr);
    if (!flg) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Error with Nest\n");CHKERRQ(ierr);
      ierr = MatMatMult(NA,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&T);CHKERRQ(ierr);
      ierr = MatAXPY(T,-1.0,X,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
      ierr = MatView(T,NULL);CHKERRQ(ierr);
      ierr = MatDestroy(&T);CHKERRQ(ierr);
    }
    ierr = MatDestroy(&NA);CHKERRQ(ierr);
  }

  if (testtranspose) { /* test with Transpose */
    Mat TA;

    ierr = MatCreateHermitianTranspose(A,&TA);CHKERRQ(ierr);
    ierr = MatMatMult(TA,X,MAT_REUSE_MATRIX,PETSC_DEFAULT,&B);CHKERRQ(ierr);
    ierr = CheckLocal(B,X,aB,aX);CHKERRQ(ierr);
    ierr = MatMatMultEqual(TA,X,B,10,&flg);CHKERRQ(ierr);
    if (!flg) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Error with Transpose\n");CHKERRQ(ierr);
      ierr = MatMatMult(TA,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&T);CHKERRQ(ierr);
      ierr = MatAXPY(T,-1.0,X,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
      ierr = MatView(T,NULL);CHKERRQ(ierr);
      ierr = MatDestroy(&T);CHKERRQ(ierr);
    }
    ierr = MatDestroy(&TA);CHKERRQ(ierr);
  }

  if (testtt) { /* test with Transpose(Transpose) */
    Mat TA, TTA;

    ierr = MatCreateHermitianTranspose(A,&TA);CHKERRQ(ierr);
    ierr = MatCreateHermitianTranspose(TA,&TTA);CHKERRQ(ierr);
    ierr = MatMatMult(TTA,B,MAT_REUSE_MATRIX,PETSC_DEFAULT,&X);CHKERRQ(ierr);
    ierr = CheckLocal(B,X,aB,aX);CHKERRQ(ierr);
    ierr = MatMatMultEqual(TTA,B,X,10,&flg);CHKERRQ(ierr);
    if (!flg) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Error with Transpose(Transpose)\n");CHKERRQ(ierr);
      ierr = MatMatMult(TTA,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&T);CHKERRQ(ierr);
      ierr = MatAXPY(T,-1.0,X,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
      ierr = MatView(T,NULL);CHKERRQ(ierr);
      ierr = MatDestroy(&T);CHKERRQ(ierr);
    }
    ierr = MatDestroy(&TA);CHKERRQ(ierr);
    ierr = MatDestroy(&TTA);CHKERRQ(ierr);
  }

  if (testcircular) { /* test circular */
    Mat AB;

    ierr = MatMatMult(A,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&AB);CHKERRQ(ierr);
    ierr = MatMatMult(A,B,MAT_REUSE_MATRIX,PETSC_DEFAULT,&X);CHKERRQ(ierr);
    ierr = CheckLocal(B,X,aB,aX);CHKERRQ(ierr);
    if (M == N && N == K) {
      ierr = MatMatMult(A,X,MAT_REUSE_MATRIX,PETSC_DEFAULT,&B);CHKERRQ(ierr);
    } else {
      ierr = MatTransposeMatMult(A,X,MAT_REUSE_MATRIX,PETSC_DEFAULT,&B);CHKERRQ(ierr);
    }
    ierr = CheckLocal(B,X,aB,aX);CHKERRQ(ierr);
    ierr = MatDestroy(&AB);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&X);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = PetscFree(dataX);CHKERRQ(ierr);
  ierr = PetscFree(dataB);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

  test:
    suffix: 1
    args: -local {{0 1}}

  test:
    output_file: output/ex70_1.out
    nsize: 2
    suffix: 1_par
    args: -testtranspose 0 -local {{0 1}}

  test:
    TODO: MPIAIJ x MPIDENSE broken for MatTransposeMatMult
    output_file: output/ex70_1.out
    nsize: 2
    suffix: 1_par_broken
    args: -testtranspose -local {{0 1}}

  test:
    output_file: output/ex70_1.out
    suffix: 2
    nsize: 1
    args: -M {{7 11}} -N {{12 9}} -K {{1 3}} -local {{0 1}}

  test:
    output_file: output/ex70_1.out
    suffix: 2_par
    nsize: 2
    args: -M {{7 11}} -N {{12 9}} -K {{1 3}} -local {{0 1}} -testcircular 0

  test:
    TODO: MatTransposeMatMultSymbolic_SeqAIJ_SeqDense plays with the destroy routine
    output_file: output/ex70_1.out
    suffix: 2_broken
    nsize: 2
    args: -M {{7 11}} -N {{12 9}} -K {{1 3}} -local {{0 1}} -testcircular 1

  test:
    output_file: output/ex70_1.out
    suffix: 3
    nsize: 1
    args: -M 13 -N 13 -K {{1 3}} -local {{0 1}} -mat_type sbaij -symm -testtranspose 0

  test:
    output_file: output/ex70_1.out
    suffix: 4
    nsize: 1
    args: -M 3 -N 3 -K 3 -local {{0 1}} -testcircular

  test:
    output_file: output/ex70_1.out
    suffix: 5
    nsize: {{2 4}}
    args: -M 3 -N 3 -K 3 -local 1 -testcircular -testtranspose 0

  test:
    TODO: MatCreateDense broken with some processors not having local rows
    output_file: output/ex70_1.out
    suffix: 5_broken
    nsize: {{2 4}}
    args: -M 3 -N 3 -K 3 -local 0 -testcircular -testtranspose 0

  test:
    output_file: output/ex70_1.out
    suffix: 6
    nsize: 1
    args: -M {{1 3}} -N {{2 5}} -K {{1 2}} -local {{0 1}} -testcircular -testtranspose 0

  test:
    TODO: MatTransposeMatMultSymbolic_SeqAIJ_SeqDense plays with the destroy routine
    output_file: output/ex70_1.out
    suffix: 6_broken
    nsize: 1
    args: -M {{1 3}} -N {{2 5}} -K {{1 2}} -local {{0 1}} -testcircular -testtranspose

  test:
    output_file: output/ex70_1.out
    suffix: 7
    nsize: 1
    args: -M 13 -N 13 -K {{1 3}} -local {{0 1}} -mat_type dense -testnest 0 -testcircular

  test:
    TODO: NEST x DENSE with dense nested matrices seems broken in this case
    output_file: output/ex70_1.out
    suffix: 7_broken
    nsize: 1
    args: -M 13 -N 13 -K {{1 3}} -local {{0 1}} -mat_type dense -testnest -testcircular
TEST*/
