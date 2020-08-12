
static char help[] = "Test Matrix products for AIJ matrices\n\
Input arguments are:\n\
  -fA <input_file> -fB <input_file> -fC <input_file>: file to load\n\n";
/* Example of usage:
   ./ex62 -fA <A_binary> -fB <B_binary>
   mpiexec -n 3 ./ex62 -fA medium -fB medium
*/

#include <petscmat.h>

/*
     B = A - B
     norm = norm(B)
*/
PetscErrorCode MatNormDifference(Mat A,Mat B,PetscReal *norm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatAXPY(B,-1.0,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatNorm(B,NORM_FROBENIUS,norm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **args)
{
  Mat            A,A_save,B,C,P,C1,R;
  PetscViewer    viewer;
  PetscErrorCode ierr;
  PetscMPIInt    size,rank;
  PetscInt       i,j,*idxn,M,N,nzp,PN,rstart,rend;
  PetscReal      norm;
  PetscRandom    rdm;
  char           file[2][128];
  PetscScalar    *a,rval,alpha;
  PetscBool      Test_MatMatMult=PETSC_TRUE,Test_MatTrMat=PETSC_TRUE,Test_MatMatTr=PETSC_TRUE;
  PetscBool      Test_MatPtAP=PETSC_TRUE,Test_MatRARt=PETSC_TRUE,flg,seqaij;
  MatInfo        info;
  MatType        mattype;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  /*  Load the matrices A_save and B */
  ierr = PetscOptionsGetString(NULL,NULL,"-fA",file[0],sizeof(file[0]),&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate a file name for small matrix A with the -fA option.");
  ierr = PetscOptionsGetString(NULL,NULL,"-fB",file[1],sizeof(file[1]),&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate a file name for small matrix B with the -fB option.");

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[0],FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A_save);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A_save);CHKERRQ(ierr);
  ierr = MatLoad(A_save,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[1],FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetFromOptions(B);CHKERRQ(ierr);
  ierr = MatLoad(B,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = MatGetType(B,&mattype);CHKERRQ(ierr);

  ierr = MatGetSize(B,&M,&N);CHKERRQ(ierr);
  nzp  = PetscMax((PetscInt)(0.1*M),5);
  ierr = PetscMalloc((nzp+1)*(sizeof(PetscInt)+sizeof(PetscScalar)),&idxn);CHKERRQ(ierr);
  a    = (PetscScalar*)(idxn + nzp);

  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rdm);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(rdm);CHKERRQ(ierr);

  /* 1) MatMatMult() */
  /* ----------------*/
  if (Test_MatMatMult) {
    ierr = MatDuplicate(A_save,MAT_COPY_VALUES,&A);CHKERRQ(ierr);

    /* (1.1) Test developer API */
    ierr = MatProductCreate(A,B,NULL,&C);CHKERRQ(ierr);
    ierr = MatSetOptionsPrefix(C,"AB_");CHKERRQ(ierr);
    ierr = MatProductSetType(C,MATPRODUCT_AB);CHKERRQ(ierr);
    ierr = MatProductSetAlgorithm(C,"default");CHKERRQ(ierr);
    ierr = MatProductSetFill(C,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = MatProductSetFromOptions(C);CHKERRQ(ierr);
    ierr = MatProductSymbolic(C);CHKERRQ(ierr);
    ierr = MatProductNumeric(C);CHKERRQ(ierr);

    /* Test reuse symbolic C */
    alpha = 0.9;
    ierr = MatScale(A,alpha);CHKERRQ(ierr);
    ierr = MatProductNumeric(C);CHKERRQ(ierr);

    ierr = MatMatMultEqual(A,B,C,10,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Error in C=A*B");
    ierr = MatDestroy(&C);CHKERRQ(ierr);

    /* (1.2) Test user driver */
    ierr = MatMatMult(A,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);

    /* Test MAT_REUSE_MATRIX - reuse symbolic C */
    alpha = 1.0;
    for (i=0; i<2; i++) {
      alpha -= 0.1;
      ierr   = MatScale(A,alpha);CHKERRQ(ierr);
      ierr   = MatMatMult(A,B,MAT_REUSE_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);
    }
    ierr = MatMatMultEqual(A,B,C,10,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Error: MatMatMult()");
    ierr = MatDestroy(&A);CHKERRQ(ierr);

    /* Test MatProductClear() */
    ierr = MatProductClear(C);CHKERRQ(ierr);
    ierr = MatDestroy(&C);CHKERRQ(ierr);

    /* Test MatMatMult() for dense and aij matrices */
    ierr = MatConvert(A_save,MATDENSE,MAT_INITIAL_MATRIX,&A);CHKERRQ(ierr);
    ierr = MatMatMult(A,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);
    ierr = MatDestroy(&C);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);

  }

  /* Create P and R = P^T  */
  /* --------------------- */
  PN   = M/2;
  nzp  = 5; /* num of nonzeros in each row of P */
  ierr = MatCreate(PETSC_COMM_WORLD,&P);CHKERRQ(ierr);
  ierr = MatSetSizes(P,PETSC_DECIDE,PETSC_DECIDE,M,PN);CHKERRQ(ierr);
  ierr = MatSetType(P,mattype);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(P,nzp,NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(P,nzp,NULL,nzp,NULL);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(P,&rstart,&rend);CHKERRQ(ierr);
  for (i=0; i<nzp; i++) {
    ierr = PetscRandomGetValue(rdm,&a[i]);CHKERRQ(ierr);
  }
  for (i=rstart; i<rend; i++) {
    for (j=0; j<nzp; j++) {
      ierr    = PetscRandomGetValue(rdm,&rval);CHKERRQ(ierr);
      idxn[j] = (PetscInt)(PetscRealPart(rval)*PN);
    }
    ierr = MatSetValues(P,1,&i,nzp,idxn,a,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatTranspose(P,MAT_INITIAL_MATRIX,&R);CHKERRQ(ierr);

  /* 2) MatTransposeMatMult() */
  /* ------------------------ */
  if (Test_MatTrMat) {
    /* (2.1) Test developer driver C = P^T*B */
    ierr = MatProductCreate(P,B,NULL,&C);CHKERRQ(ierr);
    ierr = MatSetOptionsPrefix(C,"AtB_");CHKERRQ(ierr);
    ierr = MatProductSetType(C,MATPRODUCT_AtB);CHKERRQ(ierr);
    ierr = MatProductSetAlgorithm(C,"default");CHKERRQ(ierr);
    ierr = MatProductSetFill(C,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = MatProductSetFromOptions(C);CHKERRQ(ierr);
    ierr = MatProductSymbolic(C);CHKERRQ(ierr); /* equivalent to MatSetUp() */
    ierr = MatSetOption(C,MAT_USE_INODES,PETSC_FALSE);CHKERRQ(ierr); /* illustrate how to call MatSetOption() */
    ierr = MatProductNumeric(C);CHKERRQ(ierr);
    ierr = MatProductNumeric(C);CHKERRQ(ierr); /* test reuse symbolic C */

    ierr = MatTransposeMatMultEqual(P,B,C,10,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Error: developer driver C = P^T*B");
    ierr = MatDestroy(&C);CHKERRQ(ierr);

    /* (2.2) Test user driver C = P^T*B */
    ierr = MatTransposeMatMult(P,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);
    ierr = MatTransposeMatMult(P,B,MAT_REUSE_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);
    ierr = MatGetInfo(C,MAT_GLOBAL_SUM,&info);CHKERRQ(ierr);
    ierr = MatFreeIntermediateDataStructures(C);CHKERRQ(ierr);

    /* Compare P^T*B and R*B */
    ierr = MatMatMult(R,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C1);CHKERRQ(ierr);
    ierr = MatNormDifference(C,C1,&norm);CHKERRQ(ierr);
    if (norm > PETSC_SMALL) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Error in MatTransposeMatMult(): %g",(double)norm);
    ierr = MatDestroy(&C1);CHKERRQ(ierr);

    /* Test MatDuplicate() of C=P^T*B */
    ierr = MatDuplicate(C,MAT_COPY_VALUES,&C1);CHKERRQ(ierr);
    ierr = MatDestroy(&C1);CHKERRQ(ierr);
    ierr = MatDestroy(&C);CHKERRQ(ierr);
  }

  /* 3) MatMatTransposeMult() */
  /* ------------------------ */
  if (Test_MatMatTr) {
    /* C = B*R^T */
    ierr = PetscObjectTypeCompare((PetscObject)B,MATSEQAIJ,&seqaij);CHKERRQ(ierr);
    if (size == 1 && seqaij) {
      ierr = MatMatTransposeMult(B,R,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);
      ierr = MatSetOptionsPrefix(C,"ABt_");CHKERRQ(ierr); /* enable '-ABt_' for matrix C */
      ierr = MatGetInfo(C,MAT_GLOBAL_SUM,&info);CHKERRQ(ierr);

      /* Test MAT_REUSE_MATRIX - reuse symbolic C */
      ierr = MatMatTransposeMult(B,R,MAT_REUSE_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);

      /* Check */
      ierr = MatMatMult(B,P,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C1);CHKERRQ(ierr);
      ierr = MatNormDifference(C,C1,&norm);CHKERRQ(ierr);
      if (norm > PETSC_SMALL) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Error in MatMatTransposeMult() %g",(double)norm);
      ierr = MatDestroy(&C1);CHKERRQ(ierr);
      ierr = MatDestroy(&C);CHKERRQ(ierr);
    }
  }

  /* 4) Test MatPtAP() */
  /*-------------------*/
  if (Test_MatPtAP) {
    ierr = MatDuplicate(A_save,MAT_COPY_VALUES,&A);CHKERRQ(ierr);

    /* (4.1) Test developer API */
    ierr = MatProductCreate(A,P,NULL,&C);CHKERRQ(ierr);
    ierr = MatSetOptionsPrefix(C,"PtAP_");CHKERRQ(ierr);
    ierr = MatProductSetType(C,MATPRODUCT_PtAP);CHKERRQ(ierr);
    ierr = MatProductSetAlgorithm(C,"default");CHKERRQ(ierr);
    ierr = MatProductSetFill(C,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = MatProductSetFromOptions(C);CHKERRQ(ierr);
    ierr = MatProductSymbolic(C);CHKERRQ(ierr);
    ierr = MatProductNumeric(C);CHKERRQ(ierr);
    ierr = MatProductNumeric(C);CHKERRQ(ierr); /* reuse symbolic C */

    ierr = MatPtAPMultEqual(A,P,C,10,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Error in MatProduct_PtAP");
    ierr = MatDestroy(&C);CHKERRQ(ierr);

    /* (4.2) Test user driver */
    ierr = MatPtAP(A,P,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);

    /* Test MAT_REUSE_MATRIX - reuse symbolic C */
    alpha=1.0;
    for (i=0; i<2; i++) {
      alpha -= 0.1;
      ierr   = MatScale(A,alpha);CHKERRQ(ierr);
      ierr   = MatPtAP(A,P,MAT_REUSE_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);
    }

    ierr = MatPtAPMultEqual(A,P,C,10,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Error in MatPtAP");

    /* 5) Test MatRARt() */
    /* ----------------- */
    if (Test_MatRARt) {
      Mat RARt;
      ierr = MatTranspose(P,MAT_REUSE_MATRIX,&R);CHKERRQ(ierr);

      /* (5.1) Test developer driver RARt = R*A*Rt */
      ierr = MatProductCreate(A,R,NULL,&RARt);CHKERRQ(ierr);
      ierr = MatSetOptionsPrefix(RARt,"RARt_");CHKERRQ(ierr);
      ierr = MatProductSetType(RARt,MATPRODUCT_RARt);CHKERRQ(ierr);
      ierr = MatProductSetAlgorithm(RARt,"default");CHKERRQ(ierr);
      ierr = MatProductSetFill(RARt,PETSC_DEFAULT);CHKERRQ(ierr);
      ierr = MatProductSetFromOptions(RARt);CHKERRQ(ierr);
      ierr = MatProductSymbolic(RARt);CHKERRQ(ierr); /* equivalent to MatSetUp() */
      ierr = MatSetOption(RARt,MAT_USE_INODES,PETSC_FALSE);CHKERRQ(ierr); /* illustrate how to call MatSetOption() */
      ierr = MatProductNumeric(RARt);CHKERRQ(ierr);
      ierr = MatProductNumeric(RARt);CHKERRQ(ierr); /* test reuse symbolic RARt */
      ierr = MatDestroy(&RARt);CHKERRQ(ierr);

      /* (2.2) Test user driver RARt = R*A*Rt */
      ierr = MatRARt(A,R,MAT_INITIAL_MATRIX,2.0,&RARt);CHKERRQ(ierr);
      ierr = MatRARt(A,R,MAT_REUSE_MATRIX,2.0,&RARt);CHKERRQ(ierr);

      ierr = MatNormDifference(C,RARt,&norm);CHKERRQ(ierr);
      if (norm > PETSC_SMALL) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_PLIB,"|PtAP - RARt| = %g",(double)norm);
      ierr = MatDestroy(&RARt);CHKERRQ(ierr);
    }

    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = MatDestroy(&C);CHKERRQ(ierr);
  }

  /* Destroy objects */
  ierr = PetscRandomDestroy(&rdm);CHKERRQ(ierr);
  ierr = PetscFree(idxn);CHKERRQ(ierr);

  ierr = MatDestroy(&A_save);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = MatDestroy(&P);CHKERRQ(ierr);
  ierr = MatDestroy(&R);CHKERRQ(ierr);

  PetscFinalize();
  return ierr;
}

/*TEST
   test:
     suffix: 1
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium
     output_file: output/ex62_1.out

   test:
     suffix: 2_ab_scalable
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via scalable -matmatmult_via scalable -AtB_matproduct_atb_via outerproduct -mattransposematmult_via outerproduct
     output_file: output/ex62_1.out

   test:
     suffix: 3_ab_scalable_fast
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via scalable_fast -matmatmult_via scalable_fast -matmattransmult_via color
     output_file: output/ex62_1.out

   test:
     suffix: 4_ab_heap
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via heap -matmatmult_via heap -PtAP_matproduct_ptap_via rap -matptap_via rap
     output_file: output/ex62_1.out

   test:
     suffix: 5_ab_btheap
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via btheap -matmatmult_via btheap -matrart_via r*art
     output_file: output/ex62_1.out

   test:
     suffix: 6_ab_llcondensed
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via llcondensed -matmatmult_via llcondensed -matrart_via coloring_rart
     output_file: output/ex62_1.out

   test:
     suffix: 7_ab_rowmerge
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via rowmerge -matmatmult_via rowmerge
     output_file: output/ex62_1.out

   test:
     suffix: 8_ab_hypre
     requires: hypre datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via hypre -matmatmult_via hypre -PtAP_matproduct_ptap_via hypre -matptap_via hypre
     output_file: output/ex62_1.out

   test:
     suffix: 10
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     nsize: 3
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium
     output_file: output/ex62_1.out

   test:
     suffix: 11_ab_scalable
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     nsize: 3
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via scalable -matmatmult_via scalable -AtB_matproduct_atb_via scalable -mattransposematmult_via scalable
     output_file: output/ex62_1.out

   test:
     suffix: 12_ab_seqmpi
     requires: datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     nsize: 3
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via seqmpi -matmatmult_via seqmpi -AtB_matproduct_atb_via at*b -mattransposematmult_via at*b
     output_file: output/ex62_1.out

   test:
     suffix: 13_ab_hypre
     requires: hypre datafilespath !complex double !define(PETSC_USE_64BIT_INDICES)
     nsize: 3
     args: -fA ${DATAFILESPATH}/matrices/medium -fB ${DATAFILESPATH}/matrices/medium -AB_matproduct_ab_via hypre -matmatmult_via hypre -PtAP_matproduct_ptap_via hypre -matptap_via hypre
     output_file: output/ex62_1.out

TEST*/
