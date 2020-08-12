static char help[] = "Testing PtAP for SeqMAIJ matrix, P, with SeqAIJ matrix, A.\n\n";

#include <petscmat.h>

int main(int argc,char **argv)
{
  Mat            pA,P,aijP;
  PetscScalar    pa[]={1.,-1.,0.,0.,1.,-1.,0.,0.,1.};
  PetscInt       i,pij[]={0,1,2};
  PetscInt       aij[3][3]={{0,1,2},{3,4,5},{6,7,8}};
  Mat            A,mC,C;
  PetscScalar    one=1.;
  PetscErrorCode ierr;
  PetscMPIInt    size;
  PetscBool      flg;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This is a uniprocessor example only!");

  /* Create MAIJ matrix, P */
  ierr = MatCreate(PETSC_COMM_SELF,&pA);CHKERRQ(ierr);
  ierr = MatSetSizes(pA,3,3,3,3);CHKERRQ(ierr);
  ierr = MatSetType(pA,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(pA);CHKERRQ(ierr);
  ierr = MatSetOption(pA,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatSetValues(pA,3,pij,3,pij,pa,ADD_VALUES);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(pA,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(pA,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateMAIJ(pA,3,&P);CHKERRQ(ierr);
  ierr = MatDestroy(&pA);CHKERRQ(ierr);

  /* Create AIJ equivalent matrix, aijP, for comparison testing */
  ierr = MatConvert(P,MATSEQAIJ,MAT_INITIAL_MATRIX,&aijP);CHKERRQ(ierr);

  /* Create AIJ matrix A */
  ierr = MatCreate(PETSC_COMM_SELF,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,9,9,9,9);CHKERRQ(ierr);
  ierr = MatSetType(A,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatSetValues(A,3,aij[0],3,aij[0],pa,ADD_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(A,3,aij[1],3,aij[1],pa,ADD_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(A,3,aij[2],3,aij[2],pa,ADD_VALUES);CHKERRQ(ierr);
  for (i=0; i<9; i++) {
    ierr = MatSetValue(A,i,i,one,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Perform PtAP_SeqAIJ_SeqAIJ for comparison testing */
  ierr = MatPtAP(A,aijP,MAT_INITIAL_MATRIX,1.,&C);CHKERRQ(ierr);

  /* Perform PtAP_SeqAIJ_SeqMAIJ */
  /* Developer API */
  ierr = MatProductCreate(A,P,NULL,&mC);CHKERRQ(ierr);
  ierr = MatProductSetType(mC,MATPRODUCT_PtAP);CHKERRQ(ierr);
  ierr = MatProductSetAlgorithm(mC,"default");CHKERRQ(ierr);
  ierr = MatProductSetFill(mC,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = MatProductSetFromOptions(mC);CHKERRQ(ierr);
  ierr = MatProductSymbolic(mC);CHKERRQ(ierr);
  ierr = MatProductNumeric(mC);CHKERRQ(ierr);
  ierr = MatProductNumeric(mC);CHKERRQ(ierr);

  /* Check mC = C */
  ierr = MatEqual(C,mC,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"MatProduct C != mC");
  ierr = MatDestroy(&mC);CHKERRQ(ierr);

  /* User API */
  ierr = MatPtAP(A,P,MAT_INITIAL_MATRIX,1.,&mC);CHKERRQ(ierr);
  ierr = MatPtAP(A,P,MAT_REUSE_MATRIX,1.,&mC);CHKERRQ(ierr);

  /* Check mC = C */
  ierr = MatEqual(C,mC,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"MatPtAP C != mC");
  ierr = MatDestroy(&mC);CHKERRQ(ierr);

  /* Cleanup */
  ierr = MatDestroy(&P);CHKERRQ(ierr);
  ierr = MatDestroy(&aijP);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&C);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

   test:
      output_file: output/ex101.out

TEST*/
