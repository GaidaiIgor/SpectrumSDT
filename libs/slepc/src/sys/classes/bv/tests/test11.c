/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test BV block orthogonalization.\n\n";

#include <slepcbv.h>

/*
   Compute the Frobenius norm ||A(l:k,l:k)-diag||_F
 */
PetscErrorCode MyMatNorm(Mat A,PetscInt lda,PetscInt l,PetscInt k,PetscScalar diag,PetscReal *norm)
{
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscScalar    *pA;
  PetscReal      s,val;

  PetscFunctionBeginUser;
  ierr = MatDenseGetArray(A,&pA);CHKERRQ(ierr);
  s = 0.0;
  for (i=l;i<k;i++) {
    for (j=l;j<k;j++) {
      val = (i==j)? PetscAbsScalar(pA[i+j*lda]-diag): PetscAbsScalar(pA[i+j*lda]);
      s += val*val;
    }
  }
  *norm = PetscSqrtReal(s);
  ierr = MatDenseRestoreArray(A,&pA);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  BV             X,Y,Z,cached;
  Mat            B=NULL,M,R=NULL;
  Vec            v,t;
  PetscInt       i,j,n=20,l=2,k=8,Istart,Iend;
  PetscViewer    view;
  PetscBool      withb,resid,rand,verbose;
  PetscReal      norm;
  PetscScalar    alpha;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-l",&l,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-withb",&withb);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-resid",&resid);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-rand",&rand);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test BV block orthogonalization (length %D, l=%D, k=%D)%s.\n",n,l,k,withb?" with non-standard inner product":"");CHKERRQ(ierr);

  /* Create template vector */
  ierr = VecCreate(PETSC_COMM_WORLD,&t);CHKERRQ(ierr);
  ierr = VecSetSizes(t,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(t);CHKERRQ(ierr);

  /* Create BV object X */
  ierr = BVCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)X,"X");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(X,t,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(X);CHKERRQ(ierr);

  /* Set up viewer */
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&view);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscViewerPushFormat(view,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }

  /* Fill X entries */
  if (rand) {
    ierr = BVSetRandom(X);CHKERRQ(ierr);
  } else {
    for (j=0;j<k;j++) {
      ierr = BVGetColumn(X,j,&v);CHKERRQ(ierr);
      ierr = VecSet(v,0.0);CHKERRQ(ierr);
      for (i=0;i<=n/2;i++) {
        if (i+j<n) {
          alpha = (3.0*i+j-2)/(2*(i+j+1));
          ierr = VecSetValue(v,i+j,alpha,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
      ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(X,j,&v);CHKERRQ(ierr);
    }
  }
  if (verbose) {
    ierr = BVView(X,view);CHKERRQ(ierr);
  }

  if (withb) {
    /* Create inner product matrix */
    ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
    ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(B);CHKERRQ(ierr);
    ierr = MatSetUp(B);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)B,"B");CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(B,&Istart,&Iend);CHKERRQ(ierr);
    for (i=Istart;i<Iend;i++) {
      if (i>0) { ierr = MatSetValue(B,i,i-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
      if (i<n-1) { ierr = MatSetValue(B,i,i+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
      ierr = MatSetValue(B,i,i,2.0,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    if (verbose) {
      ierr = MatView(B,view);CHKERRQ(ierr);
    }
    ierr = BVSetMatrix(X,B,PETSC_FALSE);CHKERRQ(ierr);
  }

  /* Create copy on Y */
  ierr = BVDuplicate(X,&Y);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)Y,"Y");CHKERRQ(ierr);
  ierr = BVCopy(X,Y);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&M);CHKERRQ(ierr);

  if (resid) {
    /* Create matrix R to store triangular factor */
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,k,k,NULL,&R);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)R,"R");CHKERRQ(ierr);
  }

  if (l>0) {
    /* First orthogonalize leading columns */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Orthogonalizing leading columns\n");CHKERRQ(ierr);
    ierr = BVSetActiveColumns(Y,0,l);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(X,0,l);CHKERRQ(ierr);
    ierr = BVOrthogonalize(Y,R);CHKERRQ(ierr);
    if (verbose) {
      ierr = BVView(Y,view);CHKERRQ(ierr);
      if (resid) { ierr = MatView(R,view);CHKERRQ(ierr); }
    }

    if (withb) {
      /* Extract cached BV and check it is equal to B*X */
      ierr = BVGetCachedBV(Y,&cached);CHKERRQ(ierr);
      ierr = BVDuplicate(X,&Z);CHKERRQ(ierr);
      ierr = BVSetMatrix(Z,NULL,PETSC_FALSE);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(Z,0,l);CHKERRQ(ierr);
      ierr = BVCopy(X,Z);CHKERRQ(ierr);
      ierr = BVMatMult(X,B,Z);CHKERRQ(ierr);
      ierr = BVMult(Z,-1.0,1.0,cached,NULL);CHKERRQ(ierr);
      ierr = BVNorm(Z,NORM_FROBENIUS,&norm);CHKERRQ(ierr);
      if (norm<100*PETSC_MACHINE_EPSILON) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Difference ||cached-BX|| < 100*eps\n");CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Difference ||cached-BX||: %g\n",(double)norm);CHKERRQ(ierr);
      }
      ierr = BVDestroy(&Z);CHKERRQ(ierr);
    }

    /* Check orthogonality */
    ierr = BVDot(Y,Y,M);CHKERRQ(ierr);
    ierr = MyMatNorm(M,k,0,l,1.0,&norm);CHKERRQ(ierr);
    if (norm<100*PETSC_MACHINE_EPSILON) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q1 < 100*eps\n");CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q1: %g\n",(double)norm);CHKERRQ(ierr);
    }

    if (resid) {
      /* Check residual */
      ierr = BVDuplicate(X,&Z);CHKERRQ(ierr);
      ierr = BVSetMatrix(Z,NULL,PETSC_FALSE);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(Z,0,l);CHKERRQ(ierr);
      ierr = BVCopy(X,Z);CHKERRQ(ierr);
      ierr = BVMult(Z,-1.0,1.0,Y,R);CHKERRQ(ierr);
      ierr = BVNorm(Z,NORM_FROBENIUS,&norm);CHKERRQ(ierr);
      if (norm<100*PETSC_MACHINE_EPSILON) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual ||X1-Q1*R11|| < 100*eps\n");CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual ||X1-Q1*R11||: %g\n",(double)norm);CHKERRQ(ierr);
      }
      ierr = BVDestroy(&Z);CHKERRQ(ierr);
    }

  }

  /* Now orthogonalize the rest of columns */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Orthogonalizing active columns\n");CHKERRQ(ierr);
  ierr = BVSetActiveColumns(Y,l,k);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(X,l,k);CHKERRQ(ierr);
  ierr = BVOrthogonalize(Y,R);CHKERRQ(ierr);
  if (verbose) {
    ierr = BVView(Y,view);CHKERRQ(ierr);
    if (resid) { ierr = MatView(R,view);CHKERRQ(ierr); }
  }

  if (l>0) {
    /* Check orthogonality */
    ierr = BVDot(Y,Y,M);CHKERRQ(ierr);
    ierr = MyMatNorm(M,k,l,k,1.0,&norm);CHKERRQ(ierr);
    if (norm<100*PETSC_MACHINE_EPSILON) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q2 < 100*eps\n");CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q2: %g\n",(double)norm);CHKERRQ(ierr);
    }
  }

  /* Check the complete decomposition */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Overall decomposition\n");CHKERRQ(ierr);
  ierr = BVSetActiveColumns(Y,0,k);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(X,0,k);CHKERRQ(ierr);

  /* Check orthogonality */
  ierr = BVDot(Y,Y,M);CHKERRQ(ierr);
  ierr = MyMatNorm(M,k,0,k,1.0,&norm);CHKERRQ(ierr);
  if (norm<100*PETSC_MACHINE_EPSILON) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q < 100*eps\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of orthogonality of Q: %g\n",(double)norm);CHKERRQ(ierr);
  }

  if (resid) {
    /* Check residual */
    ierr = BVMult(X,-1.0,1.0,Y,R);CHKERRQ(ierr);
    ierr = BVSetMatrix(X,NULL,PETSC_FALSE);CHKERRQ(ierr);
    ierr = BVNorm(X,NORM_FROBENIUS,&norm);CHKERRQ(ierr);
    if (norm<100*PETSC_MACHINE_EPSILON) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual ||X-Q*R|| < 100*eps\n");CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Residual ||X-Q*R||: %g\n",(double)norm);CHKERRQ(ierr);
    }
    ierr = MatDestroy(&R);CHKERRQ(ierr);
  }

  if (B) { ierr = MatDestroy(&B);CHKERRQ(ierr); }
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = BVDestroy(&Y);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 2
      args: -bv_orthog_block gs -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_1.out

   test:
      suffix: 1_cuda
      nsize: 2
      args: -bv_orthog_block gs -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_1.out

   test:
      suffix: 2
      nsize: 2
      args: -bv_orthog_block chol -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_1.out

   test:
      suffix: 2_cuda
      nsize: 2
      args: -bv_orthog_block chol -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_1.out

   test:
      suffix: 3
      nsize: 2
      args: -bv_orthog_block tsqr -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_1.out

   test:
      suffix: 3_cuda
      nsize: 2
      args: -bv_orthog_block tsqr -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_1.out

   test:
      suffix: 4
      nsize: 2
      args: -withb -bv_orthog_block gs -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_4.out

   test:
      suffix: 4_cuda
      nsize: 2
      args: -withb -bv_orthog_block gs -bv_type svec -vec_type cuda -mat_type aijcusparse
      requires: cuda
      output_file: output/test11_4.out

   test:
      suffix: 5
      nsize: 2
      args: -withb -bv_orthog_block chol -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_4.out

   test:
      suffix: 5_cuda
      nsize: 2
      args: -withb -bv_orthog_block chol -bv_type svec -vec_type cuda -mat_type aijcusparse
      requires: cuda
      output_file: output/test11_4.out

   test:
      suffix: 6
      nsize: 2
      args: -resid -bv_orthog_block gs -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_6.out

   test:
      suffix: 6_cuda
      nsize: 2
      args: -resid -bv_orthog_block gs -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_6.out

   test:
      suffix: 7
      nsize: 2
      args: -resid -bv_orthog_block chol -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_6.out

   test:
      suffix: 7_cuda
      nsize: 2
      args: -resid -bv_orthog_block chol -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_6.out

   test:
      suffix: 8
      nsize: 2
      args: -resid -bv_orthog_block tsqr -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_6.out

   test:
      suffix: 8_cuda
      nsize: 2
      args: -resid -bv_orthog_block tsqr -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_6.out

   test:
      suffix: 9
      nsize: 2
      args: -resid -withb -bv_orthog_block gs -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_9.out

   test:
      suffix: 9_cuda
      nsize: 2
      args: -resid -withb -bv_orthog_block gs -bv_type svec -vec_type cuda -mat_type aijcusparse
      requires: cuda
      output_file: output/test11_9.out

   test:
      suffix: 10
      nsize: 2
      args: -resid -withb -bv_orthog_block chol -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_9.out

   test:
      suffix: 10_cuda
      nsize: 2
      args: -resid -withb -bv_orthog_block chol -bv_type svec -vec_type cuda -mat_type aijcusparse
      requires: cuda
      output_file: output/test11_9.out

   test:
      suffix: 11
      nsize: 7
      args: -bv_orthog_block tsqr -bv_type {{vecs contiguous svec mat}shared output}
      requires: !valgrind
      output_file: output/test11_1.out

   test:
      suffix: 11_cuda
      nsize: 7
      args: -bv_orthog_block tsqr -bv_type svec -vec_type cuda
      requires: cuda !valgrind
      output_file: output/test11_1.out

   test:
      suffix: 12
      nsize: 9
      args: -resid -n 180 -l 0 -k 7 -bv_orthog_block tsqr -bv_type {{vecs contiguous svec mat}shared output}
      requires: !single !valgrind
      output_file: output/test11_12.out

   test:
      suffix: 12_cuda
      nsize: 9
      args: -resid -n 180 -l 0 -k 7 -bv_orthog_block tsqr -bv_type svec -vec_type cuda
      requires: !single !valgrind cuda
      output_file: output/test11_12.out

   test:
      suffix: 13
      nsize: 2
      args: -bv_orthog_block tsqrchol -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_1.out

   test:
      suffix: 13_cuda
      nsize: 2
      args: -bv_orthog_block tsqrchol -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_1.out

   test:
      suffix: 14
      nsize: 2
      args: -resid -bv_orthog_block tsqrchol -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_6.out

   test:
      suffix: 14_cuda
      nsize: 2
      args: -resid -bv_orthog_block tsqrchol -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_6.out

   test:
      suffix: 15
      nsize: 2
      args: -bv_orthog_block svqb -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_1.out

   test:
      suffix: 15_cuda
      nsize: 2
      args: -bv_orthog_block svqb -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_1.out

   test:
      suffix: 16
      nsize: 2
      args: -withb -bv_orthog_block svqb -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_4.out

   test:
      suffix: 16_cuda
      nsize: 2
      args: -withb -bv_orthog_block svqb -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_4.out

   test:
      suffix: 17
      nsize: 2
      args: -resid -bv_orthog_block svqb -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_6.out

   test:
      suffix: 17_cuda
      nsize: 2
      args: -resid -bv_orthog_block svqb -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_6.out

   test:
      suffix: 18
      nsize: 2
      args: -resid -withb -bv_orthog_block svqb -bv_type {{vecs contiguous svec mat}shared output}
      output_file: output/test11_9.out

   test:
      suffix: 18_cuda
      nsize: 2
      args: -resid -withb -bv_orthog_block svqb -bv_type svec -vec_type cuda
      requires: cuda
      output_file: output/test11_9.out

TEST*/
