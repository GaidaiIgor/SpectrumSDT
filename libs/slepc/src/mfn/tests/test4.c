/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates use of MFNSetBV().\n\n"
  "The command line options are:\n"
  "  -t <sval>, where <sval> = scalar value that multiplies the argument.\n"
  "  -file <filename>, where <filename> = matrix file in PETSc binary form.\n\n";

#include <slepcmfn.h>

int main(int argc,char **argv)
{
  Mat                A;           /* problem matrix */
  MFN                mfn;
  FN                 f;
  BV                 V;
  PetscInt           k=8;
  PetscReal          norm;
  PetscScalar        t=2.0;
  Vec                v,y;
  PetscErrorCode     ierr;
  PetscViewer        viewer;
  PetscBool          flg;
  char               filename[PETSC_MAX_PATH_LEN];
  MFNConvergedReason reason;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetScalar(NULL,NULL,"-t",&t,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-k",&k,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nMatrix exponential y=exp(t*A)*e, loaded from file\n\n");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Load matrix A from binary file
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PetscOptionsGetString(NULL,NULL,"-file",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate a file name with the -file option");

#if defined(PETSC_USE_COMPLEX)
  ierr = PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrix from a binary file...\n");CHKERRQ(ierr);
#else
  ierr = PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrix from a binary file...\n");CHKERRQ(ierr);
#endif
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatLoad(A,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* set v = ones(n,1) */
  ierr = MatCreateVecs(A,NULL,&y);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&v);CHKERRQ(ierr);
  ierr = VecSet(v,1.0);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Create the BV object
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = BVCreate(PETSC_COMM_WORLD,&V);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)V,"V");CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(V,v,k);CHKERRQ(ierr);
  ierr = BVSetFromOptions(V);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MFNCreate(PETSC_COMM_WORLD,&mfn);CHKERRQ(ierr);
  ierr = MFNSetOperator(mfn,A);CHKERRQ(ierr);
  ierr = MFNSetBV(mfn,V);CHKERRQ(ierr);
  ierr = MFNGetFN(mfn,&f);CHKERRQ(ierr);
  ierr = FNSetType(f,FNEXP);CHKERRQ(ierr);
  ierr = FNSetScale(f,t,1.0);CHKERRQ(ierr);
  ierr = MFNSetDimensions(mfn,k);CHKERRQ(ierr);
  ierr = MFNSetFromOptions(mfn);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the problem, y=exp(t*A)*v
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MFNSolve(mfn,v,y);CHKERRQ(ierr);
  ierr = MFNGetConvergedReason(mfn,&reason);CHKERRQ(ierr);
  if (reason<0) SETERRQ(PETSC_COMM_WORLD,1,"Solver did not converge");
  ierr = VecNorm(y,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Computed vector at time t=%.4g has norm %g\n\n",(double)PetscRealPart(t),(double)norm);CHKERRQ(ierr);

  /*
     Free work space
  */
  ierr = BVDestroy(&V);CHKERRQ(ierr);
  ierr = MFNDestroy(&mfn);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      args: -file ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -k 12
      requires: double !complex !define(PETSC_USE_64BIT_INDICES)

   test:
      suffix: 2
      args: -file ${DATAFILESPATH}/matrices/complex/qc324.petsc -k 12
      requires: double complex datafilespath !define(PETSC_USE_64BIT_INDICES)

TEST*/
