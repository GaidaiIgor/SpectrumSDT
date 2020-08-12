/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test MFN interface functions.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in x dimension.\n\n";

#include <slepcmfn.h>

int main(int argc,char **argv)
{
  Mat                  A,B;
  MFN                  mfn;
  FN                   f;
  MFNConvergedReason   reason;
  MFNType              type;
  PetscReal            norm,tol;
  Vec                  v,y;
  PetscInt             N,n=4,Istart,Iend,i,j,II,ncv,its,maxit;
  PetscBool            flg,testprefix=PETSC_FALSE;
  const char           *prefix;
  PetscErrorCode       ierr;
  PetscViewerAndFormat *vf;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  N = n*n;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSquare root of Laplacian y=sqrt(A)*e_1, N=%D (%Dx%D grid)\n\n",N,n,n);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-test_prefix",&testprefix,NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 Compute the discrete 2-D Laplacian, A
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (II=Istart;II<Iend;II++) {
    i = II/n; j = II-i*n;
    if (i>0) { ierr = MatSetValue(A,II,II-n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A,II,II+n,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j>0) { ierr = MatSetValue(A,II,II-1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    if (j<n-1) { ierr = MatSetValue(A,II,II+1,-1.0,INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,II,II,4.0,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);

  ierr = MatCreateVecs(A,NULL,&v);CHKERRQ(ierr);
  ierr = VecSetValue(v,0,1.0,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
  ierr = VecDuplicate(v,&y);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             Create the solver, set the matrix and the function
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MFNCreate(PETSC_COMM_WORLD,&mfn);CHKERRQ(ierr);
  ierr = MFNSetOperator(mfn,A);CHKERRQ(ierr);
  ierr = MFNGetFN(mfn,&f);CHKERRQ(ierr);
  ierr = FNSetType(f,FNSQRT);CHKERRQ(ierr);

  ierr = MFNSetType(mfn,MFNKRYLOV);CHKERRQ(ierr);
  ierr = MFNGetType(mfn,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Type set to %s\n",type);CHKERRQ(ierr);

  /* test prefix usage */
  if (testprefix) {
    ierr = MFNSetOptionsPrefix(mfn,"check_");CHKERRQ(ierr);
    ierr = MFNAppendOptionsPrefix(mfn,"myprefix_");CHKERRQ(ierr);
    ierr = MFNGetOptionsPrefix(mfn,&prefix);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," MFN prefix is currently: %s\n",prefix);CHKERRQ(ierr);
  }

  /* test some interface functions */
  ierr = MFNGetOperator(mfn,&B);CHKERRQ(ierr);
  ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = MFNSetTolerances(mfn,1e-4,500);CHKERRQ(ierr);
  ierr = MFNSetDimensions(mfn,6);CHKERRQ(ierr);
  ierr = MFNSetErrorIfNotConverged(mfn,PETSC_TRUE);CHKERRQ(ierr);
  /* test monitors */
  ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);CHKERRQ(ierr);
  ierr = MFNMonitorSet(mfn,(PetscErrorCode (*)(MFN,PetscInt,PetscReal,void*))MFNMonitorDefault,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);
  /* ierr = MFNMonitorCancel(mfn);CHKERRQ(ierr); */
  ierr = MFNSetFromOptions(mfn);CHKERRQ(ierr);

  /* query properties and print them */
  ierr = MFNGetTolerances(mfn,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Tolerance: %g, max iterations: %D\n",(double)tol,maxit);CHKERRQ(ierr);
  ierr = MFNGetDimensions(mfn,&ncv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Subspace dimension: %D\n",ncv);CHKERRQ(ierr);
  ierr = MFNGetErrorIfNotConverged(mfn,&flg);CHKERRQ(ierr);
  if (flg) { ierr = PetscPrintf(PETSC_COMM_WORLD," Erroring out if convergence fails\n");CHKERRQ(ierr); }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                           Solve  y=sqrt(A)*v
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MFNSolve(mfn,v,y);CHKERRQ(ierr);
  ierr = MFNGetConvergedReason(mfn,&reason);CHKERRQ(ierr);
  ierr = MFNGetIterationNumber(mfn,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Finished - converged reason = %d\n",(int)reason);CHKERRQ(ierr);
  /* ierr = PetscPrintf(PETSC_COMM_WORLD," its = %D\n",its);CHKERRQ(ierr); */
  ierr = VecNorm(y,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," sqrt(A)*v has norm %g\n",(double)norm);CHKERRQ(ierr);

  /*
     Free work space
  */
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
      args: -mfn_monitor_cancel -mfn_converged_reason -mfn_view -log_exclude mfn,bv,fn

   test:
      suffix: 2
      args: -test_prefix -check_myprefix_mfn_monitor
      filter: sed -e "s/estimate [0-9]\.[0-9]*e[+-]\([0-9]*\)/estimate (removed)/g" | sed -e "s/4.0[0-9]*e-10/4.03e-10/"

TEST*/
