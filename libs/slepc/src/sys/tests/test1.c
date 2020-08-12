/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Demonstrates SlepcInitializeNoArguments()
*/

#include <slepcsys.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscBool      isInitialized;

  ierr = SlepcInitialized(&isInitialized);if (ierr) return ierr;
  if (!isInitialized) {
    ierr = SlepcInitializeNoArguments();if (ierr) return ierr;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Initialize SLEPc.\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"SLEPc was already initialized.\n");CHKERRQ(ierr);
  }
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1

TEST*/
