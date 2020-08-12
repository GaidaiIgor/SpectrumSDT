/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Test setting FN parameters from the command line.\n\n";

#include <slepcfn.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  FN             fn1,fn2;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = FNCreate(PETSC_COMM_WORLD,&fn1);CHKERRQ(ierr);
  ierr = FNSetOptionsPrefix(fn1,"f1_");CHKERRQ(ierr);
  ierr = FNSetFromOptions(fn1);CHKERRQ(ierr);
  ierr = FNView(fn1,NULL);CHKERRQ(ierr);
  ierr = FNDestroy(&fn1);CHKERRQ(ierr);
  ierr = FNCreate(PETSC_COMM_WORLD,&fn2);CHKERRQ(ierr);
  ierr = FNSetOptionsPrefix(fn2,"f2_");CHKERRQ(ierr);
  ierr = FNSetFromOptions(fn2);CHKERRQ(ierr);
  ierr = FNView(fn2,NULL);CHKERRQ(ierr);
  ierr = FNDestroy(&fn2);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1
      nsize: 1
      args: -f1_fn_type exp -f1_fn_scale -2.5 -f2_fn_type rational -f2_fn_rational_numerator -1,1 -f2_fn_rational_denominator 1,-6,4

TEST*/
