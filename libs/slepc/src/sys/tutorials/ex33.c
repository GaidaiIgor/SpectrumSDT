/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Demonstrates SlepcGetVersionNumber().\n\n";

#include <slepcsys.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  char           version[128];
  PetscInt       major,minor,subminor;
  PetscBool      verbose;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Checking SLEPc version.\n");CHKERRQ(ierr);

  ierr = SlepcGetVersion(version,sizeof(version));CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-verbose",&verbose);CHKERRQ(ierr);
  if (verbose) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Version information:\n%s\n",version);CHKERRQ(ierr);
  }

  ierr = SlepcGetVersionNumber(&major,&minor,&subminor,NULL);CHKERRQ(ierr);
  if (major != SLEPC_VERSION_MAJOR) SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Library major %d does not equal include %d",(int)major,SLEPC_VERSION_MAJOR);
  if (minor != SLEPC_VERSION_MINOR) SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Library minor %d does not equal include %d",(int)minor,SLEPC_VERSION_MINOR);
  if (subminor != SLEPC_VERSION_SUBMINOR) SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"Library subminor %d does not equal include %d",(int)subminor,SLEPC_VERSION_SUBMINOR);

  ierr = SlepcFinalize();
  return ierr;
}

/*TEST

   test:
      suffix: 1

TEST*/
