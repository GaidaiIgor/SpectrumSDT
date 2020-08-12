/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/svdimpl.h>

static PetscBool SVDPackageInitialized = PETSC_FALSE;

const char *SVDErrorTypes[] = {"ABSOLUTE","RELATIVE","SVDErrorType","SVD_ERROR_",0};
const char *SVDPRIMMEMethods[] = {"","HYBRID","NORMALEQUATIONS","AUGMENTED","SVDPRIMMEMethod","SVD_PRIMME_",0};
const char *const SVDConvergedReasons_Shifted[] = {"","","DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","CONVERGED_USER","SVDConvergedReason","SVD_",0};
const char *const*SVDConvergedReasons = SVDConvergedReasons_Shifted + 4;

/*@C
   SVDFinalizePackage - This function destroys everything in the Slepc interface
   to the SVD package. It is called from SlepcFinalize().

   Level: developer

.seealso: SlepcFinalize()
@*/
PetscErrorCode SVDFinalizePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListDestroy(&SVDList);CHKERRQ(ierr);
  SVDPackageInitialized = PETSC_FALSE;
  SVDRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
   SVDInitializePackage - This function initializes everything in the SVD package.
   It is called from PetscDLLibraryRegister() when using dynamic libraries, and
   on the first call to SVDCreate() when using static libraries.

   Level: developer

.seealso: SlepcInitialize()
@*/
PetscErrorCode SVDInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (SVDPackageInitialized) PetscFunctionReturn(0);
  SVDPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  ierr = PetscClassIdRegister("SVD Solver",&SVD_CLASSID);CHKERRQ(ierr);
  /* Register Constructors */
  ierr = SVDRegisterAll();CHKERRQ(ierr);
  /* Register Events */
  ierr = PetscLogEventRegister("SVDSetUp",SVD_CLASSID,&SVD_SetUp);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("SVDSolve",SVD_CLASSID,&SVD_Solve);CHKERRQ(ierr);
  /* Process Info */
  classids[0] = SVD_CLASSID;
  ierr = PetscInfoProcessClass("svd",1,&classids[0]);CHKERRQ(ierr);
  /* Process summary exclusions */
  ierr = PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt);CHKERRQ(ierr);
  if (opt) {
    ierr = PetscStrInList("svd",logList,',',&pkg);CHKERRQ(ierr);
    if (pkg) { ierr = PetscLogEventDeactivateClass(SVD_CLASSID);CHKERRQ(ierr); }
  }
  /* Register package finalizer */
  ierr = PetscRegisterFinalize(SVDFinalizePackage);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the SVD methods that are in the basic SLEPc libslepcsvd
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcsvd()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SVDInitializePackage();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */

