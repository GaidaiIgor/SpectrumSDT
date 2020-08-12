/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/nepimpl.h>

static PetscBool NEPPackageInitialized = PETSC_FALSE;

const char *NEPErrorTypes[] = {"ABSOLUTE","RELATIVE","BACKWARD","NEPErrorType","NEP_ERROR_",0};
const char *NEPRefineTypes[] = {"NONE","SIMPLE","MULTIPLE","NEPRefine","NEP_REFINE_",0};
const char *NEPRefineSchemes[] = {"","SCHUR","MBE","EXPLICIT","NEPRefineScheme","NEP_REFINE_SCHEME_",0};
const char *const NEPConvergedReasons_Shifted[] = {"DIVERGED_SUBSPACE_EXHAUSTED","DIVERGED_LINEAR_SOLVE","","DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","CONVERGED_USER","NEPConvergedReason","NEP_",0};
const char *const*NEPConvergedReasons = NEPConvergedReasons_Shifted + 5;

/*@C
   NEPFinalizePackage - This function destroys everything in the Slepc interface
   to the NEP package. It is called from SlepcFinalize().

   Level: developer

.seealso: SlepcFinalize()
@*/
PetscErrorCode NEPFinalizePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListDestroy(&NEPList);CHKERRQ(ierr);
  NEPPackageInitialized = PETSC_FALSE;
  NEPRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
   NEPInitializePackage - This function initializes everything in the NEP package.
   It is called from PetscDLLibraryRegister() when using dynamic libraries, and
   on the first call to NEPCreate() when using static libraries.

   Level: developer

.seealso: SlepcInitialize()
@*/
PetscErrorCode NEPInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (NEPPackageInitialized) PetscFunctionReturn(0);
  NEPPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  ierr = PetscClassIdRegister("NEP Solver",&NEP_CLASSID);CHKERRQ(ierr);
  /* Register Constructors */
  ierr = NEPRegisterAll();CHKERRQ(ierr);
  /* Register Events */
  ierr = PetscLogEventRegister("NEPSetUp",NEP_CLASSID,&NEP_SetUp);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("NEPSolve",NEP_CLASSID,&NEP_Solve);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("NEPRefine",NEP_CLASSID,&NEP_Refine);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("NEPFunctionEval",NEP_CLASSID,&NEP_FunctionEval);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("NEPJacobianEval",NEP_CLASSID,&NEP_JacobianEval);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("NEPResolvent",NEP_CLASSID,&NEP_Resolvent);CHKERRQ(ierr);
  /* Process Info */
  classids[0] = NEP_CLASSID;
  ierr = PetscInfoProcessClass("nep",1,&classids[0]);CHKERRQ(ierr);
  /* Process summary exclusions */
  ierr = PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt);CHKERRQ(ierr);
  if (opt) {
    ierr = PetscStrInList("nep",logList,',',&pkg);CHKERRQ(ierr);
    if (pkg) { ierr = PetscLogEventDeactivateClass(NEP_CLASSID);CHKERRQ(ierr); }
  }
  /* Register package finalizer */
  ierr = PetscRegisterFinalize(NEPFinalizePackage);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the NEP methods that are in the basic SLEPc libslepcnep
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcnep()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = NEPInitializePackage();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */

