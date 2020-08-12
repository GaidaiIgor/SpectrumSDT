/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/pepimpl.h>

static PetscBool PEPPackageInitialized = PETSC_FALSE;

const char *PEPBasisTypes[] = {"MONOMIAL","CHEBYSHEV1","CHEBYSHEV2","LEGENDRE","LAGUERRE","HERMITE","PEPBasis","PEP_BASIS_",0};
const char *PEPScaleTypes[] = {"NONE","SCALAR","DIAGONAL","BOTH","PEPScale","PEP_SCALE_",0};
const char *PEPRefineTypes[] = {"NONE","SIMPLE","MULTIPLE","PEPRefine","PEP_REFINE_",0};
const char *PEPRefineSchemes[] = {"","SCHUR","MBE","EXPLICIT","PEPRefineScheme","PEP_REFINE_SCHEME_",0};
const char *PEPExtractTypes[] = {"","NONE","NORM","RESIDUAL","STRUCTURED","PEPExtract","PEP_EXTRACT_",0};
const char *PEPErrorTypes[] = {"ABSOLUTE","RELATIVE","BACKWARD","PEPErrorType","PEP_ERROR_",0};
const char *const PEPConvergedReasons_Shifted[] = {"","DIVERGED_SYMMETRY_LOST","DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","CONVERGED_USER","PEPConvergedReason","PEP_",0};
const char *const*PEPConvergedReasons = PEPConvergedReasons_Shifted + 4;
const char *PEPJDProjectionTypes[] = {"HARMONIC","ORTHOGONAL","PEPJDProjection","PEP_JD_PROJECTION_",0};

/*@C
   PEPFinalizePackage - This function destroys everything in the Slepc interface
   to the PEP package. It is called from SlepcFinalize().

   Level: developer

.seealso: SlepcFinalize()
@*/
PetscErrorCode PEPFinalizePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListDestroy(&PEPList);CHKERRQ(ierr);
  PEPPackageInitialized = PETSC_FALSE;
  PEPRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
   PEPInitializePackage - This function initializes everything in the PEP package.
   It is called from PetscDLLibraryRegister() when using dynamic libraries, and
   on the first call to PEPCreate() when using static libraries.

   Level: developer

.seealso: SlepcInitialize()
@*/
PetscErrorCode PEPInitializePackage(void)
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (PEPPackageInitialized) PetscFunctionReturn(0);
  PEPPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  ierr = PetscClassIdRegister("PEP Solver",&PEP_CLASSID);CHKERRQ(ierr);
  /* Register Constructors */
  ierr = PEPRegisterAll();CHKERRQ(ierr);
  /* Register Events */
  ierr = PetscLogEventRegister("PEPSetUp",PEP_CLASSID,&PEP_SetUp);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("PEPSolve",PEP_CLASSID,&PEP_Solve);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("PEPRefine",PEP_CLASSID,&PEP_Refine);CHKERRQ(ierr);
  /* Process Info */
  classids[0] = PEP_CLASSID;
  ierr = PetscInfoProcessClass("pep",1,&classids[0]);CHKERRQ(ierr);
  /* Process summary exclusions */
  ierr = PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt);CHKERRQ(ierr);
  if (opt) {
    ierr = PetscStrInList("pep",logList,',',&pkg);CHKERRQ(ierr);
    if (pkg) { ierr = PetscLogEventDeactivateClass(PEP_CLASSID);CHKERRQ(ierr); }
  }
  /* Register package finalizer */
  ierr = PetscRegisterFinalize(PEPFinalizePackage);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the PEP methods that are in the basic SLEPc libslepcpep
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcpep()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PEPInitializePackage();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */

