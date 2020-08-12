/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/epsimpl.h>

static PetscBool EPSPackageInitialized = PETSC_FALSE;

const char *EPSBalanceTypes[] = {"NONE","ONESIDE","TWOSIDE","USER","EPSBalance","EPS_BALANCE_",0};
const char *EPSErrorTypes[] = {"ABSOLUTE","RELATIVE","BACKWARD","EPSErrorType","EPS_ERROR_",0};
const char *EPSPowerShiftTypes[] = {"CONSTANT","RAYLEIGH","WILKINSON","EPSPowerShiftType","EPS_POWER_SHIFT_",0};
const char *EPSLanczosReorthogTypes[] = {"LOCAL","FULL","SELECTIVE","PERIODIC","PARTIAL","DELAYED","EPSLanczosReorthogType","EPS_LANCZOS_REORTHOG_",0};
const char *EPSPRIMMEMethods[] = {"","DYNAMIC","DEFAULT_MIN_TIME","DEFAULT_MIN_MATVECS","ARNOLDI","GD","GD_PLUSK","GD_OLSEN_PLUSK","JD_OLSEN_PLUSK","RQI","JDQR","JDQMR","JDQMR_ETOL","SUBSPACE_ITERATION","LOBPCG_ORTHOBASIS","LOBPCG_ORTHOBASISW","EPSPRIMMEMethod","EPS_PRIMME_",0};
const char *EPSCISSQuadRules[] = {"(not set yet)","TRAPEZOIDAL","CHEBYSHEV","EPSCISSQuadRule","EPS_CISS_QUADRULE_",0};
const char *EPSCISSExtractions[] = {"RITZ","HANKEL","EPSCISSExtraction","EPS_CISS_EXTRACTION_",0};
const char *const EPSConvergedReasons_Shifted[] = {"","DIVERGED_SYMMETRY_LOST","DIVERGED_BREAKDOWN","DIVERGED_ITS","CONVERGED_ITERATING","CONVERGED_TOL","CONVERGED_USER","EPSConvergedReason","EPS_",0};
const char *const*EPSConvergedReasons = EPSConvergedReasons_Shifted + 4;

/*@C
  EPSFinalizePackage - This function destroys everything in the SLEPc interface
  to the EPS package. It is called from SlepcFinalize().

  Level: developer

.seealso: SlepcFinalize()
@*/
PetscErrorCode EPSFinalizePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListDestroy(&EPSList);CHKERRQ(ierr);
  EPSPackageInitialized = PETSC_FALSE;
  EPSRegisterAllCalled  = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
  EPSInitializePackage - This function initializes everything in the EPS package.
  It is called from PetscDLLibraryRegister() when using dynamic libraries, and
  on the first call to EPSCreate() when using static libraries.

  Level: developer

.seealso: SlepcInitialize()
@*/
PetscErrorCode EPSInitializePackage()
{
  char           logList[256];
  PetscBool      opt,pkg;
  PetscClassId   classids[1];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (EPSPackageInitialized) PetscFunctionReturn(0);
  EPSPackageInitialized = PETSC_TRUE;
  /* Register Classes */
  ierr = PetscClassIdRegister("EPS Solver",&EPS_CLASSID);CHKERRQ(ierr);
  /* Register Constructors */
  ierr = EPSRegisterAll();CHKERRQ(ierr);
  /* Register Events */
  ierr = PetscLogEventRegister("EPSSetUp",EPS_CLASSID,&EPS_SetUp);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("EPSSolve",EPS_CLASSID,&EPS_Solve);CHKERRQ(ierr);
  /* Process Info */
  classids[0] = EPS_CLASSID;
  ierr = PetscInfoProcessClass("eps",1,&classids[0]);CHKERRQ(ierr);
  /* Process summary exclusions */
  ierr = PetscOptionsGetString(NULL,NULL,"-log_exclude",logList,sizeof(logList),&opt);CHKERRQ(ierr);
  if (opt) {
    ierr = PetscStrInList("eps",logList,',',&pkg);CHKERRQ(ierr);
    if (pkg) { ierr = PetscLogEventDeactivateClass(EPS_CLASSID);CHKERRQ(ierr); }
  }
  /* Register package finalizer */
  ierr = PetscRegisterFinalize(EPSFinalizePackage);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the EPS methods that are in the basic SLEPc libslepceps
  library.
 */
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepceps()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = EPSInitializePackage();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */

