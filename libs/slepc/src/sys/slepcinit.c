/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/slepcimpl.h>           /*I "slepcsys.h" I*/

#if defined(SLEPC_HAVE_HPDDM)
#include <petscksp.h>
SLEPC_EXTERN PetscErrorCode KSPCreate_HPDDM(KSP);
SLEPC_EXTERN PetscErrorCode PCCreate_HPDDM(PC);
#endif

/*@C
    SlepcGetVersion - Gets the SLEPc version information in a string.

    Not collective

    Input Parameter:
.   len - length of the string

    Output Parameter:
.   version - version string

    Level: developer

.seealso: SlepcGetVersionNumber()
@*/
PetscErrorCode SlepcGetVersion(char version[],size_t len)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
#if (SLEPC_VERSION_RELEASE == 1)
  ierr = PetscSNPrintf(version,len,"SLEPc Release Version %d.%d.%d, %s",SLEPC_VERSION_MAJOR,SLEPC_VERSION_MINOR,SLEPC_VERSION_SUBMINOR,SLEPC_VERSION_DATE);CHKERRQ(ierr);
#else
  ierr = PetscSNPrintf(version,len,"SLEPc Development GIT revision: %s  GIT Date: %s",SLEPC_VERSION_GIT,SLEPC_VERSION_DATE_GIT);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

/*@C
    SlepcGetVersionNumber - Gets the SLEPc version information from the library.

    Not collective

    Output Parameters:
+   major    - the major version
.   minor    - the minor version
.   subminor - the subminor version (patch number)
-   release  - indicates the library is from a release

    Notes:
    Pass NULL in any argument that is not requested.

    The C macros SLEPC_VERSION_MAJOR, SLEPC_VERSION_MINOR, SLEPC_VERSION_SUBMINOR,
    SLEPC_VERSION_RELEASE provide the information at compile time. This can be used to confirm
    that the shared library being loaded at runtime has the appropriate version updates.

    This function can be called before SlepcInitialize().

    Level: developer

.seealso: SlepcGetVersion(), SlepcInitialize()
@*/
PetscErrorCode SlepcGetVersionNumber(PetscInt *major,PetscInt *minor,PetscInt *subminor,PetscInt *release)
{
  if (major)    *major    = SLEPC_VERSION_MAJOR;
  if (minor)    *minor    = SLEPC_VERSION_MINOR;
  if (subminor) *subminor = SLEPC_VERSION_SUBMINOR;
  if (release)  *release  = SLEPC_VERSION_RELEASE;
  return 0;
}

/*
   SlepcPrintVersion - Prints SLEPc version info.

   Collective
*/
static PetscErrorCode SlepcPrintVersion(MPI_Comm comm)
{
  PetscErrorCode ierr;
  char           version[256];

  PetscFunctionBegin;
  ierr = SlepcGetVersion(version,256);CHKERRQ(ierr);
  ierr = (*PetscHelpPrintf)(comm,"%s\n",version);CHKERRQ(ierr);
  ierr = (*PetscHelpPrintf)(comm,SLEPC_AUTHOR_INFO);CHKERRQ(ierr);
  ierr = (*PetscHelpPrintf)(comm,"See docs/manual.html for help.\n");CHKERRQ(ierr);
  ierr = (*PetscHelpPrintf)(comm,"SLEPc libraries linked from %s\n",SLEPC_LIB_DIR);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   SlepcPrintHelpIntro - Prints introductory SLEPc help info.

   Collective
*/
static PetscErrorCode SlepcPrintHelpIntro(MPI_Comm comm)
{
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = (*PetscHelpPrintf)(comm,"SLEPc help information includes that for the PETSc libraries, which provide\n");CHKERRQ(ierr);
  ierr = (*PetscHelpPrintf)(comm,"low-level system infrastructure and linear algebra tools.\n");CHKERRQ(ierr);
  ierr = (*PetscHelpPrintf)(comm,"----------------------------------------\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------Nasty global variables -------------------------------*/
/*
   Indicates whether SLEPc started PETSc, or whether it was
   already started before SLEPc was initialized.
*/
PetscBool SlepcBeganPetsc = PETSC_FALSE;
PetscBool SlepcInitializeCalled = PETSC_FALSE;

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES) && defined(PETSC_USE_SHARED_LIBRARIES)
static PetscErrorCode SlepcLoadDynamicLibrary(const char *name,PetscBool *found)
{
  char           libs[PETSC_MAX_PATH_LEN],dlib[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscStrcpy(libs,SLEPC_LIB_DIR);CHKERRQ(ierr);
  ierr = PetscStrcat(libs,"/libslepc");CHKERRQ(ierr);
  ierr = PetscStrcat(libs,name);CHKERRQ(ierr);
  ierr = PetscDLLibraryRetrieve(PETSC_COMM_WORLD,libs,dlib,1024,found);CHKERRQ(ierr);
  if (*found) {
    ierr = PetscDLLibraryAppend(PETSC_COMM_WORLD,&PetscDLLibrariesLoaded,dlib);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
#endif

#if defined(PETSC_HAVE_THREADSAFETY)
SLEPC_EXTERN PetscErrorCode STInitializePackage(void);
SLEPC_EXTERN PetscErrorCode DSInitializePackage(void);
SLEPC_EXTERN PetscErrorCode FNInitializePackage(void);
SLEPC_EXTERN PetscErrorCode BVInitializePackage(void);
SLEPC_EXTERN PetscErrorCode RGInitializePackage(void);
SLEPC_EXTERN PetscErrorCode EPSInitializePackage(void);
SLEPC_EXTERN PetscErrorCode SVDInitializePackage(void);
SLEPC_EXTERN PetscErrorCode PEPInitializePackage(void);
SLEPC_EXTERN PetscErrorCode NEPInitializePackage(void);
SLEPC_EXTERN PetscErrorCode MFNInitializePackage(void);
SLEPC_EXTERN PetscErrorCode LMEInitializePackage(void);
#endif

/*
    SlepcInitialize_DynamicLibraries - Adds the default dynamic link libraries to the
    search path.
*/
PetscErrorCode SlepcInitialize_DynamicLibraries(void)
{
#if (defined(PETSC_HAVE_DYNAMIC_LIBRARIES) && defined(PETSC_USE_SHARED_LIBRARIES)) || defined(PETSC_HAVE_THREADSAFETY)
  PetscErrorCode ierr;
#endif
#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES) && defined(PETSC_USE_SHARED_LIBRARIES)
  PetscBool      found,preload;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES) && defined(PETSC_USE_SHARED_LIBRARIES)
  preload = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-dynamic_library_preload",&preload,NULL);CHKERRQ(ierr);
  if (preload) {
#if defined(PETSC_USE_SINGLE_LIBRARY)
    ierr = SlepcLoadDynamicLibrary("",&found);CHKERRQ(ierr);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to locate SLEPc dynamic library\nYou cannot move the dynamic libraries!");
#else
    ierr = SlepcLoadDynamicLibrary("sys",&found);CHKERRQ(ierr);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to locate SLEPc sys dynamic library\nYou cannot move the dynamic libraries!");
    ierr = SlepcLoadDynamicLibrary("eps",&found);CHKERRQ(ierr);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to locate SLEPc EPS dynamic library\nYou cannot move the dynamic libraries!");
    ierr = SlepcLoadDynamicLibrary("pep",&found);CHKERRQ(ierr);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to locate SLEPc PEP dynamic library\nYou cannot move the dynamic libraries!");
    ierr = SlepcLoadDynamicLibrary("nep",&found);CHKERRQ(ierr);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to locate SLEPc NEP dynamic library\nYou cannot move the dynamic libraries!");
    ierr = SlepcLoadDynamicLibrary("svd",&found);CHKERRQ(ierr);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to locate SLEPc SVD dynamic library\nYou cannot move the dynamic libraries!");
    ierr = SlepcLoadDynamicLibrary("mfn",&found);CHKERRQ(ierr);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to locate SLEPc MFN dynamic library\nYou cannot move the dynamic libraries!");
    ierr = SlepcLoadDynamicLibrary("lme",&found);CHKERRQ(ierr);
    if (!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to locate SLEPc LME dynamic library\nYou cannot move the dynamic libraries!");
#endif
  }
#endif

#if defined(PETSC_HAVE_THREADSAFETY)
  ierr = STInitializePackage();CHKERRQ(ierr);
  ierr = DSInitializePackage();CHKERRQ(ierr);
  ierr = FNInitializePackage();CHKERRQ(ierr);
  ierr = BVInitializePackage();CHKERRQ(ierr);
  ierr = RGInitializePackage();CHKERRQ(ierr);
  ierr = EPSInitializePackage();CHKERRQ(ierr);
  ierr = SVDInitializePackage();CHKERRQ(ierr);
  ierr = PEPInitializePackage();CHKERRQ(ierr);
  ierr = NEPInitializePackage();CHKERRQ(ierr);
  ierr = MFNInitializePackage();CHKERRQ(ierr);
  ierr = LMEInitializePackage();CHKERRQ(ierr);
#endif

#if defined(SLEPC_HAVE_HPDDM)
  ierr = KSPRegister(KSPHPDDM,KSPCreate_HPDDM);CHKERRQ(ierr);
  ierr = PCRegister(PCHPDDM,PCCreate_HPDDM);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode SlepcCitationsInitialize()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscCitationsRegister("@Article{slepc-toms,\n"
    "   author = \"Vicente Hernandez and Jose E. Roman and Vicente Vidal\",\n"
    "   title = \"{SLEPc}: A Scalable and Flexible Toolkit for the Solution of Eigenvalue Problems\",\n"
    "   journal = \"{ACM} Trans. Math. Software\",\n"
    "   volume = \"31\",\n"
    "   number = \"3\",\n"
    "   pages = \"351--362\",\n"
    "   year = \"2005,\"\n"
    "   doi = \"https://doi.org/10.1145/1089014.1089019\"\n"
    "}\n",NULL);CHKERRQ(ierr);
  ierr = PetscCitationsRegister("@TechReport{slepc-manual,\n"
    "   author = \"J. E. Roman and C. Campos and E. Romero and A. Tomas\",\n"
    "   title = \"{SLEPc} Users Manual\",\n"
    "   number = \"DSIC-II/24/02 - Revision 3.13\",\n"
    "   institution = \"D. Sistemes Inform\\`atics i Computaci\\'o, Universitat Polit\\`ecnica de Val\\`encia\",\n"
    "   year = \"2020\"\n"
    "}\n",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   SlepcInitialize - Initializes the SLEPc library. SlepcInitialize() calls
   PetscInitialize() if that has not been called yet, so this routine should
   always be called near the beginning of your program.

   Collective on MPI_COMM_WORLD or PETSC_COMM_WORLD if it has been set

   Input Parameters:
+  argc - count of number of command line arguments
.  args - the command line arguments
.  file - [optional] PETSc database file, defaults to ~username/.petscrc
          (use NULL for default)
-  help - [optional] Help message to print, use NULL for no message

   Fortran Note:
   Fortran syntax is very similar to that of PetscInitialize()

   Level: beginner

.seealso: SlepcFinalize(), PetscInitialize()
@*/
PetscErrorCode SlepcInitialize(int *argc,char ***args,const char file[],const char help[])
{
  PetscErrorCode ierr;
  PetscBool      flg;

  PetscFunctionBegin;
  if (SlepcInitializeCalled) PetscFunctionReturn(0);
  ierr = PetscSetHelpVersionFunctions(SlepcPrintHelpIntro,SlepcPrintVersion);CHKERRQ(ierr);
  ierr = PetscInitialized(&flg);CHKERRQ(ierr);
  if (!flg) {
    ierr = PetscInitialize(argc,args,file,help);CHKERRQ(ierr);
    SlepcBeganPetsc = PETSC_TRUE;
  }

  ierr = SlepcCitationsInitialize();CHKERRQ(ierr);

  /* Load the dynamic libraries (on machines that support them), this registers all the solvers etc. */
  ierr = SlepcInitialize_DynamicLibraries();CHKERRQ(ierr);

#if defined(PETSC_HAVE_DRAND48)
  /* work-around for Cygwin drand48() initialization bug */
  srand48(0);
#endif

  SlepcInitializeCalled = PETSC_TRUE;
  ierr = PetscInfo(0,"SLEPc successfully started\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   SlepcFinalize - Checks for options to be called at the conclusion
   of the SLEPc program and calls PetscFinalize().

   Collective on PETSC_COMM_WORLD

   Level: beginner

.seealso: SlepcInitialize(), PetscFinalize()
@*/
PetscErrorCode SlepcFinalize(void)
{
  PetscErrorCode ierr = 0;

  PetscFunctionBegin;
  ierr = PetscInfo(0,"SlepcFinalize() called\n");CHKERRQ(ierr);
  if (SlepcBeganPetsc) {
    ierr = PetscFinalize();
  }
  SlepcInitializeCalled = PETSC_FALSE;
  PetscFunctionReturn(ierr);
}

/*@C
   SlepcInitializeNoArguments - Calls SlepcInitialize() from C/C++ without
   the command line arguments.

   Collective

   Level: advanced

.seealso: SlepcInitialize(), SlepcInitializeFortran()
@*/
PetscErrorCode SlepcInitializeNoArguments(void)
{
  PetscErrorCode ierr;
  int            argc = 0;
  char           **args = 0;

  PetscFunctionBegin;
  ierr = SlepcInitialize(&argc,&args,NULL,NULL);
  PetscFunctionReturn(ierr);
}

/*@
   SlepcInitialized - Determine whether SLEPc is initialized.

   Level: beginner

.seealso: SlepcInitialize(), SlepcInitializeFortran()
@*/
PetscErrorCode SlepcInitialized(PetscBool *isInitialized)
{
  PetscFunctionBegin;
  PetscValidBoolPointer(isInitialized,1);
  *isInitialized = SlepcInitializeCalled;
  PetscFunctionReturn(0);
}

PETSC_EXTERN PetscBool PetscBeganMPI;

/*
   SlepcInitializeNoPointers - Calls SlepcInitialize() from C/C++ without the pointers
   to argc and args (analogue to PetscInitializeNoPointers).

   Collective

   Level: advanced

.seealso: SlepcInitialize()
*/
PetscErrorCode SlepcInitializeNoPointers(int argc,char **args,const char *filename,const char *help)
{
  PetscErrorCode ierr;
  int            myargc = argc;
  char           **myargs = args;

  PetscFunctionBegin;
  ierr = SlepcInitialize(&myargc,&myargs,filename,help);CHKERRQ(ierr);
  ierr = PetscPopSignalHandler();CHKERRQ(ierr);
  PetscBeganMPI = PETSC_FALSE;
  PetscFunctionReturn(ierr);
}

