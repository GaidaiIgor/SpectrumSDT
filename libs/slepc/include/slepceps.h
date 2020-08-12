/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the SLEPc linear eigenvalue solvers
*/

#if !defined(SLEPCEPS_H)
#define SLEPCEPS_H
#include <slepcst.h>
#include <slepcbv.h>
#include <slepcds.h>
#include <slepcrg.h>
#include <slepclme.h>
#include <petscsnes.h>

SLEPC_EXTERN PetscErrorCode EPSInitializePackage(void);

/*S
    EPS - Abstract SLEPc object that manages all the eigenvalue
    problem solvers.

    Level: beginner

.seealso:  EPSCreate(), ST
S*/
typedef struct _p_EPS* EPS;

/*J
    EPSType - String with the name of a SLEPc eigensolver

    Level: beginner

.seealso: EPSSetType(), EPS
J*/
typedef const char* EPSType;
#define EPSPOWER       "power"
#define EPSSUBSPACE    "subspace"
#define EPSARNOLDI     "arnoldi"
#define EPSLANCZOS     "lanczos"
#define EPSKRYLOVSCHUR "krylovschur"
#define EPSGD          "gd"
#define EPSJD          "jd"
#define EPSRQCG        "rqcg"
#define EPSLOBPCG      "lobpcg"
#define EPSCISS        "ciss"
#define EPSLYAPII      "lyapii"
#define EPSLAPACK      "lapack"
#define EPSARPACK      "arpack"
#define EPSBLZPACK     "blzpack"
#define EPSTRLAN       "trlan"
#define EPSBLOPEX      "blopex"
#define EPSPRIMME      "primme"

/* Logging support */
SLEPC_EXTERN PetscClassId EPS_CLASSID;

/*E
    EPSProblemType - Determines the type of eigenvalue problem

    Level: beginner

.seealso: EPSSetProblemType(), EPSGetProblemType()
E*/
typedef enum { EPS_HEP=1,
               EPS_GHEP,
               EPS_NHEP,
               EPS_GNHEP,
               EPS_PGNHEP,
               EPS_GHIEP } EPSProblemType;

/*E
    EPSExtraction - Determines the type of extraction technique employed
    by the eigensolver

    Level: advanced

.seealso: EPSSetExtraction(), EPSGetExtraction()
E*/
typedef enum { EPS_RITZ,
               EPS_HARMONIC,
               EPS_HARMONIC_RELATIVE,
               EPS_HARMONIC_RIGHT,
               EPS_HARMONIC_LARGEST,
               EPS_REFINED,
               EPS_REFINED_HARMONIC } EPSExtraction;

/*E
    EPSWhich - Determines which part of the spectrum is requested

    Level: intermediate

.seealso: EPSSetWhichEigenpairs(), EPSGetWhichEigenpairs()
E*/
typedef enum { EPS_LARGEST_MAGNITUDE=1,
               EPS_SMALLEST_MAGNITUDE,
               EPS_LARGEST_REAL,
               EPS_SMALLEST_REAL,
               EPS_LARGEST_IMAGINARY,
               EPS_SMALLEST_IMAGINARY,
               EPS_TARGET_MAGNITUDE,
               EPS_TARGET_REAL,
               EPS_TARGET_IMAGINARY,
               EPS_ALL,
               EPS_WHICH_USER } EPSWhich;

/*E
    EPSBalance - The type of balancing used for non-Hermitian problems

    Level: intermediate

.seealso: EPSSetBalance()
E*/
typedef enum { EPS_BALANCE_NONE,
               EPS_BALANCE_ONESIDE,
               EPS_BALANCE_TWOSIDE,
               EPS_BALANCE_USER } EPSBalance;
SLEPC_EXTERN const char *EPSBalanceTypes[];

/*E
    EPSErrorType - The error type used to assess accuracy of computed solutions

    Level: intermediate

.seealso: EPSComputeError()
E*/
typedef enum { EPS_ERROR_ABSOLUTE,
               EPS_ERROR_RELATIVE,
               EPS_ERROR_BACKWARD } EPSErrorType;
SLEPC_EXTERN const char *EPSErrorTypes[];

/*E
    EPSConv - Determines the convergence test

    Level: intermediate

.seealso: EPSSetConvergenceTest(), EPSSetConvergenceTestFunction()
E*/
typedef enum { EPS_CONV_ABS,
               EPS_CONV_REL,
               EPS_CONV_NORM,
               EPS_CONV_USER } EPSConv;

/*E
    EPSStop - Determines the stopping test

    Level: advanced

.seealso: EPSSetStoppingTest(), EPSSetStoppingTestFunction()
E*/
typedef enum { EPS_STOP_BASIC,
               EPS_STOP_USER } EPSStop;

/*E
    EPSConvergedReason - Reason an eigensolver was said to
         have converged or diverged

    Level: intermediate

.seealso: EPSSolve(), EPSGetConvergedReason(), EPSSetTolerances()
E*/
typedef enum {/* converged */
              EPS_CONVERGED_TOL                =  1,
              EPS_CONVERGED_USER               =  2,
              /* diverged */
              EPS_DIVERGED_ITS                 = -1,
              EPS_DIVERGED_BREAKDOWN           = -2,
              EPS_DIVERGED_SYMMETRY_LOST       = -3,
              EPS_CONVERGED_ITERATING          =  0} EPSConvergedReason;
SLEPC_EXTERN const char *const*EPSConvergedReasons;

SLEPC_EXTERN PetscErrorCode EPSCreate(MPI_Comm,EPS*);
SLEPC_EXTERN PetscErrorCode EPSDestroy(EPS*);
SLEPC_EXTERN PetscErrorCode EPSReset(EPS);
SLEPC_EXTERN PetscErrorCode EPSSetType(EPS,EPSType);
SLEPC_EXTERN PetscErrorCode EPSGetType(EPS,EPSType*);
SLEPC_EXTERN PetscErrorCode EPSSetProblemType(EPS,EPSProblemType);
SLEPC_EXTERN PetscErrorCode EPSGetProblemType(EPS,EPSProblemType*);
SLEPC_EXTERN PetscErrorCode EPSSetExtraction(EPS,EPSExtraction);
SLEPC_EXTERN PetscErrorCode EPSGetExtraction(EPS,EPSExtraction*);
SLEPC_EXTERN PetscErrorCode EPSSetBalance(EPS,EPSBalance,PetscInt,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSGetBalance(EPS,EPSBalance*,PetscInt*,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSSetOperators(EPS,Mat,Mat);
SLEPC_EXTERN PetscErrorCode EPSGetOperators(EPS,Mat*,Mat*);
SLEPC_EXTERN PetscErrorCode EPSSetFromOptions(EPS);
SLEPC_EXTERN PetscErrorCode EPSSetUp(EPS);
SLEPC_EXTERN PetscErrorCode EPSSolve(EPS);
SLEPC_EXTERN PetscErrorCode EPSView(EPS,PetscViewer);
SLEPC_EXTERN PetscErrorCode EPSViewFromOptions(EPS,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode EPSErrorView(EPS,EPSErrorType,PetscViewer);
PETSC_DEPRECATED_FUNCTION("Use EPSErrorView()") PETSC_STATIC_INLINE PetscErrorCode EPSPrintSolution(EPS eps,PetscViewer v) {return EPSErrorView(eps,EPS_ERROR_RELATIVE,v);}
SLEPC_EXTERN PetscErrorCode EPSErrorViewFromOptions(EPS);
SLEPC_EXTERN PetscErrorCode EPSReasonView(EPS,PetscViewer);
SLEPC_EXTERN PetscErrorCode EPSReasonViewFromOptions(EPS);
SLEPC_EXTERN PetscErrorCode EPSValuesView(EPS,PetscViewer);
SLEPC_EXTERN PetscErrorCode EPSValuesViewFromOptions(EPS);
SLEPC_EXTERN PetscErrorCode EPSVectorsView(EPS,PetscViewer);
SLEPC_EXTERN PetscErrorCode EPSVectorsViewFromOptions(EPS);

SLEPC_EXTERN PetscErrorCode EPSSetTarget(EPS,PetscScalar);
SLEPC_EXTERN PetscErrorCode EPSGetTarget(EPS,PetscScalar*);
SLEPC_EXTERN PetscErrorCode EPSSetInterval(EPS,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSGetInterval(EPS,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSSetST(EPS,ST);
SLEPC_EXTERN PetscErrorCode EPSGetST(EPS,ST*);
SLEPC_EXTERN PetscErrorCode EPSSetBV(EPS,BV);
SLEPC_EXTERN PetscErrorCode EPSGetBV(EPS,BV*);
SLEPC_EXTERN PetscErrorCode EPSSetRG(EPS,RG);
SLEPC_EXTERN PetscErrorCode EPSGetRG(EPS,RG*);
SLEPC_EXTERN PetscErrorCode EPSSetDS(EPS,DS);
SLEPC_EXTERN PetscErrorCode EPSGetDS(EPS,DS*);
SLEPC_EXTERN PetscErrorCode EPSSetTolerances(EPS,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGetTolerances(EPS,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSSetConvergenceTestFunction(EPS,PetscErrorCode (*)(EPS,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*),void*,PetscErrorCode (*)(void*));
SLEPC_EXTERN PetscErrorCode EPSSetConvergenceTest(EPS,EPSConv);
SLEPC_EXTERN PetscErrorCode EPSGetConvergenceTest(EPS,EPSConv*);
SLEPC_EXTERN PetscErrorCode EPSConvergedAbsolute(EPS,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode EPSConvergedRelative(EPS,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode EPSConvergedNorm(EPS,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode EPSSetStoppingTestFunction(EPS,PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscInt,PetscInt,EPSConvergedReason*,void*),void*,PetscErrorCode (*)(void*));
SLEPC_EXTERN PetscErrorCode EPSSetStoppingTest(EPS,EPSStop);
SLEPC_EXTERN PetscErrorCode EPSGetStoppingTest(EPS,EPSStop*);
SLEPC_EXTERN PetscErrorCode EPSStoppingBasic(EPS,PetscInt,PetscInt,PetscInt,PetscInt,EPSConvergedReason*,void*);
SLEPC_EXTERN PetscErrorCode EPSGetConvergedReason(EPS,EPSConvergedReason*);

SLEPC_EXTERN PetscErrorCode EPSSetDimensions(EPS,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGetDimensions(EPS,PetscInt*,PetscInt*,PetscInt*);

SLEPC_EXTERN PetscErrorCode EPSGetConverged(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSGetEigenpair(EPS,PetscInt,PetscScalar*,PetscScalar*,Vec,Vec);
SLEPC_EXTERN PetscErrorCode EPSGetEigenvalue(EPS,PetscInt,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode EPSGetEigenvector(EPS,PetscInt,Vec,Vec);
SLEPC_EXTERN PetscErrorCode EPSGetLeftEigenvector(EPS,PetscInt,Vec,Vec);

SLEPC_EXTERN PetscErrorCode EPSComputeError(EPS,PetscInt,EPSErrorType,PetscReal*);
PETSC_DEPRECATED_FUNCTION("Use EPSComputeError()") PETSC_STATIC_INLINE PetscErrorCode EPSComputeRelativeError(EPS eps,PetscInt i,PetscReal *r) {return EPSComputeError(eps,i,EPS_ERROR_RELATIVE,r);}
PETSC_DEPRECATED_FUNCTION("Use EPSComputeError() with EPS_ERROR_ABSOLUTE") PETSC_STATIC_INLINE PetscErrorCode EPSComputeResidualNorm(EPS eps,PetscInt i,PetscReal *r) {return EPSComputeError(eps,i,EPS_ERROR_ABSOLUTE,r);}
SLEPC_EXTERN PetscErrorCode EPSGetInvariantSubspace(EPS,Vec[]);
SLEPC_EXTERN PetscErrorCode EPSGetErrorEstimate(EPS,PetscInt,PetscReal*);

SLEPC_EXTERN PetscErrorCode EPSMonitor(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSMonitorSet(EPS,PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*),void*,PetscErrorCode (*)(void**));
SLEPC_EXTERN PetscErrorCode EPSMonitorSetFromOptions(EPS,const char*,const char*,const char*,PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscBool);
SLEPC_EXTERN PetscErrorCode EPSConvMonitorSetFromOptions(EPS,const char*,const char*,const char*,PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,SlepcConvMonitor));
SLEPC_EXTERN PetscErrorCode EPSMonitorCancel(EPS);
SLEPC_EXTERN PetscErrorCode EPSGetMonitorContext(EPS,void**);
SLEPC_EXTERN PetscErrorCode EPSGetIterationNumber(EPS,PetscInt*);

SLEPC_EXTERN PetscErrorCode EPSSetWhichEigenpairs(EPS,EPSWhich);
SLEPC_EXTERN PetscErrorCode EPSGetWhichEigenpairs(EPS,EPSWhich*);
SLEPC_EXTERN PetscErrorCode EPSSetTwoSided(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGetTwoSided(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSSetTrueResidual(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGetTrueResidual(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSSetPurify(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGetPurify(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSSetEigenvalueComparison(EPS,PetscErrorCode (*func)(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*),void*);
SLEPC_EXTERN PetscErrorCode EPSSetArbitrarySelection(EPS,PetscErrorCode (*func)(PetscScalar,PetscScalar,Vec,Vec,PetscScalar*,PetscScalar*,void*),void*);
SLEPC_EXTERN PetscErrorCode EPSIsGeneralized(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSIsHermitian(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSIsPositive(EPS,PetscBool*);

SLEPC_EXTERN PetscErrorCode EPSMonitorFirst(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode EPSMonitorAll(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode EPSMonitorConverged(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,SlepcConvMonitor);
SLEPC_EXTERN PetscErrorCode EPSMonitorLGCreate(MPI_Comm,const char[],const char[],int,int,int,int,PetscDrawLG*);
SLEPC_EXTERN PetscErrorCode EPSMonitorLG(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
SLEPC_EXTERN PetscErrorCode EPSMonitorLGAll(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);

SLEPC_EXTERN PetscErrorCode EPSSetTrackAll(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGetTrackAll(EPS,PetscBool*);

SLEPC_EXTERN PetscErrorCode EPSSetDeflationSpace(EPS,PetscInt,Vec[]);
SLEPC_EXTERN PetscErrorCode EPSSetInitialSpace(EPS,PetscInt,Vec[]);
SLEPC_EXTERN PetscErrorCode EPSSetLeftInitialSpace(EPS,PetscInt,Vec[]);

SLEPC_EXTERN PetscErrorCode EPSSetOptionsPrefix(EPS,const char*);
SLEPC_EXTERN PetscErrorCode EPSAppendOptionsPrefix(EPS,const char*);
SLEPC_EXTERN PetscErrorCode EPSGetOptionsPrefix(EPS,const char*[]);

SLEPC_EXTERN PetscFunctionList EPSList;
SLEPC_EXTERN PetscErrorCode EPSRegister(const char[],PetscErrorCode(*)(EPS));

SLEPC_EXTERN PetscErrorCode EPSSetWorkVecs(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSAllocateSolution(EPS,PetscInt);

/* --------- options specific to particular eigensolvers -------- */

/*E
    EPSPowerShiftType - determines the type of shift used in the Power iteration

    Level: advanced

.seealso: EPSPowerSetShiftType(), EPSPowerGetShiftType()
E*/
typedef enum { EPS_POWER_SHIFT_CONSTANT,
               EPS_POWER_SHIFT_RAYLEIGH,
               EPS_POWER_SHIFT_WILKINSON } EPSPowerShiftType;
SLEPC_EXTERN const char *EPSPowerShiftTypes[];

SLEPC_EXTERN PetscErrorCode EPSPowerSetShiftType(EPS,EPSPowerShiftType);
SLEPC_EXTERN PetscErrorCode EPSPowerGetShiftType(EPS,EPSPowerShiftType*);
SLEPC_EXTERN PetscErrorCode EPSPowerSetNonlinear(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSPowerGetNonlinear(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSPowerSetUpdate(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSPowerGetUpdate(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSPowerSetSNES(EPS,SNES);
SLEPC_EXTERN PetscErrorCode EPSPowerGetSNES(EPS,SNES*);

SLEPC_EXTERN PetscErrorCode EPSArnoldiSetDelayed(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSArnoldiGetDelayed(EPS,PetscBool*);

SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetRestart(EPS,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetRestart(EPS,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetLocking(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetLocking(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetPartitions(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetPartitions(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetDetectZeros(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetDetectZeros(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetDimensions(EPS,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetDimensions(EPS,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurSetSubintervals(EPS,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetSubintervals(EPS,PetscReal**);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetInertias(EPS,PetscInt*,PetscReal**,PetscInt**);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetSubcommInfo(EPS,PetscInt*,PetscInt*,Vec*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetSubcommPairs(EPS,PetscInt,PetscScalar*,Vec);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurGetSubcommMats(EPS,Mat*,Mat*);
SLEPC_EXTERN PetscErrorCode EPSKrylovSchurUpdateSubcommMats(EPS,PetscScalar,PetscScalar,Mat,PetscScalar,PetscScalar, Mat,MatStructure,PetscBool);

/*E
    EPSLanczosReorthogType - determines the type of reorthogonalization
    used in the Lanczos method

    Level: advanced

.seealso: EPSLanczosSetReorthog(), EPSLanczosGetReorthog()
E*/
typedef enum { EPS_LANCZOS_REORTHOG_LOCAL,
               EPS_LANCZOS_REORTHOG_FULL,
               EPS_LANCZOS_REORTHOG_SELECTIVE,
               EPS_LANCZOS_REORTHOG_PERIODIC,
               EPS_LANCZOS_REORTHOG_PARTIAL,
               EPS_LANCZOS_REORTHOG_DELAYED } EPSLanczosReorthogType;
SLEPC_EXTERN const char *EPSLanczosReorthogTypes[];

SLEPC_EXTERN PetscErrorCode EPSLanczosSetReorthog(EPS,EPSLanczosReorthogType);
SLEPC_EXTERN PetscErrorCode EPSLanczosGetReorthog(EPS,EPSLanczosReorthogType*);

SLEPC_EXTERN PetscErrorCode EPSBlzpackSetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSBlzpackGetBlockSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSBlzpackSetNSteps(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSBlzpackGetNSteps(EPS,PetscInt*);

/*E
    EPSPRIMMEMethod - determines the method selected in the PRIMME library

    Level: advanced

.seealso: EPSPRIMMESetMethod(), EPSPRIMMEGetMethod()
E*/
typedef enum { EPS_PRIMME_DYNAMIC=1,
               EPS_PRIMME_DEFAULT_MIN_TIME,
               EPS_PRIMME_DEFAULT_MIN_MATVECS,
               EPS_PRIMME_ARNOLDI,
               EPS_PRIMME_GD,
               EPS_PRIMME_GD_PLUSK,
               EPS_PRIMME_GD_OLSEN_PLUSK,
               EPS_PRIMME_JD_OLSEN_PLUSK,
               EPS_PRIMME_RQI,
               EPS_PRIMME_JDQR,
               EPS_PRIMME_JDQMR,
               EPS_PRIMME_JDQMR_ETOL,
               EPS_PRIMME_SUBSPACE_ITERATION,
               EPS_PRIMME_LOBPCG_ORTHOBASIS,
               EPS_PRIMME_LOBPCG_ORTHOBASISW } EPSPRIMMEMethod;
SLEPC_EXTERN const char *EPSPRIMMEMethods[];

SLEPC_EXTERN PetscErrorCode EPSPRIMMESetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSPRIMMEGetBlockSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSPRIMMESetMethod(EPS,EPSPRIMMEMethod);
SLEPC_EXTERN PetscErrorCode EPSPRIMMEGetMethod(EPS,EPSPRIMMEMethod*);

SLEPC_EXTERN PetscErrorCode EPSGDSetKrylovStart(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGDGetKrylovStart(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSGDSetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGDGetBlockSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSGDSetRestart(EPS,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGDGetRestart(EPS,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSGDSetInitialSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSGDGetInitialSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSGDSetBOrth(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGDGetBOrth(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSGDSetDoubleExpansion(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSGDGetDoubleExpansion(EPS,PetscBool*);

SLEPC_EXTERN PetscErrorCode EPSJDSetKrylovStart(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSJDGetKrylovStart(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSJDSetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSJDGetBlockSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSJDSetRestart(EPS,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSJDGetRestart(EPS,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSJDSetInitialSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSJDGetInitialSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSJDSetFix(EPS,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSJDGetFix(EPS,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSJDSetConstCorrectionTol(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSJDGetConstCorrectionTol(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSJDSetBOrth(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSJDGetBOrth(EPS,PetscBool*);

SLEPC_EXTERN PetscErrorCode EPSRQCGSetReset(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSRQCGGetReset(EPS,PetscInt*);

SLEPC_EXTERN PetscErrorCode EPSLOBPCGSetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGGetBlockSize(EPS,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGSetRestart(EPS,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGGetRestart(EPS,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGSetLocking(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSLOBPCGGetLocking(EPS,PetscBool*);

/*E
    EPSCISSQuadRule - determines the quadrature rule in the CISS solver

    Level: advanced

.seealso: EPSCISSSetQuadRule(), EPSCISSGetQuadRule()
E*/
typedef enum { EPS_CISS_QUADRULE_TRAPEZOIDAL=1,
               EPS_CISS_QUADRULE_CHEBYSHEV } EPSCISSQuadRule;
SLEPC_EXTERN const char *EPSCISSQuadRules[];

/*E
    EPSCISSExtraction - determines the extraction technique in the CISS solver

    Level: advanced

.seealso: EPSCISSSetExtraction(), EPSCISSGetExtraction()
E*/
typedef enum { EPS_CISS_EXTRACTION_RITZ,
               EPS_CISS_EXTRACTION_HANKEL } EPSCISSExtraction;
SLEPC_EXTERN const char *EPSCISSExtractions[];

SLEPC_EXTERN PetscErrorCode EPSCISSSetExtraction(EPS,EPSCISSExtraction);
SLEPC_EXTERN PetscErrorCode EPSCISSGetExtraction(EPS,EPSCISSExtraction*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetQuadRule(EPS,EPSCISSQuadRule);
SLEPC_EXTERN PetscErrorCode EPSCISSGetQuadRule(EPS,EPSCISSQuadRule*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetSizes(EPS,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSCISSGetSizes(EPS,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetThreshold(EPS,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode EPSCISSGetThreshold(EPS,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetRefinement(EPS,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSCISSGetRefinement(EPS,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode EPSCISSSetUseST(EPS,PetscBool);
SLEPC_EXTERN PetscErrorCode EPSCISSGetUseST(EPS,PetscBool*);
SLEPC_EXTERN PetscErrorCode EPSCISSGetKSPs(EPS,PetscInt*,KSP**);

SLEPC_EXTERN PetscErrorCode EPSLyapIISetLME(EPS,LME);
SLEPC_EXTERN PetscErrorCode EPSLyapIIGetLME(EPS,LME*);
SLEPC_EXTERN PetscErrorCode EPSLyapIISetRanks(EPS,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSLyapIIGetRanks(EPS,PetscInt*,PetscInt*);

SLEPC_EXTERN PetscErrorCode EPSBLOPEXSetBlockSize(EPS,PetscInt);
SLEPC_EXTERN PetscErrorCode EPSBLOPEXGetBlockSize(EPS,PetscInt*);

#endif

