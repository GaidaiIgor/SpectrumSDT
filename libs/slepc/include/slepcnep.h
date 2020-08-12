/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for SLEPc's nonlinear eigenvalue solvers
*/

#if !defined(SLEPCNEP_H)
#define SLEPCNEP_H
#include <slepceps.h>
#include <slepcpep.h>
#include <slepcfn.h>

SLEPC_EXTERN PetscErrorCode NEPInitializePackage(void);

/*S
     NEP - Abstract SLEPc object that manages all solvers for
     nonlinear eigenvalue problems.

   Level: beginner

.seealso:  NEPCreate()
S*/
typedef struct _p_NEP* NEP;

/*J
    NEPType - String with the name of a nonlinear eigensolver

   Level: beginner

.seealso: NEPSetType(), NEP
J*/
typedef const char* NEPType;
#define NEPRII       "rii"
#define NEPSLP       "slp"
#define NEPNARNOLDI  "narnoldi"
#define NEPCISS      "ciss"
#define NEPINTERPOL  "interpol"
#define NEPNLEIGS    "nleigs"

/* Logging support */
SLEPC_EXTERN PetscClassId NEP_CLASSID;

/*E
    NEPProblemType - Determines the type of the nonlinear eigenproblem

    Level: intermediate

.seealso: NEPSetProblemType(), NEPGetProblemType()
E*/
typedef enum { NEP_GENERAL=1,
               NEP_RATIONAL     /* NEP defined in split form with all f_i rational */
             } NEPProblemType;

/*E
    NEPWhich - Determines which part of the spectrum is requested

    Level: intermediate

.seealso: NEPSetWhichEigenpairs(), NEPGetWhichEigenpairs()
E*/
typedef enum { NEP_LARGEST_MAGNITUDE=1,
               NEP_SMALLEST_MAGNITUDE,
               NEP_LARGEST_REAL,
               NEP_SMALLEST_REAL,
               NEP_LARGEST_IMAGINARY,
               NEP_SMALLEST_IMAGINARY,
               NEP_TARGET_MAGNITUDE,
               NEP_TARGET_REAL,
               NEP_TARGET_IMAGINARY,
               NEP_ALL,
               NEP_WHICH_USER } NEPWhich;

/*E
    NEPErrorType - The error type used to assess accuracy of computed solutions

    Level: intermediate

.seealso: NEPComputeError()
E*/
typedef enum { NEP_ERROR_ABSOLUTE,
               NEP_ERROR_RELATIVE,
               NEP_ERROR_BACKWARD } NEPErrorType;
SLEPC_EXTERN const char *NEPErrorTypes[];

/*E
    NEPRefine - The refinement type

    Level: intermediate

.seealso: NEPSetRefine()
E*/
typedef enum { NEP_REFINE_NONE,
               NEP_REFINE_SIMPLE,
               NEP_REFINE_MULTIPLE } NEPRefine;
SLEPC_EXTERN const char *NEPRefineTypes[];

/*E
    NEPRefineScheme - The scheme used for solving linear systems during iterative refinement

    Level: intermediate

.seealso: NEPSetRefine()
E*/
typedef enum { NEP_REFINE_SCHEME_SCHUR=1,
               NEP_REFINE_SCHEME_MBE,
               NEP_REFINE_SCHEME_EXPLICIT } NEPRefineScheme;
SLEPC_EXTERN const char *NEPRefineSchemes[];

/*E
    NEPConv - Determines the convergence test

    Level: intermediate

.seealso: NEPSetConvergenceTest(), NEPSetConvergenceTestFunction()
E*/
typedef enum { NEP_CONV_ABS,
               NEP_CONV_REL,
               NEP_CONV_NORM,
               NEP_CONV_USER } NEPConv;

/*E
    NEPStop - Determines the stopping test

    Level: advanced

.seealso: NEPSetStoppingTest(), NEPSetStoppingTestFunction()
E*/
typedef enum { NEP_STOP_BASIC,
               NEP_STOP_USER } NEPStop;

/*E
    NEPConvergedReason - Reason a nonlinear eigensolver was said to
         have converged or diverged

    Level: intermediate

.seealso: NEPSolve(), NEPGetConvergedReason(), NEPSetTolerances()
E*/
typedef enum {/* converged */
              NEP_CONVERGED_TOL                =  1,
              NEP_CONVERGED_USER               =  2,
              /* diverged */
              NEP_DIVERGED_ITS                 = -1,
              NEP_DIVERGED_BREAKDOWN           = -2,
                    /* unused                  = -3 */
              NEP_DIVERGED_LINEAR_SOLVE        = -4,
              NEP_DIVERGED_SUBSPACE_EXHAUSTED  = -5,
              NEP_CONVERGED_ITERATING          =  0} NEPConvergedReason;
SLEPC_EXTERN const char *const*NEPConvergedReasons;

SLEPC_EXTERN PetscErrorCode NEPCreate(MPI_Comm,NEP*);
SLEPC_EXTERN PetscErrorCode NEPDestroy(NEP*);
SLEPC_EXTERN PetscErrorCode NEPReset(NEP);
SLEPC_EXTERN PetscErrorCode NEPSetType(NEP,NEPType);
SLEPC_EXTERN PetscErrorCode NEPGetType(NEP,NEPType*);
SLEPC_EXTERN PetscErrorCode NEPSetProblemType(NEP,NEPProblemType);
SLEPC_EXTERN PetscErrorCode NEPGetProblemType(NEP,NEPProblemType*);
SLEPC_EXTERN PetscErrorCode NEPSetTarget(NEP,PetscScalar);
SLEPC_EXTERN PetscErrorCode NEPGetTarget(NEP,PetscScalar*);
SLEPC_EXTERN PetscErrorCode NEPSetFromOptions(NEP);
SLEPC_EXTERN PetscErrorCode NEPSetUp(NEP);
SLEPC_EXTERN PetscErrorCode NEPSolve(NEP);
SLEPC_EXTERN PetscErrorCode NEPView(NEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPViewFromOptions(NEP,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode NEPErrorView(NEP,NEPErrorType,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPErrorViewFromOptions(NEP);
SLEPC_EXTERN PetscErrorCode NEPReasonView(NEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPReasonViewFromOptions(NEP);
SLEPC_EXTERN PetscErrorCode NEPValuesView(NEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPValuesViewFromOptions(NEP);
SLEPC_EXTERN PetscErrorCode NEPVectorsView(NEP,PetscViewer);
SLEPC_EXTERN PetscErrorCode NEPVectorsViewFromOptions(NEP);

SLEPC_EXTERN PetscErrorCode NEPSetFunction(NEP,Mat,Mat,PetscErrorCode (*)(NEP,PetscScalar,Mat,Mat,void*),void*);
SLEPC_EXTERN PetscErrorCode NEPGetFunction(NEP,Mat*,Mat*,PetscErrorCode (**)(NEP,PetscScalar,Mat,Mat,void*),void**);
SLEPC_EXTERN PetscErrorCode NEPSetJacobian(NEP,Mat,PetscErrorCode (*)(NEP,PetscScalar,Mat,void*),void*);
SLEPC_EXTERN PetscErrorCode NEPGetJacobian(NEP,Mat*,PetscErrorCode (**)(NEP,PetscScalar,Mat,void*),void**);
PETSC_DEPRECATED_FUNCTION("Use NEPSetFunction() and NEPSetJacobian") PETSC_STATIC_INLINE PetscErrorCode NEPSetDerivatives(NEP nep,Mat A,PetscErrorCode (*fun)(NEP,PetscScalar,PetscInt,Mat,void*),void *ctx) {(void)A;(void)fun;(void)ctx;SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Not implemented in this version");}
PETSC_DEPRECATED_FUNCTION("Use NEPGetFunction() and NEPGetJacobian") PETSC_STATIC_INLINE PetscErrorCode NEPGetDerivatives(NEP nep,Mat *A,PetscErrorCode (**fun)(NEP,PetscScalar,PetscInt,Mat,void*),void **ctx) {(void)A;(void)fun;(void)ctx;SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Not implemented in this version");}
SLEPC_EXTERN PetscErrorCode NEPSetSplitOperator(NEP,PetscInt,Mat[],FN[],MatStructure);
SLEPC_EXTERN PetscErrorCode NEPGetSplitOperatorTerm(NEP,PetscInt,Mat*,FN*);
SLEPC_EXTERN PetscErrorCode NEPGetSplitOperatorInfo(NEP,PetscInt*,MatStructure*);

SLEPC_EXTERN PetscErrorCode NEPSetBV(NEP,BV);
SLEPC_EXTERN PetscErrorCode NEPGetBV(NEP,BV*);
SLEPC_EXTERN PetscErrorCode NEPSetRG(NEP,RG);
SLEPC_EXTERN PetscErrorCode NEPGetRG(NEP,RG*);
SLEPC_EXTERN PetscErrorCode NEPSetDS(NEP,DS);
SLEPC_EXTERN PetscErrorCode NEPGetDS(NEP,DS*);
SLEPC_EXTERN PetscErrorCode NEPRefineGetKSP(NEP,KSP*);
SLEPC_EXTERN PetscErrorCode NEPSetTolerances(NEP,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPGetTolerances(NEP,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPSetConvergenceTestFunction(NEP,PetscErrorCode (*)(NEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*),void*,PetscErrorCode (*)(void*));
SLEPC_EXTERN PetscErrorCode NEPSetConvergenceTest(NEP,NEPConv);
SLEPC_EXTERN PetscErrorCode NEPGetConvergenceTest(NEP,NEPConv*);
SLEPC_EXTERN PetscErrorCode NEPConvergedAbsolute(NEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode NEPConvergedRelative(NEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode NEPConvergedNorm(NEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode NEPSetStoppingTestFunction(NEP,PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscInt,PetscInt,NEPConvergedReason*,void*),void*,PetscErrorCode (*)(void*));
SLEPC_EXTERN PetscErrorCode NEPSetStoppingTest(NEP,NEPStop);
SLEPC_EXTERN PetscErrorCode NEPGetStoppingTest(NEP,NEPStop*);
SLEPC_EXTERN PetscErrorCode NEPStoppingBasic(NEP,PetscInt,PetscInt,PetscInt,PetscInt,NEPConvergedReason*,void*);
SLEPC_EXTERN PetscErrorCode NEPSetDimensions(NEP,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPGetDimensions(NEP,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPSetRefine(NEP,NEPRefine,PetscInt,PetscReal,PetscInt,NEPRefineScheme);
SLEPC_EXTERN PetscErrorCode NEPGetRefine(NEP,NEPRefine*,PetscInt*,PetscReal*,PetscInt*,NEPRefineScheme*);

SLEPC_EXTERN PetscErrorCode NEPGetConverged(NEP,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPGetEigenpair(NEP,PetscInt,PetscScalar*,PetscScalar*,Vec,Vec);
SLEPC_EXTERN PetscErrorCode NEPGetLeftEigenvector(NEP,PetscInt,Vec,Vec);

SLEPC_EXTERN PetscErrorCode NEPComputeError(NEP,PetscInt,NEPErrorType,PetscReal*);
PETSC_DEPRECATED_FUNCTION("Use NEPComputeError()") PETSC_STATIC_INLINE PetscErrorCode NEPComputeRelativeError(NEP nep,PetscInt i,PetscReal *r) {return NEPComputeError(nep,i,NEP_ERROR_RELATIVE,r);}
PETSC_DEPRECATED_FUNCTION("Use NEPComputeError() with NEP_ERROR_ABSOLUTE") PETSC_STATIC_INLINE PetscErrorCode NEPComputeResidualNorm(NEP nep,PetscInt i,PetscReal *r) {return NEPComputeError(nep,i,NEP_ERROR_ABSOLUTE,r);}
SLEPC_EXTERN PetscErrorCode NEPGetErrorEstimate(NEP,PetscInt,PetscReal*);

SLEPC_EXTERN PetscErrorCode NEPComputeFunction(NEP,PetscScalar,Mat,Mat);
SLEPC_EXTERN PetscErrorCode NEPComputeJacobian(NEP,PetscScalar,Mat);
SLEPC_EXTERN PetscErrorCode NEPApplyFunction(NEP,PetscScalar,Vec,Vec,Vec,Mat,Mat);
SLEPC_EXTERN PetscErrorCode NEPApplyAdjoint(NEP,PetscScalar,Vec,Vec,Vec,Mat,Mat);
SLEPC_EXTERN PetscErrorCode NEPApplyJacobian(NEP,PetscScalar,Vec,Vec,Vec,Mat);
SLEPC_EXTERN PetscErrorCode NEPProjectOperator(NEP,PetscInt,PetscInt);

SLEPC_EXTERN PetscErrorCode NEPMonitor(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPMonitorSet(NEP,PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*),void*,PetscErrorCode (*)(void**));
SLEPC_EXTERN PetscErrorCode NEPMonitorSetFromOptions(NEP,const char*,const char*,const char*,PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscBool);
SLEPC_EXTERN PetscErrorCode NEPConvMonitorSetFromOptions(NEP,const char*,const char*,const char*,PetscErrorCode (*)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,SlepcConvMonitor));
SLEPC_EXTERN PetscErrorCode NEPMonitorCancel(NEP);
SLEPC_EXTERN PetscErrorCode NEPGetMonitorContext(NEP,void **);
SLEPC_EXTERN PetscErrorCode NEPGetIterationNumber(NEP,PetscInt*);

SLEPC_EXTERN PetscErrorCode NEPSetInitialSpace(NEP,PetscInt,Vec[]);
SLEPC_EXTERN PetscErrorCode NEPSetWhichEigenpairs(NEP,NEPWhich);
SLEPC_EXTERN PetscErrorCode NEPGetWhichEigenpairs(NEP,NEPWhich*);
SLEPC_EXTERN PetscErrorCode NEPSetTwoSided(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPGetTwoSided(NEP,PetscBool*);

SLEPC_EXTERN PetscErrorCode NEPApplyResolvent(NEP,RG,PetscScalar,Vec,Vec);

SLEPC_EXTERN PetscErrorCode NEPSetEigenvalueComparison(NEP,PetscErrorCode (*func)(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*),void*);

SLEPC_EXTERN PetscErrorCode NEPMonitorAll(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode NEPMonitorFirst(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode NEPMonitorConverged(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,SlepcConvMonitor);
SLEPC_EXTERN PetscErrorCode NEPMonitorLGCreate(MPI_Comm,const char[],const char[],int,int,int,int,PetscDrawLG*);
SLEPC_EXTERN PetscErrorCode NEPMonitorLG(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
SLEPC_EXTERN PetscErrorCode NEPMonitorLGAll(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);

SLEPC_EXTERN PetscErrorCode NEPSetTrackAll(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPGetTrackAll(NEP,PetscBool*);

SLEPC_EXTERN PetscErrorCode NEPSetOptionsPrefix(NEP,const char*);
SLEPC_EXTERN PetscErrorCode NEPAppendOptionsPrefix(NEP,const char*);
SLEPC_EXTERN PetscErrorCode NEPGetOptionsPrefix(NEP,const char*[]);

SLEPC_EXTERN PetscErrorCode NEPGetConvergedReason(NEP,NEPConvergedReason*);

SLEPC_EXTERN PetscFunctionList NEPList;
SLEPC_EXTERN PetscErrorCode NEPRegister(const char[],PetscErrorCode(*)(NEP));

SLEPC_EXTERN PetscErrorCode NEPSetWorkVecs(NEP,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPAllocateSolution(NEP,PetscInt);

/* --------- options specific to particular eigensolvers -------- */

SLEPC_EXTERN PetscErrorCode NEPRIISetMaximumIterations(NEP,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPRIIGetMaximumIterations(NEP,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPRIISetLagPreconditioner(NEP,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPRIIGetLagPreconditioner(NEP,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPRIISetConstCorrectionTol(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPRIIGetConstCorrectionTol(NEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPRIISetHermitian(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPRIIGetHermitian(NEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPRIISetDeflationThreshold(NEP,PetscReal);
SLEPC_EXTERN PetscErrorCode NEPRIIGetDeflationThreshold(NEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode NEPRIISetKSP(NEP,KSP);
SLEPC_EXTERN PetscErrorCode NEPRIIGetKSP(NEP,KSP*);

SLEPC_EXTERN PetscErrorCode NEPSLPSetDeflationThreshold(NEP,PetscReal);
SLEPC_EXTERN PetscErrorCode NEPSLPGetDeflationThreshold(NEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode NEPSLPSetEPS(NEP,EPS);
SLEPC_EXTERN PetscErrorCode NEPSLPGetEPS(NEP,EPS*);
SLEPC_EXTERN PetscErrorCode NEPSLPSetEPSLeft(NEP,EPS);
SLEPC_EXTERN PetscErrorCode NEPSLPGetEPSLeft(NEP,EPS*);
SLEPC_EXTERN PetscErrorCode NEPSLPSetKSP(NEP,KSP);
SLEPC_EXTERN PetscErrorCode NEPSLPGetKSP(NEP,KSP*);

SLEPC_EXTERN PetscErrorCode NEPNArnoldiSetKSP(NEP,KSP);
SLEPC_EXTERN PetscErrorCode NEPNArnoldiGetKSP(NEP,KSP*);
SLEPC_EXTERN PetscErrorCode NEPNArnoldiSetLagPreconditioner(NEP,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPNArnoldiGetLagPreconditioner(NEP,PetscInt*);

#if defined(PETSC_USE_COMPLEX)
SLEPC_EXTERN PetscErrorCode NEPCISSSetSizes(NEP,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPCISSGetSizes(NEP,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPCISSSetThreshold(NEP,PetscReal,PetscReal);
SLEPC_EXTERN PetscErrorCode NEPCISSGetThreshold(NEP,PetscReal*,PetscReal*);
SLEPC_EXTERN PetscErrorCode NEPCISSSetRefinement(NEP,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPCISSGetRefinement(NEP,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPCISSGetKSPs(NEP,PetscInt*,KSP**);
#else
#define SlepcNEPCISSUnavailable(nep) do { \
    PetscFunctionBegin; \
    SETERRQ1(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"%s() not available with real scalars",PETSC_FUNCTION_NAME); \
    } while (0)
PETSC_STATIC_INLINE PetscErrorCode NEPCISSSetSizes(NEP nep,PETSC_UNUSED PetscInt ip,PETSC_UNUSED PetscInt bs,PETSC_UNUSED PetscInt ms,PETSC_UNUSED PetscInt npart,PETSC_UNUSED PetscInt bsmax,PETSC_UNUSED PetscBool realmats) {SlepcNEPCISSUnavailable(nep);}
PETSC_STATIC_INLINE PetscErrorCode NEPCISSGetSizes(NEP nep,PETSC_UNUSED PetscInt *ip,PETSC_UNUSED PetscInt *bs,PETSC_UNUSED PetscInt *ms,PETSC_UNUSED PetscInt *npart,PETSC_UNUSED PetscInt *bsmak,PETSC_UNUSED PetscBool *realmats) {SlepcNEPCISSUnavailable(nep);}
PETSC_STATIC_INLINE PetscErrorCode NEPCISSSetThreshold(NEP nep,PETSC_UNUSED PetscReal delta,PETSC_UNUSED PetscReal spur) {SlepcNEPCISSUnavailable(nep);}
PETSC_STATIC_INLINE PetscErrorCode NEPCISSGetThreshold(NEP nep,PETSC_UNUSED PetscReal *delta,PETSC_UNUSED PetscReal *spur) {SlepcNEPCISSUnavailable(nep);}
PETSC_STATIC_INLINE PetscErrorCode NEPCISSSetRefinement(NEP nep,PETSC_UNUSED PetscInt inner,PETSC_UNUSED PetscInt blsize) {SlepcNEPCISSUnavailable(nep);}
PETSC_STATIC_INLINE PetscErrorCode NEPCISSGetRefinement(NEP nep,PETSC_UNUSED PetscInt *inner,PETSC_UNUSED PetscInt *blsize) {SlepcNEPCISSUnavailable(nep);}
PETSC_STATIC_INLINE PetscErrorCode NEPCISSGetKSPs(NEP nep,PETSC_UNUSED PetscInt *nsolve,PETSC_UNUSED KSP **ksp) {SlepcNEPCISSUnavailable(nep);}
#undef SlepcNEPCISSUnavailable
#endif

SLEPC_EXTERN PetscErrorCode NEPInterpolSetPEP(NEP,PEP);
SLEPC_EXTERN PetscErrorCode NEPInterpolGetPEP(NEP,PEP*);
SLEPC_EXTERN PetscErrorCode NEPInterpolSetInterpolation(NEP,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPInterpolGetInterpolation(NEP,PetscReal*,PetscInt*);

SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetSingularitiesFunction(NEP,PetscErrorCode (*)(NEP,PetscInt*,PetscScalar*,void*),void*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetSingularitiesFunction(NEP,PetscErrorCode (**)(NEP,PetscInt*,PetscScalar*,void*),void **);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetRestart(NEP,PetscReal);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetRestart(NEP,PetscReal*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetLocking(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetLocking(NEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetInterpolation(NEP,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetInterpolation(NEP,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetRKShifts(NEP,PetscInt,PetscScalar[]);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetRKShifts(NEP,PetscInt*,PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetKSPs(NEP,PetscInt*,KSP**);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetFullBasis(NEP,PetscBool);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetFullBasis(NEP,PetscBool*);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSSetEPS(NEP,EPS);
SLEPC_EXTERN PetscErrorCode NEPNLEIGSGetEPS(NEP,EPS*);

#endif

