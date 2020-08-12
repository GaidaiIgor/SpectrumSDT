/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for SLEPc's singular value solvers
*/

#if !defined(SLEPCSVD_H)
#define SLEPCSVD_H
#include <slepceps.h>
#include <slepcbv.h>
#include <slepcds.h>

SLEPC_EXTERN PetscErrorCode SVDInitializePackage(void);

/*S
     SVD - Abstract SLEPc object that manages all the singular value
     problem solvers.

   Level: beginner

.seealso:  SVDCreate()
S*/
typedef struct _p_SVD* SVD;

/*J
    SVDType - String with the name of a SLEPc singular value solver

   Level: beginner

.seealso: SVDSetType(), SVD
J*/
typedef const char* SVDType;
#define SVDCROSS       "cross"
#define SVDCYCLIC      "cyclic"
#define SVDLAPACK      "lapack"
#define SVDLANCZOS     "lanczos"
#define SVDTRLANCZOS   "trlanczos"
#define SVDPRIMME      "primme"

/* Logging support */
SLEPC_EXTERN PetscClassId SVD_CLASSID;

/*E
    SVDWhich - Determines whether largest or smallest singular triplets
    are to be computed

    Level: intermediate

.seealso: SVDSetWhichSingularTriplets(), SVDGetWhichSingularTriplets()
E*/
typedef enum { SVD_LARGEST,
               SVD_SMALLEST } SVDWhich;

/*E
    SVDErrorType - The error type used to assess accuracy of computed solutions

    Level: intermediate

.seealso: SVDComputeError()
E*/
typedef enum { SVD_ERROR_ABSOLUTE,
               SVD_ERROR_RELATIVE } SVDErrorType;
SLEPC_EXTERN const char *SVDErrorTypes[];

/*E
    SVDConv - Determines the convergence test

    Level: intermediate

.seealso: SVDSetConvergenceTest(), SVDSetConvergenceTestFunction()
E*/
typedef enum { SVD_CONV_ABS,
               SVD_CONV_REL,
               SVD_CONV_USER } SVDConv;

/*E
    SVDStop - Determines the stopping test

    Level: advanced

.seealso: SVDSetStoppingTest(), SVDSetStoppingTestFunction()
E*/
typedef enum { SVD_STOP_BASIC,
               SVD_STOP_USER } SVDStop;

/*E
    SVDConvergedReason - Reason a singular value solver was said to
         have converged or diverged

   Level: intermediate

.seealso: SVDSolve(), SVDGetConvergedReason(), SVDSetTolerances()
E*/
typedef enum {/* converged */
              SVD_CONVERGED_TOL                =  1,
              SVD_CONVERGED_USER               =  2,
              /* diverged */
              SVD_DIVERGED_ITS                 = -1,
              SVD_DIVERGED_BREAKDOWN           = -2,
              SVD_CONVERGED_ITERATING          =  0 } SVDConvergedReason;
SLEPC_EXTERN const char *const*SVDConvergedReasons;

SLEPC_EXTERN PetscErrorCode SVDCreate(MPI_Comm,SVD*);
SLEPC_EXTERN PetscErrorCode SVDSetBV(SVD,BV,BV);
SLEPC_EXTERN PetscErrorCode SVDGetBV(SVD,BV*,BV*);
SLEPC_EXTERN PetscErrorCode SVDSetDS(SVD,DS);
SLEPC_EXTERN PetscErrorCode SVDGetDS(SVD,DS*);
SLEPC_EXTERN PetscErrorCode SVDSetType(SVD,SVDType);
SLEPC_EXTERN PetscErrorCode SVDGetType(SVD,SVDType*);
SLEPC_EXTERN PetscErrorCode SVDSetOperator(SVD,Mat);
SLEPC_EXTERN PetscErrorCode SVDGetOperator(SVD,Mat*);
SLEPC_EXTERN PetscErrorCode SVDSetInitialSpaces(SVD,PetscInt,Vec[],PetscInt,Vec[]);
PETSC_DEPRECATED_FUNCTION("Use SVDSetInitialSpaces()") PETSC_STATIC_INLINE PetscErrorCode SVDSetInitialSpace(SVD svd,PetscInt nr,Vec *isr) {return SVDSetInitialSpaces(svd,nr,isr,0,NULL);}
PETSC_DEPRECATED_FUNCTION("Use SVDSetInitialSpaces()") PETSC_STATIC_INLINE PetscErrorCode SVDSetInitialSpaceLeft(SVD svd,PetscInt nl,Vec *isl) {return SVDSetInitialSpaces(svd,0,NULL,nl,isl);}
SLEPC_EXTERN PetscErrorCode SVDSetImplicitTranspose(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDGetImplicitTranspose(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDSetDimensions(SVD,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDGetDimensions(SVD,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDSetTolerances(SVD,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDGetTolerances(SVD,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDSetWhichSingularTriplets(SVD,SVDWhich);
SLEPC_EXTERN PetscErrorCode SVDGetWhichSingularTriplets(SVD,SVDWhich*);
SLEPC_EXTERN PetscErrorCode SVDSetFromOptions(SVD);
SLEPC_EXTERN PetscErrorCode SVDSetOptionsPrefix(SVD,const char*);
SLEPC_EXTERN PetscErrorCode SVDAppendOptionsPrefix(SVD,const char*);
SLEPC_EXTERN PetscErrorCode SVDGetOptionsPrefix(SVD,const char*[]);
SLEPC_EXTERN PetscErrorCode SVDSetUp(SVD);
SLEPC_EXTERN PetscErrorCode SVDSolve(SVD);
SLEPC_EXTERN PetscErrorCode SVDGetIterationNumber(SVD,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDSetConvergenceTestFunction(SVD,PetscErrorCode (*)(SVD,PetscReal,PetscReal,PetscReal*,void*),void*,PetscErrorCode (*)(void*));
SLEPC_EXTERN PetscErrorCode SVDSetConvergenceTest(SVD,SVDConv);
SLEPC_EXTERN PetscErrorCode SVDGetConvergenceTest(SVD,SVDConv*);
SLEPC_EXTERN PetscErrorCode SVDConvergedAbsolute(SVD,PetscReal,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode SVDConvergedRelative(SVD,PetscReal,PetscReal,PetscReal*,void*);
SLEPC_EXTERN PetscErrorCode SVDSetStoppingTestFunction(SVD,PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscInt,PetscInt,SVDConvergedReason*,void*),void*,PetscErrorCode (*)(void*));
SLEPC_EXTERN PetscErrorCode SVDSetStoppingTest(SVD,SVDStop);
SLEPC_EXTERN PetscErrorCode SVDGetStoppingTest(SVD,SVDStop*);
SLEPC_EXTERN PetscErrorCode SVDStoppingBasic(SVD,PetscInt,PetscInt,PetscInt,PetscInt,SVDConvergedReason*,void*);
SLEPC_EXTERN PetscErrorCode SVDGetConvergedReason(SVD,SVDConvergedReason*);
SLEPC_EXTERN PetscErrorCode SVDGetConverged(SVD,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDGetSingularTriplet(SVD,PetscInt,PetscReal*,Vec,Vec);
SLEPC_EXTERN PetscErrorCode SVDComputeError(SVD,PetscInt,SVDErrorType,PetscReal*);
PETSC_DEPRECATED_FUNCTION("Use SVDComputeError()") PETSC_STATIC_INLINE PetscErrorCode SVDComputeRelativeError(SVD svd,PetscInt i,PetscReal *r) {return SVDComputeError(svd,i,SVD_ERROR_RELATIVE,r);}
PETSC_DEPRECATED_FUNCTION("Use SVDComputeError() with SVD_ERROR_ABSOLUTE") PETSC_STATIC_INLINE PetscErrorCode SVDComputeResidualNorms(SVD svd,PetscInt i,PetscReal *r1,PETSC_UNUSED PetscReal *r2) {return SVDComputeError(svd,i,SVD_ERROR_ABSOLUTE,r1);}
SLEPC_EXTERN PetscErrorCode SVDView(SVD,PetscViewer);
SLEPC_EXTERN PetscErrorCode SVDViewFromOptions(SVD,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode SVDErrorView(SVD,SVDErrorType,PetscViewer);
PETSC_DEPRECATED_FUNCTION("Use SVDErrorView()") PETSC_STATIC_INLINE PetscErrorCode SVDPrintSolution(SVD svd,PetscViewer v) {return SVDErrorView(svd,SVD_ERROR_RELATIVE,v);}
SLEPC_EXTERN PetscErrorCode SVDErrorViewFromOptions(SVD);
SLEPC_EXTERN PetscErrorCode SVDReasonView(SVD,PetscViewer);
SLEPC_EXTERN PetscErrorCode SVDReasonViewFromOptions(SVD);
SLEPC_EXTERN PetscErrorCode SVDValuesView(SVD,PetscViewer);
SLEPC_EXTERN PetscErrorCode SVDValuesViewFromOptions(SVD);
SLEPC_EXTERN PetscErrorCode SVDVectorsView(SVD,PetscViewer);
SLEPC_EXTERN PetscErrorCode SVDVectorsViewFromOptions(SVD);
SLEPC_EXTERN PetscErrorCode SVDDestroy(SVD*);
SLEPC_EXTERN PetscErrorCode SVDReset(SVD);

SLEPC_EXTERN PetscErrorCode SVDMonitor(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDMonitorSet(SVD,PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*),void*,PetscErrorCode (*)(void**));
SLEPC_EXTERN PetscErrorCode SVDMonitorSetFromOptions(SVD,const char*,const char*,const char*,PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscBool);
SLEPC_EXTERN PetscErrorCode SVDConvMonitorSetFromOptions(SVD,const char*,const char*,const char*,PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,SlepcConvMonitor));
SLEPC_EXTERN PetscErrorCode SVDMonitorCancel(SVD);
SLEPC_EXTERN PetscErrorCode SVDGetMonitorContext(SVD,void **);
SLEPC_EXTERN PetscErrorCode SVDMonitorAll(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode SVDMonitorFirst(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode SVDMonitorConverged(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,SlepcConvMonitor);
SLEPC_EXTERN PetscErrorCode SVDMonitorLGCreate(MPI_Comm,const char[],const char[],int,int,int,int,PetscDrawLG*);
SLEPC_EXTERN PetscErrorCode SVDMonitorLG(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);
SLEPC_EXTERN PetscErrorCode SVDMonitorLGAll(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);

SLEPC_EXTERN PetscErrorCode SVDSetTrackAll(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDGetTrackAll(SVD,PetscBool*);

SLEPC_EXTERN PetscFunctionList SVDList;
SLEPC_EXTERN PetscErrorCode SVDRegister(const char[],PetscErrorCode(*)(SVD));

SLEPC_EXTERN PetscErrorCode SVDAllocateSolution(SVD,PetscInt);

/* --------- options specific to particular solvers -------- */

SLEPC_EXTERN PetscErrorCode SVDCrossSetExplicitMatrix(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDCrossGetExplicitMatrix(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDCrossSetEPS(SVD,EPS);
SLEPC_EXTERN PetscErrorCode SVDCrossGetEPS(SVD,EPS*);

SLEPC_EXTERN PetscErrorCode SVDCyclicSetExplicitMatrix(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDCyclicGetExplicitMatrix(SVD,PetscBool*);
SLEPC_EXTERN PetscErrorCode SVDCyclicSetEPS(SVD,EPS);
SLEPC_EXTERN PetscErrorCode SVDCyclicGetEPS(SVD,EPS*);

SLEPC_EXTERN PetscErrorCode SVDLanczosSetOneSide(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDLanczosGetOneSide(SVD,PetscBool*);

SLEPC_EXTERN PetscErrorCode SVDTRLanczosSetOneSide(SVD,PetscBool);
SLEPC_EXTERN PetscErrorCode SVDTRLanczosGetOneSide(SVD,PetscBool*);

/*E
    SVDPRIMMEMethod - determines the SVD method selected in the PRIMME library

    Level: advanced

.seealso: SVDPRIMMESetMethod(), SVDPRIMMEGetMethod()
E*/
typedef enum { SVD_PRIMME_HYBRID=1,
               SVD_PRIMME_NORMALEQUATIONS,
               SVD_PRIMME_AUGMENTED } SVDPRIMMEMethod;
SLEPC_EXTERN const char *SVDPRIMMEMethods[];

SLEPC_EXTERN PetscErrorCode SVDPRIMMESetBlockSize(SVD,PetscInt);
SLEPC_EXTERN PetscErrorCode SVDPRIMMEGetBlockSize(SVD,PetscInt*);
SLEPC_EXTERN PetscErrorCode SVDPRIMMESetMethod(SVD,SVDPRIMMEMethod);
SLEPC_EXTERN PetscErrorCode SVDPRIMMEGetMethod(SVD,SVDPRIMMEMethod*);

#endif
