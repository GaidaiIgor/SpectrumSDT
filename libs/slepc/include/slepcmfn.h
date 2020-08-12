/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the SLEPc matrix function solver object
*/

#if !defined(SLEPCMFN_H)
#define SLEPCMFN_H
#include <slepcbv.h>
#include <slepcfn.h>

SLEPC_EXTERN PetscErrorCode MFNInitializePackage(void);

/*S
    MFN - SLEPc object that encapsulates functionality for matrix functions.

    Level: beginner

.seealso:  MFNCreate()
S*/
typedef struct _p_MFN* MFN;

/*J
    MFNType - String with the name of a method for computing matrix functions.

    Level: beginner

.seealso: MFNSetType(), MFN
J*/
typedef const char* MFNType;
#define MFNKRYLOV   "krylov"
#define MFNEXPOKIT  "expokit"

/* Logging support */
SLEPC_EXTERN PetscClassId MFN_CLASSID;

SLEPC_EXTERN PetscErrorCode MFNCreate(MPI_Comm,MFN *);
SLEPC_EXTERN PetscErrorCode MFNDestroy(MFN*);
SLEPC_EXTERN PetscErrorCode MFNReset(MFN);
SLEPC_EXTERN PetscErrorCode MFNSetType(MFN,MFNType);
SLEPC_EXTERN PetscErrorCode MFNGetType(MFN,MFNType*);
SLEPC_EXTERN PetscErrorCode MFNSetOperator(MFN,Mat);
SLEPC_EXTERN PetscErrorCode MFNGetOperator(MFN,Mat*);
SLEPC_EXTERN PetscErrorCode MFNSetFromOptions(MFN);
SLEPC_EXTERN PetscErrorCode MFNSetUp(MFN);
SLEPC_EXTERN PetscErrorCode MFNSolve(MFN,Vec,Vec);
SLEPC_EXTERN PetscErrorCode MFNSolveTranspose(MFN,Vec,Vec);
SLEPC_EXTERN PetscErrorCode MFNView(MFN,PetscViewer);
SLEPC_EXTERN PetscErrorCode MFNViewFromOptions(MFN,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode MFNReasonView(MFN,PetscViewer);
SLEPC_EXTERN PetscErrorCode MFNReasonViewFromOptions(MFN);

SLEPC_EXTERN PetscErrorCode MFNSetBV(MFN,BV);
SLEPC_EXTERN PetscErrorCode MFNGetBV(MFN,BV*);
SLEPC_EXTERN PetscErrorCode MFNSetFN(MFN,FN);
SLEPC_EXTERN PetscErrorCode MFNGetFN(MFN,FN*);
SLEPC_EXTERN PetscErrorCode MFNSetTolerances(MFN,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode MFNGetTolerances(MFN,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode MFNSetDimensions(MFN,PetscInt);
SLEPC_EXTERN PetscErrorCode MFNGetDimensions(MFN,PetscInt*);

SLEPC_EXTERN PetscErrorCode MFNMonitor(MFN,PetscInt,PetscReal);
SLEPC_EXTERN PetscErrorCode MFNMonitorSet(MFN,PetscErrorCode (*)(MFN,PetscInt,PetscReal,void*),void*,PetscErrorCode (*)(void**));
SLEPC_EXTERN PetscErrorCode MFNMonitorSetFromOptions(MFN,const char*,const char*,const char*,PetscErrorCode (*)(MFN,PetscInt,PetscReal,PetscViewerAndFormat*));
SLEPC_EXTERN PetscErrorCode MFNMonitorCancel(MFN);
SLEPC_EXTERN PetscErrorCode MFNGetMonitorContext(MFN,void **);
SLEPC_EXTERN PetscErrorCode MFNGetIterationNumber(MFN,PetscInt*);

SLEPC_EXTERN PetscErrorCode MFNSetErrorIfNotConverged(MFN,PetscBool);
SLEPC_EXTERN PetscErrorCode MFNGetErrorIfNotConverged(MFN,PetscBool*);

SLEPC_EXTERN PetscErrorCode MFNMonitorDefault(MFN,PetscInt,PetscReal,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode MFNMonitorLGCreate(MPI_Comm,const char[],const char[],int,int,int,int,PetscDrawLG*);
SLEPC_EXTERN PetscErrorCode MFNMonitorLG(MFN,PetscInt,PetscReal,void*);

SLEPC_EXTERN PetscErrorCode MFNSetOptionsPrefix(MFN,const char*);
SLEPC_EXTERN PetscErrorCode MFNAppendOptionsPrefix(MFN,const char*);
SLEPC_EXTERN PetscErrorCode MFNGetOptionsPrefix(MFN,const char*[]);

/*E
    MFNConvergedReason - reason a matrix function iteration was said to
         have converged or diverged

    Level: intermediate

.seealso: MFNSolve(), MFNGetConvergedReason(), MFNSetTolerances()
E*/
typedef enum {/* converged */
              MFN_CONVERGED_TOL                =  1,
              MFN_CONVERGED_ITS                =  2,
              /* diverged */
              MFN_DIVERGED_ITS                 = -1,
              MFN_DIVERGED_BREAKDOWN           = -2,
              MFN_CONVERGED_ITERATING          =  0} MFNConvergedReason;
SLEPC_EXTERN const char *const*MFNConvergedReasons;

SLEPC_EXTERN PetscErrorCode MFNGetConvergedReason(MFN,MFNConvergedReason *);

SLEPC_EXTERN PetscFunctionList MFNList;
SLEPC_EXTERN PetscErrorCode MFNRegister(const char[],PetscErrorCode(*)(MFN));

SLEPC_EXTERN PetscErrorCode MFNAllocateSolution(MFN,PetscInt);

#endif

