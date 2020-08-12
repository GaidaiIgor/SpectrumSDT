/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the SLEPc object for solving linear matrix equations
*/

#if !defined(SLEPCLME_H)
#define SLEPCLME_H
#include <slepcbv.h>

SLEPC_EXTERN PetscErrorCode LMEInitializePackage(void);

/*S
    LME - SLEPc object that encapsulates functionality for linear matrix equations

    Level: beginner

.seealso:  LMECreate()
S*/
typedef struct _p_LME* LME;

/*J
    LMEType - String with the name of a method for solving linear matrix equations

    Level: beginner

.seealso: LMESetType(), LME
J*/
typedef const char* LMEType;
#define LMEKRYLOV   "krylov"

/* Logging support */
SLEPC_EXTERN PetscClassId LME_CLASSID;

/*E
    LMEProblemType - Determines the type of linear matrix equation

    Level: beginner

.seealso: LMESetProblemType(), LMEGetProblemType()
E*/
typedef enum { LME_LYAPUNOV,
               LME_SYLVESTER,
               LME_GEN_LYAPUNOV,
               LME_GEN_SYLVESTER,
               LME_DT_LYAPUNOV ,
               LME_STEIN} LMEProblemType;
SLEPC_EXTERN const char *LMEProblemTypes[];

SLEPC_EXTERN PetscErrorCode LMECreate(MPI_Comm,LME *);
SLEPC_EXTERN PetscErrorCode LMEDestroy(LME*);
SLEPC_EXTERN PetscErrorCode LMEReset(LME);
SLEPC_EXTERN PetscErrorCode LMESetType(LME,LMEType);
SLEPC_EXTERN PetscErrorCode LMEGetType(LME,LMEType*);
SLEPC_EXTERN PetscErrorCode LMESetProblemType(LME,LMEProblemType);
SLEPC_EXTERN PetscErrorCode LMEGetProblemType(LME,LMEProblemType*);
SLEPC_EXTERN PetscErrorCode LMESetCoefficients(LME,Mat,Mat,Mat,Mat);
SLEPC_EXTERN PetscErrorCode LMEGetCoefficients(LME,Mat*,Mat*,Mat*,Mat*);
SLEPC_EXTERN PetscErrorCode LMESetRHS(LME,Mat);
SLEPC_EXTERN PetscErrorCode LMEGetRHS(LME,Mat*);
SLEPC_EXTERN PetscErrorCode LMESetSolution(LME,Mat);
SLEPC_EXTERN PetscErrorCode LMEGetSolution(LME,Mat*);
SLEPC_EXTERN PetscErrorCode LMESetFromOptions(LME);
SLEPC_EXTERN PetscErrorCode LMESetUp(LME);
SLEPC_EXTERN PetscErrorCode LMESolve(LME);
SLEPC_EXTERN PetscErrorCode LMEView(LME,PetscViewer);
SLEPC_EXTERN PetscErrorCode LMEViewFromOptions(LME,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode LMEReasonView(LME,PetscViewer);
SLEPC_EXTERN PetscErrorCode LMEReasonViewFromOptions(LME);

SLEPC_EXTERN PetscErrorCode LMESetBV(LME,BV);
SLEPC_EXTERN PetscErrorCode LMEGetBV(LME,BV*);
SLEPC_EXTERN PetscErrorCode LMESetTolerances(LME,PetscReal,PetscInt);
SLEPC_EXTERN PetscErrorCode LMEGetTolerances(LME,PetscReal*,PetscInt*);
SLEPC_EXTERN PetscErrorCode LMESetDimensions(LME,PetscInt);
SLEPC_EXTERN PetscErrorCode LMEGetDimensions(LME,PetscInt*);

SLEPC_EXTERN PetscErrorCode LMEMonitor(LME,PetscInt,PetscReal);
SLEPC_EXTERN PetscErrorCode LMEMonitorSet(LME,PetscErrorCode (*)(LME,PetscInt,PetscReal,void*),void*,PetscErrorCode (*)(void**));
SLEPC_EXTERN PetscErrorCode LMEMonitorSetFromOptions(LME,const char*,const char*,const char*,PetscErrorCode (*)(LME,PetscInt,PetscReal,PetscViewerAndFormat*));
SLEPC_EXTERN PetscErrorCode LMEMonitorCancel(LME);
SLEPC_EXTERN PetscErrorCode LMEGetMonitorContext(LME,void **);
SLEPC_EXTERN PetscErrorCode LMEGetIterationNumber(LME,PetscInt*);

SLEPC_EXTERN PetscErrorCode LMEGetErrorEstimate(LME,PetscReal*);
SLEPC_EXTERN PetscErrorCode LMEComputeError(LME,PetscReal*);
SLEPC_EXTERN PetscErrorCode LMESetErrorIfNotConverged(LME,PetscBool);
SLEPC_EXTERN PetscErrorCode LMEGetErrorIfNotConverged(LME,PetscBool*);

SLEPC_EXTERN PetscErrorCode LMEDenseLyapunov(LME,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt);
SLEPC_EXTERN PetscErrorCode LMEDenseHessLyapunovChol(LME,PetscInt,PetscScalar*,PetscInt,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscReal*);
PETSC_DEPRECATED_FUNCTION("Use LMEDenseHessLyapunovChol()") PETSC_STATIC_INLINE PetscErrorCode LMEDenseLyapunovChol(LME lme,PetscScalar *H,PetscInt m,PetscInt ldh,PetscScalar *r,PetscScalar *L,PetscInt ldl,PetscReal *res) {return LMEDenseHessLyapunovChol(lme,m,H,ldh,1,r,m,L,ldl,res);}

SLEPC_EXTERN PetscErrorCode LMEMonitorDefault(LME,PetscInt,PetscReal,PetscViewerAndFormat*);
SLEPC_EXTERN PetscErrorCode LMEMonitorLGCreate(MPI_Comm,const char[],const char[],int,int,int,int,PetscDrawLG*);
SLEPC_EXTERN PetscErrorCode LMEMonitorLG(LME,PetscInt,PetscReal,void*);

SLEPC_EXTERN PetscErrorCode LMESetOptionsPrefix(LME,const char*);
SLEPC_EXTERN PetscErrorCode LMEAppendOptionsPrefix(LME,const char*);
SLEPC_EXTERN PetscErrorCode LMEGetOptionsPrefix(LME,const char*[]);

/*E
    LMEConvergedReason - reason a matrix function iteration was said to
         have converged or diverged

    Level: intermediate

.seealso: LMESolve(), LMEGetConvergedReason(), LMESetTolerances()
E*/
typedef enum {/* converged */
              LME_CONVERGED_TOL                =  1,
              /* diverged */
              LME_DIVERGED_ITS                 = -1,
              LME_DIVERGED_BREAKDOWN           = -2,
              LME_CONVERGED_ITERATING          =  0} LMEConvergedReason;
SLEPC_EXTERN const char *const*LMEConvergedReasons;

SLEPC_EXTERN PetscErrorCode LMEGetConvergedReason(LME,LMEConvergedReason *);

SLEPC_EXTERN PetscFunctionList LMEList;
SLEPC_EXTERN PetscErrorCode LMERegister(const char[],PetscErrorCode(*)(LME));

SLEPC_EXTERN PetscErrorCode LMEAllocateSolution(LME,PetscInt);

#endif

