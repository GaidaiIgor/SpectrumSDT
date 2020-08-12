/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for the direct solver object in SLEPc
*/

#if !defined(SLEPCDS_H)
#define SLEPCDS_H
#include <slepcsc.h>
#include <slepcfn.h>

#define DS_MAX_SOLVE 6

SLEPC_EXTERN PetscErrorCode DSInitializePackage(void);
/*S
    DS - Direct solver (or dense system), to represent low-dimensional
    eigenproblems that must be solved within iterative solvers. This is an
    auxiliary object and is not normally needed by application programmers.

    Level: beginner

.seealso:  DSCreate()
S*/
typedef struct _p_DS* DS;

/*J
    DSType - String with the name of the type of direct solver. Roughly,
    there are as many types as problem types are available within SLEPc.

    Level: advanced

.seealso: DSSetType(), DS
J*/
typedef const char* DSType;
#define DSHEP    "hep"
#define DSNHEP   "nhep"
#define DSGHEP   "ghep"
#define DSGHIEP  "ghiep"
#define DSGNHEP  "gnhep"
#define DSSVD    "svd"
#define DSPEP    "pep"
#define DSNEP    "nep"

/* Logging support */
SLEPC_EXTERN PetscClassId DS_CLASSID;

/*E
    DSStateType - Indicates in which state the direct solver is

    Level: advanced

.seealso: DSSetState()
E*/
typedef enum { DS_STATE_RAW,
               DS_STATE_INTERMEDIATE,
               DS_STATE_CONDENSED,
               DS_STATE_TRUNCATED } DSStateType;
SLEPC_EXTERN const char *DSStateTypes[];

/*E
    DSMatType - Used to refer to one of the matrices stored internally in DS

    Notes:
    The matrices preferently refer to
+   DS_MAT_A  - first matrix of eigenproblem/singular value problem
.   DS_MAT_B  - second matrix of a generalized eigenproblem
.   DS_MAT_C  - third matrix of a quadratic eigenproblem
.   DS_MAT_T  - tridiagonal matrix
.   DS_MAT_D  - diagonal matrix
.   DS_MAT_Q  - orthogonal matrix of (right) Schur vectors
.   DS_MAT_Z  - orthogonal matrix of left Schur vectors
.   DS_MAT_X  - right eigenvectors
.   DS_MAT_Y  - left eigenvectors
.   DS_MAT_U  - left singular vectors
.   DS_MAT_VT - right singular vectors
.   DS_MAT_W  - workspace matrix
-   DS_MAT_Ex - extra matrices (x=0,..,9)

    All matrices can have space to hold ld x ld elements, except for
    DS_MAT_T that has space for 3 x ld elements (ld = leading dimension)
    and DS_MAT_D that has space for just ld elements.

    In DSPEP problems, matrices A, B, W can have space for d*ld x d*ld,
    where d is the polynomial degree, and X can have ld x d*ld.

    Level: advanced

.seealso: DSAllocate(), DSGetArray(), DSGetArrayReal(), DSVectors()
E*/
typedef enum { DS_MAT_A,
               DS_MAT_B,
               DS_MAT_C,
               DS_MAT_T,
               DS_MAT_D,
               DS_MAT_Q,
               DS_MAT_Z,
               DS_MAT_X,
               DS_MAT_Y,
               DS_MAT_U,
               DS_MAT_VT,
               DS_MAT_W,
               DS_MAT_E0,
               DS_MAT_E1,
               DS_MAT_E2,
               DS_MAT_E3,
               DS_MAT_E4,
               DS_MAT_E5,
               DS_MAT_E6,
               DS_MAT_E7,
               DS_MAT_E8,
               DS_MAT_E9,
               DS_NUM_MAT } DSMatType;

/* Convenience for indexing extra matrices */
SLEPC_EXTERN DSMatType DSMatExtra[];
#define DS_NUM_EXTRA  10

/*E
    DSParallelType - Indicates the parallel mode that the direct solver will use

    Level: advanced

.seealso: DSSetParallel()
E*/
typedef enum { DS_PARALLEL_REDUNDANT,
               DS_PARALLEL_SYNCHRONIZED } DSParallelType;
SLEPC_EXTERN const char *DSParallelTypes[];

SLEPC_EXTERN PetscErrorCode DSCreate(MPI_Comm,DS*);
SLEPC_EXTERN PetscErrorCode DSSetType(DS,DSType);
SLEPC_EXTERN PetscErrorCode DSGetType(DS,DSType*);
SLEPC_EXTERN PetscErrorCode DSSetOptionsPrefix(DS,const char *);
SLEPC_EXTERN PetscErrorCode DSAppendOptionsPrefix(DS,const char *);
SLEPC_EXTERN PetscErrorCode DSGetOptionsPrefix(DS,const char *[]);
SLEPC_EXTERN PetscErrorCode DSSetFromOptions(DS);
SLEPC_EXTERN PetscErrorCode DSView(DS,PetscViewer);
SLEPC_EXTERN PetscErrorCode DSViewFromOptions(DS,PetscObject,const char[]);
SLEPC_EXTERN PetscErrorCode DSViewMat(DS,PetscViewer,DSMatType);
SLEPC_EXTERN PetscErrorCode DSDestroy(DS*);
SLEPC_EXTERN PetscErrorCode DSReset(DS);
SLEPC_EXTERN PetscErrorCode DSDuplicate(DS,DS*);

SLEPC_EXTERN PetscErrorCode DSAllocate(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSGetLeadingDimension(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSSetState(DS,DSStateType);
SLEPC_EXTERN PetscErrorCode DSGetState(DS,DSStateType*);
SLEPC_EXTERN PetscErrorCode DSSetDimensions(DS,PetscInt,PetscInt,PetscInt,PetscInt);
SLEPC_EXTERN PetscErrorCode DSGetDimensions(DS,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSSetBlockSize(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSGetBlockSize(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSTruncate(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSSetIdentity(DS,DSMatType);
SLEPC_EXTERN PetscErrorCode DSSetMethod(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSGetMethod(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSSetParallel(DS,DSParallelType);
SLEPC_EXTERN PetscErrorCode DSGetParallel(DS,DSParallelType*);
SLEPC_EXTERN PetscErrorCode DSSetCompact(DS,PetscBool);
SLEPC_EXTERN PetscErrorCode DSGetCompact(DS,PetscBool*);
SLEPC_EXTERN PetscErrorCode DSSetExtraRow(DS,PetscBool);
SLEPC_EXTERN PetscErrorCode DSGetExtraRow(DS,PetscBool*);
SLEPC_EXTERN PetscErrorCode DSSetRefined(DS,PetscBool);
SLEPC_EXTERN PetscErrorCode DSGetRefined(DS,PetscBool*);
SLEPC_EXTERN PetscErrorCode DSGetMat(DS,DSMatType,Mat*);
SLEPC_EXTERN PetscErrorCode DSRestoreMat(DS,DSMatType,Mat*);
SLEPC_EXTERN PetscErrorCode DSGetArray(DS,DSMatType,PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode DSRestoreArray(DS,DSMatType,PetscScalar*[]);
SLEPC_EXTERN PetscErrorCode DSGetArrayReal(DS,DSMatType,PetscReal*[]);
SLEPC_EXTERN PetscErrorCode DSRestoreArrayReal(DS,DSMatType,PetscReal*[]);
SLEPC_EXTERN PetscErrorCode DSVectors(DS,DSMatType,PetscInt*,PetscReal*);
SLEPC_EXTERN PetscErrorCode DSSolve(DS,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode DSSort(DS,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSSortWithPermutation(DS,PetscInt*,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode DSSynchronize(DS,PetscScalar*,PetscScalar*);
SLEPC_EXTERN PetscErrorCode DSCopyMat(DS,DSMatType,PetscInt,PetscInt,Mat,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool);
SLEPC_EXTERN PetscErrorCode DSMatGetSize(DS,DSMatType,PetscInt*,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSMatIsHermitian(DS,DSMatType,PetscBool*);
SLEPC_EXTERN PetscErrorCode DSSetSlepcSC(DS,SlepcSC);
SLEPC_EXTERN PetscErrorCode DSGetSlepcSC(DS,SlepcSC*);
SLEPC_EXTERN PetscErrorCode DSUpdateExtraRow(DS);
SLEPC_EXTERN PetscErrorCode DSCond(DS,PetscReal*);
SLEPC_EXTERN PetscErrorCode DSTranslateHarmonic(DS,PetscScalar,PetscReal,PetscBool,PetscScalar*,PetscReal*);
SLEPC_EXTERN PetscErrorCode DSTranslateRKS(DS,PetscScalar);
SLEPC_EXTERN PetscErrorCode DSOrthogonalize(DS,DSMatType,PetscInt,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSPseudoOrthogonalize(DS,DSMatType,PetscInt,PetscReal*,PetscInt*,PetscReal*);

/* --------- options specific to particular solvers -------- */

SLEPC_EXTERN PetscErrorCode DSPEPSetDegree(DS,PetscInt);
SLEPC_EXTERN PetscErrorCode DSPEPGetDegree(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSPEPSetCoefficients(DS,PetscReal*);
SLEPC_EXTERN PetscErrorCode DSPEPGetCoefficients(DS,PetscReal**);

SLEPC_EXTERN PetscErrorCode DSNEPSetFN(DS,PetscInt,FN*);
SLEPC_EXTERN PetscErrorCode DSNEPGetFN(DS,PetscInt,FN*);
SLEPC_EXTERN PetscErrorCode DSNEPGetNumFN(DS,PetscInt*);
SLEPC_EXTERN PetscErrorCode DSNEPSetComputeMatrixFunction(DS,PetscErrorCode (*fun)(DS,PetscScalar,PetscBool,DSMatType,void*),void *ctx);

SLEPC_EXTERN PetscFunctionList DSList;
SLEPC_EXTERN PetscErrorCode DSRegister(const char[],PetscErrorCode(*)(DS));

#endif
