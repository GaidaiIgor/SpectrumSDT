/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if !defined(SLEPCSTIMPL_H)
#define SLEPCSTIMPL_H

#include <slepcst.h>
#include <slepc/private/slepcimpl.h>

SLEPC_EXTERN PetscBool STRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode STRegisterAll(void);
SLEPC_EXTERN PetscLogEvent ST_SetUp,ST_ComputeOperator,ST_Apply,ST_ApplyTranspose,ST_MatSetUp,ST_MatMult,ST_MatMultTranspose,ST_MatSolve,ST_MatSolveTranspose;

typedef struct _STOps *STOps;

struct _STOps {
  PetscErrorCode (*apply)(ST,Vec,Vec);
  PetscErrorCode (*applytrans)(ST,Vec,Vec);
  PetscErrorCode (*backtransform)(ST,PetscInt,PetscScalar*,PetscScalar*);
  PetscErrorCode (*setshift)(ST,PetscScalar);
  PetscErrorCode (*getbilinearform)(ST,Mat*);
  PetscErrorCode (*setup)(ST);
  PetscErrorCode (*computeoperator)(ST);
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,ST);
  PetscErrorCode (*postsolve)(ST);
  PetscErrorCode (*destroy)(ST);
  PetscErrorCode (*reset)(ST);
  PetscErrorCode (*view)(ST,PetscViewer);
  PetscErrorCode (*checknullspace)(ST,BV);
  PetscErrorCode (*setdefaultksp)(ST);
};

/*
     'Updated' state means STSetUp must be called because matrices have been
     modified, but the pattern is the same (hence reuse symbolic factorization)
*/
typedef enum { ST_STATE_INITIAL,
               ST_STATE_SETUP,
               ST_STATE_UPDATED } STStateType;

struct _p_ST {
  PETSCHEADER(struct _STOps);
  /*------------------------- User parameters --------------------------*/
  Mat              *A;               /* matrices that define the eigensystem */
  PetscInt         nmat;             /* number of user-provided matrices */
  PetscScalar      sigma;            /* value of the shift */
  PetscScalar      defsigma;         /* default value of the shift */
  STMatMode        matmode;          /* how the transformation matrix is handled */
  MatStructure     str;              /* whether matrices have the same pattern or not */
  PetscBool        transform;        /* whether transformed matrices are computed */
  Vec              D;                /* diagonal matrix for balancing */

  /*------------------------- Misc data --------------------------*/
  KSP              ksp;              /* linear solver used in some ST's */
  PetscBool        usesksp;          /* whether the KSP object is used or not */
  PetscInt         nwork;            /* number of work vectors */
  Vec              *work;            /* work vectors */
  Vec              wb;               /* balancing requires an extra work vector */
  Vec              wht;              /* extra work vector for hermitian transpose apply */
  STStateType      state;            /* initial -> setup -> with updated matrices */
  PetscObjectState *Astate;          /* matrix state (to identify the original matrices) */
  Mat              *T;               /* matrices resulting from transformation */
  Mat              Op;               /* shell matrix for operator = alpha*D*inv(P)*M*inv(D) */
  PetscBool        opseized;         /* whether Op has been seized by user */
  PetscBool        opready;          /* whether Op is up-to-date or need be computed  */
  Mat              P;                /* matrix from which preconditioner is built */
  Mat              M;                /* matrix corresponding to the non-inverted part of the operator */
  PetscBool        sigma_set;        /* whether the user provided the shift or not */
  PetscBool        asymm;            /* the user matrices are all symmetric */
  PetscBool        aherm;            /* the user matrices are all hermitian */
  void             *data;
};

/*
    Macros to test valid ST arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define STCheckMatrices(h,arg) do {} while (0)
#define STCheckNotSeized(h,arg) do {} while (0)

#else

#define STCheckMatrices(h,arg) \
  do { \
    if (!(h)->A) SETERRQ1(PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"ST matrices have not been set: Parameter #%d",arg); \
  } while (0)
#define STCheckNotSeized(h,arg) \
  do { \
    if (h->opseized) SETERRQ1(PetscObjectComm((PetscObject)h),PETSC_ERR_ARG_WRONGSTATE,"Must call STRestoreOperator() first: Parameter #%d",arg); \
  } while (0)

#endif

SLEPC_INTERN PetscErrorCode STGetBilinearForm_Default(ST,Mat*);
SLEPC_INTERN PetscErrorCode STCheckNullSpace_Default(ST,BV);
SLEPC_INTERN PetscErrorCode STMatShellCreate(ST,PetscScalar,PetscInt,PetscInt*,PetscScalar*,Mat*);
SLEPC_INTERN PetscErrorCode STMatShellShift(Mat,PetscScalar);
SLEPC_INTERN PetscErrorCode STCheckFactorPackage(ST);
SLEPC_INTERN PetscErrorCode STMatMAXPY_Private(ST,PetscScalar,PetscScalar,PetscInt,PetscScalar*,PetscBool,Mat*);
SLEPC_INTERN PetscErrorCode STCoeffs_Monomial(ST,PetscScalar*);
SLEPC_INTERN PetscErrorCode STSetDefaultKSP(ST);
SLEPC_INTERN PetscErrorCode STSetDefaultKSP_Default(ST);
SLEPC_INTERN PetscErrorCode STIsInjective_Shell(ST,PetscBool*);
SLEPC_INTERN PetscErrorCode STComputeOperator(ST);
SLEPC_INTERN PetscErrorCode STGetOperator_Private(ST,Mat*);
SLEPC_INTERN PetscErrorCode STApply_Generic(ST,Vec,Vec);
SLEPC_INTERN PetscErrorCode STApplyTranspose_Generic(ST,Vec,Vec);

#endif
