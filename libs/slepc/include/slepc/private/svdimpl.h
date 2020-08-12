/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if !defined(SLEPCSVDIMPL_H)
#define SLEPCSVDIMPL_H

#include <slepcsvd.h>
#include <slepc/private/slepcimpl.h>

SLEPC_EXTERN PetscBool SVDRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode SVDRegisterAll(void);
SLEPC_EXTERN PetscLogEvent SVD_SetUp,SVD_Solve;

typedef struct _SVDOps *SVDOps;

struct _SVDOps {
  PetscErrorCode (*solve)(SVD);
  PetscErrorCode (*setup)(SVD);
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,SVD);
  PetscErrorCode (*publishoptions)(SVD);
  PetscErrorCode (*destroy)(SVD);
  PetscErrorCode (*reset)(SVD);
  PetscErrorCode (*view)(SVD,PetscViewer);
};

/*
     Maximum number of monitors you can run with a single SVD
*/
#define MAXSVDMONITORS 5

typedef enum { SVD_STATE_INITIAL,
               SVD_STATE_SETUP,
               SVD_STATE_SOLVED,
               SVD_STATE_VECTORS } SVDStateType;

/*
   Defines the SVD data structure.
*/
struct _p_SVD {
  PETSCHEADER(struct _SVDOps);
  /*------------------------- User parameters ---------------------------*/
  Mat            OP;               /* problem matrix */
  PetscInt       max_it;           /* max iterations */
  PetscInt       nsv;              /* number of requested values */
  PetscInt       ncv;              /* basis size */
  PetscInt       mpd;              /* maximum dimension of projected problem */
  PetscInt       nini,ninil;       /* number of initial vecs (negative means not copied yet) */
  PetscReal      tol;              /* tolerance */
  SVDConv        conv;             /* convergence test */
  SVDStop        stop;             /* stopping test */
  SVDWhich       which;            /* which singular values are computed */
  PetscBool      impltrans;        /* implicit transpose mode */
  PetscBool      trackall;         /* whether all the residuals must be computed */

  /*-------------- User-provided functions and contexts -----------------*/
  PetscErrorCode (*converged)(SVD,PetscReal,PetscReal,PetscReal*,void*);
  PetscErrorCode (*convergeduser)(SVD,PetscReal,PetscReal,PetscReal*,void*);
  PetscErrorCode (*convergeddestroy)(void*);
  PetscErrorCode (*stopping)(SVD,PetscInt,PetscInt,PetscInt,PetscInt,SVDConvergedReason*,void*);
  PetscErrorCode (*stoppinguser)(SVD,PetscInt,PetscInt,PetscInt,PetscInt,SVDConvergedReason*,void*);
  PetscErrorCode (*stoppingdestroy)(void*);
  void           *convergedctx;
  void           *stoppingctx;
  PetscErrorCode (*monitor[MAXSVDMONITORS])(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*);
  PetscErrorCode (*monitordestroy[MAXSVDMONITORS])(void**);
  void           *monitorcontext[MAXSVDMONITORS];
  PetscInt       numbermonitors;

  /*----------------- Child objects and working data -------------------*/
  DS             ds;               /* direct solver object */
  BV             U,V;              /* left and right singular vectors */
  SlepcSC        sc;               /* sorting criterion data */
  Mat            A;                /* problem matrix (m>n) */
  Mat            AT;               /* transposed matrix */
  Vec            *IS,*ISL;         /* placeholder for references to user initial space */
  PetscReal      *sigma;           /* singular values */
  PetscInt       *perm;            /* permutation for singular value ordering */
  PetscReal      *errest;          /* error estimates */
  void           *data;            /* placeholder for solver-specific stuff */

  /* ----------------------- Status variables -------------------------- */
  SVDStateType   state;            /* initial -> setup -> solved -> vectors */
  PetscInt       nconv;            /* number of converged values */
  PetscInt       its;              /* iteration counter */
  PetscBool      leftbasis;        /* if U is filled by the solver */
  SVDConvergedReason reason;
};

/*
    Macros to test valid SVD arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define SVDCheckSolved(h,arg) do {} while (0)

#else

#define SVDCheckSolved(h,arg) \
  do { \
    if ((h)->state<SVD_STATE_SOLVED) SETERRQ1(PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"Must call SVDSolve() first: Parameter #%d",arg); \
  } while (0)

#endif

PETSC_STATIC_INLINE PetscErrorCode SVDMatMult(SVD svd,PetscBool trans,Vec x,Vec y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (trans) {
    if (svd->AT) {
      ierr = MatMult(svd->AT,x,y);CHKERRQ(ierr);
    } else {
#if defined(PETSC_USE_COMPLEX)
      ierr = MatMultHermitianTranspose(svd->A,x,y);CHKERRQ(ierr);
#else
      ierr = MatMultTranspose(svd->A,x,y);CHKERRQ(ierr);
#endif
    }
  } else {
    if (svd->A) {
      ierr = MatMult(svd->A,x,y);CHKERRQ(ierr);
    } else {
#if defined(PETSC_USE_COMPLEX)
      ierr = MatMultHermitianTranspose(svd->AT,x,y);CHKERRQ(ierr);
#else
      ierr = MatMultTranspose(svd->AT,x,y);CHKERRQ(ierr);
#endif
    }
  }
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode SVDMatCreateVecs(SVD svd,Vec *x,Vec *y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (svd->A) {
    ierr = MatCreateVecs(svd->A,x,y);CHKERRQ(ierr);
  } else {
    ierr = MatCreateVecs(svd->AT,y,x);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode SVDMatCreateVecsEmpty(SVD svd,Vec *x,Vec *y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (svd->A) {
    ierr = MatCreateVecsEmpty(svd->A,x,y);CHKERRQ(ierr);
  } else {
    ierr = MatCreateVecsEmpty(svd->AT,y,x);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode SVDMatGetSize(SVD svd,PetscInt *m,PetscInt *n)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (svd->A) {
    ierr = MatGetSize(svd->A,m,n);CHKERRQ(ierr);
  } else {
    ierr = MatGetSize(svd->AT,n,m);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode SVDMatGetLocalSize(SVD svd,PetscInt *m,PetscInt *n)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (svd->A) {
    ierr = MatGetLocalSize(svd->A,m,n);CHKERRQ(ierr);
  } else {
    ierr = MatGetLocalSize(svd->AT,n,m);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

SLEPC_INTERN PetscErrorCode SVDTwoSideLanczos(SVD,PetscReal*,PetscReal*,BV,BV,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode SVDSetDimensions_Default(SVD);
SLEPC_INTERN PetscErrorCode SVDComputeVectors(SVD);

#endif
