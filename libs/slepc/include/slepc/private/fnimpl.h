/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if !defined(SLEPCFNIMPL_H)
#define SLEPCFNIMPL_H

#include <slepcfn.h>
#include <slepc/private/slepcimpl.h>

SLEPC_EXTERN PetscBool FNRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode FNRegisterAll(void);
SLEPC_EXTERN PetscLogEvent FN_Evaluate;

typedef struct _FNOps *FNOps;

struct _FNOps {
  PetscErrorCode (*evaluatefunction)(FN,PetscScalar,PetscScalar*);
  PetscErrorCode (*evaluatederivative)(FN,PetscScalar,PetscScalar*);
  PetscErrorCode (*evaluatefunctionmat[FN_MAX_SOLVE])(FN,Mat,Mat);
  PetscErrorCode (*evaluatefunctionmatvec[FN_MAX_SOLVE])(FN,Mat,Vec);
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,FN);
  PetscErrorCode (*view)(FN,PetscViewer);
  PetscErrorCode (*duplicate)(FN,MPI_Comm,FN*);
  PetscErrorCode (*destroy)(FN);
};

#define FN_MAX_W 6

struct _p_FN {
  PETSCHEADER(struct _FNOps);
  /*------------------------- User parameters --------------------------*/
  PetscScalar    alpha;          /* inner scaling (argument) */
  PetscScalar    beta;           /* outer scaling (result) */
  PetscInt       method;         /* the method to compute matrix functions */
  FNParallelType pmode;          /* parallel mode (redundant or synchronized) */

  /*---------------------- Cached data and workspace -------------------*/
  Mat            W[FN_MAX_W];    /* workspace matrices */
  PetscInt       nw;             /* number of allocated W matrices */
  PetscInt       cw;             /* current W matrix */
  void           *data;
};

/*
  FN_AllocateWorkMat - Allocate a work Mat of the same dimension of A and copy
  its contents. The work matrix is returned in M and should be freed with
  FN_FreeWorkMat().
*/
PETSC_STATIC_INLINE PetscErrorCode FN_AllocateWorkMat(FN fn,Mat A,Mat *M)
{
  PetscErrorCode ierr;
  PetscInt       n,na;
  PetscBool      create=PETSC_FALSE;

  PetscFunctionBegin;
  *M = NULL;
  if (fn->cw==FN_MAX_W) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Too many requested work matrices %D",fn->cw);
  if (fn->nw<=fn->cw) {
    create=PETSC_TRUE;
    fn->nw++;
  } else {
    ierr = MatGetSize(fn->W[fn->cw],&n,NULL);CHKERRQ(ierr);
    ierr = MatGetSize(A,&na,NULL);CHKERRQ(ierr);
    if (n!=na) {
      ierr = MatDestroy(&fn->W[fn->cw]);CHKERRQ(ierr);
      create=PETSC_TRUE;
    }
  }
  if (create) {
    ierr = MatDuplicate(A,MAT_COPY_VALUES,&fn->W[fn->cw]);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)fn,(PetscObject)fn->W[fn->cw]);CHKERRQ(ierr);
  } else {
    ierr = MatCopy(A,fn->W[fn->cw],SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  }
  *M = fn->W[fn->cw];
  fn->cw++;
  PetscFunctionReturn(0);
}

/*
  FN_FreeWorkMat - Release a work matrix created with FN_AllocateWorkMat().
*/
PETSC_STATIC_INLINE PetscErrorCode FN_FreeWorkMat(FN fn,Mat *M)
{
  PetscFunctionBegin;
  if (!fn->cw) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"There are no work matrices");
  fn->cw--;
  if (fn->W[fn->cw]!=*M) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Work matrices must be freed in the reverse order of their creation");
  *M = NULL;
  PetscFunctionReturn(0);
}

SLEPC_INTERN PetscErrorCode SlepcMatDenseSqrt(PetscBLASInt,PetscScalar*,PetscBLASInt);
SLEPC_INTERN PetscErrorCode SlepcSqrtmSchur(PetscBLASInt,PetscScalar*,PetscBLASInt,PetscBool);
SLEPC_INTERN PetscErrorCode SlepcSqrtmDenmanBeavers(PetscBLASInt,PetscScalar*,PetscBLASInt,PetscBool);
SLEPC_INTERN PetscErrorCode SlepcSqrtmNewtonSchulz(PetscBLASInt,PetscScalar*,PetscBLASInt,PetscBool);
SLEPC_INTERN PetscErrorCode SlepcNormAm(PetscBLASInt,PetscScalar*,PetscInt,PetscScalar*,PetscRandom,PetscReal*);
SLEPC_INTERN PetscErrorCode FNEvaluateFunctionMat_Private(FN,Mat,Mat,PetscBool);
SLEPC_INTERN PetscErrorCode FNEvaluateFunctionMatVec_Private(FN,Mat,Vec,PetscBool);
SLEPC_INTERN PetscErrorCode FNEvaluateFunctionMat_Exp_Higham(FN,Mat,Mat); /* used in FNPHI */

#endif
