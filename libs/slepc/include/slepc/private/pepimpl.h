/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if !defined(SLEPCPEPIMPL_H)
#define SLEPCPEPIMPL_H

#include <slepcpep.h>
#include <slepc/private/slepcimpl.h>

SLEPC_EXTERN PetscBool PEPRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode PEPRegisterAll(void);
SLEPC_EXTERN PetscLogEvent PEP_SetUp,PEP_Solve,PEP_Refine;

typedef struct _PEPOps *PEPOps;

struct _PEPOps {
  PetscErrorCode (*solve)(PEP);
  PetscErrorCode (*setup)(PEP);
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,PEP);
  PetscErrorCode (*publishoptions)(PEP);
  PetscErrorCode (*destroy)(PEP);
  PetscErrorCode (*reset)(PEP);
  PetscErrorCode (*view)(PEP,PetscViewer);
  PetscErrorCode (*backtransform)(PEP);
  PetscErrorCode (*computevectors)(PEP);
  PetscErrorCode (*extractvectors)(PEP);
  PetscErrorCode (*setdefaultst)(PEP);
};

/*
     Maximum number of monitors you can run with a single PEP
*/
#define MAXPEPMONITORS 5

typedef enum { PEP_STATE_INITIAL,
               PEP_STATE_SETUP,
               PEP_STATE_SOLVED,
               PEP_STATE_EIGENVECTORS } PEPStateType;

/*
   Defines the PEP data structure.
*/
struct _p_PEP {
  PETSCHEADER(struct _PEPOps);
  /*------------------------- User parameters ---------------------------*/
  PetscInt       max_it;           /* maximum number of iterations */
  PetscInt       nev;              /* number of eigenvalues to compute */
  PetscInt       ncv;              /* number of basis vectors */
  PetscInt       mpd;              /* maximum dimension of projected problem */
  PetscInt       nini;             /* number of initial vectors (negative means not copied yet) */
  PetscScalar    target;           /* target value */
  PetscReal      tol;              /* tolerance */
  PEPConv        conv;             /* convergence test */
  PEPStop        stop;             /* stopping test */
  PEPWhich       which;            /* which part of the spectrum to be sought */
  PetscReal      inta,intb;        /* interval [a,b] for spectrum slicing */
  PEPBasis       basis;            /* polynomial basis used to represent the problem */
  PEPProblemType problem_type;     /* which kind of problem to be solved */
  PEPScale       scale;            /* scaling strategy to be used */
  PetscReal      sfactor,dsfactor; /* scaling factors */
  PetscInt       sits;             /* number of iterations of the scaling method */
  PetscReal      slambda;          /* norm eigenvalue approximation for scaling */
  PEPRefine      refine;           /* type of refinement to be applied after solve */
  PetscInt       npart;            /* number of partitions of the communicator */
  PetscReal      rtol;             /* tolerance for refinement */
  PetscInt       rits;             /* number of iterations of the refinement method */
  PEPRefineScheme scheme;          /* scheme for solving linear systems within refinement */
  PEPExtract     extract;          /* type of extraction used */
  PetscBool      trackall;         /* whether all the residuals must be computed */

  /*-------------- User-provided functions and contexts -----------------*/
  PetscErrorCode (*converged)(PEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
  PetscErrorCode (*convergeduser)(PEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*);
  PetscErrorCode (*convergeddestroy)(void*);
  PetscErrorCode (*stopping)(PEP,PetscInt,PetscInt,PetscInt,PetscInt,PEPConvergedReason*,void*);
  PetscErrorCode (*stoppinguser)(PEP,PetscInt,PetscInt,PetscInt,PetscInt,PEPConvergedReason*,void*);
  PetscErrorCode (*stoppingdestroy)(void*);
  void           *convergedctx;
  void           *stoppingctx;
  PetscErrorCode (*monitor[MAXPEPMONITORS])(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*);
  PetscErrorCode (*monitordestroy[MAXPEPMONITORS])(void**);
  void           *monitorcontext[MAXPEPMONITORS];
  PetscInt        numbermonitors;

  /*----------------- Child objects and working data -------------------*/
  ST             st;               /* spectral transformation object */
  DS             ds;               /* direct solver object */
  BV             V;                /* set of basis vectors and computed eigenvectors */
  RG             rg;               /* optional region for filtering */
  SlepcSC        sc;               /* sorting criterion data */
  Mat            *A;               /* coefficient matrices of the polynomial */
  PetscInt       nmat;             /* number of matrices */
  Vec            Dl,Dr;            /* diagonal matrices for balancing */
  Vec            *IS;              /* references to user-provided initial space */
  PetscScalar    *eigr,*eigi;      /* real and imaginary parts of eigenvalues */
  PetscReal      *errest;          /* error estimates */
  PetscInt       *perm;            /* permutation for eigenvalue ordering */
  PetscReal      *pbc;             /* coefficients defining the polynomial basis */
  PetscScalar    *solvematcoeffs;  /* coefficients to compute the matrix to be inverted */
  PetscInt       nwork;            /* number of work vectors */
  Vec            *work;            /* work vectors */
  KSP            refineksp;        /* ksp used in refinement */
  PetscSubcomm   refinesubc;       /* context for sub-communicators */
  void           *data;            /* placeholder for solver-specific stuff */

  /* ----------------------- Status variables --------------------------*/
  PEPStateType   state;            /* initial -> setup -> solved -> eigenvectors */
  PetscInt       nconv;            /* number of converged eigenvalues */
  PetscInt       its;              /* number of iterations so far computed */
  PetscInt       n,nloc;           /* problem dimensions (global, local) */
  PetscReal      *nrma;            /* computed matrix norms */
  PetscReal      nrml[2];          /* computed matrix norms for the linearization */
  PetscBool      sfactor_set;      /* flag to indicate the user gave sfactor */
  PetscBool      lineariz;         /* current solver is based on linearization */
  PEPConvergedReason reason;
};

/*
    Macros to test valid PEP arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define PEPCheckSolved(h,arg) do {} while (0)

#else

#define PEPCheckSolved(h,arg) \
  do { \
    if ((h)->state<PEP_STATE_SOLVED) SETERRQ1(PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"Must call PEPSolve() first: Parameter #%d",arg); \
  } while (0)

#endif

SLEPC_INTERN PetscErrorCode PEPSetDimensions_Default(PEP,PetscInt,PetscInt*,PetscInt*);
SLEPC_INTERN PetscErrorCode PEPExtractVectors(PEP);
SLEPC_INTERN PetscErrorCode PEPBackTransform_Default(PEP);
SLEPC_INTERN PetscErrorCode PEPComputeVectors(PEP);
SLEPC_INTERN PetscErrorCode PEPComputeVectors_Default(PEP);
SLEPC_INTERN PetscErrorCode PEPComputeVectors_Indefinite(PEP);
SLEPC_INTERN PetscErrorCode PEPComputeResidualNorm_Private(PEP,PetscScalar,PetscScalar,Vec,Vec,Vec*,PetscReal*);
SLEPC_INTERN PetscErrorCode PEPKrylovConvergence(PEP,PetscBool,PetscInt,PetscInt,PetscReal,PetscInt*);
SLEPC_INTERN PetscErrorCode PEPComputeScaleFactor(PEP);
SLEPC_INTERN PetscErrorCode PEPBuildDiagonalScaling(PEP);
SLEPC_INTERN PetscErrorCode PEPBasisCoefficients(PEP,PetscReal*);
SLEPC_INTERN PetscErrorCode PEPEvaluateBasis(PEP,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode PEPEvaluateBasisDerivative(PEP,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode PEPEvaluateBasisMat(PEP,PetscInt,PetscScalar*,PetscInt,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode PEPNewtonRefinement_TOAR(PEP,PetscScalar,PetscInt*,PetscReal*,PetscInt,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode PEPNewtonRefinementSimple(PEP,PetscInt*,PetscReal,PetscInt);
SLEPC_INTERN PetscErrorCode PEPSetDefaultST(PEP);
SLEPC_INTERN PetscErrorCode PEPSetDefaultST_Transform(PEP);

#endif
