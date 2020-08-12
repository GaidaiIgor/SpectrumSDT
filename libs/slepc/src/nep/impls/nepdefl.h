/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Private header for Deflation in NEP
*/
#if !defined(SLEPC_NEPDEFL_H)
#define SLEPC_NEPDEFL_H

# define MAX_MINIDX 1

typedef struct _n_nep_ext_op *NEP_EXT_OP;
typedef struct _n_nep_def_fun_solve *NEP_DEF_FUN_SOLVE;
typedef struct _n_nep_def_project *NEP_DEF_PROJECT;

/* Structure characterizing a deflation context */
struct _n_nep_ext_op {
  NEP               nep;
  PetscScalar       *H;     /* invariant pair (X,H) */
  BV                X;      /* locked eigenvectors */
  PetscScalar       *bc;    /* polinomial basis roots */
  RG                rg;
  PetscInt          midx;   /* minimality index */
  PetscInt          max_midx;
  PetscInt          szd;    /* maximum size for deflation */
  PetscInt          n;      /* invariant pair size */
  PetscBool         ref;    /* refine with original size */
  Mat               MF;     /* function shell matrix */
  Mat               MJ;     /* Jacobian shell matrix */
  PetscBool         simpU;  /* the way U is computed */
  NEP_DEF_FUN_SOLVE solve;  /* MatSolve context for the operator */
  NEP_DEF_PROJECT   proj;   /* context for the projected eigenproblem */
  /* auxiliary computations */
  BV                W;
  PetscScalar       *Hj;    /* matrix containing the powers of the invariant pair matrix */
  PetscScalar       *XpX;   /* X^*X */
  DS                ds;
  Vec               w;
};

struct _n_nep_def_fun_solve {
  KSP          ksp;   /* */
  PetscBool    sincf;
  Mat          T;
  PetscScalar  theta;
  PetscInt     n;
  PetscScalar  *M;
  PetscScalar  *work;
  Vec          w[2];
  BV           T_1U;
  NEP_EXT_OP   extop;
};

typedef struct {
  NEP          nep;
  Mat          T;
  BV           U;
  PetscScalar  *A;
  PetscScalar  *B;
  PetscScalar  theta;
  PetscInt     n;
  NEP_EXT_OP   extop;
  PetscBool    jacob;
  Vec          w[2];
  PetscScalar  *work;
  PetscScalar  *hfj;
  PetscScalar  *hfjp;
  PetscBool    hfjset;
} NEP_DEF_MATSHELL;

struct _n_nep_def_project {
  Mat          *V1pApX;
  Mat          XpV1;
  PetscScalar  *V2;
  Vec          w;
  BV           V1;
  PetscInt     dim;
  PetscScalar  *work;
  PetscInt     lwork;
  NEP_EXT_OP   extop;
};

#if 0
typedef struct {
  PC          pc;      /* basic preconditioner */
  PetscScalar *M;
  PetscScalar *ps;
  PetscInt    ld;
  Vec         *work;
  BV          X;
  PetscInt    n;
} NEP_DEF_PCSHELL;
#endif
#endif

SLEPC_INTERN PetscErrorCode NEPDeflationCopyToExtendedVec(NEP_EXT_OP,Vec,PetscScalar*,Vec,PetscBool);
SLEPC_INTERN PetscErrorCode NEPDeflationReset(NEP_EXT_OP);
SLEPC_INTERN PetscErrorCode NEPDeflationInitialize(NEP,BV,KSP,PetscBool,PetscInt,NEP_EXT_OP*);
SLEPC_INTERN PetscErrorCode NEPDeflationCreateVec(NEP_EXT_OP,Vec*);
SLEPC_INTERN PetscErrorCode NEPDeflationComputeFunction(NEP_EXT_OP,PetscScalar,Mat*);
SLEPC_INTERN PetscErrorCode NEPDeflationComputeJacobian(NEP_EXT_OP,PetscScalar,Mat*);
SLEPC_INTERN PetscErrorCode NEPDeflationSolveSetUp(NEP_EXT_OP,PetscScalar);
SLEPC_INTERN PetscErrorCode NEPDeflationFunctionSolve(NEP_EXT_OP,Vec,Vec);
SLEPC_INTERN PetscErrorCode NEPDeflationGetInvariantPair(NEP_EXT_OP,BV*,Mat*);
SLEPC_INTERN PetscErrorCode NEPDeflationLocking(NEP_EXT_OP,Vec,PetscScalar);
SLEPC_INTERN PetscErrorCode NEPDeflationSetRandomVec(NEP_EXT_OP,Vec);
SLEPC_INTERN PetscErrorCode NEPDeflationProjectOperator(NEP_EXT_OP,BV,DS,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode NEPDeflationCreateBV(NEP_EXT_OP,PetscInt,BV*);
SLEPC_INTERN PetscErrorCode NEPDeflationSetRefine(NEP_EXT_OP,PetscBool);
SLEPC_INTERN PetscErrorCode NEPDeflationExtractEigenpair(NEP_EXT_OP,PetscInt,Vec,PetscScalar,DS);
