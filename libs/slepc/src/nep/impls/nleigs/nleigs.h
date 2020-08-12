/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Private header for NLEIGS
*/

#if !defined(SLEPC_NLEIGS_H)
#define SLEPC_NLEIGS_H


#define  LBPOINTS  100   /* default value of the maximum number of Leja-Bagby points */
#define  NDPOINTS  1e4   /* number of discretization points */

typedef struct {
  BV             V;         /* tensor vector basis for the linearization */
  BV             W;         /* tensor vector basis for the linearization */
  PetscInt       nmat;      /* number of interpolation points */
  PetscScalar    *s,*xi;    /* Leja-Bagby points */
  PetscScalar    *beta;     /* scaling factors */
  Mat            *D;        /* divided difference matrices */
  PetscScalar    *coeffD;   /* coefficients for divided differences in split form */
  PetscInt       nshifts;   /* provided number of shifts */
  PetscScalar    *shifts;   /* user-provided shifts for the Rational Krylov variant */
  PetscInt       nshiftsw;  /* actual number of shifts (1 if Krylov-Schur) */
  PetscReal      ddtol;     /* tolerance for divided difference convergence */
  PetscInt       ddmaxit;   /* maximum number of divided difference terms */
  PetscReal      keep;      /* restart parameter */
  PetscBool      lock;      /* locking/non-locking variant */
  PetscInt       idxrk;     /* index of next shift to use */
  KSP            *ksp;      /* ksp array for storing shift factorizations */
  Vec            vrn;       /* random vector with normally distributed value */
  PetscBool      fullbasis; /* use full Krylov basis instead of TOAR basis */
  EPS            eps;       /* eigensolver used in the full basis variant */
  Mat            A;         /* shell matrix used for the eps in full basis */
  Vec            w[6];      /* work vectors */
  void           *singularitiesctx;
  PetscErrorCode (*computesingularities)(NEP,PetscInt*,PetscScalar*,void*);
} NEP_NLEIGS;

typedef struct {
  PetscInt    nmat,maxnmat;
  PetscScalar *coeff;
  Mat         *A;
  Vec         t;
} NEP_NLEIGS_MATSHELL;

PETSC_STATIC_INLINE PetscErrorCode NEPNLEIGSSetShifts(NEP nep,PetscInt *nshiftsw)
{
  NEP_NLEIGS *ctx = (NEP_NLEIGS*)nep->data;

  PetscFunctionBegin;
  if (!ctx->nshifts) {
    ctx->shifts = &nep->target;
    *nshiftsw = 1;
  } else *nshiftsw = ctx->nshifts;
  PetscFunctionReturn(0);
}

SLEPC_INTERN PetscErrorCode NEPSetUp_NLEIGS_FullBasis(NEP);
SLEPC_INTERN PetscErrorCode NEPNLEIGSSetEPS_NLEIGS(NEP,EPS);
SLEPC_INTERN PetscErrorCode NEPNLEIGSGetEPS_NLEIGS(NEP,EPS*);
SLEPC_INTERN PetscErrorCode NEPNLEIGSBackTransform(PetscObject,PetscInt,PetscScalar*,PetscScalar *vali);
SLEPC_INTERN PetscErrorCode NEPNLEIGSEvalNRTFunct(NEP,PetscInt,PetscScalar,PetscScalar*);
SLEPC_INTERN PetscErrorCode NEPSolve_NLEIGS_FullBasis(NEP);
SLEPC_INTERN PetscErrorCode NEPSolve_NLEIGS(NEP);

#endif
