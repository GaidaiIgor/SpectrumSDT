/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Private header for TOAR and STOAR
*/

#if !defined(SLEPC_PEPKRYLOV_H)
#define SLEPC_PEPKRYLOV_H

SLEPC_INTERN PetscErrorCode PEPExtractVectors_TOAR(PEP);
SLEPC_INTERN PetscErrorCode PEPSTOARSetUpInnerMatrix(PEP,Mat*);
SLEPC_INTERN PetscErrorCode PEPSolve_STOAR(PEP);
SLEPC_INTERN PetscErrorCode PEPSolve_STOAR_QSlice(PEP);
SLEPC_INTERN PetscErrorCode PEPSetUp_STOAR_QSlice(PEP);
SLEPC_INTERN PetscErrorCode PEPReset_STOAR_QSlice(PEP);

typedef struct {
  PetscReal     keep;         /* restart parameter */
  PetscBool     lock;         /* locking/non-locking variant */
  BV            V;            /* tensor basis vectors object for the linearization */
} PEP_TOAR;

/* Structure characterizing a shift in spectrum slicing */
typedef struct _n_shift *PEP_shift;
struct _n_shift {
  PetscReal     value;
  PetscInt      inertia;
  PetscBool     comp[2];      /* Shows completion of subintervals (left and right) */
  PEP_shift     neighb[2];    /* Adjacent shifts */
  PetscInt      index;        /* Index in eig where found values are stored */
  PetscInt      neigs;        /* Number of values found */
  PetscReal     ext[2];       /* Limits for accepted values */
  PetscInt      nsch[2];      /* Number of missing values for each subinterval */
  PetscInt      nconv[2];     /* Converged on each side (accepted or not) */
};

/* Identifies the TOAR vectors for each eigenvector in the global array */
typedef struct {
  PetscInt      nq;
  PetscInt      *q;
} PEP_QInfo;

/* Structure for storing the state of spectrum slicing */
struct _n_SR {
  PetscReal     int0,int1;         /* Extremes of the interval */
  PetscInt      dir;               /* Determines the order of values in eig (+1 incr, -1 decr) */
  PetscBool     hasEnd;            /* Tells whether the interval has an end */
  PetscBool     dirch;             /* Tells if dir has been changed */
  PetscInt      inertia0,inertia1;
  PetscScalar   *back;
  PetscInt      numEigs;           /* Number of eigenvalues in the interval */
  PetscInt      indexEig;
  PEP_shift     sPres;             /* Present shift */
  PEP_shift     *pending;          /* Pending shifts array */
  PetscInt      nPend;             /* Number of pending shifts */
  PetscInt      maxPend;           /* Size of "pending" array */
  PetscInt      *idxDef0,*idxDef1; /* For deflation */
  PetscInt      ndef0,ndef1;       /* Index in deflation arrays */
  PetscInt      nMAXCompl;
  PetscInt      iterCompl;
  PetscInt      itsKs;             /* Krylovschur restarts */
  PetscInt      nleap;
  PEP_shift     s0;                /* Initial shift */
  PEP_shift     sPrev;
  PetscInt      nv;                /* position of restart vector */
  BV            V;                 /* full TOAR basis */
  PetscScalar   *S;                /* TOAR coefficients */
  PetscInt      ld;                /* Leading dimension for each block of S */
  BV            Vnext;             /* temporary working basis during change of shift */
  PetscScalar   *eigr,*eigi;       /* eigenvalues */
  PetscReal     *errest;           /* error estimates */
  PetscInt      *perm;             /* permutation */
  PEP_QInfo     *qinfo;            /* TOAR vectors for each eigenvector */
  PetscInt      intcorr;           /* Global inertia correction */
  PetscInt      type;              /* Global type of eigenvalues in general case */
  PetscInt      symmlost;          /* Counter for symmetry lost */
  Vec           v[3];
  EPS           eps;
  PetscReal     mu;
};
typedef struct _n_SR *PEP_SR;

typedef struct {
  PetscReal     keep;           /* restart parameter */
  PetscBool     lock;           /* locking/non-locking variant */
  BV            V;              /* tensor basis vectors object for the linearization */
  PEP_SR        sr;             /* spectrum slicing context */
  PetscReal     *shifts;        /* array containing global shifts */
  PetscInt      *inertias;      /* array containing global inertias */
  PetscInt      nshifts;        /* elements in the arrays of shifts and inertias */
  PetscInt      nev;            /* number of eigenvalues to compute */
  PetscInt      ncv;            /* number of basis vectors */
  PetscInt      mpd;            /* maximum dimension of projected problem */
  PetscBool     detect;         /* check for zeros during factorizations */
  PetscBool     hyperbolic;     /* hyperbolic problem flag */
  PetscReal     alpha,beta;     /* coefficients defining the linearization */
  PetscBool     checket;        /* check eigenvalue type during spectrum slicing */
} PEP_STOAR;

#endif

