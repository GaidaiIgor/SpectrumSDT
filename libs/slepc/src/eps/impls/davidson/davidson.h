/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Method: General Davidson Method (includes GD and JD)

   References:
    - Ernest R. Davidson. Super-matrix methods. Computer Physics Communications,
      53:49-60, May 1989.
*/

#include <slepc/private/epsimpl.h>
#include <slepc/private/vecimplslepc.h>

typedef PetscInt MatType_t;
#define DVD_MAT_HERMITIAN (1<<1)
#define DVD_MAT_NEG_DEF (1<<2)
#define DVD_MAT_POS_DEF (1<<3)
#define DVD_MAT_SINGULAR (1<<4)
#define DVD_MAT_COMPLEX (1<<5)
#define DVD_MAT_IMPLICIT (1<<6)
#define DVD_MAT_IDENTITY (1<<7)
#define DVD_MAT_DIAG (1<<8)
#define DVD_MAT_TRIANG (1<<9)
#define DVD_MAT_UTRIANG (1<<9)
#define DVD_MAT_LTRIANG (1<<10)
#define DVD_MAT_UNITARY (1<<11)

typedef PetscInt EPType_t;
#define DVD_EP_STD (1<<1)
#define DVD_EP_HERMITIAN (1<<2)
#define DVD_EP_INDEFINITE (1<<3)

#define DVD_IS(T,P) ((T) & (P))
#define DVD_ISNOT(T,P) (((T) & (P)) ^ (P))

struct _dvdDashboard;
typedef PetscErrorCode (*dvdCallback)(struct _dvdDashboard*);
typedef struct _dvdFunctionList {
  dvdCallback f;
  struct _dvdFunctionList *next;
} dvdFunctionList;

typedef enum {
  DVD_HARM_NONE,
  DVD_HARM_RR,
  DVD_HARM_RRR,
  DVD_HARM_REIGS,
  DVD_HARM_LEIGS
} HarmType_t;

typedef enum {
  DVD_INITV_CLASSIC,
  DVD_INITV_KRYLOV
} InitType_t;

/*
   Dashboard struct: contains the methods that will be employed and the tuning
   options.
*/

typedef struct _dvdDashboard {
  /**** Function steps ****/
  /* Initialize V */
  PetscErrorCode (*initV)(struct _dvdDashboard*);
  void *initV_data;

  /* Find the approximate eigenpairs from V */
  PetscErrorCode (*calcPairs)(struct _dvdDashboard*);
  void *calcPairs_data;

  /* Eigenpair test for convergence */
  PetscBool (*testConv)(struct _dvdDashboard*,PetscScalar,PetscScalar,PetscReal,PetscReal*);
  void *testConv_data;

  /* Improve the selected eigenpairs */
  PetscErrorCode (*improveX)(struct _dvdDashboard*,PetscInt,PetscInt,PetscInt*);
  void *improveX_data;

  /* Check for restarting */
  PetscErrorCode (*isRestarting)(struct _dvdDashboard*,PetscBool*);
  void *isRestarting_data;

  /* Perform restarting */
  PetscErrorCode (*restartV)(struct _dvdDashboard*);
  void *restartV_data;

  /* Update V */
  PetscErrorCode (*updateV)(struct _dvdDashboard*);
  void *updateV_data;

  /**** Problem specification ****/
  Mat         A,B;            /* problem matrices */
  MatType_t   sA,sB;          /* matrix specifications */
  EPType_t    sEP;            /* problem specifications */
  PetscInt    nev;            /* number of eigenpairs */
  EPSWhich    which;          /* spectrum selection */
  PetscBool   withTarget;     /* if there is a target */
  PetscScalar target[2];      /* target value */
  PetscReal   tol;            /* tolerance */
  PetscBool   correctXnorm;   /* if true, norm of X are computed */

  /**** Subspaces specification ****/
  PetscInt nconv;             /* number of converged eigenpairs */
  PetscInt npreconv;          /* number of pairs ready to converge */

  BV       W;                 /* left basis for harmonic case */
  BV       AX;                /* A*V */
  BV       BX;                /* B*V */
  PetscInt size_D;            /* active vectors */
  PetscInt max_size_proj;     /* max size projected problem */
  PetscInt max_size_P;        /* max unconverged vectors in the projector */
  PetscInt bs;                /* max vectors that expands the subspace every iteration */
  EPS      eps;               /* connection to SLEPc */

  /**** Auxiliary space ****/
  VecPool auxV;               /* auxiliary vectors */
  BV      auxBV;              /* auxiliary vectors */

  /**** Eigenvalues and errors ****/
  PetscScalar *ceigr,*ceigi;  /* converged eigenvalues */
  PetscScalar *eigr,*eigi;    /* current eigenvalues */
  PetscReal   *nR;            /* residual norm */
  PetscReal   *real_nR;       /* original nR */
  PetscReal   *nX;            /* X norm */
  PetscReal   *real_nX;       /* original nX */
  PetscReal   *errest;        /* relative error eigenpairs */
  PetscReal   *nBds;          /* B-norms of projected problem  */

  /**** Shared function and variables ****/
  PetscErrorCode (*e_Vchanged)(struct _dvdDashboard*,PetscInt,PetscInt,PetscInt,PetscInt);
  void *e_Vchanged_data;
  PetscErrorCode (*calcpairs_residual)(struct _dvdDashboard*,PetscInt,PetscInt);
  PetscErrorCode (*calcpairs_selectPairs)(struct _dvdDashboard*,PetscInt);
  void *calcpairs_residual_data;
  PetscErrorCode (*improvex_precond)(struct _dvdDashboard*,PetscInt,Vec,Vec);
  void *improvex_precond_data;
  PetscErrorCode (*improvex_jd_proj_uv)(struct _dvdDashboard*,PetscInt,PetscInt,Vec*,Vec*,Vec*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt);
  PetscErrorCode (*improvex_jd_lit)(struct _dvdDashboard*,PetscInt,PetscScalar*,PetscScalar*,PetscInt*,PetscReal*);
  PetscErrorCode (*calcpairs_W)(struct _dvdDashboard*);
  void *calcpairs_W_data;
  PetscErrorCode (*calcpairs_proj_trans)(struct _dvdDashboard*);
  PetscErrorCode (*calcpairs_eigs_trans)(struct _dvdDashboard*);
  PetscErrorCode (*calcpairs_eig_backtrans)(struct _dvdDashboard*,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*);
  PetscErrorCode (*calcpairs_proj_res)(struct _dvdDashboard*,PetscInt,PetscInt,Vec*);
  PetscErrorCode (*preTestConv)(struct _dvdDashboard*,PetscInt,PetscInt,PetscInt,PetscInt*);
  PetscErrorCode (*e_newIteration)(struct _dvdDashboard*);
  void *e_newIteration_data;

  dvdFunctionList *startList;  /* starting list */
  dvdFunctionList *endList;    /* ending list */
  dvdFunctionList *destroyList;/* destructor list */

  Mat       H,G;               /* projected problem matrices */
  Mat       auxM;              /* auxiliary dense matrix */
  PetscInt  size_MT;           /* rows in MT */

  PetscInt  V_tra_s;
  PetscInt  V_tra_e;       /* cX <- [cX V*MT(0:V_tra_s-1)], V <- V*MT(V_tra_s:V_tra_e) */
  PetscInt  V_new_s;
  PetscInt  V_new_e;           /* added to V the columns V_new_s:V_new_e */
  PetscBool W_shift;           /* if true W is shifted when vectors converge */
} dvdDashboard;

typedef struct {
  /*------------------------- User parameters ---------------------------*/
  PetscInt  blocksize;     /* block size */
  PetscInt  initialsize;   /* initial size of V */
  PetscInt  minv;          /* size of V after restarting */
  PetscInt  plusk;         /* keep plusk eigenvectors from the last iteration */
  PetscBool ipB;           /* true if B-ortho is used */
  PetscReal fix;           /* the fix parameter */
  PetscBool krylovstart;   /* true if the starting subspace is a Krylov basis */
  PetscBool dynamic;       /* true if dynamic stopping criterion is used */
  PetscBool doubleexp;     /* double expansion in GD (GD2) */

  /*----------------- Child objects and working data -------------------*/
  dvdDashboard ddb;
} EPS_DAVIDSON;

PETSC_STATIC_INLINE PetscErrorCode EPSDavidsonFLAdd(dvdFunctionList **fl,dvdCallback f)
{
  PetscErrorCode ierr;
  dvdFunctionList *l;

  PetscFunctionBegin;
  ierr = PetscNew(&l);CHKERRQ(ierr);
  l->f = f;
  l->next = *fl;
  *fl = l;
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode EPSDavidsonFLCall(dvdFunctionList *fl,dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdFunctionList *l;

  PetscFunctionBegin;
  for (l=fl;l;l=l->next) { ierr = (l->f)(d);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode EPSDavidsonFLDestroy(dvdFunctionList **fl)
{
  PetscErrorCode  ierr;
  dvdFunctionList *l,*l0;

  PetscFunctionBegin;
  for (l=*fl;l;l=l0) {
    l0 = l->next;
    ierr = PetscFree(l);CHKERRQ(ierr);
  }
  *fl = NULL;
  PetscFunctionReturn(0);
}

/*
  The blackboard configuration structure: saves information about the memory
  and other requirements.

  The starting memory structure:

  V           W?          AV          BV?          tKZ
  |-----------|-----------|-----------|------------|-------|
  nev+mpd     nev+mpd     scP+mpd     nev?+mpd     sP+scP
              scP+mpd                 scP+mpd

  The final memory structure considering W_shift:

  cX  V       cY?  W?     cAV AV      BcX? BV?     KZ  tKZ
  |---|-------|----|------|---|-------|----|-------|---|---|
  nev mpd     nev  mpd    scP mpd     nev  mpd     scP sP    <- shift
              scP                     scP                    <- !shift
*/
typedef struct {
  PetscInt max_size_V;         /* max size of the searching subspace (mpd) */
  PetscInt max_size_X;         /* max size of X (bs) */
  PetscInt size_V;             /* real size of V (nev+size_P+mpd) */
  PetscInt max_size_oldX;      /* max size of oldX */
  PetscInt max_nev;            /* max number of converged pairs */
  PetscInt max_size_P;         /* number of computed vectors for the projector */
  PetscInt max_size_cP;        /* number of converged vectors in the projectors */
  PetscInt max_size_proj;      /* max size projected problem */
  PetscInt max_size_cX_proj;   /* max converged vectors in the projected problem */
  PetscInt state;              /* method states:
                                   0: preconfiguring
                                   1: configuring
                                   2: running */
} dvdBlackboard;

#define DVD_STATE_PRECONF 0
#define DVD_STATE_CONF 1
#define DVD_STATE_RUN 2

/* Prototypes of non-static auxiliary functions */
SLEPC_INTERN PetscErrorCode dvd_calcpairs_qz(dvdDashboard*,dvdBlackboard*,PetscBool,PetscBool);
SLEPC_INTERN PetscErrorCode dvd_improvex_gd2(dvdDashboard*,dvdBlackboard*,KSP,PetscInt);
SLEPC_INTERN PetscErrorCode dvd_improvex_jd(dvdDashboard*,dvdBlackboard*,KSP,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode dvd_improvex_jd_proj_uv(dvdDashboard*,dvdBlackboard*);
SLEPC_INTERN PetscErrorCode dvd_improvex_jd_lit_const(dvdDashboard*,dvdBlackboard*,PetscInt,PetscReal,PetscReal);
SLEPC_INTERN PetscErrorCode dvd_improvex_compute_X(dvdDashboard*,PetscInt,PetscInt,Vec*,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode dvd_initV(dvdDashboard*,dvdBlackboard*,PetscInt,PetscInt,PetscBool);
SLEPC_INTERN PetscErrorCode dvd_orthV(BV,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode dvd_schm_basic_preconf(dvdDashboard*,dvdBlackboard*,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,HarmType_t,KSP,InitType_t,PetscBool,PetscBool,PetscBool);
SLEPC_INTERN PetscErrorCode dvd_schm_basic_conf(dvdDashboard*,dvdBlackboard*,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,HarmType_t,PetscBool,PetscScalar,KSP,PetscReal,InitType_t,PetscBool,PetscBool,PetscBool,PetscBool);
SLEPC_INTERN PetscErrorCode dvd_testconv_slepc(dvdDashboard*,dvdBlackboard*);
SLEPC_INTERN PetscErrorCode dvd_managementV_basic(dvdDashboard*,dvdBlackboard*,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool,PetscBool);
SLEPC_INTERN PetscErrorCode dvd_static_precond_PC(dvdDashboard*,dvdBlackboard*,PC);
SLEPC_INTERN PetscErrorCode dvd_harm_updateproj(dvdDashboard*);
SLEPC_INTERN PetscErrorCode dvd_harm_conf(dvdDashboard*,dvdBlackboard*,HarmType_t,PetscBool,PetscScalar);

/* Internal interface routines */
SLEPC_INTERN PetscErrorCode EPSReset_XD(EPS);
SLEPC_INTERN PetscErrorCode EPSSetUp_XD(EPS);
SLEPC_INTERN PetscErrorCode EPSSolve_XD(EPS);
SLEPC_INTERN PetscErrorCode EPSComputeVectors_XD(EPS);
SLEPC_INTERN PetscErrorCode EPSXDSetKrylovStart_XD(EPS,PetscBool);
SLEPC_INTERN PetscErrorCode EPSXDGetKrylovStart_XD(EPS,PetscBool*);
SLEPC_INTERN PetscErrorCode EPSXDSetBlockSize_XD(EPS,PetscInt);
SLEPC_INTERN PetscErrorCode EPSXDGetBlockSize_XD(EPS,PetscInt*);
SLEPC_INTERN PetscErrorCode EPSXDSetRestart_XD(EPS,PetscInt,PetscInt);
SLEPC_INTERN PetscErrorCode EPSXDGetRestart_XD(EPS,PetscInt*,PetscInt*);
SLEPC_INTERN PetscErrorCode EPSXDGetInitialSize_XD(EPS,PetscInt*);
SLEPC_INTERN PetscErrorCode EPSXDSetInitialSize_XD(EPS,PetscInt);
SLEPC_INTERN PetscErrorCode EPSXDSetBOrth_XD(EPS,PetscBool);
SLEPC_INTERN PetscErrorCode EPSXDGetBOrth_XD(EPS,PetscBool*);
SLEPC_INTERN PetscErrorCode EPSJDGetFix_JD(EPS,PetscReal*);
SLEPC_INTERN PetscErrorCode EPSJDGetConstCorrectionTol_JD(EPS,PetscBool*);

