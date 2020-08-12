/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "lanczos"

   Method: Explicitly Restarted Symmetric/Hermitian Lanczos

   Algorithm:

       Lanczos method for symmetric (Hermitian) problems, with explicit
       restart and deflation. Several reorthogonalization strategies can
       be selected.

   References:

       [1] "Lanczos Methods in SLEPc", SLEPc Technical Report STR-5,
           available at https://slepc.upv.es.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include <slepcblaslapack.h>

typedef struct {
  EPSLanczosReorthogType reorthog;      /* user-provided reorthogonalization parameter */
  PetscInt               allocsize;     /* number of columns of work BV's allocated at setup */
  BV                     AV;            /* work BV used in selective reorthogonalization */
} EPS_LANCZOS;

PetscErrorCode EPSSetUp_Lanczos(EPS eps)
{
  EPS_LANCZOS        *lanczos = (EPS_LANCZOS*)eps->data;
  BVOrthogRefineType refine;
  BVOrthogBlockType  btype;
  PetscReal          eta;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  ierr = EPSSetDimensions_Default(eps,eps->nev,&eps->ncv,&eps->mpd);CHKERRQ(ierr);
  if (eps->ncv>eps->nev+eps->mpd) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv must not be larger than nev+mpd");
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(100,2*eps->n/eps->ncv);
  if (!eps->which) { ierr = EPSSetWhichEigenpairs_Default(eps);CHKERRQ(ierr); }
  switch (eps->which) {
    case EPS_LARGEST_IMAGINARY:
    case EPS_SMALLEST_IMAGINARY:
    case EPS_TARGET_IMAGINARY:
    case EPS_ALL:
      SETERRQ(PetscObjectComm((PetscObject)eps),1,"Wrong value of eps->which");
    default: ; /* default case to remove warning */
  }
  if (!eps->extraction) {
    ierr = EPSSetExtraction(eps,EPS_RITZ);CHKERRQ(ierr);
  } else if (eps->extraction!=EPS_RITZ) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");

  ierr = EPSAllocateSolution(eps,1);CHKERRQ(ierr);
  ierr = EPS_SetInnerProduct(eps);CHKERRQ(ierr);
  if (lanczos->reorthog != EPS_LANCZOS_REORTHOG_FULL) {
    ierr = BVGetOrthogonalization(eps->V,NULL,&refine,&eta,&btype);CHKERRQ(ierr);
    ierr = BVSetOrthogonalization(eps->V,BV_ORTHOG_MGS,refine,eta,btype);CHKERRQ(ierr);
    ierr = PetscInfo(eps,"Switching to MGS orthogonalization\n");CHKERRQ(ierr);
  }
  if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_SELECTIVE) {
    if (!lanczos->allocsize) {
      ierr = BVDuplicate(eps->V,&lanczos->AV);CHKERRQ(ierr);
      ierr = BVGetSizes(lanczos->AV,NULL,NULL,&lanczos->allocsize);CHKERRQ(ierr);
    } else { /* make sure V and AV have the same size */
      ierr = BVGetSizes(eps->V,NULL,NULL,&lanczos->allocsize);CHKERRQ(ierr);
      ierr = BVResize(lanczos->AV,lanczos->allocsize,PETSC_FALSE);CHKERRQ(ierr);
    }
  }

  ierr = DSSetType(eps->ds,DSHEP);CHKERRQ(ierr);
  ierr = DSSetCompact(eps->ds,PETSC_TRUE);CHKERRQ(ierr);
  ierr = DSAllocate(eps->ds,eps->ncv+1);CHKERRQ(ierr);
  if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_LOCAL) {
    ierr = EPSSetWorkVecs(eps,1);CHKERRQ(ierr);
  }

  if (!eps->ishermitian) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Requested method is only available for Hermitian problems");
  if (eps->isgeneralized && eps->ishermitian && !eps->ispositive) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Requested method does not work for indefinite problems");
  PetscFunctionReturn(0);
}

/*
   EPSLocalLanczos - Local reorthogonalization.

   This is the simplest variant. At each Lanczos step, the corresponding Lanczos vector
   is orthogonalized with respect to the two previous Lanczos vectors, according to
   the three term Lanczos recurrence. WARNING: This variant does not track the loss of
   orthogonality that occurs in finite-precision arithmetic and, therefore, the
   generated vectors are not guaranteed to be (semi-)orthogonal.
*/
static PetscErrorCode EPSLocalLanczos(EPS eps,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *M,PetscBool *breakdown)
{
  PetscErrorCode ierr;
  PetscInt       i,j,m = *M;
  Mat            Op;
  PetscBool      *which,lwhich[100];
  PetscScalar    *hwork,lhwork[100];

  PetscFunctionBegin;
  if (m > 100) {
    ierr = PetscMalloc2(m,&which,m,&hwork);CHKERRQ(ierr);
  } else {
    which = lwhich;
    hwork = lhwork;
  }
  for (i=0;i<k;i++) which[i] = PETSC_TRUE;

  ierr = BVSetActiveColumns(eps->V,0,m);CHKERRQ(ierr);
  ierr = STGetOperator(eps->st,&Op);CHKERRQ(ierr);
  for (j=k;j<m;j++) {
    ierr = BVMatMultColumn(eps->V,Op,j);CHKERRQ(ierr);
    which[j] = PETSC_TRUE;
    if (j-2>=k) which[j-2] = PETSC_FALSE;
    ierr = BVOrthogonalizeSomeColumn(eps->V,j+1,which,hwork,beta+j,breakdown);CHKERRQ(ierr);
    alpha[j] = PetscRealPart(hwork[j]);
    if (*breakdown) {
      *M = j+1;
      break;
    } else {
      ierr = BVScaleColumn(eps->V,j+1,1/beta[j]);CHKERRQ(ierr);
    }
  }
  ierr = STRestoreOperator(eps->st,&Op);CHKERRQ(ierr);
  if (m > 100) {
    ierr = PetscFree2(which,hwork);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   DenseTridiagonal - Solves a real tridiagonal Hermitian Eigenvalue Problem.

   Input Parameters:
+  n   - dimension of the eigenproblem
.  D   - pointer to the array containing the diagonal elements
-  E   - pointer to the array containing the off-diagonal elements

   Output Parameters:
+  w  - pointer to the array to store the computed eigenvalues
-  V  - pointer to the array to store the eigenvectors

   Notes:
   If V is NULL then the eigenvectors are not computed.

   This routine use LAPACK routines xSTEVR.
*/
static PetscErrorCode DenseTridiagonal(PetscInt n_,PetscReal *D,PetscReal *E,PetscReal *w,PetscScalar *V)
{
  PetscErrorCode ierr;
  PetscReal      abstol = 0.0,vl,vu,*work;
  PetscBLASInt   il,iu,m,*isuppz,n,lwork,*iwork,liwork,info;
  const char     *jobz;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       i,j;
  PetscReal      *VV=NULL;
#endif

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(20*n_,&lwork);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(10*n_,&liwork);CHKERRQ(ierr);
  if (V) {
    jobz = "V";
#if defined(PETSC_USE_COMPLEX)
    ierr = PetscMalloc1(n*n,&VV);CHKERRQ(ierr);
#endif
  } else jobz = "N";
  ierr = PetscMalloc3(2*n,&isuppz,lwork,&work,liwork,&iwork);CHKERRQ(ierr);
  ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKstevr",LAPACKstevr_(jobz,"A",&n,D,E,&vl,&vu,&il,&iu,&abstol,&m,w,VV,&n,isuppz,work,&lwork,iwork,&liwork,&info));
#else
  PetscStackCallBLAS("LAPACKstevr",LAPACKstevr_(jobz,"A",&n,D,E,&vl,&vu,&il,&iu,&abstol,&m,w,V,&n,isuppz,work,&lwork,iwork,&liwork,&info));
#endif
  ierr = PetscFPTrapPop();CHKERRQ(ierr);
  SlepcCheckLapackInfo("stevr",info);
#if defined(PETSC_USE_COMPLEX)
  if (V) {
    for (i=0;i<n;i++)
      for (j=0;j<n;j++)
        V[i*n+j] = VV[i*n+j];
    ierr = PetscFree(VV);CHKERRQ(ierr);
  }
#endif
  ierr = PetscFree3(isuppz,work,iwork);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   EPSSelectiveLanczos - Selective reorthogonalization.
*/
static PetscErrorCode EPSSelectiveLanczos(EPS eps,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *M,PetscBool *breakdown,PetscReal anorm)
{
  PetscErrorCode ierr;
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;
  PetscInt       i,j,m = *M,n,nritz=0,nritzo;
  Vec            vj1,av;
  Mat            Op;
  PetscReal      *d,*e,*ritz,norm;
  PetscScalar    *Y,*hwork;
  PetscBool      *which;

  PetscFunctionBegin;
  ierr = PetscCalloc6(m+1,&d,m,&e,m,&ritz,m*m,&Y,m,&which,m,&hwork);CHKERRQ(ierr);
  for (i=0;i<k;i++) which[i] = PETSC_TRUE;
  ierr = STGetOperator(eps->st,&Op);CHKERRQ(ierr);

  for (j=k;j<m;j++) {
    ierr = BVSetActiveColumns(eps->V,0,m);CHKERRQ(ierr);

    /* Lanczos step */
    ierr = BVMatMultColumn(eps->V,Op,j);CHKERRQ(ierr);
    which[j] = PETSC_TRUE;
    if (j-2>=k) which[j-2] = PETSC_FALSE;
    ierr = BVOrthogonalizeSomeColumn(eps->V,j+1,which,hwork,&norm,breakdown);CHKERRQ(ierr);
    alpha[j] = PetscRealPart(hwork[j]);
    beta[j] = norm;
    if (*breakdown) {
      *M = j+1;
      break;
    }

    /* Compute eigenvalues and eigenvectors Y of the tridiagonal block */
    n = j-k+1;
    for (i=0;i<n;i++) {
      d[i] = alpha[i+k];
      e[i] = beta[i+k];
    }
    ierr = DenseTridiagonal(n,d,e,ritz,Y);CHKERRQ(ierr);

    /* Estimate ||A|| */
    for (i=0;i<n;i++)
      if (PetscAbsReal(ritz[i]) > anorm) anorm = PetscAbsReal(ritz[i]);

    /* Compute nearly converged Ritz vectors */
    nritzo = 0;
    for (i=0;i<n;i++) {
      if (norm*PetscAbsScalar(Y[i*n+n-1]) < PETSC_SQRT_MACHINE_EPSILON*anorm) nritzo++;
    }
    if (nritzo>nritz) {
      nritz = 0;
      for (i=0;i<n;i++) {
        if (norm*PetscAbsScalar(Y[i*n+n-1]) < PETSC_SQRT_MACHINE_EPSILON*anorm) {
          ierr = BVSetActiveColumns(eps->V,k,k+n);CHKERRQ(ierr);
          ierr = BVGetColumn(lanczos->AV,nritz,&av);CHKERRQ(ierr);
          ierr = BVMultVec(eps->V,1.0,0.0,av,Y+i*n);CHKERRQ(ierr);
          ierr = BVRestoreColumn(lanczos->AV,nritz,&av);CHKERRQ(ierr);
          nritz++;
        }
      }
    }
    if (nritz > 0) {
      ierr = BVGetColumn(eps->V,j+1,&vj1);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(lanczos->AV,0,nritz);CHKERRQ(ierr);
      ierr = BVOrthogonalizeVec(lanczos->AV,vj1,hwork,&norm,breakdown);CHKERRQ(ierr);
      ierr = BVRestoreColumn(eps->V,j+1,&vj1);CHKERRQ(ierr);
      if (*breakdown) {
        *M = j+1;
        break;
      }
    }
    ierr = BVScaleColumn(eps->V,j+1,1.0/norm);CHKERRQ(ierr);
  }

  ierr = STRestoreOperator(eps->st,&Op);CHKERRQ(ierr);
  ierr = PetscFree6(d,e,ritz,Y,which,hwork);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static void update_omega(PetscReal *omega,PetscReal *omega_old,PetscInt j,PetscReal *alpha,PetscReal *beta,PetscReal eps1,PetscReal anorm)
{
  PetscInt  k;
  PetscReal T,binv;

  PetscFunctionBegin;
  /* Estimate of contribution to roundoff errors from A*v
       fl(A*v) = A*v + f,
     where ||f|| \approx eps1*||A||.
     For a full matrix A, a rule-of-thumb estimate is eps1 = sqrt(n)*eps */
  T = eps1*anorm;
  binv = 1.0/beta[j+1];

  /* Update omega(1) using omega(0)==0 */
  omega_old[0]= beta[1]*omega[1] + (alpha[0]-alpha[j])*omega[0] - beta[j]*omega_old[0];
  if (omega_old[0] > 0) omega_old[0] = binv*(omega_old[0] + T);
  else omega_old[0] = binv*(omega_old[0] - T);

  /* Update remaining components */
  for (k=1;k<j-1;k++) {
    omega_old[k] = beta[k+1]*omega[k+1] + (alpha[k]-alpha[j])*omega[k] + beta[k]*omega[k-1] - beta[j]*omega_old[k];
    if (omega_old[k] > 0) omega_old[k] = binv*(omega_old[k] + T);
    else omega_old[k] = binv*(omega_old[k] - T);
  }
  omega_old[j-1] = binv*T;

  /* Swap omega and omega_old */
  for (k=0;k<j;k++) {
    omega[k] = omega_old[k];
    omega_old[k] = omega[k];
  }
  omega[j] = eps1;
  PetscFunctionReturnVoid();
}

static void compute_int(PetscBool *which,PetscReal *mu,PetscInt j,PetscReal delta,PetscReal eta)
{
  PetscInt  i,k,maxpos;
  PetscReal max;
  PetscBool found;

  PetscFunctionBegin;
  /* initialize which */
  found = PETSC_FALSE;
  maxpos = 0;
  max = 0.0;
  for (i=0;i<j;i++) {
    if (PetscAbsReal(mu[i]) >= delta) {
      which[i] = PETSC_TRUE;
      found = PETSC_TRUE;
    } else which[i] = PETSC_FALSE;
    if (PetscAbsReal(mu[i]) > max) {
      maxpos = i;
      max = PetscAbsReal(mu[i]);
    }
  }
  if (!found) which[maxpos] = PETSC_TRUE;

  for (i=0;i<j;i++) {
    if (which[i]) {
      /* find left interval */
      for (k=i;k>=0;k--) {
        if (PetscAbsReal(mu[k])<eta || which[k]) break;
        else which[k] = PETSC_TRUE;
      }
      /* find right interval */
      for (k=i+1;k<j;k++) {
        if (PetscAbsReal(mu[k])<eta || which[k]) break;
        else which[k] = PETSC_TRUE;
      }
    }
  }
  PetscFunctionReturnVoid();
}

/*
   EPSPartialLanczos - Partial reorthogonalization.
*/
static PetscErrorCode EPSPartialLanczos(EPS eps,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *M,PetscBool *breakdown,PetscReal anorm)
{
  PetscErrorCode ierr;
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;
  PetscInt       i,j,m = *M;
  Mat            Op;
  PetscReal      norm,*omega,lomega[100],*omega_old,lomega_old[100],eps1,delta,eta;
  PetscBool      *which,lwhich[100],*which2,lwhich2[100];
  PetscBool      reorth = PETSC_FALSE,force_reorth = PETSC_FALSE;
  PetscBool      fro = PETSC_FALSE,estimate_anorm = PETSC_FALSE;
  PetscScalar    *hwork,lhwork[100];

  PetscFunctionBegin;
  if (m>100) {
    ierr = PetscMalloc5(m,&omega,m,&omega_old,m,&which,m,&which2,m,&hwork);CHKERRQ(ierr);
  } else {
    omega     = lomega;
    omega_old = lomega_old;
    which     = lwhich;
    which2    = lwhich2;
    hwork     = lhwork;
  }

  eps1 = PetscSqrtReal((PetscReal)eps->n)*PETSC_MACHINE_EPSILON/2;
  delta = PETSC_SQRT_MACHINE_EPSILON/PetscSqrtReal((PetscReal)eps->ncv);
  eta = PetscPowReal(PETSC_MACHINE_EPSILON,3.0/4.0)/PetscSqrtReal((PetscReal)eps->ncv);
  if (anorm < 0.0) {
    anorm = 1.0;
    estimate_anorm = PETSC_TRUE;
  }
  for (i=0;i<PetscMax(100,m);i++) omega[i] = omega_old[i] = 0.0;
  for (i=0;i<k;i++) which[i] = PETSC_TRUE;

  ierr = BVSetActiveColumns(eps->V,0,m);CHKERRQ(ierr);
  ierr = STGetOperator(eps->st,&Op);CHKERRQ(ierr);
  for (j=k;j<m;j++) {
    ierr = BVMatMultColumn(eps->V,Op,j);CHKERRQ(ierr);
    if (fro) {
      /* Lanczos step with full reorthogonalization */
      ierr = BVOrthogonalizeColumn(eps->V,j+1,hwork,&norm,breakdown);CHKERRQ(ierr);
      alpha[j] = PetscRealPart(hwork[j]);
    } else {
      /* Lanczos step */
      which[j] = PETSC_TRUE;
      if (j-2>=k) which[j-2] = PETSC_FALSE;
      ierr = BVOrthogonalizeSomeColumn(eps->V,j+1,which,hwork,&norm,breakdown);CHKERRQ(ierr);
      alpha[j] = PetscRealPart(hwork[j]);
      beta[j] = norm;

      /* Estimate ||A|| if needed */
      if (estimate_anorm) {
        if (j>k) anorm = PetscMax(anorm,PetscAbsReal(alpha[j])+norm+beta[j-1]);
        else anorm = PetscMax(anorm,PetscAbsReal(alpha[j])+norm);
      }

      /* Check if reorthogonalization is needed */
      reorth = PETSC_FALSE;
      if (j>k) {
        update_omega(omega,omega_old,j,alpha,beta-1,eps1,anorm);
        for (i=0;i<j-k;i++) {
          if (PetscAbsReal(omega[i]) > delta) reorth = PETSC_TRUE;
        }
      }
      if (reorth || force_reorth) {
        for (i=0;i<k;i++) which2[i] = PETSC_FALSE;
        for (i=k;i<=j;i++) which2[i] = PETSC_TRUE;
        if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_PERIODIC) {
          /* Periodic reorthogonalization */
          if (force_reorth) force_reorth = PETSC_FALSE;
          else force_reorth = PETSC_TRUE;
          for (i=0;i<j-k;i++) omega[i] = eps1;
        } else {
          /* Partial reorthogonalization */
          if (force_reorth) force_reorth = PETSC_FALSE;
          else {
            force_reorth = PETSC_TRUE;
            compute_int(which2+k,omega,j-k,delta,eta);
            for (i=0;i<j-k;i++) {
              if (which2[i+k]) omega[i] = eps1;
            }
          }
        }
        ierr = BVOrthogonalizeSomeColumn(eps->V,j+1,which2,hwork,&norm,breakdown);CHKERRQ(ierr);
      }
    }

    if (*breakdown || norm < eps->n*anorm*PETSC_MACHINE_EPSILON) {
      *M = j+1;
      break;
    }
    if (!fro && norm*delta < anorm*eps1) {
      fro = PETSC_TRUE;
      ierr = PetscInfo1(eps,"Switching to full reorthogonalization at iteration %D\n",eps->its);CHKERRQ(ierr);
    }
    beta[j] = norm;
    ierr = BVScaleColumn(eps->V,j+1,1.0/norm);CHKERRQ(ierr);
  }

  ierr = STRestoreOperator(eps->st,&Op);CHKERRQ(ierr);
  if (m>100) {
    ierr = PetscFree5(omega,omega_old,which,which2,hwork);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   EPSBasicLanczos - Computes an m-step Lanczos factorization. The first k
   columns are assumed to be locked and therefore they are not modified. On
   exit, the following relation is satisfied:

                    OP * V - V * T = f * e_m^T

   where the columns of V are the Lanczos vectors, T is a tridiagonal matrix,
   f is the residual vector and e_m is the m-th vector of the canonical basis.
   The Lanczos vectors (together with vector f) are B-orthogonal (to working
   accuracy) if full reorthogonalization is being used, otherwise they are
   (B-)semi-orthogonal. On exit, beta contains the B-norm of f and the next
   Lanczos vector can be computed as v_{m+1} = f / beta.

   This function simply calls another function which depends on the selected
   reorthogonalization strategy.
*/
static PetscErrorCode EPSBasicLanczos(EPS eps,PetscReal *alpha,PetscReal *beta,PetscInt k,PetscInt *m,PetscBool *breakdown,PetscReal anorm)
{
  PetscErrorCode     ierr;
  EPS_LANCZOS        *lanczos = (EPS_LANCZOS*)eps->data;
  PetscScalar        *T;
  PetscInt           i,n=*m;
  PetscReal          betam;
  BVOrthogRefineType orthog_ref;
  Mat                Op;

  PetscFunctionBegin;
  switch (lanczos->reorthog) {
    case EPS_LANCZOS_REORTHOG_LOCAL:
      ierr = EPSLocalLanczos(eps,alpha,beta,k,m,breakdown);CHKERRQ(ierr);
      break;
    case EPS_LANCZOS_REORTHOG_FULL:
      ierr = STGetOperator(eps->st,&Op);CHKERRQ(ierr);
      ierr = BVMatLanczos(eps->V,Op,alpha,beta,k,m,breakdown);CHKERRQ(ierr);
      ierr = STRestoreOperator(eps->st,&Op);CHKERRQ(ierr);
      break;
    case EPS_LANCZOS_REORTHOG_SELECTIVE:
      ierr = EPSSelectiveLanczos(eps,alpha,beta,k,m,breakdown,anorm);CHKERRQ(ierr);
      break;
    case EPS_LANCZOS_REORTHOG_PERIODIC:
    case EPS_LANCZOS_REORTHOG_PARTIAL:
      ierr = EPSPartialLanczos(eps,alpha,beta,k,m,breakdown,anorm);CHKERRQ(ierr);
      break;
    case EPS_LANCZOS_REORTHOG_DELAYED:
      ierr = PetscMalloc1(n*n,&T);CHKERRQ(ierr);
      ierr = BVGetOrthogonalization(eps->V,NULL,&orthog_ref,NULL,NULL);CHKERRQ(ierr);
      if (orthog_ref == BV_ORTHOG_REFINE_NEVER) {
        ierr = EPSDelayedArnoldi1(eps,T,n,k,m,&betam,breakdown);CHKERRQ(ierr);
      } else {
        ierr = EPSDelayedArnoldi(eps,T,n,k,m,&betam,breakdown);CHKERRQ(ierr);
      }
      for (i=k;i<n-1;i++) {
        alpha[i] = PetscRealPart(T[n*i+i]);
        beta[i] = PetscRealPart(T[n*i+i+1]);
      }
      alpha[n-1] = PetscRealPart(T[n*(n-1)+n-1]);
      beta[n-1] = betam;
      ierr = PetscFree(T);CHKERRQ(ierr);
      break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_Lanczos(EPS eps)
{
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;
  PetscErrorCode ierr;
  PetscInt       nconv,i,j,k,l,x,n,*perm,restart,ncv=eps->ncv,r,ld;
  Vec            vi,vj,w;
  Mat            U;
  PetscScalar    *Y,*ritz,stmp;
  PetscReal      *d,*e,*bnd,anorm,beta,norm,rtmp,resnorm;
  PetscBool      breakdown;
  char           *conv,ctmp;

  PetscFunctionBegin;
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  ierr = PetscMalloc4(ncv,&ritz,ncv,&bnd,ncv,&perm,ncv,&conv);CHKERRQ(ierr);

  /* The first Lanczos vector is the normalized initial vector */
  ierr = EPSGetStartVector(eps,0,NULL);CHKERRQ(ierr);

  anorm = -1.0;
  nconv = 0;

  /* Restart loop */
  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;

    /* Compute an ncv-step Lanczos factorization */
    n = PetscMin(nconv+eps->mpd,ncv);
    ierr = DSGetArrayReal(eps->ds,DS_MAT_T,&d);CHKERRQ(ierr);
    e = d + ld;
    ierr = EPSBasicLanczos(eps,d,e,nconv,&n,&breakdown,anorm);CHKERRQ(ierr);
    beta = e[n-1];
    ierr = DSRestoreArrayReal(eps->ds,DS_MAT_T,&d);CHKERRQ(ierr);
    ierr = DSSetDimensions(eps->ds,n,0,nconv,0);CHKERRQ(ierr);
    ierr = DSSetState(eps->ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(eps->V,nconv,n);CHKERRQ(ierr);

    /* Solve projected problem */
    ierr = DSSolve(eps->ds,ritz,NULL);CHKERRQ(ierr);
    ierr = DSSort(eps->ds,ritz,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DSSynchronize(eps->ds,ritz,NULL);CHKERRQ(ierr);

    /* Estimate ||A|| */
    for (i=nconv;i<n;i++)
      anorm = PetscMax(anorm,PetscAbsReal(PetscRealPart(ritz[i])));

    /* Compute residual norm estimates as beta*abs(Y(m,:)) + eps*||A|| */
    ierr = DSGetArray(eps->ds,DS_MAT_Q,&Y);CHKERRQ(ierr);
    for (i=nconv;i<n;i++) {
      resnorm = beta*PetscAbsScalar(Y[n-1+i*ld]) + PETSC_MACHINE_EPSILON*anorm;
      ierr = (*eps->converged)(eps,ritz[i],eps->eigi[i],resnorm,&bnd[i],eps->convergedctx);CHKERRQ(ierr);
      if (bnd[i]<eps->tol) conv[i] = 'C';
      else conv[i] = 'N';
    }
    ierr = DSRestoreArray(eps->ds,DS_MAT_Q,&Y);CHKERRQ(ierr);

    /* purge repeated ritz values */
    if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_LOCAL) {
      for (i=nconv+1;i<n;i++) {
        if (conv[i] == 'C' && PetscAbsScalar((ritz[i]-ritz[i-1])/ritz[i]) < eps->tol) conv[i] = 'R';
      }
    }

    /* Compute restart vector */
    if (breakdown) {
      ierr = PetscInfo2(eps,"Breakdown in Lanczos method (it=%D norm=%g)\n",eps->its,(double)beta);CHKERRQ(ierr);
    } else {
      restart = nconv;
      while (restart<n && conv[restart] != 'N') restart++;
      if (restart >= n) {
        breakdown = PETSC_TRUE;
      } else {
        for (i=restart+1;i<n;i++) {
          if (conv[i] == 'N') {
            ierr = SlepcSCCompare(eps->sc,ritz[restart],0.0,ritz[i],0.0,&r);CHKERRQ(ierr);
            if (r>0) restart = i;
          }
        }
        ierr = DSGetArray(eps->ds,DS_MAT_Q,&Y);CHKERRQ(ierr);
        ierr = BVMultColumn(eps->V,1.0,0.0,n,Y+restart*ld+nconv);CHKERRQ(ierr);
        ierr = DSRestoreArray(eps->ds,DS_MAT_Q,&Y);CHKERRQ(ierr);
      }
    }

    /* Count and put converged eigenvalues first */
    for (i=nconv;i<n;i++) perm[i] = i;
    for (k=nconv;k<n;k++) {
      if (conv[perm[k]] != 'C') {
        j = k + 1;
        while (j<n && conv[perm[j]] != 'C') j++;
        if (j>=n) break;
        l = perm[k]; perm[k] = perm[j]; perm[j] = l;
      }
    }

    /* Sort eigenvectors according to permutation */
    ierr = DSGetArray(eps->ds,DS_MAT_Q,&Y);CHKERRQ(ierr);
    for (i=nconv;i<k;i++) {
      x = perm[i];
      if (x != i) {
        j = i + 1;
        while (perm[j] != i) j++;
        /* swap eigenvalues i and j */
        stmp = ritz[x]; ritz[x] = ritz[i]; ritz[i] = stmp;
        rtmp = bnd[x]; bnd[x] = bnd[i]; bnd[i] = rtmp;
        ctmp = conv[x]; conv[x] = conv[i]; conv[i] = ctmp;
        perm[j] = x; perm[i] = i;
        /* swap eigenvectors i and j */
        for (l=0;l<n;l++) {
          stmp = Y[l+x*ld]; Y[l+x*ld] = Y[l+i*ld]; Y[l+i*ld] = stmp;
        }
      }
    }
    ierr = DSRestoreArray(eps->ds,DS_MAT_Q,&Y);CHKERRQ(ierr);

    /* compute converged eigenvectors */
    ierr = DSGetMat(eps->ds,DS_MAT_Q,&U);CHKERRQ(ierr);
    ierr = BVMultInPlace(eps->V,U,nconv,k);CHKERRQ(ierr);
    ierr = MatDestroy(&U);CHKERRQ(ierr);

    /* purge spurious ritz values */
    if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_LOCAL) {
      for (i=nconv;i<k;i++) {
        ierr = BVGetColumn(eps->V,i,&vi);CHKERRQ(ierr);
        ierr = VecNorm(vi,NORM_2,&norm);CHKERRQ(ierr);
        ierr = VecScale(vi,1.0/norm);CHKERRQ(ierr);
        w = eps->work[0];
        ierr = STApply(eps->st,vi,w);CHKERRQ(ierr);
        ierr = VecAXPY(w,-ritz[i],vi);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->V,i,&vi);CHKERRQ(ierr);
        ierr = VecNorm(w,NORM_2,&norm);CHKERRQ(ierr);
        ierr = (*eps->converged)(eps,ritz[i],eps->eigi[i],norm,&bnd[i],eps->convergedctx);CHKERRQ(ierr);
        if (bnd[i]>=eps->tol) conv[i] = 'S';
      }
      for (i=nconv;i<k;i++) {
        if (conv[i] != 'C') {
          j = i + 1;
          while (j<k && conv[j] != 'C') j++;
          if (j>=k) break;
          /* swap eigenvalues i and j */
          stmp = ritz[j]; ritz[j] = ritz[i]; ritz[i] = stmp;
          rtmp = bnd[j]; bnd[j] = bnd[i]; bnd[i] = rtmp;
          ctmp = conv[j]; conv[j] = conv[i]; conv[i] = ctmp;
          /* swap eigenvectors i and j */
          ierr = BVGetColumn(eps->V,i,&vi);CHKERRQ(ierr);
          ierr = BVGetColumn(eps->V,j,&vj);CHKERRQ(ierr);
          ierr = VecSwap(vi,vj);CHKERRQ(ierr);
          ierr = BVRestoreColumn(eps->V,i,&vi);CHKERRQ(ierr);
          ierr = BVRestoreColumn(eps->V,j,&vj);CHKERRQ(ierr);
        }
      }
      k = i;
    }

    /* store ritz values and estimated errors */
    for (i=nconv;i<n;i++) {
      eps->eigr[i] = ritz[i];
      eps->errest[i] = bnd[i];
    }
    nconv = k;
    ierr = EPSMonitor(eps,eps->its,nconv,eps->eigr,eps->eigi,eps->errest,n);CHKERRQ(ierr);
    ierr = (*eps->stopping)(eps,eps->its,eps->max_it,nconv,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);

    if (eps->reason == EPS_CONVERGED_ITERATING) { /* copy restart vector */
      ierr = BVCopyColumn(eps->V,n,nconv);CHKERRQ(ierr);
      if (lanczos->reorthog == EPS_LANCZOS_REORTHOG_LOCAL && !breakdown) {
        /* Reorthonormalize restart vector */
        ierr = BVOrthonormalizeColumn(eps->V,nconv,PETSC_FALSE,NULL,&breakdown);CHKERRQ(ierr);
      }
      if (breakdown) {
        /* Use random vector for restarting */
        ierr = PetscInfo(eps,"Using random vector for restart\n");CHKERRQ(ierr);
        ierr = EPSGetStartVector(eps,nconv,&breakdown);CHKERRQ(ierr);
      }
      if (breakdown) { /* give up */
        eps->reason = EPS_DIVERGED_BREAKDOWN;
        ierr = PetscInfo(eps,"Unable to generate more start vectors\n");CHKERRQ(ierr);
      }
    }
  }
  eps->nconv = nconv;

  ierr = PetscFree4(ritz,bnd,perm,conv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_Lanczos(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode         ierr;
  EPS_LANCZOS            *lanczos = (EPS_LANCZOS*)eps->data;
  PetscBool              flg;
  EPSLanczosReorthogType reorthog;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS Lanczos Options");CHKERRQ(ierr);

    ierr = PetscOptionsEnum("-eps_lanczos_reorthog","Lanczos reorthogonalization","EPSLanczosSetReorthog",EPSLanczosReorthogTypes,(PetscEnum)lanczos->reorthog,(PetscEnum*)&reorthog,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSLanczosSetReorthog(eps,reorthog);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLanczosSetReorthog_Lanczos(EPS eps,EPSLanczosReorthogType reorthog)
{
  EPS_LANCZOS *lanczos = (EPS_LANCZOS*)eps->data;

  PetscFunctionBegin;
  switch (reorthog) {
    case EPS_LANCZOS_REORTHOG_LOCAL:
    case EPS_LANCZOS_REORTHOG_FULL:
    case EPS_LANCZOS_REORTHOG_DELAYED:
    case EPS_LANCZOS_REORTHOG_SELECTIVE:
    case EPS_LANCZOS_REORTHOG_PERIODIC:
    case EPS_LANCZOS_REORTHOG_PARTIAL:
      if (lanczos->reorthog != reorthog) {
        lanczos->reorthog = reorthog;
        eps->state = EPS_STATE_INITIAL;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid reorthogonalization type");
  }
  PetscFunctionReturn(0);
}

/*@
   EPSLanczosSetReorthog - Sets the type of reorthogonalization used during the Lanczos
   iteration.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  reorthog - the type of reorthogonalization

   Options Database Key:
.  -eps_lanczos_reorthog - Sets the reorthogonalization type (either 'local', 'selective',
                         'periodic', 'partial', 'full' or 'delayed')

   Level: advanced

.seealso: EPSLanczosGetReorthog(), EPSLanczosReorthogType
@*/
PetscErrorCode EPSLanczosSetReorthog(EPS eps,EPSLanczosReorthogType reorthog)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,reorthog,2);
  ierr = PetscTryMethod(eps,"EPSLanczosSetReorthog_C",(EPS,EPSLanczosReorthogType),(eps,reorthog));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSLanczosGetReorthog_Lanczos(EPS eps,EPSLanczosReorthogType *reorthog)
{
  EPS_LANCZOS *lanczos = (EPS_LANCZOS*)eps->data;

  PetscFunctionBegin;
  *reorthog = lanczos->reorthog;
  PetscFunctionReturn(0);
}

/*@
   EPSLanczosGetReorthog - Gets the type of reorthogonalization used during
   the Lanczos iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  reorthog - the type of reorthogonalization

   Level: advanced

.seealso: EPSLanczosSetReorthog(), EPSLanczosReorthogType
@*/
PetscErrorCode EPSLanczosGetReorthog(EPS eps,EPSLanczosReorthogType *reorthog)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(reorthog,2);
  ierr = PetscUseMethod(eps,"EPSLanczosGetReorthog_C",(EPS,EPSLanczosReorthogType*),(eps,reorthog));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_Lanczos(EPS eps)
{
  PetscErrorCode ierr;
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;

  PetscFunctionBegin;
  ierr = BVDestroy(&lanczos->AV);CHKERRQ(ierr);
  lanczos->allocsize = 0;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_Lanczos(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLanczosSetReorthog_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLanczosGetReorthog_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_Lanczos(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  EPS_LANCZOS    *lanczos = (EPS_LANCZOS*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  %s reorthogonalization\n",EPSLanczosReorthogTypes[lanczos->reorthog]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_Lanczos(EPS eps)
{
  EPS_LANCZOS    *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&ctx);CHKERRQ(ierr);
  eps->data = (void*)ctx;

  eps->useds = PETSC_TRUE;

  eps->ops->solve          = EPSSolve_Lanczos;
  eps->ops->setup          = EPSSetUp_Lanczos;
  eps->ops->setfromoptions = EPSSetFromOptions_Lanczos;
  eps->ops->destroy        = EPSDestroy_Lanczos;
  eps->ops->reset          = EPSReset_Lanczos;
  eps->ops->view           = EPSView_Lanczos;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_Hermitian;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLanczosSetReorthog_C",EPSLanczosSetReorthog_Lanczos);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSLanczosGetReorthog_C",EPSLanczosGetReorthog_Lanczos);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

