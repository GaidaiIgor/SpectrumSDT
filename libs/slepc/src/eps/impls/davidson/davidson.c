/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Skeleton of Davidson solver. Actual solvers are GD and JD.
*/

#include "davidson.h"

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-davidson,\n"
  "   author = \"E. Romero and J. E. Roman\",\n"
  "   title = \"A parallel implementation of {Davidson} methods for large-scale eigenvalue problems in {SLEPc}\",\n"
  "   journal = \"{ACM} Trans. Math. Software\",\n"
  "   volume = \"40\",\n"
  "   number = \"2\",\n"
  "   pages = \"13:1--13:29\",\n"
  "   year = \"2014,\"\n"
  "   doi = \"https://doi.org/10.1145/2543696\"\n"
  "}\n";

PetscErrorCode EPSSetUp_XD(EPS eps)
{
  PetscErrorCode ierr;
  EPS_DAVIDSON   *data = (EPS_DAVIDSON*)eps->data;
  dvdDashboard   *dvd = &data->ddb;
  dvdBlackboard  b;
  PetscInt       min_size_V,bs,initv,nmat;
  Mat            A,B;
  KSP            ksp;
  PetscBool      ipB,ispositive;
  HarmType_t     harm;
  InitType_t     init;
  PetscScalar    target;

  PetscFunctionBegin;
  /* Setup EPS options and get the problem specification */
  bs = data->blocksize;
  if (bs <= 0) bs = 1;
  if (eps->ncv!=PETSC_DEFAULT) {
    if (eps->ncv<eps->nev) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The value of ncv must be at least nev");
  } else if (eps->mpd!=PETSC_DEFAULT) eps->ncv = eps->mpd + eps->nev + bs;
  else if (eps->nev<500) eps->ncv = PetscMin(eps->n-bs,PetscMax(2*eps->nev,eps->nev+15))+bs;
  else eps->ncv = PetscMin(eps->n-bs,eps->nev+500)+bs;
  if (eps->mpd==PETSC_DEFAULT) eps->mpd = eps->ncv;
  if (eps->mpd > eps->ncv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The mpd has to be less or equal than ncv");
  if (eps->mpd < 2) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The mpd has to be greater than 2");
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(100*eps->ncv,2*eps->n);
  if (!eps->which) eps->which = EPS_LARGEST_MAGNITUDE;
  if (eps->ishermitian && (eps->which==EPS_LARGEST_IMAGINARY || eps->which==EPS_SMALLEST_IMAGINARY)) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Wrong value of eps->which");
  if (!(eps->nev + bs <= eps->ncv)) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The ncv has to be greater than nev plus blocksize");
  if (eps->trueres) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"-eps_true_residual is temporally disable in this solver.");

  ierr = EPSXDSetRestart_XD(eps,data->minv,data->plusk);CHKERRQ(ierr);
  min_size_V = data->minv;
  if (!min_size_V) min_size_V = PetscMin(PetscMax(bs,5),eps->mpd/2);
  if (!(min_size_V+bs <= eps->mpd)) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The value of minv must be less than mpd minus blocksize");
  initv = data->initialsize;
  if (eps->mpd < initv) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"The initv has to be less or equal than mpd");

  /* Change the default sigma to inf if necessary */
  if (eps->which == EPS_LARGEST_MAGNITUDE || eps->which == EPS_LARGEST_REAL || eps->which == EPS_LARGEST_IMAGINARY) {
    ierr = STSetDefaultShift(eps->st,PETSC_MAX_REAL);CHKERRQ(ierr);
  }

  /* Set up preconditioner */
  ierr = STSetUp(eps->st);CHKERRQ(ierr);

  /* Setup problem specification in dvd */
  ierr = STGetNumMatrices(eps->st,&nmat);CHKERRQ(ierr);
  ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
  if (nmat>1) { ierr = STGetMatrix(eps->st,1,&B);CHKERRQ(ierr); }
  ierr = EPSReset_XD(eps);CHKERRQ(ierr);
  ierr = PetscMemzero(dvd,sizeof(dvdDashboard));CHKERRQ(ierr);
  dvd->A = A; dvd->B = eps->isgeneralized? B: NULL;
  ispositive = eps->ispositive;
  dvd->sA = DVD_MAT_IMPLICIT | (eps->ishermitian? DVD_MAT_HERMITIAN: 0) | ((ispositive && !eps->isgeneralized) ? DVD_MAT_POS_DEF: 0);
  /* Asume -eps_hermitian means hermitian-definite in generalized problems */
  if (!ispositive && !eps->isgeneralized && eps->ishermitian) ispositive = PETSC_TRUE;
  if (!eps->isgeneralized) dvd->sB = DVD_MAT_IMPLICIT | DVD_MAT_HERMITIAN | DVD_MAT_IDENTITY | DVD_MAT_UNITARY | DVD_MAT_POS_DEF;
  else dvd->sB = DVD_MAT_IMPLICIT | (eps->ishermitian? DVD_MAT_HERMITIAN: 0) | (ispositive? DVD_MAT_POS_DEF: 0);
  ipB = (dvd->B && data->ipB && DVD_IS(dvd->sB,DVD_MAT_HERMITIAN))?PETSC_TRUE:PETSC_FALSE;
  dvd->sEP = ((!eps->isgeneralized || (eps->isgeneralized && ipB))? DVD_EP_STD: 0) | (ispositive? DVD_EP_HERMITIAN: 0) | ((eps->problem_type == EPS_GHIEP && ipB) ? DVD_EP_INDEFINITE : 0);
  if (data->ipB && !ipB) data->ipB = PETSC_FALSE;
  dvd->correctXnorm = (dvd->B && (DVD_IS(dvd->sB,DVD_MAT_HERMITIAN)||DVD_IS(dvd->sEP,DVD_EP_INDEFINITE)))?PETSC_TRUE:PETSC_FALSE;
  dvd->nev        = eps->nev;
  dvd->which      = eps->which;
  dvd->withTarget = PETSC_TRUE;
  switch (eps->which) {
    case EPS_TARGET_MAGNITUDE:
    case EPS_TARGET_IMAGINARY:
      dvd->target[0] = target = eps->target;
      dvd->target[1] = 1.0;
      break;
    case EPS_TARGET_REAL:
      dvd->target[0] = PetscRealPart(target = eps->target);
      dvd->target[1] = 1.0;
      break;
    case EPS_LARGEST_REAL:
    case EPS_LARGEST_MAGNITUDE:
    case EPS_LARGEST_IMAGINARY: /* TODO: think about this case */
      dvd->target[0] = 1.0;
      dvd->target[1] = target = 0.0;
      break;
    case EPS_SMALLEST_MAGNITUDE:
    case EPS_SMALLEST_REAL:
    case EPS_SMALLEST_IMAGINARY: /* TODO: think about this case */
      dvd->target[0] = target = 0.0;
      dvd->target[1] = 1.0;
      break;
    case EPS_WHICH_USER:
      ierr = STGetShift(eps->st,&target);CHKERRQ(ierr);
      dvd->target[0] = target;
      dvd->target[1] = 1.0;
      break;
    case EPS_ALL:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported option: which == EPS_ALL");
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported value of option 'which'");
  }
  dvd->tol = (eps->tol==PETSC_DEFAULT)? SLEPC_DEFAULT_TOL: eps->tol;
  dvd->eps = eps;

  /* Setup the extraction technique */
  if (!eps->extraction) {
    if (ipB || ispositive) eps->extraction = EPS_RITZ;
    else {
      switch (eps->which) {
        case EPS_TARGET_REAL:
        case EPS_TARGET_MAGNITUDE:
        case EPS_TARGET_IMAGINARY:
        case EPS_SMALLEST_MAGNITUDE:
        case EPS_SMALLEST_REAL:
        case EPS_SMALLEST_IMAGINARY:
          eps->extraction = EPS_HARMONIC;
          break;
        case EPS_LARGEST_REAL:
        case EPS_LARGEST_MAGNITUDE:
        case EPS_LARGEST_IMAGINARY:
          eps->extraction = EPS_HARMONIC_LARGEST;
          break;
        default:
          eps->extraction = EPS_RITZ;
      }
    }
  }
  switch (eps->extraction) {
    case EPS_RITZ:              harm = DVD_HARM_NONE; break;
    case EPS_HARMONIC:          harm = DVD_HARM_RR; break;
    case EPS_HARMONIC_RELATIVE: harm = DVD_HARM_RRR; break;
    case EPS_HARMONIC_RIGHT:    harm = DVD_HARM_REIGS; break;
    case EPS_HARMONIC_LARGEST:  harm = DVD_HARM_LEIGS; break;
    default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Unsupported extraction type");
  }

  /* Setup the type of starting subspace */
  init = data->krylovstart? DVD_INITV_KRYLOV: DVD_INITV_CLASSIC;

  /* Preconfigure dvd */
  ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr);
  ierr = dvd_schm_basic_preconf(dvd,&b,eps->mpd,min_size_V,bs,initv,PetscAbs(eps->nini),data->plusk,harm,ksp,init,eps->trackall,data->ipB,data->doubleexp);CHKERRQ(ierr);

  /* Allocate memory */
  ierr = EPSAllocateSolution(eps,0);CHKERRQ(ierr);

  /* Setup orthogonalization */
  ierr = EPS_SetInnerProduct(eps);CHKERRQ(ierr);
  if (!(ipB && dvd->B)) {
    ierr = BVSetMatrix(eps->V,NULL,PETSC_FALSE);CHKERRQ(ierr);
  }

  /* Configure dvd for a basic GD */
  ierr = dvd_schm_basic_conf(dvd,&b,eps->mpd,min_size_V,bs,initv,PetscAbs(eps->nini),data->plusk,harm,dvd->withTarget,target,ksp,data->fix,init,eps->trackall,data->ipB,data->dynamic,data->doubleexp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_XD(EPS eps)
{
  EPS_DAVIDSON   *data = (EPS_DAVIDSON*)eps->data;
  dvdDashboard   *d = &data->ddb;
  PetscInt       l,k;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscCitationsRegister(citation,&cited);CHKERRQ(ierr);
  /* Call the starting routines */
  ierr = EPSDavidsonFLCall(d->startList,d);CHKERRQ(ierr);

  while (eps->reason == EPS_CONVERGED_ITERATING) {

    /* Initialize V, if it is needed */
    ierr = BVGetActiveColumns(d->eps->V,&l,&k);CHKERRQ(ierr);
    if (l == k) { ierr = d->initV(d);CHKERRQ(ierr); }

    /* Find the best approximated eigenpairs in V, X */
    ierr = d->calcPairs(d);CHKERRQ(ierr);

    /* Test for convergence */
    ierr = (*eps->stopping)(eps,eps->its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);
    if (eps->reason != EPS_CONVERGED_ITERATING) break;

    /* Expand the subspace */
    ierr = d->updateV(d);CHKERRQ(ierr);

    /* Monitor progress */
    eps->nconv = d->nconv;
    eps->its++;
    ierr = BVGetActiveColumns(d->eps->V,NULL,&k);CHKERRQ(ierr);
    ierr = EPSMonitor(eps,eps->its,eps->nconv+d->npreconv,eps->eigr,eps->eigi,eps->errest,PetscMin(k,eps->nev));CHKERRQ(ierr);
  }

  /* Call the ending routines */
  ierr = EPSDavidsonFLCall(d->endList,d);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_XD(EPS eps)
{
  EPS_DAVIDSON   *data = (EPS_DAVIDSON*)eps->data;
  dvdDashboard   *dvd = &data->ddb;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Call step destructors and destroys the list */
  ierr = EPSDavidsonFLCall(dvd->destroyList,dvd);CHKERRQ(ierr);
  ierr = EPSDavidsonFLDestroy(&dvd->destroyList);CHKERRQ(ierr);
  ierr = EPSDavidsonFLDestroy(&dvd->startList);CHKERRQ(ierr);
  ierr = EPSDavidsonFLDestroy(&dvd->endList);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDSetKrylovStart_XD(EPS eps,PetscBool krylovstart)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  data->krylovstart = krylovstart;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDGetKrylovStart_XD(EPS eps,PetscBool *krylovstart)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *krylovstart = data->krylovstart;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDSetBlockSize_XD(EPS eps,PetscInt blocksize)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (blocksize == PETSC_DEFAULT || blocksize == PETSC_DECIDE) blocksize = 1;
  if (blocksize <= 0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid blocksize value");
  data->blocksize = blocksize;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDGetBlockSize_XD(EPS eps,PetscInt *blocksize)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *blocksize = data->blocksize;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDSetRestart_XD(EPS eps,PetscInt minv,PetscInt plusk)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (minv == PETSC_DEFAULT || minv == PETSC_DECIDE) minv = 5;
  if (minv <= 0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid minv value");
  if (plusk == PETSC_DEFAULT || plusk == PETSC_DECIDE) plusk = eps->problem_type == EPS_GHIEP?0:1;
  if (plusk < 0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid plusk value");
  data->minv = minv;
  data->plusk = plusk;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDGetRestart_XD(EPS eps,PetscInt *minv,PetscInt *plusk)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (minv) *minv = data->minv;
  if (plusk) *plusk = data->plusk;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDGetInitialSize_XD(EPS eps,PetscInt *initialsize)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *initialsize = data->initialsize;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDSetInitialSize_XD(EPS eps,PetscInt initialsize)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (initialsize == PETSC_DEFAULT || initialsize == PETSC_DECIDE) initialsize = 5;
  if (initialsize <= 0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid initial size value");
  data->initialsize = initialsize;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDSetBOrth_XD(EPS eps,PetscBool borth)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  data->ipB = borth;
  PetscFunctionReturn(0);
}

PetscErrorCode EPSXDGetBOrth_XD(EPS eps,PetscBool *borth)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *borth = data->ipB;
  PetscFunctionReturn(0);
}

/*
  EPSComputeVectors_XD - Compute eigenvectors from the vectors
  provided by the eigensolver. This version is intended for solvers
  that provide Schur vectors from the QZ decomposition. Given the partial
  Schur decomposition OP*V=V*T, the following steps are performed:
      1) compute eigenvectors of (S,T): S*Z=T*Z*D
      2) compute eigenvectors of OP: X=V*Z
 */
PetscErrorCode EPSComputeVectors_XD(EPS eps)
{
  PetscErrorCode ierr;
  Mat            X;
  PetscBool      symm;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)eps->ds,DSHEP,&symm);CHKERRQ(ierr);
  if (symm) PetscFunctionReturn(0);
  ierr = DSVectors(eps->ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);

  /* V <- V * X */
  ierr = DSGetMat(eps->ds,DS_MAT_X,&X);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,0,eps->nconv);CHKERRQ(ierr);
  ierr = BVMultInPlace(eps->V,X,0,eps->nconv);CHKERRQ(ierr);
  ierr = MatDestroy(&X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
