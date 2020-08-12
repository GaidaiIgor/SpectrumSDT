/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This file implements a wrapper to the ARPACK package
*/

#include <slepc/private/epsimpl.h>
#include "arpack.h"

PetscErrorCode EPSSetUp_ARPACK(EPS eps)
{
  PetscErrorCode ierr;
  PetscInt       ncv;
  PetscBool      istrivial;
  EPS_ARPACK     *ar = (EPS_ARPACK*)eps->data;

  PetscFunctionBegin;
  if (eps->ncv!=PETSC_DEFAULT) {
    if (eps->ncv<eps->nev+2) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"The value of ncv must be at least nev+2");
  } else eps->ncv = PetscMin(PetscMax(20,2*eps->nev+1),eps->n); /* set default value of ncv */
  if (eps->mpd!=PETSC_DEFAULT) { ierr = PetscInfo(eps,"Warning: parameter mpd ignored\n");CHKERRQ(ierr); }
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(300,(PetscInt)(2*eps->n/eps->ncv));
  if (!eps->which) eps->which = EPS_LARGEST_MAGNITUDE;

  ncv = eps->ncv;
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscFree(ar->rwork);CHKERRQ(ierr);
  ierr = PetscMalloc1(ncv,&ar->rwork);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,ncv*sizeof(PetscReal));CHKERRQ(ierr);
  ar->lworkl = 3*ncv*ncv+5*ncv;
  ierr = PetscFree(ar->workev);CHKERRQ(ierr);
  ierr = PetscMalloc1(3*ncv,&ar->workev);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,3*ncv*sizeof(PetscScalar));CHKERRQ(ierr);
#else
  if (eps->ishermitian) {
    ar->lworkl = ncv*(ncv+8);
  } else {
    ar->lworkl = 3*ncv*ncv+6*ncv;
    ierr = PetscFree(ar->workev);CHKERRQ(ierr);
    ierr = PetscMalloc1(3*ncv,&ar->workev);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)eps,3*ncv*sizeof(PetscScalar));CHKERRQ(ierr);
  }
#endif
  ierr = PetscFree(ar->workl);CHKERRQ(ierr);
  ierr = PetscMalloc1(ar->lworkl,&ar->workl);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,ar->lworkl*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscFree(ar->select);CHKERRQ(ierr);
  ierr = PetscMalloc1(ncv,&ar->select);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,ncv*sizeof(PetscInt));CHKERRQ(ierr);
  ierr = PetscFree(ar->workd);CHKERRQ(ierr);
  ierr = PetscMalloc1(3*eps->nloc,&ar->workd);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)eps,3*eps->nloc*sizeof(PetscScalar));CHKERRQ(ierr);

  if (eps->extraction) { ierr = PetscInfo(eps,"Warning: extraction type ignored\n");CHKERRQ(ierr); }

  if (eps->balance!=EPS_BALANCE_NONE) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Balancing not supported in the Arpack interface");
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");
  if (eps->stopping!=EPSStoppingBasic) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"External packages do not support user-defined stopping test");

  ierr = EPSAllocateSolution(eps,0);CHKERRQ(ierr);
  ierr = EPS_SetInnerProduct(eps);CHKERRQ(ierr);
  ierr = EPSSetWorkVecs(eps,2);CHKERRQ(ierr);

  ierr = RGIsTrivial(eps->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support region filtering");
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_ARPACK(EPS eps)
{
  PetscErrorCode ierr;
  EPS_ARPACK     *ar = (EPS_ARPACK*)eps->data;
  char           bmat[1],howmny[] = "A";
  const char     *which;
  PetscInt       n,iparam[11],ipntr[14],ido,info,nev,ncv,rvec;
#if !defined(PETSC_HAVE_MPIUNI)
  MPI_Fint       fcomm;
#endif
  PetscScalar    sigmar,*pV,*resid;
  Vec            x,y,w = eps->work[0];
  Mat            A;
  PetscBool      isSinv,isShift;
#if !defined(PETSC_USE_COMPLEX)
  PetscScalar    sigmai = 0.0;
#endif

  PetscFunctionBegin;
  nev = eps->nev;
  ncv = eps->ncv;
#if !defined(PETSC_HAVE_MPIUNI)
  fcomm = MPI_Comm_c2f(PetscObjectComm((PetscObject)eps));
#endif
  n = eps->nloc;
  ierr = EPSGetStartVector(eps,0,NULL);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(eps->V,0,0);CHKERRQ(ierr);  /* just for deflation space */
  ierr = BVCopyVec(eps->V,0,eps->work[1]);CHKERRQ(ierr);
  ierr = BVGetArray(eps->V,&pV);CHKERRQ(ierr);
  ierr = VecGetArray(eps->work[1],&resid);CHKERRQ(ierr);

  ido  = 0;            /* first call to reverse communication interface */
  info = 1;            /* indicates an initial vector is provided */
  iparam[0] = 1;       /* use exact shifts */
  iparam[2] = eps->max_it;  /* max Arnoldi iterations */
  iparam[3] = 1;       /* blocksize */
  iparam[4] = 0;       /* number of converged Ritz values */

  /*
     Computational modes ([]=not supported):
            symmetric    non-symmetric    complex
        1     1  'I'        1  'I'         1  'I'
        2     3  'I'        3  'I'         3  'I'
        3     2  'G'        2  'G'         2  'G'
        4     3  'G'        3  'G'         3  'G'
        5   [ 4  'G' ]    [ 3  'G' ]
        6   [ 5  'G' ]    [ 4  'G' ]
   */
  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&isSinv);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&isShift);CHKERRQ(ierr);
  ierr = STGetShift(eps->st,&sigmar);CHKERRQ(ierr);
  ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
  ierr = MatCreateVecsEmpty(A,&x,&y);CHKERRQ(ierr);

  if (isSinv) {
    /* shift-and-invert mode */
    iparam[6] = 3;
    if (eps->ispositive) bmat[0] = 'G';
    else bmat[0] = 'I';
  } else if (isShift && eps->ispositive) {
    /* generalized shift mode with B positive definite */
    iparam[6] = 2;
    bmat[0] = 'G';
  } else {
    /* regular mode */
    if (eps->ishermitian && eps->isgeneralized) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Spectral transformation not supported by ARPACK hermitian solver");
    iparam[6] = 1;
    bmat[0] = 'I';
  }

#if !defined(PETSC_USE_COMPLEX)
    if (eps->ishermitian) {
      switch (eps->which) {
        case EPS_TARGET_MAGNITUDE:
        case EPS_LARGEST_MAGNITUDE:  which = "LM"; break;
        case EPS_SMALLEST_MAGNITUDE: which = "SM"; break;
        case EPS_TARGET_REAL:
        case EPS_LARGEST_REAL:       which = "LA"; break;
        case EPS_SMALLEST_REAL:      which = "SA"; break;
        default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Wrong value of eps->which");
      }
    } else {
#endif
      switch (eps->which) {
        case EPS_TARGET_MAGNITUDE:
        case EPS_LARGEST_MAGNITUDE:  which = "LM"; break;
        case EPS_SMALLEST_MAGNITUDE: which = "SM"; break;
        case EPS_TARGET_REAL:
        case EPS_LARGEST_REAL:       which = "LR"; break;
        case EPS_SMALLEST_REAL:      which = "SR"; break;
        case EPS_TARGET_IMAGINARY:
        case EPS_LARGEST_IMAGINARY:  which = "LI"; break;
        case EPS_SMALLEST_IMAGINARY: which = "SI"; break;
        default: SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Wrong value of eps->which");
      }
#if !defined(PETSC_USE_COMPLEX)
    }
#endif

  do {

#if !defined(PETSC_USE_COMPLEX)
    if (eps->ishermitian) {
      PetscStackCall("ARPACKsaupd",ARPACKsaupd_(&fcomm,&ido,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&n,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,&info));
    } else {
      PetscStackCall("ARPACKnaupd",ARPACKnaupd_(&fcomm,&ido,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&n,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,&info));
    }
#else
    PetscStackCall("ARPACKnaupd",ARPACKnaupd_(&fcomm,&ido,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&n,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,ar->rwork,&info));
#endif

    if (ido == -1 || ido == 1 || ido == 2) {
      if (ido == 1 && iparam[6] == 3 && bmat[0] == 'G') {
        /* special case for shift-and-invert with B semi-positive definite*/
        ierr = VecPlaceArray(x,&ar->workd[ipntr[2]-1]);CHKERRQ(ierr);
      } else {
        ierr = VecPlaceArray(x,&ar->workd[ipntr[0]-1]);CHKERRQ(ierr);
      }
      ierr = VecPlaceArray(y,&ar->workd[ipntr[1]-1]);CHKERRQ(ierr);

      if (ido == -1) {
        /* Y = OP * X for for the initialization phase to
           force the starting vector into the range of OP */
        ierr = STApply(eps->st,x,y);CHKERRQ(ierr);
      } else if (ido == 2) {
        /* Y = B * X */
        ierr = BVApplyMatrix(eps->V,x,y);CHKERRQ(ierr);
      } else { /* ido == 1 */
        if (iparam[6] == 3 && bmat[0] == 'G') {
          /* Y = OP * X for shift-and-invert with B semi-positive definite */
          ierr = STMatSolve(eps->st,x,y);CHKERRQ(ierr);
        } else if (iparam[6] == 2) {
          /* X=A*X Y=B^-1*X for shift with B positive definite */
          ierr = MatMult(A,x,y);CHKERRQ(ierr);
          if (sigmar != 0.0) {
            ierr = BVApplyMatrix(eps->V,x,w);CHKERRQ(ierr);
            ierr = VecAXPY(y,sigmar,w);CHKERRQ(ierr);
          }
          ierr = VecCopy(y,x);CHKERRQ(ierr);
          ierr = STMatSolve(eps->st,x,y);CHKERRQ(ierr);
        } else {
          /* Y = OP * X */
          ierr = STApply(eps->st,x,y);CHKERRQ(ierr);
        }
        ierr = BVOrthogonalizeVec(eps->V,y,NULL,NULL,NULL);CHKERRQ(ierr);
      }

      ierr = VecResetArray(x);CHKERRQ(ierr);
      ierr = VecResetArray(y);CHKERRQ(ierr);
    } else if (ido != 99) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Internal error in ARPACK reverse comunication interface (ido=%d)",ido);

  } while (ido != 99);

  eps->nconv = iparam[4];
  eps->its = iparam[2];

  if (info==3) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"No shift could be applied in xxAUPD.\nTry increasing the size of NCV relative to NEV");
  else if (info!=0 && info!=1) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Error reported by ARPACK subroutine xxAUPD (%d)",(int)info);

  rvec = PETSC_TRUE;

  if (eps->nconv > 0) {
#if !defined(PETSC_USE_COMPLEX)
    if (eps->ishermitian) {
      PetscStackCall("ARPACKseupd",ARPACKseupd_(&fcomm,&rvec,howmny,ar->select,eps->eigr,pV,&n,&sigmar,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&n,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,&info));
    } else {
      PetscStackCall("ARPACKneupd",ARPACKneupd_(&fcomm,&rvec,howmny,ar->select,eps->eigr,eps->eigi,pV,&n,&sigmar,&sigmai,ar->workev,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&n,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,&info));
    }
#else
    PetscStackCall("ARPACKneupd",ARPACKneupd_(&fcomm,&rvec,howmny,ar->select,eps->eigr,pV,&n,&sigmar,ar->workev,bmat,&n,which,&nev,&eps->tol,resid,&ncv,pV,&n,iparam,ipntr,ar->workd,ar->workl,&ar->lworkl,ar->rwork,&info));
#endif
    if (info!=0) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_LIB,"Error reported by ARPACK subroutine xxEUPD (%d)",(int)info);
  }

  ierr = BVRestoreArray(eps->V,&pV);CHKERRQ(ierr);
  ierr = VecRestoreArray(eps->work[1],&resid);CHKERRQ(ierr);
  if (eps->nconv >= eps->nev) eps->reason = EPS_CONVERGED_TOL;
  else eps->reason = EPS_DIVERGED_ITS;

  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSBackTransform_ARPACK(EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      isSinv;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&isSinv);CHKERRQ(ierr);
  if (!isSinv) {
    ierr = EPSBackTransform_Default(eps);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_ARPACK(EPS eps)
{
  PetscErrorCode ierr;
  EPS_ARPACK     *ar = (EPS_ARPACK*)eps->data;

  PetscFunctionBegin;
  ierr = PetscFree(ar->workev);CHKERRQ(ierr);
  ierr = PetscFree(ar->workl);CHKERRQ(ierr);
  ierr = PetscFree(ar->select);CHKERRQ(ierr);
  ierr = PetscFree(ar->workd);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscFree(ar->rwork);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_ARPACK(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_ARPACK(EPS eps)
{
  EPS_ARPACK     *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&ctx);CHKERRQ(ierr);
  eps->data = (void*)ctx;

  eps->ops->solve          = EPSSolve_ARPACK;
  eps->ops->setup          = EPSSetUp_ARPACK;
  eps->ops->destroy        = EPSDestroy_ARPACK;
  eps->ops->reset          = EPSReset_ARPACK;
  eps->ops->backtransform  = EPSBackTransform_ARPACK;
  PetscFunctionReturn(0);
}

