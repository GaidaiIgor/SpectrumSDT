/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "davidson"

   Step: test for convergence
*/

#include "davidson.h"

static PetscBool dvd_testconv_slepc_0(dvdDashboard *d,PetscScalar eigvr,PetscScalar eigvi,PetscReal r,PetscReal *err)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*d->eps->converged)(d->eps,eigvr,eigvi,r,err,d->eps->convergedctx);CHKERRABORT(PetscObjectComm((PetscObject)d->eps),ierr);
  PetscFunctionReturn(PetscNot(*err>=d->eps->tol));
}

PetscErrorCode dvd_testconv_slepc(dvdDashboard *d, dvdBlackboard *b)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    ierr = PetscFree(d->testConv_data);CHKERRQ(ierr);
    d->testConv = dvd_testconv_slepc_0;
  }
  PetscFunctionReturn(0);
}

