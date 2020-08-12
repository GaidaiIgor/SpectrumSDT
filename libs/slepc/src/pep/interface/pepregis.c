/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/pepimpl.h>      /*I "slepcpep.h" I*/

SLEPC_EXTERN PetscErrorCode PEPCreate_Linear(PEP);
SLEPC_EXTERN PetscErrorCode PEPCreate_QArnoldi(PEP);
SLEPC_EXTERN PetscErrorCode PEPCreate_TOAR(PEP);
SLEPC_EXTERN PetscErrorCode PEPCreate_STOAR(PEP);
SLEPC_EXTERN PetscErrorCode PEPCreate_JD(PEP);

/*@C
   PEPRegisterAll - Registers all the solvers in the PEP package.

   Not Collective

   Level: advanced

.seealso:  PEPRegister()
@*/
PetscErrorCode PEPRegisterAll(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (PEPRegisterAllCalled) PetscFunctionReturn(0);
  PEPRegisterAllCalled = PETSC_TRUE;
  ierr = PEPRegister(PEPLINEAR,PEPCreate_Linear);CHKERRQ(ierr);
  ierr = PEPRegister(PEPQARNOLDI,PEPCreate_QArnoldi);CHKERRQ(ierr);
  ierr = PEPRegister(PEPTOAR,PEPCreate_TOAR);CHKERRQ(ierr);
  ierr = PEPRegister(PEPSTOAR,PEPCreate_STOAR);CHKERRQ(ierr);
  ierr = PEPRegister(PEPJD,PEPCreate_JD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

