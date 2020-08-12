/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/rgimpl.h>      /*I "slepcrg.h" I*/

SLEPC_EXTERN PetscErrorCode RGCreate_Interval(RG);
SLEPC_EXTERN PetscErrorCode RGCreate_Ellipse(RG);
SLEPC_EXTERN PetscErrorCode RGCreate_Ring(RG);
SLEPC_EXTERN PetscErrorCode RGCreate_Polygon(RG);

/*@C
   RGRegisterAll - Registers all of the regions in the RG package.

   Not Collective

   Level: advanced
@*/
PetscErrorCode RGRegisterAll(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (RGRegisterAllCalled) PetscFunctionReturn(0);
  RGRegisterAllCalled = PETSC_TRUE;
  ierr = RGRegister(RGINTERVAL,RGCreate_Interval);CHKERRQ(ierr);
  ierr = RGRegister(RGELLIPSE,RGCreate_Ellipse);CHKERRQ(ierr);
  ierr = RGRegister(RGRING,RGCreate_Ring);CHKERRQ(ierr);
  ierr = RGRegister(RGPOLYGON,RGCreate_Polygon);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

