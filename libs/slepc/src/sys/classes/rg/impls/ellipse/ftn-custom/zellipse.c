/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepcrg.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define rgellipsegetparameters000_ RGELLIPSEGETPARAMETERS000
#define rgellipsegetparameters100_ RGELLIPSEGETPARAMETERS100
#define rgellipsegetparameters010_ RGELLIPSEGETPARAMETERS010
#define rgellipsegetparameters110_ RGELLIPSEGETPARAMETERS110
#define rgellipsegetparameters001_ RGELLIPSEGETPARAMETERS001
#define rgellipsegetparameters101_ RGELLIPSEGETPARAMETERS101
#define rgellipsegetparameters011_ RGELLIPSEGETPARAMETERS011
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define rgellipsegetparameters000_ rgellipsegetparameters000
#define rgellipsegetparameters100_ rgellipsegetparameters100
#define rgellipsegetparameters010_ rgellipsegetparameters010
#define rgellipsegetparameters110_ rgellipsegetparameters110
#define rgellipsegetparameters001_ rgellipsegetparameters001
#define rgellipsegetparameters101_ rgellipsegetparameters101
#define rgellipsegetparameters011_ rgellipsegetparameters011
#endif

SLEPC_EXTERN void rgellipsegetparameters_(RG *rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscErrorCode *ierr)
{
  CHKFORTRANNULLSCALAR(center);
  CHKFORTRANNULLREAL(radius);
  CHKFORTRANNULLREAL(vscale);
  *ierr = RGEllipseGetParameters(*rg,center,radius,vscale);
}

SLEPC_EXTERN void rgellipsegetparameters000_(RG *rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscErrorCode *ierr)
{
  rgellipsegetparameters_(rg,center,radius,vscale,ierr);
}

SLEPC_EXTERN void rgellipsegetparameters100_(RG *rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscErrorCode *ierr)
{
  rgellipsegetparameters_(rg,center,radius,vscale,ierr);
}

SLEPC_EXTERN void rgellipsegetparameters010_(RG *rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscErrorCode *ierr)
{
  rgellipsegetparameters_(rg,center,radius,vscale,ierr);
}

SLEPC_EXTERN void rgellipsegetparameters110_(RG *rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscErrorCode *ierr)
{
  rgellipsegetparameters_(rg,center,radius,vscale,ierr);
}

SLEPC_EXTERN void rgellipsegetparameters001_(RG *rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscErrorCode *ierr)
{
  rgellipsegetparameters_(rg,center,radius,vscale,ierr);
}

SLEPC_EXTERN void rgellipsegetparameters101_(RG *rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscErrorCode *ierr)
{
  rgellipsegetparameters_(rg,center,radius,vscale,ierr);
}

SLEPC_EXTERN void rgellipsegetparameters011_(RG *rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscErrorCode *ierr)
{
  rgellipsegetparameters_(rg,center,radius,vscale,ierr);
}

