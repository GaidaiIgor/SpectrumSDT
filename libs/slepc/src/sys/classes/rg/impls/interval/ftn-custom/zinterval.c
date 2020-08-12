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
#define rgintervalgetendpoints0000_ RGINTERVALGETENDPOINTS0000
#define rgintervalgetendpoints1000_ RGINTERVALGETENDPOINTS1000
#define rgintervalgetendpoints0100_ RGINTERVALGETENDPOINTS0100
#define rgintervalgetendpoints1100_ RGINTERVALGETENDPOINTS1100
#define rgintervalgetendpoints0010_ RGINTERVALGETENDPOINTS0010
#define rgintervalgetendpoints1010_ RGINTERVALGETENDPOINTS1010
#define rgintervalgetendpoints0110_ RGINTERVALGETENDPOINTS0110
#define rgintervalgetendpoints1110_ RGINTERVALGETENDPOINTS1110
#define rgintervalgetendpoints0001_ RGINTERVALGETENDPOINTS0001
#define rgintervalgetendpoints1001_ RGINTERVALGETENDPOINTS1001
#define rgintervalgetendpoints0101_ RGINTERVALGETENDPOINTS0101
#define rgintervalgetendpoints1101_ RGINTERVALGETENDPOINTS1101
#define rgintervalgetendpoints0011_ RGINTERVALGETENDPOINTS0011
#define rgintervalgetendpoints1011_ RGINTERVALGETENDPOINTS1011
#define rgintervalgetendpoints0111_ RGINTERVALGETENDPOINTS0111
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define rgintervalgetendpoints0000_ rgintervalgetendpoints0000
#define rgintervalgetendpoints1000_ rgintervalgetendpoints1000
#define rgintervalgetendpoints0100_ rgintervalgetendpoints0100
#define rgintervalgetendpoints1100_ rgintervalgetendpoints1100
#define rgintervalgetendpoints0010_ rgintervalgetendpoints0010
#define rgintervalgetendpoints1010_ rgintervalgetendpoints1010
#define rgintervalgetendpoints0110_ rgintervalgetendpoints0110
#define rgintervalgetendpoints1110_ rgintervalgetendpoints1110
#define rgintervalgetendpoints0001_ rgintervalgetendpoints0001
#define rgintervalgetendpoints1001_ rgintervalgetendpoints1001
#define rgintervalgetendpoints0101_ rgintervalgetendpoints0101
#define rgintervalgetendpoints1101_ rgintervalgetendpoints1101
#define rgintervalgetendpoints0011_ rgintervalgetendpoints0011
#define rgintervalgetendpoints1011_ rgintervalgetendpoints1011
#define rgintervalgetendpoints0111_ rgintervalgetendpoints0111
#endif

SLEPC_EXTERN void rgintervalgetendpoints_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  CHKFORTRANNULLREAL(a);
  CHKFORTRANNULLREAL(b);
  CHKFORTRANNULLREAL(c);
  CHKFORTRANNULLREAL(d);
  *ierr = RGIntervalGetEndpoints(*rg,a,b,c,d);
}

SLEPC_EXTERN void rgintervalgetendpoints0000_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints1000_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints0100_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints1100_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints0010_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints1010_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints0110_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints1110_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints0001_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints1001_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints0101_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints1101_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints0011_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints1011_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

SLEPC_EXTERN void rgintervalgetendpoints0111_(RG *rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d,PetscErrorCode *ierr)
{
  rgintervalgetendpoints_(rg,a,b,c,d,ierr);
}

