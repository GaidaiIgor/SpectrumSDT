/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepcsys.h>
#include <petsc/private/fortranimpl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define slepcinitializefortran_     SLEPCINITIALIZEFORTRAN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define slepcinitializefortran_     slepcinitializefortran
#endif

/*@C
   SlepcInitializeFortran - Routine that should be called from C after
   the call to SlepcInitialize() if one is using a C main program
   that calls Fortran routines that in turn call SLEPc routines.

   Collective on PETSC_COMM_WORLD

   Level: beginner

   Notes:
   SlepcInitializeFortran() initializes some of the default SLEPc variables
   for use in Fortran if a user's main program is written in C.
   SlepcInitializeFortran() is NOT needed if a user's main
   program is written in Fortran; in this case, just calling
   SlepcInitialize() in the main (Fortran) program is sufficient.

.seealso:  SlepcInitialize()

@*/

PetscErrorCode SlepcInitializeFortran(void)
{
  PetscInitializeFortran();
  return 0;
}

SLEPC_EXTERN void slepcinitializefortran_(PetscErrorCode *info)
{
  *info = SlepcInitializeFortran();
}

