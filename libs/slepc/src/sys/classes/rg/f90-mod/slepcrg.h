!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the RG object in SLEPc
!
#include "slepc/finclude/slepcrg.h"

      type tRG
        PetscFortranAddr:: v PETSC_FORTRAN_TYPE_INITIALIZE
      end type tRG

      RG, parameter :: SLEPC_NULL_RG = tRG(0)

!
!  End of Fortran include file for the RG package in SLEPc
!
