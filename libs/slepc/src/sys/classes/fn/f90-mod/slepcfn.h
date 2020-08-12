!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the FN object in SLEPc
!
#include "slepc/finclude/slepcfn.h"

      type tFN
        PetscFortranAddr:: v PETSC_FORTRAN_TYPE_INITIALIZE
      end type tFN

      FN, parameter :: SLEPC_NULL_FN = tFN(0)

      PetscEnum, parameter :: FN_COMBINE_ADD           =  0
      PetscEnum, parameter :: FN_COMBINE_MULTIPLY      =  1
      PetscEnum, parameter :: FN_COMBINE_DIVIDE        =  2
      PetscEnum, parameter :: FN_COMBINE_COMPOSE       =  3

!
!  End of Fortran include file for the FN package in SLEPc
!
