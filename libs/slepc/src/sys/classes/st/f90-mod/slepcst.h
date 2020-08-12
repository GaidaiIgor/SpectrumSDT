!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the ST object in SLEPc
!
#include "slepc/finclude/slepcst.h"

      type tST
        PetscFortranAddr:: v PETSC_FORTRAN_TYPE_INITIALIZE
      end type tST

      ST, parameter :: SLEPC_NULL_ST = tST(0)

      PetscEnum, parameter :: ST_MATMODE_COPY          =  0
      PetscEnum, parameter :: ST_MATMODE_INPLACE       =  1
      PetscEnum, parameter :: ST_MATMODE_SHELL         =  2

!
!  End of Fortran include file for the ST package in SLEPc
!
