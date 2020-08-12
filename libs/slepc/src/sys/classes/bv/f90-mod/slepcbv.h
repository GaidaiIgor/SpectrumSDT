!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the BV object in SLEPc
!
#include "slepc/finclude/slepcbv.h"

      type tBV
        PetscFortranAddr:: v PETSC_FORTRAN_TYPE_INITIALIZE
      end type tBV

      BV, parameter :: SLEPC_NULL_BV = tBV(0)

      PetscEnum, parameter :: BV_ORTHOG_CGS             =  0
      PetscEnum, parameter :: BV_ORTHOG_MGS             =  1

      PetscEnum, parameter :: BV_ORTHOG_REFINE_IFNEEDED =  0
      PetscEnum, parameter :: BV_ORTHOG_REFINE_NEVER    =  1
      PetscEnum, parameter :: BV_ORTHOG_REFINE_ALWAYS   =  2

      PetscEnum, parameter :: BV_ORTHOG_BLOCK_GS        =  0
      PetscEnum, parameter :: BV_ORTHOG_BLOCK_CHOL      =  1
      PetscEnum, parameter :: BV_ORTHOG_BLOCK_TSQR      =  2
      PetscEnum, parameter :: BV_ORTHOG_BLOCK_TSQRCHOL  =  3
      PetscEnum, parameter :: BV_ORTHOG_BLOCK_SVQB      =  4

      PetscEnum, parameter :: BV_MATMULT_VECS           =  0
      PetscEnum, parameter :: BV_MATMULT_MAT            =  1
      PetscEnum, parameter :: BV_MATMULT_MAT_SAVE       =  2

!
!  End of Fortran include file for the BV package in SLEPc
!
