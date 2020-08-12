!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the SVD object in SLEPc
!
#include "slepc/finclude/slepcsvd.h"

      type tSVD
        PetscFortranAddr:: v PETSC_FORTRAN_TYPE_INITIALIZE
      end type tSVD

      SVD, parameter :: SLEPC_NULL_SVD = tSVD(0)

      PetscEnum, parameter :: SVD_CONVERGED_TOL          =  1
      PetscEnum, parameter :: SVD_CONVERGED_USER         =  2
      PetscEnum, parameter :: SVD_DIVERGED_ITS           = -1
      PetscEnum, parameter :: SVD_DIVERGED_BREAKDOWN     = -2
      PetscEnum, parameter :: SVD_CONVERGED_ITERATING    =  0

      PetscEnum, parameter :: SVD_LARGEST                =  0
      PetscEnum, parameter :: SVD_SMALLEST               =  1

      PetscEnum, parameter :: SVD_ERROR_ABSOLUTE         =  0
      PetscEnum, parameter :: SVD_ERROR_RELATIVE         =  1

      PetscEnum, parameter :: SVD_CONV_ABS               =  0
      PetscEnum, parameter :: SVD_CONV_REL               =  1
      PetscEnum, parameter :: SVD_CONV_USER              =  2

      PetscEnum, parameter :: SVD_STOP_BASIC             =  0
      PetscEnum, parameter :: SVD_STOP_USER              =  1

      PetscEnum, parameter :: SVD_PRIMME_HYBRID          =  1
      PetscEnum, parameter :: SVD_PRIMME_NORMALEQUATIONS =  2
      PetscEnum, parameter :: SVD_PRIMME_AUGMENTED       =  3

!
!   Possible arguments to SVDMonitorSet()
!
      external SVDMONITORALL
      external SVDMONITORLG
      external SVDMONITORLGALL
      external SVDMONITORCONVERGED
      external SVDMONITORFIRST

!
!  End of Fortran include file for the SVD package in SLEPc
!
