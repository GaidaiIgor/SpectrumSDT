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
#if !defined(SLEPCSVDDEF_H)
#define SLEPCSVDDEF_H

#include "slepc/finclude/slepcbv.h"
#include "slepc/finclude/slepcds.h"
#include "slepc/finclude/slepceps.h"

#define SVD type(tSVD)

#define SVDType            character*(80)
#define SVDConvergedReason PetscEnum
#define SVDErrorType       PetscEnum
#define SVDWhich           PetscEnum
#define SVDConv            PetscEnum
#define SVDStop            PetscEnum
#define SVDPRIMMEMethod    PetscEnum

#define SVDCROSS     'cross'
#define SVDCYCLIC    'cyclic'
#define SVDLAPACK    'lapack'
#define SVDLANCZOS   'lanczos'
#define SVDTRLANCZOS 'trlanczos'
#define SVDPRIMME    'primme'

#endif

