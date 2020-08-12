!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the PEP object in SLEPc
!
#if !defined(SLEPCPEPDEF_H)
#define SLEPCPEPDEF_H

#include "slepc/finclude/slepcbv.h"
#include "slepc/finclude/slepcst.h"
#include "slepc/finclude/slepcds.h"
#include "slepc/finclude/slepcrg.h"
#include "slepc/finclude/slepceps.h"

#define PEP type(tPEP)

#define PEPType            character*(80)
#define PEPProblemType     PetscEnum
#define PEPWhich           PetscEnum
#define PEPBasis           PetscEnum
#define PEPScale           PetscEnum
#define PEPRefine          PetscEnum
#define PEPRefineScheme    PetscEnum
#define PEPExtract         PetscEnum
#define PEPConv            PetscEnum
#define PEPStop            PetscEnum
#define PEPErrorType       PetscEnum
#define PEPConvergedReason PetscEnum
#define PEPJDProjection    PetscEnum

#define PEPLINEAR    'linear'
#define PEPQARNOLDI  'qarnoldi'
#define PEPTOAR      'toar'
#define PEPSTOAR     'stoar'
#define PEPJD        'jd'

#endif

