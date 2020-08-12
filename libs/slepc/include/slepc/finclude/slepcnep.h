!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the NEP object in SLEPc
!
#if !defined(SLEPCNEPDEF_H)
#define SLEPCNEPDEF_H

#include "slepc/finclude/slepcbv.h"
#include "slepc/finclude/slepcds.h"
#include "slepc/finclude/slepcrg.h"
#include "slepc/finclude/slepcfn.h"
#include "slepc/finclude/slepceps.h"
#include "slepc/finclude/slepcpep.h"

#define NEP type(tNEP)

#define NEPType            character*(80)
#define NEPProblemType     PetscEnum
#define NEPConvergedReason PetscEnum
#define NEPErrorType       PetscEnum
#define NEPWhich           PetscEnum
#define NEPConv            PetscEnum
#define NEPStop            PetscEnum
#define NEPRefine          PetscEnum
#define NEPRefineScheme    PetscEnum

#define NEPRII       'rii'
#define NEPSLP       'slp'
#define NEPNARNOLDI  'narnoldi'
#define NEPCISS      'ciss'
#define NEPINTERPOL  'interpol'
#define NEPNLEIGS    'nleigs'

#endif

