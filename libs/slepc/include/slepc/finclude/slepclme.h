!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the LME object in SLEPc
!
#if !defined(SLEPCLMEDEF_H)
#define SLEPCLMEDEF_H

#include "petsc/finclude/petscmat.h"
#include "slepc/finclude/slepcbv.h"

#define LME type(tLME)

#define LMEType            character*(80)
#define LMEConvergedReason PetscEnum
#define LMEProblemType     PetscEnum

#define LMEKRYLOV      'krylov'

#endif
