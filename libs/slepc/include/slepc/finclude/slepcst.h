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
#if !defined(SLEPCSTDEF_H)
#define SLEPCSTDEF_H

#include "petsc/finclude/petscksp.h"
#include "slepc/finclude/slepcbv.h"

#define ST type(tST)

#define STType     character*(80)
#define STMatMode  PetscEnum

#define STSHELL    'shell'
#define STSHIFT    'shift'
#define STSINVERT  'sinvert'
#define STCAYLEY   'cayley'
#define STPRECOND  'precond'
#define STFILTER   'filter'

#endif

