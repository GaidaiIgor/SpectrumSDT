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
#if !defined(SLEPCRGDEF_H)
#define SLEPCRGDEF_H

#include "petsc/finclude/petscsys.h"

#define RG type(tRG)

#define RGType      character*(80)

#define RGINTERVAL  'interval'
#define RGPOLYGON   'polygon'
#define RGELLIPSE   'ellipse'
#define RGRING      'ring'

#endif

