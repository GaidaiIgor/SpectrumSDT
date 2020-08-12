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
#if !defined(SLEPCFNDEF_H)
#define SLEPCFNDEF_H

#include "petsc/finclude/petscmat.h"

#define FN type(tFN)

#define FNType         character*(80)
#define FNCombineType  PetscEnum
#define FNParallelType PetscEnum

#define FNCOMBINE  'combine'
#define FNRATIONAL 'rational'
#define FNEXP      'exp'
#define FNLOG      'log'
#define FNPHI      'phi'
#define FNSQRT     'sqrt'
#define FNINVSQRT  'invsqrt'

#endif

