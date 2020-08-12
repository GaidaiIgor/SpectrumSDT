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
#if !defined(SLEPCBVDEF_H)
#define SLEPCBVDEF_H

#include "petsc/finclude/petscmat.h"

#define BV type(tBV)

#define BVType             character*(80)
#define BVOrthogType       PetscEnum
#define BVOrthogRefineType PetscEnum
#define BVOrthogBlockType  PetscEnum
#define BVMatMultType      PetscEnum

#define BVMAT        'mat'
#define BVSVEC       'svec'
#define BVVECS       'vecs'
#define BVCONTIGUOUS 'contiguous'
#define BVTENSOR     'tensor'

#endif
