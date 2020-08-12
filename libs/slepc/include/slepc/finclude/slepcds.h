!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Include file for Fortran use of the DS object in SLEPc
!
#if !defined(SLEPCDSDEF_H)
#define SLEPCDSDEF_H

#include "petsc/finclude/petscmat.h"
#include "slepc/finclude/slepcfn.h"

#define DS type(tDS)

#define DSType         character*(80)
#define DSStateType    PetscEnum
#define DSMatType      PetscEnum
#define DSParallelType PetscEnum

#define DSHEP       'hep'
#define DSNHEP      'nhep'
#define DSGHEP      'ghep'
#define DSGHIEP     'ghiep'
#define DSGNHEP     'gnhep'
#define DSSVD       'svd'
#define DSPEP       'pep'
#define DSNEP       'nep'

#endif

