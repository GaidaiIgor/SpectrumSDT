!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcsysdef
        use petscsysdef
        use petscmatdef
        use petscsys
#include <../src/sys/f90-mod/slepcsys.h>
        end module

        module slepcsys
        use slepcsysdef
#include <../src/sys/f90-mod/slepcsys.h90>
        interface
#include <../src/sys/f90-mod/ftn-auto-interfaces/slepcsys.h90>
        end interface
        end module

