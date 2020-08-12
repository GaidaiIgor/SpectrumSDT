!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcrgdef
        use slepcsys
#include <../src/sys/classes/rg/f90-mod/slepcrg.h>
        end module

        module slepcrg
        use slepcrgdef
#include <../src/sys/classes/rg/f90-mod/slepcrg.h90>
        interface
#include <../src/sys/classes/rg/f90-mod/ftn-auto-interfaces/slepcrg.h90>
        end interface
        end module

