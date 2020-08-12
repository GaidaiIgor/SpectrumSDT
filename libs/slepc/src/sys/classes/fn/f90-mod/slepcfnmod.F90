!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcfndef
        use slepcsys
#include <../src/sys/classes/fn/f90-mod/slepcfn.h>
        end module

        module slepcfn
        use slepcfndef
#include <../src/sys/classes/fn/f90-mod/slepcfn.h90>
        interface
#include <../src/sys/classes/fn/f90-mod/ftn-auto-interfaces/slepcfn.h90>
        end interface
        end module

