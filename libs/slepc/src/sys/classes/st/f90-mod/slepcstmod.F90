!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcstdef
        use petscksp
        use slepcsys
        use slepcbvdef
#include <../src/sys/classes/st/f90-mod/slepcst.h>
        end module

        module slepcst
        use slepcstdef
#include <../src/sys/classes/st/f90-mod/slepcst.h90>
        interface
#include <../src/sys/classes/st/f90-mod/ftn-auto-interfaces/slepcst.h90>
        end interface
        end module

