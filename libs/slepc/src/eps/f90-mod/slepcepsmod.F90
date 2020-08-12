!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcepsdef
        use slepcsys
        use slepcstdef
        use slepcbvdef
        use slepcrgdef
        use slepcdsdef
        use slepclmedef
        use petscsnesdef
#include <../src/eps/f90-mod/slepceps.h>
        end module

        module slepceps
        use slepcepsdef
#include <../src/eps/f90-mod/slepceps.h90>
        interface
#include <../src/eps/f90-mod/ftn-auto-interfaces/slepceps.h90>
        end interface
        end module

