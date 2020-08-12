!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepclmedef
        use slepcbvdef
        use petscvecdef
        use petscmatdef
#include <../src/lme/f90-mod/slepclme.h>
        end module

        module slepclme
        use slepclmedef
        use slepcbv
        use petscvec
        use petscmat
#include <../src/lme/f90-mod/slepclme.h90>
        interface
#include <../src/lme/f90-mod/ftn-auto-interfaces/slepclme.h90>
        end interface
        end module

