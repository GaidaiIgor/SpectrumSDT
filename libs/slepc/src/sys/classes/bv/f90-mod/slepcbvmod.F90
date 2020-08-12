!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcbvdef
        use petscmat
        use petscvec
        use slepcsys
#include <../src/sys/classes/bv/f90-mod/slepcbv.h>
        end module

        module slepcbv
        use slepcbvdef
#include <../src/sys/classes/bv/f90-mod/slepcbv.h90>
        interface
#include <../src/sys/classes/bv/f90-mod/ftn-auto-interfaces/slepcbv.h90>
        end interface
        end module

