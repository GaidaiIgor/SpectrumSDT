!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcsvddef
        use slepcbvdef
        use slepcepsdef
#include <../src/svd/f90-mod/slepcsvd.h>
        end module

        module slepcsvd
        use slepcsvddef
#include <../src/svd/f90-mod/slepcsvd.h90>
        interface
#include <../src/svd/f90-mod/ftn-auto-interfaces/slepcsvd.h90>
        end interface
        end module

