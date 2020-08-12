!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        module slepcnepdef
        use slepcbvdef
        use slepcrgdef
        use slepcdsdef
        use slepcfndef
        use slepcepsdef
        use slepcpepdef
#include <../src/nep/f90-mod/slepcnep.h>
        end module

        module slepcnep
        use slepcnepdef
#include <../src/nep/f90-mod/slepcnep.h90>
        interface
#include <../src/nep/f90-mod/ftn-auto-interfaces/slepcnep.h90>
        end interface
        end module

