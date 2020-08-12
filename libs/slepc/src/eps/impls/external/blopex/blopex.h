/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if !defined(SLEPC_BLOPEX_H)
#define SLEPC_BLOPEX_H

#include <lobpcg.h>
#include "petsc-interface.h"

SLEPC_INTERN PetscInt slepc_blopex_useconstr;

extern int
SLEPCSetupInterpreter(mv_InterfaceInterpreter *ii);

#endif

