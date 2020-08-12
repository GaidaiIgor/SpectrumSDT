#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

ALL: lib

#requiresdefine  'SLEPC_HAVE_BLOPEX'

CFLAGS   = ${BLOPEX_INCLUDE} -DBlopexInt=PetscInt
FFLAGS   =
SOURCEC  = blopex.c slepc-interface.c petsc-interface.c
SOURCEF  =
SOURCEH  = blopex.h petsc-interface.h
LIBBASE  = libslepceps
DIRS     =
MANSEC   = EPS
LOCDIR   = src/eps/impls/external/blopex/

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common


