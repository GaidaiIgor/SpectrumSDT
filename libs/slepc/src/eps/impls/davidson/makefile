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

CFLAGS   =
FFLAGS   =
SOURCEC  = davidson.c dvdcalcpairs.c dvdimprovex.c dvdinitv.c \
           dvdschm.c dvdtestconv.c dvdupdatev.c dvdutils.c dvdgd2.c
SOURCEF  =
SOURCEH  = davidson.h
LIBBASE  = libslepceps
DIRS     = gd jd
LOCDIR   = src/eps/impls/davidson/
MANSEC   = EPS

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

