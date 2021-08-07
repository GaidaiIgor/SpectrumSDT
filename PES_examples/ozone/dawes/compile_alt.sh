#!/bin/sh
DAWES_SRC=../dawes/src/*
mpifort -c ../../../path_resolution.f90 -ffree-line-length-0 -g
mpifort -c ../dawes/src/splin.f90 -ffree-line-length-0 -g
mpifort -c ../dawes/src/dynamic_parameters_3D.f90 -ffree-line-length-0 -g
mpifort -c ../../../../src/base/coordinate_conversion.f90 -ffree-line-length-0 -g
mpifort -c ../../../pes_utils.f90 -ffree-line-length-0 -g
mpifort -c ../../ozone_utils.f90 -ffree-line-length-0 -g
mpifort path_resolution.o coordinate_conversion.o pes_utils.o ozone_utils.o ../ozone_pes_alt.f90 ${DAWES_SRC} -ffree-line-length-0 -g -w -o ozone_pes_alt

