#!/bin/sh
DAWES_SRC=../dawes/src/*
mpifort -c ../path_resolution.f90 -g -ffree-line-length-0
mpifort -c ../dawes/src/splin.f90 -g -ffree-line-length-0
mpifort -c ../dawes/src/dynamic_parameters_3D.f90 -g -ffree-line-length-0
mpifort -c ../../../src/base/coordinate_conversion.f90 -g -ffree-line-length-0
mpifort ../ozone_pes_alt.f90 path_resolution.o coordinate_conversion.o ${DAWES_SRC} -ffree-line-length-0 -g -w -o ozone_pes_alt
