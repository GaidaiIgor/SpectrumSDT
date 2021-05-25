#!/bin/sh
DAWES_SRC=../dawes/src/*
ftn -c ../path_resolution.f90 -g
ftn -c ../dawes/src/splin.f90 -g
ftn -c ../dawes/src/dynamic_parameters_3D.f90 -g
ftn -c ../../../../src/base/coordinate_conversion.f90 -g
ftn -c ../../../pes_utils.f90 -g
ftn -c ../../ozone_utils.f90 -g
ftn ../ozone_pes.f90 path_resolution.o coordinate_conversion.o pes_utils.o ozone_utils.o ${DAWES_SRC} -g -w -o ozone_pes
