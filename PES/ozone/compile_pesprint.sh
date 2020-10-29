#!/bin/sh

DAWES_SRC=../dawes/src/*
mpifort -c ../dawes/src/path_resolution.f90 -g -ffree-line-length-0
mpifort -c ../dawes/src/splin.f90 -g -ffree-line-length-0
mpifort -c ../dawes/src/dynamic_parameters_3D.f90 -g -ffree-line-length-0
mpifort ../pesprint.f90 ${DAWES_SRC} -ffree-line-length-0 -g -w -o pesprint
