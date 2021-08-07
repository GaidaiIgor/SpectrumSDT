#!/bin/sh
mpifort -c ../../pes_utils.f90 -ffree-line-length-0 -g
mpifort pes_utils.o ../water_pes.f90 ../partridge/h2opes.f ../shirin/potb.f -ffree-line-length-0 -g -w -o water_pes

