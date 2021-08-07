#!/bin/sh
ftn -c ../../pes_utils.f90 -g
ftn pes_utils.o ../water_pes.f90 ../partridge/h2opes.f ../shirin/potb.f -g -w -o water_pes

