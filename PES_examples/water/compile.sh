#!/bin/sh
mpifort ../water_pes.f90 ../partridge/h2opes.f ../shirin/potb.f  -ffree-line-length-0 -g -w -o water_pes
