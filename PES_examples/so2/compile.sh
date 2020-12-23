#!/bin/sh
mpifort -c ../path_resolution.f90 -g -ffree-line-length-0
mpifort ../so2_pes.f90 ../klos/so2-X-mrcif12.f path_resolution.o -ffree-line-length-0 -g -w -o so2_pes
