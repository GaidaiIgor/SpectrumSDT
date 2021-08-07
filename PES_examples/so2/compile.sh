#!/bin/sh
mpifort -c ../../path_resolution.f90 -ffree-line-length-0 -g
mpifort -c ../../pes_utils.f90 -ffree-line-length-0 -g
mpifort path_resolution.o pes_utils.o ../so2_pes.f90 ../klos/so2-X-mrcif12.f -ffree-line-length-0 -g -w -o so2_pes

