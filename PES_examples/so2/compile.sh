#!/bin/sh
ftn -c ../../path_resolution.f90 -g
ftn -c ../../pes_utils.f90 -g
ftn path_resolution.o pes_utils.o ../so2_pes.f90 ../klos/so2-X-mrcif12.f -g -w -o so2_pes

