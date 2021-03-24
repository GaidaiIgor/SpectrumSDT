#!/bin/sh
ftn -c ../../../pes_utils.f90 -g
ftn -c ../ozone_nr_pes_2013_hyps.f -g
ftn ../ozone_pes.f90 ozone_nr_pes_2013_hyps.o pes_utils.o -g -w -o ozone_pes
