DAWES_SRC=../dawes/src/*
mpifort -c -g ../dawes/src/path_resolution.f90
mpifort -c -g ../dawes/src/splin.f90
mpifort ../pesprint.f90 ${DAWES_SRC} -ffree-line-length-0 -g -w -o pesprint
