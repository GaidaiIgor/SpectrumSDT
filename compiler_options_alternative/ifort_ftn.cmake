# An example of compiler options for Intel's ifort (called through Cray's ftn wrapper)
set(CMAKE_Fortran_COMPILER ftn)
add_compile_options(-g -traceback -xHost -O3 -warn all -diag-disable 7712)
