# An example of compiler options for Intel's ifort (called through Cray's ftn wrapper)
set(CMAKE_Fortran_COMPILER ftn)
add_compile_options(-g -traceback)
# Makes MKL implementation of blas and lapack available for linker
set(CMAKE_EXE_LINKER_FLAGS "-mkl=cluster")
