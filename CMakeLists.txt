cmake_minimum_required(VERSION 3.5)

set(CMAKE_Fortran_COMPILER gfortran)
add_compile_options(-g -fbacktrace -ffree-line-length-0 -Wall -Wextra -Wno-uninitialized -Wno-unused-dummy-argument -std=f2018)

# set(CMAKE_Fortran_COMPILER ifort)
# add_compile_options(-g -traceback -warn all -diag-disable 7712)

project(SpectrumSDT)
enable_language(Fortran)

set(prefix ${CMAKE_CURRENT_LIST_DIR})
# Specify where to look for external module files (-I option)
include_directories(${prefix}/src/include ${prefix}/libs/fdict ${prefix}/libs/petsc/include ${prefix}/libs/petsc/debug_64/include ${prefix}/libs/slepc/include ${prefix}/libs/slepc/debug_64/include)
# Specify where to look for libraries (-L option)
link_directories(${prefix}/libs/fdict ${prefix}/libs/petsc/debug_64/lib ${prefix}/libs/slepc/debug_64/lib)
include(${prefix}/source_dependencies.cmake)
# Link libraries (-llibname option)
target_link_libraries(spectrumsdt fdict blas lapack petsc slepc)
