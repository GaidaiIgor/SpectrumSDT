cmake_minimum_required(VERSION 3.5)
enable_language(Fortran)

# Define sets of source files.
set(general_utils src/utils/general/general.F90)
set(constants src/utils/constants.f90)
set(general_vars src/tools/general_vars.f90)
set(debug_tools src/tools/debug_tools.f90)
set(input_params src/tools/input_params.f90)
set(rovib_utils ${input_params} src/utils/rovib_utils.f90)
set(mkl src/tools/deriv/mkl.f90)
set(index_conversion src/tools/index_conversion.f90)
set(parabola src/tools/parabola.f90)

set(algorithms ${general_utils} src/utils/algorithms/algorithms.f90)
set(parallel_utils ${algorithms} ${general_utils} src/utils/parallel/parallel.f90)
set(path_utils ${general_utils} ${input_params} src/utils/path/path.f90)
set(string ${general_utils} src/utils/string/string.f90)

set(dawes_dynamic_parameters external/pes/dawes/src/dynamic_parameters_3D.f90)
set(dawes_linkables external/pes/dawes/src/2D-spline.f external/pes/dawes/src/VLGRoo.f external/pes/dawes/src/CJe.f external/pes/dawes/src/C6Js.f external/pes/dawes/src/cleb.f external/pes/dawes/src/Celr.f external/pes/dawes/src/C6r.f external/pes/dawes/src/faclog.f external/pes/dawes/src/index.f90)
set(dawes_ssplin ${path_utils} external/pes/dawes/src/splin.f90)
set(dawes ${dawes_linkables} ${path_utils} ${dawes_ssplin} ${dawes_dynamic_parameters} external/pes/dawes/src/PES_3D_subroutine_ozone_L3.f90)

set(vector_integer ${general_utils} src/utils/vector/vector_integer.F90)
set(vector_real ${general_utils} src/utils/vector/vector_real.F90)
set(vector_complex ${general_utils} src/utils/vector/vector_complex.F90)
set(vector_string ${general_utils} ${string} src/utils/vector/vector_string.F90)
set(vector_vector_integer ${general_utils} ${vector_integer} src/utils/vector/vector_vector_integer.F90)
set(vector_vector_real ${general_utils} ${vector_real} src/utils/vector/vector_vector_real.F90)
set(array_1d_integer ${vector_integer} src/utils/array_1d/array_1d_integer.F90)
set(array_1d_real ${vector_real} ${vector_integer} src/utils/array_1d/array_1d_real.F90)
set(array_1d_complex ${vector_complex} ${vector_integer} src/utils/array_1d/array_1d_complex.F90)
set(vector_array_1d_integer ${general_utils} ${array_1d_integer} src/utils/vector/vector_array_1d_integer.F90)
set(vector_array_1d_real ${general_utils} ${array_1d_real} src/utils/vector/vector_array_1d_real.F90)
set(vector_array_1d_complex ${general_utils} ${array_1d_complex} src/utils/vector/vector_array_1d_complex.F90)
set(vector ${vector_integer} ${vector_real} ${vector_complex} ${vector_string} ${vector_vector_integer} ${vector_vector_real} ${vector_array_1d_integer} ${vector_array_1d_real} ${vector_array_1d_complex} src/utils/vector/vector.f90)

set(array_1d ${array_1d_real} ${array_1d_integer} ${array_1d_complex} src/utils/array_1d/array_1d.F90)
set(array_2d_real src/utils/array_2d/array_2d_real.F90)
set(array_2d_complex src/utils/array_2d/array_2d_complex.F90)
set(array_2d ${array_2d_real} ${array_2d_complex} src/utils/array_2d/array_2d.F90)

set(dict_integer src/utils/dict/dict_integer.F90)
set(dict_integer_array src/utils/dict/dict_integer_array.F90)
set(dict_char_str src/utils/dict/dict_char_str.F90)
set(dict_utils ${dict_integer} ${dict_integer_array} ${dict_char_str} ${string} ${vector_string} src/utils/dict/dict.F90)

set(string_utils ${vector_string} src/utils/string/string_utils.f90)

set(io_base ${general_utils} src/utils/io/io_base.f90)
set(io_integer ${io_base} ${string} ${string_utils} ${vector} src/utils/io/io_integer.F90)
set(io_real ${io_base} ${string} ${string_utils} ${vector} src/utils/io/io_real.F90)
set(io_complex ${io_base} ${string} ${string_utils} ${vector} src/utils/io/io_complex.F90)
set(io_utils ${io_base} ${io_integer} ${io_real} ${io_complex} src/utils/io/io.f90)

set(formulas ${general_utils} ${constants} src/tools/formulas.f90)
set(config ${general_utils} ${parallel_utils} ${path_utils} ${constants} ${rovib_utils} ${input_params} ${dict_utils} ${string_utils} ${string} ${io_utils} src/tools/config.f90)
set(parpack ${parallel_utils} ${algorithms} ${index_conversion} src/tools/eicalc/parpack.f90)
set(scalapack ${algorithms} ${index_conversion} src/tools/eicalc/scalapack.f90)
set(pesgeneral ${formulas} ${general_utils} ${input_params} ${constants} src/tools/pes/pesgeneral.f90)
set(pes ${input_params} ${constants} ${pesgeneral} src/tools/pes/pesinterface.f90)
set(sdt ${debug_tools} ${formulas} ${arnoldi_operator} ${io_utils} ${parabola} ${scalapack} ${algorithms} ${constant_params} ${block_info} ${path} ${mkl} ${parpack} ${pesgeneral} ${general_vars} ${constants} src/tools/sdt.f90)

set(block_borders src/rovib/block_borders.f90)
set(matrix_block_info ${block_borders} src/rovib/matrix_block_info.f90)
set(k_block_info ${matrix_block_info} ${io_utils} ${algorithms} src/rovib/k_block_info.f90)
set(rovib_io ${array_1d} ${array_2d} ${io_utils} src/rovib/rovib_io.f90)
set(overlaps_extra ${rovib_io} ${formulas} ${array_2d} ${array_1d} ${path_utils} ${parallel_utils} ${input_params} src/rovib/overlaps_extra.f90)
set(distributed_rovib_hamiltonian ${block_borders} ${rovib_utils} ${formulas} ${general_utils} ${input_params} ${k_block_info} ${matrix_block_info} ${parallel_utils} ${parpack} ${path_utils} ${debug_tools} src/rovib/distributed_rovib_hamiltonian.f90)
set(state_properties ${array_1d} ${array_2d} ${constants} ${rovib_io} ${rovib_utils} ${input_params} ${general_utils} ${parallel_utils} ${path_utils} src/rovib/state_properties.f90)
set(arnoldi_operator ${distributed_rovib_hamiltonian} ${matrix_block_info} src/tools/eicalc/arnoldi_operator.f90)
set(slepc_solver ${arnoldi_operator} ${parallel_utils} ${debug_tools} ${constants} src/tools/eicalc/slepc_solver.F90)

set(numerical_recipies src/tools/numerical_recipies.f90)
set(optgrid_tools ${general_utils} ${numerical_recipies} src/tools/optgrid_tools.f90)

# Specify where to look for external module files (-I option)
include_directories(external/extras/mklmods external/include external/fdict/include src/include external/ftl_lib/include/ftl external/petsc/include external/petsc/arch-linux-c-debug/include external/slepc-3.12.2/include external/slepc-3.12.2/arch-linux-c-debug/include)
# Specify where to look for libraries (-L option)
link_directories(external/lib external/fdict/lib external/ftl_lib/lib external/petsc/arch-linux-c-debug/lib external/slepc-3.12.2/arch-linux-c-debug/lib)
# Add common compile flag for all sources. -u forces implicitnone for all functions.
add_compile_options(-mkl=cluster -assume bscc -g -traceback -u)
# Add common linker flags for all targets
set(CMAKE_EXE_LINKER_FLAGS "-dynamic -mkl=cluster")

add_executable(optgrid ${constants} ${optgrid_tools} src/programs/optgrid.f90)
target_link_libraries(optgrid parpack arpack mkl_blas95_lp64 mkl_lapack95_lp64)

add_executable(pesprint ${config} ${input_params} ${pes} ${dawes} src/programs/pesprint.f90)
target_link_libraries(pesprint fdict ftl pcreposix parpack arpack mkl_blas95_lp64 mkl_lapack95_lp64)

add_executable(spectrumsdt ${arnoldi_operator} ${config} ${constants} ${debug_tools} ${dict_utils} ${distributed_rovib_hamiltonian} ${general_vars} ${index_conversion} ${input_params} ${io_utils} ${overlaps_extra} ${parallel_utils} ${path_utils} ${pesgeneral} ${sdt} ${slepc_solver} ${state_properties} src/programs/spectrumsdt.f90)
target_link_libraries(spectrumsdt parpack arpack fdict ftl mkl_blas95_lp64 mkl_lapack95_lp64 pcreposix petsc slepc)

add_executable(kinetics ${pesgeneral} src/programs/kinetics.f90)
# target_link_libraries(pesprint fdict parpack arpack mkl_blas95_lp64 mkl_lapack95_lp64)

add_executable(test_prg ${parallel_utils} src/test/test_prg.F90)
target_link_libraries(test_prg fdict ftl parpack arpack mkl_blas95_lp64 mkl_lapack95_lp64 petsc slepc)

add_executable(dawes_de ${utils} ${tools} ${dawes} src/programs/misc/dawes_de.f90)
target_link_libraries(dawes_de parpack arpack mkl_blas95_lp64 mkl_lapack95_lp64)

# add_executable(zpeO2 ${utils} ${tools} ${dawes} src/programs/zpeO2.f)
# target_link_libraries(zpeO2 parpack arpack mkl_blas95_lp64 mkl_lapack95_lp64)
