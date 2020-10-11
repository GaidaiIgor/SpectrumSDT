cmake_minimum_required(VERSION 3.5)
set(prefix ${CMAKE_CURRENT_LIST_DIR})

if (EXISTS ${prefix}/compiler_options.cmake)
  include(${prefix}/compiler_options.cmake)
else()
  include(${prefix}/compiler_options_default.cmake)
endif()

project(SpectrumSDT)
enable_language(Fortran)

# Specify where to look for external module files (-I option)
include_directories(${prefix}/include ${prefix}/src/include ${prefix}/libs/fdict ${prefix}/libs/ftl/build/include/ftl ${prefix}/libs/petsc/include ${prefix}/libs/petsc/debug_64/include ${prefix}/libs/slepc/include ${prefix}/libs/slepc/debug_64/include)
# Specify where to look for libraries (-L option)
link_directories(${prefix}/libs/fdict ${prefix}/libs/ftl/build/lib ${prefix}/libs/petsc/debug_64/lib ${prefix}/libs/slepc/debug_64/lib)

# External sources. Dependency of Dawes' pes on path_utils is factored into pesinterface to avoid compiling path_utils with Dawes' flags
set(dawes_pes ${prefix}/PES/dawes/src/PES_3D_subroutine_ozone_L3.f90 ${prefix}/PES/dawes/src/splin.f90 ${prefix}/PES/dawes/src/dynamic_parameters_3D.f90 ${prefix}/PES/dawes/src/2D-spline.f ${prefix}/PES/dawes/src/VLGRoo.f ${prefix}/PES/dawes/src/CJe.f ${prefix}/PES/dawes/src/C6Js.f ${prefix}/PES/dawes/src/cleb.f ${prefix}/PES/dawes/src/Celr.f ${prefix}/PES/dawes/src/C6r.f ${prefix}/PES/dawes/src/faclog.f ${prefix}/PES/dawes/src/index.f90)
set(numerical_recipies ${prefix}/src/tools/numerical_recipies.f90)
set(external_sources ${dawes_pes} ${numerical_recipies})

# Source file dependencies
set(constants ${prefix}/src/utils/constants.f90)
set(input_params ${prefix}/src/input_params/input_params.f90)

set(general_char_str ${prefix}/src/utils/general/general_char_str.F90)
set(general_integer ${general_char_str} ${prefix}/src/utils/general/general_integer.F90)
set(general_real ${general_char_str} ${prefix}/src/utils/general/general_real.F90)
set(general_real_array ${prefix}/src/utils/general/general_real_array.F90)
set(general_utils ${general_integer} ${general_real} ${general_real_array} ${general_char_str} ${prefix}/src/utils/general/general.F90)

set(debug_tools ${general_utils} ${input_params} ${prefix}/src/tools/debug_tools.f90)
set(fourier_transform ${constants} ${general_utils} ${prefix}/src/utils/fourier_transform.f90)
set(formulas ${constants} ${general_utils} ${prefix}/src/tools/formulas.f90)
set(rovib_utils ${input_params} ${prefix}/src/utils/rovib_utils.f90)
set(algorithms ${general_utils} ${prefix}/src/utils/algorithms.f90)
set(parallel_utils ${algorithms} ${general_utils} ${prefix}/src/utils/parallel.f90)
set(cap ${constants} ${general_vars} ${input_params} ${prefix}/src/tools/cap.f90)
set(string ${general_utils} ${prefix}/src/utils/string/string.f90)

set(vector_integer ${general_utils} ${prefix}/src/utils/vector/vector_integer.F90)
set(vector_real ${general_utils} ${prefix}/src/utils/vector/vector_real.F90)
set(vector_complex ${general_utils} ${prefix}/src/utils/vector/vector_complex.F90)
set(vector_string ${general_utils} ${string} ${prefix}/src/utils/vector/vector_string.F90)
set(vector_vector_integer ${general_utils} ${vector_integer} ${prefix}/src/utils/vector/vector_vector_integer.F90)
set(vector_vector_real ${general_utils} ${vector_real} ${prefix}/src/utils/vector/vector_vector_real.F90)
set(array_1d_integer ${vector_integer} ${prefix}/src/utils/array_1d/array_1d_integer.F90)
set(array_1d_real ${vector_real} ${vector_integer} ${prefix}/src/utils/array_1d/array_1d_real.F90)
set(array_1d_complex ${vector_complex} ${vector_integer} ${prefix}/src/utils/array_1d/array_1d_complex.F90)
set(vector_array_1d_integer ${general_utils} ${array_1d_integer} ${prefix}/src/utils/vector/vector_array_1d_integer.F90)
set(vector_array_1d_real ${general_utils} ${array_1d_real} ${prefix}/src/utils/vector/vector_array_1d_real.F90)
set(vector_array_1d_complex ${general_utils} ${array_1d_complex} ${prefix}/src/utils/vector/vector_array_1d_complex.F90)
set(vector ${vector_integer} ${vector_real} ${vector_complex} ${vector_string} ${vector_vector_integer} ${vector_vector_real} ${vector_array_1d_integer} ${vector_array_1d_real} ${vector_array_1d_complex} ${prefix}/src/utils/vector/vector.f90)

set(array_1d ${array_1d_real} ${array_1d_integer} ${array_1d_complex} ${prefix}/src/utils/array_1d/array_1d.F90)
set(array_2d_real ${prefix}/src/utils/array_2d/array_2d_real.F90)
set(array_2d_complex ${prefix}/src/utils/array_2d/array_2d_complex.F90)
set(array_2d ${array_2d_real} ${array_2d_complex} ${prefix}/src/utils/array_2d/array_2d.F90)

set(dict_integer ${prefix}/src/utils/dict/dict_integer.F90)
set(dict_integer_array ${prefix}/src/utils/dict/dict_integer_array.F90)
set(dict_char_str ${prefix}/src/utils/dict/dict_char_str.F90)
set(dict_utils ${dict_integer} ${dict_integer_array} ${dict_char_str} ${string} ${vector_string} ${prefix}/src/utils/dict/dict.F90)

set(stage_metadata ${dict_utils} ${string} ${prefix}/src/input_params/stage_metadata.F90)
set(string_utils ${vector_string} ${prefix}/src/utils/string/string_utils.f90)

set(io_base ${general_utils} ${prefix}/src/utils/io/io_base.f90)
set(io_integer ${io_base} ${string} ${string_utils} ${vector} ${prefix}/src/utils/io/io_integer.F90)
set(io_real ${io_base} ${string} ${string_utils} ${vector} ${prefix}/src/utils/io/io_real.F90)
set(io_complex ${io_base} ${string} ${string_utils} ${vector} ${prefix}/src/utils/io/io_complex.F90)
set(io_utils ${io_base} ${io_integer} ${io_real} ${io_complex} ${prefix}/src/utils/io/io.f90)

set(system ${general_utils} ${io_utils} ${prefix}/src/utils/system.f90)
set(path_utils ${general_utils} ${parallel_utils} ${system} ${prefix}/src/utils/path.f90)
set(spectrumsdt_paths ${general_utils} ${input_params} ${path_utils} ${prefix}/src/base/spectrumsdt_paths.F90)
set(general_vars ${constants} ${input_params} ${path_utils} ${spectrumsdt_paths} ${prefix}/src/tools/general_vars.f90)
set(config ${constants} ${dict_utils} ${general_utils} ${input_params} ${io_utils} ${parallel_utils} ${rovib_utils} ${string} ${stage_metadata} ${string_utils} ${prefix}/src/input_params/config.f90)
set(potential ${formulas} ${input_params} ${general_vars} ${spectrumsdt_paths} ${prefix}/src/base/potential.F90)
set(pesinterface ${constants} ${dawes_pes} ${path_utils} ${parallel_utils} ${prefix}/src/tools/pes/pesinterface.f90)

set(block_borders ${prefix}/src/base/hamiltonian/block_borders.F90)
set(matrix_block_info ${block_borders} ${general_utils} ${prefix}/src/base/hamiltonian/matrix_block_info.F90)
set(k_block_info ${algorithms} ${io_utils} ${matrix_block_info} ${prefix}/src/base/hamiltonian/k_block_info.F90)
set(rovib_io ${array_1d} ${array_2d} ${constants} ${io_utils} ${input_params} ${parallel_utils} ${rovib_utils} ${spectrumsdt_paths} ${string} ${prefix}/src/rovib/rovib_io.f90)
set(overlaps_extra ${array_1d} ${array_2d} ${formulas} ${input_params} ${k_block_info} ${path_utils} ${parallel_utils} ${rovib_io} ${prefix}/src/rovib/overlaps_extra.f90)
set(distributed_rovib_hamiltonian ${block_borders} ${general_utils} ${formulas} ${input_params} ${k_block_info} ${matrix_block_info} ${parallel_utils} ${rovib_utils} ${spectrumsdt_paths} ${prefix}/src/base/hamiltonian/distributed_rovib_hamiltonian.F90)
set(matmul_operator ${distributed_rovib_hamiltonian} ${input_params} ${matrix_block_info} ${prefix}/src/base/eigencalc/matmul_operator.F90)
set(state_properties ${algorithms} ${array_1d} ${array_2d} ${constants} ${general_utils} ${input_params} ${parallel_utils} ${path_utils} ${rovib_io} ${rovib_utils} ${prefix}/src/base/state_properties.f90)
set(slepc_solver ${general_utils} ${matmul_operator} ${parallel_utils} ${prefix}/src/base/eigencalc/slepc_solver.F90)
set(lapack_interface ${general_utils} ${prefix}/src/interface/lapack.f90)
set(sdt ${constants} ${fourier_transform} ${general_vars} ${input_params} ${lapack_interface} ${parallel_utils} ${potential} ${spectrumsdt_paths} ${prefix}/src/base/sdt.f90)
set(spectrum ${cap} ${constants} ${general_utils} ${general_vars} ${input_params} ${matmul_operator} ${parallel_utils} ${path_utils} ${sdt} ${slepc_solver} ${spectrumsdt_paths} ${prefix}/src/base/eigencalc/spectrum.F90)
set(grids ${config} ${constants} ${general_utils} ${input_params} ${numerical_recipies} ${path_utils} ${vector} ${prefix}/src/base/grids.F90)

set(pesprint ${parallel_utils} ${pesinterface} ${prefix}/src/programs/pesprint.f90)

set(spectrumsdt ${cap} ${config} ${debug_tools} ${general_vars} ${input_params} ${grids} ${overlaps_extra} ${parallel_utils} ${potential} ${sdt} ${spectrum} ${state_properties} ${prefix}/src/programs/spectrumsdt.F90)

# Disable all warnings and standard conformation checks for external sources
set_source_files_properties(${external_sources} PROPERTIES COMPILE_FLAGS "-w -std=gnu")

add_executable(pesprint ${pesprint})
target_link_libraries(pesprint fdict ftl blas lapack)

add_executable(spectrumsdt ${spectrumsdt})
target_link_libraries(spectrumsdt fdict ftl blas lapack petsc slepc)
