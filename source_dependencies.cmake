set(numerical_recipies ${prefix}/src/base/numerical_recipies.f90)
set(external_sources ${numerical_recipies})

# Source file dependencies
set(grid_info ${prefix}/src/base/grid_info.f90)
set(constants ${prefix}/src/base/constants.f90)
set(coordinate_conversion ${prefix}/src/base/coordinate_conversion.f90)
set(rovib_utils_base ${prefix}/src/utils/rovib_utils_base.f90)
set(debug_tools_base ${prefix}/src/base/debug_tools_base.f90)

set(general_char_str ${prefix}/src/utils/general/general_char_str.F90)
set(general_integer ${general_char_str} ${prefix}/src/utils/general/general_integer.F90)
set(general_real ${general_char_str} ${prefix}/src/utils/general/general_real.F90)
set(general_real_array ${prefix}/src/utils/general/general_real_array.F90)
set(general_utils ${general_integer} ${general_real} ${general_real_array} ${general_char_str} ${prefix}/src/utils/general/general.F90)

set(path_utils ${general_utils} ${prefix}/src/utils/path.f90)
set(fourier_transform ${constants} ${general_utils} ${prefix}/src/utils/fourier_transform.f90)
set(formulas ${constants} ${general_utils} ${prefix}/src/base/formulas.f90)
set(parallel_utils ${algorithms} ${general_utils} ${prefix}/src/utils/parallel.f90)
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

set(algorithms_integer ${general_utils} ${vector} ${prefix}/src/utils/algorithms/algorithms_integer.F90)
set(algorithms_real ${general_utils} ${vector} ${prefix}/src/utils/algorithms/algorithms_real.F90)
set(algorithms ${algorithms_integer} ${algorithms_real} ${prefix}/src/utils/algorithms/algorithms.f90)

set(array_1d ${array_1d_real} ${array_1d_integer} ${array_1d_complex} ${prefix}/src/utils/array_1d/array_1d.F90)
set(array_2d_real ${prefix}/src/utils/array_2d/array_2d_real.F90)
set(array_2d_complex ${prefix}/src/utils/array_2d/array_2d_complex.F90)
set(array_2d ${array_2d_real} ${array_2d_complex} ${prefix}/src/utils/array_2d/array_2d.F90)

set(dict_integer ${prefix}/src/interface/fdict/dict_integer.F90)
set(dict_integer_array ${prefix}/src/interface/fdict/dict_integer_array.F90)
set(dict_char_str ${prefix}/src/interface/fdict/dict_char_str.F90)
set(dict_utils ${dict_integer} ${dict_integer_array} ${dict_char_str} ${string} ${vector_string} ${prefix}/src/interface/fdict/dict.F90)

set(string_utils ${vector_string} ${prefix}/src/utils/string/string_utils.f90)

set(io_base ${general_utils} ${prefix}/src/utils/io/io_base.f90)
set(io_integer ${io_base} ${string} ${string_utils} ${vector} ${prefix}/src/utils/io/io_integer.F90)
set(io_real ${io_base} ${string} ${string_utils} ${vector} ${prefix}/src/utils/io/io_real.F90)
set(io_complex ${io_base} ${string} ${string_utils} ${vector} ${prefix}/src/utils/io/io_complex.F90)
set(io_utils ${io_base} ${io_integer} ${io_real} ${io_complex} ${prefix}/src/utils/io/io.f90)

set(config ${dict_utils} ${general_utils} ${io_utils} ${parallel_utils} ${string} ${string_utils} ${prefix}/src/base/input_params/config.f90)
set(grid_params ${config} ${dict_utils} ${general_utils} ${string} ${prefix}/src/base/input_params/grid_params.f90)
set(optgrid_params ${config} ${dict_utils} ${grid_params} ${parallel_utils} ${string} ${prefix}/src/base/input_params/optgrid_params.f90)
set(wf_section_params ${config} ${constants} ${dict_utils} ${general_utils} ${string} ${string_utils} ${prefix}/src/base/input_params/wf_section_params.f90)
set(cap_params ${config} ${constants} ${dict_utils} ${general_utils} ${string} ${prefix}/src/base/input_params/cap_params.f90)
set(eigensolve_params ${config} ${constants} ${dict_utils} ${general_utils} ${parallel_utils} ${string} ${prefix}/src/base/input_params/eigensolve_params.f90)
set(fixed_basis_params ${config} ${constants} ${dict_utils} ${general_utils} ${string} ${prefix}/src/base/input_params/fixed_basis_params.f90)
set(basis_params ${config} ${constants} ${dict_utils} ${fixed_basis_params} ${general_utils} ${string} ${prefix}/src/base/input_params/basis_params.f90)
set(input_params ${algorithms} ${basis_params} ${cap_params} ${config} ${constants} ${dict_utils} ${eigensolve_params} ${general_utils} ${grid_info} ${grid_params} ${optgrid_params} ${parallel_utils} ${rovib_utils_base} ${string} ${string_utils} ${wf_section_params} ${prefix}/src/base/input_params/input_params.f90)

set(debug_tools ${debug_tools_base} ${general_utils} ${input_params} ${prefix}/src/base/debug_tools.f90)
set(rovib_utils ${input_params} ${rovib_utils_base} ${prefix}/src/utils/rovib_utils.f90)
set(spectrumsdt_paths ${general_utils} ${input_params} ${path_utils} ${rovib_utils_base} ${prefix}/src/base/spectrumsdt_paths.f90)
set(cap ${constants} ${formulas} ${grid_info} ${input_params} ${prefix}/src/base/cap.f90)
set(potential ${formulas} ${input_params} ${spectrumsdt_paths} ${prefix}/src/base/potential.f90)

set(block_borders ${prefix}/src/base/hamiltonian/block_borders.f90)
set(matrix_block_info ${block_borders} ${general_utils} ${prefix}/src/base/hamiltonian/matrix_block_info.f90)
set(rovib_io ${array_1d} ${array_2d} ${constants} ${io_utils} ${input_params} ${parallel_utils} ${rovib_utils} ${spectrumsdt_paths} ${prefix}/src/rovib/rovib_io.f90)
set(k_block_info ${algorithms} ${io_utils} ${matrix_block_info} ${rovib_io} ${prefix}/src/base/hamiltonian/k_block_info.f90)
set(overlaps_extra ${array_1d} ${array_2d} ${formulas} ${input_params} ${k_block_info} ${path_utils} ${parallel_utils} ${rovib_io} ${prefix}/src/rovib/overlaps_extra.f90)
set(distributed_rovib_hamiltonian ${block_borders} ${general_utils} ${formulas} ${input_params} ${k_block_info} ${matrix_block_info} ${parallel_utils} ${rovib_utils} ${spectrumsdt_paths} ${prefix}/src/base/hamiltonian/distributed_rovib_hamiltonian.f90)
set(matmul_operator ${distributed_rovib_hamiltonian} ${input_params} ${matrix_block_info} ${prefix}/src/base/eigensolve/matmul_operator.f90)
set(state_properties ${algorithms} ${array_1d} ${array_2d} ${constants} ${general_utils} ${grid_info} ${input_params} ${parallel_utils} ${path_utils} ${rovib_io} ${rovib_utils} ${vector} ${prefix}/src/base/state_properties.f90)
set(slepc_solver ${general_utils} ${matmul_operator} ${parallel_utils} ${prefix}/src/base/eigensolve/slepc_solver.F90)
set(lapack_interface ${general_utils} ${prefix}/src/interface/lapack.f90)
set(sdt ${algorithms} ${array_1d} ${array_2d} ${constants} ${formulas} ${fourier_transform} ${general_utils} ${input_params} ${lapack_interface} ${parallel_utils} ${rovib_io} ${spectrumsdt_paths} ${prefix}/src/base/sdt.f90)
set(spectrum ${cap} ${constants} ${formulas} ${general_utils} ${grid_info} ${input_params} ${io_utils} ${matmul_operator} ${parallel_utils} ${path_utils} ${sdt} ${slepc_solver} ${spectrumsdt_paths} ${prefix}/src/base/eigensolve/spectrum.f90)
set(grids ${config} ${constants} ${coordinate_conversion} ${formulas} ${general_utils} ${grid_info} ${input_params} ${io_utils} ${numerical_recipies} ${path_utils} ${spectrumsdt_paths} ${vector} ${prefix}/src/base/grids.f90)

set(spectrumsdt ${cap} ${config} ${debug_tools} ${dict_utils} ${formulas} ${input_params} ${grid_info} ${grids} ${overlaps_extra} ${parallel_utils} ${potential} ${sdt} ${spectrum} ${state_properties} ${prefix}/src/base/spectrumsdt.f90)

# Disable all warnings and standard conformation checks for external sources
set_source_files_properties(${external_sources} PROPERTIES COMPILE_FLAGS "-w -std=gnu")

add_executable(spectrumsdt ${spectrumsdt})
