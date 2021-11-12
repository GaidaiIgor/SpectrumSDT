set(numerical_recipies ${prefix}/src/base/numerical_recipies.f90)
set(external_sources ${numerical_recipies})

# Source file dependencies
set(grid_info ${prefix}/src/base/grid_info.f90)
set(constants ${prefix}/src/base/constants.f90)
set(coordinate_conversion ${prefix}/src/base/coordinate_conversion.f90)
set(rovib_utils_base ${prefix}/src/rovib/rovib_utils_base.f90)
set(debug_tools_base ${prefix}/src/base/debug/debug_tools_base.f90)
set(block_borders ${prefix}/src/base/hamiltonian/block_borders.f90)

set(general_utils_char_str ${prefix}/src/utils/general/general_utils_char_str.F90)
set(general_utils_integer ${general_utils_char_str} ${prefix}/src/utils/general/general_utils_integer.F90)
set(general_utils_real ${general_utils_char_str} ${prefix}/src/utils/general/general_utils_real.F90)
set(general_utils_real_array ${prefix}/src/utils/general/general_utils_real_array.F90)
set(general_utils ${general_utils_integer} ${general_utils_real} ${general_utils_real_array} ${general_utils_char_str} ${prefix}/src/utils/general/general_utils.f90)

set(path_utils ${general_utils} ${prefix}/src/utils/path_utils.f90)
set(fourier_transform ${constants} ${general_utils} ${prefix}/src/utils/fourier_transform.f90)
set(formulas ${constants} ${general_utils} ${prefix}/src/base/formulas.f90)
set(string ${general_utils} ${prefix}/src/utils/string/string.f90)
set(lapack_interface ${general_utils} ${prefix}/src/interface/lapack.f90)
set(matrix_block_info ${block_borders} ${general_utils} ${prefix}/src/base/hamiltonian/matrix_block_info.f90)

set(general_string ${string} ${prefix}/src/utils/general/general_string.F90)
set(general_utils_ext ${general_utils} ${general_string} ${prefix}/src/utils/general/general_utils_ext.f90)

set(parallel_utils_base ${algorithms} ${general_utils} ${prefix}/src/utils/parallel/parallel_utils_base.f90)
set(parallel_utils_integer ${parallel_utils_base} ${prefix}/src/utils/parallel/parallel_utils_integer.F90)
set(parallel_utils_real ${parallel_utils_base} ${prefix}/src/utils/parallel/parallel_utils_real.F90)
set(parallel_utils ${parallel_utils_integer} ${parallel_utils_real} ${prefix}/src/utils/parallel/parallel_utils.f90)

set(vector_integer ${general_utils} ${prefix}/src/utils/vector/vector_integer.F90)
set(vector_real ${general_utils} ${prefix}/src/utils/vector/vector_real.F90)
set(vector_complex ${general_utils} ${prefix}/src/utils/vector/vector_complex.F90)
set(vector_string ${general_utils} ${string} ${prefix}/src/utils/vector/vector_string.F90)
set(vector_vector_integer ${general_utils} ${vector_integer} ${prefix}/src/utils/vector/vector_vector_integer.F90)
set(vector_vector_real ${general_utils} ${vector_real} ${prefix}/src/utils/vector/vector_vector_real.F90)

set(array_1d_integer ${vector_integer} ${prefix}/src/utils/array_1d/array_1d_integer.F90)
set(array_1d_real ${vector_real} ${vector_integer} ${prefix}/src/utils/array_1d/array_1d_real.F90)
set(array_1d_complex ${vector_complex} ${vector_integer} ${prefix}/src/utils/array_1d/array_1d_complex.F90)
set(array_1d ${array_1d_integer} ${array_1d_real} ${array_1d_complex} ${prefix}/src/utils/array_1d/array_1d.f90)

set(vector_array_1d_integer ${general_utils} ${array_1d_integer} ${prefix}/src/utils/vector/vector_array_1d_integer.F90)
set(vector_array_1d_real ${general_utils} ${array_1d_real} ${prefix}/src/utils/vector/vector_array_1d_real.F90)
set(vector_array_1d_complex ${general_utils} ${array_1d_complex} ${prefix}/src/utils/vector/vector_array_1d_complex.F90)
set(vector ${vector_integer} ${vector_real} ${vector_complex} ${vector_string} ${vector_vector_integer} ${vector_vector_real} ${vector_array_1d_integer} ${vector_array_1d_real} ${vector_array_1d_complex} ${prefix}/src/utils/vector/vector.f90)

set(array_2d_real ${prefix}/src/utils/array_2d/array_2d_real.F90)
set(array_2d_complex ${prefix}/src/utils/array_2d/array_2d_complex.F90)
set(array_2d ${array_2d_real} ${array_2d_complex} ${prefix}/src/utils/array_2d/array_2d.f90)

set(array_3d_real ${prefix}/src/utils/array_3d/array_3d_real.F90)
set(array_3d ${array_3d_real} ${prefix}/src/utils/array_3d/array_3d.f90)

set(string_utils ${vector_string} ${prefix}/src/utils/string/string_utils.f90)

set(algorithms_base ${general_utils} ${vector} ${prefix}/src/utils/algorithms/algorithms_base.f90)
set(algorithms_integer ${algorithms_base} ${prefix}/src/utils/algorithms/algorithms_integer.F90)
set(algorithms_real ${algorithms_base} ${prefix}/src/utils/algorithms/algorithms_real.F90)
set(algorithms ${algorithms_integer} ${algorithms_real} ${prefix}/src/utils/algorithms/algorithms.f90)

set(dict_utils_base ${string} ${vector_string} ${prefix}/src/interface/fdict/dict_utils_base.f90)
set(dict_utils_integer ${dict_utils_base} ${prefix}/src/interface/fdict/dict_utils_integer.F90)
set(dict_utils_integer_array ${dict_utils_base} ${prefix}/src/interface/fdict/dict_utils_integer_array.F90)
set(dict_utils_char_str ${dict_utils_base} ${prefix}/src/interface/fdict/dict_utils_char_str.F90)
set(dict_utils ${dict_utils_integer} ${dict_utils_integer_array} ${dict_utils_char_str} ${prefix}/src/interface/fdict/dict_utils.f90)

set(io_utils_base ${general_utils} ${string} ${string_utils} ${vector} ${prefix}/src/utils/io/io_utils_base.f90)
set(io_utils_integer ${io_utils_base} ${prefix}/src/utils/io/io_utils_integer.F90)
set(io_utils_real ${io_utils_base} ${prefix}/src/utils/io/io_utils_real.F90)
set(io_utils_complex ${io_utils_base} ${io_utils_real} ${prefix}/src/utils/io/io_utils_complex.F90)
set(io_utils ${io_utils_integer} ${io_utils_real} ${io_utils_complex} ${prefix}/src/utils/io/io_utils.f90)

set(config ${dict_utils} ${general_utils} ${io_utils} ${parallel_utils} ${string} ${string_utils} ${prefix}/src/base/input_params/config.f90)
set(grid_params ${config} ${dict_utils} ${general_utils} ${string} ${prefix}/src/base/input_params/grid_params.f90)
set(optgrid_params ${config} ${dict_utils} ${grid_params} ${parallel_utils} ${string} ${prefix}/src/base/input_params/optgrid_params.f90)
set(wf_section_params ${config} ${constants} ${dict_utils} ${general_utils} ${string} ${string_utils} ${prefix}/src/base/input_params/wf_section_params.f90)
set(cap_params ${config} ${constants} ${dict_utils} ${general_utils} ${string} ${prefix}/src/base/input_params/cap_params.f90)
set(eigensolve_params ${config} ${constants} ${dict_utils} ${general_utils} ${parallel_utils} ${string} ${prefix}/src/base/input_params/eigensolve_params.f90)
set(fixed_basis_params ${config} ${constants} ${dict_utils} ${general_utils} ${string} ${prefix}/src/base/input_params/fixed_basis_params.f90)
set(basis_params ${config} ${constants} ${dict_utils} ${fixed_basis_params} ${general_utils} ${string} ${prefix}/src/base/input_params/basis_params.f90)
set(debug_params ${config} ${dict_utils} ${general_utils} ${string} ${string_utils} ${prefix}/src/base/input_params/debug_params.f90)
set(spectrumsdt_utils ${constants} ${general_utils} ${prefix}/src/base/spectrumsdt_utils.f90)
set(input_params ${algorithms} ${basis_params} ${cap_params} ${config} ${constants} ${debug_params} ${dict_utils} ${eigensolve_params} ${general_utils_ext} ${grid_info} ${grid_params} ${optgrid_params} ${parallel_utils} ${rovib_utils_base} ${spectrumsdt_utils} ${string} ${string_utils} ${wf_section_params} ${prefix}/src/base/input_params/input_params.f90)

set(spectrumsdt_paths ${general_utils} ${input_params} ${path_utils} ${rovib_utils} ${prefix}/src/base/spectrumsdt_paths.f90)
set(debug_tools ${debug_tools_base} ${input_params} ${io_utils} ${prefix}/src/base/debug/debug_tools.f90)
set(rovib_utils ${general_utils} ${input_params} ${rovib_utils_base} ${prefix}/src/rovib/rovib_utils.f90)
set(spectrumsdt_utils_ext ${input_params} ${rovib_utils} ${spectrumsdt_utils} ${prefix}/src/base/spectrumsdt_utils_ext.f90)

set(cap ${constants} ${grid_info} ${input_params} ${spectrumsdt_utils} ${prefix}/src/base/cap.f90)
set(potential ${input_params} ${spectrumsdt_paths} ${spectrumsdt_utils} ${prefix}/src/base/potential.f90)

set(spectrumsdt_io_base ${array_1d} ${array_2d} ${constants} ${general_utils} ${io_utils} ${input_params} ${parallel_utils} ${rovib_utils} ${spectrumsdt_paths} ${spectrumsdt_utils} ${prefix}/src/base/io/spectrumsdt_io_base.f90)
set(spectrumsdt_io_real ${spectrumsdt_io_base} ${prefix}/src/base/io/spectrumsdt_io_real.F90)
set(spectrumsdt_io_complex ${spectrumsdt_io_base} ${prefix}/src/base/io/spectrumsdt_io_complex.F90)
set(spectrumsdt_io ${spectrumsdt_io_real} ${spectrumsdt_io_complex} ${prefix}/src/base/io/spectrumsdt_io.f90)

set(k_block_info ${algorithms} ${io_utils} ${matrix_block_info} ${spectrumsdt_io} ${prefix}/src/base/hamiltonian/k_block_info.f90)

set(grids ${config} ${constants} ${coordinate_conversion} ${general_utils} ${grid_info} ${input_params} ${io_utils} ${numerical_recipies} ${path_utils} ${spectrumsdt_paths} ${spectrumsdt_utils} ${vector} ${prefix}/src/base/grids.f90)

set(basis_base ${algorithms} ${array_1d} ${array_2d} ${constants} ${fourier_transform} ${general_utils} ${input_params} ${lapack_interface} ${parallel_utils} ${spectrumsdt_paths} ${spectrumsdt_utils_ext} ${prefix}/src/base/basis/basis_base.f90)
set(basis_real ${basis_base} ${prefix}/src/base/basis/basis_real.F90)
set(basis_complex ${basis_base} ${prefix}/src/base/basis/basis_complex.F90)
set(basis ${basis_real} ${basis_complex} ${prefix}/src/base/basis/basis.f90)

set(overlaps_base ${array_1d} ${array_2d} ${constants} ${general_utils} ${input_params} ${parallel_utils} ${spectrumsdt_io} ${spectrumsdt_paths} ${prefix}/src/base/overlaps/overlaps_base.f90)
set(overlaps_real ${overlaps_base} ${prefix}/src/base/overlaps/overlaps_real.F90)
set(overlaps_complex ${overlaps_base} ${prefix}/src/base/overlaps/overlaps_complex.F90)
set(overlaps ${overlaps_real} ${overlaps_complex} ${prefix}/src/base/overlaps/overlaps.f90)

set(wf_print_base ${array_1d} ${array_2d} ${basis} ${constants} ${general_utils} ${input_params} ${io_utils} ${overlaps} ${spectrumsdt_io} ${spectrumsdt_paths} ${prefix}/src/base/debug/wf_print/wf_print_base.f90)
set(wf_print_real ${wf_print_base} ${prefix}/src/base/debug/wf_print/wf_print_real.F90)
set(wf_print_complex ${wf_print_base} ${prefix}/src/base/debug/wf_print/wf_print_complex.F90)
set(wf_print ${wf_print_real} ${wf_print_complex} ${prefix}/src/base/debug/wf_print/wf_print.f90)

set(overlaps_extra_base ${array_1d} ${array_2d} ${input_params} ${k_block_info} ${overlaps} ${parallel_utils} ${path_utils} ${spectrumsdt_io} ${spectrumsdt_utils} ${prefix}/src/rovib/overlaps_extra/overlaps_extra_base.f90)
set(overlaps_extra_real ${overlaps_extra_base} ${prefix}/src/rovib/overlaps_extra/overlaps_extra_real.F90)
set(overlaps_extra_complex ${overlaps_extra_base} ${prefix}/src/rovib/overlaps_extra/overlaps_extra_complex.F90)
set(overlaps_extra ${overlaps_extra_real} ${overlaps_extra_complex} ${prefix}/src/rovib/overlaps_extra/overlaps_extra.f90)

set(distributed_rovib_hamiltonian ${block_borders} ${general_utils} ${input_params} ${k_block_info} ${matrix_block_info} ${parallel_utils} ${rovib_utils} ${spectrumsdt_io} ${spectrumsdt_paths} ${prefix}/src/base/hamiltonian/distributed_rovib_hamiltonian.f90)
set(matmul_operator ${distributed_rovib_hamiltonian} ${input_params} ${matrix_block_info} ${prefix}/src/base/eigensolve/matmul_operator.f90)
set(slepc_solver ${general_utils} ${matmul_operator} ${parallel_utils} ${prefix}/src/base/eigensolve/slepc_solver.F90)
set(eigensolve ${basis} ${cap} ${constants} ${general_utils} ${grid_info} ${input_params} ${io_utils} ${matmul_operator} ${parallel_utils} ${path_utils} ${slepc_solver} ${spectrumsdt_paths} ${spectrumsdt_utils} ${prefix}/src/base/eigensolve/eigensolve.f90)

set(properties_base ${algorithms} ${array_1d} ${array_2d} ${array_3d} ${constants} ${general_utils} ${grid_info} ${input_params} ${parallel_utils} ${path_utils} ${spectrumsdt_io} ${spectrumsdt_utils} ${rovib_utils} ${vector} ${prefix}/src/base/properties/properties_base.f90)
set(properties_real ${properties_base} ${prefix}/src/base/properties/properties_real.F90)
set(properties_complex ${properties_base} ${prefix}/src/base/properties/properties_complex.F90)
set(properties ${properties_real} ${properties_complex} ${prefix}/src/base/properties/properties.f90)

set(spectrumsdt ${basis} ${cap} ${config} ${debug_tools} ${dict_utils} ${eigensolve} ${input_params} ${grid_info} ${grids} ${overlaps} ${overlaps_extra} ${parallel_utils} ${potential} ${properties} ${spectrumsdt_utils} ${wf_print} ${prefix}/src/base/spectrumsdt.f90)

# Disable all warnings and standard conformation checks for external sources
set_source_files_properties(${external_sources} PROPERTIES COMPILE_FLAGS "-w -std=gnu")

add_executable(spectrumsdt ${spectrumsdt})
