!-------------------------------------------------------------------------------------------------------------------------------------------
! Module for printing wave functions on grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
module wf_print_mod
  use array_1d_mod
  use array_2d_mod
  use basis_mod
  use general_utils_mod
  use input_params_mod
  use io_utils_mod
  use iso_fortran_env, only: real64
  use overlaps_mod
  use spectrumsdt_io_mod
  use spectrumsdt_paths_mod
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts a given set of 1D solutions over FBR to solutions over grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_wf_grid_1d(params, solutions_1d, grid_phi) result(solutions_grid)
    class(input_params), intent(in) :: params
    complex(real64), intent(in) :: solutions_1d(:, :)
    real(real64), intent(in) :: grid_phi(:)
    real(real64), allocatable :: fbr_grid(:, :)
    complex(real64), allocatable :: solutions_grid(:, :)

    fbr_grid = get_phi_basis_grid(params, grid_phi)
    solutions_grid = matmul(fbr_grid, solutions_1d)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Reads stored FBR 1D solutions for given values of rho and theta, converts them to grid solutions and prints to a file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine reprint_wf_grid_1d(params, rho_ind, theta_ind, grid_phi)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind, theta_ind
    real(real64), intent(in) :: grid_phi(:)
    integer, allocatable :: num_solutions_1d(:)
    complex(real64), allocatable :: solutions_grid(:, :)
    character(:), allocatable :: sym_path, out_file_name
    type(array_1d_real), allocatable :: energies_1d(:)
    type(array_2d_complex), allocatable :: solutions_1d(:)

    sym_path = get_sym_path(params)
    call load_solutions_1D(get_solutions_1d_path(sym_path, rho_ind), num_solutions_1d, energies_1d, solutions_1d)
    solutions_grid = get_wf_grid_1d(params, solutions_1d(theta_ind) % p, grid_phi)
    out_file_name = 'wfs_1d_grid.' // num2str(rho_ind) // '.' // num2str(theta_ind)
    call write_matrix_complex_2(solutions_grid, out_file_name)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints specified 2D wave function on grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_wf_grid_2d(params, rho_ind, solution_ind, grid_phi)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind, solution_ind
    real(real64), intent(in) :: grid_phi(:)
    integer, allocatable :: num_solutions_1d(:)
    real(real64), allocatable :: energies_2d(:), selected_solution_fbr(:)
    real(real64), allocatable :: solutions_2d(:, :), shaped_selected_solution_fbr(:, :), fbr_grid(:, :), selected_solution_grid(:, :)
    character(:), allocatable :: sym_path, out_file_name
    type(array_1d_real), allocatable :: energies_1d(:)
    type(array_2d_real), allocatable :: solutions_1d(:)

    sym_path = get_sym_path(params)
    call load_solutions_1D(get_solutions_1d_path(sym_path, rho_ind), num_solutions_1d, energies_1d, solutions_1d)
    call load_solutions_2D(get_solutions_2d_path(sym_path, rho_ind), energies_2d, solutions_2d)
    selected_solution_fbr = transform_basis_1d_to_fbr(num_solutions_1d, solutions_1d, solutions_2d(:, solution_ind))
    shaped_selected_solution_fbr = reshape(selected_solution_fbr, [size(solutions_1d(1) % p, 1), size(solutions_1d)])
    fbr_grid = get_phi_basis_grid(params, grid_phi)
    selected_solution_grid = matmul(fbr_grid, shaped_selected_solution_fbr)

    out_file_name = 'wf_2d.' // num2str(rho_ind) // '.' // num2str(solution_ind)
    call write_matrix(selected_solution_grid, out_file_name)
  end subroutine

end module

