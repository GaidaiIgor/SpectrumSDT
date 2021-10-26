#include "funcs.macro"
use wf_print_base_mod

  interface convert_matrix_fbr_to_grid
    module procedure :: CONCAT2(convert_matrix_fbr_to_grid_,TEMPLATE_TYPE_NAME)
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Changes the basis of matrix columns from FBR to grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(convert_matrix_fbr_to_grid_,TEMPLATE_TYPE_NAME)(params, fbr_columns, grid_phi) result(grid_columns)
    class(input_params), intent(in) :: params
    TEMPLATE_TYPE, intent(in) :: fbr_columns(:, :)
    real(real64), intent(in) :: grid_phi(:)
    real(real64), allocatable :: fbr_to_grid(:, :)
    TEMPLATE_TYPE, allocatable :: grid_columns(:, :)

    fbr_to_grid = get_phi_basis_grid(params, grid_phi)
    grid_columns = matmul(fbr_to_grid, fbr_columns)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Reads stored energies and FBR 1D solutions for given values of rho and theta, converts them to grid solutions and prints to a file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(reprint_wf_grid_1d_,TEMPLATE_TYPE_NAME)(params, rho_ind, theta_ind, grid_phi)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind, theta_ind
    real(real64), intent(in) :: grid_phi(:)
    integer, allocatable :: num_solutions_1d(:)
    TEMPLATE_TYPE, allocatable :: solutions_grid(:, :)
    character(:), allocatable :: sym_path, out_file_name
    type(array_1d_real), allocatable :: energies_1d(:)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable :: solutions_1d(:)

    sym_path = get_sym_path(params)
    call load_solutions_1D(get_solutions_1d_path(sym_path, rho_ind), num_solutions_1d, energies_1d, solutions_1d)

    call write_array(energies_1d(theta_ind) % p * au_to_wn, 'vals')
    solutions_grid = convert_matrix_fbr_to_grid(params, solutions_1d(theta_ind) % p, grid_phi)

    ! out_file_name = 'wfs_1d_grid.' // num2str(rho_ind) // '.' // num2str(theta_ind)
    out_file_name = 'wfs_grid.fwc'
#if TYPE_ID == COMPLEX_ID
    call write_matrix_complex_2(solutions_grid, out_file_name)
#else
    call write_matrix(solutions_grid, out_file_name)
#endif
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Reads stored 2D solutions at given *rho_ind*, converts *solution_ind* to grid representation and prints it to a file.
! File format: PxT matrix, where P is the size of phi grid, and T is the size of theta grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(reprint_wf_grid_2d_,TEMPLATE_TYPE_NAME)(params, rho_ind, solution_ind, grid_phi)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind, solution_ind
    real(real64), intent(in) :: grid_phi(:)
    integer, allocatable :: num_solutions_1d(:)
    real(real64), allocatable :: energies_2d(:)
    TEMPLATE_TYPE, allocatable :: selected_solution_fbr(:)
    TEMPLATE_TYPE, allocatable :: solutions_2d(:, :), shaped_selected_solution_fbr(:, :), selected_solution_grid(:, :)
    character(:), allocatable :: sym_path, out_file_name
    type(array_1d_real), allocatable :: energies_1d(:)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable :: solutions_1d(:)

    sym_path = get_sym_path(params)
    call load_solutions_1D(get_solutions_1d_path(sym_path, rho_ind), num_solutions_1d, energies_1d, solutions_1d)
    call load_solutions_2D(get_solutions_2d_path(sym_path, rho_ind), energies_2d, solutions_2d)
    call write_array(energies_2d * au_to_wn, 'vals_2d')

    selected_solution_fbr = transform_basis_1d_to_fbr(num_solutions_1d, solutions_1d, solutions_2d(:, solution_ind))
    shaped_selected_solution_fbr = reshape(selected_solution_fbr, [size(solutions_1d(1) % p, 1), size(solutions_1d)])
    selected_solution_grid = convert_matrix_fbr_to_grid(params, shaped_selected_solution_fbr, grid_phi)

    out_file_name = 'wf_2d.' // num2str(solution_ind)
#if TYPE_ID == COMPLEX_ID
    call write_matrix_complex_2(selected_solution_grid, out_file_name)
#else
    call write_matrix(selected_solution_grid, out_file_name)
#endif
  end subroutine

