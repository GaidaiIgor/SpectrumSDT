!-------------------------------------------------------------------------------------------------------------------------------------------
! Module for printing wave functions on grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
module wf_print_mod
  use basis_mod
  use general_utils_mod
  use input_params_mod
  use overlaps_mod
  use spectrumsdt_io_mod
  use spectrumsdt_paths_mod
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints specified 2D wave function on grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_wf_grid_2d(params, rho_ind, wf_ind, num_points_theta)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind, wf_ind, num_points_theta
    character(:), allocatable :: sym_path, out_file_name

    sym_path = get_sym_path(params)
    call load_solutions_1D(get_solutions_1d_path(sym_path, rho_ind), num_points_theta, params % basis % num_funcs_phi_per_sym, num_solutions_1d_col, energies_1d_col, solutions_1d_col)
    call load_solutions_2D(get_solutions_2d_path(sym_path, rho_ind), energies_2d_col, solutions_2d_col)

    out_file_name = 'wf_2d.' // num2str(rho_ind) // '.' // num2str(wf_ind)
  end subroutine

end module

