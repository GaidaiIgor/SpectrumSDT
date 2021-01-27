program spectrumsdt
  use cap_mod, only: calc_real_cap
  use config_mod, only: read_config_dict
  use debug_tools
  use dict_utils, only: extract_string, item_or_default
  use dictionary
  use formulas_mod, only: get_reduced_mass
  use input_params_mod, only: input_params
  use iso_fortran_env, only: real64
  use mpi
  use grid_info_mod
  use grids_mod, only: generate_grids, load_grids
  use overlaps_extra_mod, only: calculate_overlaps_extra
  use parallel_utils, only: print_parallel
  use potential_mod, only: load_potential
  use sdt, only: calculate_basis, calculate_overlaps
  use spectrum_mod, only: calculate_states
  use state_properties_mod, only: calculate_state_properties
  implicit none

  integer :: ierr
  type(dictionary_t) :: config_dict, auxiliary_info
  type(input_params) :: params
  type(grid_info) :: rho_info, theta_info, phi_info

  call read_config_dict('spectrumsdt.config', config_dict, auxiliary_info)
  ! We have to decide whether to enable MPI before we create an instance of input_params, otherwise we cannot restrict printing during instantiation to 0th processor
  if (extract_string(config_dict, 'stage') /= 'grids' .and. item_or_default(config_dict, 'use_parallel', '1') == '1') then
    call MPI_Init(ierr)
  end if
  call params % checked_init(config_dict, auxiliary_info)
  call init_debug(params)

  if (params % stage /= 'grids') then
    call load_grids(params, rho_info, theta_info, phi_info)
    call params % check_resolve_grids(rho_info, theta_info, phi_info)
  end if
  call process_stage(params, rho_info, theta_info, phi_info)

  if (params % use_parallel == 1) then
    call MPI_Finalize(ierr)
  end if

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Selects appropriate procedures for the current stage.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine process_stage(params, rho_info, theta_info, phi_info)
    type(input_params), intent(in) :: params
    class(grid_info), intent(in) :: rho_info, theta_info, phi_info
    real(real64), allocatable :: cap(:)
    real(real64), allocatable :: potential(:, :, :)

    call print_parallel('Stage: ' // params % stage)
    select case (params % stage)
      case ('grids')
        call generate_grids(params)
      case ('basis')
        potential = load_potential(params, rho_info % points, theta_info % points, phi_info % points)
        call calculate_basis(params, rho_info % points, theta_info % points, phi_info % points, potential)
      case ('overlaps')
        call calculate_overlaps(params, size(theta_info % points))
        call calculate_overlaps_extra(params, get_reduced_mass(params % mass), rho_info % points, theta_info % points)
      case ('eigensolve')
        call calculate_states(params, rho_info)
      case ('properties')
        if (params % cap % type /= 'none') then
          cap = calc_real_cap(params, rho_info)
          call calculate_state_properties(params, rho_info, theta_info, cap)
        else
          call calculate_state_properties(params, rho_info, theta_info)
        end if
    end select
    call print_parallel('Done')
  end subroutine

end program
