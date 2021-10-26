program spectrumsdt
  use basis_mod
  use cap_mod, only: calc_real_cap
  use config_mod, only: read_config_dict
  use debug_tools_mod
  use dict_utils_mod, only: extract_string, item_or_default
  use dictionary
  use eigensolve_mod, only: calculate_states
  use input_params_mod, only: input_params
  use iso_fortran_env, only: real64
  use mpi
  use grid_info_mod
  use grids_mod, only: generate_grids, load_grids
  use overlaps_mod
  use overlaps_extra_mod
  use parallel_utils_mod, only: print_parallel
  use potential_mod, only: load_potential
  use properties_mod
  use spectrumsdt_utils_mod, only: get_reduced_mass

  use wf_print_mod
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

    ! Debug
    real(real64), allocatable :: wf_grid(:)
    integer :: rho_ind, theta_ind, solution_ind

    call print_parallel('Stage: ' // params % stage)
    select case (params % stage)
      case ('grids')
        call generate_grids(params)

      case ('basis')
        potential = load_potential(params, rho_info % points, theta_info % points, phi_info % points)

        ! Debug
        if (debug_mode == 'print_wf' .or. &
            debug_mode == 'half_arg' .and. allocated(debug_ints) .and. debug_ints(1) == 1 .or. &
            debug_mode == '3m' .and. allocated(debug_ints) .and. debug_ints(1) == 1) then

          ! Well
          rho_ind = 13
          theta_ind = 42
          wf_grid = linspace(0d0, 2*pi, 1440)

          ! Channel
          ! rho_ind = 93 ! 46
          ! theta_ind = 64
          ! wf_grid = linspace(-pi, pi, 1440)

          solution_ind = 5
          call write_array(wf_grid, 'wf_grid')

          if (params % use_geometric_phase == 1) then
            ! ! 1D
            ! call reprint_wf_grid_1d_complex(params, rho_ind, theta_ind, wf_grid)
            ! call write_array(potential(:, theta_ind, rho_ind) * au_to_wn, 'pot')

            ! 2D
            call reprint_wf_grid_2d_complex(params, rho_ind, solution_ind, wf_grid)
            call write_matrix(potential(:, :, rho_ind) * au_to_wn, 'pot')
          else
            ! 1D
            call reprint_wf_grid_1d_real(params, rho_ind, theta_ind, wf_grid)
            call write_array(potential(:, theta_ind, rho_ind) * au_to_wn, 'pot')

            ! 2D
            ! call reprint_wf_grid_2d_real(params, rho_ind, solution_ind, wf_grid)
          end if

          stop 'print_wf'
        end if

        if (params % use_geometric_phase == 0) then
          call calculate_basis_real(params, rho_info % points, theta_info % points, phi_info % points, potential)
        else
          call calculate_basis_complex(params, rho_info % points, theta_info % points, phi_info % points, potential)
        end if

      case ('overlaps')
        if (params % use_geometric_phase == 0) then
          call calculate_overlaps_real(params)
          call calculate_overlaps_extra_real(params, get_reduced_mass(params % mass), rho_info % points, theta_info % points)
        else
          call calculate_overlaps_complex(params)
          call calculate_overlaps_extra_complex(params, get_reduced_mass(params % mass), rho_info % points, theta_info % points)
        end if

      case ('eigensolve')
        call calculate_states(params, rho_info)

      case ('properties')
        cap = calc_real_cap(params, rho_info) ! Returns 0s if no cap
        if (params % use_geometric_phase == 0) then
          call calculate_state_properties_real(params, rho_info, theta_info, cap)
        else
          call calculate_state_properties_complex(params, rho_info, theta_info, cap)
        end if
    end select
    call print_parallel('Done')
  end subroutine

end program
