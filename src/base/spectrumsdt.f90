program spectrumsdt
  use cap_mod, only: init_caps, get_real_cap
  use config_mod, only: read_config_dict
  use debug_tools
  use dict_utils, only: extract_string, item_or_default
  use dictionary
  use general_vars, only: mu, g1, g2, init_masses, init_grids
  use input_params_mod, only: input_params
  use iso_fortran_env, only: real64
  use mpi
  use grids_mod, only: generate_grids
  use overlaps_extra_mod, only: calculate_overlaps_extra
  use parallel_utils, only: print_parallel
  use potential_mod, only: init_potential
  use sdt, only: init_sdt, init_output_directories, calc_basis, calc_overlap
  use spectrum_mod, only: calculate_states
  use state_properties_mod, only: calculate_state_properties
  implicit none

  integer :: ierr
  type(dictionary_t) :: config_dict
  type(input_params) :: params

  config_dict = read_config_dict('spectrumsdt.config')
  ! We have to decide whether to enable MPI before we create an instance of input_params, otherwise we cannot restrict printing during instantiation to 0th processor
  if (extract_string(config_dict, 'stage') /= 'grids' .and. item_or_default(config_dict, 'use_parallel', '1') == '1') then
    call MPI_Init(ierr)
  end if
  call params % checked_init(config_dict)

  call init_debug(params)
  call init_masses(params)
  if (params % stage /= 'grids') then
    call init_grids(params)
  end if
  if (params % stage == 'basis') then
    call init_potential(params)
  end if
  if (params % stage /= 'grids') then
    call init_sdt(params)
  end if
  if (params % stage == 'eigencalc' .or. params % stage == 'properties') then
    call init_caps(params)
  end if
  call init_output_directories(params)
  call process_stage(params)

  if (params % use_parallel == 1) then
    call MPI_Finalize(ierr)
  end if

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Selects appropriate procedures for the current stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine process_stage(params)
    type(input_params), intent(in) :: params

    ! Write information
    call print_parallel('Stage: ' // params % stage)
    ! Call appropriate subroutine
    select case (params % stage)
      case ('grids')
        call generate_grids(params)

      case ('basis')
        call calc_basis(params)

      case ('overlaps')
        call calc_overlap(params)
        call calculate_overlaps_extra(params, mu, g1, g2)

      case ('eigencalc')
        call calculate_states(params)

      case ('properties')
        if (params % cap_type /= 'none') then
          call calculate_state_properties(params, g1, size(g2), get_real_cap())
        else
          call calculate_state_properties(params, g1, size(g2))
        end if
    end select
    call print_parallel('Done')
  end subroutine

end program
