program spectrumsdt
  use cap_mod, only: init_caps, get_real_cap
  use config_mod, only: read_config_dict
  use debug_tools
  use dict_utils, only: extract_string, item_or_default
  use dictionary
  use formulas_mod, only: get_reduced_mass
  use general_vars, only: g1, g2, g3, init_grids
  use input_params_mod, only: input_params
  use iso_fortran_env, only: real64
  use mpi
  use grids_mod, only: generate_grids
  use overlaps_extra_mod, only: calculate_overlaps_extra
  use parallel_utils, only: print_parallel
  use potential_mod, only: init_potential
  use sdt, only: calculate_basis, calculate_overlaps
  use spectrum_mod, only: calculate_states
  use state_properties_mod, only: calculate_state_properties
  implicit none

  integer :: ierr
  type(dictionary_t) :: config_dict, auxiliary_info
  type(input_params) :: params

  call read_config_dict('spectrumsdt.config', config_dict, auxiliary_info)
  ! We have to decide whether to enable MPI before we create an instance of input_params, otherwise we cannot restrict printing during instantiation to 0th processor
  if (extract_string(config_dict, 'stage') /= 'grids' .and. item_or_default(config_dict, 'use_parallel', '1') == '1') then
    call MPI_Init(ierr)
  end if
  call params % checked_init(config_dict, auxiliary_info)

  call init_all(params)
  call process_stage(params)

  if (params % use_parallel == 1) then
    call MPI_Finalize(ierr)
  end if

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Inits all necessary information.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_all(params)
    class(input_params), intent(inout) :: params
    call init_debug(params)
    if (params % stage == 'basis' .or. params % stage == 'overlaps' .or. params % stage == 'eigencalc' .or. params % stage == 'properties') then
      call init_grids(params)
      call params % check_resolve_grids(g1, g2, g3)
    end if
    if (params % stage == 'basis') then
      call init_potential(params)
    end if
    if (params % stage == 'eigencalc' .or. params % stage == 'properties') then
      call init_caps(params, g1)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Selects appropriate procedures for the current stage.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine process_stage(params)
    type(input_params), intent(in) :: params
    real(real64) :: mu

    call print_parallel('Stage: ' // params % stage)
    mu = get_reduced_mass(params % mass)
    select case (params % stage)
      case ('grids')
        call generate_grids(params)
      case ('basis')
        call calculate_basis(params)
      case ('overlaps')
        call calculate_overlaps(params)
        call calculate_overlaps_extra(params, mu, g1, g2)
      case ('eigencalc')
        call calculate_states(params)
      case ('properties')
        if (params % cap % type /= 'none') then
          call calculate_state_properties(params, g1, g2, get_real_cap())
        else
          call calculate_state_properties(params, g1, g2)
        end if
    end select
    call print_parallel('Done')
  end subroutine

end program
