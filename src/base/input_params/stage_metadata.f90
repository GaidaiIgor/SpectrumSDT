!-------------------------------------------------------------------------------------------------------------------------------------------
! Handles parameter properties (mandatory, optional, etc.) for each stage
!-------------------------------------------------------------------------------------------------------------------------------------------
module stage_metadata_mod
  use dict_utils
  use dictionary
  use string_mod
  implicit none

  private
  public :: get_mandatory_keys, get_optional_keys, get_all_keys

contains
!-------------------------------------------------------------------------------------------------------------------------------------------
! Adds the first key from *keys* present in *config_dict* to *key_set*
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine add_first_present(key_set, config_dict, keys)
    class(dictionary_t), intent(inout) :: key_set
    class(dictionary_t) :: config_dict ! intent(in)
    class(string), intent(in) :: keys(:)
    integer :: i
    character(:), allocatable :: next_key

    do i = 1, size(keys)
      next_key = keys(i) % to_char_str()
      if (next_key .in. config_dict) then
        call add_if_absent(key_set, next_key, 'set')
        exit
      end if
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns mandatory keys at grids stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_grids(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    character(:), allocatable :: use_optimized_grid_rho

    call add_if_absent(res, 'stage', 'set')
    call add_if_absent(res, 'grid_rho_from', 'set')
    call add_if_absent(res, 'grid_rho_to', 'set')
    call add_if_absent(res, 'grid_theta_from', 'set')
    call add_if_absent(res, 'grid_theta_to', 'set')
    call add_if_absent(res, 'grid_phi_from', 'set')
    call add_if_absent(res, 'grid_phi_to', 'set')

    call add_first_present(res, config_dict, string([character(100) :: 'grid_rho_npoints', 'grid_rho_step']))
    call add_first_present(res, config_dict, string([character(100) :: 'grid_theta_npoints', 'grid_theta_step']))
    call add_first_present(res, config_dict, string([character(100) :: 'grid_phi_npoints', 'grid_phi_step']))

    use_optimized_grid_rho = item_or_default(config_dict, 'use_optimized_grid_rho', '0')
    if (use_optimized_grid_rho == '1') then
      call add_if_absent(res, 'mass_central', 'set')
      call add_if_absent(res, 'mass_terminal1', 'set')
      call add_if_absent(res, 'mass_terminal2', 'set')
      call add_if_absent(res, 'envelope_rho_path', 'set')
      call add_if_absent(res, 'envelope_rho_max_energy', 'set')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns mandatory keys at basis stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_basis(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res

    call add_if_absent(res, 'stage', 'set')
    call add_if_absent(res, 'mass_central', 'set')
    call add_if_absent(res, 'mass_terminal1', 'set')
    call add_if_absent(res, 'mass_terminal2', 'set')
    call add_if_absent(res, 'J', 'set')
    call add_if_absent(res, 'K', 'set')
    call add_if_absent(res, 'symmetry', 'set')
    call add_if_absent(res, 'basis_size_phi', 'set')
    call add_if_absent(res, 'cutoff_energy', 'set')
    call add_if_absent(res, 'grid_path', 'set')
    call add_if_absent(res, 'root_path', 'set')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns mandatory keys at overlaps stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_overlaps(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res

    call add_if_absent(res, 'stage', 'set')
    call add_if_absent(res, 'use_rovib_coupling', 'set')
    call add_if_absent(res, 'use_fix_basis_jk', 'set')
    call add_if_absent(res, 'mass_central', 'set')
    call add_if_absent(res, 'mass_terminal1', 'set')
    call add_if_absent(res, 'mass_terminal2', 'set')
    call add_if_absent(res, 'K', 'set')
    call add_if_absent(res, 'symmetry', 'set')
    call add_if_absent(res, 'basis_size_phi', 'set')
    call add_if_absent(res, 'cutoff_energy', 'set')
    call add_if_absent(res, 'grid_path', 'set')
    call add_if_absent(res, 'root_path', 'set')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns mandatory keys at eigencalc stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_eigencalc(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    character(:), allocatable :: use_rovib_coupling, use_fix_basis_jk

    call add_if_absent(res, 'stage', 'set')
    call add_if_absent(res, 'use_rovib_coupling', 'set')
    call add_if_absent(res, 'use_fix_basis_jk', 'set')
    call add_if_absent(res, 'mass_central', 'set')
    call add_if_absent(res, 'mass_terminal1', 'set')
    call add_if_absent(res, 'mass_terminal2', 'set')
    call add_if_absent(res, 'J', 'set')
    call add_if_absent(res, 'symmetry', 'set')
    call add_if_absent(res, 'num_states', 'set')
    call add_if_absent(res, 'grid_path', 'set')
    call add_if_absent(res, 'root_path', 'set')

    use_rovib_coupling = item_or_default(res, 'use_rovib_coupling', 'unset')
    use_fix_basis_jk = item_or_default(res, 'use_fix_basis_jk', 'unset')
    if (use_rovib_coupling == '1') then
      call add_if_absent(res, 'parity', 'set')
      if (use_fix_basis_jk == '1') then
        call add_if_absent(res, 'basis_root_path', 'set')
        call add_if_absent(res, 'basis_J', 'set')
        call add_if_absent(res, 'basis_K', 'set')
      end if
    else
      call add_if_absent(res, 'K', 'set')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns mandatory keys at properties stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_properties(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    character(:), allocatable :: use_rovib_coupling, use_fix_basis_jk

    call add_if_absent(res, 'stage', 'set')
    call add_if_absent(res, 'use_rovib_coupling', 'set')
    call add_if_absent(res, 'use_fix_basis_jk', 'set')
    call add_if_absent(res, 'mass_central', 'set')
    call add_if_absent(res, 'mass_terminal1', 'set')
    call add_if_absent(res, 'mass_terminal2', 'set')
    call add_if_absent(res, 'J', 'set')
    call add_if_absent(res, 'symmetry', 'set')
    call add_if_absent(res, 'basis_size_phi', 'set')
    call add_if_absent(res, 'num_states', 'set')
    call add_if_absent(res, 'grid_path', 'set')
    call add_if_absent(res, 'root_path', 'set')
    call add_if_absent(res, 'channels_root', 'set')

    use_rovib_coupling = item_or_default(res, 'use_rovib_coupling', 'unset')
    use_fix_basis_jk = item_or_default(res, 'use_fix_basis_jk', 'unset')
    if (use_rovib_coupling == '1') then
      call add_if_absent(res, 'parity', 'set')
      if (use_fix_basis_jk == '1') then
        call add_if_absent(res, 'basis_root_path', 'set')
        call add_if_absent(res, 'basis_J', 'set')
        call add_if_absent(res, 'basis_K', 'set')
      end if
    else
      call add_if_absent(res, 'K', 'set')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns mandatory keys for this stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    character(:), allocatable :: stage

    stage = item_or_default(config_dict, 'stage', 'unset')
    if (stage == 'grids') then
      res = get_mandatory_keys_grids(config_dict)
    else if (stage == 'basis') then
      res = get_mandatory_keys_basis(config_dict)
    else if (stage == 'overlaps') then
      res = get_mandatory_keys_overlaps(config_dict)
    else if (stage == 'eigencalc') then
      res = get_mandatory_keys_eigencalc(config_dict)
    else if (stage == 'properties') then
      res = get_mandatory_keys_properties(config_dict)
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns optional keys at grids stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys_grids(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    character(:), allocatable :: use_optimized_grid_rho

    call add_if_absent(res, 'use_optimized_grid_rho', 'set')
    use_optimized_grid_rho = item_or_default(config_dict, 'use_optimized_grid_rho', '0')
    if (use_optimized_grid_rho == '1') then
      call add_if_absent(res, 'optimized_grid_solver_steps', 'set')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns optional keys at basis stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys_basis(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    call add_if_absent(res, 'sequential', 'set')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns optional keys at overlaps stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys_overlaps(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    call add_if_absent(res, 'sequential', 'set')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns optional keys at eigencalc stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys_eigencalc(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    character(:), allocatable :: use_rovib_coupling

    call add_if_absent(res, 'cap_type', 'set')
    call add_if_absent(res, 'ncv', 'set')
    call add_if_absent(res, 'mpd', 'set')
    call add_if_absent(res, 'max_iterations', 'set')
    call add_if_absent(res, 'sequential', 'set')

    use_rovib_coupling = item_or_default(res, 'use_rovib_coupling', 'unset')
    if (use_rovib_coupling == '1') then
      call add_if_absent(res, 'K', 'set')
      call add_if_absent(res, 'enable_terms', 'set')
      call add_if_absent(res, 'optimized_mult', 'set')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns optional keys at properties stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys_properties(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    character(:), allocatable :: use_rovib_coupling

    call add_if_absent(res, 'cap_type', 'set')
    call add_if_absent(res, 'sequential', 'set')

    use_rovib_coupling = item_or_default(res, 'use_rovib_coupling', 'unset')
    if (use_rovib_coupling == '1') then
      call add_if_absent(res, 'K', 'set')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns optional keys for this stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res
    character(:), allocatable :: stage

    stage = item_or_default(config_dict, 'stage', 'unset')
    if (stage == 'grids') then
      res = get_optional_keys_grids(config_dict) 
    else if (stage == 'basis') then
      res = get_optional_keys_basis(config_dict)
    else if (stage == 'overlaps') then
      res = get_optional_keys_overlaps(config_dict)
    else if (stage == 'eigencalc') then
      res = get_optional_keys_eigencalc(config_dict)
    else if (stage == 'properties') then
      res = get_optional_keys_properties(config_dict)
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns all known keys (reflection of input_params type)
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys() result(res)
    type(dictionary_t) :: res

    call add_if_absent(res, 'stage', 'set')
    call add_if_absent(res, 'use_optimized_grid_rho', 'set')
    call add_if_absent(res, 'use_rovib_coupling', 'set')
    call add_if_absent(res, 'use_fix_basis_jk', 'set')
    call add_if_absent(res, 'cap_type', 'set')

    call add_if_absent(res, 'grid_rho_from', 'set')
    call add_if_absent(res, 'grid_rho_to', 'set')
    call add_if_absent(res, 'grid_rho_npoints', 'set')
    call add_if_absent(res, 'grid_rho_step', 'set')
    call add_if_absent(res, 'envelope_rho_path', 'set')
    call add_if_absent(res, 'envelope_rho_max_energy', 'set')
    call add_if_absent(res, 'optimized_grid_solver_steps', 'set')

    call add_if_absent(res, 'grid_theta_from', 'set')
    call add_if_absent(res, 'grid_theta_to', 'set')
    call add_if_absent(res, 'grid_theta_npoints', 'set')
    call add_if_absent(res, 'grid_theta_step', 'set')

    call add_if_absent(res, 'grid_phi_from', 'set')
    call add_if_absent(res, 'grid_phi_to', 'set')
    call add_if_absent(res, 'grid_phi_npoints', 'set')
    call add_if_absent(res, 'grid_phi_step', 'set')

    call add_if_absent(res, 'mass_central', 'set')
    call add_if_absent(res, 'mass_terminal1', 'set')
    call add_if_absent(res, 'mass_terminal2', 'set')
    call add_if_absent(res, 'J', 'set')
    call add_if_absent(res, 'K', 'set')
    call add_if_absent(res, 'parity', 'set')
    call add_if_absent(res, 'symmetry', 'set')

    call add_if_absent(res, 'basis_size_phi', 'set')
    call add_if_absent(res, 'cutoff_energy', 'set')
    call add_if_absent(res, 'basis_root_path', 'set')
    call add_if_absent(res, 'basis_J', 'set')
    call add_if_absent(res, 'basis_K', 'set')

    call add_if_absent(res, 'num_states', 'set')
    call add_if_absent(res, 'ncv', 'set')
    call add_if_absent(res, 'mpd', 'set')
    call add_if_absent(res, 'max_iterations', 'set')
    
    call add_if_absent(res, 'grid_path', 'set')
    call add_if_absent(res, 'root_path', 'set')
    call add_if_absent(res, 'channels_root', 'set')

    call add_if_absent(res, 'sequential', 'set')
    call add_if_absent(res, 'enable_terms', 'set')
    call add_if_absent(res, 'optimized_mult', 'set')
    call add_if_absent(res, 'debug_mode', 'set')
    call add_if_absent(res, 'test_mode', 'set')
    call add_if_absent(res, 'debug_param_1', 'set')
  end function

end module
