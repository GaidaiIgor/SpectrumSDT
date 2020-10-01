!-------------------------------------------------------------------------------------------------------------------------------------------
! Handles parameter properties (mandatory, optional, etc.) for each stage
!-------------------------------------------------------------------------------------------------------------------------------------------
module stage_metadata_mod
  use dict_utils
  use dictionary
  implicit none

  private
  public :: get_mandatory_keys, get_optional_keys, get_all_keys

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns mandatory keys at basis stage
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_basis(config_dict) result(res)
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: res

    call add_if_absent(res, 'molecule', 'set')
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

    call add_if_absent(res, 'molecule', 'set')
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
    character(:), allocatable :: rovib_coupling, fix_basis_jk

    call add_if_absent(res, 'molecule', 'set')
    call add_if_absent(res, 'J', 'set')
    call add_if_absent(res, 'symmetry', 'set')
    call add_if_absent(res, 'num_states', 'set')
    call add_if_absent(res, 'grid_path', 'set')
    call add_if_absent(res, 'root_path', 'set')

    rovib_coupling = item_or_default(res, 'rovib_coupling', 'unset')
    fix_basis_jk = item_or_default(res, 'fix_basis_jk', 'unset')
    if (rovib_coupling == '1') then
      call add_if_absent(res, 'parity', 'set')
      if (fix_basis_jk == '1') then
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
    character(:), allocatable :: rovib_coupling, fix_basis_jk

    call add_if_absent(res, 'molecule', 'set')
    call add_if_absent(res, 'J', 'set')
    call add_if_absent(res, 'symmetry', 'set')
    call add_if_absent(res, 'basis_size_phi', 'set')
    call add_if_absent(res, 'num_states', 'set')
    call add_if_absent(res, 'grid_path', 'set')
    call add_if_absent(res, 'root_path', 'set')
    call add_if_absent(res, 'channels_root', 'set')

    rovib_coupling = item_or_default(res, 'rovib_coupling', 'unset')
    fix_basis_jk = item_or_default(res, 'fix_basis_jk', 'unset')
    if (rovib_coupling == '1') then
      call add_if_absent(res, 'parity', 'set')
      if (fix_basis_jk == '1') then
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
    if (stage == 'basis') then
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
    character(:), allocatable :: rovib_coupling

    call add_if_absent(res, 'ncv', 'set')
    call add_if_absent(res, 'mpd', 'set')
    call add_if_absent(res, 'max_iterations', 'set')
    call add_if_absent(res, 'cap_type', 'set')
    call add_if_absent(res, 'sequential', 'set')

    rovib_coupling = item_or_default(res, 'rovib_coupling', 'unset')
    if (rovib_coupling == '1') then
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
    character(:), allocatable :: rovib_coupling

    call add_if_absent(res, 'cap_type', 'set')
    call add_if_absent(res, 'sequential', 'set')

    rovib_coupling = item_or_default(res, 'rovib_coupling', 'unset')
    if (rovib_coupling == '1') then
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
    if (stage == 'basis') then
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
    call add_if_absent(res, 'rovib_coupling', 'set')
    call add_if_absent(res, 'fix_basis_jk', 'set')

    call add_if_absent(res, 'molecule', 'set')
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
    
    call add_if_absent(res, 'cap_type', 'set')
    
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
