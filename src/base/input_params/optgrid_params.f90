module optgrid_params_mod
  use config_mod
  use dictionary
  use dict_utils
  use grid_params_mod
  use iso_fortran_env, only: real64
  use string_mod
  implicit none

  type, extends(grid_params) :: optgrid_params
    integer :: optimized = -1 ! optimized distribution of grid points
    character(:), allocatable :: envelope_path
    real(real64) :: max_energy = -1
    integer :: solver_steps = -1

  contains
    procedure :: assign_dict => assign_dict_optgrid_params
    procedure :: get_mandatory_keys => get_mandatory_keys_optgrid_params
    procedure :: get_optional_keys => get_optional_keys_optgrid_params
    procedure :: get_all_keys => get_all_keys_optgrid_params
    procedure :: check_values => check_values_optgrid_params
    procedure :: checked_init => checked_init_optgrid_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of grid_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_optgrid_params(this, config_dict)
    class(optgrid_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    integer :: i
    character(:), allocatable :: next_key, next_value
    type(string), allocatable :: key_set(:)

    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      next_value = extract_string(config_dict, next_key)
      select case (next_key)
        case ('optimized')
          this % optimized = str2int(next_value)
        case ('envelope_path')
          this % envelope_path = next_value
        case ('max_energy')
          this % max_energy = str2real(next_value)
        case ('solver_steps')
          this % solver_steps = str2int(next_value)
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_optgrid_params(this) result(keys)
    class(optgrid_params), intent(in) :: this
    type(dictionary_t) :: keys

    keys = this % grid_params % get_mandatory_keys()
    if (this % optimized == 1) then
      call put_string(keys, 'envelope_path')
      call put_string(keys, 'max_energy')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of optional keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys_optgrid_params(this) result(keys)
    class(optgrid_params), intent(in) :: this
    type(dictionary_t) :: keys

    keys = this % grid_params % get_optional_keys()
    call put_string(keys, 'optimized', '0')
    if (this % optimized == 1) then
      call put_string(keys, 'solver_steps', '2048')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_optgrid_params(this) result(keys)
    class(optgrid_params), intent(in) :: this
    type(dictionary_t) :: keys

    keys = this % grid_params % get_all_keys()
    call put_string(keys, 'optimized')
    call put_string(keys, 'envelope_path')
    call put_string(keys, 'max_energy')
    call put_string(keys, 'solver_steps')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_optgrid_params(this) 
    class(optgrid_params), intent(in) :: this
    call assert(any(this % optimized == [-1, 0, 1]), 'Error: optimized should be 0 or 1')
    call assert(this % solver_steps == -1 .or. this % solver_steps > 0, 'Error: solver_steps should be > 0')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of grid_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_optgrid_params(this, config_dict, check_extra)
    class(optgrid_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    integer, optional, intent(in) :: check_extra
    integer :: check_extra_act
    type(dictionary_t) :: mandatory_keys, optional_keys, optional_nonset_keys, all_keys

    call this % grid_params % checked_init(config_dict, check_extra = 0) ! First, check init parental type
    call this % assign_dict(config_dict) ! Init type as is to ease branching in getting keys
    mandatory_keys = this % get_mandatory_keys() ! Save mandatory keys because we will reuse them later for unused key check
    call check_mandatory_keys(config_dict, mandatory_keys) ! Make sure mandatory keys are set to ease checking values

    check_extra_act = arg_or_default(check_extra, 1)
    if (check_extra_act == 1) then
      all_keys = this % get_all_keys() ! Generally optional part comes last
      call check_extra_keys(config_dict, all_keys) ! Checks that unknown keys were not specified
      optional_keys = this % get_optional_keys()
      call check_unused_keys(config_dict, mandatory_keys, optional_keys) ! Prints info if one of the known keys will not be used with current settings
    end if

    if (len(optional_keys) > 0) then
      optional_nonset_keys = set_difference(optional_keys, config_dict)
      call this % assign_dict(optional_nonset_keys)
      call announce_defaults(optional_nonset_keys)
    end if
    call this % check_values() ! Make sure the settings have valid values
  end subroutine

end module
