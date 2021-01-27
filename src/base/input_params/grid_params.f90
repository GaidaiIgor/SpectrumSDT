module grid_params_mod
  use config_mod
  use dictionary
  use dict_utils
  use general_utils
  use iso_fortran_env, only: real64
  use string_mod
  implicit none

  private
  public :: grid_params

  type :: grid_params
    character(:), allocatable :: prefix
    real(real64) :: from = -1
    real(real64) :: to = -1
    real(real64) :: step = -1
    integer :: num_points = -1

  contains
    procedure :: assign_dict => assign_dict_grid_params
    procedure :: get_mandatory_keys => get_mandatory_keys_grid_params
    procedure :: get_all_keys => get_all_keys_grid_params
    procedure :: set_defaults => set_defaults_grid_params
    procedure :: check_values => check_values_grid_params
    procedure :: checked_init => checked_init_grid_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of grid_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_grid_params(this, config_dict, auxiliary_info)
    class(grid_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    integer :: i
    character(:), allocatable :: next_key, full_key, next_value
    type(string), allocatable :: key_set(:)

    this % prefix = extract_string(auxiliary_info, 'prefix')
    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      full_key = this % prefix // next_key
      next_value = extract_string(config_dict, next_key)
      select case (next_key)
        case ('from')
          this % from = str2real_config(next_value, full_key)
        case ('to')
          this % to = str2real_config(next_value, full_key)
        case ('step')
          this % step = str2real_config(next_value, full_key)
        case ('num_points')
          this % num_points = str2int_config(next_value, full_key)
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_grid_params(this) result(keys)
    class(grid_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'from')
    call put_string(keys, 'to')
    if (this % num_points == -1) then
      call put_string(keys, 'step')
    else
      call put_string(keys, 'num_points')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_grid_params(this) result(keys)
    class(grid_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'from')
    call put_string(keys, 'to')
    call put_string(keys, 'step')
    call put_string(keys, 'num_points')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets defaults.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_defaults_grid_params(this, config_dict)
    class(grid_params), intent(inout) :: this
    class(dictionary_t), intent(in) :: config_dict
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_grid_params(this)
    class(grid_params), intent(in) :: this
    call assert(this % from .age. 0d0, 'Error: ' // this % prefix // 'from should be >= 0')
    call assert(this % to > 0, 'Error: ' // this % prefix // 'to should be > 0')
    call assert(this % from < this % to, 'Error: ' // this % prefix // 'from should be < to')
    call assert(this % num_points == -1 .or. this % num_points > 0, 'Error: ' // this % prefix // 'npoints should be > 0')
    call assert((this % step .aeq. -1d0) .or. this % step > 0, 'Error: '// this % prefix // 'step should be > 0')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of grid_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_grid_params(this, config_dict, auxiliary_info, check_extra, check_values)
    class(grid_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    integer, optional, intent(in) :: check_extra, check_values
    integer :: check_extra_act, check_values_act
    type(dictionary_t) :: mandatory_keys, all_keys

    ! Skippable to enable calls from derived types since their parameters will be unknown here
    check_extra_act = arg_or_default(check_extra, 1)
    if (check_extra_act == 1) then
      all_keys = this % get_all_keys()
      call check_extra_keys(config_dict, all_keys, this % prefix)
    end if

    call check_key_types(config_dict, auxiliary_info, 'string')
    call this % assign_dict(config_dict, auxiliary_info)
    ! These two can never be set together, so it's an error if they are
    ! Has to be checked before mandatory keys since error message is more specific
    call check_only_one_set(config_dict, string([character(100) :: 'step', 'num_points']), this % prefix)
    mandatory_keys = this % get_mandatory_keys()
    call check_mandatory_keys(config_dict, mandatory_keys, this % prefix)
    call this % set_defaults(config_dict)

    check_values_act = arg_or_default(check_values, 1)
    if (check_values_act == 1) then
      call this % check_values()
    end if
  end subroutine

end module
