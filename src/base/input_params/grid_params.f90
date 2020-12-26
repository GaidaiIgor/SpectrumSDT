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
    procedure :: get_optional_keys => get_optional_keys_grid_params
    procedure :: get_all_keys => get_all_keys_grid_params
    procedure :: check_values => check_values_grid_params
    procedure :: checked_init => checked_init_grid_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of grid_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_grid_params(this, config_dict)
    class(grid_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    integer :: i
    character(:), allocatable :: next_key, next_value
    type(string), allocatable :: key_set(:)

    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      next_value = extract_string(config_dict, next_key)
      select case (next_key)
        case ('prefix')
          this % prefix = next_value
        case ('from')
          this % from = str2real(next_value)
        case ('to')
          this % to = str2real(next_value)
        case ('step')
          this % step = str2real(next_value)
        case ('num_points')
          this % num_points = str2int(next_value)
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_grid_params(this) result(keys)
    class(grid_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'prefix')
    call put_string(keys, 'from')
    call put_string(keys, 'to')
    if (this % num_points == -1) then
      call put_string(keys, 'step')
    else
      call put_string(keys, 'num_points')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of optional keys. A placeholder method for calls from optgrid_params and other derived classes.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys_grid_params(this) result(keys)
    class(grid_params), intent(in) :: this
    type(dictionary_t) :: keys
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_grid_params(this) result(keys)
    class(grid_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'prefix')
    call put_string(keys, 'from')
    call put_string(keys, 'to')
    call put_string(keys, 'step')
    call put_string(keys, 'num_points')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_grid_params(this)
    class(grid_params), intent(in) :: this
    call assert(this % to > 0, 'Error: ' // this % prefix // 'to should be > 0')
    call assert(this % from < this % to, 'Error: ' // this % prefix // 'from should be < to')
    call assert(this % num_points == -1 .or. this % num_points > 0, 'Error: ' // this % prefix // 'npoints should be > 0')
    call assert((this % step .aeq. -1d0) .or. this % step > 0, 'Error: '// this % prefix // 'step should be > 0')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of grid_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_grid_params(this, config_dict, check_extra)
    class(grid_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    integer, optional, intent(in) :: check_extra
    integer :: check_extra_act
    type(dictionary_t) :: mandatory_keys, all_keys

    call this % assign_dict(config_dict)
    ! These two can never be set together, so it's an error if they are
    ! Has to be checked before mandatory keys since error message is more specific
    call check_only_one_set(config_dict, string([character(100) :: 'step', 'num_points']), this % prefix)
    mandatory_keys = this % get_mandatory_keys()
    call check_mandatory_keys(config_dict, mandatory_keys, this % prefix)
    call this % check_values()

    ! Skippable to enable calls from derived types since their parameters will be erroneously considered unknown here
    check_extra_act = arg_or_default(check_extra, 1)
    if (check_extra_act == 1) then
      all_keys = this % get_all_keys()
      call check_extra_keys(config_dict, all_keys, this % prefix)
    end if
  end subroutine

end module
