module eigencalc_params_mod
  use config_mod
  use constants
  use dictionary
  use dict_utils
  use general_utils
  use iso_fortran_env, only: real64
  use parallel_utils
  use string_mod
  implicit none

  private
  public :: eigencalc_params

  type :: eigencalc_params
    character(:), allocatable :: prefix
    integer :: num_states = -1
    integer :: ncv = -1
    integer :: mpd = -1
    integer :: max_iterations = 10000

  contains
    procedure :: assign_dict => assign_dict_eigencalc_params
    procedure :: get_mandatory_keys => get_mandatory_keys_eigencalc_params
    procedure :: get_all_keys => get_all_keys_eigencalc_params
    procedure :: set_defaults => set_defaults_eigencalc_params
    procedure :: check_values => check_values_eigencalc_params
    procedure :: checked_init => checked_init_eigencalc_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of eigencalc_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_eigencalc_params(this, config_dict, auxiliary_info)
    class(eigencalc_params), intent(inout) :: this
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
        case ('num_states')
          this % num_states = str2int_config(next_value, full_key)
        case ('ncv')
          this % ncv = str2int_config(next_value, full_key)
        case ('mpd')
          this % mpd = str2int_config(next_value, full_key)
        case ('max_iterations')
          this % max_iterations = str2int_config(next_value, full_key)
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_eigencalc_params(this) result(keys)
    class(eigencalc_params), intent(in) :: this
    type(dictionary_t) :: keys
    call put_string(keys, 'num_states')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_eigencalc_params(this) result(keys)
    class(eigencalc_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'num_states')
    call put_string(keys, 'ncv')
    call put_string(keys, 'mpd')
    call put_string(keys, 'max_iterations')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets defaults.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_defaults_eigencalc_params(this, config_dict, auxiliary_info) 
    class(eigencalc_params), intent(in) :: this
    class(dictionary_t), intent(in) :: config_dict, auxiliary_info
    character(:), allocatable :: prefix

    prefix = extract_string(auxiliary_info, 'prefix')
    if (.not. ('mpd' .in. config_dict)) then
      call print_parallel(prefix // 'mpd is not specified. Its value will be determined by SLEPc')
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_eigencalc_params(this) 
    class(eigencalc_params), intent(in) :: this
    call assert(this % num_states > 0, 'Error: ' // this % prefix // 'num_states should be >= 0')
    call assert(this % ncv == -1 .or. this % ncv > this % num_states, 'Error: ' // this % prefix // 'ncv should be > ' // this % prefix // ' num_states')
    call assert(this % mpd == -1 .or. this % mpd > 1, 'Error: ' // this % prefix // 'mpd should be > 1')
    call assert(this % max_iterations > 0, 'Error: ' // this % prefix // 'max_iterations should be > 0')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of eigencalc_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_eigencalc_params(this, config_dict, auxiliary_info)
    class(eigencalc_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    type(dictionary_t) :: mandatory_keys, all_keys

    call check_key_types(config_dict, auxiliary_info, 'string')
    call this % assign_dict(config_dict, auxiliary_info)
    mandatory_keys = this % get_mandatory_keys()
    call check_mandatory_keys(config_dict, mandatory_keys, this % prefix)

    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys, this % prefix)
    call this % set_defaults(config_dict, auxiliary_info)
    call this % check_values()
  end subroutine

end module
