module debug_params_mod
  use config_mod
  use dictionary
  use dict_utils
  use general_utils
  use iso_fortran_env, only: real64
  use string_mod
  implicit none

  private
  public :: debug_params

  type :: debug_params
    character(:), allocatable :: prefix
    integer :: enable_terms(2) = 1 ! 1st digit - coriolis, 2nd - asymmetric
    integer :: optimized_mult = 1 ! disables matrix-vector multiplication optimizations
    integer :: treat_tp_as_xy = 0 ! treats theta and phi grids as x and y grids for aph plots
    character(:), allocatable :: mode
    integer :: int_param = 0
    real(real64) :: real_param = 0d0

  contains
    procedure :: assign_dict => assign_dict_debug_params
    procedure :: get_all_keys => get_all_keys_debug_params
    procedure :: check_values => check_values_debug_params
    procedure :: checked_init => checked_init_debug_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses user-provided value of `enable_terms`.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_enable_terms(enable_terms_str) result(enable_terms)
    character(*), intent(in) :: enable_terms_str
    integer :: enable_terms(2)
    
    enable_terms(1) = str2int(enable_terms_str(1:1))
    enable_terms(2) = str2int(enable_terms_str(2:2))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of debug_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_debug_params(this, config_dict)
    class(debug_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    integer :: i
    character(:), allocatable :: next_key, next_value
    type(string), allocatable :: key_set(:)

    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      next_value = extract_string(config_dict, next_key)

      select case (next_key)
        case ('enable_terms')
          this % enable_terms = parse_enable_terms(next_value)
        case ('optimized_mult')
          this % optimized_mult = str2int(next_value)
        case ('treat_tp_as_xy')
          this % treat_tp_as_xy = str2int(next_value)
        case ('mode')
          this % mode = next_value
        case ('int_param')
          this % int_param = str2int(next_value)
        case ('real_param')
          this % real_param = str2real(next_value)
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_debug_params(this) result(keys)
    class(debug_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'enable_terms')
    call put_string(keys, 'optimized_mult')
    call put_string(keys, 'treat_tp_as_xy')
    call put_string(keys, 'mode')
    call put_string(keys, 'int_param')
    call put_string(keys, 'real_param')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_debug_params(this) 
    class(debug_params), intent(in) :: this
    call assert(any(this % enable_terms(1) == [0, 1]), 'Error: ' // this % prefix // 'enable_terms(1) should be 0 or 1')
    call assert(any(this % enable_terms(2) == [0, 1]), 'Error: ' // this % prefix // 'enable_terms(2) should be 0 or 1')
    call assert(any(this % optimized_mult == [0, 1]), 'Error: ' // this % prefix // 'optimized_mult should be 0 or 1')
    call assert(any(this % treat_tp_as_xy == [0, 1]), 'Error: ' // this % prefix // 'treat_tp_as_xy should be 0 or 1')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of debug_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_debug_params(this, config_dict, auxiliary_info)
    class(debug_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    type(dictionary_t) :: mandatory_keys, all_keys

    this % prefix = extract_string(auxiliary_info, 'prefix')
    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys, this % prefix)
    call check_key_types(config_dict, auxiliary_info, 'string')
    call this % assign_dict(config_dict)
    call this % check_values()
  end subroutine

end module
