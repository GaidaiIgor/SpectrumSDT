module cap_params_mod
  use config_mod
  use dictionary
  use dict_utils
  use general_utils
  use iso_fortran_env, only: real64
  use string_mod
  implicit none

  private
  public :: cap_params

  type :: cap_params
    character(:), allocatable :: prefix
    character(:), allocatable :: type
    real(real64) :: emin = -1

  contains
    procedure :: assign_dict => assign_dict_cap_params
    procedure :: get_mandatory_keys => get_mandatory_keys_cap_params
    procedure :: get_optional_keys => get_optional_keys_cap_params
    procedure :: get_all_keys => get_all_keys_cap_params
    procedure :: check_values => check_values_cap_params
    procedure :: checked_init => checked_init_cap_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of cap_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_cap_params(this, config_dict)
    class(cap_params), intent(inout) :: this
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
        case ('type')
          this % type = next_value
        case ('emin')
          this % emin = str2real(next_value)
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_cap_params(this) result(keys)
    class(cap_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'prefix')
    if (this % type == 'Manolopoulos') then
      call put_string(keys, 'emin')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of optional keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys_cap_params(this) result(keys)
    class(cap_params), intent(in) :: this
    type(dictionary_t) :: keys
    call put_string(keys, 'type', 'none')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_cap_params(this) result(keys)
    class(cap_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'prefix')
    call put_string(keys, 'type')
    call put_string(keys, 'emin')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_cap_params(this) 
    class(cap_params), intent(in) :: this
    call assert(any(this % type == [character(100) :: 'none', 'Manolopoulos']), 'Error: ' // this % prefix // 'type can be "none" or "Manolopoulos"')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of cap_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_cap_params(this, config_dict)
    class(cap_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: mandatory_keys, optional_keys, optional_nonset_keys, all_keys

    call this % assign_dict(config_dict)
    mandatory_keys = this % get_mandatory_keys()
    call check_mandatory_keys(config_dict, mandatory_keys, this % prefix)

    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys, this % prefix)
    optional_keys = this % get_optional_keys()
    call check_unused_keys(config_dict, mandatory_keys, optional_keys, this % prefix)

    if (len(optional_keys) > 0) then
      optional_nonset_keys = set_difference(optional_keys, config_dict)
      if (len(optional_nonset_keys) > 0) then
        call this % assign_dict(optional_nonset_keys)
        call announce_defaults(optional_nonset_keys, prefix = this % prefix)
      end if
    end if
    call this % check_values()
  end subroutine

end module
