module cap_params_mod
  use config_mod
  use constants
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
    real(real64) :: min_absorbed_energy = -1

  contains
    procedure :: init_default => init_default_cap_params
    procedure :: assign_dict => assign_dict_cap_params
    procedure :: get_mandatory_keys => get_mandatory_keys_cap_params
    procedure :: get_all_keys => get_all_keys_cap_params
    procedure :: check_values => check_values_cap_params
    procedure :: checked_init => checked_init_cap_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Inits default values that cannot be statically initialized.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_default_cap_params(this)
    class(cap_params), intent(inout) :: this
    this % type = 'none'
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of cap_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_cap_params(this, config_dict, auxiliary_info)
    class(cap_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    integer :: i
    character(:), allocatable :: next_key, next_value
    type(string), allocatable :: key_set(:)

    this % prefix = extract_string(auxiliary_info, 'prefix')
    this % type = 'Manolopoulos'
    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      next_value = extract_string(config_dict, next_key)
      select case (next_key)
        case ('min_absorbed_energy')
          this % min_absorbed_energy = str2real(next_value) / au_to_wn
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_cap_params(this) result(keys)
    class(cap_params), intent(in) :: this
    type(dictionary_t) :: keys
    call put_string(keys, 'min_absorbed_energy')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_cap_params(this) result(keys)
    class(cap_params), intent(in) :: this
    type(dictionary_t) :: keys
    call put_string(keys, 'min_absorbed_energy')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_cap_params(this) 
    class(cap_params), intent(in) :: this
    call assert(this % min_absorbed_energy > 0, 'Error: ' // this % prefix // 'min_absorbed_energy should be > 0')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of cap_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_cap_params(this, config_dict, auxiliary_info)
    class(cap_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    type(dictionary_t) :: mandatory_keys, all_keys

    call check_key_types(config_dict, auxiliary_info, 'string')
    call this % assign_dict(config_dict, auxiliary_info)
    mandatory_keys = this % get_mandatory_keys()
    call check_mandatory_keys(config_dict, mandatory_keys, this % prefix)

    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys, this % prefix)
    call this % check_values()
  end subroutine

end module
