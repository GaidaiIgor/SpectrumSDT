module fixed_basis_params_mod
  use config_mod
  use constants
  use dictionary
  use dict_utils
  use general_utils
  use iso_fortran_env, only: real64
  use string_mod
  implicit none

  private
  public :: fixed_basis_params

  type :: fixed_basis_params
    character(:), allocatable :: prefix
    integer :: enabled = 0
    integer :: J = -1
    integer :: K = -1
    character(:), allocatable :: root_path

  contains
    procedure :: assign_dict => assign_dict_fixed_basis_params
    procedure :: get_mandatory_keys => get_mandatory_keys_fixed_basis_params
    procedure :: get_all_keys => get_all_keys_fixed_basis_params
    procedure :: check_values => check_values_fixed_basis_params
    procedure :: checked_init => checked_init_fixed_basis_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of fixed_basis_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_fixed_basis_params(this, config_dict)
    class(fixed_basis_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    integer :: i
    character(:), allocatable :: next_key, next_value
    type(string), allocatable :: key_set(:)

    this % enabled = 1
    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      next_value = extract_string(config_dict, next_key)
      select case (next_key)
        case ('prefix')
          this % prefix = next_value
        case ('J')
          this % J = str2int(next_value)
        case ('K')
          this % K = str2int(next_value)
        case ('root_path')
          this % root_path = next_value
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_fixed_basis_params(this) result(keys)
    class(fixed_basis_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'prefix')
    call put_string(keys, 'J')
    call put_string(keys, 'K')
    call put_string(keys, 'root_path')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_fixed_basis_params(this) result(keys)
    class(fixed_basis_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'prefix')
    call put_string(keys, 'J')
    call put_string(keys, 'K')
    call put_string(keys, 'root_path')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_fixed_basis_params(this) 
    class(fixed_basis_params), intent(in) :: this
    call assert(this % J >= 0, 'Error: ' // this % prefix // 'J should be >= 0')
    call assert(this % K >= 0, 'Error: ' // this % prefix // 'K should be >= 0')
    call assert(this % K <= this % J, 'Error: ' // this % prefix // 'K should be <= ' // this % prefix // 'J')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of fixed_basis_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_fixed_basis_params(this, config_dict)
    class(fixed_basis_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: mandatory_keys, all_keys

    call this % assign_dict(config_dict)
    mandatory_keys = this % get_mandatory_keys()
    call check_mandatory_keys(config_dict, mandatory_keys, this % prefix)

    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys, this % prefix)
    call this % check_values()
  end subroutine

end module
