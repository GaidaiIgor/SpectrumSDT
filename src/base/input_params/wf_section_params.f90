module wf_section_params_mod
  use config_mod
  use constants
  use dictionary
  use dict_utils
  use general_utils
  use iso_fortran_env, only: real64
  use string_mod
  use string_utils
  implicit none

  private
  public :: wf_section_params

  type :: wf_section_params
    character(:), allocatable :: name
    real(real64) :: rho(2) = -1
    real(real64) :: theta(2) = -1
    real(real64) :: phi(2) = -1
    integer :: K(2) = -1

  contains
    procedure :: assign_dict => assign_dict_wf_section_params
    procedure :: get_optional_keys => get_optional_keys_wf_section_params
    procedure :: get_all_keys => get_all_keys_wf_section_params
    procedure :: check_values => check_values_wf_section_params
    procedure :: checked_init => checked_init_wf_section_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses range string into tokens.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_range(range_str) result(tokens)
    character(*), intent(in) :: range_str
    type(string) :: tokens(2)

    tokens = strsplit(range_str, '..')
    ! Transform keywords to placeholder values
    if (tokens(1) == 'start') then
      tokens(1) = '-2'
    end if
    if (tokens(2) == 'end') then
      tokens(2) = '-2'
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Copies whatever settings it can recognize from the given dict. No checks other than parsing errors.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_wf_section_params(this, config_dict)
    class(wf_section_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    integer :: i
    character(:), allocatable :: next_key, next_value
    type(string) :: range_tokens(2)
    type(string), allocatable :: key_set(:)

    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      next_value = extract_string(config_dict, next_key)
      select case (next_key)
        case ('name')
          this % name = next_value
        case ('rho')
          range_tokens = parse_range(next_value)
          this % rho(1) = str2real(range_tokens(1) % to_char_str())
          this % rho(2) = str2real(range_tokens(2) % to_char_str())
        case ('theta')
          range_tokens = parse_range(next_value)
          this % theta(1) = str2real(range_tokens(1) % to_char_str())
          this % theta(2) = str2real(range_tokens(2) % to_char_str())
        case ('phi')
          range_tokens = parse_range(next_value)
          this % phi(1) = str2real(range_tokens(1) % to_char_str())
          this % phi(2) = str2real(range_tokens(2) % to_char_str())
        case ('K')
          range_tokens = parse_range(next_value)
          this % K(1) = str2int(range_tokens(1) % to_char_str())
          this % K(2) = str2int(range_tokens(2) % to_char_str())
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of optional keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_optional_keys_wf_section_params(this) result(keys)
    class(wf_section_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'name', 'use_key_name')
    call put_string(keys, 'rho', 'start .. end')
    call put_string(keys, 'theta', 'start .. end')
    call put_string(keys, 'phi', 'start .. end')
    call put_string(keys, 'K', 'start .. end')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_wf_section_params(this) result(keys)
    class(wf_section_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'name')
    call put_string(keys, 'rho')
    call put_string(keys, 'theta')
    call put_string(keys, 'phi')
    call put_string(keys, 'K')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values. 
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_wf_section_params(this) 
    class(wf_section_params), intent(in) :: this
    ! Check of individual boundaries are performed later, when grid and K information is available
    call assert(this % rho(1) < this % rho(2), 'Error: rho(1) should be < rho(2)')
    call assert(this % theta(1) < this % theta(2), 'Error: theta(1) should be < theta(2)')
    call assert(this % phi(1) < this % phi(2), 'Error: phi(1) should be < phi(2)')
    call assert(this % K(1) <= this % K(2), 'Error: K(1) should be <= K(2)')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of wf_section_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_wf_section_params(this, config_dict)
    class(wf_section_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: mandatory_keys, optional_keys, optional_nonset_keys, all_keys

    call this % assign_dict(config_dict)
    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys)
    optional_keys = this % get_optional_keys()
    call check_unused_keys(config_dict, mandatory_keys, optional_keys)

    if (len(optional_keys) > 0) then
      optional_nonset_keys = set_difference(optional_keys, config_dict)
      call this % assign_dict(optional_nonset_keys)
    end if
    call this % check_values()
  end subroutine

end module
