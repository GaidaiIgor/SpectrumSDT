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
    character(:), allocatable :: stat
    integer :: K(2) = -1
    real(real64) :: rho(2) = -1
    real(real64) :: theta(2) = -1
    real(real64) :: phi(2) = -1

  contains
    procedure :: assign_dict => assign_dict_wf_section_params
    procedure :: get_optional_keys => get_optional_keys_wf_section_params
    procedure :: get_all_keys => get_all_keys_wf_section_params
    procedure :: check_values => check_values_wf_section_params
    procedure :: checked_init => checked_init_wf_section_params
    procedure :: checked_resolve_Ks => checked_resolve_Ks_wf_section_params
    procedure :: checked_resolve_rho_bounds => checked_resolve_rho_bounds_wf_section_params
    procedure :: checked_resolve_theta_bounds => checked_resolve_theta_bounds_wf_section_params
    procedure :: checked_resolve_phi_bounds => checked_resolve_phi_bounds_wf_section_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses range string into tokens.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_range(range_str) result(tokens)
    character(*), intent(in) :: range_str
    type(string) :: tokens(2)

    tokens = strsplit(range_str, '..')
    call tokens % trim()
    ! Transform keywords to placeholder values to be resolved later
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
    integer :: i, j
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
        case ('stat')
          this % stat = next_value
        case ('K')
          range_tokens = parse_range(next_value)
          this % K(1) = str2int(range_tokens(1) % to_char_str())
          this % K(2) = str2int(range_tokens(2) % to_char_str())
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
          do j = 1, 2
            this % phi(j) = str2real(range_tokens(j) % to_char_str())
            if (.not. (this % phi(j) .aeq. -2d0)) then
              this % phi(j) = this % phi(j) / 180 * pi
            end if
          end do
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
    call put_string(keys, 'stat', 'probability')
    call put_string(keys, 'K', 'start .. end')
    call put_string(keys, 'rho', 'start .. end')
    call put_string(keys, 'theta', 'start .. end')
    call put_string(keys, 'phi', 'start .. end')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_wf_section_params(this) result(keys)
    class(wf_section_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'name')
    call put_string(keys, 'stat')
    call put_string(keys, 'K')
    call put_string(keys, 'rho')
    call put_string(keys, 'theta')
    call put_string(keys, 'phi')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values. 
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_wf_section_params(this) 
    class(wf_section_params), intent(in) :: this
    call assert(any(this % stat == [character(100) :: 'probability', 'gamma']), 'Error: stat can be probability or gamma')
    ! Check of individual boundaries are performed later, when grid and K information is available
    call assert(any(this % K == -2) .or. this % K(1) <= this % K(2), 'Error: K(1) should be <= K(2)')
    call assert(any(this % rho .aeq. -2d0) .or. this % rho(1) < this % rho(2), 'Error: rho(1) should be < rho(2)')
    call assert(any(this % theta .aeq. -2d0) .or. this % theta(1) < this % theta(2), 'Error: theta(1) should be < theta(2)')
    call assert(any(this % phi .aeq. -2d0) .or. this % phi(1) < this % phi(2), 'Error: phi(1) should be < phi(2)')
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

!-------------------------------------------------------------------------------------------------------------------------------------------
! Replaces -2 placeholders with given values of *K_min* and *K_max* and checks that set values are within limits.
!-------------------------------------------------------------------------------------------------------------------------------------------
  impure elemental subroutine checked_resolve_Ks_wf_section_params(this, K_min, K_max)
    class(wf_section_params), intent(inout) :: this
    integer, intent(in) :: K_min, K_max

    if (this % K(1) == -2) then
      this % K(1) = K_min
    else
      call assert(this % K(1) >= K_min, 'Error: K(1) should be >= K_min')
    end if

    if (this % K(2) == -2) then
      this % K(2) = K_max
    else
      call assert(this % K(2) <= K_max, 'Error: K(2) should be <= K_max')
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Replaces -2 placeholders with given values of *rho_min* and *rho_max* and checks that set values are within limits.
!-------------------------------------------------------------------------------------------------------------------------------------------
  impure elemental subroutine checked_resolve_rho_bounds_wf_section_params(this, rho_min, rho_max)
    class(wf_section_params), intent(inout) :: this
    real(real64), intent(in) :: rho_min, rho_max

    if (this % rho(1) .aeq. -2d0) then
      this % rho(1) = rho_min
    else
      call assert(this % rho(1) .age. rho_min, 'Error: rho(1) should be >= rho_min')
    end if

    if (this % rho(2) .aeq. -2d0) then
      this % rho(2) = rho_max
    else
      call assert(this % rho(2) .ale. rho_max, 'Error: rho(2) should be <= rho_max')
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Replaces -2 placeholders with given values of *theta_min* and *theta_max* and checks that set values are within limits.
!-------------------------------------------------------------------------------------------------------------------------------------------
  impure elemental subroutine checked_resolve_theta_bounds_wf_section_params(this, theta_min, theta_max)
    class(wf_section_params), intent(inout) :: this
    real(real64), intent(in) :: theta_min, theta_max

    if (this % theta(1) .aeq. -2d0) then
      this % theta(1) = theta_min
    else
      call assert(this % theta(1) .age. theta_min, 'Error: theta(1) should be >= theta_min')
    end if

    if (this % theta(2) .aeq. -2d0) then
      this % theta(2) = theta_max
    else
      call assert(this % theta(2) .ale. theta_max, 'Error: theta(2) should be <= theta_max')
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Replaces -2 placeholders with given values of *phi_min* and *phi_max* and checks that set values are within limits.
!-------------------------------------------------------------------------------------------------------------------------------------------
  impure elemental subroutine checked_resolve_phi_bounds_wf_section_params(this, phi_min, phi_max)
    class(wf_section_params), intent(inout) :: this
    real(real64), intent(in) :: phi_min, phi_max

    if (this % phi(1) .aeq. -2d0) then
      this % phi(1) = phi_min
    else
      call assert(this % phi(1) .age. phi_min, 'Error: phi(1) should be >= phi_min')
    end if

    if (this % phi(2) .aeq. -2d0) then
      this % phi(2) = phi_max
    else
      call assert(this % phi(2) .ale. phi_max, 'Error: phi(2) should be <= phi_max')
    end if
  end subroutine

end module
