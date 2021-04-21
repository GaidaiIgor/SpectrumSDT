module basis_params_mod
  use config_mod
  use constants, only: au_to_wn
  use dictionary
  use dict_utils
  use fixed_basis_params_mod
  use general_utils
  use iso_fortran_env, only: real64
  use string_mod
  implicit none

  private
  public :: basis_params

  type :: basis_params
    character(:), allocatable :: prefix
    integer :: num_functions_phi = -1 ! number of sines or cosines in basis of 1D problem
    integer :: symmetry = -1 ! 0 (even, cos, +) or 1 (odd, sin, -). In case of coupled hamiltonian means symmetry of K=0, even when K=0 is not included.
    real(real64) :: cutoff_energy_1d = 0 ! 1D solutions with energies higher than this are discarded from basis
    real(real64) :: cutoff_energy_2d = 0 ! 2D solutions with energies higher than this are discarded from basis
    integer :: min_solutions_1d = 1 ! minimum number of 1D solutions in each slice, kept even if energies are higher than cutoff_energy_1d
    integer :: min_solutions_2d = 1 ! minimum number of 2D solutions in each slice, kept even if energies are higher than cutoff_energy_2d
    integer :: print_energies_1d = 0 ! whether or not to create a file with formatted 1D energies
    integer :: print_energies_2d = 0 ! whether or not to create a file with formatted 2D energies
    type(fixed_basis_params) :: fixed ! parameters describing rotational state of fixed basis, if used

  contains
    procedure, nopass :: check_key_types_basis_params
    procedure :: assign_dict => assign_dict_basis_params
    procedure :: get_mandatory_keys => get_mandatory_keys_basis_params
    procedure :: get_all_keys => get_all_keys_basis_params
    procedure :: set_defaults => set_defaults_basis_params
    procedure :: check_values => check_values_basis_params
    procedure :: checked_init => checked_init_basis_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Makes sure all keys have correct types.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_key_types_basis_params(config_dict, auxiliary_info)
    class(dictionary_t), intent(in) :: config_dict, auxiliary_info
    type(dictionary_t) :: dict_types

    call put_string(dict_types, 'fixed', 'dict')
    call check_key_types(config_dict, auxiliary_info, 'string', dict_types)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of basis_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_basis_params(this, config_dict, auxiliary_info)
    class(basis_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    integer :: i
    character(:), allocatable :: next_key, full_key, key_type, next_value
    type(string), allocatable :: key_set(:)
    type(dictionary_t) :: subdict, auxiliary_subdict

    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      full_key = this % prefix // next_key
      key_type = extract_string(auxiliary_info, next_key // '_type')

      select case (key_type)
        case ('string')
          next_value = extract_string(config_dict, next_key)
        case ('dict')
          call associate(subdict, config_dict, next_key)
          call associate(auxiliary_subdict, auxiliary_info, next_key)
        case default
          stop 'Error: wrong key type'
      end select

      select case (next_key)
        case ('num_functions_phi')
          this % num_functions_phi = str2int_config(next_value, full_key)
        case ('symmetry')
          this % symmetry = str2int_config(next_value, full_key)
        case ('cutoff_energy_1d')
          this % cutoff_energy_1d = str2real_config(next_value, full_key) / au_to_wn
        case ('cutoff_energy_2d')
          this % cutoff_energy_2d = str2real_config(next_value, full_key) / au_to_wn
        case ('min_solutions_1d')
          this % min_solutions_1d = str2int_config(next_value, full_key)
        case ('min_solutions_2d')
          this % min_solutions_2d = str2int_config(next_value, full_key)
        case ('print_energies_1d')
          this % print_energies_1d = str2int_config(next_value, full_key)
        case ('print_energies_2d')
          this % print_energies_2d = str2int_config(next_value, full_key)
        case ('fixed')
          call this % fixed % checked_init(subdict, auxiliary_subdict)
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_basis_params(this) result(keys)
    class(basis_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'num_functions_phi')
    call put_string(keys, 'symmetry')
    call put_string(keys, 'min_solutions_2d')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_basis_params(this) result(keys)
    class(basis_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'num_functions_phi')
    call put_string(keys, 'symmetry')
    call put_string(keys, 'cutoff_energy_1d')
    call put_string(keys, 'cutoff_energy_2d')
    call put_string(keys, 'min_solutions_1d')
    call put_string(keys, 'min_solutions_2d')
    call put_string(keys, 'print_energies_1d')
    call put_string(keys, 'print_energies_2d')
    call put_string(keys, 'fixed')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets default values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_defaults_basis_params(this, config_dict, stage)
    class(basis_params), intent(inout) :: this
    class(dictionary_t), intent(in) :: config_dict
    character(*), intent(in) :: stage

    if (.not. ('min_solutions_1d' .in. config_dict) .and. stage == 'basis') then
      this % min_solutions_1d = this % min_solutions_2d
      call print_parallel(this % prefix // 'min_functions_1d is not specified. Setting it equal to ' // this % prefix // 'min_functions_2d (likely an overkill).')
    end if
    if (.not. ('cutoff_energy_1d' .in. config_dict) .and. stage == 'basis') then
      call print_parallel(this % prefix // 'cutoff_energy_1d is not specified. Assuming a default value of ' // num2str(this % cutoff_energy_1d) // '.')
    end if
    if (.not. ('cutoff_energy_2d' .in. config_dict) .and. stage == 'basis') then
      call print_parallel(this % prefix // 'cutoff_energy_2d is not specified. Assuming a default value of ' // num2str(this % cutoff_energy_2d) // '.')
    end if
    if (.not. ('fixed' .in. config_dict) .and. any(stage == [character(100) :: 'overlaps', 'eigensolve', 'properties'])) then
      call print_parallel(this % prefix // 'fixed is not specified. Assuming adiabatic basis.')
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_basis_params(this) 
    class(basis_params), intent(in) :: this
    call assert(this % num_functions_phi > 0, 'Error: ' // this % prefix // 'num_functions_phi should be > 0')
    call assert(any(this % symmetry == [0, 1]), 'Error: ' // this % prefix // 'symmery should be 0 or 1')
    call assert(this % min_solutions_1d > 0, 'Error: ' // this % prefix // 'min_solutions_1d should be > 0')
    call assert(this % min_solutions_2d > 0, 'Error: ' // this % prefix // 'min_solutions_2d should be > 0')
    call assert(any(this % print_energies_1d == [0, 1]), 'Error: ' // this % prefix // 'print_energies_1d should be 0 or 1')
    call assert(any(this % print_energies_2d == [0, 1]), 'Error: ' // this % prefix // 'print_energies_2d should be 0 or 1')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of basis_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_basis_params(this, config_dict, auxiliary_info, stage)
    class(basis_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    character(*), intent(in) :: stage
    type(dictionary_t) :: mandatory_keys, all_keys

    this % prefix = extract_string(auxiliary_info, 'prefix')
    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys, this % prefix)
    call check_key_types_basis_params(config_dict, auxiliary_info)
    call this % assign_dict(config_dict, auxiliary_info)
    mandatory_keys = this % get_mandatory_keys()
    call check_mandatory_keys(config_dict, mandatory_keys, this % prefix)
    call this % set_defaults(config_dict, stage)
    call this % check_values()
  end subroutine

end module
