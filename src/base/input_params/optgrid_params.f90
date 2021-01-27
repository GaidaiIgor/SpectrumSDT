module optgrid_params_mod
  use config_mod
  use dictionary
  use dict_utils
  use grid_params_mod
  use iso_fortran_env, only: real64
  use parallel_utils
  use string_mod
  implicit none

  type, extends(grid_params) :: optgrid_params
    integer :: optimized = 0 ! optimized distribution of grid points
    character(:), allocatable :: envelope_path
    real(real64) :: max_energy = -1
    integer :: solver_steps = 2048

  contains
    procedure :: assign_dict => assign_dict_optgrid_params
    procedure :: get_mandatory_keys => get_mandatory_keys_optgrid_params
    procedure :: get_all_keys => get_all_keys_optgrid_params
    procedure :: set_defaults => set_defaults_optgrid_params
    procedure :: check_values => check_values_optgrid_params
    procedure :: checked_init => checked_init_optgrid_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of grid_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_optgrid_params(this, config_dict, auxiliary_info)
    class(optgrid_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    integer :: i
    character(:), allocatable :: next_key, full_key, next_value
    type(string), allocatable :: key_set(:)

    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      full_key = this % prefix // next_key
      next_value = extract_string(config_dict, next_key)
      select case (next_key)
        case ('optimized')
          this % optimized = str2int_config(next_value, full_key)
        case ('envelope_path')
          this % envelope_path = next_value
        case ('max_energy')
          this % max_energy = str2real_config(next_value, full_key)
        case ('solver_steps')
          this % solver_steps = str2int_config(next_value, full_key)
      end select
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_optgrid_params(this) result(keys)
    class(optgrid_params), intent(in) :: this
    type(dictionary_t) :: keys

    keys = this % grid_params % get_mandatory_keys()
    if (this % optimized == 1) then
      call put_string(keys, 'envelope_path')
      call put_string(keys, 'max_energy')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_optgrid_params(this) result(keys)
    class(optgrid_params), intent(in) :: this
    type(dictionary_t) :: keys

    keys = this % grid_params % get_all_keys()
    call put_string(keys, 'optimized')
    call put_string(keys, 'envelope_path')
    call put_string(keys, 'max_energy')
    call put_string(keys, 'solver_steps')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets defaults.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_defaults_optgrid_params(this, config_dict)
    class(optgrid_params), intent(inout) :: this
    class(dictionary_t), intent(in) :: config_dict

    call this % grid_params % set_defaults(config_dict)
    if (.not. ('optimized' .in. config_dict)) then
      call print_parallel(this % prefix // 'optimized is not specified, assuming equidistant grid')
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_optgrid_params(this)
    class(optgrid_params), intent(in) :: this
    call assert(any(this % optimized == [0, 1]), 'Error: ' // this % prefix // 'optimized should be 0 or 1')
    call assert(this % solver_steps > 0, 'Error: ' // this % prefix // 'solver_steps should be > 0')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of grid_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_optgrid_params(this, config_dict, auxiliary_info, check_extra, check_values)
    class(optgrid_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    integer, optional, intent(in) :: check_extra, check_values ! have to keep the same signature in child
    type(dictionary_t) :: mandatory_keys, all_keys

    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys, this % prefix) ! Checks that unknown keys were not specified
    call check_key_types(config_dict, auxiliary_info, 'string')
    call this % grid_params % checked_init(config_dict, auxiliary_info, check_extra = 0) ! First, check init parental type
    call this % assign_dict(config_dict, auxiliary_info) ! Init type as is to ease branching in getting keys
    mandatory_keys = this % get_mandatory_keys() ! Save mandatory keys because we will reuse them later for unused key check
    call check_mandatory_keys(config_dict, mandatory_keys, this % prefix) ! Make sure mandatory keys are set to ease checking values
    call this % set_defaults(config_dict)
    call this % check_values() ! Make sure the settings have valid values
  end subroutine

end module
