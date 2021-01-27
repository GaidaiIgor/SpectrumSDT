module input_params_mod
  use algorithms_mod
  use cap_params_mod
  use config_mod
  use constants
  use dictionary
  use dict_utils
  use eigensolve_params_mod
  use fixed_basis_params_mod
  use general_utils
  use grid_info_mod
  use grid_params_mod
  use optgrid_params_mod
  use iso_fortran_env, only: real64
  use parallel_utils
  use rovib_utils_base_mod
  use string_mod
  use string_utils
  use wf_section_params_mod
  implicit none

  private
  public :: input_params

  type :: input_params
    character(:), allocatable :: prefix

    ! Behavior control
    character(:), allocatable :: stage ! grids, basis, overlaps, eigensolve or properties
    integer :: use_rovib_coupling = -1 ! enables/disables use_rovib_coupling coupling

    ! Grids
    type(optgrid_params) :: grid_rho
    type(grid_params) :: grid_theta
    integer :: num_points_phi = -1
    character(:), allocatable :: output_coordinate_system

    ! System
    real(real64) :: mass(3) = -1
    integer :: J = -1 ! total angular momentum quantum number
    integer :: K(2) = -1 ! boundaries of K-range for use_rovib_coupling calculation. In symmetric top rotor both values should be the same
    integer :: parity = -1 ! 0(+) or 1(-)
    integer :: symmetry = -1 ! 0 (even, cos, +) or 1 (odd, sin, -). In case of coupled hamiltonian means symmetry of K=0, even when K=0 is not included.

    ! Basis
    integer :: basis_size_phi = -1 ! number of sines or cosines for 1D step
    real(real64) :: cutoff_energy = -1 ! solutions with energies higher than this are discarded from basis
    type(fixed_basis_params) :: fixed_basis

    ! Eigensolve
    type(eigensolve_params) :: eigensolve ! Eigensolver parameters
    type(cap_params) :: cap ! Complex Absorbin Potential parameters

    ! Properties
    type(wf_section_params), allocatable :: wf_sections(:) ! wave function sections for integration on the properties stage
    
    ! Paths
    character(:), allocatable :: grid_path ! path to folder with grid calculations
    character(:), allocatable :: root_path ! path to root folder for main calculations

    ! Debug
    integer :: use_parallel = -1 ! parallel execution
    integer :: enable_terms(2) = 1 ! 1st digit - coriolis, 2nd - asymmetric
    integer :: optimized_mult = 1 ! disables matrix-vector multiplication optimizations
    integer :: treat_tp_as_xy = 0 ! treats theta and phi grids as x and y grids for aph plots
    character(:), allocatable :: debug_mode
    character(:), allocatable :: debug_param_1

  contains
    procedure, nopass :: check_key_types_input_params
    procedure :: assign_dict => assign_dict_input_params
    procedure :: get_mandatory_keys => get_mandatory_keys_input_params
    procedure :: get_all_keys => get_all_keys_input_params
    procedure :: set_defaults => set_defaults_input_params
    procedure :: check_values => check_values_input_params
    procedure :: checked_init => checked_init_input_params
    procedure :: check_resolve_grids => check_resolve_grids_input_params
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Makes sure all keys have correct types.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_key_types_input_params(config_dict, auxiliary_info)
    class(dictionary_t), intent(in) :: config_dict, auxiliary_info
    type(dictionary_t) :: dict_types

    call put_string(dict_types, 'grid_rho', 'dict')
    call put_string(dict_types, 'grid_theta', 'dict')
    call put_string(dict_types, 'fixed_basis', 'dict')
    call put_string(dict_types, 'eigensolve', 'dict')
    call put_string(dict_types, 'cap', 'dict')
    call put_string(dict_types, 'wf_sections', 'dict')
    call check_key_types(config_dict, auxiliary_info, 'string', dict_types)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts all fields in *params* from degrees to radians.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine convert_grid_params_deg_to_rad(params)
    class(grid_params), intent(inout) :: params
    params % from = params % from / rad_to_deg
    params % to = params % to / rad_to_deg
    if (.not. (params % step .aeq. -1d0)) then
      params % step = params % step / rad_to_deg
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses user-provded value of mass.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_mass(mass_str) result(mass)
    character(*), intent(in) :: mass_str
    real(real64) :: mass(3)
    integer :: i
    character(:), allocatable :: next_atom
    type(string), allocatable :: tokens(:)

    tokens = strsplit(mass_str, delim = ',')
    call assert(size(tokens) == 3, 'Error: mass shoud specify 3 comma separated values')
    do i = 1, size(mass)
      next_atom = adjustl(trim(tokens(i) % to_char_str()))
      select case (next_atom)
        case ('H1')
          mass(i) = hydrogen_masses(1)
        case ('H2')
          mass(i) = hydrogen_masses(2)
        case ('H3')
          mass(i) = hydrogen_masses(3)
        case ('He3')
          mass(i) = helium_masses(1)
        case ('He4')
          mass(i) = helium_masses(2)
        case ('Li6')
          mass(i) = lithium_masses(1)
        case ('Li7')
          mass(i) = lithium_masses(2)
        case ('Be9')
          mass(i) = beryllium_masses(1)
        case ('B10')
          mass(i) = boron_masses(1)
        case ('B11')
          mass(i) = boron_masses(2)
        case ('C11')
          mass(i) = carbon_masses(1)
        case ('C12')
          mass(i) = carbon_masses(2)
        case ('C13')
          mass(i) = carbon_masses(3)
        case ('C14')
          mass(i) = carbon_masses(4)
        case ('N14')
          mass(i) = nitrogen_masses(1)
        case ('N15')
          mass(i) = nitrogen_masses(2)
        case ('O16')
          mass(i) = oxygen_masses(1)
        case ('O17')
          mass(i) = oxygen_masses(2)
        case ('O18')
          mass(i) = oxygen_masses(3)
        case ('F18')
          mass(i) = fluorine_masses(1)
        case ('F19')
          mass(i) = fluorine_masses(2)
        case ('Ne20')
          mass(i) = neon_masses(1)
        case ('Ne21')
          mass(i) = neon_masses(2)
        case ('Ne22')
          mass(i) = neon_masses(3)
        case ('Na22')
          mass(i) = sodium_masses(1)
        case ('Na23')
          mass(i) = sodium_masses(2)
        case ('Na24')
          mass(i) = sodium_masses(3)
        case ('Mg24')
          mass(i) = magnesium_masses(1)
        case ('Mg25')
          mass(i) = magnesium_masses(2)
        case ('Mg26')
          mass(i) = magnesium_masses(3)
        case ('Al27')
          mass(i) = aluminium_masses(1)
        case ('Si28')
          mass(i) = silicon_masses(1)
        case ('Si29')
          mass(i) = silicon_masses(2)
        case ('Si30')
          mass(i) = silicon_masses(3)
        case ('P31')
          mass(i) = phosphorus_masses(1)
        case ('P32')
          mass(i) = phosphorus_masses(2)
        case ('S32')
          mass(i) = sulfur_masses(1)
        case ('S33')
          mass(i) = sulfur_masses(2)
        case ('S34')
          mass(i) = sulfur_masses(3)
        case ('S35')
          mass(i) = sulfur_masses(4)
        case ('S36')
          mass(i) = sulfur_masses(5)
        case ('Cl35')
          mass(i) = chlorine_masses(1)
        case ('Cl37')
          mass(i) = chlorine_masses(2)
        case ('Ar36')
          mass(i) = argon_masses(1)
        case ('Ar38')
          mass(i) = argon_masses(2)
        case ('Ar40')
          mass(i) = argon_masses(3)
        case default
          call print_parallel('Info: treating mass as a number')
          mass(i) = str2real_config(next_atom, 'mass' // '(' // num2str(i) // ')') * amu_to_aum
      end select
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses given *dict* into wf_sections. Any item in the *dict* has to be a dict representing an individual section.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_wf_sections(dict, auxiliary_info) result(wf_sections)
    class(dictionary_t) :: dict, auxiliary_info ! intent(in)
    type(wf_section_params), allocatable :: wf_sections(:)
    integer :: i
    integer, allocatable :: section_indices(:), key_order(:)
    character(:), allocatable :: next_key
    type(string), allocatable :: key_set(:)
    type(dictionary_t) :: subdict, auxiliary_subdict

    call check_key_types(dict, auxiliary_info, 'dict')
    key_set = get_key_set(dict)
    allocate(section_indices(size(key_set)))
    ! We need to iterate the sections in order they were given.
    ! Order of keys in dict can be arbitrary, so we need to use dict index stored from parsing.
    do i = 1, size(key_set)
      call associate(auxiliary_subdict, auxiliary_info, key_set(i) % to_char_str())
      call assign(section_indices(i), auxiliary_subdict, 'dict_index')
    end do
    section_indices = bubble_sort(section_indices, index_permutation = key_order)

    allocate(wf_sections(size(key_set)))
    do i = 1, size(key_order)
      next_key = key_set(key_order(i)) % to_char_str()
      call associate(subdict, dict, next_key)
      call associate(auxiliary_subdict, auxiliary_info, next_key)

      ! use key name if name attribute is not provided
      if (.not. ('name' .in. subdict)) then
        call put_string(subdict, 'name', next_key)
        call put_string(auxiliary_subdict, 'name_type', 'string')
      end if
      call wf_sections(i) % checked_init(subdict, auxiliary_subdict)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses user-provided value of `enable_terms`.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_enable_terms(enable_terms_str) result(enable_terms)
    character(*), intent(in) :: enable_terms_str
    integer :: enable_terms(2)
    
    enable_terms(1) = str2int_config(enable_terms_str(1:1), 'enable_terms(1)')
    enable_terms(2) = str2int_config(enable_terms_str(2:2), 'enable_terms(2)')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses user-provided value of `K`.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_K(K_str, J, parity) result(K)
    character(*), intent(in) :: K_str
    integer, intent(in) :: J, parity
    integer :: K(2)
    integer :: pos
    type(string), allocatable :: tokens(:)

    if (K_str == 'all') then
      K(1) = get_k_start(J, parity)
      K(2) = J
    else
      pos = index(K_str, '..')
      if (pos == 0) then
        ! Single value of K
        K(1) = str2int_config(K_str, 'K')
        K(2) = K(1)
      else
        ! K range
        tokens = strsplit(K_str, '..')
        K(1) = str2int_config(tokens(1) % to_char_str(), 'K(1)')
        K(2) = str2int_config(tokens(2) % to_char_str(), 'K(2)')
      end if
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of input_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_input_params(this, config_dict, auxiliary_info)
    class(input_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    logical :: wf_sections_provided
    integer :: i, check_values
    character(:), allocatable :: next_key, key_type, next_value, K_str
    type(string), allocatable :: key_set(:)
    type(dictionary_t) :: subdict, auxiliary_subdict

    wf_sections_provided = .false.
    this % prefix = extract_string(auxiliary_info, 'prefix')
    ! Iterate over the keys given by user
    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      key_type = extract_string(auxiliary_info, next_key // '_type')

      if (key_type == 'string') then
        next_value = extract_string(config_dict, next_key)
      else if (key_type == 'dict') then
        call associate(subdict, config_dict, next_key)
        call associate(auxiliary_subdict, auxiliary_info, next_key)
      else
        stop 'Error: wrong key type'
      end if

      select case (next_key)
        case ('stage')
          this % stage = next_value
        case ('use_rovib_coupling')
          this % use_rovib_coupling = str2int_config(next_value, next_key)
        case ('grid_rho')
          call this % grid_rho % checked_init(subdict, auxiliary_subdict)
        case ('grid_theta')
          check_values = iff('treat_tp_as_xy' .in. config_dict, 0, 1)
          call this % grid_theta % checked_init(subdict, auxiliary_subdict, check_values = check_values)
          if (check_values == 1) then
            call convert_grid_params_deg_to_rad(this % grid_theta)
          end if
        case ('num_points_phi')
          this % num_points_phi = str2int_config(next_value, next_key)
        case ('output_coordinate_system')
          this % output_coordinate_system = next_value
        case ('mass')
          this % mass = parse_mass(next_value)
        case ('J')
          this % J = str2int_config(next_value, next_key)
        case ('K')
          K_str = next_value
        case ('parity')
          this % parity = str2int_config(next_value, next_key)
        case ('symmetry')
          this % symmetry = str2int_config(next_value, next_key)
        case ('basis_size_phi')
          this % basis_size_phi = str2int_config(next_value, next_key)
        case ('cutoff_energy')
          this % cutoff_energy = str2real_config(next_value, next_key) / au_to_wn
        case ('fixed_basis')
          call this % fixed_basis % checked_init(subdict, auxiliary_subdict)
        case ('eigensolve')
          call this % eigensolve % checked_init(subdict, auxiliary_subdict)
        case ('cap')
          call this % cap % checked_init(subdict, auxiliary_subdict)
        case ('wf_sections')
          wf_sections_provided = .true.
          this % wf_sections = parse_wf_sections(subdict, auxiliary_subdict)
        case ('grid_path')
          this % grid_path = next_value
        case ('root_path')
          this % root_path = next_value
        case ('use_parallel')
          this % use_parallel = str2int_config(next_value, next_key)
        case ('enable_terms')
          this % enable_terms = parse_enable_terms(next_value)
        case ('optimized_mult')
          this % optimized_mult = str2int_config(next_value, next_key)
        case ('treat_tp_as_xy')
          this % treat_tp_as_xy = str2int_config(next_value, next_key)
        case ('debug_mode')
          this % debug_mode = next_value
        case ('debug_param_1')
          this % debug_param_1 = next_value
      end select
    end do

    ! Values whose interpretation depends on other values are initialized after the main set loop
    if (allocated(K_str)) then
      this % K = parse_K(K_str, this % J, this % parity)
    end if
    if (wf_sections_provided) then
      call this % wf_sections % checked_resolve_Ks(this % K(1), this % K(2))
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of mandatory keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mandatory_keys_input_params(this) result(keys)
    class(input_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'stage')
    if (this % stage == 'grids') then
      call put_string(keys, 'grid_rho')
      call put_string(keys, 'grid_theta')
      call put_string(keys, 'num_points_phi')
    end if

    if (this % stage == 'grids' .and. (this % grid_rho % optimized == 1 .or. any(this % output_coordinate_system == [character(100) :: 'jacobi', 'cartesian', 'all bonds', 'internal'])) .or. &
        this % stage == 'basis' .or. this % stage == 'overlaps' .or. this % stage == 'eigensolve' .or. this % stage == 'properties') then
      call put_string(keys, 'mass')
    end if

    if (any(this % stage == [character(100) :: 'basis', 'overlaps', 'eigensolve', 'properties'])) then
      call put_string(keys, 'use_rovib_coupling')
      call put_string(keys, 'J')
      call put_string(keys, 'K')
      call put_string(keys, 'symmetry')
      if (this % use_rovib_coupling == 1) then
        call put_string(keys, 'parity')
      end if

      call put_string(keys, 'cutoff_energy')
      call put_string(keys, 'basis_size_phi')
      call put_string(keys, 'grid_path')
      call put_string(keys, 'root_path')
    end if

    if (any(this % stage == [character(100) :: 'eigensolve', 'properties'])) then
      call put_string(keys, 'eigensolve')
    end if

    if (this % stage == 'properties') then
      call put_string(keys, 'wf_sections')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_input_params(this) result(keys)
    class(input_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'stage')
    call put_string(keys, 'use_rovib_coupling')
    call put_string(keys, 'grid_rho')
    call put_string(keys, 'grid_theta')
    call put_string(keys, 'num_points_phi')
    call put_string(keys, 'output_coordinate_system')
    call put_string(keys, 'mass')
    call put_string(keys, 'J')
    call put_string(keys, 'K')
    call put_string(keys, 'parity')
    call put_string(keys, 'symmetry')
    call put_string(keys, 'basis_size_phi')
    call put_string(keys, 'cutoff_energy')
    call put_string(keys, 'fixed_basis')
    call put_string(keys, 'eigensolve')
    call put_string(keys, 'cap')
    call put_string(keys, 'wf_sections')
    call put_string(keys, 'grid_path')
    call put_string(keys, 'root_path')
    call put_string(keys, 'use_parallel')
    call put_string(keys, 'enable_terms')
    call put_string(keys, 'optimized_mult')
    call put_string(keys, 'treat_tp_as_xy')
    call put_string(keys, 'debug_mode')
    call put_string(keys, 'debug_param_1')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets default values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_defaults_input_params(this, config_dict)
    class(input_params), intent(inout) :: this
    class(dictionary_t), intent(in) :: config_dict

    if (.not. ('output_coordinate_system' .in. config_dict)) then
      this % output_coordinate_system = 'aph'
      if (this % stage == 'grids') then
        call print_parallel('output_coordinate_system is not specified, assuming aph')
      end if
    end if

    if (.not. ('use_parallel' .in. config_dict)) then
      this % use_parallel = iff(this % stage == 'grids', 0, 1)
    end if

    if (.not. ('cap' .in. config_dict)) then
      call this % cap % init_default()
      if (this % stage == 'eigensolve') then
        call print_parallel('cap is not specified, assuming no cap')
      end if
    end if

    if (.not. ('fixed_basis' .in. config_dict)) then
      if (any(this % stage == [character(100) :: 'overlaps', 'eigensolve', 'properties'])) then
        call print_parallel('fixed_basis is not specified, assuming adiabatic basis')
      end if
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values. -1 means the value was not provided.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_input_params(this) 
    class(input_params), intent(in) :: this
    call assert(any(this % stage == [character(100) :: 'grids', 'basis', 'overlaps', 'eigensolve', 'properties']), 'Error: stage can be "grids", "basis", "overlaps", "eigensolve" or "properties"')
    call assert(any(this % use_rovib_coupling == [-1, 0, 1]), 'Error: use_rovib_coupling can be 0 or 1')
    if (this % treat_tp_as_xy /= 1) then
      call assert((this % grid_theta % from .aeq. -1d0) .or. (this % grid_theta % from .ale. pi/2), 'Error: grid_theta % from should be <= pi/2')
      call assert((this % grid_theta % to .aeq. -1d0) .or. (this % grid_theta % to .ale. pi/2), 'Error: grid_theta % to should be <= pi/2')
    end if
    call assert(this % num_points_phi == -1 .or. this % num_points_phi > 1, 'Error: num_points_phi should be > 1')
    call assert(any(this % output_coordinate_system == [character(100) :: 'aph', 'mass jacobi', 'jacobi', 'cartesian', 'all bonds', 'internal']), &
        'Error: output_coordinate_system can be "aph", "mass jacobi", "jacobi", "cartesian", "all bonds" or "internal"')
    call assert(all(this % mass .aeq. -1d0) .or. all(this % mass > 0), 'Error: all mass should be > 0')
    call assert(this % J >= -1, 'Error: J should be >= 0')
    call assert(any(this % parity == [-1, 0, 1]), 'Error: parity value can be 0 or 1')
    if (any(this % stage == [character(len = 100) :: 'eigensolve', 'properties']) .and. this % use_rovib_coupling == 1) then
      call assert(all(this % K == -1) .or. all(this % K >= get_k_start(this % J, this % parity)), 'Error: K should be >= mod(J+p, 2)')
    else
      call assert(all(this % K >= -1), 'Error: K should be >= 0')
    end if
    call assert(all(this % K <= this % J), 'Error: K should be <= J')
    call assert(any(this % symmetry == [-1, 0, 1]), 'Error: symmery value can be 0 or 1')
    call assert(this % basis_size_phi == -1 .or. this % basis_size_phi > 0, 'Error: basis_size_phi should be > 0')
    call assert(any(this % use_parallel == [-1, 0, 1]), 'Error: use_parallel can be 0 or 1')
    call assert(any(this % enable_terms(1) == [-1, 0, 1]), 'Error: enable_terms(1) can be 0 or 1')
    call assert(any(this % enable_terms(2) == [-1, 0, 1]), 'Error: enable_terms(2) can be 0 or 1')
    call assert(any(this % optimized_mult == [-1, 0, 1]), 'Error: optimized_mult can be 0 or 1')
    call assert(any(this % treat_tp_as_xy == [-1, 0, 1]), 'Error: treat_tp_as_xy can be 0 or 1')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of input_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_input_params(this, config_dict, auxiliary_info)
    class(input_params), intent(inout) :: this
    class(dictionary_t) :: config_dict, auxiliary_info ! intent(in)
    type(dictionary_t) :: mandatory_keys, all_keys

    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys)
    call check_key_types_input_params(config_dict, auxiliary_info)
    call this % assign_dict(config_dict, auxiliary_info)
    mandatory_keys = this % get_mandatory_keys()
    call check_mandatory_keys(config_dict, mandatory_keys)
    call this % set_defaults(config_dict)
    call this % check_values()
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Resolves parameters that require grid information and performs grid related checks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_resolve_grids_input_params(this, rho_info, theta_info, phi_info)
    class(input_params), intent(inout) :: this
    class(grid_info), intent(in) :: rho_info, theta_info, phi_info

    call assert(this % basis_size_phi == -1 .or. this % basis_size_phi <= size(phi_info % points) / 2, 'Error: basis_size_phi should be <= num_points_phi / 2')
    if (this % stage == 'properties') then
      call this % wf_sections % checked_resolve_rho_bounds(rho_info % from, rho_info % to)
      call this % wf_sections % checked_resolve_theta_bounds(theta_info % from, theta_info % to)
      call this % wf_sections % checked_resolve_phi_bounds(0d0, 2*pi)
    end if
  end subroutine

end module
