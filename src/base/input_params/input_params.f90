module input_params_mod
  use algorithms_mod
  use config_mod
  use constants
  use dictionary
  use dict_utils
  use general_utils
  use grid_params_mod
  use optgrid_params_mod
  use iso_fortran_env, only: real64
  use rovib_utils_base_mod
  use string_mod
  use string_utils
  use wf_section_params_mod
  implicit none

  private
  public :: input_params

  type :: input_params
    ! Behavior control
    character(:), allocatable :: stage ! grids, basis, overlaps, eigencalc or properties
    integer :: use_rovib_coupling = -1 ! enables/disables use_rovib_coupling coupling
    integer :: use_fixed_basis_JK = -1 ! use basis set with the same fixed values of J and K for all calculations
    character(:), allocatable :: cap_type ! type of Complex Absorbing Potential

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
    character(:), allocatable :: basis_root_path ! path to root folder for basis calculations
    integer :: basis_J = -1 ! J of basis set
    integer :: basis_K = -1 ! K of basis set

    ! Eigencalc
    integer :: num_states = -1 ! desired number of eigenstates from eigencalc
    integer :: ncv = -1 ! number of vectors in arnoldi basis during eigencalc
    integer :: mpd = -1 ! maximum projected dimension, slepc only
    integer :: max_iterations = -1 ! maximum number of iterations during eigencalc

    ! Properties
    type(wf_section_params), allocatable :: wf_sections(:) ! wave function sections for integration on the properties stage
    
    ! Misc paths
    character(:), allocatable :: grid_path ! path to folder with grid calculations
    character(:), allocatable :: root_path ! path to root folder for main calculations

    ! Debug
    integer :: use_parallel = -1 ! parallel execution
    integer :: enable_terms(2) = -1 ! 1st digit - coriolis, 2nd - asymmetric
    integer :: optimized_mult = -1 ! disables matrix-vector multiplication optimizations
    character(:), allocatable :: debug_mode
    character(:), allocatable :: debug_param_1

  contains
    procedure :: assign_dict => assign_dict_input_params
    procedure :: get_mandatory_keys => get_mandatory_keys_input_params
    procedure :: get_optional_keys => get_optional_keys_input_params
    procedure :: get_all_keys => get_all_keys_input_params
    procedure :: check_values => check_values_input_params
    procedure :: checked_init => checked_init_input_params
    procedure :: check_resolve_grids => check_resolve_grids_input_params
  end type

contains

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
! Parses user-provded value of mass
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
        case ('o16')
          mass(i) = oxygen_masses(1)
        case ('o17')
          mass(i) = oxygen_masses(2)
        case ('o18')
          mass(i) = oxygen_masses(3)
        case default
          mass(i) = str2real(next_atom) * amu_to_aum
      end select
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses given *dict* into wf_sections. Any item in the *dict* has to be a dict representing an individual section.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_wf_sections(dict) result(wf_sections)
    class(dictionary_t) :: dict ! intent(in)
    type(wf_section_params), allocatable :: wf_sections(:)
    integer :: i
    integer, allocatable :: section_indices(:), key_order(:)
    character(:), allocatable :: next_key
    type(string), allocatable :: key_set(:)
    type(dictionary_t) :: subdict

    key_set = get_key_set(dict)
    allocate(section_indices(size(key_set)))
    ! We need to iterate the sections in order they were given.
    ! Order of keys in dict can be arbitrary, so we need to use dict index stored from parsing.
    do i = 1, size(key_set)
      call associate(subdict, dict, key_set(i) % to_char_str())
      call assign(section_indices(i), subdict, 'dict_index')
      call nullify(subdict, 'dict_index')
    end do
    section_indices = bubble_sort(section_indices, index_permutation = key_order)

    allocate(wf_sections(size(key_set)))
    do i = 1, size(key_order)
      next_key = key_set(key_order(i)) % to_char_str()
      call associate(subdict, dict, next_key)
      call wf_sections(i) % checked_init(subdict)
      if (wf_sections(i) % name == 'use_key_name') then
        wf_sections(i) % name = next_key
      end if
    end do
  end function

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
        K(1) = str2int(K_str)
        K(2) = K(1)
      else
        ! K range
        tokens = strsplit(K_str, '..')
        K(1) = str2int(tokens(1) % to_char_str())
        K(2) = str2int(tokens(2) % to_char_str())
      end if
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of input_params from a given *config_dict* with user set key-value parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_dict_input_params(this, config_dict)
    class(input_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    logical :: wf_sections_provided
    integer :: i
    character(:), allocatable :: next_key, K_str
    type(string), allocatable :: key_set(:)
    type(dictionary_t) :: subdict

    wf_sections_provided = .false.
    ! Iterate over the keys given by user
    key_set = get_key_set(config_dict)
    do i = 1, size(key_set)
      next_key = key_set(i) % to_char_str()
      select case (next_key)
        case ('stage')
          this % stage = extract_string(config_dict, next_key)
        case ('use_rovib_coupling')
          this % use_rovib_coupling = str2int(extract_string(config_dict, next_key))
        case ('use_fixed_basis_JK')
          this % use_fixed_basis_JK = str2int(extract_string(config_dict, next_key))
        case ('cap_type')
          this % cap_type = extract_string(config_dict, next_key)
        case ('grid_rho')
          call associate(subdict, config_dict, next_key)
          call this % grid_rho % checked_init(subdict)
        case ('grid_theta')
          call associate(subdict, config_dict, next_key)
          call this % grid_theta % checked_init(subdict)
          call convert_grid_params_deg_to_rad(this % grid_theta)
        case ('num_points_phi')
          this % num_points_phi = str2int(extract_string(config_dict, next_key))
        case ('output_coordinate_system')
          this % output_coordinate_system = extract_string(config_dict, next_key)
        case ('mass')
          this % mass = parse_mass(extract_string(config_dict, next_key))
        case ('J')
          this % J = str2int(extract_string(config_dict, next_key))
        case ('K')
          K_str = extract_string(config_dict, next_key)
        case ('parity')
          this % parity = str2int(extract_string(config_dict, next_key))
        case ('symmetry')
          this % symmetry = str2int(extract_string(config_dict, next_key))
        case ('basis_size_phi')
          this % basis_size_phi = str2int(extract_string(config_dict, next_key))
        case ('cutoff_energy')
          this % cutoff_energy = str2real(extract_string(config_dict, next_key)) / au_to_wn
        case ('basis_root_path')
          this % basis_root_path = extract_string(config_dict, next_key)
        case ('basis_J')
          this % basis_J = str2int(extract_string(config_dict, next_key))
        case ('basis_K')
          this % basis_K = str2int(extract_string(config_dict, next_key))
        case ('num_states')
          this % num_states = str2int(extract_string(config_dict, next_key))
        case ('ncv')
          this % ncv = str2int(extract_string(config_dict, next_key))
        case ('mpd')
          this % mpd = str2int(extract_string(config_dict, next_key))
        case ('max_iterations')
          this % max_iterations = str2int(extract_string(config_dict, next_key))
        case ('wf_sections')
          wf_sections_provided = .true.
          call associate(subdict, config_dict, next_key)
          this % wf_sections = parse_wf_sections(subdict)
        case ('grid_path')
          this % grid_path = extract_string(config_dict, next_key)
        case ('root_path')
          this % root_path = extract_string(config_dict, next_key)
        case ('use_parallel')
          this % use_parallel = str2int(extract_string(config_dict, next_key))
        case ('enable_terms')
          this % enable_terms = parse_enable_terms(extract_string(config_dict, next_key))
        case ('optimized_mult')
          this % optimized_mult = str2int(extract_string(config_dict, next_key))
        case ('debug_mode')
          this % debug_mode = extract_string(config_dict, next_key)
        case ('debug_param_1')
          this % debug_param_1 = extract_string(config_dict, next_key)
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
      if (this % basis_size_phi == -1) then
        call put_string(keys, 'num_points_phi')
      end if
    end if

    if (this % stage == 'grids' .and. (this % grid_rho % optimized == 1 .or. any(this % output_coordinate_system == [character(100) :: 'jacobi', 'cartesian', 'all bonds', 'internal'])) .or. &
        this % stage == 'basis' .or. this % stage == 'overlaps' .or. this % stage == 'eigencalc' .or. this % stage == 'properties') then
      call put_string(keys, 'mass')
    end if

    if (this % stage == 'basis' .or. this % stage == 'overlaps') then
      call put_string(keys, 'cutoff_energy')
    end if

    if (this % stage == 'basis' .or. this % stage == 'overlaps' .or. this % stage == 'eigencalc' .or. this % stage == 'properties') then
      call put_string(keys, 'symmetry')
      call put_string(keys, 'grid_path')
      call put_string(keys, 'root_path')
    end if

    if (this % stage == 'basis' .or. this % stage == 'overlaps' .or. (this % stage == 'eigencalc' .or. this % stage == 'properties') .and. this % use_rovib_coupling == 0) then
      call put_string(keys, 'K')
    end if

    if (this % stage == 'grids' .and. this % basis_size_phi /= -1 .and. this % num_points_phi == -1 .or. &
        this % stage == 'basis' .or. this % stage == 'overlaps' .or. this % stage == 'properties') then
      call put_string(keys, 'basis_size_phi')
    end if

    if (this % stage == 'basis' .or. this % stage == 'eigencalc' .or. this % stage == 'properties') then
      call put_string(keys, 'J')
    end if

    if (this % stage == 'overlaps' .or. this % stage == 'eigencalc' .or. this % stage == 'properties') then
      call put_string(keys, 'use_rovib_coupling')
      call put_string(keys, 'use_fixed_basis_JK')
    end if

    if (this % stage == 'eigencalc' .or. this % stage == 'properties') then
      call put_string(keys, 'num_states')
    end if

    if ((this % stage == 'eigencalc' .or. this % stage == 'properties') .and. this % use_rovib_coupling == 1) then
      call put_string(keys, 'parity')
    end if

    if ((this % stage == 'eigencalc' .or. this % stage == 'properties') .and. this % use_fixed_basis_JK == 1) then
      call put_string(keys, 'basis_root_path')
      call put_string(keys, 'basis_J')
      call put_string(keys, 'basis_K')
    end if

    if (this % stage == 'properties') then
      call put_string(keys, 'wf_sections')
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of optional keys and custom messages for default announcements.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_optional_keys_input_params(this, keys, messages)
    class(input_params), intent(in) :: this
    type(dictionary_t), intent(out) :: keys, messages

    if (this % stage == 'grids') then
      call put_string(keys, 'num_points_phi', num2str(2 * this % basis_size_phi))
      call put_string(messages, 'num_points_phi', 'Assuming default = 2 * basis_size_phi')
      call put_string(keys, 'output_coordinate_system', 'aph')
    end if

    if (this % stage == 'basis' .or. this % stage == 'overlaps' .or. this % stage == 'eigencalc' .or. this % stage == 'properties') then
      call put_string(keys, 'use_parallel', '1')
      call put_string(messages, 'use_parallel', '')
    end if

    if (this % stage == 'eigencalc') then
      call put_string(keys, 'ncv', '-1')
      call put_string(messages, 'ncv', 'Its value will be determined by SLEPc')
      call put_string(keys, 'mpd', '-1')
      call put_string(messages, 'mpd', 'Its value will be determined by SLEPc')
      call put_string(keys, 'max_iterations', '10000')
    end if

    if (this % stage == 'eigencalc' .and. this % use_rovib_coupling == 1) then
      call put_string(keys, 'enable_terms', '11')
      call put_string(messages, 'enable_terms', '')
      call put_string(keys, 'optimized_mult', '1')
      call put_string(messages, 'optimized_mult', '')
    end if

    if (this % stage == 'eigencalc' .or. this % stage == 'properties') then
      call put_string(keys, 'cap_type', 'none')
    end if

    if ((this % stage == 'eigencalc' .or. this % stage == 'properties') .and. this % use_rovib_coupling == 1) then
      call put_string(keys, 'K', 'all')
    end if
end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a set of all known keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_all_keys_input_params(this) result(keys)
    class(input_params), intent(in) :: this
    type(dictionary_t) :: keys

    call put_string(keys, 'stage')
    call put_string(keys, 'use_rovib_coupling')
    call put_string(keys, 'use_fixed_basis_JK')
    call put_string(keys, 'cap_type')
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
    call put_string(keys, 'basis_root_path')
    call put_string(keys, 'basis_J')
    call put_string(keys, 'basis_K')
    call put_string(keys, 'num_states')
    call put_string(keys, 'ncv')
    call put_string(keys, 'mpd')
    call put_string(keys, 'max_iterations')
    call put_string(keys, 'wf_sections')
    call put_string(keys, 'grid_path')
    call put_string(keys, 'root_path')
    call put_string(keys, 'use_parallel')
    call put_string(keys, 'enable_terms')
    call put_string(keys, 'optimized_mult')
    call put_string(keys, 'debug_mode')
    call put_string(keys, 'debug_param_1')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of provided values. -1 means the value was not provided.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_values_input_params(this) 
    class(input_params), intent(in) :: this
    call assert(any(this % stage == [character(100) :: 'grids', 'basis', 'overlaps', 'eigencalc', 'properties']), 'Error: stage can be "grids", "basis", "overlaps", "eigencalc" or "properties"')
    call assert(any(this % use_rovib_coupling == [-1, 0, 1]), 'Error: use_rovib_coupling can be 0 or 1')
    call assert(any(this % use_fixed_basis_JK == [-1, 0, 1]), 'Error: use_fixed_basis_JK can be 0 or 1')
    if (allocated(this % cap_type)) then
      call assert(any(this % cap_type == [character(100) :: 'none', 'Manolopoulos']), 'Error: cap_type can be "none" or "Manolopoulos"')
    end if
    call assert((this % grid_theta % from .aeq. -1d0) .or. this % grid_theta % from < pi, 'Error: grid_theta % from should be < pi')
    call assert((this % grid_theta % to .aeq. -1d0) .or. this % grid_theta % to < pi, 'Error: grid_theta % to should be < pi')
    call assert(this % num_points_phi == -1 .or. this % num_points_phi > 1, 'Error: num_points_phi should be > 1')
    call assert(.not. allocated(this % output_coordinate_system) .or. &
        any(this % output_coordinate_system == [character(100) :: 'aph', 'mass jacobi', 'jacobi', 'cartesian', 'all bonds', 'internal']), &
        'Error: output_coordinate_system can be "aph", "mass jacobi", "jacobi", "cartesian", "all bonds" or "internal"')
    call assert(all(this % mass .aeq. -1d0) .or. all(this % mass > 0), 'Error: all mass should be > 0')
    call assert(this % J >= -1, 'Error: J should be >= 0')
    call assert(any(this % parity == [-1, 0, 1]), 'Error: parity value can be 0 or 1')
    if (any(this % stage == [character(len = 100) :: 'eigencalc', 'properties']) .and. this % use_rovib_coupling == 1) then
      call assert(all(this % K == -1) .or. all(this % K >= get_k_start(this % J, this % parity)), 'Error: K should be >= mod(J+p, 2)')
    else
      call assert(all(this % K >= -1), 'Error: K should be >= 0')
    end if
    call assert(all(this % K <= this % J), 'Error: K should be <= J')
    call assert(any(this % symmetry == [-1, 0, 1]), 'Error: symmery value can be 0 or 1')
    call assert(this % basis_size_phi == -1 .or. this % basis_size_phi > 0, 'Error: basis_size_phi should be > 0')
    call assert(this % basis_J >= -1, 'Error: basis_J should be >= 0')
    call assert(this % basis_K >= -1, 'Error: basis_K should be >= 0')
    call assert(this % basis_K <= this % basis_J, 'Error: basis_K should be <= basis_J')
    call assert(this % num_states == -1 .or. this % num_states > 0, 'Error: num_states should be > 0')
    call assert(this % ncv == -1 .or. this % ncv > this % num_states, 'Error: ncv should be > num_states')
    call assert(this % mpd == -1 .or. this % mpd > 1, 'Error: mpd should be > 1')
    call assert(this % max_iterations == -1 .or. this % max_iterations > 0, 'Error: max_iterations should be > 0')
    call assert(any(this % use_parallel == [-1, 0, 1]), 'Error: use_parallel can be 0 or 1')
    call assert(any(this % enable_terms(1) == [-1, 0, 1]), 'Error: enable_terms(1) can be 0 or 1')
    call assert(any(this % enable_terms(2) == [-1, 0, 1]), 'Error: enable_terms(2) can be 0 or 1')
    call assert(any(this % optimized_mult == [-1, 0, 1]), 'Error: optimized_mult can be 0 or 1')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes an instance of input_params from a given *config_dict* with user set key-value parameters.
! Validates created instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine checked_init_input_params(this, config_dict)
    class(input_params), intent(inout) :: this
    class(dictionary_t) :: config_dict ! intent(in)
    type(dictionary_t) :: mandatory_keys, optional_keys, messages, optional_nonset_keys, all_keys

    call this % assign_dict(config_dict)
    mandatory_keys = this % get_mandatory_keys()
    call check_mandatory_keys(config_dict, mandatory_keys)

    all_keys = this % get_all_keys()
    call check_extra_keys(config_dict, all_keys)
    call this % get_optional_keys(optional_keys, messages)
    call check_unused_keys(config_dict, mandatory_keys, optional_keys)

    if (len(optional_keys) > 0) then
      optional_nonset_keys = set_difference(optional_keys, config_dict)
      if (len(optional_nonset_keys) > 0) then
        call this % assign_dict(optional_nonset_keys)
        call announce_defaults(optional_nonset_keys, messages)
      end if
    end if
    call this % check_values()
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Resolves parameters that require grid information and performs grid related checks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_resolve_grids_input_params(this, grid_rho, grid_theta, grid_phi)
    class(input_params), intent(inout) :: this
    real(real64), intent(in) :: grid_rho(:), grid_theta(:), grid_phi(:)

    call assert(this % basis_size_phi == -1 .or. this % basis_size_phi <= size(grid_phi) / 2, 'Error: basis_size_phi should be <= size(grid_phi) / 2')
    if (this % stage == 'properties') then
      call this % wf_sections % checked_resolve_rho_bounds(grid_rho(1), grid_rho(size(grid_rho)))
      call this % wf_sections % checked_resolve_theta_bounds(grid_theta(1), grid_theta(size(grid_theta)))
      call this % wf_sections % checked_resolve_phi_bounds(0d0, 2*pi)
    end if
  end subroutine

end module
