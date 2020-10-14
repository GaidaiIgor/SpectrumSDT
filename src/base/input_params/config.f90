module config_mod
  use constants
  use dict_utils
  use dictionary
  use general_utils
  use input_params_mod
  use io_utils
  use iso_fortran_env, only: real64
  use parallel_utils
  use rovib_utils_mod
  use stage_metadata_mod
  use string_mod
  use string_utils
  implicit none

  private
  public :: process_user_settings
  
contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Reads config file (a set of key-value pairs)
!-------------------------------------------------------------------------------------------------------------------------------------------
  function read_config_dict(config_path) result(config)
    character(*), intent(in) :: config_path
    type(dictionary_t) :: config
    integer :: file_unit, iostat, line_num
    character(:), allocatable :: line, key, value
    type(string), allocatable :: line_tokens(:)
    
    open(newunit = file_unit, file = config_path)
    line_num = 0
    do
      line_num = line_num + 1
      line = read_line(file_unit, iostat)
      if (is_iostat_end(iostat)) then
        exit ! EOF
      end if 
      
      line_tokens = strsplit(line, '!') ! separate inline comments
      line = trim(line_tokens(1) % to_char_str())
      if (line == '') then
        cycle ! skip empty/commented lines
      end if
      
      line_tokens = strsplit(line, '=') ! separate key value pair
      call assert(size(line_tokens) == 2, 'Config error on line ' // num2str(line_num) // '. Only key-value assignments and comments are allowed.')
      
      key = trim(adjustl(line_tokens(1) % to_char_str()))
      value = trim(adjustl(line_tokens(2) % to_char_str()))
      call put_string(config, key, value)
    end do
    close(file_unit)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! If *key* is not present in *config_dict*, it adds *key*/*value* pair to both *config_dict* and *set_keys*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_default(key, value, config_dict, set_keys)
    character(*), intent(in) :: key, value
    class(dictionary_t), intent(inout) :: config_dict, set_keys

    if (key .in. config_dict) then
      return
    end if
    call add_if_absent(config_dict, key, value)
    call add_if_absent(set_keys, key, value)
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets default values in *config_dict* and returns actually set keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function set_defaults(config_dict) result(set_keys)
    class(dictionary_t), intent(inout) :: config_dict
    type(dictionary_t) :: set_keys
    character(:), allocatable :: stage, sequential

    call set_default('use_optimized_grid_rho', '0', config_dict, set_keys)
    call set_default('optimized_grid_solver_steps', '2048', config_dict, set_keys)
    call set_default('cap_type', 'none', config_dict, set_keys)
    call set_default('ncv', '-1', config_dict, set_keys)
    call set_default('mpd', '-1', config_dict, set_keys)

    ! These parameters should not normally be changed and their defaults are not announced
    call set_default('max_iterations', '10000', config_dict, set_keys)
    call set_default('enable_terms', '11', config_dict, set_keys)
    call set_default('optimized_mult', '1', config_dict, set_keys)

    ! Conditional defaults
    stage = item_or_default(config_dict, 'stage', 'unset')
    sequential = iff(stage /= 'grids', '0', '1')
    call set_default('sequential', sequential, config_dict, set_keys)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Makes sure *key* is present in *config_dict* and the corresponding value is set to one of the allowed values
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_value_one_of(config_dict, key, allowed_values)
    class(dictionary_t) :: config_dict ! intent(in)
    character(*), intent(in) :: key
    class(string), intent(in) :: allowed_values(:)
    character(:), allocatable :: value

    value = item_or_default(config_dict, key, 'unset')
    call assert(value /= 'unset', 'Error: the following key has to be specified: ' // key)
    call assert(any(allowed_values == value), 'Error: allowed values of "' // key // '" are: ' // string_arr_to_char_str(allowed_values))
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks whether behavior control variables have legit values. This has to be done ahead of normal checking since these values alter
! checking behavior itself.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_behavior_control_values(config_dict)
    class(dictionary_t) :: config_dict ! intent(in)
    character(:), allocatable :: stage

    call check_value_one_of(config_dict, 'stage', string([character(100) :: 'grids', 'basis', 'overlaps', 'eigencalc', 'properties']))
    stage = item_or_default(config_dict, 'stage', 'unset')

    if (stage == 'grids') then
      call check_value_one_of(config_dict, 'use_optimized_grid_rho', string(['0', '1']))
    end if

    if (stage == 'overlaps') then
      call check_value_one_of(config_dict, 'use_rovib_coupling', string(['0', '1']))
      call check_value_one_of(config_dict, 'use_fix_basis_jk', string(['0', '1']))
    end if

    if (stage == 'eigencalc' .or. stage == 'properties') then
      call check_value_one_of(config_dict, 'use_rovib_coupling', string(['0', '1']))
      call check_value_one_of(config_dict, 'use_fix_basis_jk', string(['0', '1']))
      call check_value_one_of(config_dict, 'cap_type', string([character(100) :: 'none', 'Manolopoulos']))
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Makes sure only one of the parameters in the group is set
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_only_one_set(config_dict, group)
    class(dictionary_t) :: config_dict ! intent(in)
    class(string), intent(in) :: group(:)
    logical :: is_set
    integer :: i
    character(:), allocatable :: next_key

    is_set = .false.
    do i = 1, size(group)
      next_key = group(i) % to_char_str()
      if (next_key .in. config_dict) then
        call assert(.not. is_set, 'Only one of the following keys can be set: ' // string_arr_to_char_str(group))
        is_set = .true.
      end if
    end do
    call assert(is_set, 'One of the following keys has be to specified: ' // string_arr_to_char_str(group))
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks group conditions
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_parameter_groups(config_dict)
    class(dictionary_t) :: config_dict ! intent(in)
    character(:), allocatable :: stage

    stage = item_or_default(config_dict, 'stage', 'unset')
    if (stage == 'grids') then
      call check_only_one_set(config_dict, string([character(100) :: 'grid_rho_npoints', 'grid_rho_step']))
      call check_only_one_set(config_dict, string([character(100) :: 'grid_theta_npoints', 'grid_theta_step']))
      call check_only_one_set(config_dict, string([character(100) :: 'grid_phi_npoints', 'grid_phi_step']))
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks whether all mandatory keys of this stage are specified in config
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_mandatory_keys(config_dict, mandatory_keys)
    class(dictionary_t) :: config_dict, mandatory_keys ! intent(in)
    integer :: i
    character(:), allocatable :: next_key
    type(string), allocatable :: mandatory_keys_plain(:)

    mandatory_keys_plain = key_set(mandatory_keys)
    do i = 1, size(mandatory_keys_plain)
      next_key = mandatory_keys_plain(i) % to_char_str()
      call assert(next_key .in. config_dict, 'Error: the following key has to be specified: ' // next_key)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks whether some of the specified keys will be unused
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_extra_keys(config_dict, all_keys)
    class(dictionary_t) :: config_dict, all_keys ! intent(in)
    integer :: i
    character(:), allocatable :: next_key
    type(string), allocatable :: config_keys(:)

    config_keys = key_set(config_dict)
    do i = 1, size(config_keys)
      next_key = config_keys(i) % to_char_str()
      call assert(next_key .in. all_keys, 'Error: the following key is not recognized: ' // next_key)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks whether some of the specified keys will be unused
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_unused_keys(config_dict, mandatory_keys, optional_keys, set_default_keys)
    class(dictionary_t) :: config_dict, mandatory_keys, optional_keys, set_default_keys ! intent(in)
    integer :: i
    character(:), allocatable :: next_key
    type(string), allocatable :: config_keys(:)

    config_keys = key_set(config_dict)
    do i = 1, size(config_keys)
      next_key = config_keys(i) % to_char_str()
      if (.not. ((next_key .in. mandatory_keys) .or. (next_key .in. optional_keys) .or. (next_key .in. set_default_keys))) then
        call print_parallel('Info: the following key is not used: ' // next_key)
      end if
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses user-provded value of mass
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_mass(mass_str) result(mass)
    character(*), intent(in) :: mass_str
    real(real64) :: mass

    if (mass_str == 'o16') then
      mass = oxygen_masses(1)
    else if (mass_str == 'o17') then
      mass = oxygen_masses(2)
    else if (mass_str == 'o18') then
      mass = oxygen_masses(3)
    else
      ! Plain value
      mass = str2real(mass_str) * amu_to_aum
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses user-provided value of `K`
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
! Parses user-provided value of `enable_terms`
!-------------------------------------------------------------------------------------------------------------------------------------------
  function parse_enable_terms(enable_terms_str) result(enable_terms)
    character(*), intent(in) :: enable_terms_str
    integer :: enable_terms(2)
    
    if (enable_terms_str == '-1') then
      ! Mark as not set
      enable_terms(1) = -1
      enable_terms(2) = -1
    else
      enable_terms(1) = str2int(enable_terms_str(1:1))
      enable_terms(2) = str2int(enable_terms_str(2:2))
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Interprets the values set by user in *config_dict* and fills out *config* structure. This procedure does not make assumptions about
! default values of not specified optional arguments.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function process_raw_config(config_dict) result(config)
    class(dictionary_t) :: config_dict ! intent(in)
    type(input_params) :: config
    integer :: use_optimized_grid_rho, use_rovib_coupling, use_fix_basis_jk, grid_rho_npoints, optimized_grid_solver_steps, grid_theta_npoints, &
        grid_phi_npoints, J, parity, symmetry, basis_size_phi, basis_J, basis_K, ncv, mpd, num_states, max_iterations, sequential, optimized_mult
    integer :: K(2), enable_terms(2)
    real(real64) :: grid_rho_from, grid_rho_to, grid_rho_step, envelope_rho_max_energy, grid_theta_from, grid_theta_to, grid_theta_step, grid_phi_from, grid_phi_to, grid_phi_step, &
        mass_central, mass_terminal1, mass_terminal2, cutoff_energy
    character(:), allocatable :: stage, envelope_rho_path, basis_root_path, cap_type, grid_path, root_path, channels_root, debug_mode, test_mode, debug_param_1

    ! Extract parameters from dictionary. The default values here are only assigned to unused parameters or parameters with non-constant defaults
    stage = item_or_default(config_dict, 'stage', '-1')
    use_optimized_grid_rho = str2int(item_or_default(config_dict, 'use_optimized_grid_rho', '-1'))
    use_rovib_coupling = str2int(item_or_default(config_dict, 'use_rovib_coupling', '-1'))
    use_fix_basis_jk = str2int(item_or_default(config_dict, 'use_fix_basis_jk', '-1'))
    cap_type = item_or_default(config_dict, 'cap_type', '-1')

    grid_rho_from = str2real(item_or_default(config_dict, 'grid_rho_from', '-1'))
    grid_rho_to = str2real(item_or_default(config_dict, 'grid_rho_to', '-1'))
    grid_rho_npoints = str2int(item_or_default(config_dict, 'grid_rho_npoints', '-1'))
    grid_rho_step = str2real(item_or_default(config_dict, 'grid_rho_step', '-1'))
    envelope_rho_path = item_or_default(config_dict, 'envelope_rho_path', '-1')
    envelope_rho_max_energy = str2real(item_or_default(config_dict, 'envelope_rho_max_energy', '-1'))
    optimized_grid_solver_steps = str2int(item_or_default(config_dict, 'optimized_grid_solver_steps', '-1'))

    grid_theta_from = str2real(item_or_default(config_dict, 'grid_theta_from', '-1'))
    grid_theta_to = str2real(item_or_default(config_dict, 'grid_theta_to', '-1'))
    grid_theta_npoints = str2int(item_or_default(config_dict, 'grid_theta_npoints', '-1'))
    grid_theta_step = str2real(item_or_default(config_dict, 'grid_theta_step', '-1'))

    grid_phi_from = str2real(item_or_default(config_dict, 'grid_phi_from', '-1'))
    grid_phi_to = str2real(item_or_default(config_dict, 'grid_phi_to', '-1'))
    grid_phi_npoints = str2int(item_or_default(config_dict, 'grid_phi_npoints', '-1'))
    grid_phi_step = str2real(item_or_default(config_dict, 'grid_phi_step', '-1'))

    mass_central = parse_mass(item_or_default(config_dict, 'mass_central', '-1'))
    mass_terminal1 = parse_mass(item_or_default(config_dict, 'mass_terminal1', '-1'))
    mass_terminal2 = parse_mass(item_or_default(config_dict, 'mass_terminal2', '-1'))
    J = str2int(item_or_default(config_dict, 'J', '-1'))
    parity = str2int(item_or_default(config_dict, 'parity', '-1'))
    K = parse_K(item_or_default(config_dict, 'K', '-1'), J, parity)
    symmetry = str2int(item_or_default(config_dict, 'symmetry', '-1'))

    basis_size_phi = str2int(item_or_default(config_dict, 'basis_size_phi', '-1'))
    cutoff_energy = str2real(item_or_default(config_dict, 'cutoff_energy', '-1')) / au_to_wn
    basis_root_path = item_or_default(config_dict, 'basis_root_path', '-1')
    basis_J = str2int(item_or_default(config_dict, 'basis_J', '-1'))
    basis_K = str2int(item_or_default(config_dict, 'basis_K', '-1'))

    num_states = str2int(item_or_default(config_dict, 'num_states', '-1'))
    ncv = str2int(item_or_default(config_dict, 'ncv', '-1'))
    mpd = str2int(item_or_default(config_dict, 'mpd', '-1'))
    max_iterations = str2int(item_or_default(config_dict, 'max_iterations', '-1'))

    grid_path = item_or_default(config_dict, 'grid_path', '-1')
    root_path = item_or_default(config_dict, 'root_path', '-1')
    channels_root = item_or_default(config_dict, 'channels_root', '-1')

    sequential = str2int(item_or_default(config_dict, 'sequential', '-1'))
    enable_terms = parse_enable_terms(item_or_default(config_dict, 'enable_terms', '-1'))
    optimized_mult = str2int(item_or_default(config_dict, 'optimized_mult', '-1'))
    debug_mode = item_or_default(config_dict, 'debug_mode', '-1')
    test_mode = item_or_default(config_dict, 'test_mode', '-1')
    debug_param_1 = item_or_default(config_dict, 'debug_param_1', '-1')

    config = input_params(stage, use_optimized_grid_rho, use_rovib_coupling, use_fix_basis_jk, cap_type, grid_rho_from, grid_rho_to, grid_rho_npoints, grid_rho_step, &
        envelope_rho_path, envelope_rho_max_energy, optimized_grid_solver_steps, grid_theta_from, grid_theta_to, grid_theta_npoints, grid_theta_step, &
        grid_phi_from, grid_phi_to, grid_phi_npoints, grid_phi_step, &
        mass_central, mass_terminal1, mass_terminal2, J, K, parity, symmetry, basis_size_phi, cutoff_energy, basis_root_path, basis_J, basis_K, num_states, ncv, &
        mpd, max_iterations, grid_path, root_path, channels_root, sequential, enable_terms, optimized_mult, debug_mode, test_mode, debug_param_1)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of values of parameters
! A value of -1 means the value was unassigned, which means it is unused in this mode.
! Therefore -1 is always considered a valid value (but that is not displayed in the error messages)
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_params_values(params)
    class(input_params), intent(in) :: params
    integer :: K_min

    call assert((params % grid_rho_from .aeq. -1d0) .or. params % grid_rho_from >= 0, 'Error: grid_rho_from should be >= 0')
    call assert((params % grid_rho_to .aeq. -1d0) .or. params % grid_rho_to > 0, 'Error: grid_rho_to should be > 0')
    if (.not. (params % grid_rho_from .aeq. -1d0) .and. .not. (params % grid_rho_to .aeq. -1d0)) then
      call assert(params % grid_rho_from < params % grid_rho_to, 'Error: grid_rho_from should be < grid_rho_to')
    end if
    call assert(params % grid_rho_npoints == -1 .or. params % grid_rho_npoints > 0, 'Error: grid_rho_npoints should be > 0')
    call assert((params % grid_rho_step .aeq. -1d0) .or. params % grid_rho_step > 0, 'Error: grid_rho_step should be > 0')
    call assert(params % optimized_grid_solver_steps == -1 .or. params % optimized_grid_solver_steps > 0, 'Error: optimized_grid_solver_steps should be > 0')

    call assert((params % grid_theta_from .aeq. -1d0) .or. params % grid_theta_from >= 0, 'Error: grid_theta_from should be >= 0')
    call assert((params % grid_theta_from .aeq. -1d0) .or. params % grid_theta_from < pi, 'Error: grid_theta_from should be < pi')
    call assert((params % grid_theta_to .aeq. -1d0) .or. params % grid_theta_to > 0, 'Error: grid_theta_to should be > 0')
    call assert((params % grid_theta_to .aeq. -1d0) .or. params % grid_theta_to < pi, 'Error: grid_theta_to should be < pi')
    if (.not. (params % grid_theta_from .aeq. -1d0) .and. .not. (params % grid_theta_to .aeq. -1d0)) then
      call assert(params % grid_theta_from < params % grid_theta_to, 'Error: grid_theta_from should be < grid_theta_to')
    end if
    call assert(params % grid_theta_npoints == -1 .or. params % grid_theta_npoints > 0, 'Error: grid_theta_npoints should be > 0')
    call assert((params % grid_theta_step .aeq. -1d0) .or. params % grid_theta_step > 0, 'Error: grid_theta_step should be > 0')

    call assert((params % grid_phi_from .aeq. -1d0) .or. params % grid_phi_from >= 0, 'Error: grid_phi_from should be >= 0')
    call assert((params % grid_phi_from .aeq. -1d0) .or. params % grid_phi_from < 2*pi, 'Error: grid_phi_from should be < 2*pi')
    call assert((params % grid_phi_to .aeq. -1d0) .or. params % grid_phi_to > 0, 'Error: grid_phi_to should be > 0')
    call assert((params % grid_phi_to .aeq. -1d0) .or. (params % grid_phi_to .ale. 2*pi), 'Error: grid_phi_to should be <= 2*pi')
    if (.not. (params % grid_phi_from .aeq. -1d0) .and. .not. (params % grid_phi_to .aeq. -1d0)) then
      call assert(params % grid_phi_from < params % grid_phi_to, 'Error: grid_phi_from should be < grid_phi_to')
    end if
    call assert(params % grid_phi_npoints == -1 .or. params % grid_phi_npoints > 0, 'Error: grid_phi_npoints should be > 0')
    call assert((params % grid_phi_step .aeq. -1d0) .or. params % grid_phi_step > 0, 'Error: grid_phi_step should be > 0')

    call assert((params % mass_central .aeq. -1d0) .or. params % mass_central > 0, 'Error: mass_central should be > 0')
    call assert((params % mass_terminal1 .aeq. -1d0) .or. params % mass_terminal1 > 0, 'Error: mass_terminal1 should be > 0')
    call assert((params % mass_terminal2 .aeq. -1d0) .or. params % mass_terminal2 > 0, 'Error: mass_terminal2 should be > 0')
    call assert(params % J >= -1, 'Error: J should be >= 0')
    call assert(any(params % parity == [-1, 0, 1]), 'Error: parity value can be 0 or 1')
    if (any(params % stage == [character(len = 100) :: 'eigencalc', 'properties']) .and. params % use_rovib_coupling == 1) then
      K_min = get_k_start(params % J, params % parity)
      call assert(all(params % K == -1) .or. all(params % K >= K_min), 'Error: K(1:2) should be >= mod(J+p, 2)')
    else
      call assert(all(params % K >= -1), 'Error: K(1:2) should be >= 0')
    end if
    call assert(all(params % K <= params % J), 'Error: K(1:2) should be <= J')
    call assert(any(params % symmetry == [-1, 0, 1]), 'Error: symmery value can be 0 or 1')

    call assert(params % basis_size_phi == -1 .or. params % basis_size_phi > 0, 'Error: basis_size_phi should be > 0')
    call assert(params % basis_J >= -1, 'Error: basis_J should be >= 0')
    call assert(params % basis_K >= -1, 'Error: basis_K should be >= 0')
    call assert(params % basis_K <= params % basis_J, 'Error: basis_K should be <= basis_J')

    call assert(params % num_states == -1 .or. params % num_states > 0, 'Error: num_states should be > 0')
    call assert(params % ncv == -1 .or. params % ncv > 0, 'Error: ncv should be > 0')
    call assert(params % mpd == -1 .or. params % mpd > 0, 'Error: mpd should be > 0')
    call assert(params % max_iterations == -1 .or. params % max_iterations > 0, 'Error: max_iterations should be > 0')

    call assert(any(params % sequential == [-1, 0, 1]), 'Error: sequential can be 0 or 1')
    call assert(any(params % enable_terms(1) == [-1, 0, 1]), 'Error: enable_terms(1) can be 0 or 1')
    call assert(any(params % enable_terms(2) == [-1, 0, 1]), 'Error: enable_terms(2) can be 0 or 1')
    call assert(any(params % optimized_mult == [-1, 0, 1]), 'Error: optimized_mult can be 0 or 1')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! If *key* was default-set and is actually used with current settings, prints a standard notification message or custom *message*, if given.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine announce_default(key, set_default_keys, optional_keys, message)
    character(*), intent(in) :: key
    class(dictionary_t) :: set_default_keys, optional_keys ! intent(in)
    character(*), optional, intent(in) :: message
    character(:), allocatable :: default_value

    if ((key .in. set_default_keys) .and. (key .in. optional_keys)) then
      if (present(message)) then
        call print_parallel(message)
      else
        default_value = item_or_default(set_default_keys, key, 'unset')
        call print_parallel(key // ' is not specified. Assuming default = ' // default_value)
      end if
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Announces previously set default if its value is going to be used, as specified by *optional_keys*
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine announce_defaults(set_default_keys, optional_keys)
    class(dictionary_t) :: set_default_keys, optional_keys ! intent(in)
    call announce_default('use_optimized_grid_rho', set_default_keys, optional_keys)
    call announce_default('cap_type', set_default_keys, optional_keys)
    call announce_default('ncv', set_default_keys, optional_keys, message = 'ncv is not specified. Its value will be determined by SLEPc')
    call announce_default('mpd', set_default_keys, optional_keys, message = 'mpd is not specified. Its value will be determined by SLEPc')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Assigns default values to some parameters if unset by user. Not all parameters have reasonable defaults.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function process_user_settings(settings_path) result(params)
    character(*), intent(in) :: settings_path ! path to the file with settings
    type(dictionary_t) :: config_dict, set_default_keys, mandatory_keys, optional_keys, all_keys
    type(input_params) :: params
    integer :: ierr
    
    ! Parse the raw file and structure its content into a dictionary
    config_dict = read_config_dict(settings_path)
    set_default_keys = set_defaults(config_dict)
    ! Pull out sequential ahead of time as we need to decide if parallel mode should be initialized first
    if (item_or_default(config_dict, 'sequential', '0') == '0') then
      call MPI_Init(ierr)
    end if

    call check_behavior_control_values(config_dict)
    call check_parameter_groups(config_dict)
    mandatory_keys = get_mandatory_keys(config_dict)
    optional_keys = get_optional_keys(config_dict)
    all_keys = get_all_keys()

    call check_mandatory_keys(config_dict, mandatory_keys)
    call check_extra_keys(config_dict, all_keys)
    call check_unused_keys(config_dict, mandatory_keys, optional_keys, set_default_keys)

    params = process_raw_config(config_dict)
    call check_params_values(params)
    call announce_defaults(set_default_keys, optional_keys)
  end function
end module
