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
! Checks whether behavior control variables have legit values. This has to be done ahead of normal checking since these values alter
! checking behavior.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_behavior_control_values(config_dict)
    class(dictionary_t) :: config_dict ! intent(in)
    character(:), allocatable :: stage, rovib_coupling, fix_basis_jk

    stage = item_or_default(config_dict, 'stage', 'unset')
    rovib_coupling = item_or_default(config_dict, 'rovib_coupling', 'unset')
    fix_basis_jk = item_or_default(config_dict, 'fix_basis_jk', 'unset')

    call assert(stage /= 'unset', 'Error: stage has to be specified')
    call assert(any(stage == [character(100) :: 'grids', 'basis', 'overlaps', 'eigencalc', 'properties']), &
        'Error: allowed values of "stage" are "grids", "basis", "overlaps", "eigencalc" or "properties"')

    if (any(stage == [character(100) :: 'overlaps', 'eigencalc', 'properties'])) then
      call assert(rovib_coupling /= 'unset', 'Error: rovib_coupling has to be specified')
      call assert(fix_basis_jk /= 'unset', 'Error: fix_basis_jk has to be specified')
      call assert(rovib_coupling == '0' .or. rovib_coupling == '1', 'Error: allowed values of "rovib_coupling" are 0 or 1')
      call assert(fix_basis_jk == '0' .or. fix_basis_jk == '1', 'Error: allowed values of "fix_basis_jk" are 0 or 1')
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
    call check_only_one_set(config_dict, to_string_char_str_arr_trim([character(100) :: 'grid_rho_npoints', 'grid_rho_step']))
    call check_only_one_set(config_dict, to_string_char_str_arr_trim([character(100) :: 'grid_theta_npoints', 'grid_theta_step']))
    call check_only_one_set(config_dict, to_string_char_str_arr_trim([character(100) :: 'grid_phi_npoints', 'grid_phi_step']))
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
      call assert(next_key .in. config_dict, 'Error: the following key has to be specified at this stage: ' // next_key)
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
  subroutine check_unused_keys(config_dict, mandatory_keys, optional_keys)
    class(dictionary_t) :: config_dict, mandatory_keys, optional_keys ! intent(in)
    integer :: i
    character(:), allocatable :: next_key
    type(string), allocatable :: config_keys(:)

    config_keys = key_set(config_dict)
    do i = 1, size(config_keys)
      next_key = config_keys(i) % to_char_str()
      if (.not. ((next_key .in. mandatory_keys) .or. (next_key .in. optional_keys))) then
        call print_parallel('Info: the following key is not used at this stage: ' // next_key)
      end if
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Interprets the values set by user in *config_dict* and fills out *config* structure. This procedure does not make assumptions about
! default values of not specified optional arguments.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function process_raw_config(config_dict) result(config)
    class(dictionary_t) :: config_dict ! intent(in)
    type(input_params) :: config
    integer :: rovib, fix_basis_jk, grid_rho_npoints, grid_theta_npoints, grid_phi_npoints, J, parity, symmetry, basis_size_phi, basis_J, basis_K, ncv, mpd, num_states, &
        max_iterations, sequential, optimized_mult
    integer :: pos
    integer :: K(2), enable_terms(2)
    real(real64) :: grid_rho_from, grid_rho_to, grid_rho_step, grid_theta_from, grid_theta_to, grid_theta_step, grid_phi_from, grid_phi_to, grid_phi_step, cutoff_energy
    character(:), allocatable :: stage, molecule, K_str, basis_root_path, cap_type, grid_path, root_path, channels_root, enable_terms_str, debug_mode, test_mode, debug_param_1
    type(string), allocatable :: tokens(:)

    ! Extract parameters from dictionary. The default values here are only assigned to unused parameters or parameters with non-constant defaults
    stage = item_or_default(config_dict, 'stage', '-1')
    rovib = str2int(item_or_default(config_dict, 'rovib_coupling', '-1'))
    fix_basis_jk = str2int(item_or_default(config_dict, 'fix_basis_jk', '-1'))

    grid_rho_from = str2real(item_or_default(config_dict, 'grid_rho_from', '-1'))
    grid_rho_to = str2real(item_or_default(config_dict, 'grid_rho_to', '-1'))
    grid_rho_npoints = str2int(item_or_default(config_dict, 'grid_rho_npoints', '-1'))
    grid_rho_step = str2real(item_or_default(config_dict, 'grid_rho_step', '-1'))

    grid_theta_from = str2real(item_or_default(config_dict, 'grid_theta_from', '-1'))
    grid_theta_to = str2real(item_or_default(config_dict, 'grid_theta_to', '-1'))
    grid_theta_npoints = str2int(item_or_default(config_dict, 'grid_theta_npoints', '-1'))
    grid_theta_step = str2real(item_or_default(config_dict, 'grid_theta_step', '-1'))

    grid_phi_from = str2real(item_or_default(config_dict, 'grid_phi_from', '-1'))
    grid_phi_to = str2real(item_or_default(config_dict, 'grid_phi_to', '-1'))
    grid_phi_npoints = str2int(item_or_default(config_dict, 'grid_phi_npoints', '-1'))
    grid_phi_step = str2real(item_or_default(config_dict, 'grid_phi_step', '-1'))

    molecule = item_or_default(config_dict, 'molecule', '-1')
    J = str2int(item_or_default(config_dict, 'J', '-1'))
    parity = str2int(item_or_default(config_dict, 'parity', '-1'))
    K_str = item_or_default(config_dict, 'K', '-1')
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
    symmetry = str2int(item_or_default(config_dict, 'symmetry', '-1'))

    basis_size_phi = str2int(item_or_default(config_dict, 'basis_size_phi', '-1'))
    cutoff_energy = str2real(item_or_default(config_dict, 'cutoff_energy', '-1')) / autown
    basis_root_path = item_or_default(config_dict, 'basis_root_path', '-1')
    basis_J = str2int(item_or_default(config_dict, 'basis_J', '-1'))
    basis_K = str2int(item_or_default(config_dict, 'basis_K', '-1'))

    num_states = str2int(item_or_default(config_dict, 'num_states', '-1'))
    ncv = str2int(item_or_default(config_dict, 'ncv', '-1'))
    mpd = str2int(item_or_default(config_dict, 'mpd', '-1'))
    max_iterations = str2int(item_or_default(config_dict, 'max_iterations', '-1'))

    cap_type = item_or_default(config_dict, 'cap_type', '-1')

    grid_path = item_or_default(config_dict, 'grid_path', '-1')
    root_path = item_or_default(config_dict, 'root_path', '-1')
    channels_root = item_or_default(config_dict, 'channels_root', '-1')

    sequential = str2int(item_or_default(config_dict, 'sequential', '-1'))
    enable_terms_str = item_or_default(config_dict, 'enable_terms', '-1')
    if (enable_terms_str == '-1') then
      ! Mark as not set
      enable_terms(1) = -1
      enable_terms(2) = -1
    else
      enable_terms(1) = str2int(enable_terms_str(1:1))
      enable_terms(2) = str2int(enable_terms_str(2:2))
    end if
    optimized_mult = str2int(item_or_default(config_dict, 'optimized_mult', '-1'))
    debug_mode = item_or_default(config_dict, 'debug_mode', '-1')
    test_mode = item_or_default(config_dict, 'test_mode', '-1')
    debug_param_1 = item_or_default(config_dict, 'debug_param_1', '-1')

    config = input_params(stage, rovib, fix_basis_jk, grid_rho_from, grid_rho_to, grid_rho_npoints, grid_rho_step, grid_theta_from, grid_theta_to, grid_theta_npoints, &
        grid_theta_step, grid_phi_from, grid_phi_to, grid_phi_npoints, grid_phi_step, molecule, J, K, parity, symmetry, basis_size_phi, cutoff_energy, basis_root_path, basis_J, basis_K, &
        num_states, ncv, mpd, max_iterations, cap_type, grid_path, root_path, channels_root, sequential, enable_terms, optimized_mult, debug_mode, test_mode, debug_param_1)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets default values
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_defaults(params, optional_keys)
    class(input_params), intent(inout) :: params
    class(dictionary_t) :: optional_keys ! intent(in)
    integer :: sequential_default

    if (params % cap_type == '-1' .and. ('cap_type' .in. optional_keys)) then
      params % cap_type = 'none'
      call print_parallel('cap_type is not specified. Assuming default = none')
    end if
    if (params % ncv == -1 .and. ('ncv' .in. optional_keys)) then
      call print_parallel('ncv is not specified. Its value will be determined by SLEPc')
    end if
    if (params % mpd == -1 .and. ('mpd' .in. optional_keys)) then
      call print_parallel('mpd is not specified. Its value will be determined by SLEPc')
    end if

    ! Silent defaults. These parameters should not normally be changed.
    params % max_iterations = iff(params % max_iterations == -1, 10000, params % max_iterations)
    sequential_default = iff(params % stage /= 'grids', 0, 1)
    params % sequential = iff(params % sequential == -1, sequential_default, params % sequential)
    params % enable_terms(1) = iff(params % enable_terms(1) == -1, 1, params % enable_terms(1))
    params % enable_terms(2) = iff(params % enable_terms(2) == -1, 1, params % enable_terms(2))
    params % optimized_mult = iff(params % optimized_mult == -1, 1, params % optimized_mult)
  end subroutine

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

    call assert(any(params % molecule == [character(len = 100) :: '-1', '686', '868']), 'Error: molecule can be "686" or "868"')
    call assert(params % J >= -1, 'Error: J should be >= 0')
    call assert(any(params % parity == [-1, 0, 1]), 'Error: parity value can be 0 or 1')
    if (any(params % stage == [character(len = 100) :: 'eigencalc', 'properties']) .and. params % rovib_coupling == 1) then
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

    call assert(any(params % cap_type == [character(len = 100) :: '-1', 'none', 'Manolopoulos']), 'Error: cap_type can be "none" or "Manolopoulos"')
    call assert(any(params % sequential == [-1, 0, 1]), 'Error: sequential can be 0 or 1')
    call assert(any(params % enable_terms(1) == [-1, 0, 1]), 'Error: enable_terms(1) can be 0 or 1')
    call assert(any(params % enable_terms(2) == [-1, 0, 1]), 'Error: enable_terms(2) can be 0 or 1')
    call assert(any(params % optimized_mult == [-1, 0, 1]), 'Error: optimized_mult can be 0 or 1')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Assigns default values to some parameters if unset by user. Not all parameters have reasonable defaults.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function process_user_settings(settings_path) result(params)
    character(*), intent(in) :: settings_path ! path to the file with settings
    character(:), allocatable :: stage, sequential_default, sequential
    type(dictionary_t) :: raw_config, mandatory_keys, optional_keys, all_keys
    type(input_params) :: params
    integer :: ierr
    
    ! Parse the raw file and structure its content into a dictionary
    raw_config = read_config_dict(settings_path)
    ! Pull out sequential ahead of time as we need to decide if parallel mode should be initialized first
    stage = item_or_default(raw_config, 'stage', 'unset')
    sequential_default = iff(stage /= 'grids', '0', '1')
    sequential = item_or_default(raw_config, 'sequential', sequential_default) ! Sequential execution is disabled by default
    if (sequential == '0') then
      call MPI_Init(ierr)
    end if

    call check_behavior_control_values(raw_config)
    call check_parameter_groups(raw_config)
    mandatory_keys = get_mandatory_keys(raw_config)
    optional_keys = get_optional_keys(raw_config)
    all_keys = get_all_keys()
    call check_mandatory_keys(raw_config, mandatory_keys)
    call check_extra_keys(raw_config, all_keys)
    call check_unused_keys(raw_config, mandatory_keys, optional_keys)

    params = process_raw_config(raw_config)
    call set_defaults(params, optional_keys)
    call check_params_values(params)
  end function
end module
