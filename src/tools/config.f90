module config
  use constants
  use dict_utils
  use dictionary
  use ftlRegexModule
  use general_utils
  use input_params_mod
  use io_utils
  use parallel_utils
  use rovib_utils_mod
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
      if (iostat == -1) then
        exit ! EOF
      end if 
      
      line_tokens = strsplit(line, '!') ! separate inline comments
      line = trim(line_tokens(1) % to_char_str())
      if (line == '') then
        cycle ! skip empty/commented lines
      end if
      
      line_tokens = strsplit(line, '=') ! separate key value pair
      if (size(line_tokens) /= 2) then
        print '(A)', 'An error in config format on line ' // num2str(line_num) // '. Only key-value assignments and comments are allowed.'
        stop 
      end if
      
      key = trim(adjustl(line_tokens(1) % to_char_str()))
      value = trim(adjustl(line_tokens(2) % to_char_str()))
      call put_string(config, key, value)
    end do
    close(file_unit)
  end function
  
! This function can be uncommented if there is ever a need to introduce parameter groups, i.e. a set of settings such that the settings 
! within the group are either all set or all unset
!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks if all the *keys* are simultaneously present or absent in *config_dict*
!-------------------------------------------------------------------------------------------------------------------------------------------
  ! subroutine check_setting_group(config_dict, keys)
  !   class(dictionary_t) :: config_dict ! intent(in)
  !   class(string), intent(in) :: keys(:)
  !   integer :: i
  !   integer :: status(size(keys))
  !   character(:), allocatable :: key
  !
  !   status = 0
  !   do i = 1, size(keys)
  !     key = keys(i) % to_char_str()
  !     if (key .in. config_dict) then
  !       status(i) = 1
  !     end if
  !   end do
  !
  !   if (any(status == 1) .and. .not. all(status == 1)) then
  !     call print_parallel('Error: the following keys must all be set if any of them is set: ' // string_arr_to_char_str(keys))
  !     stop
  !   end if
  ! end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks that all setting groups are set/unset
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_setting_groups(config_dict)
    class(dictionary_t) :: config_dict ! intent(in)
    ! call check_setting_group(config_dict, string(['basis_root_path', 'basis_J', 'basis_K']))
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Fills out 2 dictionaries. A key is a name of a config parameter. The corresponding value is a regexp of all mode ids where this parameter
! is mandatory (1st dict) or optional (2nd dict). See get_mode_id for description of id string
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine populate_mode_settings(mandatory_modes, optional_modes)
    class(dictionary_t), intent(out) :: mandatory_modes, optional_modes

    call put_string(mandatory_modes, 'mode', '.*')
    call put_string(mandatory_modes, 'rovib_coupling', '^[2-4].*')
    call put_string(mandatory_modes, 'fix_basis_jk', '^[2-4].*')

    call put_string(mandatory_modes, 'molecule', '^([0-4]).*')
    call put_string(mandatory_modes, 'J', '^(1|[3-4]).*')
    call put_string(mandatory_modes, 'K', '^([1-2]|[3-4]0).*')
    call put_string(mandatory_modes, 'parity', '^[3-4]1.*')
    call put_string(mandatory_modes, 'symmetry', '^[1-4].*')
    call put_string(mandatory_modes, 'basis_size_phi', '^([1-2|4]).*')
    call put_string(mandatory_modes, 'cutoff_energy', '^[1-2].*')
    call put_string(mandatory_modes, 'basis_root_path', '^[3-4]11.*')
    call put_string(mandatory_modes, 'basis_J', '^[3-4]11.*')
    call put_string(mandatory_modes, 'basis_K', '^[3-4]11.*')
    call put_string(mandatory_modes, 'num_states', '^[3-4].*')
    call put_string(mandatory_modes, 'grid_path', '^[1-4].*')
    call put_string(mandatory_modes, 'root_path', '^[1-4].*')
    call put_string(mandatory_modes, 'channels_root', '^4.*')

    call put_string(optional_modes, 'K', '^[3-4]1.*')
    call put_string(optional_modes, 'solver', '^3.*')
    call put_string(optional_modes, 'ncv', '^3.*')
    call put_string(optional_modes, 'mpd', '^3.*')
    call put_string(optional_modes, 'max_iterations', '^3.*')
    call put_string(optional_modes, 'cap_type', '^(3|4).*')
    call put_string(optional_modes, 'print_potential', '^0.*')
    call put_string(optional_modes, 'sequential', '^[1-4].*')
    call put_string(optional_modes, 'enable_terms', '^31.*')
    call put_string(optional_modes, 'optimized_mult', '^31.*')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Builds a string that incorporates settings relevant for mode determination
! 1st digit - mode (basis, overlaps, diagonalization, etc.)
! 2nd digit - whether rovib coupling is enabled
! 3rd digit - whether fixed basis is enabled
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mode_id(config_dict) result(mode_id)
    class(dictionary_t) :: config_dict ! intent(in)
    character(3) :: mode_id
    character(:), allocatable :: mode, rovib, fix_basis_jk

    mode_id = repeat('x', len(mode_id))
    mode = item_or_default(config_dict, 'mode', 'unset')
    call assert(mode /= 'unset', 'Error: mode has to be specified')
    if (mode == 'pesprint') then
      mode_id(1:1) = '0'
    else if (mode == 'basis') then
      mode_id(1:1) = '1'
    else if (mode == 'overlaps') then
      mode_id(1:1) = '2'
    else if (mode == 'diagonalization') then
      mode_id(1:1) = '3'
    else if (mode == 'properties') then
      mode_id(1:1) = '4'
    else
      call print_parallel('Error: unknown mode')
      stop
    end if

    mode_id(2:2) = '0'
    rovib = item_or_default(config_dict, 'rovib_coupling', 'unset')
    if (mode == 'overlaps' .or. mode == 'diagonalization' .or. mode == 'properties') then
      call assert(rovib /= 'unset', 'Error: rovib_coupling has to be specified in this mode')
      call assert(rovib == '0' .or. rovib == '1', 'Error: rovib_coupling can be 0 or 1')
      mode_id(2:2) = rovib(1:1)
    else
      if (rovib /= 'unset') then
        call print_parallel('Warning: rovib_coupling parameter has no effect in this mode')
      end if
    end if

    mode_id(3:3) = '0'
    fix_basis_jk = item_or_default(config_dict, 'fix_basis_jk', 'unset')
    if (mode == 'overlaps' .or. mode == 'diagonalization' .or. mode == 'properties') then
      call assert(fix_basis_jk /= 'unset', 'Error: fix_basis_jk has to be specified in this mode')
      call assert(fix_basis_jk == '0' .or. fix_basis_jk == '1', 'Error: fix_basis_jk can be 0 or 1')
      mode_id(3:3) = fix_basis_jk(1:1)
    else
      if (fix_basis_jk /= 'unset') then
        call print_parallel('Warning: fix_basis_jk parameter has no effect in this mode')
      end if
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks whether all mandatory keys of this mode are specified in config
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_mandatory_keys(mode_id, config_dict, mandatory_modes)
    character(*), intent(in) :: mode_id
    class(dictionary_t) :: config_dict, mandatory_modes ! intent(in)
    integer :: i
    character(:), allocatable :: key, mode_pattern
    type(string), allocatable :: keys(:)
    type(ftlRegex) :: re
    type(ftlRegexMatch), allocatable :: m(:)

    keys = key_set(mandatory_modes)
    do i = 1, size(keys)
      key = keys(i) % to_char_str()
      mode_pattern = extract_string(mandatory_modes, key)
      call re % New(mode_pattern)
      m = re % Match(mode_id)
      if (size(m) > 0 .and. .not. (key .in. config_dict)) then
        call print_parallel('Error: ' // key // ' value must be set in the current mode')
        stop
      end if
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks if some of the specified keys will be unused
!-------------------------------------------------------------------------------------------------------------------------------------------
  function check_mode_match(modes, key, mode_id) result(is_match)
    class(dictionary_t) :: modes ! intent(in)
    character(*), intent(in) :: key, mode_id
    logical :: is_match
    character(:), allocatable :: mode_pattern
    type(ftlRegex) :: re
    type(ftlRegexMatch), allocatable :: matches(:)

    is_match = .false.
    mode_pattern = item_or_default(modes, key, '')
    if (len(mode_pattern) > 0) then
      call re % New(mode_pattern)
      matches = re % Match(mode_id)
      is_match = size(matches) > 0
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks whether some of the specified keys will be unused
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_extra_keys(mode_id, config_dict, mandatory_modes, optional_modes)
    character(*), intent(in) :: mode_id
    class(dictionary_t) :: config_dict, mandatory_modes, optional_modes ! intent(in)
    integer :: i
    character(:), allocatable :: key
    type(string), allocatable :: keys(:)

    keys = key_set(config_dict)
    do i = 1, size(keys)
      key = keys(i) % to_char_str()
      if (.not. (key .in. mandatory_modes) .and. .not. (key .in. optional_modes)) then
        call print_parallel('Error: "' // key // '" keyword is not recognized')
        stop
      else if (.not. check_mode_match(mandatory_modes, key, mode_id) .and. .not. check_mode_match(optional_modes, key, mode_id)) then
        call print_parallel('Warning: ' // key // ' value is not used in this mode')
      end if
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Extracts info from *config_dict* and fills out *config* structure
!-------------------------------------------------------------------------------------------------------------------------------------------
  function process_raw_config(config_dict) result(config)
    class(dictionary_t) :: config_dict ! intent(in)
    type(input_params) :: config
    character(:), allocatable :: mode, molecule, K_str, basis_root_path, solver, cap_type, grid_path, root_path, channels_root, enable_terms_str, debug_mode, test_mode, &
      debug_param_1
    integer :: rovib, fix_basis_jk, J, parity, symmetry, basis_size_phi, basis_J, basis_K, ncv, mpd, num_states, print_potential, max_iterations, sequential, &
      optimized_mult
    integer :: pos
    integer :: K(2), enable_terms(2)
    real*8 :: cutoff_energy
    type(string), allocatable :: tokens(:)

    ! Extract parameters from dictionary. The default values here are only assigned to unused parameters or parameters with non-constant defaults
    mode = item_or_default(config_dict, 'mode', '-1')
    rovib = str2int(item_or_default(config_dict, 'rovib_coupling', '-1'))
    fix_basis_jk = str2int(item_or_default(config_dict, 'fix_basis_jk', '-1'))

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
        K(1) = str2int(K_str)
        K(2) = K(1)
      else
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

    solver = item_or_default(config_dict, 'solver', '-1')
    num_states = str2int(item_or_default(config_dict, 'num_states', '-1'))
    ncv = str2int(item_or_default(config_dict, 'ncv', '-1'))
    mpd = str2int(item_or_default(config_dict, 'mpd', '-1'))
    max_iterations = str2int(item_or_default(config_dict, 'max_iterations', '-1'))

    cap_type = item_or_default(config_dict, 'cap_type', '-1')

    grid_path = item_or_default(config_dict, 'grid_path', '-1')
    root_path = item_or_default(config_dict, 'root_path', '-1')
    channels_root = item_or_default(config_dict, 'channels_root', '-1')

    print_potential = str2int(item_or_default(config_dict, 'print_potential', '-1'))

    sequential = str2int(item_or_default(config_dict, 'sequential', '-1'))
    enable_terms_str = item_or_default(config_dict, 'enable_terms', '-1')
    if (enable_terms_str == '-1') then
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

    config = input_params(mode, rovib, fix_basis_jk, molecule, J, K, parity, symmetry, basis_size_phi, cutoff_energy, basis_root_path, basis_J, basis_K, solver, &
        num_states, ncv, mpd, max_iterations, cap_type, grid_path, root_path, channels_root, print_potential, sequential, enable_terms, optimized_mult, debug_mode, &
        test_mode, debug_param_1)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Announces default setting of *key*, represented by *default_str*, if it is used in *mode_id*. 
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine announce_default(mode_id, key, default_str, default_comment, optional_modes, custom_message)
    character(*), intent(in) :: mode_id, key, default_str, default_comment
    character(:), allocatable :: mode_pattern
    class(dictionary_t) :: optional_modes ! intent(in)
    character(*), optional, intent(in) :: custom_message
    type(ftlRegex) :: re
    type(ftlRegexMatch), allocatable :: m(:)

    mode_pattern = item_or_default(optional_modes, key, '')
    call assert(mode_pattern /= '', 'Assertion error: attempt to set default for non-existing key')
    call re % New(mode_pattern)
    m = re % Match(mode_id)
    if (size(m) > 0) then
      if (present(custom_message)) then
        call print_parallel(custom_message)
      else
        call print_parallel(key // ' is not specified; assuming default = ' // default_str // ' ' // default_comment)
      end if
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets default values
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_defaults(mode_id, params, optional_modes)
    character(*), intent(in) :: mode_id
    class(input_params), intent(inout) :: params
    class(dictionary_t) :: optional_modes ! intent(in)

    if (params % cap_type == '-1') then
      params % cap_type = 'none'
      call announce_default(mode_id, 'cap_type', params % cap_type, '', optional_modes)
    end if

    params % solver = iff(params % solver == '-1', 'slepc', params % solver)

    if (params % ncv == -1 .and. params % solver == 'slepc') then
      call announce_default(mode_id, 'ncv', '', '', optional_modes, 'ncv is not specified; its value will be determined by SLEPc. This might affect the number of converged ' // &
         'eigenvalues')
    end if
    if (params % mpd == -1 .and. params % solver == 'slepc') then
      call announce_default(mode_id, 'mpd', '', '', optional_modes, 'mpd is not specified; its value will be determined by SLEPc. This might affect the number of converged ' // &
         'eigenvalues')
    end if

    ! Silent defaults. These parameters should not normally be changed.
    params % max_iterations = iff(params % max_iterations == -1, 10000, params % max_iterations)
    params % print_potential = iff(params % print_potential == -1, 0, params % print_potential)
    params % sequential = iff(params % sequential == -1, 0, 1)
    params % enable_terms(1) = iff(params % enable_terms(1) == -1, 1, params % enable_terms(1))
    params % enable_terms(2) = iff(params % enable_terms(2) == -1, 1, params % enable_terms(2))
    params % optimized_mult = iff(params % optimized_mult == -1, 1, params % optimized_mult)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks validity of values of parameters
! A value of -1 means the value was unassigned, which means it's unused in this mode and should not be checked
! Therefore -1 is always considered a valid value (but that is not displayed in the error messages)
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_params_values(params)
    class(input_params), intent(in) :: params
    integer :: K_min

    call assert(any(params % molecule == [character(len = 100) :: '-1', '686', '868']), 'Error: molecule can be "686" or "868"')
    call assert(params % J >= -1, 'Error: J should be >= 0')
    call assert(any(params % parity == [-1, 0, 1]), 'Error: parity value can be 0 or 1')
    if (any(params % mode == [character(len = 100) :: 'diagonalization', 'properties']) .and. params % rovib_coupling == 1) then
      K_min = get_k_start(params % J, params % parity)
      call assert(all(params % K == -1) .or. all(params % K >= K_min), 'Error: K(1:2) should be >= mod(J+p, 2)')
    else
      call assert(all(params % K >= -1), 'Error: K(1:2) should be >= 0')
    end if
    ! call assert(all(params % K <= params % J), 'Error: K(1:2) should be <= J')
    call assert(any(params % symmetry == [-1, 0, 1]), 'Error: symmery value can be 0 or 1')

    call assert(params % basis_size_phi == -1 .or. params % basis_size_phi > 0, 'Error: basis_size_phi should be > 0')
    call assert(params % basis_J >= -1, 'Error: basis_J should be >= 0')
    call assert(params % basis_K >= -1, 'Error: basis_K should be >= 0')
    ! call assert(params % basis_K <= params % basis_J, 'Error: basis_K should be <= basis_J')

    call assert(any(params % solver == [character(len = 100) :: 'slepc']), 'Error: solver can be "slepc"')
    call assert(params % num_states == -1 .or. params % num_states > 0, 'Error: num_states should be > 0')
    call assert(params % ncv == -1 .or. params % ncv > 0, 'Error: ncv should be > 0')
    call assert(params % mpd == -1 .or. params % mpd > 0, 'Error: mpd should be > 0')
    call assert(params % max_iterations == -1 .or. params % max_iterations > 0, 'Error: max_iterations should be > 0')

    call assert(any(params % cap_type == [character(len = 100) :: 'none', 'Manolopoulos']), 'Error: cap_type can be "none" or "Manolopoulos"')
    call assert(any(params % print_potential == [-1, 0, 1]), 'Error: print_potential can be 0 or 1')
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
    character(:), allocatable :: sequential
    type(dictionary_t) :: raw_config, mandatory_modes, optional_modes
    type(input_params) :: params
    integer :: ierr
    character(:), allocatable :: mode_id
    
    ! Parse the raw file and structure its content into a dictionary
    raw_config = read_config_dict(settings_path)
    ! Pull out sequential ahead of time as we need to decide if parallel mode should be initialized first
    sequential = item_or_default(raw_config, 'sequential', '0') ! Sequential execution is disabled by default
    if (sequential == '0') then
      call MPI_Init(ierr)
    end if

    call check_setting_groups(raw_config)
    call populate_mode_settings(mandatory_modes, optional_modes)
    mode_id = get_mode_id(raw_config)
    call check_mandatory_keys(mode_id, raw_config, mandatory_modes)
    call check_extra_keys(mode_id, raw_config, mandatory_modes, optional_modes)
    params = process_raw_config(raw_config)
    call set_defaults(mode_id, params, optional_modes)
    call check_params_values(params)
  end function
end module
