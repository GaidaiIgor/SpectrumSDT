!-------------------------------------------------------------------------------------------------------------------------------------------
! Common procedures for user input.
!-------------------------------------------------------------------------------------------------------------------------------------------
module config_mod
  use dictionary
  use dict_utils
  use general_utils
  use io_utils
  use iso_fortran_env, only: real64
  use parallel_utils
  use string_mod
  use string_utils
  implicit none
  
  private :: read_inner_dict

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Parses inner dictionary.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive function read_inner_dict(file_unit, recursion_level, line_num, prefix) result(dict)
    integer, intent(in) :: file_unit, recursion_level
    integer, intent(inout) :: line_num
    character(*), intent(in) :: prefix
    type(dictionary_t) :: dict
    integer :: iostat, dict_count
    character(:), allocatable :: line, key, value
    type(string), allocatable :: line_tokens(:)
    type(dictionary_t) :: inner_dict
    
    call put_string(dict, 'prefix', prefix)
    dict_count = 0
    do
      line_num = line_num + 1
      line = read_line(file_unit, iostat)
      if (is_iostat_end(iostat)) then
        call assert(recursion_level == 1, 'Config error: not enough closing brackets.')
        exit
      end if 
      
      line_tokens = strsplit(line, '!') ! separate inline comments
      line = trim(adjustl(line_tokens(1) % to_char_str()))

      if (line == '') then
        ! skip empty/commented lines
        cycle 
      end if
      if (line == '}') then
        ! finished reading inner dict
        call assert(recursion_level > 1, 'Config error at line ' // num2str(line_num) // '. Mismatched bracket.')
        exit 
      end if
      
      line_tokens = strsplit(line, '=') ! separate key value pair
      call assert(size(line_tokens) == 2, 'Config error at line ' // num2str(line_num) // '. Number of assignments is not equal to 1.')
      
      key = trim(adjustl(line_tokens(1) % to_char_str()))
      call assert(.not. (key .in. dict), 'Config error at line ' // num2str(line_num) // '. This key has already been specified.')
      value = trim(adjustl(line_tokens(2) % to_char_str()))
      if (value == '{') then
        ! new inner dict
        inner_dict = read_inner_dict(file_unit, recursion_level + 1, line_num, prefix // key // ' % ')
        dict_count = dict_count + 1
        if (recursion_level > 1) then
          ! remember subdict order in config
          call extend(inner_dict, 'dict_index' .kv. dict_count)
        end if
        call extend(dict, key .kvp. inner_dict)
      else
        call put_string(dict, key, value)
      end if
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Reads config file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function read_config_dict(config_path) result(config)
    character(*), intent(in) :: config_path
    type(dictionary_t) :: config
    integer :: file_unit, line_num
    
    open(newunit = file_unit, file = config_path)
    line_num = 0
    config = read_inner_dict(file_unit, 1, line_num, '')
    close(file_unit)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks that only one of the keys in the specified group is set.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_only_one_set(config_dict, group, prefix)
    class(dictionary_t) :: config_dict ! intent(in)
    class(string), intent(in) :: group(:)
    character(*), optional, intent(in) :: prefix
    logical :: is_set
    integer :: i
    character(:), allocatable :: next_key
    type(string), allocatable :: group_repr(:)

    group_repr = group
    if (present(prefix)) then
      do i = 1, size(group)
        group_repr(i) = prefix // group_repr(i) % to_char_str()
      end do
    end if

    is_set = .false.
    do i = 1, size(group)
      next_key = group(i) % to_char_str()
      if (next_key .in. config_dict) then
        call assert(.not. is_set, 'Only one of the following keys can be set: ' // string_arr_to_char_str(group_repr))
        is_set = .true.
      end if
    end do
    call assert(is_set, 'One of the following keys has be to specified: ' // string_arr_to_char_str(group_repr))
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks whether all mandatory keys of this stage are specified in config.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_mandatory_keys(config_dict, mandatory_keys, prefix)
    class(dictionary_t) :: config_dict, mandatory_keys ! intent(in)
    character(*), optional, intent(in) :: prefix
    logical :: all_set
    integer :: i
    character(:), allocatable :: next_key, next_key_repr
    type(string), allocatable :: mandatory_keys_plain(:)

    mandatory_keys_plain = get_key_set(mandatory_keys)
    all_set = .true.
    do i = 1, size(mandatory_keys_plain)
      next_key = mandatory_keys_plain(i) % to_char_str()
      next_key_repr = next_key
      if (present(prefix)) then
        next_key_repr = prefix // next_key_repr
      end if
      if (.not. (next_key .in. config_dict)) then
        all_set = .false.
        call print_parallel('Error: the following key has to be specified: ' // next_key_repr)
      end if
    end do
    call assert(all_set, 'Some of the mandatory keys have not been specified')
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks whether some of the specified keys will be unused.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_extra_keys(config_dict, all_keys, prefix)
    class(dictionary_t) :: config_dict, all_keys ! intent(in)
    character(*), optional, intent(in) :: prefix
    integer :: i
    character(:), allocatable :: next_key, next_key_repr
    type(string), allocatable :: config_keys(:)

    config_keys = get_key_set(config_dict)
    do i = 1, size(config_keys)
      next_key = config_keys(i) % to_char_str()
      next_key_repr = next_key
      if (present(prefix)) then
        next_key_repr = prefix // next_key_repr
      end if
      call assert(next_key .in. all_keys, 'Error: the following key is not recognized: ' // next_key_repr)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks whether some of the specified keys will be unused.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_unused_keys(config_dict, mandatory_keys, optional_keys, prefix)
    class(dictionary_t) :: config_dict, mandatory_keys, optional_keys ! intent(in)
    character(*), optional, intent(in) :: prefix
    integer :: i
    character(:), allocatable :: next_key, next_key_repr
    type(string), allocatable :: config_keys(:)

    config_keys = get_key_set(config_dict)
    do i = 1, size(config_keys)
      next_key = config_keys(i) % to_char_str()
      next_key_repr = next_key
      if (present(prefix)) then
        next_key_repr = prefix // next_key_repr
      end if
      if (.not. ((next_key .in. mandatory_keys) .or. (next_key .in. optional_keys))) then
        call print_parallel('Info: the following key is not used: ' // next_key_repr)
      end if
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Announces default values of unset keys.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine announce_defaults(optional_nonset_keys, messages, prefix)
    class(dictionary_t), intent(in) :: optional_nonset_keys
    class(dictionary_t), optional, intent(in) :: messages
    character(*), optional, intent(in) :: prefix
    integer :: i
    character(:), allocatable :: next_key, next_key_repr, message
    type(string), allocatable :: keys(:)

    keys = get_key_set(optional_nonset_keys)
    do i = 1, size(keys)
      next_key = keys(i) % to_char_str()
      next_key_repr = next_key
      if (present(prefix)) then
        next_key_repr = prefix // next_key_repr
      end if
      message = 'Assuming default = ' // extract_string(optional_nonset_keys, next_key)
      if (present(messages)) then
        if (next_key .in. messages) then
          message = extract_string(messages, next_key)
          if (message == '') then
            cycle
          end if
        end if
      end if
      call print_parallel('Info: ' // next_key_repr // ' was not specified. ' // message)
    end do
  end subroutine

end module
