module io_base_mod
  use general_utils
  use iso_fortran_env, only: output_unit
  implicit none
  
contains

!---------------------------------------------------------------------------------------------------------------------------------------------
! Opens a file and returns file unit
!---------------------------------------------------------------------------------------------------------------------------------------------
  integer function get_file_unit(file_name, position_mode)
    character(*), optional, intent(in) :: file_name
    integer, optional, intent(in) :: position_mode
    character(:), allocatable :: position_mode_string
    
    if (present(position_mode)) then
      if (position_mode == 0) then
        position_mode_string = 'asis'
      elseif (position_mode == 1) then
        position_mode_string = 'append'
      end if
    else
      position_mode_string = 'asis'
    end if

    if (present(file_name)) then
      open(newunit = get_file_unit, file = file_name, position = position_mode_string)
    else
      get_file_unit = output_unit
    end if
  end function
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Reads a line from specified unit. buf_size controls how much data is read at once. The actual line may be bigger than that.
! iostat contains operation result.
!---------------------------------------------------------------------------------------------------------------------------------------------
  function read_line(file_unit, iostat, buf_size) result(line)
    integer, intent(in) :: file_unit
    integer, optional, intent(out) :: iostat
    integer, optional, intent(in) :: buf_size
    character(:), allocatable :: line
    integer :: iostat_act, buf_size_act, chars_read
    character(:), allocatable :: buf
    
    buf_size_act = arg_or_default(buf_size, 256)
    allocate(character(buf_size_act) :: buf)
    
    line = ''
    do
      read(file_unit, '(A)', advance = 'no', iostat = iostat_act, size = chars_read) buf
      if (iostat_act > 0) then
        stop 'Error reading file'
      end if
      line = line // buf(:chars_read)
      if (iostat_act < 0) then
        if (present(iostat)) then
          iostat = iostat_act
        end if 
        return
      end if
    end do
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Reads an entire file specified by file_path
!---------------------------------------------------------------------------------------------------------------------------------------------
  function read_file(file_path, delete) result(content)
    character(*), intent(in) :: file_path
    integer, optional, intent(in) :: delete ! Deletes file after read
    character(:), allocatable :: content
    integer :: file_unit, file_size, delete_act

    open(newunit = file_unit, file = file_path, action = 'read', form = 'unformatted', access = 'stream')
    inquire(unit = file_unit, size = file_size)
    allocate(character(file_size) :: content)
    read(file_unit) content

    delete_act = arg_or_default(delete, 0)
    if (delete_act == 0) then
      close(file_unit)
    else
      close(file_unit, status = 'delete')
    end if
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Pads a string with appropriate number of spaces on both ends to center it in the given width
!---------------------------------------------------------------------------------------------------------------------------------------------
  function align_center(str, width) result(res)
    character(*), intent(in) :: str
    integer, intent(in) :: width
    character(:), allocatable :: res
    integer :: spaces, left_pad_size, right_pad_size

    spaces = width - len(str)
    left_pad_size = spaces / 2
    right_pad_size = spaces / 2 + mod(spaces, 2)
    res = repeat(' ', left_pad_size) // str // repeat(' ', right_pad_size)
  end function

end module
