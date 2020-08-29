module io_base_mod
  use general_utils
  
  interface myprint
    module procedure :: print_element, print_array, print_matrix
  end interface
  
contains

!-----------------------------------------------------------------------
! Counts number of lines in a file
!-----------------------------------------------------------------------
  integer function countLines(filePath)
    character(len = *), intent(in) :: filePath
    integer :: fileState, fileUnit

    open(newunit = fileUnit, file = filepath, status = 'old')
    countLines = 0
    do
      read(fileUnit, *, iostat = fileState)
      if (fileState < 0) then
        return
      end if
      countLines = countLines + 1
    end do
    close(fileUnit)
  end function countLines
  
!-----------------------------------------------------------------------
! Prints matrix neatly
!-----------------------------------------------------------------------
  subroutine print_real_matrix(matrix, fileName, append, transpose, print_size)
    real*8 :: matrix(:, :)
    character(len = *), optional :: fileName
    integer, optional :: append, transpose, print_size
    integer :: append_act, transpose_act, print_size_act
    integer :: fileUnit, i, j

    append_act = merge(append, 0, present(append))
    transpose_act = merge(transpose, 0, present(transpose))
    print_size_act = merge(print_size, 0, present(print_size))
    if (present(fileName)) then
      if (append_act) then
        open(newunit = fileUnit, file = fileName, position = 'append')
      else
        open(newunit = fileUnit, file = fileName)
      end if
    else
      fileUnit = output_unit
    end if

    if (print_size_act == 1) then
      write(fileUnit, '(I0,x,I0)') size(matrix, 1), size(matrix, 2)
    end if
    if (transpose_act) then
      do j = 1,size(matrix, 2)
        write(fileUnit, '(*(E26.17))') (matrix(i, j), i = 1,size(matrix, 1))
      end do
    else
      do i = 1,size(matrix, 1)
        write(fileUnit, '(*(E26.17))') (matrix(i, j), j = 1,size(matrix, 2))
      end do
    end if
    close(fileUnit)
  end subroutine

!---------------------------------------------------------------------------------------------------------------------------------------------
! Reads a real matrix printed by print_real_matrix with print_size=1
!---------------------------------------------------------------------------------------------------------------------------------------------
  function read_real_matrix(path) result(res)
    character(*), intent(in) :: path
    real*8, allocatable :: res(:, :)
    integer :: rows, cols, file_unit

    open(newunit = file_unit, file = path)
    read(file_unit, *) rows, cols
    allocate(res(cols, rows))
    read(file_unit, *) res
    close(file_unit)
    res = transpose(res)
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Prints complex matrix as A+iB, where A and B are printed separately
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_complex_matrix(matrix, fileName, append, transpose, print_size)
    complex*16 :: matrix(:, :)
    character(len = *), optional :: fileName
    integer, optional :: append, transpose, print_size

    call print_real_matrix(real(matrix), fileName // '_real', append, transpose, print_size)
    call print_real_matrix(aimag(matrix), fileName // '_imag', append, transpose, print_size)
  end subroutine

!---------------------------------------------------------------------------------------------------------------------------------------------
! Reads a complex matrix printed by print_complex_matrix with print_size=1
!---------------------------------------------------------------------------------------------------------------------------------------------
  function read_complex_matrix(path) result(res)
    character(*), intent(in) :: path
    complex*16, allocatable :: res(:, :)
    real*8, allocatable :: real_part(:, :), imag_part(:, :)

    real_part = read_real_matrix(path // '_real')
    imag_part = read_real_matrix(path // '_imag')
    res = dcmplx(real_part, imag_part)
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Core subroutine for print subroutines
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_element(smth, file_unit)
    class(*) :: smth
    integer, value :: file_unit

    select type(smth)
    type is (integer)
      write(file_unit, '(3x,I0$)') smth
    type is (real*8)
      write(file_unit, '(3x,G24.16$)') smth
    type is (complex*16)
      ! write(file_unit, '(3x,ES23.16,SP,ES23.16,A$)') real(smth), aimag(smth), 'i'
      write(file_unit, '(3x,A,G0.16,A,G0.16,A$)') '(', real(smth), ',', aimag(smth), ')'
    type is (character(*))
      write(file_unit, '(A$)') smth
    end select
  end subroutine
  
!-----------------------------------------------------------------------
! Prints matrix
!-----------------------------------------------------------------------
  subroutine print_matrix(matrix, file_name, append)
    class(*), intent(in) :: matrix(:,:)
    character(*), optional :: file_name
    integer, optional, intent(in) :: append
    integer :: file_unit, i, j

    file_unit = get_file_unit(file_name, append)
    do i = 1,size(matrix, 1)
      do j = 1,size(matrix, 2)
        call write_element(matrix(i, j), file_unit)
      end do
      write(file_unit, '(/)')
    end do
  end subroutine

!-----------------------------------------------------------------------
! Prints 1D array
!-----------------------------------------------------------------------
  subroutine print_array(array, file_name, append, vertical)
    class(*), intent(in) :: array(:)
    character(*), optional, intent(in) :: file_name
    integer, optional, intent(in) :: append, vertical
    integer :: vertical_act, file_unit, i

    vertical_act = arg_or_default(vertical, 1)
    file_unit = get_file_unit(file_name, append)
    do i = 1,size(array)
      call write_element(array(i), file_unit)
      if (vertical_act) then
        write(file_unit, '(/)')
      end if
    end do
    write(file_unit, '(/)')
    close(file_unit)
  end subroutine

!-----------------------------------------------------------------------
! Prints a single element
!-----------------------------------------------------------------------
  subroutine print_element(smth, file_name, append, new_line)
    class(*), intent(in) :: smth
    character(*), optional :: file_name
    integer, optional, intent(in) :: append, new_line
    integer :: new_line_act
    integer :: file_unit

    if (present(new_line)) then
      new_line_act = new_line
    else
      new_line_act = 1
    end if

    file_unit = get_file_unit(file_name, append)
    call write_element(smth, file_unit)
    if (new_line_act) then
      write(file_unit, '(/)')
    end if
    close(file_unit)
  end subroutine

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
! iostat contains operation result. -2 - normal read, -1 - EOF
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
  function read_file(file_path, destruct) result(content)
    character(*), intent(in) :: file_path
    integer, optional, intent(in) :: destruct ! Deletes file after read
    character(:), allocatable :: content
    integer :: file_unit, file_size, destruct_act

    open(newunit = file_unit, file = file_path, action = 'read', form = 'unformatted', access = 'stream')
    inquire(unit = file_unit, size = file_size)
    allocate(character(file_size) :: content)
    read(file_unit) content

    destruct_act = arg_or_default(destruct, 0)
    if (destruct_act == 0) then
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
