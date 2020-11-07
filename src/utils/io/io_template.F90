#include "funcs.macro"
use io_base_mod
use string_mod
use string_utils
use vector_mod
implicit none

contains

!---------------------------------------------------------------------------------------------------------------------------------------------
! Reads matrix from file.
!---------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(read_matrix_,TEMPLATE_TYPE_NAME)(file_path) result(matrix)
    character(*), intent(in) :: file_path
    TEMPLATE_TYPE_OUT, allocatable :: matrix(:, :)
    type(CONCAT2(vector_array_1d_,TEMPLATE_TYPE_NAME)) :: matrix_vector
    integer :: file_unit
    character(:), allocatable :: next_line
    type(string), allocatable :: line_tokens(:)
    TEMPLATE_TYPE_OUT, allocatable :: next_row(:)

    matrix_vector = CONCAT2(vector_array_1d_,TEMPLATE_TYPE_NAME)()
    open(newunit = file_unit, file = file_path)
    do
      next_line = read_line(file_unit)
      if (next_line == '') then
        exit ! no empty lines are allowed
      end if
      line_tokens = strsplit(next_line) ! split by spaces
      line_tokens = pack(line_tokens, line_tokens % length() > 0) ! filter out empty tokens

      next_row = line_tokens % CONCAT2(to_,TEMPLATE_TYPE_NAME)() ! convert strings to target type
      call matrix_vector % push(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)(next_row)) ! wrap into the type and store in matrix
    end do
    close(file_unit)
    matrix = CONCAT3(vector_array_1d_,TEMPLATE_TYPE_NAME,_to_2D_array)(matrix_vector)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Writes *array* to *unit*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(write_array_unit_,TEMPLATE_TYPE_NAME)(array, unit, vertical)
    TEMPLATE_TYPE, intent(in) :: array(:)
    integer, intent(in) :: unit
    integer, optional, intent(in) :: vertical
    integer :: vertical_act
    character(:), allocatable :: count_spec, elem_spec, format_spec

    vertical_act = arg_or_default(vertical, 1)
    count_spec = iff(vertical_act == 0, num2str(size(array)), '1')
    elem_spec = TEMPLATE_ELEM_SPEC
    format_spec = '(' // count_spec // '(2x,' // elem_spec // '))'
    write(unit, format_spec) array
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Writes *array* to *file_name*, or screen if not specified.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(write_array_,TEMPLATE_TYPE_NAME)(array, file_name, append, vertical)
    TEMPLATE_TYPE, intent(in) :: array(:)
    character(*), optional, intent(in) :: file_name
    integer, optional, intent(in) :: append, vertical
    integer :: unit

    unit = get_file_unit(file_name, append)
    call CONCAT2(write_array_unit_,TEMPLATE_TYPE_NAME)(array, unit, vertical)
    if (unit /= output_unit) then
      close(unit)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Writes *matrix* to *file_name*, or screen if not specified.
! If *print_size* == 1 then size of the matrix is printed on the first line.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(write_matrix_,TEMPLATE_TYPE_NAME)(matrix, file_name, append, print_size)
    TEMPLATE_TYPE, intent(in) :: matrix(:, :)
    character(*), optional :: file_name
    integer, optional, intent(in) :: append, print_size
    integer :: unit, print_size_act, i

    unit = get_file_unit(file_name, append)
    print_size_act = arg_or_default(print_size, 0)
    if (print_size_act == 1) then
      write(unit, '(2I0)') size(matrix, 1), size(matrix, 2)
    end if
    do i = 1, size(matrix, 1)
      call CONCAT2(write_array_unit_,TEMPLATE_TYPE_NAME)(matrix(i, :), unit, vertical = 0)
    end do
    if (unit /= output_unit) then
      close(unit)
    end if
  end subroutine

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes a matrix in binary form to the specified file.
!---------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(read_binary_array_,TEMPLATE_TYPE_NAME)(file_path) result(matrix)
    character(*), intent(in) :: file_path
    TEMPLATE_TYPE_OUT, allocatable :: matrix(:)
    integer :: file_unit, rows

    open(newunit = file_unit, file = file_path, form = 'unformatted')
    read(file_unit) rows
    allocate(matrix(rows))
    read(file_unit) matrix
    close(file_unit)
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes a matrix in binary form to the specified file.
!---------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(read_binary_matrix_,TEMPLATE_TYPE_NAME)(file_path) result(matrix)
    character(*), intent(in) :: file_path
    TEMPLATE_TYPE_OUT, allocatable :: matrix(:, :)
    integer :: file_unit, rows, cols

    open(newunit = file_unit, file = file_path, form = 'unformatted')
    read(file_unit) rows, cols
    allocate(matrix(rows, cols))
    read(file_unit) matrix
    close(file_unit)
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes a matrix in binary form to the specified file.
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(write_binary_array_,TEMPLATE_TYPE_NAME)(matrix, file_path)
    TEMPLATE_TYPE, intent(in) :: matrix(:)
    character(*), intent(in) :: file_path
    integer :: file_unit

    open(newunit = file_unit, file = file_path, form = 'unformatted')
    write(file_unit) size(matrix)
    write(file_unit) matrix
    close(file_unit)
  end subroutine

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes a matrix in binary form to the specified file.
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(write_binary_matrix_,TEMPLATE_TYPE_NAME)(matrix, file_path)
    TEMPLATE_TYPE, intent(in) :: matrix(:, :)
    character(*), intent(in) :: file_path
    integer :: file_unit

    open(newunit = file_unit, file = file_path, form = 'unformatted')
    write(file_unit) size(matrix, 1), size(matrix, 2)
    write(file_unit) matrix
    close(file_unit)
  end subroutine
