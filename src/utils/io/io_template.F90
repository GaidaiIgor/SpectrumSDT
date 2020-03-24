#include "funcs.macro"
use io_base_mod
use string_mod
use string_utils
use vector_mod

contains

!---------------------------------------------------------------------------------------------------------------------------------------------
! Reads matrix from file
!---------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(read_matrix_, TEMPLATE_TYPE_NAME)(file_path) result(matrix)
    character(*), intent(in) :: file_path
    TEMPLATE_TYPE, allocatable :: matrix(:, :)
    type(CONCAT2(vector_array_1d_, TEMPLATE_TYPE_NAME)) :: matrix_vector
    integer :: file_unit
    character(:), allocatable :: next_line
    type(string), allocatable :: line_tokens(:)
    TEMPLATE_TYPE, allocatable :: next_row(:)

    matrix_vector = CONCAT2(vector_array_1d_, TEMPLATE_TYPE_NAME)()
    open(newunit = file_unit, file = file_path)
    do
      next_line = read_line(file_unit)
      if (next_line == '') then
        exit ! end of file
      end if
      line_tokens = strsplit(next_line) ! split by spaces
      line_tokens = pack(line_tokens, line_tokens % length() > 0) ! filter out empty tokens

      next_row = line_tokens % CONCAT2(to_, TEMPLATE_TYPE_NAME)() ! convert strings to target type
      call matrix_vector % push(CONCAT2(array_1d_, TEMPLATE_TYPE_NAME)(next_row)) ! wrap into the type and store in matrix
    end do
    close(file_unit)
    matrix = CONCAT3(vector_array_1d_, TEMPLATE_TYPE_NAME, _to_2D_array)(matrix_vector)
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes a matrix in binary form to the specified file
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(write_binary_matrix_, TEMPLATE_TYPE_NAME)(matrix, file_path)
    TEMPLATE_TYPE, intent(in) :: matrix(:, :)
    character(*), intent(in) :: file_path
    integer :: file_unit

    open(newunit = file_unit, file = file_path, form = 'unformatted')
    write(file_unit) size(matrix, 1), size(matrix, 2)
    write(file_unit) matrix
    close(file_unit)
  end subroutine

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes a matrix in binary form to the specified file
!---------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(read_binary_matrix_, TEMPLATE_TYPE_NAME)(file_path) result(matrix)
    character(*), intent(in) :: file_path
    TEMPLATE_TYPE, allocatable :: matrix(:, :)
    integer :: file_unit, rows, cols

    open(newunit = file_unit, file = file_path, form = 'unformatted')
    read(file_unit) rows, cols
    allocate(matrix(rows, cols))
    read(file_unit) matrix
    close(file_unit)
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes a matrix in binary form to the specified file
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(write_binary_vector_, TEMPLATE_TYPE_NAME)(matrix, file_path)
    TEMPLATE_TYPE, intent(in) :: matrix(:)
    character(*), intent(in) :: file_path
    integer :: file_unit

    open(newunit = file_unit, file = file_path, form = 'unformatted')
    write(file_unit) size(matrix)
    write(file_unit) matrix
    close(file_unit)
  end subroutine

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes a matrix in binary form to the specified file
!---------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(read_binary_vector_, TEMPLATE_TYPE_NAME)(file_path) result(matrix)
    character(*), intent(in) :: file_path
    TEMPLATE_TYPE, allocatable :: matrix(:)
    integer :: file_unit, rows

    open(newunit = file_unit, file = file_path, form = 'unformatted')
    read(file_unit) rows
    allocate(matrix(rows))
    read(file_unit) matrix
    close(file_unit)
  end function
