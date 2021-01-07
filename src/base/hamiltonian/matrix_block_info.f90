module matrix_block_info_mod
  use block_borders_mod
  use general_utils
  use iso_fortran_env, only: real64
  implicit none

  private
  public :: matrix_block_info

  interface matrix_block_info
    module procedure :: construct_matrix_block_info
  end interface

  type :: matrix_block_info
    integer :: rows = -1
    integer :: columns = -1
    integer :: block_row_ind = -1 ! row index of this block in parental block
    integer :: block_col_ind = -1 ! column index if this block in parental block
    type(block_borders) :: borders = block_borders(-1, -1, -1, -1)
    type(matrix_block_info), pointer :: subblocks(:, :) => null()
  contains
    procedure :: deallocate_recursive => deallocate_recursive_matrix_block_info
    final :: finalize_matrix_block_info
    procedure :: shift_to => shift_to_matrix_block_info
    procedure :: update_block_indexes => update_block_indexes_matrix_block_info
    procedure :: compute_offdiag_subblock_sizes_diag => compute_offdiag_subblock_sizes_diag_matrix_block_info
    procedure :: remove_above_row => remove_above_row_matrix_block_info
    procedure :: remove_below_row => remove_below_row_matrix_block_info
    procedure :: cut_rows_between => cut_rows_between_matrix_block_info
    procedure :: transpose => transpose_matrix_block_info
    procedure :: extract_from_matrix => extract_from_matrix_matrix_block_info
    procedure :: is_empty => is_empty_matrix_block_info
    procedure :: calculate_compressed_columns => calculate_compressed_columns_matrix_block_info
    procedure :: copy_without_subblocks => copy_without_subblocks_matrix_block_info
    procedure :: copy => copy_matrix_block_info
    procedure :: print => print_matrix_block_info
    procedure :: print_all => print_all_matrix_block_info
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Constructor.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function construct_matrix_block_info(first_row, first_col, rows, columns) result(instance)
    integer, intent(in) :: first_row, first_col, rows, columns
    type(matrix_block_info) :: instance

    instance % rows = rows
    instance % columns = columns
    instance % block_row_ind = 1
    instance % block_col_ind = 1
    instance % borders = block_borders(first_col, first_col + columns - 1, first_row, first_row + rows - 1)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Deallocates memory recursively for all subblocks down this tree.
!-------------------------------------------------------------------------------------------------------------------------------------------
  impure elemental subroutine deallocate_recursive_matrix_block_info(this)
    class(matrix_block_info), intent(inout) :: this
    if (associated(this % subblocks)) then
      deallocate(this % subblocks)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Destructor.
! Not elemental to allow for deallocation of arrays without triggering chain deallocation of the inner layers.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine finalize_matrix_block_info(this)
    type(matrix_block_info), intent(inout) :: this
    call this % deallocate_recursive()
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Moves upper left corner of a block to given row and col position i.e. shifts all borders positions for the block and all its subblocks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine shift_to_matrix_block_info(this, row, col)
    class(matrix_block_info), intent(inout) :: this
    integer, intent(in) :: row, col
    integer :: i, j, curr_row, curr_col

    this % borders = block_borders(col, col + this % columns - 1, row, row + this % rows - 1)
    if (.not. associated(this % subblocks)) then
      return
    end if

    curr_row = row ! coordinates of upper left corner of next block
    curr_col = col
    do i = 1, size(this % subblocks, 1)
      do j = 1, size(this % subblocks, 2)
        call this % subblocks(i, j) % shift_to(curr_row, curr_col)
        curr_col = curr_col + this % subblocks(i, j) % columns
      end do
      curr_row = curr_row + this % subblocks(i, 1) % rows
      curr_col = col
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Recursively updates indexes of all nested blocks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine update_block_indexes_matrix_block_info(this)
    class(matrix_block_info), intent(inout) :: this
    integer :: i, j

    if (.not. associated(this % subblocks)) then
      return
    end if

    do i = 1, size(this % subblocks, 1)
      do j = 1, size(this % subblocks, 2)
        this % subblocks(i, j) % block_row_ind = i
        this % subblocks(i, j) % block_col_ind = j
        call this % subblocks(i, j) % update_block_indexes()
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Recursively computes sizes of all offdiagonal subblocks of a given _offdiagonal_ block located at the intersection of the two parental
! blocks. Rows are taken from parent1, columns are taken from parent2.
! Parents must be blocks of the same level (i.e. having the same number of levels of subblocks).
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine compute_offdiag_subblock_sizes_offdiag(parent1, parent2, block)
    class(matrix_block_info), intent(in) :: parent1, parent2
    type(matrix_block_info), intent(out) :: block
    integer :: i, j

    block = matrix_block_info(1, 1, parent1 % rows, parent2 % columns)
    if (.not. associated(parent1 % subblocks)) then
      return
    end if

    allocate(block % subblocks(size(parent1 % subblocks, 1), size(parent2 % subblocks, 2)))
    do i = 1, size(block % subblocks, 1)
      do j = 1, size(block % subblocks, 2)
        call compute_offdiag_subblock_sizes_offdiag(parent1 % subblocks(i, j), parent2 % subblocks(i, j), block % subblocks(i, j))
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Recursively calculates sizes of all offdiagonal subblocks of a given _diagonal_ block based on known sizes of diagonal subblocks.
! The starting block must be diagonal.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine compute_offdiag_subblock_sizes_diag_matrix_block_info(this)
    class(matrix_block_info), intent(inout) :: this
    integer :: i, j

    if (.not. associated(this % subblocks)) then
      return
    end if

    ! move down to the lowest level
    do i = 1, size(this % subblocks, 1)
      call this % subblocks(i, i) % compute_offdiag_subblock_sizes_diag()
    end do

    do i = 1, size(this % subblocks, 1)
      do j = 1, size(this % subblocks, 2)
        ! skip already processed diagonal blocks
        if (i == j) then
          cycle
        end if

        call compute_offdiag_subblock_sizes_offdiag(this % subblocks(i, i), this % subblocks(j, j), this % subblocks(i, j))
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Provides block information for matrix where rows above the specified row are removed.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine remove_above_row_matrix_block_info(this, row)
    class(matrix_block_info), intent(inout) :: this
    integer, intent(in) :: row
    integer :: j, first_non_empty_row
    type(matrix_block_info), pointer :: new_subblocks(:, :)

    ! If this block is unaffected by cutting then do nothing
    if (this % borders % top >= row) then
      return
    end if

    ! Correct upper border and block size. Size <= 0 means the block is fully eliminated.
    this % borders % top = row
    this % rows = this % borders % bottom - this % borders % top + 1

    ! If there are no subblocks then there is nothing else to do
    if (.not. associated(this % subblocks)) then
      return
    end if

    ! Find first not completely eliminated block row index
    first_non_empty_row = findloc(this % subblocks(:, 1) % borders % bottom >= row, .true., 1)
    ! If all subblocks are eliminated
    if (first_non_empty_row == 0) then
      call this % deallocate_recursive()
      return
    end if
    ! Allocate new version of this layer and copy relevant data from old version
    allocate(new_subblocks, source = this % subblocks(first_non_empty_row:, :))
    ! Deallocate cutted blocks recursively
    call this % subblocks(:first_non_empty_row - 1, :) % deallocate_recursive()
    ! Deallocate the old version of this layer
    deallocate(this % subblocks)
    ! Associate the pointer with the new version of this layer
    this % subblocks => new_subblocks

    ! The first row is partially eliminated, need to iterate recursively
    do j = 1, size(this % subblocks, 2)
      call this % subblocks(1, j) % remove_above_row(row)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Provides block information for matrix where rows below the specified row are removed.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine remove_below_row_matrix_block_info(this, row)
    class(matrix_block_info), intent(inout) :: this
    integer, intent(in) :: row
    integer :: j, last_non_empty_row
    type(matrix_block_info), pointer :: new_subblocks(:, :)
    ! Workaround for intel bug
    logical, allocatable :: temp(:)

    ! If this block is unaffected by cutting then do nothing
    if (this % borders % bottom <= row) then
      return
    end if

    ! Correct upper border and block size. Size <= 0 means the block is fully eliminated.
    this % borders % bottom = row
    this % rows = this % borders % bottom - this % borders % top + 1

    ! If there are no subblocks then there is nothing else to do
    if (.not. associated(this % subblocks)) then
      return
    end if

    ! Find last not completely eliminated block row index
    temp = this % subblocks(:, 1) % borders % top <= row ! Workaround for intel bug
    last_non_empty_row = findloc(temp, .true., 1, back = .true.)
    ! If all subblocks are eliminated
    if (last_non_empty_row == 0) then
      call this % deallocate_recursive()
      return
    end if

    ! Allocate new version of this layer and copy relevant data from old version
    allocate(new_subblocks, source = this % subblocks(:last_non_empty_row, :))
    ! Deallocate cutted blocks recursively
    call this % subblocks(last_non_empty_row + 1:, :) % deallocate_recursive()
    ! Deallocate the old version of this layer
    deallocate(this % subblocks)
    ! Associate the pointer with the new version of this layer
    this % subblocks => new_subblocks

    ! The last row is partially eliminated, need to iterate recursively
    do j = 1, size(this % subblocks, 2)
      call this % subblocks(size(this % subblocks, 1), j) % remove_below_row(row)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Cuts block section between specified rows.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine cut_rows_between_matrix_block_info(this, first_row, last_row)
    class(matrix_block_info), intent(inout) :: this
    integer, intent(in) :: first_row, last_row

    call this % remove_above_row(first_row)
    call this % remove_below_row(last_row)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transposes a given block and its subblocks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine transpose_matrix_block_info(this)
    class(matrix_block_info), intent(inout) :: this
    integer :: i, j
    type(matrix_block_info), pointer :: transposed_subblocks(:, :)

    call swap(this % rows, this % columns)
    call swap(this % block_row_ind, this % block_col_ind)
    call swap(this % borders % top, this % borders % left)
    call swap(this % borders % bottom, this % borders % right)

    if (.not. associated(this % subblocks)) then
      return
    end if

    allocate(transposed_subblocks(size(this % subblocks, 2), size(this % subblocks, 1)))
    do i = 1, size(this % subblocks, 1)
      do j = 1, size(this % subblocks, 2)
        transposed_subblocks(j, i) = this % subblocks(i, j)
        call transposed_subblocks(j, i) % transpose()
      end do
    end do
    deallocate(this % subblocks)
    this % subblocks => transposed_subblocks
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a pointer to the part of the given matrix corresponding to this block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function extract_from_matrix_matrix_block_info(this, matrix) result(block)
    class(matrix_block_info), intent(in) :: this
    complex(real64), target, intent(in) :: matrix(:, :)
    complex(real64), pointer :: block(:, :)

    call assert(.not. this % is_empty(), 'Error: attempt to extract from an empty block')
    block => matrix(this % borders % top : this % borders % bottom, this % borders % left : this % borders % right)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns part of the given matrix corresponding to this block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function is_empty_matrix_block_info(this) result(res)
    class(matrix_block_info), intent(in) :: this
    logical :: res
    res = this % rows == 0 .or. this % columns == 0
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates actual number of columns in local blocks after compression.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine calculate_compressed_columns_matrix_block_info(this)
    class(matrix_block_info), intent(inout) :: this
    integer :: j, i, block_row_ind
    integer, allocatable :: block_row_cols(:) ! columns in current block-row

    if (.not. associated(this % subblocks)) then
      return
    end if

    do j = 1, size(this % subblocks, 2)
      do i = 1, size(this % subblocks, 1)
        call this % subblocks(i, j) % calculate_compressed_columns()
      end do
    end do

    allocate(block_row_cols(size(this % subblocks, 1)))
    do block_row_ind = 1, size(this % subblocks, 1)
      block_row_cols(block_row_ind) = sum(this % subblocks(block_row_ind, :) % columns)
    end do
    this % columns = maxval(block_row_cols)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Makes a shallow copy of block_info structure without copying subblocks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine copy_without_subblocks_matrix_block_info(this, other)
    class(matrix_block_info), intent(out) :: this
    class(matrix_block_info), intent(in) :: other

    this % rows = other % rows
    this % columns = other % columns
    this % block_row_ind = other % block_row_ind
    this % block_col_ind = other % block_col_ind
    this % borders = other % borders
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Makes a deep copy.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine copy_matrix_block_info(this, other)
    class(matrix_block_info), intent(out) :: this
    class(matrix_block_info), intent(in) :: other
    integer :: i, j

    call this % copy_without_subblocks(other)
    if (.not. associated(other % subblocks)) then
      return
    end if

    allocate(this % subblocks(size(other % subblocks, 1), size(other % subblocks, 2))) ! Because mold does not work in gfortran in this case
    do j = 1, size(this % subblocks, 2)
      do i = 1, size(this % subblocks, 1)
        call this % subblocks(i, j) % copy(other % subblocks(i, j))
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints info about this block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_matrix_block_info(this)
    class(matrix_block_info), intent(in) :: this
    print *, 'Index', this % block_row_ind, this % block_col_ind
    print *, 'Size', this % rows, this % columns
    print *, 'Subblocks size', size(this % subblocks, 1), size(this % subblocks, 2)
    print *, 'Borders (t, b, l, r)', this % borders % top, this % borders % bottom, this % borders % left, this % borders % right
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints info about this block and all its subblocks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine print_all_matrix_block_info(this)
    class(matrix_block_info), intent(in) :: this
    integer :: i, j

    call this % print()
    if (.not. associated(this % subblocks)) then
      return
    end if

    do i = 1, size(this % subblocks, 1)
      do j = 1, size(this % subblocks, 2)
        print *
        print *, 'Printing subblock:', i, j
        call this % subblocks(i, j) % print_all()
      end do
    end do
  end subroutine

end module
