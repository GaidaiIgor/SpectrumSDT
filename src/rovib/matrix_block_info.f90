module matrix_block_info_mod
  use block_borders_mod
  use general_utils
  private

  public :: copy_without_subblocks

  interface matrix_block_info
    module procedure :: construct_matrix_block_info
  end interface

  type, public :: matrix_block_info
    integer :: rows
    integer :: columns
    integer :: block_row_ind ! row index of this block in parental block
    integer :: block_col_ind ! column index if this block in parental block
    type(block_borders) :: borders
    type(matrix_block_info), allocatable :: subblocks(:, :)
  contains
    procedure :: shift_to
    procedure :: update_block_indexes
    procedure :: compute_offdiag_subblock_sizes_diag
    procedure :: remove_above_row
    procedure :: remove_below_row
    procedure :: cut_rows_between
    procedure :: transpose
    procedure :: extract_from_matrix
    procedure :: update_matrix
    procedure :: is_empty
    procedure :: print
    procedure :: print_all
    procedure :: calculate_compressed_columns
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Constructor
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
! Moves upper left corner of a block to given row and col position i.e. shifts all borders positions for the block and all its subblocks
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine shift_to(this, row, col)
    class(matrix_block_info), intent(inout) :: this
    integer, intent(in) :: row, col
    integer :: i, j, curr_row, curr_col

    this % borders = block_borders(col, col + this % columns - 1, row, row + this % rows - 1)
    if (.not. allocated(this % subblocks)) then
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
! Recursively updates indexes of all nested blocks
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine update_block_indexes(this)
    class(matrix_block_info), intent(inout) :: this
    integer :: i, j

    if (.not. allocated(this % subblocks)) then
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
! Recursively computes sizes of all offdiagonal subblocks of a given _offdiagonal_ block located at the intersection of the two parental blocks
! Rows are taken from parent1, columns are from parent2
! Parents must be blocks of the same level (i.e. having the same number of levels of subblocks)
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive function compute_offdiag_subblock_sizes_offdiag(parent1, parent2) result(block)
    class(matrix_block_info), intent(inout) :: parent1, parent2
    type(matrix_block_info) :: block
    integer :: i, j

    block = matrix_block_info(1, 1, parent1 % rows, parent2 % columns)
    if (.not. allocated(parent1 % subblocks)) then
      return
    end if

    allocate(block % subblocks(size(parent1 % subblocks, 1), size(parent2 % subblocks, 2)))
    do i = 1, size(block % subblocks, 1)
      do j = 1, size(block % subblocks, 2)
        block % subblocks(i, j) = compute_offdiag_subblock_sizes_offdiag(parent1 % subblocks(i, j), parent2 % subblocks(i, j))
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Recursively calculates sizes of all offdiagonal subblocks of a given _diagonal_ block based on known sizes of diagonal subblocks
! The starting block must be diagonal
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine compute_offdiag_subblock_sizes_diag(this)
    class(matrix_block_info), intent(inout) :: this
    integer :: i, j

    if (.not. allocated(this % subblocks)) then
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

        this % subblocks(i, j) = compute_offdiag_subblock_sizes_offdiag(this % subblocks(i, i), this % subblocks(j, j))
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Provides block information for matrix where rows above the specified row are removed
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine remove_above_row(this, row)
    class(matrix_block_info), intent(inout) :: this
    integer, intent(in) :: row
    integer :: i, j, first_non_empty_row

    ! If current block is affected by cutting
    if (this % borders % top < row) then
      ! Correct upper border and block size. Size <= 0 means the block is fully eliminated.
      this % borders % top = row
      this % rows = this % borders % bottom - this % borders % top + 1
    end if

    ! If there are no subblocks or the whole block is fully eliminated then do not go into subblocks
    if (.not. allocated(this % subblocks) .or. this % rows <= 0) then
      return
    end if

    ! If block is partially eliminated - iterate through its subblocks recursively
    do i = 1, size(this % subblocks, 1)
      do j = 1, size(this % subblocks, 2)
        call this % subblocks(i, j) % remove_above_row(row)
      end do
    end do

    ! Find first not completely eliminated block row index
    first_non_empty_row = findloc(this % subblocks(:, 1) % rows > 0, .true., 1)
    this % subblocks = this % subblocks(first_non_empty_row:, :)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Provides block information for matrix where rows below the specified row are removed
! See also remove_above_row for further comments
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine remove_below_row(this, row)
    class(matrix_block_info), intent(inout) :: this
    integer, intent(in) :: row
    integer :: i, j, last_non_empty_row
    integer, allocatable :: last_non_empty_block_indexes(:)

    if (this % borders % bottom > row) then
      this % borders % bottom = row
      this % rows = this % borders % bottom - this % borders % top + 1
    end if

    if (.not. allocated(this % subblocks) .or. this % rows <= 0) then
      return
    end if

    do i = 1, size(this % subblocks, 1)
      do j = 1, size(this % subblocks, 2)
        call this % subblocks(i, j) % remove_below_row(row)
      end do
    end do

    last_non_empty_block_indexes = findloc(this % subblocks % rows > 0, .true., 1, back = .true.) ! finds last non empty block in each column of blocks
    last_non_empty_row = last_non_empty_block_indexes(1)
    this % subblocks = this % subblocks(:last_non_empty_row, :)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Cuts block section between specified rows
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine cut_rows_between(this, first_row, last_row)
    class(matrix_block_info), intent(inout) :: this
    integer, intent(in) :: first_row, last_row

    call this % remove_above_row(first_row)
    call this % remove_below_row(last_row)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transposes a given block and its subblocks
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine transpose(this)
    class(matrix_block_info), intent(inout) :: this
    integer :: i, j
    type(matrix_block_info), allocatable :: transposed_subblocks(:, :)

    call swap(this % rows, this % columns)
    call swap(this % block_row_ind, this % block_col_ind)
    call swap(this % borders % top, this % borders % left)
    call swap(this % borders % bottom, this % borders % right)

    if (.not. allocated(this % subblocks)) then
      return
    end if

    allocate(transposed_subblocks(size(this % subblocks, 2), size(this % subblocks, 1)))
    ! this % subblocks = transpose(this % subblocks)
    do i = 1, size(this % subblocks, 1)
      do j = 1, size(this % subblocks, 2)
        transposed_subblocks(j, i) = this % subblocks(i, j)
      end do
    end do
    this % subblocks = transposed_subblocks

    do i = 1, size(this % subblocks, 1)
      do j = 1, size(this % subblocks, 2)
        call this % subblocks(i, j) % transpose()
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns part of the given matrix corresponding to this block
!-------------------------------------------------------------------------------------------------------------------------------------------
  function extract_from_matrix(this, matrix) result(block)
    class(matrix_block_info), intent(in) :: this
    complex*16, intent(in) :: matrix(:, :)
    complex*16, allocatable :: block(:, :)

    call assert(.not. this % is_empty(), 'Error: attempt to extract from an empty block')
    block = matrix(this % borders % top : this % borders % bottom, this % borders % left : this % borders % right)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Updates part of the given matrix corresponding to this block with a given block
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine update_matrix(this, matrix, block)
    class(matrix_block_info), intent(in) :: this
    complex*16, intent(inout) :: matrix(:, :)
    complex*16, intent(in) :: block(:, :)

    call assert(.not. this % is_empty(), 'Error: attempt to update an empty block')
    call assert(size(block, 1) == this % rows .and. size(block, 2) == this % columns, 'Error: Supplied block has incorrect dimensions')
    matrix(this % borders % top : this % borders % bottom, this % borders % left : this % borders % right) = block
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns part of the given matrix corresponding to this block
!-------------------------------------------------------------------------------------------------------------------------------------------
  function is_empty(this) result(res)
    class(matrix_block_info), intent(in) :: this
    logical :: res
    res = this % rows == 0 .or. this % columns == 0
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! prints info about this block
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print(this)
    class(matrix_block_info), intent(in) :: this
    print *, 'Index', this % block_row_ind, this % block_col_ind
    print *, 'Size', this % rows, this % columns
    print *, 'Subblocks size', size(this % subblocks, 1), size(this % subblocks, 2)
    print *, 'Borders (t, b, l, r)', this % borders % top, this % borders % bottom, this % borders % left, this % borders % right
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! prints info about this block and all its subblocks
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine print_all(this)
    class(matrix_block_info), intent(in) :: this
    integer :: i, j

    call this % print()
    if (.not. allocated(this % subblocks)) then
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

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates actual number of columns in local blocks after compression
!-------------------------------------------------------------------------------------------------------------------------------------------
  recursive subroutine calculate_compressed_columns(this)
    class(matrix_block_info), intent(inout) :: this
    integer :: j, i, block_row_ind
    integer, allocatable :: block_row_cols(:) ! columns in current block-row

    if (.not. allocated(this % subblocks)) then
      return
    end if

    do j = 1, size(this % subblocks, 2)
      do i = 1, size(this % subblocks, 1)
        call calculate_compressed_columns(this % subblocks(i, j))
      end do
    end do

    allocate(block_row_cols(size(this % subblocks, 1)))
    do block_row_ind = 1, size(this % subblocks, 1)
      block_row_cols(block_row_ind) = sum(this % subblocks(block_row_ind, :) % columns)
    end do
    this % columns = maxval(block_row_cols)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Makes a shallow copy of block_info structure without copying subblocks
!-------------------------------------------------------------------------------------------------------------------------------------------
  function copy_without_subblocks(other) result(copy)
    class(matrix_block_info), intent(in) :: other
    type(matrix_block_info) :: copy

    copy % rows = other % rows
    copy % columns = other % columns
    copy % block_row_ind = other % block_row_ind
    copy % block_col_ind = other % block_col_ind
    copy % borders = other % borders
  end function

end module
