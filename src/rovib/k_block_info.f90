module k_block_info
  use algorithms
  use io_utils
  use matrix_block_info_mod

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
!  Initializes sizes of diagonal n blocks
!  Only block sizes are computed, not positions
!-------------------------------------------------------------------------------------------------------------------------------------------
  function init_n_block_sizes_diag(block_sizes) result(n_blocks)
    integer, intent(in) :: block_sizes(:)
    type(matrix_block_info), allocatable :: n_blocks(:, :)
    integer :: i, rows

    allocate(n_blocks(size(block_sizes), size(block_sizes)))
    do i = 1, size(block_sizes)
      rows = block_sizes(i)
      n_blocks(i, i) = matrix_block_info(1, 1, rows, rows) ! diagonal blocks are square
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
!  Loads sizes of all overlap blocks for a given k-block
!-------------------------------------------------------------------------------------------------------------------------------------------
  function load_k_subblock_sizes_diag(block_info_path) result(k_block_info)
    character(*), intent(in) :: block_info_path
    type(matrix_block_info) :: k_block_info
    integer :: rows_total
    integer, allocatable :: block_sizes(:) ! overlap subblocks of a given block

    block_sizes = load_basis_size_2d(block_info_path) ! loads num of solutions kept in each slice
    rows_total = sum(block_sizes)
    k_block_info = matrix_block_info(1, 1, rows_total, rows_total)
    k_block_info % subblocks = init_n_block_sizes_diag(block_sizes)
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
!  Loads number of 2D states kept in each slice
!-------------------------------------------------------------------------------------------------------------------------------------------
  function load_basis_size_2d(basis_size_info_path) result(basis_size_2d)
    character(*), intent(in) :: basis_size_info_path
    integer, allocatable :: basis_size_2d(:)
    integer, allocatable :: basis_size_matrix(:, :)
    
    basis_size_matrix = read_matrix_integer(basis_size_info_path)
    basis_size_2d = basis_size_matrix(:, 2)
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
!  Loads matrix block info
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_block_info_nonadiabatic(n_channels, n_rho_points, n_blocks, block_sizes, offset, n_rows)
    integer, intent(in) :: n_channels ! number of channels
    integer, intent(in) :: n_rho_points ! number of points along rho
    integer, intent(out) :: n_blocks ! number of blocks along rows or columns
    integer, allocatable, intent(out) :: block_sizes(:) ! number of solutions kept in each rho slice
    integer, allocatable, intent(out) :: offset(:) ! starting row number of each block
    integer, intent(out) :: n_rows ! total number of rows in hamiltonian

    allocate(block_sizes(n_channels))
    n_blocks = n_channels
    block_sizes = n_rho_points
    offset = prefix_sum_exclusive(block_sizes)
    n_rows = sum(block_sizes)
  end subroutine
end module