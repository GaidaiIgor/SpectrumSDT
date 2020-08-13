!-------------------------------------------------------------------------------------------------------------------------------------------
! Contains matmul operators and global variables necessary to evaluate them
! Operator signatures are restricted, so they cannot be made local
!-------------------------------------------------------------------------------------------------------------------------------------------
module matmul_operator_mod
  ! use blas95 ! matrix-vector product
  use distributed_rovib_hamiltonian_mod
  use input_params_mod
  use matrix_block_info_mod
  use mpi_f08
  implicit none

  abstract interface
    subroutine matmul_operator(proc_rows, vector, result)
      integer :: proc_rows
      complex*16 :: vector(proc_rows), result(proc_rows)
    end subroutine
    
    subroutine matmul_operator_real(proc_rows, vector, result)
      integer :: proc_rows
      real*8 :: vector(proc_rows), result(proc_rows)
    end subroutine
  end interface

  integer :: context ! BLACS context
  integer :: msize ! Hamiltonain size
  integer :: rog ! Row offset of the first row of this process
  real*8, allocatable :: hamd(:, :) ! Hamiltonian matrices
  complex*16, allocatable :: hamz(:, :)
  type(distributed_rovib_hamiltonian), target :: rovib_ham
  procedure(matmul_operator), pointer :: active_matmul_operator

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets module variables
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_matmul_variables(params)
    class(input_params), intent(in) :: params

    rog = rovib_ham % global_chunk_info % subblocks(1, 1) % borders % top - 1
    msize = rovib_ham % global_chunk_info % columns
    if (params % optimized_mult == 1) then
      active_matmul_operator => ham_mult_compressed
    else
      active_matmul_operator => ham_mult_old
    end if

    ! if (debug_mode == 'old_code_compatibility') then
    !   hamz = rovib_ham % proc_chunk
    ! end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Multiplies a single n-block of compressed hamiltonian with a given vector v
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine ham_mult_block_compressed(ham, v, local_block_info, v_start, v_proc_out)
    complex*16, intent(in) :: ham(:, :)
    complex*16, intent(in) :: v(:)
    class(matrix_block_info), intent(in) :: local_block_info
    integer, intent(in) :: v_start ! where the relevant part of v starts (block's global left border)
    complex*16, intent(inout) :: v_proc_out(:)
    integer :: row1, row2, col1, col2, i, j, v_shift

    row1 = local_block_info % borders % top
    row2 = local_block_info % borders % bottom
    col1 = local_block_info % borders % left
    col2 = local_block_info % borders % right
    v_shift = v_start - col1

    do j = col1, col2
      do i = row1, row2
        v_proc_out(i) = v_proc_out(i) + ham(i, j) * v(j + v_shift)
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Multiplies compressed Hamiltonian by a given vector parallelly
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine ham_mult_compressed(size_proc, v_proc_inp, v_proc_out)
    integer :: size_proc ! number of elements assigned to this processor
    complex*16 :: v_proc_inp(size_proc) ! This processor's chunk of the vector that we need to multiply with Hamiltonian
    complex*16 :: v_proc_out(size_proc) ! This processor's chunk of the result of multiplication
    integer :: ierr, K_row_ind, K_col_ind, n_row_ind, n_col_ind, v_start
    complex*16 :: v(msize)
    type(matrix_block_info), pointer :: local_K_block_info, global_K_block_info, local_n_block_info, global_n_block_info

    call MPI_Allgatherv(v_proc_inp, size(v_proc_inp), MPI_DOUBLE_COMPLEX, v, rovib_ham % all_counts, rovib_ham % all_shifts, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    v_proc_out = 0

    ! Iterate over all K-values of the current chunk
    do K_col_ind = 1, size(rovib_ham % local_chunk_info % subblocks, 2)
      do K_row_ind = 1, size(rovib_ham % local_chunk_info % subblocks, 1)
        local_K_block_info => rovib_ham % local_chunk_info % subblocks(K_row_ind, K_col_ind)
        global_K_block_info => rovib_ham % global_chunk_info % subblocks(K_row_ind, K_col_ind)

        ! -1 corresponds to empty K-blocks in the first and last 2 rows
        if (global_K_block_info % block_row_ind == global_K_block_info % block_col_ind .and. global_K_block_info % block_row_ind /= -1) then
          v_start = global_K_block_info % borders % left
          call ham_mult_block_compressed(rovib_ham % proc_chunk, v, local_K_block_info, v_start, v_proc_out)
        else
          do n_col_ind = 1, size(local_K_block_info % subblocks, 2)
            do n_row_ind = 1, size(local_K_block_info % subblocks, 1)
              local_n_block_info => local_K_block_info % subblocks(n_row_ind, n_col_ind)
              global_n_block_info => global_K_block_info % subblocks(n_row_ind, n_col_ind)
              v_start = global_n_block_info % borders % left
              call ham_mult_block_compressed(rovib_ham % proc_chunk, v, local_n_block_info, v_start, v_proc_out)
            end do
          end do
        end if

      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Plain old version of matrix-vector multiplication. Deprecated.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine ham_mult_old(size_proc, v_proc_inp, v_proc_out)
    integer :: size_proc ! number of elements assigned to this processor
    complex*16 :: v_proc_inp(size_proc) ! This processor's chunk of the vector that we need to multiply with Hamiltonian
    complex*16 :: v_proc_out(size_proc) ! This processor's chunk of the result of multiplication
    integer :: i, j, ierr
    complex*16 :: v(msize)

    if (size_proc == msize) then
      v = v_proc_inp
    else
      call MPI_Allgatherv(v_proc_inp, size(v_proc_inp), MPI_DOUBLE_COMPLEX, v, rovib_ham % all_counts, rovib_ham % all_shifts, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    end if

    v_proc_out = 0
    do i = 1, size(rovib_ham % proc_chunk, 1)
      do j = 1, size(rovib_ham % proc_chunk, 2)
        v_proc_out(i) = v_proc_out(i) + rovib_ham % proc_chunk(i, j) * v(j)
      end do
    end do
  end subroutine

! !-----------------------------------------------------------------------
! ! Hamiltonian operator for PARPACK. Complex version. rovib_coupling = 0 version.
! !-----------------------------------------------------------------------
!   subroutine opz(mloc, x, y)
!     integer mloc
!     complex*16 x(mloc)
!     complex*16 y(mloc)
!     complex*16 v(msize)
!
!     v = 0
!     v(rog+1:rog+mloc) = x
!     ! assembles full vector v by exchanging chunks between processes
!     call zgsum2d(context,'A',' ',msize,1,v,msize,-1,-1)
!     ! matrix-vector product: y = hamz * v
!     call gemv(hamz,v,y)
!   end subroutine

! !-----------------------------------------------------------------------
! ! Hamiltonian operator for PARPACK. Real version. rovib = 0 version. Deprecated.
! !-----------------------------------------------------------------------
!   subroutine opd(mloc, x, y)
!     integer mloc
!     real*8 x(mloc)
!     real*8 y(mloc)
!     real*8 v(msize)
!
!     v = 0
!     v(rog+1:rog+mloc) = x
!     call dgsum2d(context,'A',' ',msize,1,v,msize,-1,-1)
!     call gemv(hamd,v,y)
!   end subroutine

end module
