#include "funcs.macro"
use algorithms_mod
use mpi
use parallel_base_mod
implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Gathers local chunks of a 1D array on 0th processor.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(gather_array_1d_,TEMPLATE_TYPE_NAME)(proc_chunk, global, proc_sizes, proc_shifts)
    TEMPLATE_TYPE, intent(in) :: proc_chunk(:)
    TEMPLATE_TYPE_OUT, allocatable, intent(out) :: global(:)
    integer, allocatable, optional, intent(out) :: proc_sizes(:), proc_shifts(:)
    integer :: ierr
    integer, allocatable :: proc_sizes_act(:), proc_shifts_act(:)

    if (.not. is_mpi_enabled()) then
      global = proc_chunk
      return
    end if

    allocate(proc_sizes_act(get_num_procs()))
    call MPI_Allgather(size(proc_chunk), 1, MPI_INTEGER, proc_sizes_act, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    proc_shifts_act = prefix_sum_exclusive(proc_sizes_act)

    if (get_proc_id() == 0) then
      allocate(global(sum(proc_sizes_act)))
    else
      allocate(global(0))
    end if
    call MPI_Gatherv(proc_chunk, size(proc_chunk), MPI_TYPE_NAME, global, proc_sizes_act, proc_shifts_act, MPI_TYPE_NAME, 0, MPI_COMM_WORLD, ierr)

    if (present(proc_sizes)) then
      proc_sizes = proc_sizes_act
    end if
    if (present(proc_shifts)) then
      proc_shifts = proc_shifts_act
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Gathers local chunks of columns of a 2D array on 0th processor.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(gather_cols_array_2d_,TEMPLATE_TYPE_NAME)(proc_chunk, global)
    TEMPLATE_TYPE, intent(in) :: proc_chunk(:, :)
    TEMPLATE_TYPE_OUT, allocatable, intent(out) :: global(:, :)
    integer :: global_cols
    TEMPLATE_TYPE, allocatable :: proc_chunk_1d(:), global_1d(:)

    proc_chunk_1d = [proc_chunk]
    call CONCAT2(gather_array_1d_,TEMPLATE_TYPE_NAME)(proc_chunk_1d, global_1d)
    global_cols = size(global_1d) / size(proc_chunk, 1)
    global = reshape(global_1d, [size(proc_chunk, 1), global_cols])
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Gathers local chunks of rows of a 2D array on 0th processor.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(gather_rows_array_2d_,TEMPLATE_TYPE_NAME)(proc_chunk, global)
    TEMPLATE_TYPE, intent(in) :: proc_chunk(:, :)
    TEMPLATE_TYPE_OUT, allocatable, intent(out) :: global(:, :)

    call CONCAT2(gather_cols_array_2d_,TEMPLATE_TYPE_NAME)(transpose(proc_chunk), global)
    global = transpose(global)
  end subroutine
