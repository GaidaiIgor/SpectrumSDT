module parallel_utils
  use algorithms_mod
  use general_utils
  use iso_fortran_env, only: real64
  use mpi
  use parallel_base_mod
  use parallel_integer_mod
  use parallel_real_mod
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes range of elements to process by this processor.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_proc_elem_range(n_elems, first_elem, proc_elems, all_counts, all_shifts)
    integer, intent(in) :: n_elems ! total number of elements to be processed in parallel
    integer, intent(out) :: first_elem ! index of the first element that should be processed by this processor
    integer, intent(out) :: proc_elems ! number of elements to be processed by this processor
    integer, allocatable, optional, intent(out) :: all_counts(:) ! number of elements to be processed by all processors
    integer, allocatable, optional, intent(out) :: all_shifts(:) ! shifts of 1st element of all processors with respect to the global 1st element
    integer, allocatable :: all_counts_act(:)
    integer :: proc_id, n_procs, remaining_elems

    proc_id = get_proc_id()
    n_procs = get_num_procs()
    proc_elems = n_elems / n_procs
    remaining_elems = mod(n_elems, n_procs)
    if (proc_id < remaining_elems) then
      proc_elems = proc_elems + 1 ! distribute remaining elements
    end if
    first_elem = n_elems / n_procs * proc_id + min(proc_id, remaining_elems) + 1

    if (present(all_counts) .or. present(all_shifts)) then
      allocate(all_counts_act(n_procs))
      all_counts_act = n_elems / n_procs
      all_counts_act(1:remaining_elems) = all_counts_act(1:remaining_elems) + 1
    end if

    if (present(all_counts)) then
      all_counts = all_counts_act
    end if

    if (present(all_shifts)) then
      all_shifts = prefix_sum_exclusive(all_counts_act)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints using only 0th process.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_parallel(printable)
    class(*), intent(in) :: printable
    if (get_proc_id() == 0) then
      select type (printable)
      type is (character(*))
        print '(A)', printable
      end select
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Gathers local chunks of rows to a global 2D array on 0th processor.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine gather_rows(chunk, counts, shifts, global)
    real(real64), intent(in) :: chunk(:, :)
    integer, intent(in) :: counts(:), shifts(:)
    real(real64), allocatable, intent(out) :: global(:, :)
    integer :: ind, ierr

    ! Avoid calling MPI functions if MPI is not initialized
    if (.not. is_mpi_enabled()) then
      global = chunk
      return
    end if
    
    ! The result is collected and allocated on 0th proc only
    if (get_proc_id() == 0) then
      allocate(global(sum(counts), size(chunk, 2)))
    end if

    ! Gather results
    do ind = 1, size(chunk, 2)
      call MPI_Gatherv(chunk(1, ind), size(chunk, 1), MPI_DOUBLE_PRECISION, global(1, ind), counts, shifts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end do
  end subroutine

end module
