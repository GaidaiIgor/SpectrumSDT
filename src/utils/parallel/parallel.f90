module parallel_utils
  use algorithms
  use general_utils
  use mpi_f08

contains
!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns current process id and total number of processes
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_proc_info(proc_id, n_procs)
    integer, intent(out) :: proc_id, n_procs
    integer :: ierr
    logical :: mpi_enabled

    call MPI_Initialized(mpi_enabled, ierr)
    if (mpi_enabled) then
      call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
      call MPI_Comm_Size(MPI_COMM_WORLD, n_procs, ierr)
      ! call blacs_pinfo(proc_id, n_procs)
    else
      proc_id = 0
      n_procs = 1
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns current process id
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_proc_id() result(proc_id)
    integer :: proc_id
    integer :: n_procs
    call get_proc_info(proc_id, n_procs)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes range of elements to process by this processor
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_proc_elem_range(n_elems, first_elem, proc_elems, all_counts, all_shifts)
    integer, intent(in) :: n_elems ! total number of elements to be processed in parallel
    integer, intent(out) :: first_elem ! index of the first element that should be processed by this processor
    integer, intent(out) :: proc_elems ! number of elements to be processed by this processor
    integer, allocatable, optional, intent(out) :: all_counts(:) ! number of elements to be processed by all processors
    integer, allocatable, optional, intent(out) :: all_shifts(:) ! shifts of 1st element of all processors with respect to the global 1st element
    integer, allocatable :: all_counts_act(:)
    integer :: proc_id, n_procs, remaining_elems

    call get_proc_info(proc_id, n_procs)
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
! Prints using only 0th process
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_parallel(printable)
    class(*), intent(in) :: printable
    integer :: proc_id, n_procs
    
    call get_proc_info(proc_id, n_procs)
    if (proc_id == 0) then
      select type (printable)
      type is (character(*))
        print '(A)', printable
      end select
    end if
  end subroutine
end module
