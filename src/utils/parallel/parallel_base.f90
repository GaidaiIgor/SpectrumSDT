module parallel_base_mod
  use mpi
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Wraps MPI_Initialized into a function.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function is_mpi_enabled() result(mpi_enabled)
    logical :: mpi_enabled
    integer :: ierr
    call MPI_Initialized(mpi_enabled, ierr)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns number of processors.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_num_procs() result(num_procs)
    integer :: num_procs
    integer :: ierr

    if (is_mpi_enabled()) then
      call MPI_Comm_Size(MPI_COMM_WORLD, num_procs, ierr)
    else
      num_procs = 1
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns current process id.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_proc_id() result(proc_id)
    integer :: proc_id
    integer :: ierr

    if (is_mpi_enabled()) then
      call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
    else
      proc_id = 0
    end if
  end function

end module
