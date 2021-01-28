program so2_pes
  use iso_fortran_env, only: real64
  use mpi
  use pes_utils
  implicit none

  integer :: ierr, proc_id
  real(real64), allocatable :: pes(:)
  real(real64), allocatable :: coords(:, :)

  call MPI_Init(ierr)
  coords = load_coords('pes.in')
  call calc_pes(coords, pes)

  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
  if (proc_id == 0) then
    call print_pes(pes, 'pes.out')
    print *, 'Done'
  end if
  call MPI_Finalize(ierr)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential of Klos et al. at a configuration given by internal coordinates:
! bond O1-S (Bohr), bond S-O2 (Bohr), O1-S-O2 angle (rad).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point(coords) result(potential)
    real(real64), intent(in) :: coords(3)
    real(real64) :: potential
    real(real64), parameter :: eh_per_ev = 0.036749323814708d0
    external :: so2Xpes
    call so2Xpes(coords(1), coords(2), cos(coords(3)), potential)
    potential = potential * eh_per_ev
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential energy surface for each combination of points in *grids*. Result is defined on 0th processor only. Parallel version.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_pes(coords, pes)
    real(real64), intent(in) :: coords(:, :)
    real(real64), allocatable, intent(out) :: pes(:)
    integer :: proc_id, ierr, first_k, proc_points, proc_k, k
    integer, allocatable :: proc_counts(:), proc_shifts(:)
    real(real64), allocatable :: proc_pes(:)
    external :: potreadX

    call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
    ! Split all points between available processors
    call get_proc_elem_range(size(coords, 2), first_k, proc_points, proc_counts, proc_shifts)
    allocate(proc_pes(proc_points))

    ! Read data
    call potreadX()
    ! Calculate this proc's share
    do proc_k = 1, proc_points
      if (proc_id == 0) then
        call track_progress(proc_k * 1d0 / proc_points, 0.01d0)
      end if
      k = first_k + proc_k - 1
      proc_pes(proc_k) = calc_potential_point(coords(:, k))
    end do

    ! Allocate global storage on 0th proc
    if (proc_id == 0) then
      allocate(pes(size(coords, 2)))
    end if
    print *, 'Proc', proc_id, 'is done'
    call MPI_Gatherv(proc_pes, size(proc_pes), MPI_DOUBLE_PRECISION, pes, proc_counts, proc_shifts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine

end program
