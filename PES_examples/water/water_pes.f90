program water_pes
  use iso_fortran_env, only: real64
  use mpi
  implicit none

  integer :: ierr, proc_id
  real(real64), allocatable :: coords(:, :)
  real(real64), allocatable :: pes(:)
  character(100) :: arg_buf
  character(:), allocatable :: pes_name

  call get_command_argument(1, arg_buf)
  pes_name = trim(adjustl(arg_buf))
  if (.not. any(pes_name == [character(100) :: 'shirin', 'partridge'])) then
    stop 'Unknown PES name'
  end if

  call MPI_Init(ierr)
  coords = load_coords('pes.in')
  call calc_pes(coords, pes_name, pes)

  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
  if (proc_id == 0) then
    call print_pes(pes, 'pes.out')
    print *, 'Done'
  end if
  call MPI_Finalize(ierr)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads requested coordinates.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function load_coords(file_name) result(coords) 
    character(*), intent(in) :: file_name
    real(real64), allocatable :: coords(:, :)
    integer :: file_unit, num_points

    open(newunit = file_unit, file = file_name)
    read(file_unit, *) num_points
    allocate(coords(3, num_points))
    read(file_unit, *) ! Skip header
    read(file_unit, *) coords
    close(file_unit)
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes exclusive prefix sum, i.e. res(i) = sum( [array(1)..array(i - 1)] ).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function prefix_sum_exclusive(array) result(res)
    integer, intent(in) :: array(:)
    integer, allocatable :: res(:)
    integer :: i
    
    allocate(res(size(array)))
    res(1) = 0
    do i = 2, size(array)
      res(i) = res(i - 1) + array(i - 1)
    end do
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes range of elements to process by this processor.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_proc_elem_range(n_elems, first_elem, proc_elems, all_counts, all_shifts)
    integer, intent(in) :: n_elems ! total number of elements to be processed in parallel
    integer, intent(out) :: first_elem ! index of the first element that should be processed by this processor
    integer, intent(out) :: proc_elems ! number of elements to be processed by this processor
    integer, allocatable, intent(out) :: all_counts(:) ! number of elements to be processed by all processors
    integer, allocatable, intent(out) :: all_shifts(:) ! shifts of 1st element of all processors with respect to the global 1st element
    integer :: proc_id, n_procs, remaining_elems, ierr

    call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, n_procs, ierr)
    proc_elems = n_elems / n_procs
    remaining_elems = mod(n_elems, n_procs)
    if (proc_id < remaining_elems) then
      proc_elems = proc_elems + 1 ! distribute remaining elements
    end if
    first_elem = n_elems / n_procs * proc_id + min(proc_id, remaining_elems) + 1

    allocate(all_counts(n_procs))
    all_counts = n_elems / n_procs
    all_counts(1:remaining_elems) = all_counts(1:remaining_elems) + 1
    all_shifts = prefix_sum_exclusive(all_counts)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Compares given reals with specified precision.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function compare_reals(a, b, comp_precision) result(res)
    real(real64), intent(in) :: a, b
    real(real64), intent(in) :: comp_precision
    integer :: res
    real(real64) :: difference

    difference = a - b
    if (abs(difference) < comp_precision) then
      res = 0
    else if (difference > 0) then
      res = 1
    else if (difference < 0) then
      res = -1
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts real to integer, rounding up if within target accuracy of the next integer, otherwise down.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function real2int(a, comp_precision) result(res)
    real(real64), intent(in) :: a
    real(real64), intent(in) :: comp_precision
    integer :: res

    if (ceiling(a) - a < comp_precision) then
      res = ceiling(a)
    else
      res = int(a)
    end if
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Prints progress % at specified points.
! progress: current progress of some process (a number from 0 to 1).
! progress_step: controls how often progress should be reported.
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine track_progress(progress, progress_step)
    real(real64), intent(in) :: progress
    real(real64), intent(in) :: progress_step
    character, parameter :: carriage_return = achar(13)
    real(real64), save :: last_progress = 0

    if (compare_reals(progress, last_progress + progress_step, 1d-10) >= 0) then
      last_progress = real2int(progress / progress_step, 1d-10) * progress_step
      write(*, '(A,1x,F6.2,A)', advance = 'no') carriage_return, last_progress * 100, '% done'
      if (compare_reals(last_progress, 1d0, 1d-10) == 0) then
        print *
      end if
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential of Partridge et al. at a configuration given by internal coordinates:
! bond H1-O (Bohr), bond O-H2 (Bohr), H1-O-H2 angle (rad).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point_partridge(coords) result(potential)
    real(real64), intent(in) :: coords(3)
    real(real64) :: potential
    real(real64) :: target_coords(1, 3)
    external :: vibpot ! underlying PES procedure

    target_coords(1, :) = coords
    call vibpot(target_coords, potential, 1)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential of Shirin et al. at a configuration given by internal coordinates:
! bond H1-O (Bohr), bond O-H2 (Bohr), H1-O-H2 angle (rad).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point_shirin(coords) result(potential)
    real(real64), intent(in) :: coords(3)
    real(real64) :: potential
    real(real64) :: target_coords(3)
    external :: potv ! underlying PES procedure

    target_coords(1) = coords(1)
    target_coords(2) = coords(2)
    target_coords(3) = cos(coords(3))
    call potv(potential, target_coords(1), target_coords(2), target_coords(3))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential energy surface for each combination of points in *grids*. Result is defined on 0th processor only. Parallel version.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_pes(coords, pes_name, pes)
    real(real64), intent(in) :: coords(:, :)
    character(*), intent(in) :: pes_name
    real(real64), allocatable, intent(out) :: pes(:)
    integer :: proc_id, ierr, first_k, proc_points, proc_k, k
    integer, allocatable :: proc_counts(:), proc_shifts(:)
    real(real64), allocatable :: proc_pes(:)

    call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
    ! Split all points between available processors
    call get_proc_elem_range(size(coords, 2), first_k, proc_points, proc_counts, proc_shifts)
    allocate(proc_pes(proc_points))

    ! Calculate this proc's share
    do proc_k = 1, proc_points
      if (proc_id == 0) then
        call track_progress(proc_k * 1d0 / proc_points, 0.01d0)
      end if
      k = first_k + proc_k - 1

      select case (pes_name)
        case ('shirin')
          proc_pes(proc_k) = calc_potential_point_shirin(coords(:, k))
        case ('partridge')
          proc_pes(proc_k) = calc_potential_point_partridge(coords(:, k))
        case default
          stop 'Unknown PES name'
      end select
    end do

    ! Allocate global storage on 0th proc
    if (proc_id == 0) then
      allocate(pes(size(coords, 2)))
    end if
    print *, 'Proc', proc_id, 'is done'
    call MPI_Gatherv(proc_pes, size(proc_pes), MPI_DOUBLE_PRECISION, pes, proc_counts, proc_shifts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints PES to file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_pes(pes, file_name)
    real(real64), intent(in) :: pes(:)
    character(*), intent(in) :: file_name
    integer :: i, file_unit

    open(newunit = file_unit, file = file_name)
    do i = 1, size(pes)
      write(file_unit, '(G23.15)') pes(i)
    end do
    close(file_unit)
  end subroutine

end program
