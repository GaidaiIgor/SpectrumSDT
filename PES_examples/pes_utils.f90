module pes_utils
  use mpi
  use iso_fortran_env, only: real64
  implicit none

  real(real64), parameter :: pi = acos(-1d0)
  real(real64), parameter :: au_to_wn = 219474.6313708d0 ! cm^-1 / Eh (wavenumbers per Hartree)
  real(real64), parameter :: amu_to_kg = 1.660538921d-27 ! kg / amu (kilograms per atomic mass unit)
  real(real64), parameter :: aum_to_kg = 9.10938291d-31 ! kg / aum (kilogram per electron (atomic unit of mass))
  real(real64), parameter :: amu_to_aum = amu_to_kg / aum_to_kg ! aum / amu

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
! Loads requested coordinates.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_grids_aph(grid_rho, grid_theta, grid_phi)
    real(real64), allocatable, intent(out) :: grid_rho(:), grid_theta(:), grid_phi(:)
    integer :: file_unit, num_points, i
    real(real64) :: skip

    open(newunit = file_unit, file = 'grid_rho.dat')
    read(file_unit, *) skip, skip, skip, num_points
    allocate(grid_rho(num_points))
    do i = 1, num_points
      read(file_unit, *) grid_rho(i)
    end do
    close(file_unit)

    open(newunit = file_unit, file = 'grid_theta.dat')
    read(file_unit, *) skip, skip, skip, num_points
    allocate(grid_theta(num_points))
    do i = 1, num_points
      read(file_unit, *) grid_theta(i)
    end do
    close(file_unit)

    open(newunit = file_unit, file = 'grid_phi.dat')
    read(file_unit, *) skip, skip, skip, num_points
    allocate(grid_phi(num_points))
    do i = 1, num_points
      read(file_unit, *) grid_phi(i)
    end do
    close(file_unit)
  end subroutine
  
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
! Converts a 1D index of 1D-representation into 3 indexes corresponding to it in 3D representation.
! Assumes the 3D array was flattened using the following order of dimenisions: 3, 2, 1 (3rd coordinate is changing most frequently).
! n1, n2, n3 - sizes of the 3D array.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine convert_1d_ind_to_3d(ind_1d, n1, n2, n3, i1_3d, i2_3d, i3_3d)
    integer, intent(in) :: ind_1d, n1, n2, n3
    integer, intent(out) :: i1_3d, i2_3d, i3_3d
    integer :: ind_1d_0, i1_3d_0, i2_3d_0, i3_3d_0

    ! Convert using 0-based indexing
    ind_1d_0 = ind_1d - 1
    i1_3d_0 = ind_1d_0 / (n2 * n3)
    i2_3d_0 = (ind_1d_0 - i1_3d_0 * n2 * n3) / n3
    i3_3d_0 = ind_1d_0 - i1_3d_0 * n2 * n3 - i2_3d_0 * n3

    ! Shift back to 1-based indexing
    i1_3d = i1_3d_0 + 1
    i2_3d = i2_3d_0 + 1
    i3_3d = i3_3d_0 + 1
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints PES to file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_pes(pes, file_name)
    real(real64), intent(in) :: pes(:)
    character(*), intent(in) :: file_name
    integer :: file_unit

    open(newunit = file_unit, file = file_name)
    write(file_unit, '(G23.15)') pes
    close(file_unit)
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

end module
