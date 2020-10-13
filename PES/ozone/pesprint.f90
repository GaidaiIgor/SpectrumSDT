module pesprint_constants
  use iso_fortran_env, only: real64
  implicit none

  real(real64), parameter :: pi = acos(-1d0)
  real(real64), parameter :: amu_to_kg = 1.660538921d-27 ! kg / amu (kilograms per atomic mass unit)
  real(real64), parameter :: aum_to_kg = 9.10938291d-31 ! kg / aum (kilogram per electron (atomic unit of mass))
  real(real64), parameter :: amu_to_aum = amu_to_kg / aum_to_kg ! aum / amu
  real(real64), parameter :: ozone_isotope_masses(3) = [15.99491461956d0, 16.99913170d0, 17.9991596129d0] * amu_to_aum ! aum; masses of isotopes of oxygen 16, 17 and 18 respectively
  real(real64), parameter :: au_to_wn = 219474.6313708d0 ! cm^-1 / Eh (wavenumbers per Hartree)

  ! As computed by program zpeO2
  real(real64), parameter :: zpe_66 = 791.6373314827983d0 / au_to_wn ! J = 0
  ! real(real64), parameter :: zpe_66 = 794.5112391213112d0 / au_to_wn ! J = 1
  real(real64), parameter :: zpe_68 = 769.3702558802487d0 / au_to_wn ! J = 0
  ! real(real64), parameter :: zpe_68 = 772.0845812041375d0 / au_to_wn ! J = 1
  real(real64), parameter :: zpe_88 = 746.4311399198617d0 / au_to_wn ! J = 0
  ! real(real64), parameter :: zpe_88 = 748.9858447146220d0 / au_to_wn ! J = 1

  ! Obtained as energy of dissociated O3 on ozone PES of Dawes for r1 = 900 Bohr, angle = 90 deg (does not matter which one), r2 is optimized to 10^-15 precision (see dawes_de.f90)
  real(real64), parameter :: De = 9274.99560025014d0 / au_to_wn
end module

program pesprint
  use iso_fortran_env, only: real64
  use mpi
  use pesprint_constants
  implicit none

  integer :: ierr, proc_id
  real(real64) :: shift
  real(real64) :: mass(3)
  real(real64), allocatable :: coords(:, :)
  real(real64), allocatable :: pes(:)

  call MPI_Init(ierr)
  call init_parameters(mass, shift)
  coords = load_coords('pes.in')
  call calc_pes(coords, mass, shift, pes)

  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
  if (proc_id == 0) then
    call print_pes(pes, 'pes.out')
    print *, 'Done'
  end if
  call MPI_Finalize(ierr)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes masses and shifts. Isotopic composition of ozone has to be supplied via command line argument.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_parameters(mass, shift)
    real(real64), intent(out) :: mass(3)
    real(real64), intent(out) :: shift
    integer :: i, atom_code
    character(3) :: ozone_isotopes

    ! Read command line arguments
    call get_command_argument(1, ozone_isotopes)
    if (len(trim(ozone_isotopes)) == 0) then
      stop 'Isotopic composition has to be specified'
    end if

    do i = 1, 3
      read(ozone_isotopes(i:i), *) atom_code
      mass(i) = ozone_isotope_masses(atom_code - 5)
    end do
    
    ! mass(1) is the central atom; mass(2) and mass(3) are terminal atoms
    mass = cshift(mass, 1)

    ! select corresponding ZPE
    if (ozone_isotopes == '666') then
      shift = -zpe_66
    else if (ozone_isotopes == '686') then
      shift = -zpe_68
    else if (ozone_isotopes == '868') then
      shift = -zpe_88
    else
      stop 'Unsupported isotopic composition'
    end if
  end subroutine

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
! Computes exclusive prefix sum, i.e. res(i) = sum( [array(1)..array(i - 1)] )
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
! Compares given reals with specified precision
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
! Converts real to integer, rounding up if within target accuracy of the next integer, otherwise down
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
! progress: current progress of some process (a number from 0 to 1)
! progress_step: controls how often progress should be reported
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
! Calculates potential at a given *aph_point*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point_aph(aph_point, mass, shift) result(potential)
    real(real64), intent(in) :: aph_point(3)
    real(real64), intent(in) :: mass(3)
    real(real64), intent(in) :: shift
    real(real64) :: potential
    real(real64) :: r(3), r2(3), mass_copy(3)
    real(real64) :: cart(9)
    external :: INT_Cart, Cart_INT, IMLS ! Dawes' procedures

    r2 = aph_point
    mass_copy = mass
    call INT_Cart(r2, cart, mass_copy, 4) ! APH to Cartesian
    call Cart_INT(r2, cart, mass_copy, 2) ! Cartesian to bond-angle

    r(1) = min(r2(1), r2(2))
    r(2) = max(r2(1), r2(2))
    r(3) = 180*acos(r2(3)) / pi

    call IMLS(r, potential, 1)
    potential = potential / au_to_wn - De + shift
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential energy surface for each combination of points in *grids*. Result is defined on 0th processor only. Parallel version.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_pes(coords, mass, shift, pes)
    real(real64), intent(in) :: coords(:, :)
    real(real64), allocatable, intent(out) :: pes(:)
    real(real64) :: mass(3) ! intent(in)
    real(real64), intent(in) :: shift
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
      proc_pes(proc_k) = calc_potential_point_aph(coords(:, k), mass, shift)
    end do

    ! Allocate global storage on 0th proc
    if (proc_id == 0) then
      allocate(pes(size(coords, 2)))
      print *, 'Exchanging data...'
    end if
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
