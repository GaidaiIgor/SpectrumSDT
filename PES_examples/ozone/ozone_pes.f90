module pesprint_constants
  use iso_fortran_env, only: real64
  implicit none

  real(real64), parameter :: pi = acos(-1d0)
  real(real64), parameter :: au_to_wn = 219474.6313708d0 ! cm^-1 / Eh (wavenumbers per Hartree)
  real(real64), parameter :: amu_to_kg = 1.660538921d-27 ! kg / amu (kilograms per atomic mass unit)
  real(real64), parameter :: aum_to_kg = 9.10938291d-31 ! kg / aum (kilogram per electron (atomic unit of mass))
  real(real64), parameter :: amu_to_aum = amu_to_kg / aum_to_kg ! aum / amu
  real(real64), parameter :: ozone_isotope_masses(3) = [15.99491461956d0, 16.99913170d0, 17.9991596129d0] * amu_to_aum ! aum; masses of isotopes of oxygen 16, 17 and 18 respectively

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
  use coordinate_coversion_mod
  implicit none

  integer :: ierr, proc_id
  real(real64) :: shift
  real(real64) :: mass(3)
  real(real64), allocatable :: grid_rho(:), grid_theta(:), grid_phi(:), pes(:)

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
  call init_parameters(mass, shift)

  call load_grids(grid_rho, grid_theta, grid_phi)
  call calc_pes(grid_rho, grid_theta, grid_phi, mass, shift, pes)

  if (proc_id == 0) then
    call print_pes(pes, 'pes.out')
  end if
  call MPI_Finalize(ierr)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes masses and shifts. Isotopic composition of ozone has to be supplied via command line argument to determine appropriate shift.
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
  subroutine load_grids(grid_rho, grid_theta, grid_phi)
    real(real64), allocatable, intent(out) :: grid_rho(:), grid_theta(:), grid_phi(:)
    integer :: file_unit, num_points, i

    open(newunit = file_unit, file = 'grid_rho.dat')
    read(file_unit, *) num_points
    allocate(grid_rho(num_points))
    do i = 1, num_points
      read(file_unit, *) grid_rho(i)
    end do
    close(file_unit)

    open(newunit = file_unit, file = 'grid_theta.dat')
    read(file_unit, *) num_points
    allocate(grid_theta(num_points))
    do i = 1, num_points
      read(file_unit, *) grid_theta(i)
    end do
    close(file_unit)

    open(newunit = file_unit, file = 'grid_phi.dat')
    read(file_unit, *) num_points
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
! Calculates ozone potential at given *APH coordinates*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point(aph, mass, shift) result(potential)
    real(real64), intent(in) :: aph(3), mass(3)
    real(real64), intent(in) :: shift
    real(real64) :: potential
    real(real64) :: internal(3)
    real(real64) :: all_bonds(3, 1)
    external :: IMLS ! external PES procedure

    all_bonds = convert_aph_to_all_bonds(reshape(aph, [3, 1]), mass)
    ! convert all bonds to internal coordinates expected by IMLS
    internal(1) = all_bonds(1, 1)
    internal(2) = all_bonds(2, 1)
    internal(3) = acos((all_bonds(1, 1)**2 + all_bonds(2, 1)**2 - all_bonds(3, 1)**2) / (2 * all_bonds(1, 1) * all_bonds(2, 1))) * 180 / pi ! 866 angle

    call IMLS(internal, potential, 1)
    potential = potential / au_to_wn - De + shift
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential energy surface for each combination of points in *grids*. Result is defined on 0th processor only. Parallel version.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_pes(grid_rho, grid_theta, grid_phi, mass, shift, pes)
    real(real64), intent(in) :: grid_rho(:), grid_theta(:), grid_phi(:)
    real(real64), intent(in) :: mass(3)
    real(real64), intent(in) :: shift
    real(real64), allocatable, intent(out) :: pes(:)
    real(real64), allocatable :: proc_pes(:)
    integer :: total_points, proc_id, ierr, first_k, proc_points, proc_k, k, ind_rho, ind_theta, ind_phi
    integer, allocatable :: proc_counts(:), proc_shifts(:)

    total_points = size(grid_rho) * size(grid_theta) * size(grid_phi)
    call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
    ! Split all points between available processors
    call get_proc_elem_range(total_points, first_k, proc_points, proc_counts, proc_shifts)
    allocate(proc_pes(proc_points))

    ! Calculate this proc's share
    do proc_k = 1, proc_points
      k = first_k + proc_k - 1
      call convert_1d_ind_to_3d(k, size(grid_rho), size(grid_theta), size(grid_phi), ind_rho, ind_theta, ind_phi)
      proc_pes(proc_k) = calc_potential_point([grid_rho(ind_rho), grid_theta(ind_theta), grid_phi(ind_phi)], mass, shift)
    end do

    ! Allocate global storage on 0th proc
    if (proc_id == 0) then
      allocate(pes(total_points))
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
    integer :: file_unit

    open(newunit = file_unit, file = file_name)
    write(file_unit, '(G23.15)') pes
    close(file_unit)
  end subroutine

end program
