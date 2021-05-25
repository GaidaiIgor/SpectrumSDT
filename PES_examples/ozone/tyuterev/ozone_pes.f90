module pesprint_constants
  use iso_fortran_env, only: real64
  implicit none

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
  use pes_utils
  use pesprint_constants
  implicit none

  integer :: ierr, proc_id
  real(real64) :: shift
  real(real64), allocatable :: grid_rho(:), grid_theta(:), grid_phi(:), pes(:)

  integer :: i
  real*8 :: rho, tet, phi, potential
  external :: NR_PES_O3_2013_hyps ! external PES procedure

  ! Tyuterev's definition of phi is shifted by pi/2 relative to ours
  do i = 0, 100
    rho = 2d0 + i*0.01d0
    tet = 51 / 180d0 * pi
    phi = 90d0 / 180d0 * pi

    call NR_PES_O3_2013_hyps(rho, tet, phi, potential)
    rho = 2d0 + i*0.01d0
    print *, rho, potential
  end do
  stop 'debug'

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
  call init_parameters(shift)

  call load_grids_aph(grid_rho, grid_theta, grid_phi)
  call calc_pes(grid_rho, grid_theta, grid_phi, shift, pes)

  if (proc_id == 0) then
    call print_pes(pes, 'pes_out.txt')
  end if
  call MPI_Finalize(ierr)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes masses and shifts. Isotopic composition of ozone has to be supplied via command line argument to determine appropriate shift.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_parameters(shift)
    real(real64), intent(out) :: shift
    integer :: i, atom_code
    character(3) :: ozone_isotopes

    ! Read command line arguments
    call get_command_argument(1, ozone_isotopes)
    if (len(trim(ozone_isotopes)) == 0) then
      stop 'Isotopic composition has to be specified'
    end if

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
    
    ! temporary
    shift = 0
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates ozone potential at given *APH coordinates*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point(aph, shift) result(potential)
    real(real64), intent(in) :: aph(3)
    real(real64), intent(in) :: shift
    real(real64) :: potential
    external :: NR_PES_O3_2013_hyps ! external PES procedure

    ! Tyuterev's definition of phi is shifted by pi/2 relative to ours
    call NR_PES_O3_2013_hyps(aph(1), aph(2), aph(3) - pi/2, potential)
    potential = potential / au_to_wn + shift
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential energy surface for each combination of points in *grids*. Result is defined on 0th processor only. Parallel version.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_pes(grid_rho, grid_theta, grid_phi, shift, pes)
    real(real64), intent(in) :: grid_rho(:), grid_theta(:), grid_phi(:)
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
      proc_pes(proc_k) = calc_potential_point([grid_rho(ind_rho), grid_theta(ind_theta), grid_phi(ind_phi)], shift)
    end do

    ! Allocate global storage on 0th proc
    if (proc_id == 0) then
      allocate(pes(total_points))
    else
      allocate(pes(0))
    end if
    print *, 'Proc', proc_id, 'is done'
    call MPI_Gatherv(proc_pes, size(proc_pes), MPI_DOUBLE_PRECISION, pes, proc_counts, proc_shifts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine

end program
