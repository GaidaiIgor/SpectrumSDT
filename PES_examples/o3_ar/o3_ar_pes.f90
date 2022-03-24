!-------------------------------------------------------------------------------------------------------------------------------------------
! Reads grid files and computes corresponding ozone-argonne potential using Dawes's PES.
!-------------------------------------------------------------------------------------------------------------------------------------------
program pesprint
  use iso_fortran_env, only: real64
  use mpi
  use path_resolution_mod
  use pes_utils_mod
  implicit none
  
  integer :: ierr, proc_id
  real(real64), allocatable :: grid_rho(:), grid_theta(:), grid_phi(:), pes(:)
  real(real64), allocatable :: coords(:, :)
  character(:), allocatable :: data_path, pes_input_path_1, pes_data_path_1, pes_input_path_2, pes_data_path_2

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
  data_path = resolve_relative_exe_path('../dawes/data/')
  pes_input_path_1 = data_path // '/input-AUTOSURF-PES.dat'
  pes_data_path_1 = data_path // '/PES-O3-Ar-2112'
  pes_input_path_2 = data_path // '/input-AUTOSURF-PES2.dat'
  pes_data_path_2 = data_path // '/PES-O3-Ar-600'
  
  call load_grids_aph(grid_rho, grid_theta, grid_phi)
  coords = aph_grids_to_coord_list(grid_rho, grid_theta, grid_phi)
  call calc_pes(coords, calc_potential_point, pes = pes)

  if (proc_id == 0) then
    call print_pes(pes, 'pes_out.txt')
  end if
  call MPI_Finalize(ierr)
  
contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates O3 - Ar potential at given "aph" coordinates.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point(aph, mass, shift) result(potential)
    real(real64), intent(in) :: aph(3)
    real(real64), optional, intent(in) :: mass(3)
    real(real64), optional, intent(in) :: shift
    real(real64) :: potential
    real(real64) :: spherical(3)
    external :: PES_LR ! PES procedure

    spherical(1) = aph(1) * m_per_a0 * 1d10 ! R, in Angstroms
    spherical(2) = cos(aph(2)) ! theta, in cos(theta)
    spherical(3) = aph(3) / pi * 180 ! phi, in degrees

    call PES_LR(spherical, potential, pes_input_path_1, pes_data_path_1, pes_input_path_2, pes_data_path_2)
    potential = potential / au_to_wn
  end function
  
end program
