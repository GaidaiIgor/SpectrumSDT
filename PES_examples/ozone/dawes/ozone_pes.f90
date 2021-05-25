!-------------------------------------------------------------------------------------------------------------------------------------------
! Reads APH grid files and computes corresponding ozone potential using Dawe's PES.
!-------------------------------------------------------------------------------------------------------------------------------------------
program pesprint
  use iso_fortran_env, only: real64
  use mpi
  use ozone_utils_mod
  use pes_utils_mod
  implicit none

  integer :: ierr, proc_id
  real(real64) :: shift
  real(real64) :: mass(3)
  real(real64), allocatable :: grid_rho(:), grid_theta(:), grid_phi(:), pes(:)
  real(real64), allocatable :: coords(:, :)

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
  call init_parameters(mass, shift)

  call load_grids_aph(grid_rho, grid_theta, grid_phi)
  coords = aph_grids_to_coord_list(grid_rho, grid_theta, grid_phi)
  print *, size(coords, 1), size(coords, 2)
  call calc_pes(coords, mass, shift, calc_potential_point, pes)

  if (proc_id == 0) then
    call print_pes(pes, 'pes_out.txt')
  end if
  call MPI_Finalize(ierr)
end program
