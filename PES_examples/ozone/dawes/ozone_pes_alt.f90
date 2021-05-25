!-------------------------------------------------------------------------------------------------------------------------------------------
! This version reads pes_in.txt instead of plain grid files.
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
  real(real64), allocatable :: pes(:)
  real(real64), allocatable :: coords(:, :)

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
  call init_parameters(mass, shift)

  coords = load_coords('pes_in.txt')
  call calc_pes(coords, mass, shift, calc_potential_point, pes)

  if (proc_id == 0) then
    call print_pes(pes, 'pes_out.txt')
  end if
  call MPI_Finalize(ierr)
end program
