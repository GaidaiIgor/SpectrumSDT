program so2_pes
  use iso_fortran_env, only: real64
  use mpi
  use pes_utils_mod
  implicit none

  integer :: ierr, proc_id
  real(real64), allocatable :: pes(:)
  real(real64), allocatable :: coords(:, :)
  external :: potreadX

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)

  coords = load_coords('pes_in.txt')
  call potreadX()
  call calc_pes(coords, calc_potential_point, pes = pes)

  if (proc_id == 0) then
    call print_pes(pes, 'pes_out.txt')
  end if
  call MPI_Finalize(ierr)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential of Klos et al. at a configuration given by internal coordinates:
! bond O1-S (Bohr), bond S-O2 (Bohr), O1-S-O2 angle (rad).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point(coords, mass, shift) result(potential)
    real(real64), intent(in) :: coords(3)
    real(real64), optional, intent(in) :: mass(3)
    real(real64), optional, intent(in) :: shift
    real(real64) :: potential
    real(real64), parameter :: eh_per_ev = 0.036749323814708d0
    external :: so2Xpes

    call so2Xpes(coords(1), coords(2), cos(coords(3)), potential)
    potential = potential * eh_per_ev
  end function

end program
