program water_pes
  use iso_fortran_env, only: real64
  use mpi
  use pes_utils_mod
  implicit none

  integer :: ierr, proc_id
  real(real64), allocatable :: coords(:, :)
  real(real64), allocatable :: pes(:)
  character(100) :: arg_buf
  character(:), allocatable :: pes_name
  procedure(pes_point_calculator), pointer :: pes_procedure => null()

  call get_command_argument(1, arg_buf)
  pes_name = trim(adjustl(arg_buf))
  if (pes_name == 'shirin') then
    pes_procedure => calc_potential_point_shirin
  else if (pes_name == 'partridge') then
    pes_procedure => calc_potential_point_partridge
  else
    stop 'Unknown PES name'
  end if

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)

  coords = load_coords('pes_in.txt')
  call calc_pes(coords, pes_procedure, pes = pes)

  if (proc_id == 0) then
    call print_pes(pes, 'pes_out.txt')
  end if
  call MPI_Finalize(ierr)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential of Partridge et al. at a configuration given by internal coordinates:
! bond H1-O (Bohr), bond O-H2 (Bohr), H1-O-H2 angle (rad).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point_partridge(coords, mass, shift) result(potential)
    real(real64), intent(in) :: coords(3)
    real(real64), optional, intent(in) :: mass(3)
    real(real64), optional, intent(in) :: shift
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
  function calc_potential_point_shirin(coords, mass, shift) result(potential)
    real(real64), intent(in) :: coords(3)
    real(real64), optional, intent(in) :: mass(3)
    real(real64), optional, intent(in) :: shift
    real(real64) :: potential
    real(real64) :: target_coords(3)
    external :: potv ! underlying PES procedure

    target_coords(1) = coords(1)
    target_coords(2) = coords(2)
    target_coords(3) = cos(coords(3))
    call potv(potential, target_coords(1), target_coords(2), target_coords(3))
  end function

end program
