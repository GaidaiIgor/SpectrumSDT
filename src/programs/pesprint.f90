!-----------------------------------------------------------------------
!  Prints PES in APH coordinates.
!  Author: Alexander Teplukhin, Igor Gayday
!-----------------------------------------------------------------------
program pesprint
  use config
  use input_params_mod
  use parallel_utils
  use pesgeneral
  use pesinterface_mod

  implicit none
  integer :: ierr, file_unit
  type(input_params) :: params

  call init_parameters_pesprint(params)
  call load_optgrids()
  call calc_pots(g1, g2, g3)

  open(newunit = file_unit, file = 'potvib.dat', form = 'unformatted')
  write(file_unit) potvib
  close(file_unit)
  
  if (params % print_potential == 1 .and. get_proc_id() == 0) then
    call print_potvib()
  end if
  call MPI_Finalize(ierr)

contains

!-----------------------------------------------------------------------
!  Loads parameters.
!-----------------------------------------------------------------------
  subroutine init_parameters_pesprint(params)
    type(input_params), intent(inout) :: params
    logical :: config_exists

    inquire(file = 'spectrumsdt.config', exist = config_exists)
    call assert(config_exists, 'Error: spectrumsdt.config does not exist')
    params = process_user_settings('spectrumsdt.config')
    call init_pots(params)
  end subroutine

!-----------------------------------------------------------------------
!  Loads optimized grids in APH coordinates.
!-----------------------------------------------------------------------
  subroutine load_optgrids
    integer i
    open(1,file='grid1.dat')
    read(1,*)n1
    allocate(g1(n1))
    do i=1,n1
      read(1,*)g1(i)
    enddo
    close(1)
    open(1,file='grid2.dat')
    read(1,*)n2
    allocate(g2(n2))
    do i=1,n2
      read(1,*)g2(i)
    enddo
    close(1)
    open(1,file='grid3.dat')
    read(1,*)n3
    allocate(g3(n3))
    do i=1,n3
      read(1,*)g3(i)
    enddo
    close(1)
  end subroutine
  
end program
