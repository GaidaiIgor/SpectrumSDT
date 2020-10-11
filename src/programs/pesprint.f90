!-----------------------------------------------------------------------
!  Prints PES in APH coordinates.
!  Author: Alexander Teplukhin, Igor Gayday
!-----------------------------------------------------------------------
program pesprint
  use parallel_utils
  use pesinterface_mod

  implicit none
  integer :: ierr

  call MPI_Init(ierr)
  call init_pots()
  call load_optgrids()
  call calc_pots(g1, g2, g3)

  if (get_proc_id() == 0) then
    call print_potvib()
    print *, 'Done'
  end if
  call MPI_Finalize(ierr)

contains

!-----------------------------------------------------------------------
!  Loads optimized grids in APH coordinates.
!-----------------------------------------------------------------------
  subroutine load_optgrids()
    integer i
    open(1,file='grid_rho.dat')
    read(1,*)n1
    allocate(g1(n1))
    do i=1,n1
      read(1,*)g1(i)
    enddo
    close(1)
    open(1,file='grid_theta.dat')
    read(1,*)n2
    allocate(g2(n2))
    do i=1,n2
      read(1,*)g2(i)
    enddo
    close(1)
    open(1,file='grid_phi.dat')
    read(1,*)n3
    allocate(g3(n3))
    do i=1,n3
      read(1,*)g3(i)
    enddo
    close(1)
  end subroutine
  
  !-----------------------------------------------------------------------
  !  Prints vibrational potential on 3D grid.
  !-----------------------------------------------------------------------
  subroutine print_potvib()
    integer :: i1, i2, i3, file_unit
    open(newunit = file_unit, file = 'pes.out')
    do i1 = 1, size(potvib, 3)
      do i2 = 1, size(potvib, 2)
        do i3 = 1, size(potvib, 1)
          write(file_unit, '(G23.15)') potvib(i3, i2, i1)
        end do
      end do
    end do
    close(file_unit)
  end subroutine

end program
