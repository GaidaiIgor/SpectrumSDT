!-----------------------------------------------------------------------
!  EquGrid (equgrid.f)
!  Equidistant grid generator for ozone PES
!  Author: Alexander Teplukhin
!-----------------------------------------------------------------------
module general
  use constants
  implicit none
  real*8,allocatable::g1(:),g2(:),g3(:)
  real*8 min1,max1
  real*8 min2,max2
  real*8 min3,max3
  integer n1,n2,n3
end module

program equgrid
  use general
  implicit none
  call input_parameters
  allocate(g1(n1),g2(n2),g3(n3))
  call generate_grid(g1,n1,min1,max1)
  call generate_grid(g2,n2,min2,max2)
  call generate_grid(g3,n3,min3,max3)
  call print_grid(g1,n1,'1')
  call print_grid(g2,n2,'2')
  call print_grid(g3,n3,'3')
contains

!-----------------------------------------------------------------------
!  Input Parameters
!  Number of points and range for each coordinate.
!-----------------------------------------------------------------------
  subroutine input_parameters
    implicit none
    open(1,file='equgrid.config')
    read(1,*)n1,n2,n3
    read(1,*)min1,max1
    read(1,*)min2,max2
    read(1,*)min3,max3
    close(1)
    if(max2.eq.0)max2 = pi/2
    if(max3.eq.0)max3 = 2*pi
  end subroutine

!-----------------------------------------------------------------------
!  Generate Grid
!  Generates grid along one coordinate with required number of points
!  and within required boundaries.
!-----------------------------------------------------------------------
  subroutine generate_equidistant_grid(grid,n,minr,maxr)
    implicit none
    real*8,allocatable::grid(:)
    real*8 minr,maxr
    integer i,n
    do i=1,n
     grid(i) = minr + (maxr-minr)*(i-0.5d0)/n
    enddo
  end subroutine

!-----------------------------------------------------------------------
!  Print Grid
!  Saves grid, Jacobian to the *.dat file
!  Corresponding alpha is saved too
!-----------------------------------------------------------------------
  subroutine print_grid(c,n,name)
    implicit none
    real*8 c(n)
    integer i,n
    character name
    open(1,file='grid'//name//'.dat')
    write(1,*)n,c(2)-c(1)
    do i=1,n
     write(1,'(2F25.12)')c(i),1d0
    enddo
    close(1)
  end subroutine
  
end program