 module splin
 use path_resolution_mod
implicit none

! Global variables within this module 
integer, parameter :: nr=16      ! base points for interpolation
integer, parameter :: nlt=45  ! nombre de  coefficients ds tabCn-roo.res
double precision    coef(nlt,nr)
integer, dimension (nlt,6) ::  ll(nlt,6)
double precision xmin, xmax     ! given interval of x()
!!double precision, dimension (nr) :: xi(nr), yi(nr), bb(nr), cc(nr), dd(nr)
!!double precision, dimension (nr) ::  bq(nr), cq(nr), dq(nr)
double precision, dimension (nr) :: rro(nr), qq2(nr)
!PUBLIC :: nr, coef,  xmin,  xmax, xi, yi, b, c, d

   INTERFACE

   END INTERFACE

CONTAINS
!* This module has the following INTERANL FUNCTIONS:
! ispline
!* This module has the following INERNAL SUBROUTINES:
! inpdat, initspline,spline
!

  subroutine inpdat
  !Initialization of the input data
  implicit none
  !integer, parameter :: nr=16      ! base points for interpolation
  integer i, nl,n
  integer :: file_unit
  !double precision, dimension (nr) :: rro(nr), qq2(nr)
  double precision yy,zz
  character(len = :), allocatable :: path

  path = resolve_relative_exe_path('../dawes/data/tabCn-roo.res')
! step 1: read xi and yi from tabCn-roo.res, xmin, xmax, n
  open (newunit = file_unit, file = path, status = 'old')
  read(file_unit,*)
  do nl=1,nlt
  read(file_unit,*)(ll(nl,i),i=1,6),(coef(nl,n),n=1,nr)
  enddo
  close(file_unit)
  path = resolve_relative_exe_path('../dawes/data/PolarO2_D2h_aVQZ_roo.res')
  open(newunit = file_unit, file = path, status = 'old')
  read(file_unit,*)
  read(file_unit,*)
  do n=1,nr
  read(file_unit,*)rro(n),yy,zz,qq2(n)
  enddo
  close(file_unit)
  end subroutine inpdat
!
 subroutine initspline(np,xi,yi,b,c,d)
! Initialization of the spline
implicit none
!integer, parameter :: np=2  ! interpolation pour le npème coefficient de tabCn-roo.res
double precision, dimension (nr) :: xi(nr), yi(nr), b(nr), c(nr), d(nr)
double precision step
integer i, np

xmin = 1.8
xmax = 3.3

  step = (xmax-xmin)/(nr-1)
  do i=1,nr
  xi(i) = xmin + step*float(i-1) 
  yi(i) = coef(np,i)
!  write (*,200) xi(i), yi(i)
  end do
!  step 2: call spline to calculate spline coeficients
call spline_d (xi, yi, b, c, d, nr) 
200 format (3f12.5)
 end subroutine initspline
! 
   subroutine spline_d (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
double precision x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
double precision h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions 
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination 
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline_d

  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
!
 subroutine initspline2(xi,yi,b,c,d)
! Initialization of the spline
implicit none
!integer, parameter :: nr=16      ! base points for interpolation
double precision, dimension (nr) :: xi(nr), yi(nr), b(nr), c(nr), d(nr)
!double precision, dimension (nr) :: rro(nr),  qq2(nr)
double precision step, xmin, xmax
integer i

xmin = 1.8
xmax = 3.3

  step = (xmax-xmin)/(nr-1)
  do i=1,nr
  xi(i) = xmin + step*float(i-1) 
  if(dabs(xi(i)-rro(i)).ge.1.e-01)then
  print *,'# Erreur de distance O-O',xi(i),rro(i)
  stop
  endif
  yi(i) = qq2(i)
!  write (*,200) xi(i), yi(i)
  end do
!  step 2: call spline to calculate spline coeficients
call spline_d (xi, yi, b, c, d, nr) 
200 format (3f12.5)
 end subroutine initspline2
end module splin
