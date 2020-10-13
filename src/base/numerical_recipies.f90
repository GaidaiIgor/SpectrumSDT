module numerical_recipies
  use iso_fortran_env, only: real64
  implicit none

  ABSTRACT INTERFACE
    SUBROUTINE diff_equations_rhs(x,y,dydx)
      import
      REAL(real64), INTENT(IN) :: x
      REAL(real64), DIMENSION(:), INTENT(IN) :: y
      REAL(real64), DIMENSION(:), INTENT(OUT) :: dydx
    END SUBROUTINE
  END INTERFACE

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Solves for a vector u of size N the tridiagonal linear set using a
! serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size N,
! while a and c (off-diagonal elements) are size N − 1.
!-------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE tridag(a,b,c,r,u)
    REAL(real64), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(real64), DIMENSION(:), INTENT(OUT) :: u
    REAL(real64), DIMENSION(size(b)) :: gam ! One vector of workspace, gam is needed.
    INTEGER :: n,j
    REAL(real64) :: bet

    if (any(size(a) + 1 /= [size(b), size(c) + 1, size(r), size(u)])) then
      stop 'Size mismatch in tridag'
    end if
    n = size(a) + 1
    bet=b(1)
    if (bet == 0) then
      ! If this happens then you should rewrite your equations as a set of order N − 1, with u2 trivially eliminated.
      stop 'tridag_ser: Error at code stage 1'
    end if
    u(1)=r(1)/bet
    do j=2,n ! Decomposition and forward substitution.
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j-1)*gam(j)
      if (bet == 0) then
        stop 'tridag_ser: Error at code stage 2'
      end if
      u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1 ! Backsubstitution.
      u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridag

!-------------------------------------------------------------------------------------------------------------------------------------------
! Given arrays x and y of length N containing a tabulated function, i.e., yi = f(xi), with x1 <
! x2 < ... < xN , and given values yp1 and ypn for the first derivative of the interpolating
! function at points 1 and N, respectively, this routine returns an array y2 of length N
! that contains the second derivatives of the interpolating function at the tabulated points
! xi. If yp1 and/or ypn are equal to 10^30 or larger, the routine is signaled to set the
! corresponding boundary condition for a natural spline, with zero second derivative on that
! boundary.
!-------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE spline(x,y,yp1,ypn,y2)
    REAL(real64), DIMENSION(:), INTENT(IN) :: x,y
    REAL(real64), INTENT(IN) :: yp1,ypn
    REAL(real64), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER :: n
    REAL(real64), DIMENSION(size(x)) :: a,b,c,r

    if (size(x) /= size(y) .or. size(x) /= size(y2)) then
      stop 'Size of x, y and y2 has to be the same in spline'
    end if
    n = size(x)
    c(1:n-1)=x(2:n)-x(1:n-1) ! Set up the tridiagonal equations.
    r(1:n-1)=6*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2*(c(2:n-1)+a(2:n-1))
    b(1)=1
    b(n)=1
    if (yp1 > 0.99d30) then ! The lower boundary condition is set either to be "natural"
      r(1)=0
      c(1)=0
    else ! or else to have a specified first derivative.
      r(1)=(3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      c(1)=0.5d0
    end if ! The upper boundary condition is set either to be "natural"
    if (ypn > 0.99d30) then
      r(n)=0
      a(n)=0
    else ! or else to have a specified first derivative.
      r(n)=(-3/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
      a(n)=0.5d0
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
  END SUBROUTINE spline

!-------------------------------------------------------------------------------------------------------------------------------------------
! Given an array xx(1:N), and given a value x, returns a value j such that x is between
! xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
! j = N is returned to indicate that x is out of range.
!-------------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION locate(xx,x)
    REAL(real64), DIMENSION(:), INTENT(IN) :: xx
    REAL(real64), INTENT(IN) :: x
    INTEGER :: locate
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd

    n=size(xx)
    ascnd = (xx(n) >= xx(1)) ! True if ascending order of table, false otherwise.
    jl=0 ! Initialize lower
    ju=n+1 ! and upper limits.
    do
      if (ju-jl <= 1) then
        exit
      end if
      jm=(ju+jl)/2 ! Compute a midpoint,
      if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm ! and replace either the lower limit
      else
        ju=jm ! or the upper limit, as appropriate.
      end if
    end do
    if (x == xx(1)) then ! set the output, being careful with the endpoints.
      locate=1
    else if (x == xx(n)) then
      locate=n-1
    else
      locate=jl
    end if
  END FUNCTION locate

!-------------------------------------------------------------------------------------------------------------------------------------------
! Given the arrays xa and ya, which tabulate a function (with the xa_i’s in increasing or
! decreasing order), and given the array y2a, which is the output from spline, and
! given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
! and y2a are all of the same size.
!-------------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION splint(xa,ya,y2a,x)
    REAL(real64), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(real64), INTENT(IN) :: x
    REAL(real64) :: splint
    INTEGER :: khi,klo,n
    REAL(real64) :: a,b,h

    if (size(xa) /= size(ya) .or. size(xa) /= size(y2a)) then
      stop 'Sizes of xa, ya and y2a have to be equal in splint'
    end if
    n=size(xa)
    ! We will find the right place in the table by means of locate’s bisection algorithm. This is
    ! optimal if sequential calls to this routine are at random values of x. If sequential calls are in
    ! order, and closely spaced, one would do better to store previous values of klo and khi and
    ! test if they remain appropriate on the next call.
    klo=max(min(locate(xa,x),n-1),1)
    khi=klo+1 ! klo and khi now bracket the input value of x.
    h=xa(khi)-xa(klo)
    if (h == 0) then
      stop "Bad xa input in splint. The xa's must be distinct"
    end if
    a=(xa(khi)-x)/h ! Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6
  END FUNCTION splint

!-------------------------------------------------------------------------------------------------------------------------------------------
! Given values for the N variables y and their derivatives dydx known at x, use the fourth order Runge-Kutta method to advance the solution
! over an interval h and return the incremented variables as yout, which need not be a distinct array from y. y, dydx and yout
! are all of length N. The user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
!-------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE rk4(y,dydx,x,h,yout,derivs)
    REAL(real64), DIMENSION(:), INTENT(IN) :: y,dydx
    REAL(real64), INTENT(IN) :: x,h
    REAL(real64), DIMENSION(:), INTENT(OUT) :: yout
    PROCEDURE(diff_equations_rhs) :: derivs
    INTEGER :: ndum
    REAL(real64) :: h6,hh,xh
    REAL(real64), DIMENSION(size(y)) :: dym,dyt,yt

    if (size(y) /= size(dydx) .or. size(y) /= size(yout)) then
      stop 'Sizes of y, dydx and yout have to be the same in rk4'
    end if
    ndum=size(y)
    hh=h*0.5d0
    h6=h/6
    xh=x+hh
    yt=y+hh*dydx ! First step.
    call derivs(xh,yt,dyt) ! Second step.
    yt=y+hh*dyt
    call derivs(xh,yt,dym) ! Third step.
    yt=y+h*dym
    dym=dyt+dym
    call derivs(x+h,yt,dyt) ! Fourth step.
    yout=y+h6*(dydx+dyt+2*dym) ! Accumulate increments with proper weights.
  END SUBROUTINE rk4

!-------------------------------------------------------------------------------------------------------------------------------------------
! Starting from N initial values vstart known at x1, use fourth-order Runge-Kutta to advance nstep equal increments to x2. 
! The user-supplied subroutine derivs(x,y,dydx) evaluates derivatives. 
! The original numeric recipes version was modified to not store intermediate results in module "path" variables xx and y.
! Additionally SUBROUTINE was changed to FUNCTION that returns v (from last step). 
!-------------------------------------------------------------------------------------------------------------------------------------------
  FUNCTION rkdumb(vstart,x1,x2,nstep,derivs) result(v)
    REAL(real64), DIMENSION(:), INTENT(IN) :: vstart
    REAL(real64), INTENT(IN) :: x1,x2
    INTEGER, INTENT(IN) :: nstep
    PROCEDURE(diff_equations_rhs) :: derivs
    INTEGER :: k
    REAL(real64) :: h,x
    REAL(real64), DIMENSION(size(vstart)) :: dv,v

    v(:)=vstart(:) ! Load starting values.
    x=x1
    h=(x2-x1)/nstep
    do k=1,nstep ! Take nstep steps.
      call derivs(x,v,dv)
      call rk4(v,dv,x,h,v,derivs)
      if (x+h == x) then
        stop 'stepsize not significant in rkdumb'
      end if
      x=x+h
    end do
  END FUNCTION rkdumb

end module
