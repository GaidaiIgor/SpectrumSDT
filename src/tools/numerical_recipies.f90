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
    if (bet == 0d0) then
      ! If this happens then you should rewrite your equations as a set of order N − 1, with u2 trivially eliminated.
      stop 'tridag_ser: Error at code stage 1'
    end if
    u(1)=r(1)/bet
    do j=2,n ! Decomposition and forward substitution.
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j-1)*gam(j)
      if (bet == 0d0) then
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
    r(1:n-1)=6d0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2d0*(c(2:n-1)+a(2:n-1))
    b(1)=1d0
    b(n)=1d0
    if (yp1 > 0.99d30) then ! The lower boundary condition is set either to be "natural"
      r(1)=0d0
      c(1)=0d0
    else ! or else to have a specified first derivative.
      r(1)=(3d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      c(1)=0.5d0
    end if ! The upper boundary condition is set either to be "natural"
    if (ypn > 0.99d30) then
      r(n)=0d0
      a(n)=0d0
    else ! or else to have a specified first derivative.
      r(n)=(-3d0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
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
    if (h == 0d0) then
      stop "Bad xa input in splint. The xa's must be distinct"
    end if
    a=(xa(khi)-x)/h ! Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6d0
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
    h6=h/6d0
    xh=x+hh
    yt=y+hh*dydx ! First step.
    call derivs(xh,yt,dyt) ! Second step.
    yt=y+hh*dyt
    call derivs(xh,yt,dym) ! Third step.
    yt=y+h*dym
    dym=dyt+dym
    call derivs(x+h,yt,dyt) ! Fourth step.
    yout=y+h6*(dydx+dyt+2d0*dym) ! Accumulate increments with proper weights.
  END SUBROUTINE rk4

! !-----------------------------------------------------------------------
! ! Runge-Kutta 4th order
! ! Given values for the variables y(1:n) and their derivatives dydx(1:n) known at x, use
! ! the fourth-order Runge-Kutta method to advance the solution over an interval h and return
! ! the incremented variables as yout(1:n), which need not be a distinct array from y. The
! ! user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
! !-----------------------------------------------------------------------
!       SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
!       INTEGER n,NMAX
!       real(real64) h,x,dydx(n),y(n),yout(n)
!       EXTERNAL derivs
!       PARAMETER (NMAX=50) ! Set to the maximum number of functions.
!       INTEGER i
!       real(real64) h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
!       hh=h*0.5d0
!       h6=h/6d0
!       xh=x+hh
!       do i=1,n ! First step.
!         yt(i)=y(i)+hh*dydx(i)
!       enddo
!       call derivs(xh,yt,dyt) ! Second step.
!       do i=1,n
!         yt(i)=y(i)+hh*dyt(i)
!       enddo
!       call derivs(xh,yt,dym) ! Third step.
!       do i=1,n
!         yt(i)=y(i)+h*dym(i)
!         dym(i)=dyt(i)+dym(i)
!       enddo
!       call derivs(x+h,yt,dyt) ! Fourth step.
!       do i=1,n ! Accumulate increments with proper weights.
!         yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2d0*dym(i))
!       enddo
!       return
!       END

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

! !-----------------------------------------------------------------------
! ! USES rk4
! ! Starting from initial values vstart(1:nvar) known at x1 use fourth-order Runge-Kutta to
! ! advance nstep equal increments to x2. The user-supplied subroutine derivs(x,v,dvdx)
! ! evaluates derivatives. Results are stored in the common block path. Be sure to dimension
! ! the common block appropriately.
! !-----------------------------------------------------------------------
!       SUBROUTINE rkdumb(vstart,nvar,x1,x2,nstep,derivs,result)
!       INTEGER nstep,nvar,NMAX,NSTPMX
!       PARAMETER (NMAX=50,NSTPMX=100000)
!       ! Maximum number of functions and maximum number of values to be stored.
!       real(real64) x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX),result(nvar)
!       EXTERNAL derivs
!       ! COMMON /path/ xx,y ! Storage of results.
!       INTEGER i,k
!       real(real64) h,x,dv(NMAX),v(NMAX)
!       do i=1,nvar ! Load starting values.
!         v(i)=vstart(i)
!         y(i,1)=v(i)
!       enddo
!       xx(1)=x1
!       x=x1
!       h=(x2-x1)/nstep
!       do k=1,nstep ! Take nstep steps.
!         call derivs(x,v,dv)
!         call rk4(v,dv,nvar,x,h,v,derivs)
!         if (x+h.eq.x) stop 'stepsize not significant in rkdumb'
!         x=x+h
!         xx(k+1)=x ! Store intermediate steps.
!         do i=1,nvar
!           y(i,k+1)=v(i)
!         enddo
!       enddo
!       result = y(:nvar, nstep+1)
!       return
!       END
      
! !-----------------------------------------------------------------------
! !-----------------------------------------------------------------------
!        SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
!        implicit none
!        INTEGER n,NMAX
!        real(real64) eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
!        EXTERNAL derivs
!        PARAMETER (NMAX=50)
!        INTEGER i
!        real(real64) errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,ERRCON
!        PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89e-4)
!
!        h=htry
!
!  1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
!
!        errmax=0.d0
!
!        do 11 i=1,n
!          errmax=max(errmax,abs(yerr(i)/yscal(i)))
!  11    continue
!
!        errmax=errmax/eps
!
!        if(errmax.gt.1.d0)then
!          htemp=SAFETY*h*(errmax**PSHRNK)
!          h=sign(max(abs(htemp),0.1d0*abs(h)),h)
!          xnew=x+h
!          if(xnew.eq.x) stop 'stepsize underflow in rkqs'
!          goto 1
!        else
!          if(errmax.gt.ERRCON)then
!            hnext=SAFETY*h*(errmax**PGROW)
!          else
!            hnext=5.d0*h
!          endif
!          hdid=h
!          x=x+h
!          do 12 i=1,n
!            y(i)=ytemp(i)
!  12      continue
!
!          return
!        endif
!        END

! !-----------------------------------------------------------------------
! !-----------------------------------------------------------------------
!        SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
!        implicit none
!        INTEGER n,NMAX
!        real(real64) h,x,dydx(n),y(n),yerr(n),yout(n)
!        EXTERNAL derivs
!        PARAMETER (NMAX=50)
!        INTEGER i
!        real(real64) ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
!        PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,B63=575./13824.,B64=44275./110592., &
!        B65=253./4096.,C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)
!
!        do 11 i=1,n
!          ytemp(i)=y(i)+B21*h*dydx(i)
!  11    continue
!
!        call derivs(x+A2*h,ytemp,ak2)
!
!        do 12 i=1,n
!          ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
!  12    continue
!
!        call derivs(x+A3*h,ytemp,ak3)
!
!        do 13 i=1,n
!          ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
!  13    continue
!
!        call derivs(x+A4*h,ytemp,ak4)
!
!        do 14 i=1,n
!          ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
!  14    continue
!
!        call derivs(x+A5*h,ytemp,ak5)
!
!        do 15 i=1,n
!          ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
!  15    continue
!
!        call derivs(x+A6*h,ytemp,ak6)
!
!        do 16 i=1,n
!          yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
!  16    continue
!
!        do 17 i=1,n
!          yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
!
!  17    continue
!        return
!        END

! !-----------------------------------------------------------------------
! !-----------------------------------------------------------------------
!        SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
!        implicit none
!        INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
!        real(real64) eps,h1,hmin,x1,x2,ystart(nvar),TINY
!        EXTERNAL derivs,rkqs
!        PARAMETER (MAXSTP=100000000,NMAX=50,KMAXX=200,TINY=1.e-30)
! ! Runge-Kutta driver with adaptive stepsize control. Integrate the starting values ystart(1:nvar)
! ! from x1 to x2 with accuracy eps, storing intermediate results in the common block /path/.
! ! h1 should be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can
! ! be zero). On output nok and nbad are the number of good and bad (but retried and
! ! fixed) steps taken, and ystart is replaced by values at the end of the integration interval.
! ! derivs is the user-supplied subroutine for calculating the right-hand side derivative, while
! ! rkqs is the name of the stepper routine to be used. /path/ contains its own information
! ! about how often an intermediate value is to be stored.
!        INTEGER i,kmax,kount,nstp
!        real(real64) dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
!        ! COMMON /path/ kmax,kount,dxsav,xp,yp
!        ! User storage for intermediate results. Preset dxsav and kmax.
!        kmax=0
!        dxsav=1.d99
!
!        x=x1
!        h=sign(h1,x2-x1)
!        nok=0
!        nbad=0
!        kount=0
!        do 11 i=1,nvar
!          y(i)=ystart(i)
!  11    continue
!        if (kmax.gt.0) xsav=x-2.d0*dxsav ! Assures storage of first step.
!        do 16 nstp=1,MAXSTP ! Take at most MAXSTP steps.
!          call derivs(x,y,dydx)
!          do 12 i=1,nvar
!          ! Scaling used to monitor accuracy. This general-purpose choice can be modified if need be.
!            yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
!  12      continue
!          if(kmax.gt.0)then
!            if(abs(x-xsav).gt.abs(dxsav)) then ! Store intermediate results.
!              if(kount.lt.kmax-1)then
!                kount=kount+1
!                xp(kount)=x
!                do 13 i=1,nvar
!                  yp(i,kount)=y(i)
!  13            continue
!                xsav=x
!              endif
!            endif
!          endif
!          if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x ! If stepsize can overshoot, decrease.
!          call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
!          if(hdid.eq.h)then
!            nok=nok+1
!          else
!            nbad=nbad+1
!          endif
!          if((x-x2)*(x2-x1).ge.0.d0)then ! Are we done?
!            do 14 i=1,nvar
!              ystart(i)=y(i)
!  14        continue
!            if(kmax.ne.0)then
!              kount=kount+1 ! Save final step
!              xp(kount)=x
!              do 15 i=1,nvar
!                yp(i,kount)=y(i)
!  15          continue
!            endif
!            return ! Normal exit.
!          endif
!          if(abs(hnext).lt.hmin) stop 'stepsize smaller than minimum in odeint'
!          h=hnext
!  16    continue
!        stop 'too many steps in odeint'
!        return
!        END
end module
