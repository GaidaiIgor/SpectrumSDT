module numerical_recipies
contains

       SUBROUTINE splint(xa,ya,yp1,y2a,n,x,y)
       INTEGER n
       REAL*8 x,y,xa(n),y2a(n),ya(n),yp1
       INTEGER k,khi,klo
       REAL*8 a,b,h
       klo=1
       khi=n
 1     if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
           khi=k
         else
           klo=k
         endif
       goto 1
       endif
       h=xa(khi)-xa(klo)
       if (h.eq.0.d0) stop 658
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
       return
       END
       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
       SUBROUTINE spline(x,y,n,yp1,ypn,y2)
       INTEGER n,NMAX
       REAL*8 yp1,ypn,x(n),y(n),y2(n)
       PARAMETER (NMAX=500)
       INTEGER i,k
       REAL*8 p,qn,sig,un,u(NMAX)
       if (yp1.gt..99e30) then
         y2(1)=0.d0
         u(1)=0.d0
       else
         y2(1)=-0.5d0
         u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       endif
       do 11 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.d0)/p
         u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
 11    continue
       if (ypn.gt..99e30) then
       qn=0.d0
       un=0.d0
       else
       qn=0.5d0
       un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       endif
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
       do 12 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
 12    continue
       return
       END
       
!-----------------------------------------------------------------------
! USES rk4
! Starting from initial values vstart(1:nvar) known at x1 use fourth-order Runge-Kutta to
! advance nstep equal increments to x2. The user-supplied subroutine derivs(x,v,dvdx)
! evaluates derivatives. Results are stored in the common block path. Be sure to dimension
! the common block appropriately.
!-----------------------------------------------------------------------
      SUBROUTINE rkdumb(vstart,nvar,x1,x2,nstep,derivs,result)
      INTEGER nstep,nvar,NMAX,NSTPMX
      PARAMETER (NMAX=50,NSTPMX=100000) 
      ! Maximum number of functions and maximum number of values to be stored.
      REAL*8 x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX),result(nvar)
      EXTERNAL derivs
      ! COMMON /path/ xx,y ! Storage of results.
      INTEGER i,k
      REAL*8 h,x,dv(NMAX),v(NMAX)
      do i=1,nvar ! Load starting values.
        v(i)=vstart(i)
        y(i,1)=v(i)
      enddo
      xx(1)=x1
      x=x1
      h=(x2-x1)/nstep
      do k=1,nstep ! Take nstep steps.
        call derivs(x,v,dv)
        call rk4(v,dv,nvar,x,h,v,derivs)
        if (x+h.eq.x) stop 'stepsize not significant in rkdumb'
        x=x+h
        xx(k+1)=x ! Store intermediate steps.
        do i=1,nvar
          y(i,k+1)=v(i)
        enddo
      enddo
      result = y(:nvar, nstep+1)
      return
      END
      
!-----------------------------------------------------------------------
! Runge-Kutta 4th order
! Given values for the variables y(1:n) and their derivatives dydx(1:n) known at x, use
! the fourth-order Runge-Kutta method to advance the solution over an interval h and return
! the incremented variables as yout(1:n), which need not be a distinct array from y. The
! user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
!-----------------------------------------------------------------------
      SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
      INTEGER n,NMAX
      REAL*8 h,x,dydx(n),y(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50) ! Set to the maximum number of functions.
      INTEGER i
      REAL*8 h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
      hh=h*0.5d0
      h6=h/6d0
      xh=x+hh
      do i=1,n ! First step.
        yt(i)=y(i)+hh*dydx(i)
      enddo
      call derivs(xh,yt,dyt) ! Second step.
      do i=1,n
        yt(i)=y(i)+hh*dyt(i)
      enddo
      call derivs(xh,yt,dym) ! Third step.
      do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
      enddo
      call derivs(x+h,yt,dyt) ! Fourth step.
      do i=1,n ! Accumulate increments with proper weights.
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2d0*dym(i))
      enddo
      return
      END

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
       SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
       implicit none
       INTEGER n,NMAX
       REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
       EXTERNAL derivs
       PARAMETER (NMAX=50)
       INTEGER i
       REAL*8 errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,ERRCON
       PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89e-4)

       h=htry

 1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)

       errmax=0.d0

       do 11 i=1,n
         errmax=max(errmax,abs(yerr(i)/yscal(i)))
 11    continue

       errmax=errmax/eps

       if(errmax.gt.1.d0)then
         htemp=SAFETY*h*(errmax**PSHRNK)
         h=sign(max(abs(htemp),0.1d0*abs(h)),h)
         xnew=x+h
         if(xnew.eq.x) stop 'stepsize underflow in rkqs'  
         goto 1
       else
         if(errmax.gt.ERRCON)then
           hnext=SAFETY*h*(errmax**PGROW)
         else
           hnext=5.d0*h
         endif
         hdid=h
         x=x+h
         do 12 i=1,n
           y(i)=ytemp(i)
 12      continue

         return
       endif
       END

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
       SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
       implicit none
       INTEGER n,NMAX
       REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
       EXTERNAL derivs
       PARAMETER (NMAX=50)
       INTEGER i
       REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
       PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,B63=575./13824.,B64=44275./110592., &
       B65=253./4096.,C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)

       do 11 i=1,n
         ytemp(i)=y(i)+B21*h*dydx(i)
 11    continue

       call derivs(x+A2*h,ytemp,ak2)

       do 12 i=1,n
         ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
 12    continue

       call derivs(x+A3*h,ytemp,ak3)

       do 13 i=1,n
         ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
 13    continue

       call derivs(x+A4*h,ytemp,ak4)

       do 14 i=1,n
         ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))  
 14    continue

       call derivs(x+A5*h,ytemp,ak5)

       do 15 i=1,n
         ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
 15    continue

       call derivs(x+A6*h,ytemp,ak6)

       do 16 i=1,n
         yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
 16    continue

       do 17 i=1,n
         yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))

 17    continue
       return
       END

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
       SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
       implicit none
       INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
       REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
       EXTERNAL derivs,rkqs
       PARAMETER (MAXSTP=100000000,NMAX=50,KMAXX=200,TINY=1.e-30)
! Runge-Kutta driver with adaptive stepsize control. Integrate the starting values ystart(1:nvar)
! from x1 to x2 with accuracy eps, storing intermediate results in the common block /path/.
! h1 should be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can
! be zero). On output nok and nbad are the number of good and bad (but retried and
! fixed) steps taken, and ystart is replaced by values at the end of the integration interval.
! derivs is the user-supplied subroutine for calculating the right-hand side derivative, while
! rkqs is the name of the stepper routine to be used. /path/ contains its own information
! about how often an intermediate value is to be stored.
       INTEGER i,kmax,kount,nstp
       REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
       ! COMMON /path/ kmax,kount,dxsav,xp,yp
       ! User storage for intermediate results. Preset dxsav and kmax.
       kmax=0
       dxsav=1.d99

       x=x1
       h=sign(h1,x2-x1)
       nok=0
       nbad=0
       kount=0
       do 11 i=1,nvar
         y(i)=ystart(i)
 11    continue
       if (kmax.gt.0) xsav=x-2.d0*dxsav ! Assures storage of first step.
       do 16 nstp=1,MAXSTP ! Take at most MAXSTP steps.
         call derivs(x,y,dydx)
         do 12 i=1,nvar
         ! Scaling used to monitor accuracy. This general-purpose choice can be modified if need be.
           yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
 12      continue
         if(kmax.gt.0)then
           if(abs(x-xsav).gt.abs(dxsav)) then ! Store intermediate results.
             if(kount.lt.kmax-1)then
               kount=kount+1
               xp(kount)=x
               do 13 i=1,nvar
                 yp(i,kount)=y(i)
 13            continue
               xsav=x
             endif
           endif
         endif
         if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x ! If stepsize can overshoot, decrease.
         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
         if(hdid.eq.h)then
           nok=nok+1
         else
           nbad=nbad+1
         endif
         if((x-x2)*(x2-x1).ge.0.d0)then ! Are we done?
           do 14 i=1,nvar
             ystart(i)=y(i)
 14        continue
           if(kmax.ne.0)then
             kount=kount+1 ! Save final step
             xp(kount)=x
             do 15 i=1,nvar
               yp(i,kount)=y(i)
 15          continue
           endif
           return ! Normal exit.
         endif
         if(abs(hnext).lt.hmin) stop 'stepsize smaller than minimum in odeint'  
         h=hnext
 16    continue
       stop 'too many steps in odeint'  
       return
       END
end module
