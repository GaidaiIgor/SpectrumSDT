	subroutine dataread
        use path_resolution_mod
        implicit real*8(a-h,o-z)
        parameter (pi=3.141592653589793d0)    
        parameter (nx=100,ny=100)
        common /pest/pesmin,ydat(nx,ny),xdat(nx),ve(nx,ny),ind(nx)
        character(:), allocatable :: path
        path = resolve_relative_exe_path('../dawes/data/SOcorrdata.dat')
        open(98, file = path, status = 'old')
C        open(66,file='check.txt')
         do i=1,nx
            ind(i)=0
            do j=1,ny
             read(98,*)r1,r2,ve1

           ind(i)=ind(i)+1
           ydat(i,j)=r2
           xdat(i)=r1
           ve(i,j)=ve1

C		write (66,*) xdat(i), ydat (i,j),ve(i,j)
C		stop 
C	  write (66,*) xdat(i), ydat (i,j),ve(i,j),i,j


            enddo



C		stop 


         enddo




C          do i=1,nx

C          do k=1,ind(i)
          
C          write(66,*)xdat(i),ydat(i,k),ve(i,k)
C          enddo
C          enddo
	
	close(98)

	return
       end


        subroutine TwoDPES(r10,r20,va)
        implicit real*8(a-h,o-z)
        parameter (pi=3.141592653589793d0)
	common /ranges/xmin,ymin,xmax,ymax
	common /bvcutt/vcut,vmin

        r1=r10
        r2=r20

        CALL SPl2t(r1,r2,va)

	return
        end

       subroutine spl2t(r1,r2,v)
       implicit real*8(a-h,o-z)
       parameter (nx=100,ny=100)
	parameter (m=100,n=100,l=100)
	dimension xt(m)
       dimension y(m),ss(m),sss(m),y2(l)
       data dy1,dyn/1.0d30,1.0d30/
	common /bvcutt/vcut,vmin
        common /pest/pesmin,ydat(nx,ny),xdat(nx),ve(nx,ny),ind(nx)

        do 20 i=1,nx
	nh=0
       do 2 k=1,ind(i)
	nh=nh+1
	xt(nh)=ydat(i,k)
	y(nh)=ve(i,k)
   2   continue
	call spline_d(xt,y,nh,dy1,dyn,y2)
        call splint_d(xt,y,y2,nh,r2,y3)


        ss(i)=y3
20	continue
	nt=0
       do 5 j=1,nx
	nt=nt+1
	xt(nt)=xdat(j)
        y(nt)=ss(j)
   5   continue
	call spline_d(xt,y,nt,dy1,dyn,y2)
	call splint_d(xt,y,y2,nt,r1,yw2)


        v=yw2

	return
        end
C##################################################################
C# SPLINE ROUTINES
C#            Numerical recipes in fortran
C#            Cambrige University Press
C#            York, 2nd edition, 1992.
C##################################################################
      SUBROUTINE splint_d(xa,ya,y2a,n,x,y)
      implicit double precision  (a-h,o-z)
      DIMENSION xa(n),y2a(n),ya(n)
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.0d0) write(6,*) 'bad xa input in splint_d'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.0d0
      return
      END

C##############################################################################
      SUBROUTINE spline_d(x,y,n,yp1,ypn,y2)
      implicit double precision  (a-h,o-z)
      DIMENSION x(n),y(n),y2(n)
      PARAMETER (NMAX=100)
      DIMENSION u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.0d0
        u(1)=0.0d0
      else
        y2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0d0
        y2(i)=(sig-1.0d0)/p
        u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.0d0
        un=0.0d0
      else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
