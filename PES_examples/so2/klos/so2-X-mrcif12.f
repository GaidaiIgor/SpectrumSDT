C       Ground state adiabatic X(1^1A')  PES of SO2
C       Calculations:MRCISD-F12a+Q/AVTZ/JKFIT/MP2FIT
C       ACTIVE SPACE: 12 active orbitals
C       Author: Jacek Klos
C       Department of Chemistry 
C       University of Maryland
C       College Park, MD 20742
C*******************************************************
C       REFERENCE TO CITE PES USE
C       Jacek Klos, Millard H. Alexander,
C       Praveen Kumar, Bill Poirier,
C       Bin Jiang, and Hua Guo
C       The Journal of Chemical Physics 144, 174301 (2016)
C*******************************************************
C       Data file:X_mrcisdf12.dat
C       Firstly: call potreadX
C       For a given geometry,the potential can be obtained by
C                call so2Xpes(rso1,rso2,cth,v)
C       rso1 and rso2 are two S-O bond lengths (in bohr),
C       RANGE OF SPLINE FOR r1 and r2: 2.0 to 5.0 bohr
C       cth is the cosine of the inter-bond angle O-S-O'
C       angle runs from 0 to 180 degrees
C       v is given in eV
C       Minimum at r1(SO)=r2(SO)=2.7094 bohr and 119.25 degrees

	subroutine potreadX
c   in bondlength coordinates (theta,ros1,ros2)
        use path_resolution_mod
        implicit real*8(a-h,o-z)
        parameter (pi=3.141592653589793d0)    
        parameter (nros1=26,nth=19,nros2=26,ip=1)

!-->    minimum energy of the X state
        data vmin1/0d0/
        data dy1,dyn/1.0d30,1.0d30/
	common /bvcut/vcut
        common /pes/pesmin,ros2(nros1,nth,nros2),ros1(nros1),thth(nth),
     &  ve(ip,nros1,nth,nros2),ind(nros1,nth)
        open(98,file=resolve_relative_exe_path(
     &  '../klos/X_mrcisdf12.dat'))
        vcut=40d0
        pesmin=10000
	kstep=0
        do j=1,nth
         do i=1,nros1
           ind(i,j)=0
            do k=1,nros2
            read(98,*)th1,ra1,ra2,ed
            if (k.eq.1) then
            rra1=ra1
            rra2=ra2
            tth1=th1
            else
            dal=abs(ra1-rra1)+abs(th1-tth1)
            if (dal.gt.0.01) then
            write(*,*)' error  ',ra1,rra1,ra2,rra2
            stop
            endif
            endif
            if (ve1.lt.pesmin) then
            pesmin=ve1
            r1m=ra1
            r2m=ra2
            thm=th1
            endif
           ind(i,j)=ind(i,j)+1
            ros2(i,j,k)=ra2
            ros1(i)=ra1
            thth(j)=th1
            ve(1,i,j,k)=ed-(vmin1-2.322594366831287D-003)
            enddo
          enddo  
        enddo
c        write(55,*)r1m,r2m,thm
c        write(*,*)r1m,r2m,thm,pesmin
c        write(*,*)' th1=',(thth(i),i=1,nth)
c        write(*,*)' ros1=',(ros1(i),i=1,nros1)
          do j=1,nth
          do i=1,nros1
          do k=1,ind(i,j)
          if (ve(1,i,j,k).gt.vcut) ve(1,i,j,k)=vcut
c          write(66,1001)thth(j),ros1(i),ros2(i,j,k),ve(1,i,j,k)
c     &,ve(1,i,j,k),ve(2,i,j,k),ve(3,i,j,k),ve(4,i,j,k)
1001    format(3f10.3,6f15.8)
           enddo
         enddo
       enddo
       end

        subroutine so2Xpes(r10,r20,cthi,va)
        implicit real*8(a-h,o-z)
        parameter (rbohr=0.5291771)
        parameter (pi=3.141592653589793d0)
c        dimension va
	common /bvcut/vcut
c : r1=r(s-o1), r2=r(s-o2),r1 < r2
        r1=r10
        r2=r20
        cth=cthi
        if (r1.gt.r2) then 
        tmp=r2
        r2=r1
        r1=tmp
        endif
	if(cth.gt.1d0) cth=1d0
	if(cth.lt.-1d0) cth=-1d0
        th=dacos(cth)*180.0d0/pi
        if (r1.lt.2.0d0) r1=2.0d0
        if (r1.gt.5.0d0) r1=5.0d0
        if (r2.lt.2.0d0) r2=2.0d0
	if (r2.gt.5.0d0) r2=5.0d0
        if (th.gt.180d0) th=180d0
        if (th.lt.0d0) th=0d0
        CALL spl3d(r2,th,r1,va,1)

1002  format(1x,3f12.4,f20.8)
      end

       subroutine spl3d(r1,th,r2,v,istate)
       implicit real*8(a-h,o-z)
       parameter (nros1=26,nth=19,nros2=26,ip=1)
	parameter (m=1000,n=1000,l=1000)
	dimension xt(m)
       dimension dty(l),ddty(l),s1(l),ds1(l),dds1(l),h1(l)
       dimension dny(n),ddny(n),s2(l),ds2(l),dds2(l),h2(n)
       dimension dhy(m),ddhy(m),s3(l),ds3(l),dds3(l),h3(m)
       dimension y(m),ss(m),sss(m),y2(l)
       data dy1,dyn/1.0d30,1.0d30/
        common /bvcut/vcut
        common /pes/pesmin,ros2(nros1,nth,nros2),ros1(nros1),thth(nth),
     &  ve(ip,nros1,nth,nros2),ind(nros1,nth)
                     
c        do iop=istate,istate
        iop=istate
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE (i,j,nh,k,xt,y,y2,y3,ss,nthth,yw2)
C$OMP DO
        do 20 i=1,nros1
       do 10 j=1,nth
	nh=0
       do 2 k=1,ind(i,j)
	nh=nh+1
	xt(nh)=ros2(i,j,k)
	y(nh)=ve(iop,i,j,k)
   2   continue
  	r1a=r1
	call spline(xt,y,nh,dy1,dyn,y2)
        call splint(xt,y,y2,nh,r1a,y3)
 22	continue
	if (y3.gt.vcut) y3=vcut
        ss(j)=y3
   10   continue
	nthth=0
       do 5 j=1,nth
	nthth=nthth+1
	xt(nthth)=thth(j)
        y(nthth)=ss(j)
   5   continue
	call spline(xt,y,nthth,0.0d0,0.0d0,y2)
	call splint(xt,y,y2,nthth,th,yw2)
 33	continue
	if(yw2.gt.vcut) yw2=vcut
       sss(i)=yw2
   20    continue
C$OMP END DO NOWAIT
C$OMP END PARALLEL
	nr1=0
	nr1=nros1
c	write(*,*)' vrhh=',(sss(i),i=1,nr1)
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP&  PRIVATE (i)
C$OMP DO
	do i=1,nros1
c	nr1=nr1+1
	xt(i)=ros1(i)
	y(i)=sss(i)
  	enddo
C$OMP END DO NOWAIT
C$OMP END PARALLEL
	call spline(xt,y,nr1,dy1,dyn,y2)
	call splint(xt,y,y2,nr1,r2,yw2)
   44	continue
	if (yw2.gt.vcut) yw2=vcut
        v=yw2
c       enddo
       end


C##################################################################
C# SPLINE ROUTINES
C#            Numerical recipes in fortran
C#            Cambrige University Press
C#            York, 2nd edition, 1992.
C##################################################################
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
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
      if (h.eq.0.0d0) write(6,*) 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.0d0
      return
      END
C##############################################################################
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit double precision  (a-h,o-z)
      DIMENSION x(n),y(n),y2(n)
      PARAMETER (NMAX=1000)
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


