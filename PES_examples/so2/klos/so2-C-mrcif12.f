C
C Ground state (2^1A') C PES of SO2
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
C       Data file:C_so2_bin.dat
c      Please add the following statement in your main program:
c           call potreadc
c     use the following statement to calculate the PES when needed:
c           call so2cpes(r1,r2,cth,v)
c            th is the angle in degree between vector r1 and r2, cth=cos(th)
c       v is the potential in eV
    
       IMPLICIT REAL*8(A-H,O-Z)
        dimension vpes(180)
        pi=dacos(-1.d0)

        call potreadC

        r10=2.81758d0
        r20=3.09726d0
        cthi=dcos(105.5d0/180.d0*pi)
        call so2Cpes(r10,r20,cthi,va)
 110  FORMAT(F11.8,F13.8,F12.2,F12.6)
        print*,'SO2 C-state (2Aprime) singlet MRCISD-F12 PES'
        print*,'r(S-O1) =',r10,'bohr'
        print*,'r(S-O2) =',r20,'bohr'
        print*,'<O1-S-O2=',104.31d0,'degree'
        print*,'Calculated PES value=',va,'eV'
        print*
        print*,'A reference value:5.41173882375606 eV'
C Minimum
Cr(S-O1) =   2.82770000000000      bohr
C r(S-O2) =   3.09950000000000      bohr
C <O1-S-O2=   104.310000000000      degree
C   5.40959338653431      eV

        END

 
	subroutine potreadc
c   we use Jacobi coordinates instead of valence coordinates
c   here the coordinate denote(theta,rHH,rHBr)=(theta,r,R)
c   (nth,nso1,nso2)=(nth,nr,nrr)
        implicit real*8(a-h,o-z)
        parameter (nso1=71,nso2=71,nth=50)
        parameter (m=100,n=100,l=100)
        dimension xt(m)
       dimension y(m),ss(m),sss(m),y2(l)
        common /bvcut/vcut
      common /pes/pesmin,thth(nth,nso2,nso1),rso2(nso2,nso1),rso1(nso1),
     &ve(nth,nso2,nso1),vs(nth,nso2,nso1),ind(nso2,nso1),ind2(nso1)
       data dy1,dyn/1.0d30,1.0d30/
        open(98,file='C_so2_bin.dat')
        vcut=40.d0
	pi=dacos(-1d0)
        pesmin=1000d0
        do i=1,nso1
        ind2(i)=0
          do j=1,nso2
        ind(j,i)=0
            do k=1,nth
            read(98,*,end=100)ra1,ra2,th1,ve1
c            if (k.eq.1) then
c            rra1=ra1
c            rra2=ra2
c            tth1=th1
c            else
c            dal=abs(ra2-rra2)+abs(th1-tth1)
c            if (dal.gt.0.01) then
c            write(*,*)' error  ',ra1,rra1,ra2,rra2,th1,tth1
c            stop
c            endif
c            endif
            ind(j,i)=ind(j,i)+1
            rso2(j,i)=ra2
            rso1(i)=ra1
            thth(k,j,i)=th1
            ve(k,j,i)=ve1
            if (ve(k,j,i).lt.pesmin) then
            pesmin=ve(k,j,i)
            r1m=ra1
            r2m=ra2
            thm=th1
            endif
            if (dabs(th1-180d0).le.0.001d0) goto 11
            enddo
11          continue
C            write(777,*)j,i,ind(j,i)
            ind2(i)=ind2(i)+1
        if (dabs(rso2(j,i)-16.0d0).le.0.001d0) goto 12
            enddo
12          continue
C            write(888,*)i,ind2(i)
            enddo
100         continue
C        write(55,*)r1m,r2m,thm
C        write(*,*)r1m,r2m,thm,pesmin
C        write(*,*)' rso1=',(rso1(i),i=1,nso1)
C        write(*,*)' rso2=',(rso2(j,1),j=1,ind2(1))
C         do i=1,nso1
C          do j=1,ind2(i)
C            do k=1,ind(j,i)
C            write(55,*)' thth=',thth(k,j,i)
C        enddo;enddo;enddo
       close(98)
         do i=1,nso1
          do j=1,ind2(i)
           do k=1,ind(j,i)
          if (ve(k,j,i).ge.vcut) ve(k,j,i)=vcut
          xt(k)=thth(k,j,i)
          y(k)=ve(k,j,i)
c       write(66,1001)rso1(i),rso2(j,i),thth(k,j,i),ve(k,j,i)
 1001      format (1x,3f12.4,f15.8)
           enddo
       call spline(xt,y,ind(j,i),dy1,dyn,y2)
         do k=1,ind(j,i)
         vs(k,j,i)=y2(k)
         enddo
         enddo
       enddo
       close(66)
       end



        subroutine so2cpes(r10,r20,cthi,v)
        implicit real*8(a-h,o-z)
        parameter (rbohr=0.5291771)
        common /bvcut/vcut
c : r1=r(h-h), r2=r(br-h)
	pi=dacos(-1d0)
        r1=r10
        r2=r20
        cth=cthi
        if (r1.lt.r2) then 
        tmp=r2
        r2=r1
        r1=tmp
        endif
        if(cth.gt.1d0)cth=1d0
        if(cth.lt.-1d0)cth=-1d0
        th=dacos(cth)*180.0d0/pi
        if (r1.lt.2.0d0) r1=2.0d0
        if (r1.gt.7.0d0) r1=7.0d0
        if (r2.lt.2.0d0) r2=2.0d0
        if (r2.gt.7.0d0) r2=7.0d0
c        if (th.lt.80d0) th=80d0
c        if (th.gt.180d0) th=180d0
        CALL SPl3(r1,th,r2,va,1)
        v=va
      end




       subroutine spl3(r1,th,r2,v,iop)
       implicit real*8(a-h,o-z)
        parameter (nso1=71,nso2=71,nth=50)
	parameter (m=100,n=100,l=100)
	dimension xt(m)  
       dimension y(m),ss(m),sss(m),y2(l)
       data dy1,dyn/1.0d30,1.0d30/
       common /bvcut/vcut
      common /pes/pesmin,thth(nth,nso2,nso1),rso2(nso2,nso1),rso1(nso1),
     &ve(nth,nso2,nso1),vs(nth,nso2,nso1),ind(nso2,nso1),ind2(nso1)
                     
       do 20 i=1,nso1
       do 10 j=1,ind2(i)
       do 2 k=1,ind(j,i)
        xt(k)=thth(k,j,i)
        y(k)=ve(k,j,i)
        y2(k)=vs(k,j,i)
   2   continue
        call splint(xt,y,y2,ind(j,i),th,y3)
 22     continue
        ss(j)=y3
   10   continue
       do 5 j=1,ind2(i)
        xt(j)=rso2(j,i)
        y(j)=ss(j)
   5   continue
        call spline(xt,y,ind2(i),dy1,dyn,y2)
        call splint(xt,y,y2,ind2(i),r2,yw2)
 33     continue
       sss(i)=yw2
   20    continue
        do i=1,nso1
        xt(i)=rso1(i)
        y(i)=sss(i)
        enddo
        call spline(xt,y,nso1,dy1,dyn,y2)
        call splint(xt,y,y2,nso1,r1,yw2)
   44   continue
        v=yw2
        if(v.gt.vcut) v=vcut
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
      PARAMETER (NMAX=200)
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


