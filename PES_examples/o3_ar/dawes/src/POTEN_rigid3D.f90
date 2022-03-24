!*******************************  A U T O S U R F  *********************************
!===================================================================================
!-----------------------------------------------------------------------------------
!-                                                                                 -
!-        AUTOSURF Package: A set of programs for the automated construction       -
!-              of Potential Energy Surfaces on van der Waals systems              -
!-                                                                                 -
!-----------------------------------------------------------------------------------
!===================================================================================
!***********************************************************************************
!-       "POTEN_rigid3D": SUBROUTINE for ...                -
!-----------------------------------------------------------------------------------
!-        Input files: "input-AUTOSURF-PES.dat" & "PES-file"                       -
!-                                                                                 -
!***********************************************************************************

SUBROUTINE PES_LR(jac3, V, path_input1, path_data1, path_input2, path_data2)
implicit none
 character(*) :: path_input1, path_data1, path_input2, path_data2
 integer :: i,initflag,ncont1
 real*8 :: xi(3),jac3(3)
 real*8 :: V1,V2,V,SS



if(jac3(1)<7.5d0) then
call PES(jac3, V1, path_input1, path_data1)
V=(V1+472208.56572973d0)*349.75d0-0.65d0
!write(6,*)'1',jac3(1),V1,0d0,V
return
endif
if(jac3(1)>10.5d0) then
call PES2(jac3, V1, path_input2, path_data2)
V=(V1+472091.813122131d0)*349.75d0
!write(6,*)'2',jac3(1),0d0,V1,V
return
endif

call PES(jac3, V, path_input1, path_data1)
V1=(V+472208.56572973d0)*349.75d0-0.65d0
call PES2(jac3, V, path_input2, path_data2)
V2=(V+472091.813122131d0)*349.75d0

SS=(1d0-tanh(6d0*(jac3(1)-9d0)))/2d0

V=SS*V1+(1d0-SS)*V2
!write(6,*)'X',jac3(1),V1,V2,V


END SUBROUTINE PES_LR


SUBROUTINE PES(jac3, V, path_input, path_data)  

use dynamic_parameters
!-----------------------------------------------------------------------------------
implicit none
 character(*) :: path_input, path_data
 character (len=300) :: line
 integer :: i,initflag,ncont1
 real*8 :: xi(XDIM),jac3(XDIM)
 real*8,allocatable :: cart3(:)
 real*8 :: v,pii,temp3
 logical :: logica1
 save initflag
 data initflag /1/
!-----------------------------------------------------------------------------------
! Interface blocks
!-----------------------------------------------------------------------------------
 INTERFACE! Energy of largest basis and high-level ab initio
   FUNCTION func_actual(xi)    
     use nrtype
     USE dynamic_parameters
     IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: xi
       REAL(SP) :: func_actual
   END FUNCTION func_actual
 end interface
 INTERFACE! Energy of minimal basis and low-level ab initio
   FUNCTION func_actual_seed(xi)  
     use nrtype
     USE dynamic_parameters
     IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: xi
       REAL(SP) :: func_actual_seed
   END FUNCTION func_actual_seed
 end interface
!-----------------------------------------------------------------------------------

pii=dacos(-1d0)

!***********************************************************************************
!                            INITIALIZATION
!***********************************************************************************
IF(initflag==1)THEN   

  alpha=-1d0! coefficient in R=exp(alpha*r) coordinate
  ! basis for "minimal fit" to high level data, also used to fit low level grid.
  order_1_min=3  
  order_2_min=3
  epss=1d-14
  zz=4
  zz_low=4
  zz4=20
  W_a=0.5d0
  dist_tol=1.25d0
  allocate(jac(XDIM),jac2(XDIM))
  !# FRAGMENTS INFORMATION:
  open(unit=10,file=path_input,status='old',action='read')
  112 read(10,'(A100)')line
  ncont1=INDEX(line,'FRAGMENTS INFORMATION:')
  if (ncont1==0) goto 112
  read(10,*)
  read(10,*)nfold
  read(10,*)flip
  read(10,*)reflect
  read(10,*)natom1! number of atoms in the molecule
  natom2=1
  natom=natom1+natom2
  nbdist=natom*(natom-1)/2
  allocate(ref1(3*natom1),ref2(3*natom2),bdist(nbdist),cart(3*(natom)))
  allocate(symb(natom),mass(natom))
  do i=1,natom
    read(10,*) symb(i)! element labels of all atoms
  enddo
  do i=1,natom
    read(10,*) mass(i)! masses of all atoms
  enddo  
  do i=1,3*natom1
    read(10,*) ref1(i)! Cartesian positions for fragment 1 atoms
  enddo
  !# PES INFORMATION:
  OPEN (UNIT=652, FILE=path_data, FORM='UNFORMATTED', ACCESS='SEQUENTIAL') 
  read(652) count3
  read(652) order_1
  read(652) order_2
  read(652) order_3
  read(652) maxpoints
  read(652) mass
  read(652) rmax
  read(652) rmin
  read(652) Max_E
  read(652) low_grid
  read(652) count_seed
  call basis_size3D(order_1,order_2,basis_1)
  call basis_size3D(order_1-1,order_2-1,basis_2)  
  call basis_size3D(order_1_min,order_2_min,basis_3) 
  allocate(b2(basis_1,count3),b2_lower(basis_2,count3),b2_minimal(basis_3,count3))
  allocate(d(count3),coords(count3,XDIM)) 
  allocate(b2_seed(basis_3,count_seed),d_seed(count_seed),coords_seed(count_seed,XDIM))
  b2=0d0
  b2_lower=0d0
  b2_minimal=0d0
  b2_seed=0d0
  d=0d0
  d_seed=0d0
  coords=0d0
  coords_seed=0d0
  do i=1,count3 
     read(652) b2(:,i)       
  enddo
  do i=1,count3      
     read(652) b2_lower(:,i)       
  enddo
  do i=1,count3      
     read(652) b2_minimal(:,i)       
  enddo  
  do i=1,count3       
     read(652) d(i) 
  enddo
  do i=1,count3        
     read(652) coords(i,:)
  enddo
  if(low_grid>0)then
     read(652) Max_E_seed
     do i=1,count_seed      
        read(652) b2_seed(:,i)       
     enddo
     do i=1,count_seed 
        read(652) d_seed(i)       
     enddo
     do i=1,count_seed        
        read(652) coords_seed(i,:)       
     enddo
  endif
  close(652)
  initflag=2

ENDIF
!***********************************************************************************


allocate(cart3(3*(natom)))
xi=jac3

if(xi(1)<rmin(1))then
  v=Max_E
  !write(*,*) 'coords outside fitted range'
  !write(*,*)  xi(1),rmax(1),rmin(1)
  return
endif

! cos(TH) always from -1 to 1
if(xi(2)>1.d0)then
  xi(2)=2.d0-xi(2)
endif
if(xi(2)<-1.d0)then
  xi(2)=-2.d0-xi(2)
endif
if (flip==1) xi(2)=dabs(xi(2))

xi(3)=xi(3)*dble(nfold)
100 continue
! xi(3) always from -pi to pi
if(xi(3)>180.d0)then
  xi(3)=xi(3)-360.d0
  if (xi(3)>180.d0) goto 100
endif
if(xi(3)<-180.d0)then
  xi(3)=xi(3)+360.d0
  if (xi(3)<-180.d0) goto 100
endif
if (reflect==1) xi(3)=dabs(xi(3))
xi(3)=xi(3)*pii/180.d0

call INT_Cart3D(cart3,xi,natom1,ref1,nfold)
call cart_to_bdist_inter3D(cart3,natom1,dist_tol,dist_flag)
if(dist_flag==1) then
  v=Max_E
  !write(*,*) 'bdist less than distol'
  return
endif

if(low_grid>0)then
  temp3=func_actual_seed(xi)
  if(temp3>Max_E_seed)then
    v=Max_E
    !write(*,*)xi(1),'hit ceiling on low grid'
    return
  endif
endif

temp3=func_actual(xi)
if(temp3>Max_E)then
  temp3=Max_E
  !write(*,*)xi(1),'hit ceiling (func_actual)'
endif
V=temp3
return

END SUBROUTINE PES


SUBROUTINE PES2(jac3, V, path_input, path_data)

use dynamic_parameters_2
!-----------------------------------------------------------------------------------
implicit none
 character(*) :: path_input, path_data
 character (len=300) :: line
 integer :: i,initflag,ncont1
 real*8 :: xi(XDIM),jac3(XDIM)
 real*8,allocatable :: cart3(:)
 real*8 :: v,pii,temp3
 logical :: logica1
 save initflag
 data initflag /1/
!-----------------------------------------------------------------------------------
! Interface blocks
!-----------------------------------------------------------------------------------
 INTERFACE! Energy of largest basis and high-level ab initio
   FUNCTION func_actual_2(xi)
     use nrtype
     USE dynamic_parameters_2
     IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: xi
       REAL(SP) :: func_actual_2
   END FUNCTION func_actual_2
 end interface
 INTERFACE! Energy of minimal basis and low-level ab initio
   FUNCTION func_actual_seed_2(xi)
     use nrtype
     USE dynamic_parameters_2
     IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: xi
       REAL(SP) :: func_actual_seed_2
   END FUNCTION func_actual_seed_2
 end interface
!-----------------------------------------------------------------------------------

pii=dacos(-1d0)

!***********************************************************************************
!                            INITIALIZATION
!***********************************************************************************
IF(initflag==1)THEN

  alpha=-1d0! coefficient in R=exp(alpha*r) coordinate
  ! basis for "minimal fit" to high level data, also used to fit low level grid.
  order_1_min=3
  order_2_min=3
  epss=1d-14
  zz=4
  zz_low=4
  zz4=20
  W_a=0.5d0
  dist_tol=1.25d0
  allocate(jac(XDIM),jac2(XDIM))
  !# FRAGMENTS INFORMATION:
  open(unit=10,file=path_input,status='old',action='read')
  112 read(10,'(A100)')line
  ncont1=INDEX(line,'FRAGMENTS INFORMATION:')
  if (ncont1==0) goto 112
  read(10,*)
  read(10,*)nfold
  read(10,*)flip
  read(10,*)reflect
  read(10,*)natom1! number of atoms in the molecule
  natom2=1
  natom=natom1+natom2
  nbdist=natom*(natom-1)/2
  allocate(ref1(3*natom1),ref2(3*natom2),bdist(nbdist),cart(3*(natom)))
  allocate(symb(natom),mass(natom))
  do i=1,natom
    read(10,*) symb(i)! element labels of all atoms
  enddo
  do i=1,natom
    read(10,*) mass(i)! masses of all atoms
  enddo
  do i=1,3*natom1
    read(10,*) ref1(i)! Cartesian positions for fragment 1 atoms
  enddo
  !# PES INFORMATION:
  OPEN (UNIT=652, FILE=path_data, FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
  read(652) count3
  read(652) order_1
  read(652) order_2
  read(652) order_3
  read(652) maxpoints
  read(652) mass
  read(652) rmax
  read(652) rmin
  read(652) Max_E
  read(652) low_grid
  read(652) count_seed
  call basis_size3D(order_1,order_2,basis_1)
  call basis_size3D(order_1-1,order_2-1,basis_2)
  call basis_size3D(order_1_min,order_2_min,basis_3)
  allocate(b2(basis_1,count3),b2_lower(basis_2,count3),b2_minimal(basis_3,count3))
  allocate(d(count3),coords(count3,XDIM))
  allocate(b2_seed(basis_3,count_seed),d_seed(count_seed),coords_seed(count_seed,XDIM))
  b2=0d0
  b2_lower=0d0
  b2_minimal=0d0
  b2_seed=0d0
  d=0d0
  d_seed=0d0
  coords=0d0
  coords_seed=0d0
  do i=1,count3
     read(652) b2(:,i)
  enddo
  do i=1,count3
     read(652) b2_lower(:,i)
  enddo
  do i=1,count3
     read(652) b2_minimal(:,i)
  enddo
  do i=1,count3
     read(652) d(i)
  enddo
  do i=1,count3
     read(652) coords(i,:)
  enddo
  if(low_grid>0)then
     read(652) Max_E_seed
     do i=1,count_seed
        read(652) b2_seed(:,i)
     enddo
     do i=1,count_seed
        read(652) d_seed(i)
     enddo
     do i=1,count_seed
        read(652) coords_seed(i,:)
     enddo
  endif
  close(652)
  initflag=2

ENDIF
!***********************************************************************************


allocate(cart3(3*(natom)))
xi=jac3

if(xi(1)<rmin(1))then
  v=Max_E
  !write(*,*) 'coords outside fitted range'
  !write(*,*)  xi(1),rmax(1),rmin(1)
  return
endif

! cos(TH) always from -1 to 1
if(xi(2)>1.d0)then
  xi(2)=2.d0-xi(2)
endif
if(xi(2)<-1.d0)then
  xi(2)=-2.d0-xi(2)
endif
if (flip==1) xi(2)=dabs(xi(2))

xi(3)=xi(3)*dble(nfold)
100 continue
! xi(3) always from -pi to pi
if(xi(3)>180.d0)then
  xi(3)=xi(3)-360.d0
  if (xi(3)>180.d0) goto 100
endif
if(xi(3)<-180.d0)then
  xi(3)=xi(3)+360.d0
  if (xi(3)<-180.d0) goto 100
endif
if (reflect==1) xi(3)=dabs(xi(3))
xi(3)=xi(3)*pii/180.d0

call INT_Cart3D(cart3,xi,natom1,ref1,nfold)
call cart_to_bdist_inter3D(cart3,natom1,dist_tol,dist_flag)
if(dist_flag==1) then
  v=Max_E
  !write(*,*) 'bdist less than distol'
  return
endif

if(low_grid>0)then
  temp3=func_actual_seed_2(xi)
  if(temp3>Max_E_seed)then
    v=Max_E
    !write(*,*)xi(1),'hit ceiling on low grid'
    return
  endif
endif

temp3=func_actual_2(xi)
if(temp3>Max_E)then
  temp3=Max_E
  !write(*,*)xi(1),'hit ceiling (func_actual)'
endif
V=temp3
return

END SUBROUTINE PES2




!***********************************************************************************
!***********************************************************************************

!                   S U B R O U T I N E S   &   F U N C T I O N S                           

!***********************************************************************************
!***********************************************************************************


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      c a r t _ t o _ b d i s t _ i n t e r 3 D
! ----------------------------------------------------------------------------------
! Known the Cartesian coordinates for all atoms in the system, the internuclear 
! distance (between atoms in different frags.) is computed.
! The variable "flag" is switched to "1" if the atoms are closer than "dist_tol"
!
! *** Input ***
! X         <-- Cartesian coordinates for all atoms in the system
! natom1    <-- Number of atoms in the molecule
! dist_tol  <-- minimum non-bonded internuclear distance

subroutine cart_to_bdist_inter3D(X,natom1,dist_tol,flag)

 implicit none
 integer :: i,k,flag,natom1
 real*8 :: X(3*(natom1+1))
 real*8 :: summ,dist_tol

 flag=0
 do i=1,natom1
   summ=0d0
   do k=1,3
     summ=summ+(X(3*(i-1)+k)-X(3*natom1+k))**2
   enddo
   if(dsqrt(summ)<dist_tol)then
     flag=1
   endif
 enddo

return
end subroutine cart_to_bdist_inter3D

!***********************************************************************************
! ----------------------------------------------------------------------------------
!      B A S I S _ S I Z E 3 D
! ----------------------------------------------------------------------------------
! Calculate the size of the basis set

! *** Input ***
! order_1   <-- maximum power of R = exp(alpha*r)
! order_2   <-- maximum value of L

subroutine basis_size3D(order_1,order_2,basis)
  
 implicit none
 integer :: order_1,order_2,count,basis,l1,m
   
 ! basis calc
 count=0
  
 do l1=0,order_2
   do m=0,l1
     count=count+1
   enddo
 enddo
 basis=count*(order_1)+1

return
end subroutine basis_size3D

!***********************************************************************************
! ----------------------------------------------------------------------------------
!      I N T _ C a r t 3 D
! ----------------------------------------------------------------------------------
! Known the internal coordinates for a given configuration:
! internal(1) -> R
! internal(2) -> cos(theta)
! internal(3) -> phi
! ... the Cartesian coordinates for all atoms in the system (cart) are calculated

! *** Input ***
! internal  <-- vector containing the internal coordinates
! natom1    <-- Number of atoms in the molecule
! ref1      <-- Cartesian coord. of all the nuclei in the molecule 
!               (the origin of the frame of reference is at the CM of the molecule)
! nfold     <-- n-fold rotational symmetry of the molecule with respect to the Z axis

subroutine INT_Cart3D(cart,internal,natom1,ref1,nfold)

 implicit none
 integer,parameter :: XDIM=3
 integer :: natom1,nfold
 real*8 :: pii,sin_theta
 real*8 :: cart((natom1+1)*3),ref1(natom1*3),internal(XDIM)

 pii=dacos(-1d0)

 ! Cartesian coordinates for the atoms inside the molecule
 cart(1:natom1*3)=ref1

 ! Cartesian coordinates for the extra atom
 sin_theta=dsqrt(1.d0-internal(2)**2)
 cart(natom1*3+1)=internal(1)*sin_theta*dcos(internal(3)/nfold)! only in a reduced phi-range...
 cart(natom1*3+2)=internal(1)*sin_theta*dsin(internal(3)/nfold)! ... if n-fold rot. symm. exist
 cart(natom1*3+3)=internal(1)*internal(2)

return
end subroutine INT_Cart3D


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      d i s t _ m e t r i c 3 D
! ----------------------------------------------------------------------------------
! This subroutine computes the "distance-metric" between two given configurations 
!
! *** Input ***
! jac       <-- internal coordinates of configuration 1
! jac2      <-- internal coordinates of configuration 2
! scale     <-- scaling for R in dist metric (W_a)
!
! *** Output ***
! dist      <-- computed distance-metric


subroutine dist_metric3D(jac,jac2,scale,dist)

 implicit none
 integer,parameter :: XDIM=3
 integer :: i
 real*8 :: jac(XDIM),jac2(XDIM),temp(XDIM)
 real*8 :: scale,dist,pii,x1

 pii=acos(-1d0)

 temp(1)=((jac(1)-jac2(1))*scale)**2! dR**2 x (1/W_a)**2
 temp(2)=(dacos(jac(2))-dacos(jac2(2)))**2! dTH**2
 temp(3)=jac(3)-jac2(3)! dPHI

 if(temp(3)>pii)then
   temp(3)=temp(3)-2d0*pii
 endif
 if(temp(3)<-pii)then
   temp(3)=temp(3)+2d0*pii
 endif

 x1=(1d0-jac(2)**2)*(1d0-jac2(2)**2)! sinTH1**2 x sinTH2**2
 temp(3)=(temp(3)**2)*dsqrt(x1)! dPHI**2 x sqrt(...)
 dist=0d0
 do i=1,XDIM
   dist=dist+temp(i)
 enddo
 dist=dsqrt(dist)

return
end subroutine dist_metric3D


!!! index.f90
!c-----------------------------------------------------------------------

SUBROUTINE indexxy(n,arr,indx)

  IMPLICIT NONE
  integer,parameter :: nstack=50, m=7
  INTEGER ::indx(n),istack(nstack)
  INTEGER :: i,indxt,ir,itemp,j,jstack,k,l,n
  REAL*8 :: arr(n)
  REAL*8 :: a
 
  do j=1,n
     indx(j)=j
  enddo
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M)then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,1,-1
           if(arr(indx(i)).le.a)goto 2
           indx(i+1)=indx(i)
        enddo
        i=0
2       indx(i+1)=indxt
     enddo
     if(jstack.eq.0)return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp
     if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     endif
     if(arr(indx(l)).gt.arr(indx(ir)))then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     endif
     if(arr(indx(l+1)).gt.arr(indx(l)))then
        itemp=indx(l+1)
        indx(l+1)=indx(l)
        indx(l)=itemp
     endif
     i=l+1
     j=ir
     indxt=indx(l)
     a=arr(indxt)
3    continue
     i=i+1
     if(arr(indx(i)).lt.a)goto 3
4    continue
     j=j-1
     if(arr(indx(j)).gt.a)goto 4
     if(j.lt.i)goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3
5    indx(l)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'
     if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     endif
  endif
  goto 1
  return

end subroutine indexxy


!***********************************************************************************
!***********************************************************************************
!                                F U N C T I O N S                             
!***********************************************************************************
!***********************************************************************************


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      f u n c _ a c t u a l (xi)
! ----------------------------------------------------------------------------------
! Energy of largest basis and high-level ab initio

! *** Input ***
! xi        <-- vector containing the internal coordinates

function func_actual(xi)

use nrtype
USE dynamic_parameters
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func_actual
real*8 :: temp,weight,norm,somme,pii
real*8 :: jac3(XDIM),jac4(XDIM)
real*8,allocatable :: ind7(:),PM1(:,:),PD1(:,:)
integer :: i,ip,quitt,l1,jj,R,M,count
integer,allocatable :: ind8(:)

pii=dacos(-1d0)
allocate(ind7(count3),ind8(count3),PM1(0:order_2+1,0:order_2+1),PD1(0:order_2+1,0:order_2+1))

jac3=xi! compute and order "distance-metric" between every geometry and xi
count=0
do ip=1,count3 
  count=count+1
  Jac4=coords(ip,:)
  call dist_metric3D(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=dexp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo
call indexxy(count3,ind7,ind8)
quitt=0! number of expansions included in the interpolation
do ip=1,count3
  if(ind7(ind8(count3))/ind7(ind8(count3+1-ip))>1d11) goto 12
  quitt=quitt+1
enddo

12 jac3=xi
Jac4=jac3
jac4(1)=dexp(alpha*jac4(1))
call LPMN(order_2+1,order_2,order_2,jac4(2),PM1,PD1)

norm=0d0
temp=0d0
do i=1,quitt
  jj=ind8(count3+1-i)
!  if(pot(jj)<E_limit)then
   weight=ind7(ind8(count3+1-i)) 
   norm=norm+weight
   temp=temp+weight*b2(1,jj)
   count=1
   do R=1,order_1
    do L1=0,order_2
     do M=0,L1
       count=count+1
       temp=temp+weight*b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*dcos(dble(M)*jac4(3))
     enddo
    enddo
   enddo
!  endif
enddo

func_actual=temp/norm

return
end function func_actual
function func_actual_2(xi)

use nrtype
USE dynamic_parameters_2
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func_actual_2
real*8 :: temp,weight,norm,somme,pii
real*8 :: jac3(XDIM),jac4(XDIM)
real*8,allocatable :: ind7(:),PM1(:,:),PD1(:,:)
integer :: i,ip,quitt,l1,jj,R,M,count
integer,allocatable :: ind8(:)

pii=dacos(-1d0)
allocate(ind7(count3),ind8(count3),PM1(0:order_2+1,0:order_2+1),PD1(0:order_2+1,0:order_2+1))

jac3=xi! compute and order "distance-metric" between every geometry and xi
count=0
do ip=1,count3
  count=count+1
  Jac4=coords(ip,:)
  call dist_metric3D(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=dexp(-((somme)/d(ip)**2))/(((somme)/d(ip)**2)**(zz/2)+epss)
enddo
call indexxy(count3,ind7,ind8)
quitt=0! number of expansions included in the interpolation
do ip=1,count3
  if(ind7(ind8(count3))/ind7(ind8(count3+1-ip))>1d11) goto 12
  quitt=quitt+1
enddo

12 jac3=xi
Jac4=jac3
jac4(1)=dexp(alpha*jac4(1))
call LPMN(order_2+1,order_2,order_2,jac4(2),PM1,PD1)

norm=0d0
temp=0d0
do i=1,quitt
  jj=ind8(count3+1-i)
!  if(pot(jj)<E_limit)then
   weight=ind7(ind8(count3+1-i))
   norm=norm+weight
   temp=temp+weight*b2(1,jj)
   count=1
   do R=1,order_1
    do L1=0,order_2
     do M=0,L1
       count=count+1
       temp=temp+weight*b2(count,jj)*(jac4(1))**(R)*PM1(M,L1)*dcos(dble(M)*jac4(3))
     enddo
    enddo
   enddo
!  endif
enddo

func_actual_2=temp/norm

return
end function func_actual_2


!***********************************************************************************
! ----------------------------------------------------------------------------------
!      f u n c _ a c t u a l _ s e e d
! ----------------------------------------------------------------------------------
!  Energy of minimal basis and low-level ab initio

! *** Input ***
! xi        <-- vector containing the internal coordinates

function func_actual_seed(xi)

use nrtype
USE dynamic_parameters
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func_actual_seed
real*8 :: temp,weight,norm,somme,pii
real*8 :: jac3(XDIM),jac4(XDIM)
real*8,allocatable :: ind7(:),PM1(:,:),PD1(:,:)
integer :: i,ip,quitt,l1,jj,R,M,count
integer,allocatable :: ind8(:)

pii=dacos(-1d0)
allocate(ind7(count_seed),ind8(count_seed))
allocate(PM1(0:order_2_min+1,0:order_2_min+1),PD1(0:order_2_min+1,0:order_2_min+1))

jac3=xi! compute and order "distance-metric" between every geometry and xi
count=0
do ip=1,count_seed 
  count=count+1
  Jac4=coords_seed(ip,:)
  call dist_metric3D(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=dexp(-((somme)/d_seed(ip)**2))/(((somme)/d_seed(ip)**2)**(zz_low/2)+epss)
enddo
call indexxy(count_seed,ind7,ind8)
quitt=0! number of expansions are included in interpolation
do ip=1,count_seed
   if(ind7(ind8(count_seed))/ind7(ind8(count_seed+1-ip))>1d11) goto 12
   quitt=quitt+1
enddo
12 jac3=xi
Jac4=jac3
jac4(1)=dexp(alpha*jac4(1))
call LPMN(order_2_min+1,order_2_min,order_2_min,jac4(2),PM1,PD1)
norm=0d0
temp=0d0

do i=1,quitt
  jj=ind8(count_seed+1-i)
!  if(pot(jj)<E_limit)then
   weight=ind7(ind8(count_seed+1-i)) 
   norm=norm+weight
   temp=temp+weight*b2_seed(1,jj)
   count=1
   do R=1,order_1_min
    do L1=0,order_2_min
     do M=0,L1
       count=count+1
       temp=temp+weight*b2_seed(count,jj)*(jac4(1))**(R)*PM1(M,L1)*dcos(dble(M)*jac4(3))
     enddo
    enddo
   enddo
!  endif
enddo

func_actual_seed=temp/norm

return
end function func_actual_seed


function func_actual_seed_2(xi)

use nrtype
USE dynamic_parameters_2
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: xi
REAL(SP) :: func_actual_seed_2
real*8 :: temp,weight,norm,somme,pii
real*8 :: jac3(XDIM),jac4(XDIM)
real*8,allocatable :: ind7(:),PM1(:,:),PD1(:,:)
integer :: i,ip,quitt,l1,jj,R,M,count
integer,allocatable :: ind8(:)

pii=dacos(-1d0)
allocate(ind7(count_seed),ind8(count_seed))
allocate(PM1(0:order_2_min+1,0:order_2_min+1),PD1(0:order_2_min+1,0:order_2_min+1))

jac3=xi! compute and order "distance-metric" between every geometry and xi
count=0
do ip=1,count_seed
  count=count+1
  Jac4=coords_seed(ip,:)
  call dist_metric3D(jac3,jac4,W_a,somme)
  somme=somme**2
  ind7(count)=dexp(-((somme)/d_seed(ip)**2))/(((somme)/d_seed(ip)**2)**(zz_low/2)+epss)
enddo
call indexxy(count_seed,ind7,ind8)
quitt=0! number of expansions are included in interpolation
do ip=1,count_seed
   if(ind7(ind8(count_seed))/ind7(ind8(count_seed+1-ip))>1d11) goto 12
   quitt=quitt+1
enddo
12 jac3=xi
Jac4=jac3
jac4(1)=dexp(alpha*jac4(1))
call LPMN(order_2_min+1,order_2_min,order_2_min,jac4(2),PM1,PD1)
norm=0d0
temp=0d0

do i=1,quitt
  jj=ind8(count_seed+1-i)
!  if(pot(jj)<E_limit)then
   weight=ind7(ind8(count_seed+1-i))
   norm=norm+weight
   temp=temp+weight*b2_seed(1,jj)
   count=1
   do R=1,order_1_min
    do L1=0,order_2_min
     do M=0,L1
       count=count+1
       temp=temp+weight*b2_seed(count,jj)*(jac4(1))**(R)*PM1(M,L1)*dcos(dble(M)*jac4(3))
     enddo
    enddo
   enddo
!  endif
enddo

func_actual_seed_2=temp/norm

return
end function func_actual_seed_2


