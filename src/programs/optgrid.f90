!-----------------------------------------------------------------------
!  OptGrid (optgrid.f)
!  Optimal grid generator for ozone PES in coordinates:
!    1. APH
!    2. Valence (bonds and angle)
!  Author: Alexander Teplukhin
!-----------------------------------------------------------------------
module general
  use constants
  implicit none
  real*8, parameter ::  m0 = isomass(1), &
                        m1 = isomass(1), &
                        m2 = isomass(1), &
                        Mtot = m0+m1+m2, &
                        mu = sqrt(m0*m1*m2/Mtot), &
                        mu1 = m0*m1/(m0+m1), &
                        mu2 = m0*m2/(m0+m2), &
                        r10 = 2.410d0, &
                        r20 = 2.410d0, &
                        te0 = pi * 116.8 / 180 ! , &
                        
  real*8 :: zpe
  
  real*8,allocatable::g1(:),g2(:),g3(:)
  real*8,allocatable::jac1(:),jac2(:),jac3(:)
  real*8 a1,a2,a3
  real*8 Emax1,Emax2,Emax3
  real*8 min1,max1
  real*8 min2,max2
  real*8 min3,max3
  real*8 eps
  integer nenv1,nenv2,nenv3
  integer n1,n2,n3
  integer envtype1,envtype2,envtype3
  integer cs,rgrid,symmphi
  !--- Envelopes ---
  character*256 envpath
  real*8 env1_d1,env1_dn
  real*8 env2_d1,env2_dn
  real*8 env3_d1,env3_dn
  real*8,allocatable::env1(:),pot1(:),env1_2d(:)
  real*8,allocatable::env2(:),pot2(:),env2_2d(:)
  real*8,allocatable::env3(:),pot3(:),env3_2d(:)
  real*8 enva1,enva2,enva3
  real*8 envE1,envE2,envE3
  real*8 envm1,envm2,envm3
  real*8 energyoffset
end module

program optgrid
  use general
  use constants
  use optgrid_tools
  implicit none

  call input_parameters
  allocate(env1(nenv1), pot1(nenv1), env1_2d(nenv1), env2(nenv2), pot2(nenv2), env2_2d(nenv2), env3(nenv3), pot3(nenv3), env3_2d(nenv3), g1(n1), g2(n2), g3(n3), jac1(n1), jac2(n2), jac3(n3))

  call input_envelopes
  call find_parabola(env1,pot1,nenv1,enva1,envE1,envm1)
  call find_parabola(env2,pot2,nenv2,enva2,envE2,envm2)
  call find_parabola(env3,pot3,nenv3,enva3,envE3,envm3)
  if(symmphi==2) envm3 = pi

  ! APH coordinates: g1 - rho, g2 - tet, g3 - phi
  if(cs.eq.1)then
    select case(rgrid)
    case(1)
      call generate_grida(g2,jac2,n2,a2,min2,max2,eps,der_tet)
      call generate_grida(g3,jac3,n3,a3,min3,max3,eps,der_phi)
      a1 = (a2+a3)/2
      call generate_gridr(g1,jac1,n1,a1,min1,eps,der_rho)
    case(2)
      call generate_grida(g1,jac1,n1,a1,min1,max1,eps,der_rho)
      call generate_grida(g2,jac2,n2,a2,min2,max2,eps,der_tet)
      call generate_grida(g3,jac3,n3,a3,min3,max3,eps,der_phi)
    case default
      deallocate(g1, g2, g3, jac1, jac2, jac3)
      call generate_gridn(g1,jac1,n1,a1,min1,max1,envm1,eps,der_rho)
      ! call generate_gridn(g2,jac2,n2,a2,min2,max2,envm2,eps,der_tet)

      ! Select step/points mode
      if (n2 == 0) then
        call generate_equidistant_grid_step(min2, max2, a2, g2, jac2, n2)
      else
        call generate_equidistant_grid_points(min2, max2, n2, g2, jac2, a2)
      end if
      call generate_equidistant_grid_points(min3, max3, n3, g3, jac3, a3)
    end select
  ! Valence coordinates: g1 - R1, g2 - R2, g3 - te
  else
    if(rgrid.eq.1)then
      call generate_gridr(g1,jac1,n1,a1,min1,eps,der_r1)
      call generate_gridr(g2,jac2,n2,a2,min2,eps,der_r2)
    else
      call generate_grida(g1,jac1,n1,a1,min1,max1,eps,der_r1)
      call generate_grida(g2,jac2,n2,a2,min2,max2,eps,der_r2)
    endif
    call generate_grida(g3,jac3,n3,a3,min3,max3,eps,der_te)
  endif
  call print_grid(g1,jac1,n1,'grid1.dat',a1,potv1)
  call print_grid(g2,jac2,n2,'grid2.dat',a2,potv2)
  call print_grid(g3,jac3,n3,'grid3.dat',a3,potv3)
contains

  !-----------------------------------------------------------------------
  !  Input Parameters
  !  J vector, number of point along each coordinate, initial alpha
  !-----------------------------------------------------------------------
  subroutine input_parameters
    implicit none
    integer :: m0_code, m1_code, m2_code, total_code
    open(1,file='optgrid.config')
    read(1,'(A)')envpath
    read(1,*)n1,n2,n3
    read(1,*)nenv1,nenv2,nenv3
    read(1,*)a1,a2,a3
    read(1, *) envtype1, envtype2, envtype3
    read(1, *) min1, max1, Emax1
    read(1, *) min2, max2, Emax2
    
    if (envtype3 == 1) then
      read(1, *) min3, max3, Emax3, symmphi
    elseif (envtype3 == 2) then
      read(1, *) min3, max3, n3, symmphi
    end if
    
    read(1,*)cs,eps,rgrid
    read(1, *) m1_code, m0_code, m2_code
    if(cs.eq.2)then
      read(1,*)energyoffset
      min3 = pi * min3 / 180
      max3 = pi * max3 / 180
      energyoffset = energyoffset / autown
    endif
    close(1)
    Emax1 = Emax1 / autown
    Emax2 = Emax2 / autown
    Emax3 = Emax3 / autown
    
    total_code = m0_code + m1_code + m2_code
    if (total_code == 48) then
      zpe = zpe_66
    else if (total_code == 50) then
      zpe = zpe_68
    else if (total_code == 52) then
      zpe = zpe_88
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Input Envelopes
  !  Loads potential along MEP for each coordinate
  !-----------------------------------------------------------------------
  subroutine input_envelopes
    implicit none
    integer i
    open(1,file=trim(envpath)//'/MEP1.dat',status='old')
    do i=1,nenv1
      read(1,*) env1(i), pot1(i)
    enddo
    close(1)
    pot1 = pot1 / autown + zpe
    open(1,file=trim(envpath)//'/MEP2.dat',status='old')
    do i=1,nenv2
      read(1,*) env2(i), pot2(i)
    enddo
    close(1)
    pot2 = pot2/autown + zpe
    open(1,file=trim(envpath)//'/MEP3.dat',status='old')
    do i=1,nenv3
      read(1,*) env3(i), pot3(i)
    enddo
    close(1)
    pot3 = pot3/autown + zpe

    ! Preparation for extrapolation of envelopes
    env1_d1 = (pot1(2)-pot1(1))/(env1(2)-env1(1))
    env1_dn = (pot1(nenv1)-pot1(nenv1-1))/(env1(nenv1)-env1(nenv1-1))
    call spline(env1,pot1,nenv1,env1_d1,env1_dn,env1_2d)
    env2_d1 = (pot2(2)-pot2(1))/(env2(2)-env2(1))
    env2_dn = (pot2(nenv2)-pot2(nenv2-1))/(env2(nenv2)-env2(nenv2-1))
    call spline(env2,pot2,nenv2,env2_d1,env2_dn,env2_2d)
    env3_d1 = (pot3(2)-pot3(1))/(env3(2)-env3(1))
    env3_dn = (pot3(nenv3)-pot3(nenv3-1))/(env3(nenv3)-env3(nenv3-1))
    call spline(env3,pot3,nenv3,env3_d1,env3_dn,env3_2d)
  end subroutine

  !-----------------------------------------------------------------------
  !  Find parabola
  !-----------------------------------------------------------------------
  subroutine find_parabola(env,pot,nenv,enva,envE,envm)
    implicit none
    real*8 env(nenv),pot(nenv),enva,envE,envm
    integer nenv,i,im,il,ir
    real*8 mat(3,3),d(3)
    real*8 det0,det1,c,x
    im = 1
    do i=1,nenv
      if(pot(i).lt.pot(im))im = i
    enddo
    do il=im,1,-1
      if(pot(il)-pot(im).gt.energyoffset)exit
    enddo
    do ir=im,nenv
      if(pot(ir)-pot(im).gt.energyoffset)exit
    enddo
    mat(1,1) = env(il)**2
    mat(1,2) = env(il)
    mat(1,3) = 1
    mat(2,1) = env(im)**2
    mat(2,2) = env(im)
    mat(2,3) = 1
    mat(3,1) = env(ir)**2
    mat(3,2) = env(ir)
    mat(3,3) = 1
    d(1) = pot(il)
    d(2) = pot(im)
    d(3) = pot(ir)
    det0 = det(mat)
    mat(1,1) = d(1)
    mat(2,1) = d(2)
    mat(3,1) = d(3)
    det1 = det(mat)
    c = det1 / det0
    envE = abs(pot(im))
    enva = sqrt(envE/c)
    envm = env(im)
  end subroutine

  real*8 function det(a)
    real*8 a(3,3)
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end function
  
  !-----------------------------------------------------------------------
  !  PotVib
  !  Returns interpolated value of envelope at a given point
  !  Adds approximate ZPE value along two other coordinates
  !-----------------------------------------------------------------------
  real*8 function potv1(r)
    use general
    implicit none
    real*8 r
    if(r>env1(nenv1))then
      potv1 = pot1(nenv1)
    else
      call splint(env1,pot1,env1_d1,env1_2d,nenv1,r,potv1)
    endif
  end function

  real*8 function potv2(r)
    use general
    implicit none
    real*8 r
    if(r>env2(nenv2))then
      potv2 = pot2(nenv2)
    else
      call splint(env2,pot2,env2_d1,env2_2d,nenv2,r,potv2)
    endif
  end function

  real*8 function potv3(r)
    use general
    implicit none
    real*8 r
    if(symmphi==2.and.r<pi)then
      call splint(env3,pot3,env3_d1,env3_2d,nenv3,2*pi-r,potv3)
    else
      call splint(env3,pot3,env3_d1,env3_2d,nenv3,r,potv3)
    endif
  end function

  !-----------------------------------------------------------------------
  !  PotEnv
  !  Smooth envelopes by Eckart function
  !-----------------------------------------------------------------------
  real*8 function potenv1(r)
    use general
    implicit none
    real*8 r,x
    if(envtype1.eq.1)then
      x = r - envm1
      if(x.lt.0)then
        potenv1 = -envE1 * 4 / (exp(x/enva1) + exp(-x/enva1))**2
      else
        if(r>env1(nenv1))then
          potenv1 = pot1(nenv1)
        else
          call splint(env1,pot1,env1_d1,env1_2d,nenv1,r,potenv1)
        endif
      endif
    else
      potenv1 = - envE1
    endif
  end function

  real*8 function potenv2(r)
    use general
    implicit none
    real*8 r,x
    if(envtype2.eq.1)then
      x = r - envm2
      if(x.lt.0)then
        potenv2 = -envE2 * 4 / (exp(x/enva2) + exp(-x/enva2))**2
      else
        if(r>env2(nenv2))then
          potenv2 = pot2(nenv2)
        else
          call splint(env2,pot2,env2_d1,env2_2d,nenv2,r,potenv2)
        endif
      endif
    else
      potenv2 = - envE2
    endif
  end function

  real*8 function potenv3(r)
    use general
    implicit none
    real*8 r,x
    if(envtype3.eq.1)then
      x = r - envm3
      potenv3 = - envE3 * 4 / (exp(x/enva3) + exp(-x/enva3))**2
    else
      potenv3 = - envE3
    endif
  end function

  !-----------------------------------------------------------------------
  !  Grid derivative for rho
  !-----------------------------------------------------------------------
  subroutine der_rho(x,y,dydx)
    use general
    implicit none
    real*8 x,y,dydx
    real*8 argument
    argument = 2*mu*(Emax1-potenv1(y)-gridextra(envm1,envm2))
    if (argument.gt.0.d0) then
      dydx = pi/DSQRT(argument)
    else
      stop 'Potential is greater than Emax1'
    endif
  end subroutine

  !-----------------------------------------------------------------------
  !  Grid derivative for tet
  !-----------------------------------------------------------------------
  subroutine der_tet(x,y,dydx)
    use general
    implicit none
    real*8 x,y,dydx
    real*8 argument
    argument = mu/2*(Emax2-potenv2(y)-gridextra(envm1,envm2))
    if (argument.gt.0.d0) then
      dydx = pi/(envm1*DSQRT(argument))
    else
      stop 'Potential is greater than Emax2'
    endif
  end subroutine

  !-----------------------------------------------------------------------
  !  Grid derivative for phi
  !-----------------------------------------------------------------------
  subroutine der_phi(x,y,dydx)
    use general
    implicit none
    real*8 x,y,dydx
    real*8 argument
    argument = mu/2*(Emax3-potenv3(y)-gridextra(envm1,envm2))
    if (argument.gt.0.d0) then
      dydx = pi/(envm1*sin(envm2)*DSQRT(argument))
    else
      stop 'Potential is greater than Emax3'
    endif
  end subroutine

  !-----------------------------------------------------------------------
  !  Extra potential caused by conversion from Jacoby coordinates
  !-----------------------------------------------------------------------
  real*8 function gridextra(rho,tet)
    use general
    implicit none
    real*8 rho,tet
    gridextra = - 2/(mu*rho**2)*(sin(2*tet)**-2+0.0625d0)
  end function

  !-----------------------------------------------------------------------
  !  Grid derivative for R1
  !-----------------------------------------------------------------------
  subroutine der_r1(x,y,dydx)
    use general
    implicit none
    real*8 x,y,dydx
    real*8 argument
    argument = 2*mu1*(Emax1-potenv1(y))
    if (argument.gt.0.d0) then
      dydx = pi/DSQRT(argument)
    else
      stop 'Potential is greater than Emax1'
    endif
  end subroutine

  !-----------------------------------------------------------------------
  !  Grid derivative for R2
  !-----------------------------------------------------------------------
  subroutine der_r2(x,y,dydx)
    use general
    implicit none
    real*8 x,y,dydx
    real*8 argument
    argument = 2*mu2*(Emax2-potenv2(y))
    if (argument.gt.0.d0) then
      dydx = pi/DSQRT(argument)
    else
      stop 'Potential is greater than Emax2'
    endif
  end subroutine

  !-----------------------------------------------------------------------
  !  Grid derivative for Te
  !-----------------------------------------------------------------------
  subroutine der_te(x,y,dydx)
    use general
    implicit none
    real*8 x,y,dydx
    real*8 argument
    argument = 2/(1/(mu1*envm1**2) + 1/(mu2*envm2**2) - 2*cos(envm3)/(m0*envm1*envm2)) * (Emax3-potenv3(y))
    if (argument.gt.0.d0) then
      dydx = pi/DSQRT(argument)
    else
      stop 'Potential is greater than Emax3'
    endif
  end subroutine
end

