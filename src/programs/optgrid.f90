!-----------------------------------------------------------------------
!  OptGrid (optgrid.f)
!  Optimal grid generator for ozone PES in coordinates:
!    1. APH
!    2. Valence (bonds and angle)
!  Author: Alexander Teplukhin, Igor Gayday
!-----------------------------------------------------------------------
module global_vars
  use constants
  use iso_fortran_env, only: real64

  implicit none
  real(real64), parameter :: m0 = isomass(1), &
                        m1 = isomass(1), &
                        m2 = isomass(1), &
                        Mtot = m0+m1+m2, &
                        mu = sqrt(m0*m1*m2/Mtot), &
                        mu1 = m0*m1/(m0+m1), &
                        mu2 = m0*m2/(m0+m2), &
                        r10 = 2.410d0, &
                        r20 = 2.410d0, &
                        te0 = pi * 116.8 / 180
                        
  real(real64) :: zpe
  real(real64),allocatable :: grid_rho(:), grid_theta(:), grid_phi(:)
  real(real64),allocatable :: jac_rho(:), jac_theta(:), jac_phi(:)
  real(real64) :: alpha_rho, alpha_theta, alpha_phi
  real(real64) :: Emax1,Emax2,Emax3
  real(real64) :: rho_from, rho_to
  real(real64) :: theta_from, theta_to
  real(real64) :: phi_from, phi_to
  real(real64) :: eps
  integer :: nenv1,nenv2,nenv3
  integer :: npoints_rho, npoints_theta, npoints_phi
  integer :: envtype1,envtype2,envtype3
  integer :: rgrid,symmphi
  !--- Envelopes ---
  character(:), allocatable :: envpath
  real(real64) :: env1_d1,env1_dn
  real(real64) :: env2_d1,env2_dn
  real(real64) :: env3_d1,env3_dn
  real(real64),allocatable :: env1(:),pot1(:),env1_2d(:)
  real(real64),allocatable :: env2(:),pot2(:),env2_2d(:)
  real(real64),allocatable :: env3(:),pot3(:),env3_2d(:)
  real(real64) :: enva1,enva2,enva3
  real(real64) :: envE1,envE2,envE3
  real(real64) :: envm1,envm2,envm3
  real(real64) :: energyoffset
end module

program optgrid
  use config_mod
  use constants
  use global_vars
  use input_params_mod
  use iso_fortran_env, only: real64
  use optgrid_tools
  use path_utils

  implicit none
  type(input_params) :: params

  params = process_user_settings('spectrumsdt.config')
  call input_parameters(params)
  call generate_equidistant_grid_points(params % grid_rho_from, params % grid_rho_to, params % grid_rho_npoints, grid_rho, jac_rho, alpha_rho)
  call generate_equidistant_grid_points(params % grid_theta_from, params % grid_theta_to, params % grid_theta_npoints, grid_theta, jac_theta, alpha_theta)
  call generate_equidistant_grid_points(params % grid_phi_from, params % grid_phi_to, params % grid_phi_npoints, grid_phi, jac_phi, alpha_phi)
  call print_grid(grid_rho, jac_rho, alpha_rho, 'grid_rho.dat')
  call print_grid(grid_theta, jac_theta, alpha_theta, 'grid_theta.dat')
  call print_grid(grid_phi, jac_phi, alpha_phi, 'grid_phi.dat')
  print *, 'Done'

  ! allocate(env1(nenv1), pot1(nenv1), env1_2d(nenv1), env2(nenv2), pot2(nenv2), env2_2d(nenv2), env3(nenv3), pot3(nenv3), env3_2d(nenv3), g1(n1), g2(n2), g3(n3), jac1(n1), jac2(n2), jac3(n3))

  ! call input_envelopes()
  ! call find_parabola(env1,pot1,nenv1,enva1,envE1,envm1)
  ! call find_parabola(env2,pot2,nenv2,enva2,envE2,envm2)
  ! call find_parabola(env3,pot3,nenv3,enva3,envE3,envm3)
  ! if (symmphi == 2) envm3 = pi

  ! select case (rgrid)
  ! case(1)
  !   call generate_grida(g2,jac2,n2,a2,min2,max2,eps,der_tet)
  !   call generate_grida(g3,jac3,n3,a3,min3,max3,eps,der_phi)
  !   a1 = (a2+a3)/2
  !   call generate_gridr(g1,jac1,n1,a1,min1,eps,der_rho)
  ! case(2)
  !   call generate_grida(g1,jac1,n1,a1,min1,max1,eps,der_rho)
  !   call generate_grida(g2,jac2,n2,a2,min2,max2,eps,der_tet)
  !   call generate_grida(g3,jac3,n3,a3,min3,max3,eps,der_phi)
  ! case default
  !   deallocate(g1, g2, g3, jac1, jac2, jac3)
  !   call generate_gridn(g1,jac1,n1,a1,min1,max1,envm1,eps,der_rho)
  !   ! call generate_gridn(g2,jac2,n2,a2,min2,max2,envm2,eps,der_tet)
  !
  !   ! Select step/points mode
  !   if (n2 == 0) then
  !     call generate_equidistant_grid_step(min2, max2, a2, g2, jac2, n2)
  !   else
  !     call generate_equidistant_grid_points(min2, max2, n2, g2, jac2, a2)
  !   end if
  !   call generate_equidistant_grid_points(min3, max3, n3, g3, jac3, a3)
  ! end select

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Input parameters.
! J vector, number of point along each coordinate, initial alpha.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine input_parameters(params)
    class(input_params), intent(in) :: params
    integer :: m0_code, m1_code, m2_code, total_code

    rho_from = params % grid_rho_from
    rho_to = params % grid_rho_to
    npoints_rho = params % grid_rho_npoints
    theta_from = params % grid_theta_from
    theta_to = params % grid_theta_to
    npoints_theta = params % grid_theta_npoints
    phi_from = params % grid_phi_from
    phi_to = params % grid_phi_to
    npoints_phi = params % grid_phi_npoints

    Emax1 = 0d0
    Emax2 = 0d0
    Emax3 = 0d0
    nenv1 = 1
    nenv2 = 1
    nenv3 = 1
    alpha_rho = 0d0
    alpha_theta = 0d0
    alpha_phi = 0d0
    envtype1 = 2
    envtype2 = 2
    envtype3 = 2
    symmphi = 2
    eps = 0d0
    rgrid = 3
    m1_code = 16
    m0_code = 16
    m2_code = 16

    total_code = m0_code + m1_code + m2_code
    if (total_code == 48) then
      zpe = zpe_66
    else if (total_code == 50) then
      zpe = zpe_68
    else if (total_code == 52) then
      zpe = zpe_88
    end if
  end subroutine

  ! !-----------------------------------------------------------------------
  ! !  Input Parameters
  ! !  J vector, number of point along each coordinate, initial alpha
  ! !-----------------------------------------------------------------------
  ! subroutine input_parameters()
  !   integer :: m0_code, m1_code, m2_code, total_code
  !   open(1,file='optgrid.config')
  !   read(1,*) n1,n2,n3
  !   read(1,*) nenv1,nenv2,nenv3
  !   read(1,*) a1,a2,a3
  !   read(1, *) envtype1, envtype2, envtype3
  !   read(1, *) min1, max1, Emax1
  !   read(1, *) min2, max2, Emax2
  !
  !   envpath = resolve_relative_exe_path('../extra/spectrum/optgrid/dawes/J00K000')
  !   if (envtype3 == 1) then
  !     read(1, *) min3, max3, Emax3, symmphi
  !   elseif (envtype3 == 2) then
  !     read(1, *) min3, max3, n3, symmphi
  !   end if
  !
  !   read(1,*)cs,eps,rgrid
  !   read(1, *) m1_code, m0_code, m2_code
  !   if(cs.eq.2)then
  !     read(1,*)energyoffset
  !     min3 = pi * min3 / 180
  !     max3 = pi * max3 / 180
  !     energyoffset = energyoffset / autown
  !   endif
  !   close(1)
  !   Emax1 = Emax1 / autown
  !   Emax2 = Emax2 / autown
  !   Emax3 = Emax3 / autown
  !
  !   total_code = m0_code + m1_code + m2_code
  !   if (total_code == 48) then
  !     zpe = zpe_66
  !   else if (total_code == 50) then
  !     zpe = zpe_68
  !   else if (total_code == 52) then
  !     zpe = zpe_88
  !   end if
  ! end subroutine

  !-----------------------------------------------------------------------
  !  Input Envelopes
  !  Loads potential along MEP for each coordinate
  !-----------------------------------------------------------------------
  ! subroutine input_envelopes()
  !   integer i
  !   open(1, file = append_path_token(envpath, 'MEP1.dat'), status='old')
  !   do i=1,nenv1
  !     read(1,*) env1(i), pot1(i)
  !   enddo
  !   close(1)
  !   pot1 = pot1 / autown + zpe
  !   open(1, file = append_path_token(envpath, 'MEP2.dat'), status='old')
  !   do i=1,nenv2
  !     read(1,*) env2(i), pot2(i)
  !   enddo
  !   close(1)
  !   pot2 = pot2/autown + zpe
  !   open(1, file = append_path_token(envpath, 'MEP3.dat'), status='old')
  !   do i=1,nenv3
  !     read(1,*) env3(i), pot3(i)
  !   enddo
  !   close(1)
  !   pot3 = pot3/autown + zpe
  !
  !   ! Preparation for extrapolation of envelopes
  !   env1_d1 = (pot1(2)-pot1(1))/(env1(2)-env1(1))
  !   env1_dn = (pot1(nenv1)-pot1(nenv1-1))/(env1(nenv1)-env1(nenv1-1))
  !   call spline(env1,pot1,nenv1,env1_d1,env1_dn,env1_2d)
  !   env2_d1 = (pot2(2)-pot2(1))/(env2(2)-env2(1))
  !   env2_dn = (pot2(nenv2)-pot2(nenv2-1))/(env2(nenv2)-env2(nenv2-1))
  !   call spline(env2,pot2,nenv2,env2_d1,env2_dn,env2_2d)
  !   env3_d1 = (pot3(2)-pot3(1))/(env3(2)-env3(1))
  !   env3_dn = (pot3(nenv3)-pot3(nenv3-1))/(env3(nenv3)-env3(nenv3-1))
  !   call spline(env3,pot3,nenv3,env3_d1,env3_dn,env3_2d)
  ! end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Find parabola
!-------------------------------------------------------------------------------------------------------------------------------------------
  ! subroutine find_parabola(env,pot,nenv,enva,envE,envm)
  !   integer nenv,i,im,il,ir
  !   real(real64) env(nenv),pot(nenv),enva,envE,envm
  !   real(real64) mat(3,3),d(3)
  !   real(real64) det0,det1,c
  !   im = 1
  !   do i=1,nenv
  !     if(pot(i).lt.pot(im))im = i
  !   enddo
  !   do il=im,1,-1
  !     if(pot(il)-pot(im).gt.energyoffset)exit
  !   enddo
  !   do ir=im,nenv
  !     if(pot(ir)-pot(im).gt.energyoffset)exit
  !   enddo
  !   mat(1,1) = env(il)**2
  !   mat(1,2) = env(il)
  !   mat(1,3) = 1
  !   mat(2,1) = env(im)**2
  !   mat(2,2) = env(im)
  !   mat(2,3) = 1
  !   mat(3,1) = env(ir)**2
  !   mat(3,2) = env(ir)
  !   mat(3,3) = 1
  !   d(1) = pot(il)
  !   d(2) = pot(im)
  !   d(3) = pot(ir)
  !   det0 = det(mat)
  !   mat(1,1) = d(1)
  !   mat(2,1) = d(2)
  !   mat(3,1) = d(3)
  !   det1 = det(mat)
  !   c = det1 / det0
  !   envE = abs(pot(im))
  !   enva = sqrt(envE/c)
  !   envm = env(im)
  ! end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------------------------------------------------------------------
  real(real64) function det(a)
    real(real64) a(3,3)
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end function
  
  !-----------------------------------------------------------------------
  !  PotVib
  !  Returns interpolated value of envelope at a given point
  !  Adds approximate ZPE value along two other coordinates
  !-----------------------------------------------------------------------
  real(real64) function potv1(r)
    real(real64) r
    if(r>env1(nenv1))then
      potv1 = pot1(nenv1)
    else
      call splint(env1,pot1,env1_d1,env1_2d,nenv1,r,potv1)
    endif
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------------------------------------------------------------------
  real(real64) function potv2(r)
    real(real64) r
    if(r>env2(nenv2))then
      potv2 = pot2(nenv2)
    else
      call splint(env2,pot2,env2_d1,env2_2d,nenv2,r,potv2)
    endif
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------------------------------------------------------------------
  real(real64) function potv3(r)
    real(real64) r
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
  real(real64) function potenv1(r)
    real(real64) r,x
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

!-------------------------------------------------------------------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------------------------------------------------------------------
  real(real64) function potenv2(r)
    real(real64) r,x
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

!-------------------------------------------------------------------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------------------------------------------------------------------
  real(real64) function potenv3(r)
    real(real64) r,x
    if(envtype3.eq.1)then
      x = r - envm3
      potenv3 = - envE3 * 4 / (exp(x/enva3) + exp(-x/enva3))**2
    else
      potenv3 = - envE3
    endif
  end function

  ! !-----------------------------------------------------------------------
  ! !  Grid derivative for rho
  ! !-----------------------------------------------------------------------
  ! subroutine der_rho(x,y,dydx)
  !   real(real64) x,y,dydx
  !   real(real64) argument
  !   argument = 2*mu*(Emax1-potenv1(y)-gridextra(envm1,envm2))
  !   if (argument.gt.0.d0) then
  !     dydx = pi/DSQRT(argument)
  !   else
  !     stop 'Potential is greater than Emax1'
  !   endif
  ! end subroutine

  ! !-----------------------------------------------------------------------
  ! !  Grid derivative for tet
  ! !-----------------------------------------------------------------------
  ! subroutine der_tet(x,y,dydx)
  !   real(real64) x,y,dydx
  !   real(real64) argument
  !   argument = mu/2*(Emax2-potenv2(y)-gridextra(envm1,envm2))
  !   if (argument.gt.0.d0) then
  !     dydx = pi/(envm1*DSQRT(argument))
  !   else
  !     stop 'Potential is greater than Emax2'
  !   endif
  ! end subroutine

  ! !-----------------------------------------------------------------------
  ! !  Grid derivative for phi
  ! !-----------------------------------------------------------------------
  ! subroutine der_phi(x,y,dydx)
  !   real(real64) x,y,dydx
  !   real(real64) argument
  !   argument = mu/2*(Emax3-potenv3(y)-gridextra(envm1,envm2))
  !   if (argument.gt.0.d0) then
  !     dydx = pi/(envm1*sin(envm2)*DSQRT(argument))
  !   else
  !     stop 'Potential is greater than Emax3'
  !   endif
  ! end subroutine

  !-----------------------------------------------------------------------
  !  Extra potential caused by conversion from Jacoby coordinates
  !-----------------------------------------------------------------------
  real(real64) function gridextra(rho,tet)
    real(real64) rho,tet
    gridextra = - 2/(mu*rho**2)*(sin(2*tet)**(-2)+0.0625d0)
  end function

  ! !-----------------------------------------------------------------------
  ! !  Grid derivative for R1
  ! !-----------------------------------------------------------------------
  ! subroutine der_r1(x,y,dydx)
  !   real(real64) x,y,dydx
  !   real(real64) argument
  !   argument = 2*mu1*(Emax1-potenv1(y))
  !   if (argument.gt.0.d0) then
  !     dydx = pi/DSQRT(argument)
  !   else
  !     stop 'Potential is greater than Emax1'
  !   endif
  ! end subroutine

  ! !-----------------------------------------------------------------------
  ! !  Grid derivative for R2
  ! !-----------------------------------------------------------------------
  ! subroutine der_r2(x,y,dydx)
  !   real(real64) x,y,dydx
  !   real(real64) argument
  !   argument = 2*mu2*(Emax2-potenv2(y))
  !   if (argument.gt.0.d0) then
  !     dydx = pi/DSQRT(argument)
  !   else
  !     stop 'Potential is greater than Emax2'
  !   endif
  ! end subroutine

  ! !-----------------------------------------------------------------------
  ! !  Grid derivative for Te
  ! !-----------------------------------------------------------------------
  ! subroutine der_te(x,y,dydx)
  !   real(real64) x,y,dydx
  !   real(real64) argument
  !   argument = 2/(1/(mu1*envm1**2) + 1/(mu2*envm2**2) - 2*cos(envm3)/(m0*envm1*envm2)) * (Emax3-potenv3(y))
  !   if (argument.gt.0.d0) then
  !     dydx = pi/DSQRT(argument)
  !   else
  !     stop 'Potential is greater than Emax3'
  !   endif
  ! end subroutine
end

