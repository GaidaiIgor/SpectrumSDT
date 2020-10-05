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
                        mu = sqrt(m0*m1*m2/Mtot)
                        
  real(real64),allocatable :: grid_rho(:), grid_theta(:), grid_phi(:)
  real(real64),allocatable :: jac_rho(:), jac_theta(:), jac_phi(:)
  real(real64) :: rho_step, theta_step, phi_step
  real(real64) :: rho_from, rho_to
  real(real64) :: theta_from, theta_to
  real(real64) :: phi_from, phi_to
  real(real64) :: env_emax, eps
  integer :: env_npoints
  integer :: rho_npoints, theta_npoints, phi_npoints
  !--- Envelopes ---
  character(:), allocatable :: env_path
  real(real64), allocatable :: env_grid(:), env_values(:), spline_deriv_2nd(:)
  real(real64) :: env_fit_param, env_fit_min_y, env_fit_min_x
end module

program optgrid
  use config_mod
  use constants
  use general_utils
  use global_vars
  use input_params_mod
  use iso_fortran_env, only: real64
  use numerical_recipies
  use optgrid_tools
  use path_utils

  implicit none
  type(input_params) :: params

  params = process_user_settings('spectrumsdt.config')
  call input_parameters(params)

  if (params % optimized_grid_rho == 1) then
    allocate(env_grid(env_npoints), env_values(env_npoints), spline_deriv_2nd(env_npoints))
    call input_envelopes()
    call find_parabola(env_grid, env_values, env_npoints, env_fit_param, env_fit_min_y, env_fit_min_x)
    call generate_gridn(grid_rho, jac_rho, params % grid_rho_npoints, params % grid_rho_step, params % grid_rho_from, params % grid_rho_to, env_fit_min_x, eps, der_rho)
  else
    if (params % grid_rho_npoints == -1) then
      call generate_equidistant_grid_step(params % grid_rho_from, params % grid_rho_to, params % grid_rho_step, grid_rho, jac_rho, rho_npoints)
    else
      call generate_equidistant_grid_points(params % grid_rho_from, params % grid_rho_to, params % grid_rho_npoints, grid_rho, jac_rho, rho_step)
    end if
  end if

  if (params % grid_theta_npoints == -1) then
    call generate_equidistant_grid_step(params % grid_theta_from, params % grid_theta_to, params % grid_theta_step, grid_theta, jac_theta, theta_npoints)
  else
    call generate_equidistant_grid_points(params % grid_theta_from, params % grid_theta_to, params % grid_theta_npoints, grid_theta, jac_theta, theta_step)
  end if

  if (params % grid_phi_npoints == -1) then
    call generate_equidistant_grid_step(params % grid_phi_from, params % grid_phi_to, params % grid_phi_step, grid_phi, jac_phi, phi_npoints)
  else
    call generate_equidistant_grid_points(params % grid_phi_from, params % grid_phi_to, params % grid_phi_npoints, grid_phi, jac_phi, phi_step)
  end if

  call print_grid(grid_rho, jac_rho, rho_step, 'grid_rho.dat')
  call print_grid(grid_theta, jac_theta, theta_step, 'grid_theta.dat')
  call print_grid(grid_phi, jac_phi, phi_step, 'grid_phi.dat')
  print *, 'Done'

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
  !   call generate_gridn(g1,jac1,n1,a1,min1,max1,env_fit_min_x,eps,der_rho)
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
! Loads input parameters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine input_parameters(params)
    class(input_params), intent(in) :: params

    rho_from = params % grid_rho_from
    rho_to = params % grid_rho_to
    rho_npoints = params % grid_rho_npoints
    rho_step = params % grid_rho_step
    env_path = params % envelope_rho_path

    theta_from = params % grid_theta_from
    theta_to = params % grid_theta_to
    theta_npoints = params % grid_theta_npoints
    theta_step = params % grid_theta_step

    phi_from = params % grid_phi_from
    phi_to = params % grid_phi_to
    phi_npoints = params % grid_phi_npoints
    phi_step = params % grid_phi_step

    env_emax = 414.4466 / autown
    env_npoints = 494
    eps = 1d-12
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads rho MEP and prepares spline information.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine input_envelopes()
    integer :: i
    real(real64) :: spline_deriv_left, spline_deriv_right
    open(1, file = append_path_token(env_path, 'MEP_rho.dat'), status='old')
    do i=1,env_npoints
      read(1,*) env_grid(i), env_values(i)
    enddo
    close(1)

    ! Calculates derivatives at the endpoints for spline
    spline_deriv_left = (env_values(2) - env_values(1)) / (env_grid(2) - env_grid(1))
    spline_deriv_right = (env_values(env_npoints) - env_values(env_npoints-1)) / (env_grid(env_npoints) - env_grid(env_npoints-1))
    call spline(env_grid, env_values, spline_deriv_left, spline_deriv_right, spline_deriv_2nd)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Finds parameters of fitting parabola.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine find_parabola(env,pot,nenv,enva,envE,envm)
    integer nenv,i,im,il,ir
    real(real64) env(nenv),pot(nenv),enva,envE,envm
    real(real64) mat(3,3),d(3)
    real(real64) det0,det1,c
    im = 1
    do i=1,nenv
      if (pot(i) < pot(im)) im = i
    enddo
    il = im - 1
    ir = im + 1
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

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates 3x3 matrix determinant.
!-------------------------------------------------------------------------------------------------------------------------------------------
  real(real64) function det(a)
    real(real64) a(3,3)
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates envelope potential for optimized grid. Left side is represented by Eckart function fitted on rho MEP, right by rho MEP itself.
!-------------------------------------------------------------------------------------------------------------------------------------------
  real(real64) function envelope_potential(r)
    real(real64) r,x
    x = r - env_fit_min_x
    if (x < 0) then
      envelope_potential = -env_fit_min_y * 4 / (exp(x/env_fit_param) + exp(-x/env_fit_param))**2
    else
      if (r > env_grid(env_npoints)) then
        envelope_potential = env_values(env_npoints)
      else
        envelope_potential = splint(env_grid, env_values, spline_deriv_2nd, r)
      endif
    endif
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Right-hand side of differential equation for optimized grid generation.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine der_rho(x,y,dydx)
    real(real64) x,y,dydx
    real(real64) argument
    argument = 2 * mu * (env_emax - envelope_potential(y))
    if (argument > 0) then
      dydx = pi / sqrt(argument)
    else
      stop 'Potential is greater than env_emax'
    endif
  end subroutine

end
