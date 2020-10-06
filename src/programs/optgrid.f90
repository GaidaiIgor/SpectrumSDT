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
  real(real64) :: env_emax
  integer :: env_npoints
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
  env_emax = 414.4466 / autown
  env_npoints = 494

  if (params % optimized_grid_rho == 1) then
    allocate(env_grid(env_npoints), env_values(env_npoints), spline_deriv_2nd(env_npoints))
    call input_envelopes(params)
    call find_parabola(env_grid, env_values, env_npoints, env_fit_param, env_fit_min_y, env_fit_min_x)
    if (params % grid_rho_npoints == -1) then
      call generate_optimized_grid_step(params % grid_rho_from, params % grid_rho_to, env_fit_min_x, params % grid_rho_step, optgrid_diff_rhs, grid_rho, jac_rho, params % grid_rho_npoints)
    else
      call generate_optimized_grid_points(params % grid_rho_from, params % grid_rho_to, env_fit_min_x, params % grid_rho_npoints, optgrid_diff_rhs, grid_rho, jac_rho, params % grid_rho_step)
    end if
  else
    if (params % grid_rho_npoints == -1) then
      call generate_equidistant_grid_step(params % grid_rho_from, params % grid_rho_to, params % grid_rho_step, grid_rho, jac_rho, params % grid_rho_npoints)
    else
      call generate_equidistant_grid_points(params % grid_rho_from, params % grid_rho_to, params % grid_rho_npoints, grid_rho, jac_rho, params % grid_rho_step)
    end if
  end if

  if (params % grid_theta_npoints == -1) then
    call generate_equidistant_grid_step(params % grid_theta_from, params % grid_theta_to, params % grid_theta_step, grid_theta, jac_theta, params % grid_theta_npoints)
  else
    call generate_equidistant_grid_points(params % grid_theta_from, params % grid_theta_to, params % grid_theta_npoints, grid_theta, jac_theta, params % grid_theta_step)
  end if

  if (params % grid_phi_npoints == -1) then
    call generate_equidistant_grid_step(params % grid_phi_from, params % grid_phi_to, params % grid_phi_step, grid_phi, jac_phi, params % grid_phi_npoints)
  else
    call generate_equidistant_grid_points(params % grid_phi_from, params % grid_phi_to, params % grid_phi_npoints, grid_phi, jac_phi, params % grid_phi_step)
  end if

  call print_grid(grid_rho, jac_rho, params % grid_rho_step, 'grid_rho.dat')
  call print_grid(grid_theta, jac_theta, params % grid_theta_step, 'grid_theta.dat')
  call print_grid(grid_phi, jac_phi, params % grid_phi_step, 'grid_phi.dat')
  print *, 'Done'

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads rho MEP and prepares spline information.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine input_envelopes(params)
    class(input_params), intent(in) :: params
    integer :: i
    real(real64) :: spline_deriv_left, spline_deriv_right

    open(1, file = params % envelope_rho_path, status='old')
    do i = 1, env_npoints
      read(1, *) env_grid(i), env_values(i)
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
  impure elemental function envelope_potential(r) result(res)
    real(real64), intent(in) :: r
    real(real64) :: res
    real(real64) :: x

    x = r - env_fit_min_x
    if (x < 0) then
      res = -env_fit_min_y * 4 / (exp(x/env_fit_param) + exp(-x/env_fit_param))**2
    else
      if (r > env_grid(env_npoints)) then
        res = env_values(env_npoints)
      else
        res = splint(env_grid, env_values, spline_deriv_2nd, r)
      endif
    endif
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Right-hand side of differential equation for optimized grid generation.
! Represents a 1-equation system, so only the first element of y and dydx is used.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine optgrid_diff_rhs(x, y, dydx)
    real(real64), intent(in) :: x ! not used, but required for interface conformance
    real(real64), intent(in) :: y(:)
    real(real64), intent(out) :: dydx(:)
    real(real64), allocatable :: argument(:)

    argument = 2 * mu * (env_emax - envelope_potential(y))
    call assert(all(argument > 0), 'Error: potential is greater than env_emax')
    dydx = pi / sqrt(argument)
  end subroutine

end
