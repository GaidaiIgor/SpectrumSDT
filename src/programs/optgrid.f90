!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures related to `grids` stage
!-------------------------------------------------------------------------------------------------------------------------------------------
module optgrid_mod
  use config_mod
  use constants
  use general_utils
  use input_params_mod
  use iso_fortran_env, only: real64
  use numerical_recipies
  use path_utils
  use vector_mod
  implicit none

  private
  public :: generate_grids

  real(real64), parameter :: m0 = isomass(1), &
                        m1 = isomass(1), &
                        m2 = isomass(1), &
                        Mtot = m0+m1+m2, &
                        mu = sqrt(m0*m1*m2/Mtot)
                        
  real(real64) :: env_emax
  integer :: env_npoints
  real(real64), allocatable :: env_grid(:), env_values(:), spline_deriv_2nd(:)
  real(real64) :: env_fit_param, env_fit_min_y, env_fit_min_x

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
! Calculates 3x3 matrix determinant.
!-------------------------------------------------------------------------------------------------------------------------------------------
  real(real64) function det(a)
    real(real64) a(3,3)
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end function
  
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
! Simplified interface to rkdumb for the case of a signle equation to use scalars instead of arrays of size 1.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function rkdumb_single(vstart, x1, x2, nstep, derivs) result(v)
    real(real64) :: vstart, x1, x2
    integer :: nstep
    procedure(diff_equations_rhs) :: derivs
    real(real64) :: v
    real(real64) :: vstart_arr(1), v_arr(1)

    vstart_arr(1) = vstart
    v_arr = rkdumb(vstart_arr, x1, x2, nstep, derivs)
    v = v_arr(1)
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
! Returns the values of derivative *dy/dx* for all values of *y* and a given value of *x*.
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

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates an optimized grid from *from* to *to*.
! Step specifies step on grid pre-image. Step on the actual grid is variable. Number of points is adjusted.
! *start* is a point between *from* and *to* from where generation begins in left and right direction.
! This allows for symmetric distribution of points if *deriv* is symmetric.
! *deriv* evaluates the right-hand side of the differential equation that defines the optimized grid.
! Resulting optimized grid is returned in *grid*. The corresponding values of *deriv* at grid points are returned in *jac*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_optimized_grid_step(from, to, start, step, deriv, grid, jac, npoints)
    real(real64), intent(in) :: from, to, start, step
    procedure(diff_equations_rhs) :: deriv
    real(real64), allocatable, intent(out) :: grid(:), jac(:)
    integer, optional, intent(out) :: npoints
    real(real64) :: prev_x, next_x, next_point
    type(vector_real) :: grid_left, grid_right

    ! Generate right branch (relative to minimum at start)
    next_point = start
    prev_x = 0
    next_x = step / 2
    grid_right = vector_real()
    do
      next_point = rkdumb_single(next_point, prev_x, next_x, 2048, deriv)
      if (next_point > to) then
        exit
      end if
      call grid_right % push(next_point)
      prev_x = next_x
      next_x = next_x + step
    end do

    ! Generate left branch (relative to minimum at start)
    next_point = start
    prev_x = 0
    next_x = -step / 2
    grid_left = vector_real()
    do
      next_point = rkdumb_single(next_point, prev_x, next_x, 2048, deriv)
      ! Add an extra point on left if the total number of points is odd
      ! if (next_point < from .and. mod(grid_right % get_size() + grid_left % get_size(), 2) == 0) then
      if (next_point < from) then
        exit
      end if
      call grid_left % push(next_point)
      prev_x = next_x
      next_x = next_x - step
    end do

    ! Merge branches: grid_left is used as common storage for both branches
    call grid_left % reverse()
    call grid_left % append_vector(grid_right)
    grid = grid_left % to_array()
    allocate(jac(size(grid)))
    call deriv(0d0, grid, jac)

    if (present(npoints)) then
      npoints = size(grid)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calls `generate_optimized_grid_step` for different values of *step* until desired *number of points* is achieved.
! Final value of *step* is returned.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_optimized_grid_points(from, to, start, npoints, deriv, grid, jac, step)
    real(real64), intent(in) :: from, to, start
    integer, intent(in) :: npoints
    procedure(diff_equations_rhs) :: deriv
    real(real64), allocatable, intent(out) :: grid(:), jac(:)
    real(real64), optional, intent(out) :: step
    real(real64) :: step_power_left, step_power_right, step_power_middle ! step search boundaries

    ! Find initial left boundary
    step_power_left = -4
    do
      call generate_optimized_grid_step(from, to, start, 2 ** step_power_left, deriv, grid, jac)
      if (size(grid) > npoints) then
        exit
      end if
      step_power_left = step_power_left - 1
    end do

    ! Find initial right boundary
    step_power_right = 0
    do
      call generate_optimized_grid_step(from, to, start, 2 ** step_power_right, deriv, grid, jac)
      if (size(grid) < npoints) then
        exit
      end if
      step_power_right = step_power_right + 1
    end do

    ! Find step matching specified number of points
    do
      step_power_middle = (step_power_left + step_power_right) / 2
      call generate_optimized_grid_step(from, to, start, 2 ** step_power_middle, deriv, grid, jac)
      if (size(grid) == npoints) then
        exit
      end if
      if (size(grid) > npoints) then
        step_power_left = step_power_middle
      end if
      if (size(grid) < npoints) then
        step_power_right = step_power_middle
      end if
    end do

    if (present(step)) then
      step = 2 ** step_power_middle
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates grid with required step size and boundaries. Number of points is adjusted.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_equidistant_grid_step(from, to, step, grid, jac, npoints)
    real(real64), intent(in) :: from, to, step
    real(real64), allocatable, intent(out) :: grid(:), jac(:)
    integer, optional, intent(out) :: npoints
    
    grid = generate_real_range(from + step/2, to - step/2, step)
    allocate(jac(size(grid)))
    jac = 1
    if (present(npoints)) then
      npoints = size(grid)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates grid with required number of points and boundaries. Step size is adjusted.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_equidistant_grid_points(from, to, npoints, grid, jac, step)
    real(real64), intent(in) :: from, to
    integer, intent(in) :: npoints
    real(real64), allocatable, intent(out) :: grid(:), jac(:)
    real(real64), optional, intent(out) :: step
    real(real64) :: step_act

    step_act = (to - from) / npoints
    grid = linspace(from + step_act/2, to - step_act/2, npoints)
    allocate(jac(size(grid)))
    jac = 1
    if (present(step)) then
      step = step_act
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints grid, Jacobian and step to a file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_grid(grid, jac, step, file_name)
    real(real64), intent(in) :: grid(:), jac(:)
    character(*), intent(in) :: file_name
    real(real64), intent(in) :: step
    integer :: i, file_unit
    
    call assert(size(grid) == size(jac), 'Error: grid and jac have to have the same size')
    open(newunit = file_unit, file = file_name)
    write(file_unit, *) size(grid), step
    do i = 1, size(grid)
      write(file_unit, '(2G23.15)') grid(i), jac(i)
    enddo
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Top level procedure for grids generation
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_grids(params)
    class(input_params), intent(in) :: params
    real(real64) :: rho_step, theta_step, phi_step
    real(real64), allocatable :: grid_rho(:), grid_theta(:), grid_phi(:)
    real(real64), allocatable :: jac_rho(:), jac_theta(:), jac_phi(:)

    env_emax = 414.4466 / autown
    env_npoints = 494
    rho_step = params % grid_rho_step
    theta_step = params % grid_theta_step
    phi_step = params % grid_phi_step

    if (params % use_optimized_grid_rho == 1) then
      allocate(env_grid(env_npoints), env_values(env_npoints), spline_deriv_2nd(env_npoints))
      call input_envelopes(params)
      call find_parabola(env_grid, env_values, env_npoints, env_fit_param, env_fit_min_y, env_fit_min_x)
      if (params % grid_rho_npoints == -1) then
        call generate_optimized_grid_step(params % grid_rho_from, params % grid_rho_to, env_fit_min_x, params % grid_rho_step, optgrid_diff_rhs, grid_rho, jac_rho)
      else
        call generate_optimized_grid_points(params % grid_rho_from, params % grid_rho_to, env_fit_min_x, params % grid_rho_npoints, optgrid_diff_rhs, grid_rho, jac_rho, rho_step)
      end if
    else
      if (params % grid_rho_npoints == -1) then
        call generate_equidistant_grid_step(params % grid_rho_from, params % grid_rho_to, params % grid_rho_step, grid_rho, jac_rho)
      else
        call generate_equidistant_grid_points(params % grid_rho_from, params % grid_rho_to, params % grid_rho_npoints, grid_rho, jac_rho, rho_step)
      end if
    end if

    if (params % grid_theta_npoints == -1) then
      call generate_equidistant_grid_step(params % grid_theta_from, params % grid_theta_to, params % grid_theta_step, grid_theta, jac_theta)
    else
      call generate_equidistant_grid_points(params % grid_theta_from, params % grid_theta_to, params % grid_theta_npoints, grid_theta, jac_theta, theta_step)
    end if

    if (params % grid_phi_npoints == -1) then
      call generate_equidistant_grid_step(params % grid_phi_from, params % grid_phi_to, params % grid_phi_step, grid_phi, jac_phi)
    else
      call generate_equidistant_grid_points(params % grid_phi_from, params % grid_phi_to, params % grid_phi_npoints, grid_phi, jac_phi, phi_step)
    end if

    call print_grid(grid_rho, jac_rho, rho_step, 'grid_rho.dat')
    call print_grid(grid_theta, jac_theta, theta_step, 'grid_theta.dat')
    call print_grid(grid_phi, jac_phi, phi_step, 'grid_phi.dat')
  end subroutine
end
