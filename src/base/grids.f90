!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures and structures related to grids generation and storage.
!-------------------------------------------------------------------------------------------------------------------------------------------
module grids_mod
  use config_mod
  use constants
  use coordinate_coversion_mod
  use formulas_mod, only: get_reduced_mass
  use general_utils
  use grid_info_mod
  use input_params_mod
  use io_utils
  use iso_fortran_env, only: real64
  use numerical_recipies
  use path_utils
  use spectrumsdt_paths_mod
  use vector_mod
  implicit none

  private
  public :: generate_grids, load_grids

  real(real64) :: mu, env_emax
  real(real64), allocatable :: env_grid(:), env_values(:), spline_deriv_2nd(:)
  real(real64) :: env_fit_param, env_fit_min_abs_energy, env_fit_min_x

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads rho MEP and prepares spline information.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine input_envelopes(params)
    class(input_params), intent(in) :: params
    real(real64) :: spline_deriv_left, spline_deriv_right
    real(real64), allocatable :: envelope_matrix(:, :)

    envelope_matrix = read_matrix_real(params % grid_rho % envelope_path)
    env_grid = envelope_matrix(:, 1)
    env_values = envelope_matrix(:, 2)

    ! Calculates derivatives at the endpoints for spline
    spline_deriv_left = (env_values(2) - env_values(1)) / (env_grid(2) - env_grid(1))
    spline_deriv_right = (env_values(size(env_values)) - env_values(size(env_values) - 1)) / (env_grid(size(env_grid)) - env_grid(size(env_grid) - 1))
    allocate(spline_deriv_2nd(size(env_values)))
    call spline(env_grid, env_values, spline_deriv_left, spline_deriv_right, spline_deriv_2nd)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates 3x3 matrix determinant.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function det(a) result(res)
    real(real64) :: a(3, 3)
    real(real64) :: res
    res = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Finds parameters of fitting parabola.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine find_parabola(env, pot, enva, envE, envm)
    real(real64), intent(in) :: env(:), pot(:)
    real(real64), intent(out) :: enva, envE, envm
    integer :: im
    real(real64) :: d(3)
    real(real64) :: mat(3, 3)
    real(real64) :: det0, det1, c

    im = minloc(pot, 1)
    mat(:, 1) = env(im-1 : im+1) ** 2
    mat(:, 2) = env(im-1 : im+1)
    mat(:, 3) = 1
    d = pot(im-1 : im+1)
    det0 = det(mat)
    mat(:, 1) = d
    det1 = det(mat)
    c = det1 / det0

    envE = abs(pot(im))
    enva = sqrt(envE/c)
    envm = env(im)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates envelope potential for optimized grid. Left side is represented by Eckart function fitted on rho MEP, right by rho MEP itself.
!-------------------------------------------------------------------------------------------------------------------------------------------
  impure elemental function envelope_potential(r) result(res)
    real(real64), intent(in) :: r
    real(real64) :: res
    real(real64) :: x

    x = r - env_fit_min_x
    if (x < 0) then
      res = -env_fit_min_abs_energy * 4 / (exp(x/env_fit_param) + exp(-x/env_fit_param))**2
    else
      if (r > env_grid(size(env_grid))) then
        res = env_values(size(env_values))
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
! Simplified interface to rkdumb for the case of a signle equation to use scalars instead of arrays of size 1.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function rkdumb_single(vstart, x1, x2, rk_nsteps, derivs) result(v)
    real(real64) :: vstart, x1, x2
    integer :: rk_nsteps
    procedure(diff_equations_rhs) :: derivs
    real(real64) :: v
    real(real64) :: vstart_arr(1), v_arr(1)

    vstart_arr(1) = vstart
    v_arr = rkdumb(vstart_arr, x1, x2, rk_nsteps, derivs)
    v = v_arr(1)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates an optimized grid from *from* to *to*.
! Step specifies step on grid pre-image. Step on the actual grid is variable. Number of points is adjusted.
! *start* is a point between *from* and *to* from where generation begins in left and right direction.
! This allows for symmetric distribution of points if *deriv* is symmetric.
! *deriv* evaluates the right-hand side of the differential equation that defines the optimized grid.
! Resulting optimized grid is returned in *grid*. The corresponding values of *deriv* at grid points are returned in *jac*.
! *rk_nsteps* controls how many steps are used for differential equation propagation.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_optimized_grid_step(from, to, start, step, rk_nsteps, deriv, grid, jac, npoints)
    real(real64), intent(in) :: from, to, start, step
    integer, intent(in) :: rk_nsteps
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
      next_point = rkdumb_single(next_point, prev_x, next_x, rk_nsteps, deriv)
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
      next_point = rkdumb_single(next_point, prev_x, next_x, rk_nsteps, deriv)
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
  subroutine generate_optimized_grid_points(from, to, start, npoints, rk_nsteps, deriv, grid, jac, step)
    real(real64), intent(in) :: from, to, start
    integer, intent(in) :: npoints, rk_nsteps
    procedure(diff_equations_rhs) :: deriv
    real(real64), allocatable, intent(out) :: grid(:), jac(:)
    real(real64), optional, intent(out) :: step
    real(real64) :: step_power_left, step_power_right, step_power_middle ! step search boundaries

    ! Find initial left boundary
    step_power_left = -4
    do
      call generate_optimized_grid_step(from, to, start, 2 ** step_power_left, rk_nsteps, deriv, grid, jac)
      if (size(grid) > npoints) then
        exit
      end if
      step_power_left = step_power_left - 1
    end do

    ! Find initial right boundary
    step_power_right = 0
    do
      call generate_optimized_grid_step(from, to, start, 2 ** step_power_right, rk_nsteps, deriv, grid, jac)
      if (size(grid) < npoints) then
        exit
      end if
      step_power_right = step_power_right + 1
    end do

    ! Find step matching specified number of points
    do
      step_power_middle = (step_power_left + step_power_right) / 2
      call generate_optimized_grid_step(from, to, start, 2 ** step_power_middle, rk_nsteps, deriv, grid, jac)
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
! Prints grid, its ends, step and Jacobian to a file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_grid(grid, grid_ends, step, file_name, jac)
    real(real64), intent(in) :: grid(:)
    real(real64), intent(in) :: grid_ends(2)
    real(real64), intent(in) :: step
    character(*), intent(in) :: file_name
    real(real64), optional, intent(in) :: jac(:)
    integer :: i, file_unit
    
    open(newunit = file_unit, file = file_name)
    write(file_unit, '(3G23.15,1x,I0)') grid_ends, step, size(grid)
    if (present(jac)) then
      call assert(size(grid) == size(jac), 'Error: grid and jac have to have the same size')
      do i = 1, size(grid)
        write(file_unit, '(2G23.15)') grid(i), jac(i)
      end do
    else
      write(file_unit, '(G23.15)') grid
    end if
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Forms a list where each row is a combination of elements from the given arrays.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function product_list(array1, array2, array3) result(list)
    real(real64), intent(in) :: array1(:), array2(:), array3(:)
    real(real64) :: list(3, size(array1) * size(array2) * size(array3))
    integer :: i1, i2, i3, j

    j = 1
    do i1 = 1, size(array1)
      do i2 = 1, size(array2)
        do i3 = 1, size(array3)
          list(1, j) = array1(i1)
          list(2, j) = array2(i2)
          list(3, j) = array3(i3)
          j = j + 1
        end do
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints a file specified by *file_name* with coordinates at which the values of PES are requested.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_pes_request(grid_rho, grid_theta, grid_phi, file_name, coord_system, treat_tp_as_xy, mass)
    real(real64), intent(in) :: grid_rho(:), grid_theta(:), grid_phi(:)
    character(*), intent(in) :: file_name, coord_system
    integer, intent(in) :: treat_tp_as_xy
    real(real64), optional, intent(in) :: mass(:)
    integer, parameter :: col_width = 25
    integer :: file_unit
    real(real64), allocatable :: coord_list(:, :)
    character(:), allocatable :: col1_header, col2_header, col3_header

    call assert(any(coord_system == [character(100) :: 'jacobi', 'cartesian', 'all bonds', 'internal']) .means. present(mass), 'Error: mass has to be given for this coordinate system')
    print *, 'Converting coordinates...'
    if (treat_tp_as_xy == 1) then
      coord_list = product_list(grid_theta, grid_theta, grid_rho)
      coord_list = convert_aphxyz_to_aph(coord_list)
    else
      coord_list = product_list(grid_rho, grid_theta, grid_phi)
    end if

    select case (coord_system)
      case ('aph')
        col1_header = 'rho (Bohr)'
        col2_header = 'theta (rad)'
        col3_header = 'phi (rad)'
      case ('mass jacobi')
        col1_header = 's (Bohr)'
        col2_header = 'S (Bohr)'
        col3_header = 'Theta (rad)'
        coord_list = convert_aph_to_mass_jacobi(coord_list)
      case ('jacobi')
        col1_header = 'r (Bohr)'
        col2_header = 'R (Bohr)'
        col3_header = 'Theta (rad)'
        coord_list = convert_aph_to_jacobi(coord_list, mass)
      case ('cartesian')
        col1_header = 'x2 (Bohr)'
        col2_header = 'y2 (Bohr)'
        col3_header = 'x3 (Bohr)'
        coord_list = convert_aph_to_cartesian(coord_list, mass)
      case ('all bonds')
        col1_header = 'bond12 (Bohr)'
        col2_header = 'bond13 (Bohr)'
        col3_header = 'bond23 (Bohr)'
        coord_list = convert_aph_to_all_bonds(coord_list, mass)
      case ('internal')
        col1_header = 'bond12 (Bohr)'
        col2_header = 'bond23 (Bohr)'
        col3_header = 'angle123 (rad)'
        coord_list = convert_aph_to_internal(coord_list, mass)
      case default
        call assert(.false., 'Error: unknown coordinate system')
    end select

    print *, 'Writing coordinates file...'
    open(newunit = file_unit, file = file_name)
    write(file_unit, *) size(coord_list, 2)
    write(file_unit, '(3A' // num2str(col_width) // ')') align_center(col1_header, col_width), align_center(col2_header, col_width), align_center(col3_header, col_width)
    write(file_unit, '(3G' // num2str(col_width) // '.15)') coord_list
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Top level procedure for grids generation.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_grids(params)
    class(input_params), intent(in) :: params
    real(real64) :: rho_step, theta_step, phi_step
    real(real64), allocatable :: grid_rho(:), grid_theta(:), grid_phi(:)
    real(real64), allocatable :: jac_rho(:), jac_theta(:), jac_phi(:)

    mu = get_reduced_mass(params % mass)
    env_emax = params % grid_rho % max_energy / au_to_wn
    rho_step = params % grid_rho % step
    theta_step = params % grid_theta % step

    if (params % grid_rho % optimized == 1) then
      call input_envelopes(params)
      call find_parabola(env_grid, env_values, env_fit_param, env_fit_min_abs_energy, env_fit_min_x)
      if (params % grid_rho % num_points == -1) then
        call generate_optimized_grid_step(params % grid_rho % from, params % grid_rho % to, env_fit_min_x, params % grid_rho % step, &
            params % grid_rho % solver_steps, optgrid_diff_rhs, grid_rho, jac_rho)
      else
        call generate_optimized_grid_points(params % grid_rho % from, params % grid_rho % to, env_fit_min_x, params % grid_rho % num_points, &
            params % grid_rho % solver_steps, optgrid_diff_rhs, grid_rho, jac_rho, rho_step)
      end if
    else
      if (params % grid_rho % num_points == -1) then
        call generate_equidistant_grid_step(params % grid_rho % from, params % grid_rho % to, params % grid_rho % step, grid_rho, jac_rho)
      else
        call generate_equidistant_grid_points(params % grid_rho % from, params % grid_rho % to, params % grid_rho % num_points, grid_rho, jac_rho, rho_step)
      end if
    end if

    if (params % grid_theta % num_points == -1) then
      call generate_equidistant_grid_step(params % grid_theta % from, params % grid_theta % to, params % grid_theta % step, grid_theta, jac_theta)
    else
      call generate_equidistant_grid_points(params % grid_theta % from, params % grid_theta % to, params % grid_theta % num_points, grid_theta, jac_theta, theta_step)
    end if

    call generate_equidistant_grid_points(0d0, 2*pi, params % num_points_phi, grid_phi, jac_phi, phi_step)

    call write_grid(grid_rho, [params % grid_rho % from, params % grid_rho % to], rho_step, 'grid_rho.dat', jac = jac_rho)
    call write_grid(grid_theta, [params % grid_theta % from, params % grid_theta % to], theta_step, 'grid_theta.dat')
    call write_grid(grid_phi, [0d0, 2*pi], phi_step, 'grid_phi.dat')

    if (params % output_coordinate_system /= 'aph' .or. params % treat_tp_as_xy == 1) then
      call write_pes_request(grid_rho, grid_theta, grid_phi, 'pes.in', params % output_coordinate_system, params % treat_tp_as_xy, params % mass)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads grid information from disk.
! Jacobian in theta_info and phi_info is not allocated since they are equidistant.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_grids(params, rho_info, theta_info, phi_info)
    class(input_params), intent(in) :: params
    type(grid_info), intent(out) :: rho_info, theta_info, phi_info
    integer :: file_unit, num_points, i

    open(newunit = file_unit, file = get_rho_info_path(params))
    read(file_unit, *) rho_info % from, rho_info % to, rho_info % step, num_points
    allocate(rho_info % points(num_points), rho_info % jac(num_points))
    do i = 1, num_points
      read(file_unit, *) rho_info % points(i), rho_info % jac(i)
    end do
    close(file_unit)

    open(newunit = file_unit, file = get_theta_info_path(params))
    read(file_unit, *) theta_info % from, theta_info % to, theta_info % step, num_points
    allocate(theta_info % points(num_points))
    read(file_unit, *) theta_info % points
    close(file_unit)

    open(newunit = file_unit, file = get_phi_info_path(params))
    read(file_unit, *) phi_info % from, phi_info % to, phi_info % step, num_points
    allocate(phi_info % points(num_points))
    read(file_unit, *) phi_info % points
    close(file_unit)
  end subroutine

end module
