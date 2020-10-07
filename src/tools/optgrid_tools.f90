module optgrid_tools
  use constants
  use general_utils
  use iso_fortran_env, only: real64
  use numerical_recipies
  use vector_mod
  implicit none

contains

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
      print *, 'left npts:', size(grid), step_power_left
      if (size(grid) > npoints) then
        exit
      end if
      step_power_left = step_power_left - 1
    end do

    ! Find initial right boundary
    step_power_right = 0
    do
      call generate_optimized_grid_step(from, to, start, 2 ** step_power_right, deriv, grid, jac)
      print *, 'right npts:', size(grid), step_power_right
      if (size(grid) < npoints) then
        exit
      end if
      step_power_right = step_power_right + 1
    end do

    ! Find step matching specified number of points
    do
      step_power_middle = (step_power_left + step_power_right) / 2
      call generate_optimized_grid_step(from, to, start, 2 ** step_power_middle, deriv, grid, jac)
      print *, 'middle npts:', size(grid), step_power_middle
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

end module
