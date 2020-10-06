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

! !-------------------------------------------------------------------------------------------------------------------------------------------
! ! Generates a grid with fixed left end and free right end.
! ! Alpha is fixed.
! !-------------------------------------------------------------------------------------------------------------------------------------------
!   subroutine generate_gridr(grid,jac,n,alpha,minr,eps,der)
!     integer i,n,nok,nbad
!     real(real64) grid(n),jac(n),x(0:n)
!     real(real64) alpha,minr,eps,nextr
!     real(real64) :: nextr_ode(1)
!     external der
!
!     x(0) = 0
!     do i=1,n
!       x(i) = alpha*i - alpha/2
!     enddo
!     nextr = minr
!     do i=1,n
!       nextr_ode(1) = nextr
!       call odeint(nextr_ode,1,x(i-1),x(i),eps,alpha/10.d0,0.d0,nok,nbad,der,rkqs)
!       nextr = nextr_ode(1)
!
!       grid(i) = nextr
!       call der(0.0d0,nextr,jac(i))
!     enddo
!   end subroutine

! !-------------------------------------------------------------------------------------------------------------------------------------------
! ! Generates a grid within required boundaries.
! ! Solves from the left to the right.
! ! Alpha is adjusted.
! !-------------------------------------------------------------------------------------------------------------------------------------------
!   subroutine generate_grida(grid,jac,n,alpha,minr,maxr,eps,der)
!     integer i,n,nok,nbad
!     real(real64) grid(n),jac(n),x(0:n+1)
!     real(real64) alpha,alpha1,alpha2,minr,maxr,nextr,eps
!     real(real64) :: nextr_ode(1)
!     logical veryclose,allinrange
!     external der
!
!     ! Iteratively change alpha and regenerate the grid
!     veryclose = .false.
!     alpha1 = 0
!     do
!       ! Calculate working grid
!       x(0) = 0
!       do i=1,n
!         x(i) = alpha*i - alpha/2
!       enddo
!       x(n+1) = x(n) + alpha/2
!       ! Calculate physical grid
!       nextr = minr
!       allinrange = .false.
!       do i=1,n+1
!         nextr_ode(1) = nextr
!         call odeint(nextr_ode,1,x(i-1),x(i),eps,alpha/10.d0,0.d0,nok,nbad,der,rkqs)
!         nextr = nextr_ode(1)
!
!         if(i.eq.n+1.and.nextr.le.maxr)then
!           allinrange = .true.
!         else
!           if(nextr.gt.maxr)exit
!           grid(i) = nextr
!           call der(0.0d0,nextr,jac(i))
!         endif
!       enddo
!       if(allinrange)then
!         print *,alpha,'Before'
!       else
!         print *,alpha,'After'
!       endif
!       ! Adjust alpha
!       if(veryclose)then
!         if(allinrange)then
!           alpha1 = alpha
!         else
!           alpha2 = alpha
!         endif
!         alpha = (alpha1 + alpha2) / 2
!         if(abs(alpha1-alpha2).le.eps)exit
!       else
!         if(.not.allinrange)then
!           veryclose = .true.
!           alpha2 = alpha
!           alpha = (alpha1 + alpha2) / 2
!         else
!           alpha1 = alpha
!           alpha2 = alpha
!           alpha = alpha * 1.01
!         endif
!       endif
!     enddo
!   end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates an optimized grid from *from* to *to*.
! Step specifies step on grid pre-image. Step on the actual grid is variable. Number of points is adjusted.
! *start* is a point between *from* and *to* from where generation begins in left and right direction.
! This allows for symmetric distribution of points if *deriv* is symmetric.
! *deriv* evaluates the right-hand side of the differential equation that defines the optimized grid.
! Resulting optimized grid is returned in *grid*. The corresponding values of *deriv* at grid points are returned in *jac*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_optgrid_step(from, to, start, step, deriv, grid, jac)
    real(real64), intent(in) :: from, to, start, step
    procedure(diff_equations_rhs) :: deriv
    real(real64), allocatable, intent(out) :: grid(:), jac(:)
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
      if (next_point < from .and. mod(grid_right % get_size() + grid_left % get_size(), 2) == 0) then
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
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates grid with required number of points and boundaries. Step size is adjusted.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_equidistant_grid_points(minr, maxr, n, grid, jac, step)
    real(real64), intent(in) :: minr, maxr
    integer, intent(in) :: n
    real(real64), allocatable, intent(out) :: grid(:), jac(:)
    real(real64), intent(out) :: step
    integer :: i
    
    allocate(grid(n), jac(n))
    do i = 1, n
      grid(i) = minr + (maxr - minr) * (i - 0.5d0) / n
    end do
    jac = 1d0
    step = grid(2) - grid(1)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates grid with required step size and boundaries. Number of points is adjusted.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine generate_equidistant_grid_step(minr, maxr, step, grid, jac, n)
    real(real64), intent(in) :: minr, maxr, step
    real(real64), allocatable, intent(out) :: grid(:), jac(:)
    integer, intent(out) :: n
    
    grid = generate_real_range(minr + step / 2d0, maxr - step / 2d0, step)
    n = size(grid)
    allocate(jac(n))
    jac = 1d0
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Saves grid, Jacobian and alpha to a file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_grid(grid, jac, alpha, file_name)
    real(real64), intent(in) :: grid(:), jac(:)
    character(*), intent(in) :: file_name
    real(real64), intent(in) :: alpha
    integer :: i, file_unit
    
    call assert(size(grid) == size(jac), 'Error: grid and jac have to have the same size')
    open(newunit = file_unit, file = file_name)
    write(file_unit, *) size(grid), alpha
    do i = 1, size(grid)
      write(file_unit, '(2G23.15)') grid(i), jac(i)
    enddo
    close(file_unit)
  end subroutine

end module
