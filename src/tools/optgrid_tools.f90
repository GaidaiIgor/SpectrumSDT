module optgrid_tools
  use general_utils
  use numerical_recipies
contains
  !-----------------------------------------------------------------------
  !  OptGrid (optgrid.f)
  !  Optimal grid generator, subroutines:
  !    1. (r,radial) Left end fixed, right end free, alpha is parameter
  !    2. (a,alpha)  Both ends are fixed, alpha is adjusted
  !    3. (n,number) Both ends are fixed, number of points is adjusted
  !  Author: Alexander Teplukhin
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !  Generates a grid with fixed left end and free right end.
  !  Alpha is fixed.
  !-----------------------------------------------------------------------
  subroutine generate_gridr(grid,jac,n,alpha,minr,eps,der)
    implicit none
    real*8 grid(n),jac(n),x(0:n)
    real*8 alpha,minr,eps,nextr
    real*8 :: nextr_ode(1)
    integer i,n,nok,nbad
    external der

    x(0) = 0
    do i=1,n
      x(i) = alpha*i - alpha/2
    enddo
    nextr = minr
    do i=1,n
      nextr_ode(1) = nextr
      call odeint(nextr_ode,1,x(i-1),x(i),eps,alpha/10.d0,0.d0,nok,nbad,der,rkqs)
      nextr = nextr_ode(1)
      
      grid(i) = nextr
      call der(0.0d0,nextr,jac(i))
    enddo
  end subroutine

  !-----------------------------------------------------------------------
  !  Generates a grid within required boundaries.
  !  Solves from the left to the right.
  !  Alpha is adjusted.
  !-----------------------------------------------------------------------
  subroutine generate_grida(grid,jac,n,alpha,minr,maxr,eps,der)
    implicit none
    real*8 grid(n),jac(n),x(0:n+1)
    real*8 alpha,alpha1,alpha2,minr,maxr,nextr,eps
    real*8 :: nextr_ode(1)
    integer i,n,nok,nbad
    logical veryclose,allinrange
    external der

    ! Iteratively change alpha and regenerate the grid
    veryclose = .false.
    alpha1 = 0
    do
      ! Calculate working grid
      x(0) = 0
      do i=1,n
        x(i) = alpha*i - alpha/2
      enddo
      x(n+1) = x(n) + alpha/2
      ! Calculate physical grid
      nextr = minr
      allinrange = .false.
      do i=1,n+1
        nextr_ode(1) = nextr
        call odeint(nextr_ode,1,x(i-1),x(i),eps,alpha/10.d0,0.d0,nok,nbad,der,rkqs)
        nextr = nextr_ode(1)
        
        if(i.eq.n+1.and.nextr.le.maxr)then
          allinrange = .true.
        else
          if(nextr.gt.maxr)exit
          grid(i) = nextr
          call der(0.0d0,nextr,jac(i))
        endif
      enddo
      if(allinrange)then
        print *,alpha,'Before'
      else
        print *,alpha,'After'
      endif
      ! Adjust alpha
      if(veryclose)then
        if(allinrange)then
          alpha1 = alpha
        else
          alpha2 = alpha
        endif
        alpha = (alpha1 + alpha2) / 2
        if(abs(alpha1-alpha2).le.eps)exit
      else
        if(.not.allinrange)then
          veryclose = .true.
          alpha2 = alpha
          alpha = (alpha1 + alpha2) / 2
        else
          alpha1 = alpha
          alpha2 = alpha
          alpha = alpha * 1.01
        endif
      endif
    enddo
  end subroutine

  !-----------------------------------------------------------------------
  !  Generates a grid within required boundaries.
  !  Number of points is adjusted.
  !  Makes symmetric distribution of points for symmetric derivative.
  !-----------------------------------------------------------------------
  subroutine generate_gridn(grid,jac,n,alpha,minr,maxr,r0,eps,der)
    implicit none
    integer,parameter:: n0 = 1024
    real*8,allocatable::grid(:),jac(:)
    real*8 gridl(n0),gridr(n0),x(0:n0)
    real*8 r0,alpha,minr,maxr,nextr,eps
    real*8 :: nextr_rkdumb(1)
    integer i,n,nl,nr
    external der

    x(0) = 0
    do i=1,n0
      x(i) = alpha*i - alpha/2
    enddo
    nextr = r0
    do i=1,n0
      ! call odeint(nextr,1,x(i-1),x(i),eps,alpha/10.d0,0.d0,nok,nbad,der,rkqs)
      nextr_rkdumb(1) = nextr
      call rkdumb(nextr_rkdumb, 1, x(i-1), x(i), 2048, der, nextr_rkdumb)
      nextr = nextr_rkdumb(1)

      if(nextr>maxr)exit
      gridr(i) = nextr
    enddo
    nr = i-1
    x = - x
    nextr = r0
    do i=1,n0
      ! call odeint(nextr,1,x(i-1),x(i),eps,alpha/10.d0,0.d0,nok,nbad,der,rkqs)
      nextr_rkdumb(1) = nextr
      call rkdumb(nextr_rkdumb, 1, x(i-1), x(i), 2048, der, nextr_rkdumb)
      nextr = nextr_rkdumb(1)
      
      if(nextr<minr)then
        if(mod(nr+i-1,2)==1)then
          gridl(i) = nextr
          nl = i
        else
          nl = i - 1
        endif
        exit
      endif
      gridl(i) = nextr
    enddo
    n = nl + nr
    allocate(grid(n),jac(n))
    do i=1,n
      if(i<=nl)then
        grid(i) = gridl(nl-i+1)
      else
        grid(i) = gridr(i-nl)
      endif
      call der(0.0d0,grid(i),jac(i))
    enddo
  end subroutine

  !-----------------------------------------------------------------------
  ! Generates grid with required number of points and boundaries. Step size is adjusted
  !-----------------------------------------------------------------------
  subroutine generate_equidistant_grid_points(minr, maxr, n, grid, jac, step)
    implicit none
    real*8, intent(in) :: minr, maxr
    integer, intent(in) :: n
    real*8, allocatable, intent(out) :: grid(:), jac(:)
    real*8, intent(out) :: step
    integer :: i
    
    allocate(grid(n), jac(n))
    do i = 1, n
      grid(i) = minr + (maxr - minr) * (i - 0.5d0) / n
    end do
    jac = 1d0
    step = grid(2) - grid(1)
  end subroutine

  !-----------------------------------------------------------------------
  ! Generates grid with required step size and boundaries. Number of points is adjusted
  !-----------------------------------------------------------------------
  subroutine generate_equidistant_grid_step(minr, maxr, step, grid, jac, n)
    implicit none
    real*8, intent(in) :: minr, maxr, step
    real*8, allocatable, intent(out) :: grid(:), jac(:)
    integer, intent(out) :: n
    
    grid = generate_real_range(minr + step / 2d0, maxr - step / 2d0, step)
    n = size(grid)
    allocate(jac(n))
    jac = 1d0
  end subroutine

  !-----------------------------------------------------------------------
  !  Print Grid
  !  Saves grid, Jacobian and potential on the grid to the *.dat file
  !  Corresponding alpha is saved too
  !-----------------------------------------------------------------------
  subroutine print_grid(grid, jac, n, name, alpha, potvib)
    use constants
    implicit none
    real*8, intent(in) :: grid(n), jac(n)
    integer, intent(in) :: n
    character(*), intent(in) :: name
    real*8, intent(in) :: alpha
    real*8, external :: potvib
    integer :: i
    
    ! external potvib
    open(1,file=name)
    write(1,*)n,alpha
    do i=1,n
      write(1,'(3F30.17)')grid(i),jac(i),potvib(grid(i))*autown
    enddo
    close(1)
  end subroutine
end module
