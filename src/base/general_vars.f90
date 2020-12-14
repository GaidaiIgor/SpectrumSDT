!-------------------------------------------------------------------------------------------------------------------------------------------
! Contains historically global variables shared by some other modules.
!-------------------------------------------------------------------------------------------------------------------------------------------
module general_vars
  use constants, only: oxygen_masses
  use input_params_mod
  use iso_fortran_env, only: real64
  use path_utils
  use spectrumsdt_paths_mod
  implicit none

  ! Masses
  real(real64) :: m0          ! Mass of central atom
  real(real64) :: m1          ! Mass of 1st terminal atom
  real(real64) :: m2          ! Mass of 2nd terminal atom
  real(real64) :: mtot        ! Total mass
  real(real64) :: mu          ! Three-particle reduced mass
  real(real64) :: mu0         ! Reduced mass of m1 and m2
  real(real64) :: mu1         ! Reduced mass of m0 and m2
  real(real64) :: mu2         ! Reduced mass of m0 and m1

  ! Grids
  integer :: n1, n2, n3 ! Grid sizes
  real(real64), allocatable :: g1(:), g2(:), g3(:)
  real(real64), allocatable :: jac1(:), jac2(:), jac3(:)
  real(real64) :: alpha1, alpha2, alpha3

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes masses.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_masses(params)
    class(input_params), intent(in) :: params

    m1 = params % mass(1)
    m0 = params % mass(2)
    m2 = params % mass(3)
    mtot = m0 + m1 + m2
    mu = sqrt(m0 * m1 * m2 / mtot)
    mu0 = m1 * m2 / (m1 + m2)
    mu1 = m0 * m2 / (m0 + m2)
    mu2 = m0 * m1 / (m0 + m1)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads grids, jacobians, alphas (steps). Inits corresponding parameters in *params*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_grids(params)
    class(input_params), intent(inout) :: params
    integer :: i, file_unit
    real(real64) :: skip
    
    open(newunit = file_unit, file = get_grid_rho_path(params))
    read(file_unit, *) params % grid_rho % from, params % grid_rho % to, params % grid_rho % step, params % grid_rho % num_points
    n1 = params % grid_rho % num_points
    alpha1 = params % grid_rho % step
    allocate(g1(n1), jac1(n1))
    do i = 1, n1
      read(file_unit, *) g1(i), jac1(i)
    end do
    close(file_unit)
    
    open(newunit = file_unit, file = get_grid_theta_path(params))
    read(file_unit, *) params % grid_theta % from, params % grid_theta % to, params % grid_theta % step, params % grid_theta % num_points
    n2 = params % grid_theta % num_points
    alpha2 = params % grid_theta % step
    allocate(g2(n2), jac2(n2))
    do i = 1, n2
      read(file_unit, *) g2(i), jac2(i)
    end do
    close(file_unit)
    
    open(newunit = file_unit, file = get_grid_phi_path(params))
    read(file_unit, *) skip, skip, alpha3, params % num_points_phi
    n3 = params % num_points_phi
    allocate(g3(n3), jac3(n3))
    do i = 1, n3
      read(file_unit, *) g3(i), jac3(i)
    end do
    close(file_unit)
    
    ! Treat alpha2(3) as a grid step size
    alpha2 = alpha2 * jac2(1)
    alpha3 = alpha3 * jac3(1)
  end subroutine

end module
