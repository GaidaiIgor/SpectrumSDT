!-------------------------------------------------------------------------------------------------------------------------------------------
! Contains historically global variables shared by some other modules.
!-------------------------------------------------------------------------------------------------------------------------------------------
module general_vars
  use input_params_mod
  use iso_fortran_env, only: real64
  use path_utils
  use spectrumsdt_paths_mod
  implicit none

  real(real64), allocatable :: g1(:), g2(:), g3(:)
  real(real64), allocatable :: jac1(:), jac2(:), jac3(:)
  real(real64) :: alpha1, alpha2, alpha3

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads grids, jacobians, alphas (steps). Inits corresponding parameters in *params*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_grids(params)
    class(input_params), intent(inout) :: params
    integer :: i, file_unit
    real(real64) :: skip
    
    open(newunit = file_unit, file = get_grid_rho_path(params))
    read(file_unit, *) params % grid_rho % from, params % grid_rho % to, params % grid_rho % step, params % grid_rho % num_points
    alpha1 = params % grid_rho % step
    allocate(g1(params % grid_rho % num_points), jac1(params % grid_rho % num_points))
    do i = 1, size(g1)
      read(file_unit, *) g1(i), jac1(i)
    end do
    close(file_unit)
    
    open(newunit = file_unit, file = get_grid_theta_path(params))
    read(file_unit, *) params % grid_theta % from, params % grid_theta % to, params % grid_theta % step, params % grid_theta % num_points
    alpha2 = params % grid_theta % step
    allocate(g2(params % grid_theta % num_points), jac2(params % grid_theta % num_points))
    read(file_unit, *) g2
    jac2 = 1
    close(file_unit)
    
    open(newunit = file_unit, file = get_grid_phi_path(params))
    read(file_unit, *) skip, skip, alpha3, params % num_points_phi
    allocate(g3(params % num_points_phi), jac3(params % num_points_phi))
    read(file_unit, *) g3
    jac3 = 1
    close(file_unit)
    
    ! Treat alpha2(3) as a grid step size
    alpha2 = alpha2 * jac2(1)
    alpha3 = alpha3 * jac3(1)
  end subroutine

end module
