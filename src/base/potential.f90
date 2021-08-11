!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures related to potential.
!-------------------------------------------------------------------------------------------------------------------------------------------
module potential_mod
  use input_params_mod
  use iso_fortran_env, only: real64
  use spectrumsdt_paths_mod
  use spectrumsdt_utils_mod
  implicit none

contains

!------------------------------------------------------------------------------------------------------------------------------------------- 
! Calculates rotational potential in symmetric top approximation.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_rotational_potential(mu, rho, theta, J, K) result(res)
    real(real64), intent(in) :: mu, rho, theta
    integer, intent(in) :: J, K
    real(real64) :: res
    real(real64) :: a, b, c

    a = get_rotational_a(mu, rho, theta)
    b = get_rotational_b(mu, rho, theta)
    c = get_rotational_c(mu, rho, theta)
    res = (a + b)/2 * J*(J + 1) + (c - (a + b)/2) * K**2
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates extra potential term at a given point.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_extra_potential(mu, rho, theta) result(res)
    real(real64), intent(in) :: mu, rho, theta
    real(real64) :: res
    res = -(0.25d0 + 4/(sin(2*theta)**2)) / (2*mu*rho**2)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes total potential. Grids have to be initialized first.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function load_potential(params, grid_rho, grid_theta, grid_phi) result(potential)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: grid_rho(:), grid_theta(:), grid_phi(:)
    real(real64), allocatable :: potential(:, :, :)
    integer :: file_unit, iostat, i1, i2, i3
    real(real64) :: mu

    ! Load vibrational potential
    allocate(potential(size(grid_phi), size(grid_theta), size(grid_rho)))
    open(newunit = file_unit, file = get_pes_path(params))
    read(file_unit, *, iostat = iostat) potential
    close(file_unit)
    call assert(.not. is_iostat_end(iostat), 'Error: size of pes.out is not sufficient for the specified grids')

    ! Add rotational and extra potentials
    mu = get_reduced_mass(params % mass)
    do i1 = 1, size(grid_rho)
      do i2 = 1, size(grid_theta)
        do i3 = 1, size(grid_phi)
          potential(i3, i2, i1) = potential(i3, i2, i1) + calc_rotational_potential(mu, grid_rho(i1), grid_theta(i2), params % J, params % K(1)) + &
            calc_extra_potential(mu, grid_rho(i1), grid_theta(i2))
        end do
      end do
    end do
  end function

end module
