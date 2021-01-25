!-------------------------------------------------------------------------------------------------------------------------------------------
! Data and procedures related to complex absorbing potential (CAP).
!-------------------------------------------------------------------------------------------------------------------------------------------
module cap_mod
  use constants, only: pi
  use formulas_mod, only: get_reduced_mass
  use grid_info_mod
  use input_params_mod
  use iso_fortran_env, only: real64
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates Manolopoulos CAP/-i for a given *rho_info*.
! See D.E. Manolopoulos, J. Chem. Phys. 117, 9552 (2002) for more details.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_real_cap(params, rho_info) result(cap)
    class(input_params), intent(in) :: params
    class(grid_info), intent(in) :: rho_info
    real(real64), allocatable :: cap(:)
    integer :: i
    real(real64), parameter :: a = 0.112449d0
    real(real64), parameter :: b = 0.00828735d0
    real(real64), parameter :: c = 2.62206d0
    real(real64) :: delta, mu, kmin, absorbing_width, rho_start, rho, x

    ! Allocate array
    allocate(cap(size(rho_info % points)))
    cap = 0

    ! Setup parameters for Manolopoulos CAP
    delta = c / (4 * pi)
    mu = get_reduced_mass(params % mass)
    kmin = sqrt(2 * mu * params % cap % min_absorbed_energy)
    absorbing_width = c / (2 * delta * kmin)
    rho_start = rho_info % to - absorbing_width

    ! Fill in Manolopoulos CAP
    do i = 1, size(rho_info % points)
      rho = rho_info % points(i)
      if (rho > rho_start) then
        x = 2 * delta * kmin * (rho - rho_start)
        cap(i) = params % cap % min_absorbed_energy * (a*x - b*x**3 + 4/(c-x)**2 - 4/(c+x)**2)
      end if
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns CAP as complex(real64).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_complex_cap(params, rho_info) result(cap_complex)
    class(input_params), intent(in) :: params
    class(grid_info), intent(in) :: rho_info
    complex(real64), allocatable :: cap_complex(:)
    cap_complex = (0, -1d0) * calc_real_cap(params, rho_info)
  end function

end module
