!-------------------------------------------------------------------------------------------------------------------------------------------
! Data and procedures related to complex absorbing potential (CAP).
!-------------------------------------------------------------------------------------------------------------------------------------------
module cap_mod
  use constants
  use formulas_mod, only: get_reduced_mass
  use general_utils
  use input_params_mod
  use iso_fortran_env, only: real64
  implicit none

  private :: cap
  real(real64), allocatable :: cap(:)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates Manolopoulos CAP on a given *grid_rho*.
! See D.E. Manolopoulos, J. Chem. Phys. 117, 9552 (2002) for more details.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_cap(params, grid_rho)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: grid_rho(:)
    integer :: i
    real(real64), parameter :: a = 0.112449d0
    real(real64), parameter :: b = 0.00828735d0
    real(real64), parameter :: c = 2.62206d0
    real(real64) :: delta, mu, kmin, absorbing_width, r1, rho, x

    ! Allocate array
    allocate(cap(size(grid_rho)))
    cap = 0

    ! Setup parameters for Manolopoulos CAP
    delta = c / (4 * pi)
    mu = get_reduced_mass(params % mass)
    kmin = sqrt(2 * mu * params % cap % min_absorbed_energy)
    absorbing_width = c / (2 * delta * kmin)
    r1 = params % grid_rho % to - absorbing_width

    ! Fill in Manolopoulos CAP
    do i = 1, size(grid_rho)
      rho = grid_rho(i)
      if (rho > r1) then
        x = 2 * delta * kmin * (rho - r1)
        cap(i) = params % cap % min_absorbed_energy * (a*x - b*x**3 + 4/(c-x)**2 - 4/(c+x)**2)
      end if
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes CAPs (complex absorbing potential).
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_caps(params, grid_rho)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: grid_rho(:)

    if (params % cap % type == 'none') then
      return
    else
      call calc_cap(params, grid_rho)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns CAP/-i as real(real64).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_real_cap() result(cap_copy)
    real(real64), allocatable :: cap_copy(:)
    call assert(allocated(cap), 'Error: cap is not initialized')
    cap_copy = cap
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns CAP as complex(real64).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_complex_cap() result(cap_complex)
    complex(real64), allocatable :: cap_complex(:)
    cap_complex = (0, -1d0) * get_real_cap()
  end function

end module
