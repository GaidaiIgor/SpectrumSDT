!-------------------------------------------------------------------------------------------------------------------------------------------
! Data and procedures related to complex absorbing potential (CAP).
!-------------------------------------------------------------------------------------------------------------------------------------------
module cap_mod
  use constants
  use general_vars
  use input_params_mod
  use iso_fortran_env, only: real64
  implicit none

  private :: cap
  real(real64), allocatable :: cap(:)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates CAP on a given grid for rho.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_cap(params, grid_rho)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: grid_rho(:)
    integer :: i
    real(real64), parameter :: aM = 0.112449d0
    real(real64), parameter :: bM = 0.00828735d0
    real(real64), parameter :: cM = 2.62206d0
    real(real64) :: rc, dM, xM, rho, damp_len

    ! Allocate array
    allocate(cap(size(grid_rho)))
    cap = 0

    ! Setup parameters for Manolopoulos CAP
    dM = cM / (4 * pi)
    damp_len = cM / (2 * dM * sqrt(2 * mu * params % cap % emin))
    rc = params % grid_rho % to - damp_len

    ! Fill in Manolopoulos CAP
    do i = 1, size(grid_rho)
      rho = grid_rho(i)
      if (rho > rc) then
        xM = 2 * dM * sqrt(2 * mu * params % cap % emin) * (rho - rc)
        cap(i) = params % cap % emin * (aM * xM - bM * xM**3 + 4/(cM-xM)**2 - 4/(cM+xM)**2 )
        if (rho > rc + damp_len) then
          cap(i) = cap(i - 1)
        end if
      end if
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes CAPs (complex absorbing potential).
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_caps(params)
    class(input_params), intent(in) :: params
    if (params % cap % type == 'none') then
      return
    else
      call calc_cap(params, g1)
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
