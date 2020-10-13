!-------------------------------------------------------------------------------------------------------------------------------------------
! Contains functions that evaluate simple physical formulas
!-------------------------------------------------------------------------------------------------------------------------------------------
module formulas_mod
  use constants
  use general_utils
  use iso_fortran_env, only: real64
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates the value of rotational constant A, for given ozone reduced mass and geometry
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rotational_a(mu, rho, theta) result(a)
    real(real64), intent(in) :: mu, rho, theta
    real(real64) :: a
    a = 1d0 / (mu * rho**2 * (1+sin(theta)))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates the value of rotational constant B, for given ozone reduced mass and geometry
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rotational_b(mu, rho, theta) result(b)
    real(real64), intent(in) :: mu, rho, theta
    real(real64) :: b
    b = 1d0 / (2d0 * mu * rho**2 * sin(theta)**2)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates the value of rotational constant C, for given ozone reduced mass and geometry
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rotational_c(mu, rho, theta) result(c)
    real(real64), intent(in) :: mu, rho, theta
    real(real64) :: c
    c = 1d0 / (mu * rho**2 * (1-sin(theta)))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates the value of lambda plus function
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_lambda_plus(J, K) result(res)
    integer, intent(in) :: J, K
    real(real64) :: res
    res = sqrt(1d0 * (J + K + 1) * (J - K))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates the value of U factor for asymmetric term
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_U(J, K1, K2, parity) result(U)
    integer, intent(in) :: J, K1, K2, parity
    real(real64) :: U

    call assert(K1 + 2 == K2 .or. K1 - 2 == K2 .or. K1 == 1 .and. K2 == 1, 'Wrong input for U')
    U = 0
    if (K1 + 2 == K2) then ! above diagonal
      U = calculate_lambda_plus(J, K1) * calculate_lambda_plus(J, K1 + 1)
    else if (K1 - 2 == K2) then ! below diagonal
      U = calculate_lambda_plus(J, K1 - 1) * calculate_lambda_plus(J, K1 - 2)
    end if

    ! Extra diagonal across
    if (2 - K1 == K2) then
      U = U + (-1) ** (J + K1 + parity) * calculate_lambda_plus(J, K1 - 1) * calculate_lambda_plus(J, K1 - 2)
    end if
    U = U / 2d0

    ! Delta-term for K=0 case
    if (K1 == 0 .or. K2 == 0) then
      U = U / sqrt(2d0)
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates the value of W factor for coriolis term
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_W(J, K1, K2, parity) result(W)
    integer, intent(in) :: J, K1, K2, parity
    real(real64) :: W

    call assert(K1 + 1 == K2 .or. K1 - 1 == K2, 'Wrong input for W')
    if (K1 + 1 == K2) then ! above diagonal
      W = calculate_lambda_plus(J, K1)
    else if (K1 - 1 == K2) then ! below diagonal
      W = -calculate_lambda_plus(J, K1 - 1)
    end if

    if (1 - K1 == K2) then ! (1, 0) or (0, 1)
      W = W + (-1) ** (J + K1 + parity) * calculate_lambda_plus(J, K1 - 1)
    end if
    W = W / 2d0

    ! Delta-term for K=0 case
    if (K1 == 0 .or. K2 == 0) then
      W = W / sqrt(2d0)
    end if
  end function
end module
