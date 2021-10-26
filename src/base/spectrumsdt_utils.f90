!-------------------------------------------------------------------------------------------------------------------------------------------
! Miscellaneous procedures specific to spectrumsdt.
!-------------------------------------------------------------------------------------------------------------------------------------------
module spectrumsdt_utils_mod
  use constants_mod
  use general_utils_mod
  use input_params_mod
  use iso_fortran_env, only: real64
  implicit none

  interface get_num_funcs_phi
    module procedure :: get_num_funcs_phi_plain, get_num_funcs_phi_params
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates the value of rotational constant A, for given reduced mass and geometry.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rotational_a(mu, rho, theta) result(a)
    real(real64), intent(in) :: mu, rho, theta
    real(real64) :: a
    a = 1d0 / (mu * rho**2 * (1+sin(theta)))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates the value of rotational constant B, for given reduced mass and geometry.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rotational_b(mu, rho, theta) result(b)
    real(real64), intent(in) :: mu, rho, theta
    real(real64) :: b
    b = 1d0 / (2d0 * mu * rho**2 * sin(theta)**2)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates the value of rotational constant C, for given reduced mass and geometry.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rotational_c(mu, rho, theta) result(c)
    real(real64), intent(in) :: mu, rho, theta
    real(real64) :: c
    c = 1d0 / (mu * rho**2 * (1-sin(theta)))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns information associated with a given *m_ind*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_m_ind_info(m_ind, phi_sym, num_funcs_phi_per_sym, m, m_sym)
    integer, intent(in) :: m_ind, phi_sym
    integer, optional, intent(in) :: num_funcs_phi_per_sym
    integer, optional, intent(out) :: m, m_sym
    integer :: m_act, m_sym_act

    if (any(phi_sym == [0, 1])) then
      m_act = iff(phi_sym == 0, m_ind - 1, m_ind)
      m_sym_act = phi_sym
    else if (phi_sym == 2) then
      call assert(present(num_funcs_phi_per_sym), 'Error: num_funcs_phi_per_sym has to be supplied if phi_sym == 2')
      if (m_ind <= num_funcs_phi_per_sym) then
        m_act = m_ind - 1
        m_sym_act = 0
      else
        m_act = m_ind - num_funcs_phi_per_sym
        m_sym_act = 1
      end if
    else
      stop 'Unknown phi_sym'
    end if

    if (present(m)) then
      m = m_act
    end if
    if (present(m_sym)) then
      m_sym = m_sym_act
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates reduced mass of given 3 *masses*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_reduced_mass(masses) result(mu)
    real(real64), intent(in) :: masses(3)
    real(real64) :: mu
    mu = sqrt(product(masses) / sum(masses))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns total number of phi functions.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_num_funcs_phi_plain(num_funcs_phi_per_sym, sym) result(num_funcs_phi)
    integer, intent(in) :: num_funcs_phi_per_sym, sym
    integer :: num_funcs_phi

    if (sym == 2) then
      num_funcs_phi = 2 * num_funcs_phi_per_sym
    else
      num_funcs_phi = num_funcs_phi_per_sym
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns total number of phi functions.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_num_funcs_phi_params(params) result(num_funcs_phi)
    class(input_params), intent(in) :: params
    integer :: num_funcs_phi
    num_funcs_phi = get_num_funcs_phi_plain(params % basis % num_funcs_phi_per_sym, params % basis % symmetry)
  end function

end module
