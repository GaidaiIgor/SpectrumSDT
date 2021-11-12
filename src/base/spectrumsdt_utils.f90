!-------------------------------------------------------------------------------------------------------------------------------------------
! Miscellaneous procedures specific to spectrumsdt.
!-------------------------------------------------------------------------------------------------------------------------------------------
module spectrumsdt_utils_mod
  use constants_mod
  use general_utils_mod
  use iso_fortran_env, only: real64
  implicit none

  interface get_m_ind_info
    module procedure :: get_m_ind_info_plain
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
! Calculates reduced mass of given 3 *masses*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_reduced_mass(masses) result(mu)
    real(real64), intent(in) :: masses(3)
    real(real64) :: mu
    mu = sqrt(product(masses) / sum(masses))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns molecule type: AAA or ABA.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_molecule_type(masses) result(molecule_type)
    real(real64), intent(in) :: masses(3)
    character(3) :: molecule_type

    if ((masses(1) .aeq. masses(2)) .and. (masses(1) .aeq. masses(3))) then
      molecule_type = 'AAA'
    else
      molecule_type = 'ABA'
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a multiplier, relating number of functions per symmetry block to number of functions per basis type (sines or cosines).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_basis_type_multiplier(K0_sym, molecule_type) result(basis_type_multiplier)
    integer, intent(in) :: K0_sym
    character(*), intent(in) :: molecule_type
    integer :: basis_type_multiplier

    if (molecule_type == 'AAA') then
      if (any(K0_sym == [0, 1])) then
        basis_type_multiplier = 1
      else if (any(K0_sym == [2, 3])) then
        basis_type_multiplier = 2
      end if

    else if (molecule_type == 'ABA') then
      basis_type_multiplier = 3
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a multiplier, relating number of functions per basis type to total number of functions.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_basis_type_to_total_multiplier(use_geometric_phase) result(basis_type_to_total_multiplier)
    integer, intent(in) :: use_geometric_phase
    integer :: basis_type_to_total_multiplier
    basis_type_to_total_multiplier = iff(use_geometric_phase == 0, 1, 2)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a multiplier, relating number of functions per symmetry block to total number of functions.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_total_basis_multiplier(sym, molecule_type, use_geometric_phase) result(total_multiplier)
    integer, intent(in) :: sym, use_geometric_phase
    character(*), intent(in) :: molecule_type
    integer :: total_multiplier
    integer :: basis_type_multiplier, basis_type_to_total_multiplier

    basis_type_multiplier = get_basis_type_multiplier(sym, molecule_type)
    basis_type_to_total_multiplier = get_basis_type_to_total_multiplier(use_geometric_phase)
    total_multiplier = basis_type_multiplier * basis_type_to_total_multiplier
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns information associated with a given *m_ind*. m_type = 0 for cos, 1 for sin.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_m_ind_info_plain(m_ind, nphi_per_basis_type, K_sym, molecule_type, m, m_type)
    integer, intent(in) :: m_ind, nphi_per_basis_type, K_sym
    character(*), intent(in) :: molecule_type
    integer, optional, intent(out) :: m, m_type
    integer :: m_ind_sym, m_act, m_type_act

    m_type = iff(m_ind <= nphi_per_basis_type .and. any(K_sym == [0, 2]), 0, 1)
    m_ind_sym = iff(m_ind > nphi_per_basis_type, m_ind - nphi_per_basis_type, m_ind)
    m_ind_sym = iff(m_type == 0 .and. K_sym == 0, m_ind_sym - 1, m_ind_sym)

    if (molecule_type == 'AAA') then
      if (any(K_sym == [0, 1])) then
        m_act = 3 * m_ind_sym
      else if (any(K_sym == [2, 3])) then
        m_act = m_ind_sym + int((m_ind_sym - 1) / 2)
      end if
    else if (molecule_type == 'ABA') then
      m_act = m_ind_sym
    end if

    if (present(m)) then
      m = m_act
    end if
    if (present(m_type)) then
      m_type = m_type_act
    end if
  end subroutine

end module
