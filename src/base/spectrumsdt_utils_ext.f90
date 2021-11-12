!-------------------------------------------------------------------------------------------------------------------------------------------
! Extended version of spectrumsdt_utils with added interfaces for input_params.
!-------------------------------------------------------------------------------------------------------------------------------------------
module spectrumsdt_utils_ext_mod
  use input_params_mod
  use rovib_utils_mod
  use spectrumsdt_utils_mod
  implicit none

  interface get_m_ind_info
    module procedure :: get_m_ind_info_params
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns information associated with a given *m_ind* in first K-block. m_type = 0 for cos, 1 for sin.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_m_ind_info_params(m_ind, params, m, m_type)
    integer, intent(in) :: m_ind
    class(input_params), intent(in) :: params
    integer, optional, intent(out) :: m, m_type
    integer :: K_sym

    K_sym = get_k_symmetry(params % K(1), params)
    call get_m_ind_info(m_ind, params % get_num_funcs_phi_per_basis_type(), K_sym, get_molecule_type(params % mass), m, m_type)
  end subroutine

end module
