module rovib_utils_mod
  use input_params_mod
  use rovib_utils_base_mod
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes index of given value of K.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_ind(K, K_start) result(K_ind)
    integer, intent(in) :: K, K_start
    integer :: K_ind
    K_ind = K - K_start + 1
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes symmetry of given K block based on symmetry of K=0 block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_symmetry(K, K0_sym) result(symmetry)
    integer, intent(in) :: K, K0_sym
    integer :: symmetry
    symmetry = mod(K + K0_sym, 2)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes K based on starting K and K-block index.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function k_ind_to_k(K_ind, K_start) result(K)
    integer, intent(in) :: K_ind, K_start
    integer :: K
    K = K_start + K_ind - 1
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates common attributes of given *K* value.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_k_attributes(K, params, K_ind, K_sym, K_ind_comp)
    integer, intent(in) :: K
    class(input_params), intent(in) :: params
    integer, optional, intent(out) :: K_ind, K_sym, K_ind_comp

    K_ind = get_k_ind(K, params % K(1))
    K_sym = get_k_symmetry(K, params % symmetry)
    K_ind_comp = merge(K_sym + 1, K_ind, params % fixed_basis % enabled == 1) ! compressed index for possibly compressed structures
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes m value corresponding to given m index.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function m_ind_to_m(m_ind, K, K0_sym) result(m)
    integer, intent(in) :: m_ind, K, K0_sym
    integer :: m
    m = m_ind - 1 + get_k_symmetry(K, K0_sym)
  end function
end module
