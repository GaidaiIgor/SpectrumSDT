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

    if (any(K0_sym == [0, 1])) then
      symmetry = mod(K + K0_sym, 2)
    else if (K0_sym == 2) then
      symmetry = 2
    else
      stop 'Error: unknown symmetry'
    end if
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
    integer :: K_ind_act, K_sym_act

    K_ind_act = get_k_ind(K, params % K(1))
    if (present(K_ind)) then
      K_ind = K_ind_act
    end if

    K_sym_act = get_k_symmetry(K, params % basis % symmetry)
    if (present(K_sym)) then
      K_sym = K_sym_act
    end if

    if (present(K_ind_comp)) then
      K_ind_comp = merge(K_sym_act + 1, K_ind_act, params % basis % fixed % enabled == 1) ! compressed index for possibly compressed structures
    end if
  end subroutine

end module
