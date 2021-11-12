module rovib_utils_mod
  use general_utils_mod
  use input_params_mod
  use rovib_utils_base_mod
  implicit none

  interface get_k_symmetry
    module procedure :: get_k_symmetry_plain, get_k_symmetry_params
  end interface

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
! Computes K based on starting K and K-block index.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function k_ind_to_k(K_ind, K_start) result(K)
    integer, intent(in) :: K_ind, K_start
    integer :: K
    K = K_start + K_ind - 1
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes symmetry of given K block based on symmetry of K=0 block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_symmetry_plain(K, K0_sym, use_geometric_phase) result(symmetry)
    integer, intent(in) :: K, K0_sym, use_geometric_phase
    integer :: symmetry

    if (use_geometric_phase == 1) then
      symmetry = K0_sym
    else if (any(K0_sym == [0, 1])) then
      symmetry = mod(K + K0_sym, 2)
    else if (any(K0_sym == [2, 3])) then
      symmetry = mod(K + K0_sym - 2, 2) + 2
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes symmetry of given K block based on symmetry of K=0 block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_symmetry_params(K, params) result(symmetry)
    integer, intent(in) :: K
    class(input_params), intent(in) :: params
    integer :: symmetry
    symmetry = get_k_symmetry(K, params % basis % K0_symmetry, params % use_geometric_phase)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns regular K index for adiabatic basis and compressed K index for fixed basis.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_ind_smart(K, params) result(K_ind)
    integer, intent(in) :: K
    class(input_params), intent(in) :: params
    integer :: K_ind

    if (params % basis % fixed % enabled == 0) then
      K_ind = get_k_ind(K, params % K(1))
    else if (params % use_geometric_phase == 0) then
      K_ind = iff(mod(K - params % K(1), 2) == 0, 1, 2)
    else
      K_ind = 1
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates common attributes of given *K* value.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine get_k_attributes(K, params, K_ind, K_sym, K_ind_smart)
    integer, intent(in) :: K
    class(input_params), intent(in) :: params
    integer, optional, intent(out) :: K_ind, K_sym, K_ind_smart
    integer :: K_ind_act, K_sym_act

    K_ind_act = get_k_ind(K, params % K(1))
    if (present(K_ind)) then
      K_ind = K_ind_act
    end if

    K_sym_act = get_k_symmetry(K, params % basis % K0_symmetry, params % use_geometric_phase)
    if (present(K_sym)) then
      K_sym = K_sym_act
    end if

    if (present(K_ind_smart)) then
      K_ind_smart = get_k_ind_smart(K, params)
    end if
  end subroutine

end module
