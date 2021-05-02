!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures from rovib_utils that had to be factored out to avoid circular dependencies.
!-------------------------------------------------------------------------------------------------------------------------------------------
module rovib_utils_base_mod
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes symmetry of given K block based on symmetry of K=0 block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_start(J, parity) result(K_start)
    integer, intent(in) :: J, parity
    integer :: K_start
    K_start = mod(J + parity, 2)
  end function

end module
