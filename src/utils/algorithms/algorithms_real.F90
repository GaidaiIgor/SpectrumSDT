module algorithms_real_mod
#include "type_list.macro"
#define TYPE_ID REAL_ID
#include "type_attributes.macro"
  use iso_fortran_env, only: real64
#include "algorithms_template.F90"

!-------------------------------------------------------------------------------------------------------------------------------------------
! Finds index of all *values* in *array*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function findloc_all_real(array, values) result(inds)
    real(real64), intent(in) :: array(:), values(:)
    integer :: inds(size(values))
    integer :: i, j

    inds = -1
    do i = 1, size(values)
      do j = 1, size(array)
        if (values(i) .aeq. array(j)) then
          inds(i) = j
          exit   
        end if
      end do
    end do
  end function

end module
