module algorithms_integer_mod
#include "type_list.macro"
#define TYPE_ID INTEGER_ID
#include "type_attributes.macro"
#include "algorithms_template.F90"

!-------------------------------------------------------------------------------------------------------------------------------------------
! Finds index of all *values* in *array*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function findloc_all_integer(array, values) result(inds)
    integer, intent(in) :: array(:), values(:)
    integer :: inds(size(values))
    integer :: i

    do i = 1, size(values)
      inds(i) = findloc(array, values(i), dim = 1)
    end do
  end function

end module
