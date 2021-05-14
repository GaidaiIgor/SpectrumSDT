module algorithms_real_mod
#include "type_list.macro"
#define TYPE_ID REAL_ID
#include "type_attributes.macro"
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

!-------------------------------------------------------------------------------------------------------------------------------------------
! Integrates *array* using trapezoidal rule.
! Assumes *array* points are in the middle of intervals, therefore each value corresponds to average trapezoid height and a simple sum can be taken.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function integrate_1d(array, step) result(res)
    real(real64), intent(in) :: array(:)
    real(real64), intent(in) :: step
    real(real64) :: res
    res = sum(array) * step
  end function

end module
