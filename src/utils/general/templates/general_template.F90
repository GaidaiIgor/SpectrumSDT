#include "funcs.macro"

!-------------------------------------------------------------------------------------------------------------------------------------------
! returns arg if it's present, otherwise default
!-------------------------------------------------------------------------------------------------------------------------------------------
TEMPLATE_ELEMENTAL function CONCAT2(arg_or_default_, TEMPLATE_TYPE_NAME)(arg, default) result(res)
  TEMPLATE_TYPE, intent(in), optional :: arg
  TEMPLATE_TYPE, intent(in) :: default
  TEMPLATE_TYPE_OUT :: res
  
  if (present(arg)) then
    res = arg
  else
    res = default
  end if
end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Swaps numbers
!-------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(swap_, TEMPLATE_TYPE_NAME)(a, b)
  TEMPLATE_TYPE_OUT, intent(inout) :: a, b
  TEMPLATE_TYPE_OUT :: temp
  
  temp = a
  a = b
  b = temp
end subroutine