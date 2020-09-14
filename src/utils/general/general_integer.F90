module general_integer_mod
#include "type_list.macro"
#define TYPE_ID INTEGER_ID
#include "type_attributes.macro"

use general_char_str_mod
implicit none

contains

#include "general_template.F90"

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts integer to string
!-------------------------------------------------------------------------------------------------------------------------------------------
  function num2str_integer(num, format) result(str)
    integer, intent(in) :: num
    character(*), intent(in), optional :: format
    character(:), allocatable :: str
    character(256) :: buffer
    character(:), allocatable :: format_act

    format_act = arg_or_default_char_str(format, '(I0)')
    write(buffer, format_act) num
    str = trim(adjustl(buffer))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to integer
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function str2int(str) result(int)
    character(*), intent(in) :: str
    integer :: int
    read(str, *) int
  end function

end module
