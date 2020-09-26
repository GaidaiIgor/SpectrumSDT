module general_real_mod
#include "type_list.macro"
#define TYPE_ID REAL_ID
#include "type_attributes.macro"
  use general_char_str_mod
  use iso_fortran_env, only: real64
#include "general_template.F90"

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts real(real64) to string
!-------------------------------------------------------------------------------------------------------------------------------------------
  function num2str_real(num, format) result(str)
    real(real64), intent(in) :: num
    character(*), intent(in), optional :: format
    character(:), allocatable :: str
    character(256) :: buffer
    character(:), allocatable :: format_act

    format_act = arg_or_default_char_str(format, '(G10.2)')
    write(buffer, format_act) num
    str = trim(adjustl(buffer))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to real(real64)
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function str2real(str) result(real_num)
    character(*), intent(in) :: str
    real(real64) :: real_num
    read(str, *) real_num
  end function

end module
