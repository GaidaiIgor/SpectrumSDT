module general_utils_integer_mod
#include "type_list.macro"
#define TYPE_ID INTEGER_ID
#include "type_attributes.macro"
  use general_utils_char_str_mod
#include "general_utils_template.F90"

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts integer to string.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function num2str_integer(num, format) result(str)
    integer, intent(in) :: num
    character(*), optional, intent(in) :: format
    character(:), allocatable :: str
    character(256) :: buffer
    character(:), allocatable :: format_act

    format_act = arg_or_default_char_str(format, '(I0)')
    write(buffer, format_act) num
    str = trim(adjustl(buffer))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts integer array to string.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function arr2str_integer(arr, format) result(str)
    integer, intent(in) :: arr(:)
    character(*), optional, intent(in) :: format
    character(:), allocatable :: str
    integer :: i
    character(:), allocatable :: next_num

    str = num2str_integer(arr(1), format)
    do i = 2, size(arr)
      next_num = num2str_integer(arr(i), format)
      str = str // ', ' // next_num
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to integer.
!-------------------------------------------------------------------------------------------------------------------------------------------
  impure elemental function str2int(str, iostat) result(int)
    character(*), intent(in) :: str
    integer, optional, intent(inout) :: iostat
    integer :: int

    if (present(iostat)) then
      read(str, *, iostat = iostat) int
    else
      read(str, *) int
    end if
  end function

end module
