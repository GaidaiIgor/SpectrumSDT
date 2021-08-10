module string_mod
  use general_utils_mod
  use iso_fortran_env, only: real64
  implicit none
  
  interface string
    module procedure :: to_string_char_str
  end interface
  
  interface assignment(=)
    module procedure :: assign_char_str, assign_char_str_array
  end interface
  
  type :: string 
    character(:), allocatable :: s
  contains    
    procedure :: write => write_string
    generic :: write(formatted) => write

    procedure :: equals_char_str
    procedure :: equals_string
    generic :: operator(==) => equals_char_str, equals_string
    procedure :: not_equals_char_str
    procedure :: not_equals_string
    generic :: operator(/=) => not_equals_char_str, not_equals_string
    
    procedure :: to_char_str
    procedure :: to_integer
    procedure :: to_real
    procedure :: to_complex
    procedure :: length
    procedure :: substr
    procedure :: trim => trim_string
  end type
  
contains
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes string.
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_string(this, unit, iotype, v_list, iostat, iomsg)
    class(string), intent(in) :: this
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    write(unit, *, iostat = iostat) this % s
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Compares string with char str.
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function equals_char_str(this, char_str) result(res)
    class(string), intent(in) :: this
    character(*), intent(in) :: char_str
    logical :: res
    res = this % s == char_str
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Compares two strings.
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental function equals_string(this, other) result(res)
    class(string), intent(in) :: this, other
    logical :: res
    res = this % s == other % s
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Compares string with char str.
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function not_equals_char_str(this, char_str) result(res)
    class(string), intent(in) :: this
    character(*), intent(in) :: char_str
    logical :: res
    res = this % s /= char_str
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Compares two strings.
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental function not_equals_string(this, other) result(res)
    class(string), intent(in) :: this, other
    logical :: res
    res = this % s /= other % s
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Converts string type to char string.
!---------------------------------------------------------------------------------------------------------------------------------------------
  function to_char_str(this) result(char_str)
    class(string), intent(in) :: this
    character(:), allocatable :: char_str
    char_str = this % s
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to integer.
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function to_integer(this) result(res)
    class(string), intent(in) :: this
    integer :: res
    character(:), allocatable :: char_str
    
    char_str = this % s
    read(char_str, *) res
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to real.
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function to_real(this) result(res)
    class(string), intent(in) :: this
    real(real64) :: res
    character(:), allocatable :: char_str
    
    char_str = this % s
    read(char_str, *) res
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to complex.
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function to_complex(this) result(res)
    class(string), intent(in) :: this
    complex(real64) :: res
    character(:), allocatable :: char_str
    
    char_str = this % s
    read(char_str, *) res
  end function
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Returns length of string.
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental function length(this) result(res)
    class(string), intent(in) :: this
    integer :: res
    res = len(this % s)
  end function
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Takes substring.
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental function substr(this, start, end) result(res)
    class(string), intent(in) :: this
    integer, intent(in), optional :: start, end
    type(string) :: res
    integer :: start_act, end_act
    
    start_act = arg_or_default(start, 1)
    end_act = arg_or_default(end, this % length())
    res = this % s(start_act:end_act)
  end function
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! MODULE PROCEDURES
!---------------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------
! Constructs a string from a char string. Trims by default.
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental function to_string_char_str(char_str, trim_arg) result(str)
    character(*), intent(in) :: char_str
    integer, optional, intent(in) :: trim_arg
    type(string) :: str
    integer :: trim_arg_act
    character(:), allocatable :: char_str_act

    trim_arg_act = arg_or_default(trim_arg, 1)
    char_str_act = iff(trim_arg_act == 1, trim(char_str), char_str)
    str % s = char_str_act
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Asssigns a char string.
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental subroutine assign_char_str(this, char_str)
    type(string), intent(inout) :: this
    character(*), intent(in) :: char_str
    this = string(char_str)
  end subroutine
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Assigns a char string array.
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_char_str_array(array_lhs, array_rhs)
    type(string), allocatable, intent(out) :: array_lhs(:)
    character(*), intent(in) :: array_rhs(:)
    type(string), allocatable :: string_array_rhs(:)
    
    string_array_rhs = string(array_rhs)
    array_lhs = string_array_rhs
  end subroutine

!---------------------------------------------------------------------------------------------------------------------------------------------
! Joins strings from a string array (arr) into a single (char_str). (join_str) specifies what string separates the strings in (arr).
!---------------------------------------------------------------------------------------------------------------------------------------------
  function string_arr_to_char_str(arr, join_str) result(char_str)
    class(string), intent(in) :: arr(:)
    character(*), optional, intent(in) :: join_str
    character(:), allocatable :: char_str
    integer :: i
    character(:), allocatable :: join_str_act

    join_str_act = arg_or_default(join_str, ', ')
    char_str = arr(1) % to_char_str()
    do i = 2, size(arr)
      char_str = char_str // join_str_act // arr(i) % to_char_str()
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string array to char str array.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function string_arr_to_char_str_arr(arr) result(char_str_arr)
    class(string), intent(in) :: arr(:)
    character(:), allocatable :: char_str_arr(:)
    integer :: max_len, i

    max_len = maxval(arr % length())
    allocate(character(max_len) :: char_str_arr(size(arr)))
    do i = 1, size(arr)
      char_str_arr(i) = arr(i) % s
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Trims string, both left and right.
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental subroutine trim_string(this)
    class(string), intent(inout) :: this
    this % s = trim(adjustl(this % s))
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Finds index of *value* in *array* of strings.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function findloc_string(array, value) result(ind)
    class(string), intent(in) :: array(:)
    character(*), intent(in) :: value
    integer :: ind
    integer :: i

    ind = 0
    do i = 1, size(array)
      if (array(i) == value) then
        ind = i
        exit
      end if
    end do
  end function

end module
