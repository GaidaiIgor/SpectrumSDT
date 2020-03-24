module string_mod
  use general_utils
  
  interface string
    module procedure :: to_string_char_str
  end interface
  
  interface assignment(=)
    procedure :: assign_char_str, assign_char_str_array
  end interface
  
  type :: string 
    character(:), allocatable :: s
  contains    
    procedure :: write => write_string
    generic :: write(formatted) => write

    procedure :: compare_string
    procedure :: compare_string_arr
    generic :: operator(==) => compare_string, compare_string_arr
    
    procedure :: to_char_str
    procedure :: to_integer
    procedure :: to_real
    procedure :: to_complex
    procedure :: length
    procedure :: substr
  end type
  
contains
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_string(this, unit, iotype, v_list, iostat, iomsg)
    class(string), intent(in) :: this
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    write(unit, *, iostat = iostat) this % s // ' '
  end subroutine
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Compares two strings
!---------------------------------------------------------------------------------------------------------------------------------------------
  function compare_string(this, other) result(res)
    class(string), intent(in) :: this, other
    logical :: res
    res = this % s == other % s
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Compares a string with an array of strings
!---------------------------------------------------------------------------------------------------------------------------------------------
  function compare_string_arr(this, others) result(res)
    class(string), intent(in) :: this
    class(string), intent(in) :: others(:)
    logical, allocatable :: res(:)
    integer :: i

    allocate(res(size(others)))
    res = [(this % s == others(i) % s, i = 1, size(others))]
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Converts string type to char string
!---------------------------------------------------------------------------------------------------------------------------------------------
  function to_char_str(this) result(char_str)
    class(string), intent(in) :: this
    character(:), allocatable :: char_str
    char_str = this % s
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to int
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function to_integer(this) result(res)
    class(string), intent(in) :: this
    integer :: res
    character(:), allocatable :: char_str
    
    char_str = this % s
    read(char_str, *) res
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to real
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function to_real(this) result(res)
    class(string), intent(in) :: this
    real*8 :: res
    character(:), allocatable :: char_str
    
    char_str = this % s
    read(char_str, *) res
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to real
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function to_complex(this) result(res)
    class(string), intent(in) :: this
    complex*16 :: res
    character(:), allocatable :: char_str
    
    char_str = this % s
    read(char_str, *) res
  end function
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Length of string
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental function length(this) result(res)
    class(string), intent(in) :: this
    integer :: res
    res = len(this%s)
  end function
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Substring
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental function substr(this, start, end) result(res)
    class(string), intent(in) :: this
    integer, intent(in), optional :: start, end
    type(string) :: res
    integer :: start_act, end_act
    
    start_act = arg_or_default(start, 1)
    end_act = arg_or_default(end, this%length())
    res = this%s(start_act:end_act)
  end function
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! MODULE FUNCTIONS
!---------------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------
! Constructs a string from char string
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental function to_string_char_str(char_str) result(str)
    character(*), intent(in) :: char_str
    type(string) :: str
    str%s = char_str
  end function
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Asssigns a char string
!---------------------------------------------------------------------------------------------------------------------------------------------
  elemental subroutine assign_char_str(this, char_str)
    type(string), intent(inout) :: this
    character(*), intent(in) :: char_str
    this = string(char_str)
  end subroutine
  
!---------------------------------------------------------------------------------------------------------------------------------------------
! Assigns a char string array
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assign_char_str_array(array_lhs, array_rhs)
    type(string), allocatable, intent(out) :: array_lhs(:)
    character(*), intent(in) :: array_rhs(:)
    type(string), allocatable :: string_array_rhs(:)
    
    string_array_rhs = string(array_rhs)
    array_lhs = string_array_rhs
  end subroutine

!---------------------------------------------------------------------------------------------------------------------------------------------
! Joins strings from a string array (arr) into a single (char_str). (join_str) specifies what string separates the strings in (arr)
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

end module