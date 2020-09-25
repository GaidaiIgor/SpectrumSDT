module string_utils
  use vector_string_mod
  implicit none

  interface strsplit
    module procedure :: strsplit_char_str, strsplit_string
  end interface

contains
!---------------------------------------------------------------------------------------------------------------------------------------------
! Splits a string into tokens using specified delimiter
!---------------------------------------------------------------------------------------------------------------------------------------------
  function strsplit_char_str(str, delim) result(tokens_arr)
    character(*), intent(in) :: str
    character(*), optional, intent(in) :: delim
    type(string), allocatable :: tokens_arr(:)

    tokens_arr = strsplit_string(string(str), delim)
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Splits a string into tokens using specified delimiter
!---------------------------------------------------------------------------------------------------------------------------------------------
  function strsplit_string(str, delim) result(tokens_arr)
    class(string), intent(in) :: str
    character(*), optional, intent(in) :: delim
    type(string), allocatable :: tokens_arr(:)
    ! local
    type(vector_string) :: tokens
    character(:), allocatable :: delim_act
    integer :: delim_pos, search_start
    type(string) :: substring
    
    delim_act = arg_or_default(delim, ' ')
    tokens = vector_string()
    search_start = 1
    do
      substring = str%substr(search_start)
      delim_pos = index(substring%to_char_str(), delim_act)
      if (delim_pos == 0) then ! next delim is not found
        call tokens%push(str%substr(search_start))
        exit
      end if
      delim_pos = delim_pos + search_start - 1 ! position of delim relative to the beginning of string
      call tokens%push(str%substr(search_start, delim_pos - 1))      
      search_start = delim_pos + len(delim_act)
    end do
    tokens_arr = tokens%to_array()
  end function
end module
