module dict_utils
  use dictionary
  use dict_integer_mod
  use dict_integer_array_mod
  use dict_char_str_mod
  use string_mod
  use vector_string_mod
  
  interface item_or_default
    module procedure :: item_or_default_integer, item_or_default_integer_array, item_or_default_char_str
  end interface

  interface add_if_absent
    module procedure :: add_if_absent_integer, add_if_absent_integer_array, add_if_absent_char_str
  end interface
  
contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns all keys in a given dict
!-------------------------------------------------------------------------------------------------------------------------------------------
  function key_set(dic) result(set)
    class(dict), intent(in) :: dic
    type(string), allocatable :: set(:)
    integer :: i
    character(:), allocatable :: key
    type(vector_string) :: vec
    type(dict) :: pair
    
    vec = vector_string()
    pair = .first. dic
    key = trim(.key. pair)
    call vec % push(string(key))
    
    do i = 2, len(dic)
      pair = .next. pair
      key = trim(.key. pair)
      call vec % push(string(key))
    end do
    set = vec % to_array()
  end function
  
end module
