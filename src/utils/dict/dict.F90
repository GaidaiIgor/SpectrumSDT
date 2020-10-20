module dict_utils
  use dictionary
  use dict_integer_mod
  use dict_integer_array_mod
  use dict_char_str_mod
  use string_mod
  use vector_string_mod
  implicit none
  
  interface item_or_default
    module procedure :: item_or_default_integer, item_or_default_integer_array, item_or_default_char_str
  end interface

  interface add_if_absent
    module procedure :: add_if_absent_integer, add_if_absent_integer_array, add_if_absent_char_str
  end interface
  
contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns all keys in a given dict.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_key_set(dict) result(set)
    class(dictionary_t), intent(in) :: dict
    type(string), allocatable :: set(:)
    integer :: i
    character(:), allocatable :: key
    type(vector_string) :: vec
    type(dictionary_t) :: pair
    
    vec = vector_string()
    pair = .first. dict
    key = trim(.key. pair)
    call vec % push(string(key))
    
    do i = 2, len(dict)
      pair = .next. pair
      key = trim(.key. pair)
      call vec % push(string(key))
    end do
    set = vec % to_array()
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns dict with entries from dict1 that are not present in dict2.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function set_difference(dict1, dict2) result(diff)
    class(dictionary_t), intent(in) :: dict1, dict2 ! intent(in)
    type(dictionary_t) :: diff
    integer :: i
    character(:), allocatable :: next_key
    type(string), allocatable :: keys(:)

    diff = dict1 ! Make a copy
    keys = get_key_set(diff)
    ! Remove keys present in dict2
    do i = 1, size(keys)
      next_key = keys(i) % to_char_str()
      if (next_key .in. dict2) then
        call nullify(diff, next_key)
      end if
    end do
  end function
  
end module
