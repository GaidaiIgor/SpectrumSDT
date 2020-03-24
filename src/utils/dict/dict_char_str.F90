module dict_char_str_mod

#include "type_list.macro"
#include "reset_definitions.macro"
#define TYPE_ID CHAR_STR_ID
#include "type_attributes.macro"
#include "dict_template.F90"

!-------------------------------------------------------------------------------------------------------------------------------------------
! Forms key for size based on a given key
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_key_size_name(key) result(key_size)
    character(*), intent(in) :: key
    character(:), allocatable :: key_size
    key_size = key // '_size'
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Puts a string to dictionary
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine put_string(dic, key, str)
    type(dict), intent(inout) :: dic
    character(*), intent(in) :: key, str
    character(:), allocatable :: key_size
    type(dict) :: pair
    integer :: str_len

    ! Remove an existing value
    if (key .in. dic) then
      call nullify(dic, key)
    end if
    str_len = len(str)
    key_size = get_key_size_name(key)
    pair = (key_size .kv. str_len) // (key .kv. str)
    dic = dic // (key .kvp. pair)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Extracts a string from dictionary
!-------------------------------------------------------------------------------------------------------------------------------------------
  function extract_string(dic, key) result(str)
    type(dict), intent(inout) :: dic
    character(*), intent(in) :: key
    character(:), allocatable :: str
    integer :: str_len
    character(:), allocatable :: key_size
    type(dict) :: pair
    
    call associate(pair, dic, key)
    key_size = get_key_size_name(key)
    call assign(str_len, pair, key_size)
    allocate(character(str_len) :: str)
    call assign(str, pair, key)
  end function
end module
