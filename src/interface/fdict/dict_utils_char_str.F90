module dict_utils_char_str_mod
  use general_utils_mod

#include "type_list.macro"
#define TYPE_ID CHAR_STR_ID
#include "type_attributes.macro"
#include "dict_utils_template.F90"

!-------------------------------------------------------------------------------------------------------------------------------------------
! Forms key name for string size based on a given key.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_key_size_name(key) result(key_size)
    character(*), intent(in) :: key
    character(:), allocatable :: key_size
    key_size = key // '_size'
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Puts a string to dictionary. Updates the value if the key is already present.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine put_string(dict, key, str)
    type(dictionary_t), intent(inout) :: dict
    character(*), intent(in) :: key
    character(*), optional, intent(in) :: str
    character(:), allocatable :: key_size, str_act
    type(dictionary_t) :: pair
    integer :: str_len

    str_act = arg_or_default(str, '')
    ! Remove an existing value. fdict automatically replaces plain values, but apparently that does not extend to 'pointer' objects, so there will be a runtime error 
    ! if one tries to re-insert the same key without calling nullify first.
    if (key .in. dict) then
      call nullify(dict, key)
    end if
    str_len = len(str_act)
    key_size = get_key_size_name(key)
    pair = (key_size .kv. str_len) // (key .kv. str_act)
    ! Copying by value through .kv. does not have a defined interface, but .kvp. does.
    ! .kvp. is supposed to copy a reference, but fdict seems to be able to handle copying a reference to a local variable.
    ! I tried to save a pointer to `pair` and export it. Referencing the pointer outside of this subroutine returned invalid values, so the memory that belonged to `pair` was properly reset.
    ! And yet requesting the corresponding key from the global dict worked, which means it must've saved another copy internally, so it should be ok to reference local variables in this way...
    dict = dict // (key .kvp. pair)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Extracts a string from dictionary.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function extract_string(dict, key) result(str)
    type(dictionary_t) :: dict ! intent(in)
    character(*), intent(in) :: key
    character(:), allocatable :: str
    integer :: str_len
    character(:), allocatable :: key_size
    type(dictionary_t) :: pair

    call associate(pair, dict, key)
    key_size = get_key_size_name(key)
    call assign(str_len, pair, key_size)
    allocate(character(str_len) :: str)
    call assign(str, pair, key)
  end function

end module
