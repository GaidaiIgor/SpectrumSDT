#include "funcs.macro" 
  use dict_utils_base_mod
  implicit none

  interface item_or_default
    module procedure :: CONCAT2(item_or_default_,TEMPLATE_TYPE_NAME)
  end interface

  interface add_if_absent
    module procedure :: CONCAT2(add_if_absent_,TEMPLATE_TYPE_NAME)
  end interface
  

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns item by given key or default value if the key is not present in the dictionary
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(item_or_default_,TEMPLATE_TYPE_NAME)(dict, key, default) result(res)
    type(dictionary_t) :: dict ! intent(in)
    character(*), intent(in) :: key
    TEMPLATE_TYPE, intent(in) :: default
    TEMPLATE_TYPE_OUT :: res
    
    if (key .in. dict) then
#if TYPE_ID == CHAR_STR_ID
      res = extract_string(dict, key)
#else
      call assign(res, dict, key)
#endif
    else
      res = default
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Adds an item if it's not already present
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(add_if_absent_,TEMPLATE_TYPE_NAME)(dict, key, value)
    type(dictionary_t), intent(inout) :: dict
    character(*), intent(in) :: key
    TEMPLATE_TYPE, intent(in) :: value
    
    if (key .in. dict) then
      return
    end if
#if TYPE_ID == CHAR_STR_ID
    call put_string(dict, key, value)
#else
    call extend(dict, (key .kv. value))
#endif
  end subroutine
