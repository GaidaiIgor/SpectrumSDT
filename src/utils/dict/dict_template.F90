#include "funcs.macro" 

use dictionary
implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns item by given key or default value if the key is not present in the dictionary
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(item_or_default_,TEMPLATE_TYPE_NAME)(dic, key, default) result(res)
    type(dictionary_t), intent(inout) :: dic ! inout because assign does not declare it as in
    character(*), intent(in) :: key
    TEMPLATE_TYPE, intent(in) :: default
    TEMPLATE_TYPE_OUT :: res
    
    if (key .in. dic) then
#if TYPE_ID == CHAR_STR_ID
      res = extract_string(dic, key)
#else
      call assign(res, dic, key)
#endif
    else
      res = default
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Adds an item if it's not already present
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(add_if_absent_,TEMPLATE_TYPE_NAME)(dic, key, value)
    type(dictionary_t), intent(inout) :: dic
    character(*), intent(in) :: key
    TEMPLATE_TYPE, intent(in) :: value
    
    if (key .in. dic) then
      return
    end if
#if TYPE_ID == CHAR_STR_ID
    call put_string(dic, key, value)
#else
    call extend(dic, (key .kv. value))
#endif
  end subroutine
