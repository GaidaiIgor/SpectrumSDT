module general_string_mod
  use string_mod
#include "type_list.macro"
#define TYPE_ID STRING_ID
#include "type_attributes.macro"
#include "general_utils_template.F90"
end module
