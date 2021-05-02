module general_utils_real_array_mod
  use iso_fortran_env, only: real64
#include "type_list.macro"
#define TYPE_ID REAL_ARRAY_ID
#include "type_attributes.macro"
#include "general_utils_template.F90"
end module
