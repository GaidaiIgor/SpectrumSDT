module array_1d_real_mod
  use iso_fortran_env, only: real64
  use vector_real_mod
#include "type_list.macro"
#define TYPE_ID REAL_ID
#include "type_attributes.macro"
#include "array_1d_template.F90"
end module
