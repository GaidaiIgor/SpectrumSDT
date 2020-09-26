module array_2d_complex_mod
  use iso_fortran_env, only: real64
#include "type_list.macro"
#define TYPE_ID COMPLEX_ID
#include "type_attributes.macro"
#include "array_2d_template.F90"
end module
