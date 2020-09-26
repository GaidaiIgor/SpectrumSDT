module array_1d_complex_mod
  use iso_fortran_env, only: real64
  use vector_complex_mod
#include "type_list.macro"
#define TYPE_ID COMPLEX_ID
#include "type_attributes.macro"
#include "array_1d_template.F90"
end module
