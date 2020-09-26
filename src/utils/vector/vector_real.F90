module vector_real_mod
  use iso_fortran_env, only: real64
#include "type_list.macro"
#define TYPE_ID REAL_ID
#include "type_attributes.macro"
#include "vector_template.F90"
end module
