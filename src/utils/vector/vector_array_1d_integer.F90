module vector_array_1d_integer_mod
use array_1d_integer_mod
#include "type_list.macro"
#include "reset_definitions.macro"
#define TYPE_ID ARRAY_1D_INTEGER_ID
#include "type_attributes.macro"
#include "templates/vector_array_1d_template.F90"
end module