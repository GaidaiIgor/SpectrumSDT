module vector_array_1d_real_mod
use array_1d_real_mod
#include "type_list.macro"
#include "reset_definitions.macro"
#define TYPE_ID ARRAY_1D_REAL_ID
#include "type_attributes.macro"
#include "templates/vector_array_1d_template.F90"
end module