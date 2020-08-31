module vector_string_mod
use string_mod
#include "type_list.macro"
#include "reset_definitions.macro"
#define TYPE_ID STRING_ID
#include "type_attributes.macro"
#include "vector_template.F90"
end module
