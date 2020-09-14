module io_real_mod
#include "type_list.macro"
#include "reset_definitions.macro"
#define TYPE_ID REAL_ID
#define TEMPLATE_ELEM_SPEC 'G23.15'
#include "type_attributes.macro"
#include "io_template.F90"
end module
