module io_complex_mod
#include "type_list.macro"
#define TYPE_ID COMPLEX_ID
#define TEMPLATE_ELEM_SPEC '"(",G0.15,",",G0.15,")"'
#include "type_attributes.macro"
#include "io_template.F90"
end module
