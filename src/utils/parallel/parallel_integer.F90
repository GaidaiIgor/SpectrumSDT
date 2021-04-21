module parallel_integer_mod
#include "type_list.macro"
#define TYPE_ID INTEGER_ID
#include "type_attributes.macro"
#define MPI_TYPE_NAME MPI_INTEGER
#include "parallel_template.F90"
end module
