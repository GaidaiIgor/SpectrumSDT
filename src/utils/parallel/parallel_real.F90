module parallel_real_mod
#include "type_list.macro"
#define TYPE_ID REAL_ID
#include "type_attributes.macro"
#define MPI_TYPE_NAME MPI_DOUBLE_PRECISION
#include "parallel_template.F90"
end module
