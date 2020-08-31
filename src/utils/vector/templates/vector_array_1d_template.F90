#include "vector_template.F90"

!---------------------------------------------------------------------------------------------------------------------------------------------
! Converts ragged array into plain 2D matrix, assuming all 1D slices have the same size
!---------------------------------------------------------------------------------------------------------------------------------------------
function CONCAT3(vector_,TEMPLATE_TYPE_NAME,_to_2D_array)(vector) result(matrix)
  type(CONCAT2(vector_,TEMPLATE_TYPE_NAME)), intent(in) :: vector
  TEMPLATE_INNER_TYPE, allocatable :: matrix(:, :)
  TEMPLATE_INNER_TYPE, allocatable :: row(:)
  TEMPLATE_TYPE_OUT :: vector_elem
  integer :: n_rows, n_cols, i
  
  n_rows = vector % get_size()
  vector_elem = vector % get(1)
  n_cols = size(vector_elem % p)
  
  allocate(matrix(n_rows, n_cols))
  do i = 1, n_rows
    vector_elem = vector % get(i)
    matrix(i, :) = vector_elem % p
  end do
end function
