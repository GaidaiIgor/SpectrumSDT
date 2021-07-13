#include "funcs.macro"  
implicit none

  type CONCAT2(array_3d_,TEMPLATE_TYPE_NAME)
    TEMPLATE_TYPE_OUT, allocatable :: p(:, :, :)
  end type
