module io_utils_complex_mod
  use io_utils_real_mod
  use iso_fortran_env, only: real64
#include "type_list.macro"
#define TYPE_ID COMPLEX_ID
#define TEMPLATE_ELEM_SPEC '"(",G0.15,",",G0.15,")"'
#include "type_attributes.macro"
#include "io_utils_template.F90"

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints real and imaginary parts of a given complex *matrix* to 2 files, formed from *file_name*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_matrix_complex_2(matrix, file_name, append, print_size)
    complex(real64), intent(in) :: matrix(:, :)
    character(*), intent(in) :: file_name
    integer, optional, intent(in) :: append, print_size

    call write_matrix(real(matrix), file_name // '.real', append, print_size)
    call write_matrix(aimag(matrix), file_name // '.imag', append, print_size)
  end subroutine

end module
