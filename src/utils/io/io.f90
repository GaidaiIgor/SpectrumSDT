module io_utils
  use io_base_mod
  use io_integer_mod
  use io_real_mod
  use io_complex_mod
  implicit none

  interface write_array
    module procedure :: write_array_integer, write_array_real, write_array_complex
  end interface

  interface write_matrix
    module procedure :: write_matrix_integer, write_matrix_real, write_matrix_complex
  end interface
end module
