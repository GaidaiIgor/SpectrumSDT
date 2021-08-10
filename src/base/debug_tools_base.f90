module debug_tools_base_mod
  use iso_fortran_env, only: real64
  implicit none
  
  character(:), allocatable :: debug_mode
  integer, allocatable :: debug_ints(:)
  real(real64) :: debug_real
  integer :: signal = 0
end module
