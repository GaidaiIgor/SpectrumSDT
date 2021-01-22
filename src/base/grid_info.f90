module grid_info_mod
  use iso_fortran_env, only: real64
  implicit none

  type :: grid_info
    real(real64) :: from
    real(real64) :: to
    real(real64) :: step
    real(real64), allocatable :: points(:)
    real(real64), allocatable :: jac(:)
  end type
end module
