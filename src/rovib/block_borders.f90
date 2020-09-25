module block_borders_mod
  implicit none

  type :: block_borders
    integer :: left
    integer :: right
    integer :: top
    integer :: bottom
  end type
end module
