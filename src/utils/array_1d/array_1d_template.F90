#include "funcs.macro"
use vector_integer_mod

interface CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)
  module procedure :: CONCAT2(new_array_1d_,TEMPLATE_TYPE_NAME)
end interface

type CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)
  TEMPLATE_TYPE_OUT, allocatable :: p(:)
contains
  procedure :: write => CONCAT2(write_array_1d_,TEMPLATE_TYPE_NAME)
  generic :: write(formatted) => write
  
  procedure :: to_plain_array => CONCAT3(array_1d_,TEMPLATE_TYPE_NAME,_to_plain_array)
end type

contains

!---------------------------------------------------------------------------------------------------------------------------------------------
! Override default constructor to copy by value
!---------------------------------------------------------------------------------------------------------------------------------------------
function CONCAT2(new_array_1d_,TEMPLATE_TYPE_NAME)(p) result(new_instance)
  TEMPLATE_TYPE, intent(in) :: p(:)
  type(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)) :: new_instance
  new_instance % p = p
end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes
!---------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(write_array_1d_,TEMPLATE_TYPE_NAME)(this, unit, iotype, v_list, iostat, iomsg)
  class(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)), intent(in) :: this
  integer, intent(in) :: unit
  character(*), intent(in) :: iotype
  integer, intent(in) :: v_list(:)
  integer, intent(out) :: iostat
  character(*), intent(inout) :: iomsg
  integer :: i
  write(unit, *, iostat = iostat) this % p
end

!---------------------------------------------------------------------------------------------------------------------------------------------
! Returns inner array
!---------------------------------------------------------------------------------------------------------------------------------------------
function CONCAT3(array_1d_,TEMPLATE_TYPE_NAME,_to_plain_array)(this) result(res)
  class(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)), intent(in) :: this
  TEMPLATE_TYPE_OUT, allocatable :: res(:)
  res = this % p
end function

!-----------------------------------------------------------------------
! Appends all 1D slices to form a single 1D array
! blocks is an additional optional output argument, which saves initial slice number of each element in the resulting array
!-----------------------------------------------------------------------
function CONCAT2(flatten_1d_,TEMPLATE_TYPE_NAME)(ragged_array, blocks) result(res)
  class(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)), intent(in) :: ragged_array(:)
  TEMPLATE_TYPE_OUT, allocatable :: res(:)
  integer, allocatable, optional, intent(out) :: blocks(:)
  type(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: accumulator
  type(vector_integer) :: blocks_vec
  integer :: i, j
  accumulator = CONCAT2(vector_,TEMPLATE_TYPE_NAME)()
  blocks_vec = vector_integer()

  do i = 1,size(ragged_array)
    do j = 1,size(ragged_array(i) % p)
      call accumulator % push(ragged_array(i) % p(j))
      call blocks_vec % push(i)
    end do
  end do

  res = accumulator % to_array()
  if (present(blocks)) then
    blocks = blocks_vec % to_array()
  end if
end function
