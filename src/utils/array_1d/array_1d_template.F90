#include "funcs.macro"
use vector_integer_mod
implicit none

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
! Overrides default constructor to copy by value.
!---------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(new_array_1d_,TEMPLATE_TYPE_NAME)(p) result(new_instance)
    TEMPLATE_TYPE, intent(in) :: p(:)
    type(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)) :: new_instance
    new_instance % p = p
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Writes array_1d.
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(write_array_1d_,TEMPLATE_TYPE_NAME)(this, unit, iotype, v_list, iostat, iomsg)
    class(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)), intent(in) :: this
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    write(unit, *, iostat = iostat) this % p
  end

!---------------------------------------------------------------------------------------------------------------------------------------------
! Returns inner array.
!---------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT3(array_1d_,TEMPLATE_TYPE_NAME,_to_plain_array)(this) result(res)
    class(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)), intent(in) :: this
    TEMPLATE_TYPE_OUT, allocatable :: res(:)
    res = this % p
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Appends all elements of 1D array of 1D arrays (ragged 2D array) to form a single plain 1D array.
! *index_map* contains positions of first elements from each ragged dimension (column) of the original array in the flattened array.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(flatten_1d_array_1d_,TEMPLATE_TYPE_NAME)(ragged_array, flattened, index_map)
    class(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)), intent(in) :: ragged_array(:)
    TEMPLATE_TYPE_OUT, allocatable, intent(out) :: flattened(:)
    integer, allocatable, optional, intent(out) :: index_map(:)
    integer :: i
    integer :: index_map_act(size(ragged_array))
    type(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: accumulator

    accumulator = CONCAT2(vector_,TEMPLATE_TYPE_NAME)()
    do i = 1, size(ragged_array)
      index_map_act(i) = accumulator % get_size() + 1
      call accumulator % push_all(ragged_array(i) % p)
    end do

    flattened = accumulator % to_array()
    if (present(index_map)) then
      index_map = index_map_act
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Appends all elements of 2D array of 1D arrays (ragged 3D array) to form a single plain 1D array in column-wise order.
! *index_map* contains rows and columns of first elements from each ragged dimension of the original array in the flattened array.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(flatten_2d_array_1d_,TEMPLATE_TYPE_NAME)(ragged_array, flattened, index_map)
    class(CONCAT2(array_1d_,TEMPLATE_TYPE_NAME)), intent(in) :: ragged_array(:, :)
    TEMPLATE_TYPE_OUT, allocatable, intent(out) :: flattened(:)
    integer, allocatable, optional, intent(out) :: index_map(:, :)
    integer :: i, j
    integer :: index_map_act(size(ragged_array, 1), size(ragged_array, 2))
    type(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: accumulator

    accumulator = CONCAT2(vector_,TEMPLATE_TYPE_NAME)()
    do j = 1, size(ragged_array, 2)
      do i = 1, size(ragged_array, 1)
        index_map_act(i, j) = accumulator % get_size() + 1
        call accumulator % push_all(ragged_array(i, j) % p)
      end do
    end do

    flattened = accumulator % to_array()
    if (present(index_map)) then
      index_map = index_map_act
    end if
  end subroutine
