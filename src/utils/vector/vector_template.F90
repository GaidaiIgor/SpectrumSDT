#include "funcs.macro"
use general_utils
implicit none

interface CONCAT2(vector_,TEMPLATE_TYPE_NAME)
  module procedure :: CONCAT2(new_vector_,TEMPLATE_TYPE_NAME)
end interface

type CONCAT2(vector_,TEMPLATE_TYPE_NAME)
  integer :: size
  integer :: resize_factor
  TEMPLATE_TYPE_OUT, allocatable :: storage(:)

contains
  procedure :: push => CONCAT2(push_,TEMPLATE_TYPE_NAME)
  procedure :: push_all => CONCAT2(push_all_,TEMPLATE_TYPE_NAME)
  procedure :: append_vector => CONCAT2(append_vector_,TEMPLATE_TYPE_NAME)
  procedure :: clear => CONCAT2(clear_,TEMPLATE_TYPE_NAME)
  procedure :: get => CONCAT2(get_,TEMPLATE_TYPE_NAME)
  procedure :: top => CONCAT2(top_,TEMPLATE_TYPE_NAME)
  procedure :: set => CONCAT2(set_,TEMPLATE_TYPE_NAME)
  procedure :: get_size => CONCAT2(get_size_,TEMPLATE_TYPE_NAME)
  procedure :: reverse => CONCAT2(reverse_,TEMPLATE_TYPE_NAME)
  procedure :: to_array => CONCAT2(to_array_,TEMPLATE_TYPE_NAME)
  procedure :: to_existing_array => CONCAT2(to_existing_array_,TEMPLATE_TYPE_NAME)
  
  procedure :: write => CONCAT2(write_vector_,TEMPLATE_TYPE_NAME)
  generic :: write(formatted) => write
end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes a new instance.
!-------------------------------------------------------------------------------------------------------------------------------------------
elemental function CONCAT2(new_vector_,TEMPLATE_TYPE_NAME)(initial_capacity, resize_factor) result(new_instance)
  integer, optional, intent(in) :: initial_capacity, resize_factor
  integer :: initial_capacity_act, resize_factor_act
  type(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: new_instance
  
  initial_capacity_act = arg_or_default(initial_capacity, 256)
  resize_factor_act = arg_or_default(resize_factor, 2)
  new_instance % size = 0
  new_instance % resize_factor = resize_factor_act
  allocate(new_instance % storage(initial_capacity_act))
end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Pushes a value to vector.
!-------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(push_,TEMPLATE_TYPE_NAME)(this, value)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this
  TEMPLATE_TYPE, intent(in) :: value
  TEMPLATE_TYPE_OUT, allocatable :: temp(:)

  if (this % size == size(this % storage)) then
    allocate(temp(this % size * this % resize_factor)) ! make a new allocation
    temp(1:this % size) = this % storage ! move items to new allocation
    call move_alloc(temp, this % storage) ! swap pointers, deallocate old pointer
  end if
  
  this % size = this % size + 1
  this % storage(this % size) = value
end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calls push for all elements of *array*.
!-------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(push_all_,TEMPLATE_TYPE_NAME)(this, array)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this
  TEMPLATE_TYPE, intent(in) :: array(:)
  integer :: i

  do i = 1,size(array)
    call this % push(array(i))
  end do
end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Appends another vector.
!-------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(append_vector_,TEMPLATE_TYPE_NAME)(this, other)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this, other
  integer :: i

  do i = 1,other % get_size()
    call this % push(other % get(i))
  end do
end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Clears vector.
!-------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(clear_,TEMPLATE_TYPE_NAME)(this)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)), intent(inout) :: this
  this % size = 0
end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns a value at a given *index*.
!-------------------------------------------------------------------------------------------------------------------------------------------
function CONCAT2(get_,TEMPLATE_TYPE_NAME)(this, index) result(res)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this
  integer, intent(in) :: index
  TEMPLATE_TYPE_OUT :: res
  res = this % storage(index)
end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns the last added value.
!-------------------------------------------------------------------------------------------------------------------------------------------
function CONCAT2(top_,TEMPLATE_TYPE_NAME)(this) result(res)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this
  TEMPLATE_TYPE_OUT :: res
  res = this % storage(this % get_size())
end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets a value at a given *index* to *value*.
!-------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(set_,TEMPLATE_TYPE_NAME)(this, index, value)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this
  integer, intent(in) :: index
  TEMPLATE_TYPE, intent(in) :: value
  this % storage(index) = value
end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns number of elements in this vector.
!-------------------------------------------------------------------------------------------------------------------------------------------
function CONCAT2(get_size_,TEMPLATE_TYPE_NAME)(this) result(res)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this
  integer :: res
  res = this % size
end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Reverses the order of content of this vector.
!-------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(reverse_,TEMPLATE_TYPE_NAME)(this)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this
  this % storage(1:this % size) = this % storage(this % size:1:-1)
end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts vector into plain array.
!-------------------------------------------------------------------------------------------------------------------------------------------
function CONCAT2(to_array_,TEMPLATE_TYPE_NAME)(this) result(array)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this
  TEMPLATE_TYPE_OUT, allocatable :: array(:)

  allocate(array(this % get_size()))
  array = this % storage(1:this % get_size())
end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Moves vector content into an existing *array*.
!-------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(to_existing_array_,TEMPLATE_TYPE_NAME)(this, array)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: this
  TEMPLATE_TYPE_OUT, allocatable, intent(inout) :: array(:)

  if (size(array) < this % get_size()) then
    stop 'Target array is too small'
  end if
  array = this % storage(1:this % get_size())
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------------------
! Vector writing procedure.
!---------------------------------------------------------------------------------------------------------------------------------------------
subroutine CONCAT2(write_vector_,TEMPLATE_TYPE_NAME)(this, unit, iotype, v_list, iostat, iomsg)
  class(CONCAT2(vector_,TEMPLATE_TYPE_NAME)), intent(in) :: this
  integer, intent(in) :: unit
  character(*), intent(in) :: iotype
  integer, intent(in) :: v_list(:)
  integer, intent(out) :: iostat
  character(*), intent(inout) :: iomsg
  integer :: i
  
  do i = 1,this % get_size() - 1
    write(unit, *, iostat = iostat) this % get(i), ','
  end do
  write(unit, *, iostat = iostat) this % get(this % get_size())
end subroutine
