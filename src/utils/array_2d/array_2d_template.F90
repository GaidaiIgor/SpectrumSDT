#include "funcs.macro"  
implicit none

  interface array_2d
    module procedure :: CONCAT2(new_array_2d_,TEMPLATE_TYPE_NAME)
  end interface

  type CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)
    TEMPLATE_TYPE_OUT, allocatable :: p(:, :)
  end type

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Appends all 2D slices to form a single 2D array. Assumes the number of rows is the same in each 2D slice.
! *blocks* is an additional optional output argument, which saves initial slice number of each column in the resulting array.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(new_array_2d_,TEMPLATE_TYPE_NAME)(p) result(new_instance)
    TEMPLATE_TYPE, intent(in) :: p(:, :)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)) :: new_instance
    new_instance % p = p
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Appends all 2D slices to form a single 2D array. Assumes the number of rows is the same in each 2D slice.
! *blocks* is an additional optional output argument, which saves initial slice number of each column in the resulting array.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(flatten_array_2d_,TEMPLATE_TYPE_NAME)(array_2d, blocks) result(res)
    class(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), intent(in) :: array_2d(:)
    TEMPLATE_TYPE_OUT, allocatable :: res(:, :)
    integer, allocatable, intent(out), optional :: blocks(:)
    integer :: i, j, k, cols_total
    integer, allocatable :: blocks_temp(:)

    cols_total = CONCAT2(columns_number_,TEMPLATE_TYPE_NAME)(array_2d)
    allocate(res(size(array_2d(1) % p, 1), cols_total))
    allocate(blocks_temp(cols_total))

    k = 1
    do i = 1, size(array_2d)
      do j = 1, size(array_2d(i) % p, 2)
        res(:, k) = array_2d(i) % p(:, j)
        blocks_temp(k) = i
        k = k + 1
      end do
    end do

    if (present(blocks)) then
      blocks = blocks_temp
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes total number of columns in all slices of ragged array.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(columns_number_,TEMPLATE_TYPE_NAME)(array_2d) result(res)
    class(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), intent(in) :: array_2d(:)
    integer :: res
    integer :: i

    res = 0
    do i = 1, size(array_2d)
      res = res + size(array_2d(i) % p, 2)
    end do
  end function
