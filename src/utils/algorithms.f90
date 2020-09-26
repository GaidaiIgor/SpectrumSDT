module algorithms
  use general_utils
  use iso_fortran_env, only: real64
  implicit none

  interface prefix_sum
    module procedure :: prefix_sum_integer
  end interface

  interface prefix_sum_exclusive
    module procedure :: prefix_sum_exclusive_integer
  end interface
  
  interface bubble_sort
    module procedure :: bubble_sort_real
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes prefix sum, i.e. res(i) = sum( [array(1)..array(i)] )
!-------------------------------------------------------------------------------------------------------------------------------------------
  function prefix_sum_integer(array) result(res)
    integer, intent(in) :: array(:)
    integer, allocatable :: res(:)
    integer :: i
    
    allocate(res(size(array)))
    res(1) = array(1)
    do i = 2, size(array)
      res(i) = res(i - 1) + array(i)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes exclusive prefix sum, i.e. res(i) = sum( [array(1)..array(i - 1)] )
!-------------------------------------------------------------------------------------------------------------------------------------------
  function prefix_sum_exclusive_integer(array) result(res)
    integer, intent(in) :: array(:)
    integer, allocatable :: res(:)
    integer :: i
    
    allocate(res(size(array)))
    res(1) = 0
    do i = 2, size(array)
      res(i) = res(i - 1) + array(i - 1)
    end do
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Sorts an array using bubble sort
!-------------------------------------------------------------------------------------------------------------------------------------------
  function bubble_sort_real(array, index_permutation) result(sorted)
    real(real64), intent(in) :: array(:)
    real(real64), allocatable :: sorted(:)
    integer, allocatable, optional, intent(out) :: index_permutation(:)
    integer :: i, j
    integer, allocatable :: index_permutation_act(:)
    
    sorted = array
    index_permutation_act = [(i, i = 1, size(array))]
    ! after every major iteration, the biggest remaining element ascends to the top, so we lower the upper limit
    do i = size(sorted) - 1, 1, -1
      do j = 1, i
        if (sorted(j) > sorted(j + 1)) then
          call swap(sorted(j), sorted(j + 1))
          call swap(index_permutation_act(j), index_permutation_act(j + 1))
        end if
      end do
    end do
    
    if (present(index_permutation)) then
      index_permutation = index_permutation_act
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Estimates the value of *query_y* at *query_x* by building line equation from 2 points in *known_xs* and *known_ys*
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine linear_estimation(known_xs, known_ys, query_x, query_y)
    real(real64), intent(in) :: known_xs(2), known_ys(2)
    real(real64), intent(in) :: query_x
    real(real64), intent(out) :: query_y
    query_y = (known_ys(2) - known_ys(1)) / (known_xs(2) - known_xs(1)) * (query_x - known_xs(1)) + known_ys(1)
  end subroutine

end module
