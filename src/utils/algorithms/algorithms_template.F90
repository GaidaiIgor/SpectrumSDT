#include "funcs.macro"
  use algorithms_base_mod
  implicit none

  interface prefix_sum_exclusive
    module procedure :: CONCAT2(prefix_sum_exclusive_,TEMPLATE_TYPE_NAME)
  end interface

  interface bubble_sort
    module procedure :: CONCAT2(bubble_sort_,TEMPLATE_TYPE_NAME)
  end interface

  interface get_unique
    module procedure :: CONCAT2(get_unique_,TEMPLATE_TYPE_NAME)
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes exclusive prefix sum, i.e. res(i) = sum( [array(1)..array(i - 1)] ).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(prefix_sum_exclusive_,TEMPLATE_TYPE_NAME)(array) result(res)
    TEMPLATE_TYPE, intent(in) :: array(:)
    TEMPLATE_TYPE_OUT, allocatable :: res(:)
    integer :: i
    
    allocate(res(size(array)))
    res(1) = 0
    do i = 2, size(array)
      res(i) = res(i - 1) + array(i - 1)
    end do
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Sorts an array using bubble sort.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(bubble_sort_,TEMPLATE_TYPE_NAME)(array, index_permutation) result(sorted)
    TEMPLATE_TYPE, intent(in) :: array(:)
    TEMPLATE_TYPE_OUT, allocatable :: sorted(:)
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
! Sorts and filters out non-unique values.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(get_unique_,TEMPLATE_TYPE_NAME)(array) result(unique)
    TEMPLATE_TYPE, intent(in) :: array(:)
    TEMPLATE_TYPE_OUT, allocatable :: unique(:)
    integer :: i
    TEMPLATE_TYPE_OUT, allocatable :: sorted(:)
    type(CONCAT2(vector_,TEMPLATE_TYPE_NAME)) :: unique_vec

    ! Filter out non-unique values
    sorted = CONCAT2(bubble_sort_,TEMPLATE_TYPE_NAME)(array)
    unique_vec = CONCAT2(vector_,TEMPLATE_TYPE_NAME)()
    call unique_vec % push(sorted(1))
    do i = 2, size(sorted)
#if TYPE_ID == REAL_ID
      if (unique_vec % top() .aeq. sorted(i)) then
#else
      if (unique_vec % top() == sorted(i)) then
#endif
        cycle
      end if
      call unique_vec % push(sorted(i))
    end do
    unique = unique_vec % to_array()
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Evaluates expression of the form v^T * A * v for a vector *v* and a symmetric matrix *A*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(calculate_vT_A_v_symmetric_,TEMPLATE_TYPE_NAME)(A, v) result(res)
    TEMPLATE_TYPE, intent(in) :: A(:, :)
    complex(real64), intent(in) :: v(:)
    real(real64) :: res
    integer :: i, j
    real(real64) :: diag_sum, offdiag_sum

    diag_sum = 0
    do j = 1, size(A, 2)
      diag_sum = diag_sum + real(conjg(v(j)) * A(j, j) * v(j))
    end do

    offdiag_sum = 0
    do j = 2, size(A, 2)
      do i = 1, j - 1
        offdiag_sum = offdiag_sum + real(conjg(v(i)) * A(i, j) * v(j))
      end do
    end do

    res = diag_sum + 2 * offdiag_sum
  end function

