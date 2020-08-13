module lapack_interface_mod
  use general_utils
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Provides a more convenient interface to dsyev subroutine
! matrix - matrix, whose eigenpairs are to be found. Contains eigenvectors on exit.
! eivals - eigenvalues of matrix
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine lapack_eigensolver(matrix, eivals)
    real*8, intent(inout) :: matrix(:, :)
    real*8, intent(out) :: eivals(:)
    integer :: work_size, info
    real*8 :: work_size_query
    real*8, allocatable :: work(:)

    call assert(size(matrix, 1) == size(matrix, 2), 'Error: matrix is not square')
    call assert(size(matrix, 1) == size(eivals), 'Error: eivals array has to match matrix size')

    ! Workspace size query
    work_size = -1
    call dsyev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work_size_query, work_size, info)
    call assert(info == 0, 'Error: lapack work space query has failed, info =' // num2str(info))

    ! Actual solution
    work_size = int(work_size_query)
    allocate(work(work_size))
    call dsyev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work, work_size, info)
    call assert(info == 0, 'Error: lapack eigensolver has failed, info =' // num2str(info))
  end subroutine
end module
