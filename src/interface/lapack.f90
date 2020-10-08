module lapack_interface_mod
  use general_utils
  use iso_fortran_env, only: real64
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Provides a more convenient interface to dsyev subroutine
! matrix - matrix, whose eigenpairs are to be found. Contains eigenvectors on exit.
! eivals - eigenvalues of matrix
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine lapack_eigensolver(matrix, eivals)
    real(real64), intent(inout) :: matrix(:, :)
    real(real64), allocatable, intent(out) :: eivals(:)
    integer :: work_size, info
    real(real64) :: work_size_query(1)
    real(real64), allocatable :: work(:)

    call assert(size(matrix, 1) == size(matrix, 2), 'Error: matrix is not square in lapack_eigensolver')
    allocate(eivals(size(matrix, 1)))

    ! Workspace size query
    work_size = -1
    call dsyev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work_size_query, work_size, info)
    call assert(info == 0, 'Error: lapack work space query has failed, info =' // num2str(info))

    ! Actual solution
    work_size = int(work_size_query(1))
    allocate(work(work_size))
    call dsyev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work, work_size, info)
    call assert(info == 0, 'Error: lapack eigensolver has failed, info =' // num2str(info))
  end subroutine
end module
