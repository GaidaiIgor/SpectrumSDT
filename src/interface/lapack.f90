module lapack_interface_mod
  use general_utils_mod
  use iso_fortran_env, only: real64
  implicit none

  interface lapack_eigensolver
    module procedure :: lapack_eigensolver_real, lapack_eigensolver_complex
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Provides a more convenient interface to dsyev subroutine.
! matrix - real symmetric matrix, whose eigenpairs are to be found. Contains eigenvectors on exit.
! eivals - eigenvalues of matrix, in ascending order.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine lapack_eigensolver_real(matrix, eivals)
    real(real64), intent(inout) :: matrix(:, :)
    real(real64), allocatable, intent(out) :: eivals(:)
    integer :: work_size, info
    real(real64) :: work_size_query(1)
    real(real64), allocatable :: work(:)

    call assert(size(matrix, 1) == size(matrix, 2), 'Error: matrix is not square in lapack_eigensolver')
    allocate(eivals(size(matrix, 1)))

    ! Trivial case: return empty vector if empty matrix is given
    if (size(eivals) == 0) then
      return
    end if

    ! Workspace size query
    work_size = -1
    call dsyev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work_size_query, work_size, info)
    call assert(info == 0, 'Error: lapack workspace query has failed, info =' // num2str(info))

    ! Actual solution
    work_size = int(work_size_query(1))
    allocate(work(work_size))
    call dsyev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work, work_size, info)
    call assert(info == 0, 'Error: lapack eigensolver has failed, info =' // num2str(info))
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Provides a more convenient interface to zheev subroutine.
! matrix - Hermitian matrix, whose eigenpairs are to be found. Contains eigenvectors on exit.
! eivals - eigenvalues of matrix, in ascending order.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine lapack_eigensolver_complex(matrix, eivals)
    complex(real64), intent(inout) :: matrix(:, :)
    real(real64), allocatable, intent(out) :: eivals(:)
    integer :: work_size, info
    real(real64), allocatable :: rwork(:)
    complex(real64) :: work_size_query(1)
    complex(real64), allocatable :: work(:)

    call assert(size(matrix, 1) == size(matrix, 2), 'Error: matrix is not square in lapack_eigensolver')
    allocate(eivals(size(matrix, 1)))
    allocate(rwork(3*size(matrix, 1) - 2))

    ! Trivial case: return empty vector if empty matrix is given
    if (size(eivals) == 0) then
      return
    end if

    ! Workspace size query
    work_size = -1
    call zheev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work_size_query, work_size, rwork, info)
    call assert(info == 0, 'Error: lapack workspace query has failed, info =' // num2str(info))

    ! Actual solution
    work_size = int(work_size_query(1))
    allocate(work(work_size))
    call zheev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work, work_size, rwork, info)
    call assert(info == 0, 'Error: lapack eigensolver has failed, info =' // num2str(info))
  end subroutine

end module
