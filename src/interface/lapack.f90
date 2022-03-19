module lapack_interface_mod
  use iso_fortran_env, only: real64
  implicit none

  interface lapack_eigensolver
    module procedure :: lapack_eigensolver_real, lapack_eigensolver_complex
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks info and stops the program in case of error.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_info(info)
    integer, intent(in) :: info
    if (info /= 0) then
      print *, 'Error: lapack has failed, info =', info
      stop
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Provides a more convenient interface to dsyev subroutine.
! matrix - real square symmetric matrix, whose eigenpairs are to be found. Contains eigenvectors on exit.
! eivals - eigenvalues of matrix, in ascending order. Size has be to equal to size(matrix, 1)
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine lapack_eigensolver_real_plain(matrix, eivals)
    real(real64), intent(inout) :: matrix(:, :)
    real(real64), intent(out) :: eivals(:)
    integer :: work_size, info
    real(real64) :: work_size_query(1)
    real(real64), allocatable :: work(:)
    external :: dsyev

    ! Trivial case: return empty vector if empty matrix is given
    if (size(eivals) == 0) then
      return
    end if

    ! Workspace size query
    work_size = -1
    call dsyev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work_size_query, work_size, info)
    call check_info(info)

    ! Actual solution
    work_size = int(work_size_query(1))
    allocate(work(work_size))
    call dsyev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work, work_size, info)
    call check_info(info)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Allocatable version of the main subroutine.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine lapack_eigensolver_real(matrix, eivals)
    real(real64), intent(inout) :: matrix(:, :)
    real(real64), allocatable, intent(out) :: eivals(:)

    allocate(eivals(size(matrix, 1)))
    call lapack_eigensolver_real_plain(matrix, eivals)
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
    external :: zheev

    allocate(eivals(size(matrix, 1)))
    allocate(rwork(3*size(matrix, 1) - 2))

    ! Trivial case: return empty vector if empty matrix is given
    if (size(eivals) == 0) then
      return
    end if

    ! Workspace size query
    work_size = -1
    call zheev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work_size_query, work_size, rwork, info)
    call check_info(info)

    ! Actual solution
    work_size = int(work_size_query(1))
    allocate(work(work_size))
    call zheev('V', 'U', size(matrix, 1), matrix, size(matrix, 1), eivals, work, work_size, rwork, info)
    call check_info(info)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates LU-factorization of a matrix (dgetrf). Both L and U matrices are written into a single matrix, ignoring unit diagonal of L.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine lapack_lu_factorization(matrix, lu_factors, pivot_indices)
    real(real64), intent(in) :: matrix(:, :)
    real(real64), allocatable, intent(out) :: lu_factors(:, :)
    integer, allocatable, intent(out) :: pivot_indices(:)
    integer :: info
    external :: dgetrf

    lu_factors = matrix
    allocate(pivot_indices(min(size(lu_factors, 1), size(lu_factors, 2))))
    call dgetrf(size(lu_factors, 1), size(lu_factors, 2), lu_factors, size(lu_factors, 1), pivot_indices, info)
    call check_info(info)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates inverse of a matrix (dgetri).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function lapack_matrix_inverse(matrix) result(inverse)
    real(real64), intent(in) :: matrix(:, :)
    real(real64), allocatable :: inverse(:, :)
    integer :: work_size, info
    integer, allocatable :: pivot_indices(:)
    real(real64) :: work_size_query(1)
    real(real64), allocatable :: work(:)
    external :: dgetri

    call lapack_lu_factorization(matrix, inverse, pivot_indices)
    ! Work size query
    work_size = -1
    call dgetri(size(inverse, 1), inverse, size(inverse, 1), pivot_indices, work_size_query, work_size, info)
    call check_info(info)

    ! Actual solution
    work_size = int(work_size_query(1))
    allocate(work(work_size))
    call dgetri(size(inverse, 1), inverse, size(inverse, 1), pivot_indices, work, work_size, info)
    call check_info(info)
  end function

end module
