module slepc_solver_mod
#include <slepc/finclude/slepceps.h>
  use general_utils, only: num2str
  use matmul_operator_mod, only: active_matmul_operator, msize
  use parallel_utils, only: get_proc_id, get_proc_elem_range, print_parallel
  use slepceps

  use debug_tools

  private
  public :: find_eigenpairs_slepc

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Finds eigenpairs of a matrix currently set in matmul_operator_mod using SLEPc eigensolver
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine find_eigenpairs_slepc(n_eigs, ncv, mpd, eivals, eivecs)
    integer, intent(in) :: n_eigs
    integer, intent(in) :: ncv, mpd ! Convergence parameters. Both can be set to -1 to let SLEPc decide.
    complex*16, allocatable, intent(out) :: eivals(:)
    complex*16, allocatable, intent(out) :: eivecs(:, :)
    integer :: n_conv, proc_first_eivec, proc_eivecs, proc_first_row, proc_rows, root
    integer, allocatable :: eivec_counts(:), eivec_shifts(:), row_counts(:), row_shifts(:)
    complex*16 :: next_eivec_local(msize)
    PetscInt :: ncv_petsc, mpd_petsc, n_eigs_petsc, n_conv_petsc, i_petsc, msize_petsc
    PetscErrorCode :: ierr
    PetscScalar, pointer :: next_eivec_ptr(:)
    Mat :: A
    EPS :: eps
    DS :: ds
    Vec :: next_eivec

    msize_petsc = msize
    call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
    call MatCreateShell(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, msize_petsc, msize_petsc, PETSC_NULL_INTEGER, A, ierr)
    call MatShellSetOperation(A, MATOP_MULT, slepc_mult, ierr)

    call EPSCreate(PETSC_COMM_WORLD, eps, ierr)
    call EPSGetDS(eps, ds, ierr)
    call DSSetParallel(ds, DS_PARALLEL_SYNCHRONIZED, ierr)

    call EPSSetOperators(eps, A, PETSC_NULL_MAT, ierr)
    call EPSSetProblemType(eps, EPS_NHEP, ierr)
    call EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL, ierr)

    ncv_petsc = ncv
    mpd_petsc = mpd
    ncv_petsc = merge(PETSC_DEFAULT_INTEGER, ncv_petsc, ncv_petsc == -1)
    mpd_petsc = merge(PETSC_DEFAULT_INTEGER, mpd_petsc, mpd_petsc == -1)
    n_eigs_petsc = n_eigs
    call EPSSetDimensions(eps, n_eigs_petsc, ncv_petsc, mpd_petsc, ierr)

    call EPSSetFromOptions(eps, ierr)
    call EPSSolve(eps, ierr)
    call EPSGetConverged(eps, n_conv_petsc, ierr)
    n_conv = n_conv_petsc
    call print_parallel('Eigenpairs converged: ' // num2str(n_conv))

    ! Populate eigenvalues
    allocate(eivals(n_conv))
    do i_petsc = 0, n_conv - 1
      call EPSGetEigenpair(eps, i_petsc, eivals(i_petsc + 1), PETSC_NULL_SCALAR, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
    end do

    ! Populate eigenvectors
    call get_proc_elem_range(n_conv, proc_first_eivec, proc_eivecs, eivec_counts, eivec_shifts)
    call get_proc_elem_range(msize, proc_first_row, proc_rows, row_counts, row_shifts)
    allocate(eivecs(msize, proc_eivecs))
    call MatCreateVecs(A, next_eivec, PETSC_NULL_VEC, ierr)
    do i_petsc = 0, n_conv - 1
      call EPSGetEigenpair(eps, i_petsc, PETSC_NULL_SCALAR, PETSC_NULL_SCALAR, next_eivec, PETSC_NULL_VEC, ierr)
      call VecGetArrayReadF90(next_eivec, next_eivec_ptr, ierr)
      root = findloc(i_petsc + 1 > eivec_shifts, .true., 1, back=.true.) - 1
      call MPI_Gatherv(next_eivec_ptr, proc_rows, MPI_DOUBLE_COMPLEX, next_eivec_local, row_counts, row_shifts, MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr)
      if (get_proc_id() == root) then
        eivecs(:, i_petsc - proc_first_eivec + 2) = next_eivec_local
      end if
      call VecRestoreArrayReadF90(next_eivec, next_eivec_ptr, ierr)
    end do

    call EPSDestroy(eps, ierr)
    call MatDestroy(A, ierr)
    call VecDestroy(next_eivec, ierr)
    call SlepcFinalize(ierr)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! A wrapper around ham_mult_compressed
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine slepc_mult(A, x, y, ierr)
    Mat :: A ! intent(in)
    Vec :: x ! intent(in)
    Vec :: y ! intent(out)
    PetscErrorCode :: ierr
    integer :: size_proc
    PetscInt :: size_proc_petsc
    PetscScalar, pointer :: x_ptr(:), y_ptr(:)

    call VecGetLocalSize(x, size_proc_petsc, ierr)
    call VecGetArrayReadF90(x, x_ptr, ierr)
    call VecGetArrayF90(y, y_ptr, ierr)

    size_proc = size_proc_petsc
    call active_matmul_operator(size_proc, x_ptr, y_ptr)

    call VecRestoreArrayReadF90(x, x_ptr, ierr)
    call VecRestoreArrayF90(y, y_ptr, ierr)
  end subroutine
end module
