module slepc_solver_mod
#include <slepc/finclude/slepceps.h>
  use arnoldi_operator_mod, only: msize, ham_mult_compressed
  use parallel_utils
  use slepceps

  use constants
  use debug_tools

  private
  public :: find_eigenpairs_slepc
contains
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Finds eigenpairs of a matrix currently set in arnoldi_operator_mod using SLEPc eigensolver
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine find_eigenpairs_slepc(n_eigs, ncv, mpd, eivals, eivecs)
    integer, intent(in) :: n_eigs
    integer, intent(in) :: ncv, mpd ! Convergence parameters. Both can be set to -1 to let SLEPc decide.
    complex*16, allocatable, intent(out) :: eivals(:)
    complex*16, allocatable, intent(out) :: eivecs(:, :)
    integer :: proc_first_eivec, proc_eivecs, proc_first_row, proc_rows, root
    integer, allocatable :: eivec_counts(:), eivec_shifts(:), row_counts(:), row_shifts(:)
    complex*16 :: next_eivec_local(msize)
    PetscInt :: ncv_act, mpd_act, n_conv, i
    PetscErrorCode :: ierr
    PetscScalar, pointer :: next_eivec_ptr(:)
    Mat :: A
    EPS :: eps
    DS :: ds
    RG :: rg
    Vec :: next_eivec

    call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
    call MatCreateShell(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, msize, msize, PETSC_NULL_INTEGER, A, ierr)
    call MatShellSetOperation(A, MATOP_MULT, slepc_mult, ierr)

    call EPSCreate(PETSC_COMM_WORLD, eps, ierr)
    call EPSGetDS(eps, ds, ierr)
    call DSSetParallel(ds, DS_PARALLEL_SYNCHRONIZED, ierr)

    call EPSSetOperators(eps, A, PETSC_NULL_MAT, ierr)
    call EPSSetProblemType(eps, EPS_NHEP, ierr)

    call EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL, ierr)
    ncv_act = merge(PETSC_DEFAULT_INTEGER, ncv, ncv == -1)
    mpd_act = merge(PETSC_DEFAULT_INTEGER, mpd, mpd == -1)

    if (debug_mode == 'interval') then
      call EPSGetRG(eps, rg, ierr)
      call RGSetType(rg, RGINTERVAL, ierr)
      call RGIntervalSetEndpoints(rg, -500 / autown, 500 / autown, -1d0, 1d0, ierr)
      call EPSSetTarget(eps, (-10000, 0) / autown, ierr)
      call EPSSetWhichEigenpairs(eps, EPS_TARGET_REAL, ierr)
      ! call EPSSetDimensions(eps, 2*mpd, 3*mpd, mpd, ierr)
    else
      call EPSSetDimensions(eps, n_eigs, ncv_act, mpd_act, ierr)
    end if

    call EPSSetFromOptions(eps, ierr)
    call EPSSolve(eps, ierr)
    call EPSGetConverged(eps, n_conv, ierr)
    call print_parallel('Eigenpairs converged: ' // num2str(n_conv))

    ! Populate eigenvalues
    allocate(eivals(n_conv))
    do i = 0, n_conv - 1
      call EPSGetEigenpair(eps, i, eivals(i + 1), PETSC_NULL_SCALAR, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
    end do

    ! Populate eigenvectors
    call get_proc_elem_range(n_conv, proc_first_eivec, proc_eivecs, eivec_counts, eivec_shifts)
    call get_proc_elem_range(msize, proc_first_row, proc_rows, row_counts, row_shifts)
    allocate(eivecs(msize, proc_eivecs))
    call MatCreateVecs(A, next_eivec, PETSC_NULL_VEC, ierr)
    do i = 0, n_conv - 1
      call EPSGetEigenpair(eps, i, PETSC_NULL_SCALAR, PETSC_NULL_SCALAR, next_eivec, PETSC_NULL_VEC, ierr)
      call VecGetArrayReadF90(next_eivec, next_eivec_ptr, ierr)
      root = findloc(i + 1 > eivec_shifts, .true., 1, back=.true.) - 1
      call MPI_Gatherv(next_eivec_ptr, proc_rows, MPI_DOUBLE_COMPLEX, next_eivec_local, row_counts, row_shifts, MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr)
      if (get_proc_id() == root) then
        eivecs(:, i - proc_first_eivec + 2) = next_eivec_local
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
    Mat :: A
    Vec :: x, y
    PetscErrorCode :: ierr
    PetscInt :: size_proc
    PetscScalar, pointer :: x_ptr(:), y_ptr(:)

    call VecGetLocalSize(x, size_proc, ierr)
    call VecGetArrayReadF90(x, x_ptr, ierr)
    call VecGetArrayF90(y, y_ptr, ierr)

    call ham_mult_compressed(size_proc, x_ptr, y_ptr)

    call VecRestoreArrayReadF90(x, x_ptr, ierr)
    call VecRestoreArrayF90(y, y_ptr, ierr)
  end subroutine

end module
