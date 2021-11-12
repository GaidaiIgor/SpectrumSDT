#include "funcs.macro"
  use spectrumsdt_io_base_mod
  implicit none

  interface load_solutions_1D
    module procedure :: CONCAT2(load_solutions_1D_,TEMPLATE_TYPE_NAME)
  end interface

  interface load_1D_expansion_coefficients_K
    module procedure :: CONCAT2(load_1D_expansion_coefficients_K_,TEMPLATE_TYPE_NAME)
  end interface

  interface load_1D_expansion_coefficients
    module procedure :: CONCAT2(load_1D_expansion_coefficients_,TEMPLATE_TYPE_NAME)
  end interface

  interface load_solutions_2D
    module procedure :: CONCAT2(load_solutions_2D_,TEMPLATE_TYPE_NAME)
  end interface

  interface load_2D_expansion_coefficients_K
    module procedure :: CONCAT2(load_2D_expansion_coefficients_K_,TEMPLATE_TYPE_NAME)
  end interface

  interface load_2D_expansion_coefficients
    module procedure :: CONCAT2(load_2D_expansion_coefficients_,TEMPLATE_TYPE_NAME)
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads 1D solutions (eigenpairs).
! solutions_1d_path - path to the file with eigenvalues and eigenvectors of 1D solutions for specific rho index (n) and K.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(load_solutions_1D_,TEMPLATE_TYPE_NAME)(solutions_1d_path, num_solutions_1d, energies_1d, exp_coeffs_1d)
    character(*), intent(in) :: solutions_1d_path
    integer, allocatable, intent(out) :: num_solutions_1d(:) ! Num of 1D solutions (S_Knl) in each thread (l)
    type(array_1d_real), allocatable, intent(out) :: energies_1d(:) ! 1D solutions eigenvalues for each thread (l)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable, intent(out) :: exp_coeffs_1d(:) ! 1D solutions expansion coefficients (over sin/cos) for each thread (l)
    integer :: file_unit, theta_size, num_funcs_phi, i

    open(newunit = file_unit, file = solutions_1d_path, form = 'unformatted')
    read(file_unit) theta_size, num_funcs_phi
    allocate(num_solutions_1d(theta_size), energies_1d(theta_size), exp_coeffs_1d(theta_size))
    read(file_unit) num_solutions_1d
    do i = 1, theta_size
      allocate(energies_1d(i) % p(num_solutions_1d(i)), exp_coeffs_1d(i) % p(num_funcs_phi, num_solutions_1d(i)))
      read(file_unit) energies_1d(i) % p
      read(file_unit) exp_coeffs_1d(i) % p
    end do
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads 1D expansion coefficients of a specific K and sym stored at *sym_path* into *K_ind* of *num_solutions_1d* and *As*.
! N - basis size in rho.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(load_1D_expansion_coefficients_K_,TEMPLATE_TYPE_NAME)(sym_path, K_ind, N, num_solutions_1d, As)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: K_ind, N
    integer, intent(inout) :: num_solutions_1d(:, :, :)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), intent(inout) :: As(:, :, :)
    integer :: n_val, nphi
    integer, allocatable :: num_solutions_1d_ls(:)
    character(:), allocatable :: solutions_path ! path to file with all 1D solutions for all l-s for given K and n
    type(array_1d_real), allocatable :: energies_1d(:)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable :: As_slice(:) ! slice along l of the global array As

    do n_val = 1, N
      solutions_path = get_solutions_1d_path(sym_path, n_val)
      call load_solutions_1D(solutions_path, num_solutions_1d_ls, energies_1d, As_slice)
      num_solutions_1d(K_ind, n_val, :) = num_solutions_1d_ls
      As(K_ind, n_val, :) = As_slice
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 1D expansion coefficients in adiabatic basis mode.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(load_1D_expansion_coefficients_,TEMPLATE_TYPE_NAME)(params, N, L, As, num_solutions_1d)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N, L ! Num of points along rho and theta
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable, intent(out) :: As(:, :, :) ! K x N x L. Inner matrix is M x S_Knl
    integer, allocatable, intent(out) :: num_solutions_1d(:, :, :)
    integer :: last_K_load, total_Ks_load, K, K_ind
    character(:), allocatable :: sym_path

    last_K_load = params % get_last_k_load()
    total_Ks_load = last_K_load - params % K(1) + 1
    allocate(As(total_Ks_load, N, L), num_solutions_1d(total_Ks_load, N, L))
    do K = params % K(1), last_K_load
      K_ind = get_k_ind_smart(K, params)
      sym_path = get_sym_path_smart(params, K)
      call load_1D_expansion_coefficients_K(sym_path, K_ind, N, num_solutions_1d, As)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads 2D solutions.
! solutions_2d_path - path to the file with 2D solutions for specific K and n.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(load_solutions_2D_,TEMPLATE_TYPE_NAME)(solutions_2d_path, energies_2d, exp_coeffs_2d)
    character(*), intent(in) :: solutions_2d_path
    real(real64), allocatable, intent(out) :: energies_2d(:)
    TEMPLATE_TYPE, allocatable, intent(out) :: exp_coeffs_2d(:, :) ! 2D solutions expansion coefficients over 1D solutions
    integer :: file_unit, num_solutions, exp_size ! exp_size - number of coefficients in the basis expansion over 1D solutions

    open(newunit = file_unit, file = solutions_2d_path, form = 'unformatted')
    read(file_unit) num_solutions, exp_size
    allocate(energies_2d(num_solutions), exp_coeffs_2d(exp_size, num_solutions))
    read(file_unit) energies_2d
    read(file_unit) exp_coeffs_2d
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 2D expansion coefficients for a given K at *sym_path* into *K_ind* of corresponding arrays.
! N, L - number of points along rho and theta.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(load_2D_expansion_coefficients_K_,TEMPLATE_TYPE_NAME)(sym_path, K_ind, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: K_ind, N, L
    integer, intent(in) :: num_solutions_1d(:, :, :)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), intent(inout) :: Bs(:, :, :)
    integer, intent(inout) :: num_solutions_2d(:, :)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable, optional, intent(inout) :: Bs_plain(:, :) ! all Ls are stacked together
    integer :: n_val, l_val, next_start
    real(real64), allocatable :: energies_2d(:)
    TEMPLATE_TYPE, allocatable :: Bs_raw(:, :) ! Unarranged matrix of expansion coefficients. Single element of Bs_plain.
    character(:), allocatable :: solutions_path ! path to file with all 2D solutions for given K and n

    do n_val = 1, N
      solutions_path = get_solutions_2d_path(sym_path, n_val)
      call load_solutions_2D(solutions_path, energies_2d, Bs_raw)

      ! Not all parameter sets keep >0 solutions
      if (.not. allocated(Bs_raw)) then
        cycle
      end if

      if (present(Bs_plain)) then
        Bs_plain(K_ind, n_val) = array_2d(Bs_raw)
      end if

      num_solutions_2d(K_ind, n_val) = size(energies_2d)
      next_start = 1
      do l_val = 1, L
        Bs(K_ind, n_val, l_val) = array_2d(Bs_raw(next_start : next_start + num_solutions_1d(K_ind, n_val, l_val) - 1, :))
        next_start = next_start + num_solutions_1d(K_ind, n_val, l_val)
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 2D expansion coefficients in adiabatic basis mode.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(load_2D_expansion_coefficients_,TEMPLATE_TYPE_NAME)(params, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N, L ! Num of points along rho and theta
    integer, intent(in) :: num_solutions_1d(:, :, :)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable, intent(out) :: Bs(:, :, :) ! K x N x L
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable, optional, intent(out) :: Bs_plain(:, :) ! K x N, all Ls are stacked in columns
    integer, allocatable, intent(out) :: num_solutions_2d(:, :)
    integer :: last_K_load, total_Ks_load, K, K_ind
    character(:), allocatable :: sym_path

    last_K_load = params % get_last_k_load()
    total_Ks_load = last_k_load - params % K(1) + 1
    allocate(Bs(total_Ks_load, N, L), num_solutions_2d(total_Ks_load, N))
    num_solutions_2d = 0
    if (present(Bs_plain)) then
      allocate(Bs_plain(total_Ks_load, N))
    end if

    do K = params % K(1), last_K_load
      K_ind = get_k_ind_smart(K, params)
      sym_path = get_sym_path_smart(params, K)
      call load_2D_expansion_coefficients_K(sym_path, K_ind, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads an overlap block for given values of Ks and slice indexes.
! overlap_type: 0 - regular, 10 - symmetric J, 11 - symmetric K, 2 - coriolis, 3 - asymmetric.
! *is_file_real* is 1 if overlap file is written with real data, otherwise complex is assumed.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(load_overlap_block_,TEMPLATE_TYPE_NAME)(params, K_row, K_col, overlap_type, slice_ind_row, slice_ind_col, rows, columns) result(block)
    class(input_params), intent(in) :: params
    integer, intent(in) :: K_row, K_col, overlap_type, slice_ind_row, slice_ind_col, rows, columns
    TEMPLATE_TYPE, allocatable :: block(:, :)
    logical :: swap_condition
    integer :: file_unit, K_row_load, K_col_load, slice_ind_row_load, slice_ind_col_load, rows_load, columns_load
    character(:), allocatable :: sym_path, file_path

    ! Load identity matrix for regular diagonal overlaps
    if (K_row == K_col .and. slice_ind_row == slice_ind_col .and. overlap_type == 0) then
      call assert(rows == columns, 'Assertion error: diagonal overlap is not square')
      block = identity_matrix(rows)
      return
    end if

    ! Since the matrix is symmetric only blocks above the main diagonal are actually stored, blocks below are just conjugated transposed upper diagonal blocks
    ! Prepare to load transposed version of block if it is below main diagonal
    swap_condition = K_row > K_col .or. slice_ind_row > slice_ind_col
    K_row_load = merge(K_col, K_row, swap_condition)
    K_col_load = merge(K_row, K_col, swap_condition)
    slice_ind_row_load = merge(slice_ind_col, slice_ind_row, swap_condition)
    slice_ind_col_load = merge(slice_ind_row, slice_ind_col, swap_condition)
    rows_load = merge(columns, rows, swap_condition)
    columns_load = merge(rows, columns, swap_condition)

    sym_path = get_sym_path_smart(params, K_row_load)
    if (overlap_type == 0) then
      file_path = get_regular_overlap_path(sym_path, slice_ind_row_load, slice_ind_col_load)
    else if (overlap_type == 10) then
      file_path = get_symmetric_overlap_J_path(sym_path, slice_ind_row_load)
    else if (overlap_type == 11) then
      file_path = get_symmetric_overlap_K_path(sym_path, slice_ind_row_load)
    else if (overlap_type == 2) then
      file_path = get_coriolis_overlap_path(sym_path, slice_ind_row_load)
    else if (overlap_type == 3) then
      if (K_row == 1 .and. K_col == 1) then
        file_path = get_asymmetric_overlap_1_path(sym_path, slice_ind_row_load)
      else
        file_path = get_asymmetric_overlap_path(sym_path, slice_ind_row_load)
      end if
    end if

    open(newunit = file_unit, file = file_path, form = 'unformatted')
    allocate(block(rows_load, columns_load))
    read(file_unit) block
    close(file_unit)

    ! transpose back once loaded
    if (swap_condition) then
      block = conjg(transpose(block))
    end if
  end function

