#include "funcs.macro"
  use overlaps_base_mod
  implicit none

  interface transform_basis_1d_to_fbr
    module procedure :: CONCAT2(transform_basis_1d_to_fbr_,TEMPLATE_TYPE_NAME)
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Takes a 2D solution *vec2* expressed over 1D basis and expresses it over FBR of sin/cos, using information about 1D basis given in *nvec1* and *vec1*.
! nvec1 - Number of 1D solutions in each thread (theta-slice).
! vec1 - All 1D solutions from all threads of a given rho, expressed over sin/cos. i-th element contains 1D solutions from i-th thread.
! vec2 - 2D solution expressed over 1D solutions. Expansion coefficients over solutions in all threads are stacked together in a single vector.
! vec2_fbr - The same 2D solution expressed in the basis of sin/cos.
! vec2_fbr consists of L blocks of length M each. Each L-block contains expansion coefficient over the same sin/cos function in different ls.
! Each element of the vector is sum over i of a_nlm^i * b_nli^j (j is index of vec2 and is fixed within this procedure).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(transform_basis_1d_to_fbr_,TEMPLATE_TYPE_NAME)(nvec1, vec1, vec2) result(vec2_fbr)
    integer, intent(in) :: nvec1(:)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), intent(in) :: vec1(:)
    TEMPLATE_TYPE, intent(in) :: vec2(:)
    TEMPLATE_TYPE, allocatable :: vec2_fbr(:)
    integer :: basis_size_phi, i, j, theta_ind

    basis_size_phi = size(vec1(1) % p, 1) ! All elements of vec1 are matrices with the same number of rows (basis size phi)
    allocate(vec2_fbr(size(nvec1) * basis_size_phi))
    i = 0
    j = 0
    do theta_ind = 1, size(nvec1)
      ! The case of nvec1(theta_ind) == 0 is defined and produces a vector of 0
      vec2_fbr(j+1 : j+basis_size_phi) = matmul(vec1(theta_ind) % p, vec2(i+1 : i+nvec1(theta_ind)))
      i = i + nvec1(theta_ind)
      j = j + basis_size_phi
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(calculate_overlap_block_,TEMPLATE_TYPE_NAME)(params, rho_ind_row, rho_ind_col, num_points_theta, sym_path) result(overlap_block)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind_row, rho_ind_col, num_points_theta
    character(*), intent(in) :: sym_path
    TEMPLATE_TYPE, allocatable :: overlap_block(:, :)
    integer :: i, j
    integer, allocatable :: num_solutions_1d_col(:), num_solutions_1d_row(:)
    real(real64), allocatable :: energies_2d_col(:), energies_2d_row(:)
    TEMPLATE_TYPE, allocatable :: solution_2d_fbr_col(:), solution_2d_fbr_row(:)
    TEMPLATE_TYPE, allocatable :: solutions_2d_col(:, :), solutions_2d_row(:, :)
    type(array_1d_real), allocatable :: energies_1d_col(:), energies_1d_row(:)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable :: solutions_1d_col(:), solutions_1d_row(:)

    ! Load bases
    call load_solutions_1D(get_solutions_1d_path(sym_path, rho_ind_col), num_points_theta, params % basis % num_funcs_phi_per_sym, num_solutions_1d_col, energies_1d_col, solutions_1d_col)
    call load_solutions_2D(get_solutions_2d_path(sym_path, rho_ind_col), energies_2d_col, solutions_2d_col)
    call load_solutions_1D(get_solutions_1d_path(sym_path, rho_ind_row), num_points_theta, params % basis % num_funcs_phi_per_sym, num_solutions_1d_row, energies_1d_row, solutions_1d_row)
    call load_solutions_2D(get_solutions_2d_path(sym_path, rho_ind_row), energies_2d_row, solutions_2d_row)

    ! Calculate overlap matrix
    allocate(overlap_block(size(solutions_2d_row, 2), size(solutions_2d_col, 2)))
    do j = 1, size(overlap_block, 2)
      solution_2d_fbr_col = transform_basis_1d_to_fbr(num_solutions_1d_col, solutions_1d_col, solutions_2d_col(:, j))
      do i = 1, size(overlap_block, 1)
        solution_2d_fbr_row = transform_basis_1d_to_fbr(num_solutions_1d_row, solutions_1d_row, solutions_2d_row(:, i))
        overlap_block(i, j) = dot_product(solution_2d_fbr_row, solution_2d_fbr_col)
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates overlap blocks in upper triange of Hamiltonian and saves them on disk in binary form.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(calculate_overlaps_,TEMPLATE_TYPE_NAME)(params, num_points_theta)
    class(input_params), intent(in) :: params
    integer, intent(in) :: num_points_theta
    integer :: proc_id, num_procs, rho_ind_col, rho_ind_row, block_ind, file_unit
    integer, allocatable :: num_solutions_2d(:)
    TEMPLATE_TYPE, allocatable :: overlap_block(:, :)
    character(:), allocatable :: sym_path

    proc_id = get_proc_id()
    num_procs = get_num_procs()
    sym_path = get_sym_path(params)
    num_solutions_2d = load_basis_size_2d(get_block_info_path(sym_path))

    block_ind = 0
    do rho_ind_row = 1, size(num_solutions_2d)
      if (num_solutions_2d(rho_ind_row) == 0) then
        cycle
      end if

      do rho_ind_col = rho_ind_row + 1, size(num_solutions_2d)
        if (num_solutions_2d(rho_ind_col) == 0) then
          cycle
        end if

        ! Assign process
        block_ind = block_ind + 1
        if (mod(block_ind - 1, num_procs) /= proc_id) then
          cycle
        end if

        ! Calculate and save block
        overlap_block = CONCAT2(calculate_overlap_block_,TEMPLATE_TYPE_NAME)(params, rho_ind_row, rho_ind_col, num_points_theta, sym_path)
        open(newunit = file_unit, file = get_regular_overlap_file_path(sym_path, rho_ind_row, rho_ind_col), form = 'unformatted')
        write(file_unit) overlap_block
        close(file_unit)
      end do
    end do
  end subroutine

