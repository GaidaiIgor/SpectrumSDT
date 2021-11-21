#include "funcs.macro"
  use basis_base_mod

  use debug_tools_mod
  implicit none

  interface calc_1d
    module procedure :: CONCAT2(calc_1d_,TEMPLATE_TYPE_NAME)
  end interface

  interface build_hamiltonian_2d
    module procedure :: CONCAT2(build_hamiltonian_2d_,TEMPLATE_TYPE_NAME)
  end interface

  interface calc_2d
    module procedure :: CONCAT2(calc_2d_,TEMPLATE_TYPE_NAME)
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates 1D Hamiltonian (for phi).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(get_hamiltonian_1d_,TEMPLATE_TYPE_NAME)(params, rho_val, theta_val, grid_phi, potential_phi) result(ham)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: rho_val, theta_val
    real(real64), intent(in) :: grid_phi(:), potential_phi(:)
    TEMPLATE_TYPE, allocatable :: ham(:, :)
    integer :: nphi_per_basis_type, l_dash, m1_ind, m2_ind, m_shift
    real(real64) :: mu, step_phi, coeff, m
    real(real64), allocatable :: func(:)
    real(real64), allocatable :: basis(:, :)

    basis = get_phi_basis_grid(params, grid_phi)
    mu = params % get_reduced_mass()
    step_phi = grid_phi(2) - grid_phi(1)
    nphi_per_basis_type = params % get_num_funcs_phi_per_basis_type()
    l_dash = 3

    allocate(ham(size(basis, 2), size(basis, 2)))
    ham = 0

    ! Build potential energy matrix
    do m2_ind = 1, size(ham, 2)
      m_shift = iff(m2_ind > nphi_per_basis_type, nphi_per_basis_type, 0) ! Skip empty quadrant when both symmetries are included
      do m1_ind = 1 + m_shift, m2_ind
        func = basis(:, m1_ind) * potential_phi * basis(:, m2_ind)
        ham(m1_ind, m2_ind) = integrate_1d(func, step_phi)
        ham(m2_ind, m1_ind) = ham(m1_ind, m2_ind)
      end do
    end do

    ! Add kinetic energy matrix
    coeff = -2 / mu / rho_val**2 / sin(theta_val)**2
    do m2_ind = 1, size(ham, 2)
      call get_m_ind_info(m2_ind, params, m = m)
      ham(m2_ind, m2_ind) = ham(m2_ind, m2_ind) - coeff * m ** 2

#if TYPE_ID == COMPLEX_ID
      if (params % use_geometric_phase == 1) then
        ham(m2_ind, m2_ind) = ham(m2_ind, m2_ind) - coeff * l_dash**2 / 4
        if (m2_ind > nphi_per_basis_type .and. m2_ind < size(ham, 2)) then
          m1_ind = m + 1
          ham(m1_ind, m2_ind) = coeff * (0, 1d0) * l_dash * m
          ham(m2_ind, m1_ind) = conjg(ham(m1_ind, m2_ind))
        end if
      end if
#endif

    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Solves 1D problems for each thread in slice specified by *rho_ind*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(calc_1d_,TEMPLATE_TYPE_NAME)(params, rho_ind, rho_val, grid_theta, grid_phi, potential, nvec1, val1, vec1)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind
    real(real64), intent(in) :: rho_val
    real(real64), intent(in) :: grid_theta(:), grid_phi(:)
    real(real64), intent(in) :: potential(:, :, :)
    integer, allocatable, intent(out) :: nvec1(:) ! Num of 1D vecs in each thread
    type(array_1d_real), allocatable, intent(out) :: val1(:) ! 1D values in each thread
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable, intent(out) :: vec1(:) ! 1D vectors in each thread
    integer :: theta_ind, file_unit
    real(real64), allocatable :: val1_all(:) ! All eigenvalues
    TEMPLATE_TYPE, allocatable :: ham(:, :) ! Hamiltonian

    allocate(nvec1(size(grid_theta)), val1(size(grid_theta)), vec1(size(grid_theta)))
    ! Solve eigenvalue problem for each thread
    do theta_ind = 1, size(grid_theta)
      ham = CONCAT2(get_hamiltonian_1d_,TEMPLATE_TYPE_NAME)(params, rho_val, grid_theta(theta_ind), grid_phi, potential(:, theta_ind, rho_ind))
      call lapack_eigensolver(ham, val1_all)

      ! Save results
      nvec1(theta_ind) = min(max(findloc(val1_all < params % basis % cutoff_energy_1d, .true., dim = 1, back = .true.), params % get_min_solutions_1d_total()), size(ham, 1))
      val1(theta_ind) % p = val1_all(:nvec1(theta_ind))
      vec1(theta_ind) % p = ham(:, :nvec1(theta_ind))
    end do

    ! Write results to a binary file
    open(newunit = file_unit, file = get_solutions_1d_path(get_sym_path(params), rho_ind), form = 'unformatted')
    write(file_unit) size(nvec1), size(ham, 1)
    write(file_unit) nvec1
    do theta_ind = 1, size(grid_theta)
      write(file_unit) val1(theta_ind) % p
      write(file_unit) vec1(theta_ind) % p
    end do
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Builds a 2D Hamiltonian in 1D basis in phi.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(build_hamiltonian_2d_,TEMPLATE_TYPE_NAME)(mu, rho_val, period_theta, nvec1, val1, vec1) result(ham2)
    real(real64), intent(in) :: mu, rho_val, period_theta
    integer, intent(in) :: nvec1(:) ! Number of 1D eigenpairs in each thread
    type(array_1d_real), intent(in) :: val1(:) ! 1D eigenvalues for each thread
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), intent(in) :: vec1(:) ! 1D eigenvectors for each thread
    TEMPLATE_TYPE, allocatable :: ham2(:, :)
    integer :: ham2_size, block_col, block_row, i
    integer, allocatable :: offset(:) ! Block offsets in 2D Hamiltonian matrix
    real(real64), allocatable :: kin(:, :)
    TEMPLATE_TYPE, allocatable :: block(:, :)

    offset = prefix_sum_exclusive(nvec1)
    kin = compute_kinetic_theta(mu, rho_val, period_theta, size(nvec1))
    ham2_size = sum(nvec1)
    allocate(ham2(ham2_size, ham2_size))

    do block_col = 1, size(nvec1)
      if (nvec1(block_col) == 0) then
        cycle
      end if

      do block_row = block_col, size(nvec1)
        if (nvec1(block_row) == 0) then
          cycle
        end if

        ! Calculate overlaps
        if (block_col == block_row) then
          block = identity_matrix(nvec1(block_col)) 
        else
          block = matmul(conjg(transpose(vec1(block_row) % p)), vec1(block_col) % p)
        end if

        ! Factor in kinetic energy
        block = block * kin(block_row, block_col)
        ! Factor in 1D eigenvalues. Note: eigenvalues should not be multiplied by kinetic energy.
        if (block_col == block_row) then
          do i = 1, nvec1(block_col)
            block(i, i) = block(i, i) + val1(block_col) % p(i)
          end do
        end if

        ! Write complete block
        ham2(offset(block_row)+1 : offset(block_row)+nvec1(block_row), offset(block_col)+1 : offset(block_col)+nvec1(block_col)) = block
        if (block_col /= block_row) then
          ham2(offset(block_col)+1 : offset(block_col)+nvec1(block_col), offset(block_row)+1 : offset(block_row)+nvec1(block_row)) = conjg(transpose(block))
        end if
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Solves 2D problem for *rho_ind*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(calc_2d_,TEMPLATE_TYPE_NAME)(params, rho_ind, rho_val, period_theta, nvec1, val1, vec1, nvec2, val2)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind
    real(real64), intent(in) :: rho_val, period_theta
    integer, intent(in) :: nvec1(:) ! Num of 1D vecs in each thread
    type(array_1d_real), intent(in) :: val1(:) ! 1D eigenvalues for each thread
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), intent(in) :: vec1(:) ! 1D eigenvectors for each thread
    integer, intent(out) :: nvec2
    real(real64), allocatable, intent(out) :: val2(:)
    integer :: file_unit
    real(real64) :: mu
    real(real64), allocatable :: val2_all(:)
    TEMPLATE_TYPE, allocatable :: ham2(:, :), vec2(:, :)

    mu = params % get_reduced_mass()
    ham2 = build_hamiltonian_2d(mu, rho_val, period_theta, nvec1, val1, vec1)
    call lapack_eigensolver(ham2, val2_all)

    nvec2 = min(max(findloc(val2_all < params % basis % cutoff_energy_2d, .true., dim = 1, back = .true.), params % get_min_solutions_2d_total()), size(ham2, 1))
    val2 = val2_all(:nvec2)
    vec2 = ham2(:, :nvec2)

    open(newunit = file_unit, file = get_solutions_2d_path(get_sym_path(params), rho_ind), form = 'unformatted')
    write(file_unit) nvec2, size(ham2, 1)
    write(file_unit) val2
    write(file_unit) vec2
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates 1D and 2D basis.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(calculate_basis_,TEMPLATE_TYPE_NAME)(params, grid_rho, grid_theta, grid_phi, potential)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: grid_rho(:), grid_theta(:), grid_phi(:)
    real(real64), intent(in) :: potential(:, :, :)
    integer :: proc_first, proc_rhos, rho_ind, nvec2, file_unit
    integer, allocatable :: nvec1(:), proc_nvec2(:), all_nvec2(:)
    integer, allocatable :: proc_nvec1(:, :), all_nvec1(:, :)
    real(real64) :: step_theta
    real(real64), allocatable :: val2(:)
    type(array_1d_real), allocatable :: val1(:), proc_val2(:)
    type(array_1d_real), allocatable :: proc_val1(:, :)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable :: vec1(:)

    call get_proc_elem_range(size(grid_rho), proc_first, proc_rhos)
    allocate(proc_nvec2(proc_rhos), proc_val2(proc_rhos))
    allocate(proc_nvec1(size(grid_theta), proc_rhos), proc_val1(size(grid_theta), proc_rhos))
    step_theta = grid_theta(2) - grid_theta(1)

    do rho_ind = proc_first, proc_first + proc_rhos - 1
      call calc_1d(params, rho_ind, grid_rho(rho_ind), grid_theta, grid_phi, potential, nvec1, val1, vec1)

      proc_nvec1(:, rho_ind - proc_first + 1) = nvec1
      proc_val1(:, rho_ind - proc_first + 1) = val1

      call calc_2d(params, rho_ind, grid_rho(rho_ind), step_theta, nvec1, val1, vec1, nvec2, val2)
      proc_nvec2(rho_ind - proc_first + 1) = nvec2
      proc_val2(rho_ind - proc_first + 1) % p = val2
    end do

    if (params % basis % print_energies_1d == 1) then
      call gather_write_energies_1d(proc_val1, get_energies_1d_path(get_sym_path(params)))
    end if
    if (params % basis % print_energies_2d == 1) then
      call gather_write_energies_2d(proc_val2, get_energies_2d_path(get_sym_path(params)))
    end if

    call gather_cols_array_2d_integer(proc_nvec1, all_nvec1)
    call gather_array_1d_integer(proc_nvec2, all_nvec2)

    ! Write total number of 1D and 2D vectors
    if (get_proc_id() == 0) then
      open(newunit = file_unit, file = get_basis_1D_summary_path(get_sym_path(params)))
      write(file_unit, '(' // num2str(size(grid_theta)) // '(I5))') all_nvec1
      close(file_unit)

      open(newunit = file_unit, file = get_block_info_path(get_sym_path(params)))
      write(file_unit, '(I5)') all_nvec2
      close(file_unit)

      print *, 'Total number of 1D basis functions: ', sum(all_nvec1)
      print *, 'Maximum number of 1D basis functions per slice: ', maxval(all_nvec1)
      print *, 'Total number of 2D basis functions: ', sum(all_nvec2)
      print *, 'Maximum number of 2D basis functions per slice: ', maxval(all_nvec2)
    end if
  end subroutine

