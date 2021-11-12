#include "funcs.macro"
  use overlaps_extra_base_mod
  implicit none

  private
  public :: CONCAT2(calculate_overlaps_extra_,TEMPLATE_TYPE_NAME)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates symmetric term overlaps.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(calculate_sym_term_,TEMPLATE_TYPE_NAME)(params, rho_grid, theta_grid)
    type(input_params), intent(in) :: params
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)
    integer :: my_id, n_procs, n, l, ir, ic, num_funcs_phi, file_unit
    real(real64) :: mu ! reduced mass
    ! Arrays for data arrays length
    integer, allocatable :: num_solutions_1d(:) ! number of 1D solutions in each theta slice (for fixed rho)
    integer, allocatable :: num_solutions_2d(:) ! number of 2D solutions in each rho slice
    ! Arrays for 1D eigenstates
    type(array_1d_real), allocatable :: energies_1d(:)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable :: exp_coeffs_1d(:) ! Expansion coefficients of all 1D solutions (in a single rho-slice) over primitive basis set (sin/cos)
    ! Arrays for 2D eigenstates
    real(real64), allocatable :: energies_2d(:)
    TEMPLATE_TYPE, allocatable :: exp_coeffs_2d_prim_row(:), exp_coeffs_2d_prim_col(:) ! Expansion coefficients of a 2D solution over primitive basis set (sin/cos)
    TEMPLATE_TYPE, allocatable :: exp_coeffs_2d(:, :) ! Expansion coefficients of all 2D solutions over 1D solutions (in a single rho-slice)
    ! Calculation of sym term
    real(real64), allocatable :: sym_factors_J(:), sym_factors_K(:)
    TEMPLATE_TYPE, allocatable :: partial_sum(:)
    TEMPLATE_TYPE, allocatable :: sym_block_J(:, :), sym_block_K(:, :)
    character(:), allocatable :: sym_folder, block_info_path, sym_block_J_path, sym_block_K_path
    character(:), allocatable :: solutions_1d_path, solutions_2d_path

    call print_parallel('Calculating sym term')
    my_id = get_proc_id()
    n_procs = get_num_procs()
    num_funcs_phi = params % get_num_funcs_phi_total()
    mu = params % get_reduced_mass()
    allocate(sym_factors_J(size(theta_grid)), sym_factors_K(size(theta_grid)), partial_sum(size(theta_grid)))

    sym_folder = get_sym_path(params)
    block_info_path = get_block_info_path(sym_folder)
    num_solutions_2d = load_basis_size_2d(block_info_path)

    ! Loop over diagonal blocks
    do n = my_id + 1, size(rho_grid), n_procs
      ! Skip if no basis
      if (num_solutions_2d(n) == 0) then
        cycle
      end if

      ! Load solutions
      solutions_1d_path = get_solutions_1d_path(sym_folder, n)
      call load_solutions_1D(solutions_1d_path, num_solutions_1d, energies_1d, exp_coeffs_1d)
      solutions_2d_path = get_solutions_2d_path(sym_folder, n)
      call load_solutions_2D(solutions_2d_path, energies_2d, exp_coeffs_2d)

      ! Compute sym operator values
      do l = 1, size(theta_grid)
        ! An + Bn
        sym_factors_J(l) = get_rotational_a(mu, rho_grid(n), theta_grid(l)) + get_rotational_b(mu, rho_grid(n), theta_grid(l))
        ! Cn - (An + Bn) / 2
        sym_factors_K(l) = get_rotational_c(mu, rho_grid(n), theta_grid(l)) - sym_factors_J(l) / 2d0
      end do

      allocate(sym_block_J(num_solutions_2d(n), num_solutions_2d(n)))
      allocate(sym_block_K(num_solutions_2d(n), num_solutions_2d(n)))
      ! Compute elements of sym block
      do ic = 1, num_solutions_2d(n)
        exp_coeffs_2d_prim_col = transform_basis_1d_to_fbr(num_solutions_1d, exp_coeffs_1d, exp_coeffs_2d(:, ic))

        do ir = 1, num_solutions_2d(n)
          exp_coeffs_2d_prim_row = transform_basis_1d_to_fbr(num_solutions_1d, exp_coeffs_1d, exp_coeffs_2d(:, ir))
          ! Compute partial sums
          do l = 1, size(theta_grid)
            partial_sum(l) = dot_product(exp_coeffs_2d_prim_row((l - 1)*num_funcs_phi + 1 : l*num_funcs_phi), exp_coeffs_2d_prim_col((l - 1)*num_funcs_phi + 1 : l*num_funcs_phi))
          end do
          sym_block_J(ir, ic) = dot_product(sym_factors_J, partial_sum)
          sym_block_K(ir, ic) = dot_product(sym_factors_K, partial_sum)
        end do
      end do

      ! Save sym blocks
      sym_block_J_path = get_symmetric_overlap_J_path(sym_folder, n)
      open(newunit=file_unit, file=sym_block_J_path, form='unformatted')
      write(file_unit) sym_block_J
      close(file_unit)
      deallocate(sym_block_J)

      sym_block_K_path = get_symmetric_overlap_K_path(sym_folder, n)
      open(newunit=file_unit, file=sym_block_K_path, form='unformatted')
      write(file_unit) sym_block_K
      close(file_unit)
      deallocate(sym_block_K)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates coriolis term overlaps for upper diagonal blocks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(calculate_coriolis_term_,TEMPLATE_TYPE_NAME)(params, rho_grid, theta_grid)
    type(input_params), intent(in) :: params
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)
    integer :: my_id, n_procs, nphi_per_basis_type, nphi_total, K_row, K_col, sym_row, sym_col, m_shift, m_ind, m_sign, n, m, l, ir, ic, file_unit
    integer :: m_range_row(2), m_range_col(2), m_range_global_row(2), m_range_global_col(2)
    real(real64) :: mu ! reduced mass
    ! Arrays for data length
    integer, allocatable :: num_solutions_1d_row(:), num_solutions_1d_col(:) ! number of 1D solutions in each theta slice (for fixed rho)
    integer, allocatable :: num_solutions_2d_row(:), num_solutions_2d_col(:) ! number of 2D solutions in each rho slice
    ! Arrays for 1D eigenstates
    type(array_1d_real), allocatable :: energies_1d_row(:), energies_1d_col(:)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable :: exp_coeffs_1d_row(:), exp_coeffs_1d_col(:) ! Expansion coefficients of all 1D solutions over primitive basis set (sin/cos)
    ! Arrays for 2D eigenstates
    real(real64), allocatable :: energies_2d_row(:), energies_2d_col(:)
    TEMPLATE_TYPE, allocatable :: exp_coeffs_2d_prim_row(:), exp_coeffs_2d_prim_col(:) ! Expansion coefficients of a 2D solution over primitive basis set (sin/cos)
    TEMPLATE_TYPE, allocatable :: exp_coeffs_2d_row(:, :), exp_coeffs_2d_col(:, :) ! Expansion coefficients of all 2D solutions over 1D solutions
    ! Calculation of Coriolis term
    real(real64), allocatable :: m_values(:), cor_factors_b(:)
    TEMPLATE_TYPE, allocatable :: partial_sum(:)
    TEMPLATE_TYPE, allocatable :: cor_block(:, :)
    character(:), allocatable :: sym_folder_row, sym_folder_col, block_info_path_row, block_info_path_col, cor_block_path
    character(:), allocatable :: solutions_1d_path_row, solutions_1d_path_col, solutions_2d_path_row, solutions_2d_path_col

    call print_parallel('Calculating Coriolis term')
    my_id = get_proc_id()
    n_procs = get_num_procs()
    nphi_per_basis_type = params % get_num_funcs_phi_per_basis_type()
    nphi_total = params % get_num_funcs_phi_total()
    mu = params % get_reduced_mass()

    K_row = params % K(1)
    K_col = K_row + 1
    sym_folder_row = get_sym_path_smart(params, K_row)
    sym_folder_col = get_sym_path_smart(params, K_col)
    block_info_path_row = get_block_info_path(sym_folder_row)
    block_info_path_col = get_block_info_path(sym_folder_col)
    num_solutions_2d_row = load_basis_size_2d(block_info_path_row)
    num_solutions_2d_col = load_basis_size_2d(block_info_path_col)

    ! Cosines have m values from 0 to M-1, sines from 1 to M, so they overlap only at m from 1 to M-1. We need to shift indices accordingly.
    ! 1,2m symmetries do not require shift
    sym_row = get_k_symmetry(K_row, params)
    sym_col = get_k_symmetry(K_col, params)
    m_shift = iff(any(sym_row == [2, 3]), 0, 1)
    m_range_row(1) = iff(sym_row == 0, 2, 1)
    m_range_row(2) = iff(m_shift == 1 .and. m_range_row(1) == 1, nphi_per_basis_type - 1, nphi_per_basis_type)
    m_range_col(1) = iff(m_shift == 1 .and. m_range_row(1) == 1, 2, 1)
    m_range_col(2) = iff(m_shift == 1 .and. m_range_col(1) == 1, nphi_per_basis_type - 1, nphi_per_basis_type)
    if (params % use_geometric_phase == 1) then
      m_range_col = m_range_col + nphi_per_basis_type
    end if

    ! Prepare m-values
    allocate(m_values(nphi_per_basis_type - m_shift), cor_factors_b(size(theta_grid)), partial_sum(size(theta_grid)))
    do m_ind = m_range_row(1), m_range_row(2)
      call get_m_ind_info(m_ind, nphi_per_basis_type, sym_row, params % get_molecule_type(), m = m)
      m_values(m_ind - m_range_row(1) + 1) = m
    end do
    ! Sign of factor depends on which function we take derivative from. Always positive for the first half (cos), when geometric phase is enabled.
    m_sign = iff(any(sym_row == [0, 2]), 1, -1)
    m_values = m_values * m_sign

    ! Loop over diagonal blocks
    do n = my_id + 1, size(rho_grid), n_procs
      ! Skip if no basis
      if (num_solutions_2d_row(n) == 0 .or. num_solutions_2d_col(n) == 0) then
        cycle
      end if

      ! Load solutions
      solutions_1d_path_row = get_solutions_1d_path(sym_folder_row, n)
      call load_solutions_1D(solutions_1d_path_row, num_solutions_1d_row, energies_1d_row, exp_coeffs_1d_row)
      solutions_2d_path_row = get_solutions_2d_path(sym_folder_row, n)
      call load_solutions_2D(solutions_2d_path_row, energies_2d_row, exp_coeffs_2d_row)

      solutions_1d_path_col = get_solutions_1d_path(sym_folder_col, n)
      call load_solutions_1D(solutions_1d_path_col, num_solutions_1d_col, energies_1d_col, exp_coeffs_1d_col)
      solutions_2d_path_col = get_solutions_2d_path(sym_folder_col, n)
      call load_solutions_2D(solutions_2d_path_col, energies_2d_col, exp_coeffs_2d_col)

      ! Compute b-factors
      do l = 1, size(theta_grid)
        cor_factors_b(l) = get_rotational_b(mu, rho_grid(n), theta_grid(l)) * cos(theta_grid(l))
      end do

      allocate(cor_block(num_solutions_2d_row(n), num_solutions_2d_col(n)))
      ! Compute elements of coriolis block
      do ic = 1, num_solutions_2d_col(n)
        exp_coeffs_2d_prim_col = transform_basis_1d_to_fbr(num_solutions_1d_col, exp_coeffs_1d_col, exp_coeffs_2d_col(:, ic))
        do ir = 1, num_solutions_2d_row(n)
          exp_coeffs_2d_prim_row = transform_basis_1d_to_fbr(num_solutions_1d_row, exp_coeffs_1d_row, exp_coeffs_2d_row(:, ir))
          ! Compute partial sums
          do l = 1, size(theta_grid)
            m_range_global_row = m_range_row + (l - 1) * nphi_total
            m_range_global_col = m_range_col + (l - 1) * nphi_total
            partial_sum(l) = dot_product(exp_coeffs_2d_prim_row(m_range_global_row(1) : m_range_global_row(2)), m_values * exp_coeffs_2d_prim_col(m_range_global_col(1) : m_range_global_col(2)))
            ! Add overlap between the other halves of the basis. Sign is negative because derivative is taken from cos in this case.
            if (params % use_geometric_phase == 1) then
              partial_sum(l) = partial_sum(l) - &
                dot_product(exp_coeffs_2d_prim_row(m_range_global_col(1) : m_range_global_col(2)), m_values * exp_coeffs_2d_prim_col(m_range_global_row(1) : m_range_global_row(2)))
            end if
          end do
          cor_block(ir, ic) = dot_product(cor_factors_b, partial_sum) * 2
        end do
      end do

      ! Save cor block
      cor_block_path = get_coriolis_overlap_path(sym_folder_row, n) ! Results are saved in the folder corresponding to K_row
      open(newunit=file_unit, file=cor_block_path, form='unformatted')
      write(file_unit) cor_block
      close(file_unit)
      deallocate(cor_block)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates asymmetric term overlaps.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(calculate_asym_term_,TEMPLATE_TYPE_NAME)(params, K_row, K_col, rho_grid, theta_grid)
    type(input_params), intent(in) :: params
    integer, intent(in) :: K_row, K_col ! pair of Ks to calculate asymmetric term
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)
    integer :: my_id, n_procs, n, l, ir, ic, num_funcs_phi, file_unit
    real(real64) :: mu ! reduced mass
    ! Arrays for data arrays length
    integer, allocatable :: num_solutions_1d_row(:), num_solutions_1d_col(:) ! number of 1D solutions in each theta slice (for fixed rho)
    integer, allocatable :: num_solutions_2d_row(:), num_solutions_2d_col(:) ! number of 2D solutions in each rho slice
    ! Arrays for 1D eigenstates
    type(array_1d_real), allocatable :: energies_1d_row(:), energies_1d_col(:)
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable :: exp_coeffs_1d_row(:), exp_coeffs_1d_col(:) ! Expansion coefficients of all 1D solutions over primitive basis set (sin/cos)
    ! Arrays for 2D eigenstates
    real(real64), allocatable :: energies_2d_row(:), energies_2d_col(:)
    TEMPLATE_TYPE, allocatable :: exp_coeffs_2d_prim_row(:), exp_coeffs_2d_prim_col(:) ! Expansion coefficients of a 2D solution over primitive basis set (sin/cos)
    TEMPLATE_TYPE, allocatable :: exp_coeffs_2d_row(:, :), exp_coeffs_2d_col(:, :) ! Expansion coefficients of all 2D solutions over 1D solutions
    ! Calculation of asym term
    real(real64), allocatable :: asym_factors(:)
    TEMPLATE_TYPE, allocatable :: partial_sum(:)
    TEMPLATE_TYPE, allocatable :: asym_block(:, :)
    character(:), allocatable :: sym_folder_row, sym_folder_col, block_info_path_row, block_info_path_col, asym_block_path
    character(:), allocatable :: solutions_1d_path_row, solutions_1d_path_col, solutions_2d_path_row, solutions_2d_path_col

    call print_parallel('Calculating asym term')
    call assert(K_row + 2 == K_col .or. K_row == 1 .and. K_col == 1, 'Wrong combination of Ks for asym term')
    my_id = get_proc_id()
    n_procs = get_num_procs()
    num_funcs_phi = params % get_num_funcs_phi_total()
    mu = params % get_reduced_mass()

    sym_folder_row = get_sym_path_smart(params, K_row)
    sym_folder_col = get_sym_path_smart(params, K_col)
    block_info_path_row = get_block_info_path(sym_folder_row)
    block_info_path_col = get_block_info_path(sym_folder_col)
    num_solutions_2d_row = load_basis_size_2d(block_info_path_row)
    num_solutions_2d_col = load_basis_size_2d(block_info_path_col)

    allocate(asym_factors(size(theta_grid)), partial_sum(size(theta_grid)))
    ! Loop over diagonal blocks
    do n = my_id + 1, size(rho_grid), n_procs
      ! Skip if no basis
      if (num_solutions_2d_row(n) == 0 .or. num_solutions_2d_col(n) == 0) then
        cycle
      end if

      ! Load solutions
      solutions_1d_path_row = get_solutions_1d_path(sym_folder_row, n)
      call load_solutions_1D(solutions_1d_path_row, num_solutions_1d_row, energies_1d_row, exp_coeffs_1d_row)
      solutions_2d_path_row = get_solutions_2d_path(sym_folder_row, n)
      call load_solutions_2D(solutions_2d_path_row, energies_2d_row, exp_coeffs_2d_row)

      solutions_1d_path_col = get_solutions_1d_path(sym_folder_col, n)
      call load_solutions_1D(solutions_1d_path_col, num_solutions_1d_col, energies_1d_col, exp_coeffs_1d_col)
      solutions_2d_path_col = get_solutions_2d_path(sym_folder_col, n)
      call load_solutions_2D(solutions_2d_path_col, energies_2d_col, exp_coeffs_2d_col)

      ! Compute asym operator values
      do l = 1, size(theta_grid)
        ! An - Bn only. Division by 2 is carried out later in this function
        asym_factors(l) = get_rotational_a(mu, rho_grid(n), theta_grid(l)) - get_rotational_b(mu, rho_grid(n), theta_grid(l)) 
      end do

      allocate(asym_block(num_solutions_2d_row(n), num_solutions_2d_col(n)))
      ! Compute elements of asym block
      do ic = 1, num_solutions_2d_col(n)
        exp_coeffs_2d_prim_col = transform_basis_1d_to_fbr(num_solutions_1d_col, exp_coeffs_1d_col, exp_coeffs_2d_col(:, ic))
        do ir = 1, num_solutions_2d_row(n)
          exp_coeffs_2d_prim_row = transform_basis_1d_to_fbr(num_solutions_1d_row, exp_coeffs_1d_row, exp_coeffs_2d_row(:, ir))
          ! Compute partial sums
          do l = 1, size(theta_grid)
            partial_sum(l) = dot_product(exp_coeffs_2d_prim_row((l - 1)*num_funcs_phi + 1 : l*num_funcs_phi), exp_coeffs_2d_prim_col((l - 1)*num_funcs_phi + 1 : l*num_funcs_phi))
          end do
          asym_block(ir, ic) = dot_product(asym_factors, partial_sum) / 4
        end do
      end do

      ! Save asym block
      ! Results are saved in the folder corresponding to K_row
      if (K_row == 1 .and. K_col == 1 .and. params % basis % fixed % enabled == 0) then
        asym_block_path = get_asymmetric_overlap_1_path(sym_folder_row, n)
      else
        asym_block_path = get_asymmetric_overlap_path(sym_folder_row, n)
      end if
      open(newunit=file_unit, file=asym_block_path, form='unformatted')
      write(file_unit) asym_block
      close(file_unit)
      deallocate(asym_block)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! The main subroutine that initiates calculation of both types of extra overlaps.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(calculate_overlaps_extra_,TEMPLATE_TYPE_NAME)(params, rho_grid, theta_grid)
    type(input_params), intent(in) :: params
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)

    call check_prerequisites(params)
    if (params % basis % fixed % enabled == 1) then
      call CONCAT2(calculate_sym_term_,TEMPLATE_TYPE_NAME)(params, rho_grid, theta_grid)
    end if

    ! Further blocks are only needed if Coriolis coupling is enabled
    if (params % use_rovib_coupling == 1) then
      ! Coriolis term
      if (params % K(1) + 1 <= params % J .or. params % basis % fixed % enabled == 1) then
        call CONCAT2(calculate_coriolis_term_,TEMPLATE_TYPE_NAME)(params, rho_grid, theta_grid)
      end if

      ! Asym term for (1, 1)-block or fixed basis
      if (params % K(1) == 1 .or. params % basis % fixed % enabled == 1) then
        call CONCAT2(calculate_asym_term_,TEMPLATE_TYPE_NAME)(params, params % K(1), params % K(1), rho_grid, theta_grid)
      end if

      ! Offdiag asym term
      if (params % K(1) + 2 <= params % J .and. params % basis % fixed % enabled == 0) then
        call CONCAT2(calculate_asym_term_,TEMPLATE_TYPE_NAME)(params, params % K(1), params % K(1) + 2, rho_grid, theta_grid)
      end if
    end if
  end subroutine

