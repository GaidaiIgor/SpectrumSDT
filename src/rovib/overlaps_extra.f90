module overlaps_extra_mod
  ! Contains procedures related to calculation of extra overlap terms: asymmetric and coriolis
  use array_1d_mod
  use array_2d_mod
  use formulas_mod
  use input_params_mod
  use iso_fortran_env, only: real64
  use k_block_info
  use parallel_utils
  use path_utils
  use rovib_io_mod
  implicit none

  private
  public :: calculate_overlaps_extra

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transforms expansion coefficients of 2D functions over 1D functions to expansion coefficients over sin/cos functions.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function expansion_1D_to_primitive(num_solutions_1d, exp_coeffs_1d, exp_coeffs_2d) result(exp_coeffs_2d_prim)
    integer, intent(in) :: num_solutions_1d(:) ! number of 1D vectors in each theta slice
    type(array_2d_real), intent(in) :: exp_coeffs_1d(:) ! i-th element is a 2D array of 1D solutions in the i-th theta slice given as expansion coefficients over sin/cos basis
    real(real64), intent(in) :: exp_coeffs_2d(:) ! expansion coefficients of a given 2D state over 1D solutions
    real(real64), allocatable :: exp_coeffs_2d_prim(:) ! expansion coefficients of the same 2D state over sin/cos basis
    integer :: is, i, j ! is - theta slice index, i tracks current position in 1D expansion, j track current position in sin/cos output vector
    integer :: n_basis, n_prim_coeffs

    i = 0
    n_basis = size(exp_coeffs_1d(1) % p, 1) ! all 2D arrays have the same number of rows
    n_prim_coeffs = size(num_solutions_1d) * n_basis
    allocate(exp_coeffs_2d_prim(n_prim_coeffs))
    ! Iterate over theta slices
    do is = 1, size(num_solutions_1d)
      j = (is-1) * n_basis
      ! call gemv(exp_coeffs_1d(is) % p, exp_coeffs_2d(i+1 : i+num_solutions_1d(is)), exp_coeffs_2d_prim(j+1 : j+n_basis)) ! matrix vector product
      exp_coeffs_2d_prim(j+1 : j+n_basis) = matmul(exp_coeffs_1d(is) % p, exp_coeffs_2d(i+1 : i+num_solutions_1d(is)))
      i = i + num_solutions_1d(is)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates symmetric term overlaps.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_sym_term(params, mu, rho_grid, theta_grid)
    type(input_params), intent(in) :: params
    real(real64), intent(in) :: mu ! reduced mass
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)
    integer :: my_id, n_procs
    integer :: n, l, ir, ic, n_basis, file_unit
    ! Arrays for data arrays length
    integer, allocatable :: num_solutions_1d(:) ! number of 1D solutions in each theta slice (for fixed rho)
    integer, allocatable :: num_solutions_2d(:) ! number of 2D solutions in each rho slice
    ! Arrays for 1D eigenstates
    type(array_1d_real), allocatable :: energies_1d(:)
    type(array_2d_real), allocatable :: exp_coeffs_1d(:) ! Expansion coefficients of all 1D solutions over primitive basis set (sin/cos)
    ! Arrays for 2D eigenstates
    real(real64), allocatable :: energies_2d(:)
    real(real64), allocatable :: exp_coeffs_2d(:, :) ! Expansion coefficients of all 2D solutions over 1D solutions
    real(real64), allocatable :: exp_coeffs_2d_prim_row(:), exp_coeffs_2d_prim_col(:) ! Expansion coefficients of a 2D solution over primitive basis set (sin/cos)
    ! Calculation of asym term
    real(real64), allocatable :: sym_block_J(:, :), sym_block_K(:, :)
    real(real64), allocatable :: sym_factors_J(:), sym_factors_K(:), partial_sum(:)
    character(:), allocatable :: sym_folder, block_info_path, sym_block_J_path, sym_block_K_path
    character(:), allocatable :: solutions_1d_path, solutions_2d_path

    call print_parallel('Calculating sym term')
    my_id = get_proc_id()
    n_procs = get_num_procs()
    n_basis = params % basis % num_functions_phi
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
      call load_solutions_1D(solutions_1d_path, size(theta_grid), n_basis, num_solutions_1d, energies_1d, exp_coeffs_1d)
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
        exp_coeffs_2d_prim_col = expansion_1D_to_primitive(num_solutions_1d, exp_coeffs_1d, exp_coeffs_2d(:, ic))

        do ir = 1, num_solutions_2d(n)
          exp_coeffs_2d_prim_row = expansion_1D_to_primitive(num_solutions_1d, exp_coeffs_1d, exp_coeffs_2d(:, ir))
          ! Compute partial sums
          do l = 1, size(theta_grid)
            partial_sum(l) = dot_product(exp_coeffs_2d_prim_row((l - 1)*n_basis + 1 : l*n_basis), exp_coeffs_2d_prim_col((l - 1)*n_basis + 1 : l*n_basis))
          end do
          sym_block_J(ir, ic) = dot_product(sym_factors_J, partial_sum)
          sym_block_K(ir, ic) = dot_product(sym_factors_K, partial_sum)
        end do
      end do

      ! Save sym blocks
      sym_block_J_path = get_symmetric_overlap_J_file_path(sym_folder, n)
      open(newunit=file_unit, file=sym_block_J_path, form='unformatted')
      write(file_unit) sym_block_J
      close(file_unit)
      deallocate(sym_block_J)

      sym_block_K_path = get_symmetric_overlap_K_file_path(sym_folder, n)
      open(newunit=file_unit, file=sym_block_K_path, form='unformatted')
      write(file_unit) sym_block_K
      close(file_unit)
      deallocate(sym_block_K)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates coriolis term overlaps for upper diagonal blocks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_coriolis_term(params, mu, rho_grid, theta_grid)
    type(input_params), intent(in) :: params
    real(real64), intent(in) :: mu ! reduced mass
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)
    integer :: my_id, n_procs
    integer :: n, m, l, ir, ic, n_basis, file_unit
    integer :: K_row, K_col, m_shift_row, m_shift_col, first_m_row, last_m_row, first_m_col, last_m_col, sym_row, sym_col
    ! Arrays for data arrays length
    integer, allocatable :: num_solutions_1d_row(:), num_solutions_1d_col(:) ! number of 1D solutions in each theta slice (for fixed rho)
    integer, allocatable :: num_solutions_2d_row(:), num_solutions_2d_col(:) ! number of 2D solutions in each rho slice
    ! Arrays for 1D eigenstates
    type(array_1d_real), allocatable :: energies_1d_row(:), energies_1d_col(:)
    type(array_2d_real), allocatable :: exp_coeffs_1d_row(:), exp_coeffs_1d_col(:) ! Expansion coefficients of all 1D solutions over primitive basis set (sin/cos)
    ! Arrays for 2D eigenstates
    real(real64), allocatable :: energies_2d_row(:), energies_2d_col(:)
    real(real64), allocatable :: exp_coeffs_2d_row(:, :), exp_coeffs_2d_col(:, :) ! Expansion coefficients of all 2D solutions over 1D solutions
    real(real64), allocatable :: exp_coeffs_2d_prim_row(:), exp_coeffs_2d_prim_col(:) ! Expansion coefficients of a 2D solution over primitive basis set (sin/cos)
    ! Calculation of coriolis term
    real(real64), allocatable :: cor_block(:, :)
    real(real64), allocatable :: cor_factors_m(:), m_product(:), cor_factors_b(:), partial_sum(:)
    character(:), allocatable :: root_path, sym_folder_row, sym_folder_col, block_info_path_row, block_info_path_col, cor_block_path
    character(:), allocatable :: solutions_1d_path_row, solutions_1d_path_col, solutions_2d_path_row, solutions_2d_path_col

    call print_parallel('Calculating Coriolis term')
    my_id = get_proc_id()
    n_procs = get_num_procs()
    n_basis = params % basis % num_functions_phi
    allocate(cor_factors_m(n_basis - 1), cor_factors_b(size(theta_grid)), partial_sum(size(theta_grid)))

    root_path = params % root_path
    K_row = params % K(1)
    K_col = iff(params % basis % fixed % enabled == 1, K_row, K_row + 1)
    sym_row = params % basis % symmetry
    sym_col = 1 - sym_row ! opposite symmetry for column block

    sym_folder_row = get_sym_path_root(root_path, K_row, sym_row)
    sym_folder_col = get_sym_path_root(root_path, K_col, sym_col)
    block_info_path_row = get_block_info_path(sym_folder_row)
    block_info_path_col = get_block_info_path(sym_folder_col)
    num_solutions_2d_row = load_basis_size_2d(block_info_path_row)
    num_solutions_2d_col = load_basis_size_2d(block_info_path_col)

    ! Loop over diagonal blocks
    do n = my_id + 1, size(rho_grid), n_procs
      ! Skip if no basis
      if (num_solutions_2d_row(n) == 0 .or. num_solutions_2d_col(n) == 0) then
        cycle
      end if

      ! Load solutions
      solutions_1d_path_row = get_solutions_1d_path(sym_folder_row, n)
      call load_solutions_1D(solutions_1d_path_row, size(theta_grid), n_basis, num_solutions_1d_row, energies_1d_row, exp_coeffs_1d_row)
      solutions_2d_path_row = get_solutions_2d_path(sym_folder_row, n)
      call load_solutions_2D(solutions_2d_path_row, energies_2d_row, exp_coeffs_2d_row)

      solutions_1d_path_col = get_solutions_1d_path(sym_folder_col, n)
      call load_solutions_1D(solutions_1d_path_col, size(theta_grid), n_basis, num_solutions_1d_col, energies_1d_col, exp_coeffs_1d_col)
      solutions_2d_path_col = get_solutions_2d_path(sym_folder_col, n)
      call load_solutions_2D(solutions_2d_path_col, energies_2d_col, exp_coeffs_2d_col)

      ! Compute m-factors
      do m = 1, n_basis - 1
        ! Sign of factor depends on which function we take derivative from
        cor_factors_m(m) = (-1) ** sym_row * m
      end do

      ! Compute b-factors
      do l = 1, size(theta_grid)
        cor_factors_b(l) = get_rotational_b(mu, rho_grid(n), theta_grid(l)) * cos(theta_grid(l))
      end do

      m_shift_row = 1 - sym_row ! controls index range for "m-product"
      m_shift_col = 1 - sym_col
      allocate(cor_block(num_solutions_2d_row(n), num_solutions_2d_col(n)))
      ! Compute elements of coriolis block
      do ic = 1, num_solutions_2d_col(n)
        exp_coeffs_2d_prim_col = expansion_1D_to_primitive(num_solutions_1d_col, exp_coeffs_1d_col, exp_coeffs_2d_col(:, ic))
        do ir = 1, num_solutions_2d_row(n)
          exp_coeffs_2d_prim_row = expansion_1D_to_primitive(num_solutions_1d_row, exp_coeffs_1d_row, exp_coeffs_2d_row(:, ir))
          ! Compute partial sums
          do l = 1, size(theta_grid)
            first_m_row = (l - 1)*n_basis + 1 + m_shift_row
            last_m_row = l*n_basis - 1 + m_shift_row
            first_m_col = (l - 1)*n_basis + 1 + m_shift_col
            last_m_col = l*n_basis - 1 + m_shift_col
            m_product = exp_coeffs_2d_prim_row(first_m_row : last_m_row) * exp_coeffs_2d_prim_col(first_m_col : last_m_col)
            partial_sum(l) = dot_product(cor_factors_m, m_product)
          end do
          cor_block(ir, ic) = dot_product(cor_factors_b, partial_sum) * 4d0 ! Apply coriolis operator
        end do
      end do

      ! Save cor block
      cor_block_path = get_coriolis_overlap_file_path(sym_folder_row, n) ! Results are saved in the folder corresponding to K_row
      open(newunit=file_unit, file=cor_block_path, form='unformatted')
      write(file_unit) cor_block
      close(file_unit)
      deallocate(cor_block)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates asymmetric term overlaps.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_asym_term(params, K_row, K_col, mu, rho_grid, theta_grid)
    type(input_params), intent(in) :: params
    integer, intent(in) :: K_row, K_col ! pair of Ks to calculate asymmetric term
    real(real64), intent(in) :: mu ! reduced mass
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)
    integer :: my_id, n_procs
    integer :: n, l, ir, ic, n_basis, file_unit
    ! Arrays for data arrays length
    integer, allocatable :: num_solutions_1d_row(:), num_solutions_1d_col(:) ! number of 1D solutions in each theta slice (for fixed rho)
    integer, allocatable :: num_solutions_2d_row(:), num_solutions_2d_col(:) ! number of 2D solutions in each rho slice
    ! Arrays for 1D eigenstates
    type(array_1d_real), allocatable :: energies_1d_row(:), energies_1d_col(:)
    type(array_2d_real), allocatable :: exp_coeffs_1d_row(:), exp_coeffs_1d_col(:) ! Expansion coefficients of all 1D solutions over primitive basis set (sin/cos)
    ! Arrays for 2D eigenstates
    real(real64), allocatable :: energies_2d_row(:), energies_2d_col(:)
    real(real64), allocatable :: exp_coeffs_2d_row(:, :), exp_coeffs_2d_col(:, :) ! Expansion coefficients of all 2D solutions over 1D solutions
    real(real64), allocatable :: exp_coeffs_2d_prim_row(:), exp_coeffs_2d_prim_col(:) ! Expansion coefficients of a 2D solution over primitive basis set (sin/cos)
    ! Calculation of asym term
    real(real64), allocatable :: asym_block(:, :)
    real(real64), allocatable :: asym_factors(:), partial_sum(:)
    character(:), allocatable :: root_path, sym_folder_row, sym_folder_col, block_info_path_row, block_info_path_col, asym_block_path
    character(:), allocatable :: solutions_1d_path_row, solutions_1d_path_col, solutions_2d_path_row, solutions_2d_path_col

    call print_parallel('Calculating asym term')
    call assert(K_row + 2 == K_col .or. K_row == 1 .and. K_col == 1 .or. params % basis % fixed % enabled == 1 .and. K_row == K_col, 'Wrong combination of Ks for asym term')
    my_id = get_proc_id()
    n_procs = get_num_procs()
    n_basis = params % basis % num_functions_phi
    allocate(asym_factors(size(theta_grid)), partial_sum(size(theta_grid)))

    root_path = params % root_path
    sym_folder_row = get_sym_path_root(root_path, K_row, params % basis % symmetry)
    sym_folder_col = get_sym_path_root(root_path, K_col, params % basis % symmetry)
    block_info_path_row = get_block_info_path(sym_folder_row)
    block_info_path_col = get_block_info_path(sym_folder_col)
    num_solutions_2d_row = load_basis_size_2d(block_info_path_row)
    num_solutions_2d_col = load_basis_size_2d(block_info_path_col)

    ! Loop over diagonal blocks
    do n = my_id + 1, size(rho_grid), n_procs
      ! Skip if no basis
      if (num_solutions_2d_row(n) == 0 .or. num_solutions_2d_col(n) == 0) then
        cycle
      end if

      ! Load solutions
      solutions_1d_path_row = get_solutions_1d_path(sym_folder_row, n)
      call load_solutions_1D(solutions_1d_path_row, size(theta_grid), n_basis, num_solutions_1d_row, energies_1d_row, exp_coeffs_1d_row)
      solutions_2d_path_row = get_solutions_2d_path(sym_folder_row, n)
      call load_solutions_2D(solutions_2d_path_row, energies_2d_row, exp_coeffs_2d_row)

      solutions_1d_path_col = get_solutions_1d_path(sym_folder_col, n)
      call load_solutions_1D(solutions_1d_path_col, size(theta_grid), n_basis, num_solutions_1d_col, energies_1d_col, exp_coeffs_1d_col)
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
        exp_coeffs_2d_prim_col = expansion_1D_to_primitive(num_solutions_1d_col, exp_coeffs_1d_col, exp_coeffs_2d_col(:, ic))
        do ir = 1, num_solutions_2d_row(n)
          exp_coeffs_2d_prim_row = expansion_1D_to_primitive(num_solutions_1d_row, exp_coeffs_1d_row, exp_coeffs_2d_row(:, ir))
          ! Compute partial sums
          do l = 1, size(theta_grid)
            partial_sum(l) = dot_product(exp_coeffs_2d_prim_row((l - 1)*n_basis + 1 : l*n_basis), exp_coeffs_2d_prim_col((l - 1)*n_basis + 1 : l*n_basis))
          end do
          asym_block(ir, ic) = dot_product(asym_factors, partial_sum) / 2d0 ! Apply asym operator
        end do
      end do

      ! Save asym block
      ! Results are saved in the folder corresponding to K_row
      if (K_row == 1 .and. K_col == 1 .and. params % basis % fixed % enabled == 0) then
        asym_block_path = get_asymmetric_overlap_file_1_path(sym_folder_row, n)
      else
        asym_block_path = get_asymmetric_overlap_file_path(sym_folder_row, n)
      end if
      open(newunit=file_unit, file=asym_block_path, form='unformatted')
      write(file_unit) asym_block
      close(file_unit)
      deallocate(asym_block)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks that all required files exist.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_prerequisites(params)
    type(input_params), intent(in) :: params
    logical :: file_exists
    character(:), allocatable :: block_info_path

    block_info_path = get_block_info_path(get_sym_path_root(params % root_path, params % K(1), params % basis % symmetry))
    inquire(file = block_info_path, exist = file_exists)
    call assert(file_exists, 'Error: basis is not computed')

    ! the other symmetry is also required in this mode
    if (params % basis % fixed % enabled == 1) then
      block_info_path = get_block_info_path(get_sym_path_root(params % root_path, params % K(1), 1 - params % basis % symmetry))
      inquire(file = block_info_path, exist = file_exists)
      call assert(file_exists, 'Error: basis of both symmetries has to be computed in fixed basis mode')
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! The main subroutine that initiates calculation of both types of extra overlaps.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_overlaps_extra(params, mu, rho_grid, theta_grid)
    type(input_params), intent(in) :: params
    real(real64), intent(in) :: mu
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)

    call check_prerequisites(params)
    if (params % basis % fixed % enabled == 1) then
      call calculate_sym_term(params, mu, rho_grid, theta_grid)
    end if

    ! Further blocks are only needed if coriolis coupling is enabled
    if (params % use_rovib_coupling == 0) then
      return
    end if

    ! coriolis term
    if (params % K(1) + 1 <= params % J .or. params % basis % fixed % enabled == 1) then
      call calculate_coriolis_term(params, mu, rho_grid, theta_grid)
    end if

    ! asym term for fixed basis
    if (params % basis % fixed % enabled == 1) then
      call calculate_asym_term(params, params % K(1), params % K(1), mu, rho_grid, theta_grid)
      return
    end if

    ! Calculate extra asym for (1, 1)-block
    if (params % K(1) == 1) then
      call calculate_asym_term(params, params % K(1), params % K(1), mu, rho_grid, theta_grid)
    end if

    ! offdiag asym term
    if (params % K(1) + 2 <= params % J) then
      call calculate_asym_term(params, params % K(1), params % K(1) + 2, mu, rho_grid, theta_grid)
    end if
  end subroutine

end module
