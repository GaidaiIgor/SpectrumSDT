module state_properties_mod
  use array_1d_mod
  use array_2d_mod
  use constants
  use general_utils
  use input_params_mod
  use io_utils
  use mpi_f08
  use parallel_utils
  use path_utils
  use rovib_io_mod
  use rovib_utils_mod

  use debug_tools
  implicit none

  private
  public :: calculate_state_properties

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads 1D expansion coefficients of a specific K and sym stored at *sym_path* into *K_ind* of *num_solutions_1d* and *As*
! N, L, M - basis size in rho, theta and phi
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_1D_expansion_coefficients_K(sym_path, K_ind, N, L, M, num_solutions_1d, As)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: K_ind, N, L, M
    integer, intent(inout) :: num_solutions_1d(:, :, :)
    type(array_2d_real), intent(inout) :: As(:, :, :)
    integer :: n_val
    integer, allocatable :: num_solutions_1d_ls(:)
    character(:), allocatable :: solutions_path ! path to file with all 1D solutions for all l-s for given K and n
    type(array_1d_real), allocatable :: energies_1d(:)
    type(array_2d_real), allocatable :: As_slice(:) ! slice along l of the global array As
    
    do n_val = 1, N
      solutions_path = get_solutions_1d_path(sym_path, n_val)
      call load_solutions_1D(solutions_path, L, M, num_solutions_1d_ls, energies_1d, As_slice)
      num_solutions_1d(K_ind, n_val, :) = num_solutions_1d_ls
      As(K_ind, n_val, :) = As_slice
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 1D expansion coefficients
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_1D_expansion_coefficients(params, N, L, As, num_solutions_1d)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N, L ! Num of points along rho and theta
    type(array_2d_real), allocatable, intent(out) :: As(:, :, :) ! K x N x L. Inner matrix is M x S_Knl
    integer, allocatable, intent(out) :: num_solutions_1d(:, :, :)
    integer :: K_start, K_end, total_Ks, K, K_ind, K_sym
    character(:), allocatable :: sym_path

    K_start = params % K(1)
    K_end = params % K(2)
    total_Ks = K_end - K_start + 1
    allocate(As(total_Ks, N, L), num_solutions_1d(total_Ks, N, L))

    do K = K_start, K_end
      K_ind = get_k_ind(K, K_start)
      K_sym = get_k_symmetry(K, params % symmetry)
      sym_path = get_sym_path_root(params % root_path, K, K_sym)
      call load_1D_expansion_coefficients_K(sym_path, K_ind, N, L, params % basis_size_phi, num_solutions_1d, As)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 1D expansion coefficients for fix_basis_jk = 1
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_1D_expansion_coefficients_fixed_basis(params, N, L, As, num_solutions_1d)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N, L ! Num of points along rho and theta
    type(array_2d_real), allocatable, intent(out) :: As(:, :, :) ! first dim: 1 - even, 2 - odd. Inner matrix is M x S_Knl. 
    integer, allocatable, intent(out) :: num_solutions_1d(:, :, :)
    integer :: K_ind, K_sym
    character(:), allocatable :: sym_path

    allocate(As(2, N, L), num_solutions_1d(2, N, L))
    do K_sym = 0, 1
      K_ind = K_sym + 1
      sym_path = get_sym_path_root(params % basis_root_path, params % basis_K, K_sym)
      call load_1D_expansion_coefficients_K(sym_path, K_ind, N, L, params % basis_size_phi, num_solutions_1d, As)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 2D expansion coefficients for a given K at *sym_path* into *K_ind* of corresponding arrays
! N, L - number of points along rho and theta
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_2D_expansion_coefficients_K(sym_path, K_ind, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: K_ind, N, L
    integer, intent(in) :: num_solutions_1d(:, :, :)
    type(array_2d_real), intent(inout) :: Bs(:, :, :) 
    integer, intent(inout) :: num_solutions_2d(:, :)
    type(array_2d_real), allocatable, optional, intent(inout) :: Bs_plain(:, :) ! all Ls are stacked together
    integer :: n_val, l_val, next_start
    real*8, allocatable :: energies_2d(:)
    real*8, allocatable :: Bs_raw(:, :) ! Unarranged matrix of expansion coefficients. Single element of Bs_plain.
    character(:), allocatable :: solutions_path ! path to file with all 2D solutions for given K and n

    do n_val = 1, N
      solutions_path = get_solutions_2d_path(sym_path, n_val)
      call load_solutions_2D(solutions_path, energies_2d, Bs_raw)

      ! Not all parameter sets keep >0 solutions
      if (.not. allocated(Bs_raw)) then
        cycle
      end if

      if (present(Bs_plain)) then
        Bs_plain(K_ind, n_val) = array_2d_real(Bs_raw)
      end if

      num_solutions_2d(K_ind, n_val) = size(energies_2d)
      next_start = 1
      do l_val = 1, L
        Bs(K_ind, n_val, l_val) = array_2d_real(Bs_raw(next_start : next_start + num_solutions_1d(K_ind, n_val, l_val) - 1, :))
        next_start = next_start + num_solutions_1d(K_ind, n_val, l_val)
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 2D expansion coefficients
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_2D_expansion_coefficients(params, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N, L ! Num of points along rho and theta
    integer, intent(in) :: num_solutions_1d(:, :, :)
    type(array_2d_real), allocatable, intent(out) :: Bs(:, :, :) ! K x N x L
    type(array_2d_real), allocatable, optional, intent(out) :: Bs_plain(:, :) ! K x N, all Ls are stacked in columns
    integer, allocatable, intent(out) :: num_solutions_2d(:, :)
    integer :: K_start, K_end, total_Ks
    integer :: K, K_ind, K_sym
    character(:), allocatable :: sym_path

    K_start = params % K(1)
    K_end = params % K(2)
    total_Ks = K_end - K_start + 1
    allocate(Bs(total_Ks, N, L), num_solutions_2d(total_Ks, N))
    if (present(Bs_plain)) then
      allocate(Bs_plain(total_Ks, N))
    end if

    do K = K_start, K_end
      K_ind = get_k_ind(K, K_start)
      K_sym = get_k_symmetry(K, params % symmetry)
      sym_path = get_sym_path_root(params % root_path, K, K_sym)
      call load_2D_expansion_coefficients_K(sym_path, K_ind, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 2D expansion coefficients when fix_basis_jk = 1
! Bs - Inner dimensions: S_Knl x S_Kn. All coeffs of a single 2D solution should be collected over different Ls
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_2D_expansion_coefficients_fixed_basis(params, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N, L ! Num of points along rho and theta
    integer, intent(in) :: num_solutions_1d(:, :, :)
    type(array_2d_real), allocatable, intent(out) :: Bs(:, :, :)
    integer, allocatable, intent(out) :: num_solutions_2d(:, :)
    type(array_2d_real), allocatable, optional, intent(out) :: Bs_plain(:, :) ! K x N, all Ls are stacked together
    integer :: K_ind, K_sym
    character(:), allocatable :: sym_path

    allocate(Bs(2, N, L), num_solutions_2d(2, N))
    if (present(Bs_plain)) then
      allocate(Bs_plain(2, N))
    end if

    do K_sym = 0, 1
      K_ind = K_sym + 1
      sym_path = get_sym_path_root(params % basis_root_path, params % basis_K, K_sym)
      call load_2D_expansion_coefficients_K(sym_path, K_ind, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 3D expansion coefficients
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_3D_expansion_coefficients(params, N, num_solutions_2d, Cs, Cs_plain)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N ! Num of points along rho
    integer, intent(in) :: num_solutions_2d(:, :)
    type(array_2d_complex), allocatable, intent(out) :: Cs(:, :) ! K x N. Inner dimensions: S_Kn x S
    complex*16, allocatable, optional, intent(out) :: Cs_plain(:, :) ! Stacked over Ks and ns (in this order) version of Cs
    integer :: total_solutions_2d, first_elem, proc_elems, sln_ind, total_Ks, start_ind, K, K_ind, K_sym, K_ind_comp, n_val
    complex*16, allocatable :: Cs_raw_col(:) ! column in Cs_raw
    complex*16, allocatable :: Cs_raw(:, :) ! unarranged Cs
    character(:), allocatable :: sym_path, exp_coeffs_path

    ! Count total number of 2D solutions
    total_solutions_2d = 0
    do K = params % K(1), params % K(2)
      call get_k_attributes(K, params, K_ind, K_sym, K_ind_comp)
      total_solutions_2d = total_solutions_2d + sum(num_solutions_2d(K_ind_comp, :))
    end do

    ! Load raw matrix
    call get_proc_elem_range(params % num_states, first_elem, proc_elems)
    allocate(Cs_raw(total_solutions_2d, proc_elems))
    sym_path = get_sym_path_params(params)
    do sln_ind = first_elem, first_elem + proc_elems - 1
      exp_coeffs_path = get_solution_3d_path(sym_path, sln_ind)
      call load_solution_3D(exp_coeffs_path, total_solutions_2d, Cs_raw_col)
      Cs_raw(:, sln_ind - first_elem + 1) = Cs_raw_col
    end do

    if (present(Cs_plain)) then
      Cs_plain = Cs_raw
    end if

    ! Rearrange matrix
    total_Ks = params % K(2) - params % K(1) + 1
    allocate(Cs(total_Ks, N))
    start_ind = 1
    do K = params % K(1), params % K(2)
      call get_k_attributes(K, params, K_ind, K_sym, K_ind_comp)
      do n_val = 1, N
        Cs(K_ind, n_val) = array_2d_complex(Cs_raw(start_ind : start_ind + num_solutions_2d(K_ind_comp, n_val) - 1, :))
        start_ind = start_ind + num_solutions_2d(K_ind_comp, n_val)
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Evaluates the value of phi integral in the range [a, b] (F~_Ks(m,m') = int(F_m1^+-(phi) * F_m2^+-(phi), a, b))
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_phi_integral(m1, m2, a, b, F_sym) result(integral_value)
    integer, intent(in) :: m1, m2 ! Frequency of integrand functions
    real*8, intent(in) :: a, b ! Range of integration
    integer, intent(in) :: F_sym ! Symmetry of integrand functions: 0 or 1 (same for both)
    real*8 :: integral_value ! Result
    integer :: sign
    real*8 :: norm

    call assert(m1 == 0 .or. m2 == 0 .means. F_sym == 0, 'm can be 0 only for symmetry 0')
    call assert(m1 >= 0 .and. m2 >= 0, 'm values should be non-negative')

    sign = (-1) ** F_sym
    if (m1 == 0 .and. m2 == 0) then
      integral_value = (b - a) / (2*pi)
    elseif (m1 == m2) then
      integral_value = 1/(2*pi) * (b - a - sign*(sin(2*a*m1) - sin(2*b*m1)) / (2*m1))
    elseif (m1 /= m2) then
      norm = 1/(2*pi*sqrt(1d0*(delta(m1, 0) + 1)*(delta(m2, 0) + 1)))
      integral_value = norm * ((sin(b*(m1 - m2)) - sin(a*(m1 - m2))) / (m1 - m2) + (sign*sin(b*(m1 + m2)) - sign*sin(a*(m1 + m2))) / (m1 + m2))
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates m-sum in wf integral
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_m_sum(params, As, K_val, n, l, phi_range, j_sums) result(m_sum)
    class(input_params), intent(in) :: params
    class(array_2d_real), intent(in) :: As(:, :, :)
    integer, intent(in) :: K_val, n, l
    real*8, intent(in) :: phi_range(2)
    complex*16, intent(in) :: j_sums(:)
    real*8 :: m_sum
    integer :: K_ind, K_sym, K_ind_comp, m1, m2, m1_ind, m2_ind
    real*8 :: product_coeff, phi_integral

    call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_comp)
    m_sum = 0
    if (compare_reals(phi_range(1), 0d0) == 0 .and. compare_reals(phi_range(2), pi) == 0) then
      ! Simplified method for the case when phi range is not restricted
      do m1_ind = 1, size(As(K_ind_comp, n, l) % p, 1)
        m_sum = m_sum + real(j_sums(m1_ind) * conjg(j_sums(m1_ind)))
      end do
    else
      do m1_ind = 1, size(As(K_ind_comp, n, l) % p, 1)
        do m2_ind = m1_ind, size(As(K_ind_comp, n, l) % p, 1)
          m1 = m_ind_to_m(m1_ind, K_val, params % symmetry)
          m2 = m_ind_to_m(m2_ind, K_val, params % symmetry)
          phi_integral = calc_phi_integral(m1, m2, phi_range(1), phi_range(2), K_sym)
          product_coeff = real(j_sums(m1_ind) * conjg(j_sums(m2_ind)))
          if (m1_ind /= m2_ind) then
            product_coeff = 2 * product_coeff ! since m-integral matrix is hermitian
          end if
          m_sum = m_sum + phi_integral * product_coeff
        end do
      end do
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates <Psi|O|Psi>-like expressions in the specified region. 
! Operators O need to be supplied in form of multiplicative factors, which affect the overall sum at different stages.
! Default operator is identity i.e. the function simply evaluates total probability of proc's k-th state in the specified region of the PES.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, K_range, n_range, phi_range, n_factors) result(res)
    class(input_params), intent(in) :: params
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A (1D), S_Knl x S_Kn for B (2D)
    class(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    class(array_2d_complex), intent(in) :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    integer, intent(in) :: proc_k ! local 3D state index
    integer, optional, intent(in) :: K_range(2) ! values of K that should be considered in the sum
    integer, optional, intent(in) :: n_range(2) ! indexes of the border rho points on the grid included in the corresponding range
    real*8, optional, intent(in) :: phi_range(2) ! limits of phi integral
    real*8, optional, intent(in) :: n_factors(:) ! multiplying factors in n-sum
    real*8 :: res
    integer :: K_val, K_ind, K_ind_comp, K_sym, n, l, m_ind, j
    integer :: K_range_act(2), n_range_act(2)
    real*8 :: i_sum, l_sum
    real*8 :: phi_range_act(2)
    real*8, allocatable :: n_factors_act(:), n_factors_default(:)
    complex*16, allocatable :: j_sums(:) ! j-sums for different ms

    K_range_act = arg_or_default(K_range, [params % K(1), params % K(2)])
    n_range_act = arg_or_default(n_range, [1, size(As, 2)])
    phi_range_act = arg_or_default(phi_range, [0d0, pi])

    allocate(n_factors_default(n_range_act(2) - n_range_act(1) + 1))
    n_factors_default = 1
    n_factors_act = arg_or_default(n_factors, n_factors_default)

    call assert(size(n_factors_act) == size(n_factors_default), 'Error: size of n_factors does not match n_range')
    call assert(phi_range_act(1) >= 0 .and. phi_range_act(2) <= pi, 'Error: phi_range should be [0, pi]')

    allocate(j_sums(params % basis_size_phi))
    res = 0
    do K_val = K_range_act(1), K_range_act(2)
      call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_comp)
      do n = n_range_act(1), n_range_act(2)
        l_sum = 0
        do l = 1, size(As, 3)
          j_sums = 0
          do m_ind = 1, size(As(K_ind_comp, n, l) % p, 1)
            do j = 1, size(Bs(K_ind_comp, n, l) % p, 2)
              i_sum = sum(Bs(K_ind_comp, n, l) % p(:, j) * As(K_ind_comp, n, l) % p(m_ind, :))
              j_sums(m_ind) = j_sums(m_ind) + proc_Cs(K_ind, n) % p(j, proc_k) * i_sum
            end do
          end do
          l_sum = l_sum + calc_m_sum(params, As, K_val, n, l, phi_range_act, j_sums)
        end do
        res = res + l_sum * n_factors_act(n - n_range_act(1) + 1)
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates total probability of k-th state in symmetric molecule part of the PES (2*pi/3 <= phi <= 4*pi/3, regardless of rho)
! Asymmetric can be obtained as 1 - symmetric
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_symmetric_molecule_probability(params, As, Bs, proc_Cs, proc_k) result(sym_prob)
    class(input_params), intent(in) :: params
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A (1D), S_Knl x S_Kn for B (2D)
    type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    type(array_2d_complex), intent(in) :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    integer, intent(in) :: proc_k ! local 3D state index
    real*8 :: sym_prob
    sym_prob = 2 * calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, phi_range = [2*pi/3, pi])
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates total probability of k-th state in different parts of the PES.
! The output array contains probabilities in the following regions:
! 1 - rho in n_range_cov, phi in [2pi/3, 4pi/3]
! 2 - rho in n_range_cov, phi in [0, 2pi/3] + [4pi/3, 2pi]
! 3 - rho in n_range_vdw, phi in [2pi/3, 4pi/3]
! 4 - rho in n_range_vdw, phi in [pi/3, 2pi/3] + [4pi/3, 5pi/3]
! 5 - rho in n_range_vdw, phi in [0, pi/3] + [5pi/3, 2pi]
!-------------------------------------------------------------------------------------------------------------------------------------------
  function assign_molecule_extra(params, As, Bs, proc_Cs, n_range_cov, n_range_vdw, proc_k) result(region_probs)
    class(input_params), intent(in) :: params
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A (1D), S_Knl x S_Kn for B (2D)
    type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    type(array_2d_complex), intent(in) :: proc_Cs(:, :) ! 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    integer, intent(in) :: n_range_cov(2), n_range_vdw(2) ! indexes of the border rho points on the grid included in the corresponding range
    integer, intent(in) :: proc_k ! local 3D state index
    real*8 :: region_probs(5)

    region_probs(1) = 2 * calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, n_range = n_range_cov, phi_range = [2*pi/3, pi])
    region_probs(2) = 2 * calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, n_range = n_range_cov, phi_range = [0d0, 2*pi/3])
    region_probs(3) = 2 * calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, n_range = n_range_vdw, phi_range = [2*pi/3, pi])
    region_probs(4) = 2 * calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, n_range = n_range_vdw, phi_range = [pi/3, 2*pi/3])
    region_probs(5) = 2 * calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, n_range = n_range_vdw, phi_range = [0d0, pi/3])
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates total probability of each state in symmetric molecule part of the PES (2*pi/3 <= phi <= 4*pi/3, regardless of rho)
! Asymmetric can be obtained as 1 - symmetric
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_symmetric_molecule_probability_all(params, As, Bs, proc_Cs) result(sym_probs)
    class(input_params), intent(in) :: params
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A (1D), S_Knl x S_Kn for B (2D)
    type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    type(array_2d_complex), intent(in) :: proc_Cs(:, :) ! 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    real*8, allocatable :: sym_probs(:)
    integer :: proc_k
    ! Parallel
    integer :: proc_id, n_procs, proc_first_state, proc_states, ierr
    integer, allocatable :: recv_counts(:), recv_shifts(:)
    real*8, allocatable :: sym_probs_chunk(:)

    call get_proc_info(proc_id, n_procs)
    call get_proc_elem_range(params % num_states, proc_first_state, proc_states, recv_counts, recv_shifts)

    ! Calculate chunks
    allocate(sym_probs_chunk(proc_states))
    do proc_k = 1, proc_states
      sym_probs_chunk(proc_k) = calculate_symmetric_molecule_probability(params, As, Bs, proc_Cs, proc_k)
    end do

    ! Gather results
    if (n_procs == 1) then
      sym_probs = sym_probs_chunk
    else
      allocate(sym_probs(params % num_states))
      call MPI_Gatherv(sym_probs_chunk, proc_states, MPI_DOUBLE_PRECISION, sym_probs, recv_counts, recv_shifts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates K probability distribution for a given state
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_K_dist(params, As, Bs, proc_Cs, proc_k) result(K_dist)
    class(input_params), intent(in) :: params
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B
    type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    type(array_2d_complex), intent(in) :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    integer, intent(in) :: proc_k ! local 3D state index
    real*8, allocatable :: K_dist(:) ! Distibution over all Ks, always includes K=0
    integer :: K_val, K_sym, K_ind, K_ind_comp

    allocate(K_dist(params % J + 1))
    do K_val = params % K(1), params % K(2)
      call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_comp)
      K_dist(K_ind + params % K(1)) = calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, K_range = [K_val, K_val])
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates K probability distribution for all states
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_K_dist_all(params, As, Bs, proc_Cs) result(K_dists)
    class(input_params), intent(in) :: params
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B
    type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    type(array_2d_complex), intent(in) :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    real*8, allocatable :: K_dists(:, :)
    integer :: proc_id, n_procs, proc_first_state, proc_states, proc_k, ind, ierr
    integer, allocatable :: recv_counts(:), recv_shifts(:)
    real*8, allocatable :: K_dists_chunk(:, :)

    call get_proc_info(proc_id, n_procs)
    call get_proc_elem_range(params % num_states, proc_first_state, proc_states, recv_counts, recv_shifts)

    ! Calculate chunks
    allocate(K_dists_chunk(proc_states, params % J + 1))
    do proc_k = 1, proc_states
      K_dists_chunk(proc_k, :) = calculate_K_dist(params, As, Bs, proc_Cs, proc_k)
    end do

    ! Gather results
    if (n_procs == 1) then
      K_dists = K_dists_chunk
      return
    end if
    
    allocate(K_dists(params % num_states, params % J + 1))
    do ind = 1, size(K_dists, 2)
      call MPI_Gatherv(K_dists_chunk(1, ind), proc_states, MPI_DOUBLE_PRECISION, K_dists(1, ind), recv_counts, recv_shifts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates K probability distribution for a given state
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_gamma_channels(params, As, Bs, proc_Cs, cap, proc_k) result(gammas)
    class(input_params), intent(in) :: params
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B
    type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    type(array_2d_complex), intent(in) :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    real*8, intent(in) :: cap(:)
    integer, intent(in) :: proc_k ! local 3D state index
    real*8 :: gammas(2) ! 1 - channel A, 2 - channel B

    ! One factor of 2 is due to phi symmetry, another is due to gamma definition
    gammas(1) = 2 * 2 * calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, phi_range = [pi/3, pi], n_factors = cap)
    gammas(2) = 2 * 2 * calculate_wf_integral_region(params, As, Bs, proc_Cs, proc_k, phi_range = [0d0, pi/3], n_factors = cap)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates channel-specific gammas for all states
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_gamma_channels_all(params, As, Bs, proc_Cs, cap) result(gammas)
    class(input_params), intent(in) :: params
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B
    type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    type(array_2d_complex), intent(in) :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    real*8, intent(in) :: cap(:)
    real*8, allocatable :: gammas(:, :)
    integer :: proc_id, n_procs, proc_first_state, proc_states, proc_k, ind, ierr
    integer, allocatable :: recv_counts(:), recv_shifts(:)
    real*8, allocatable :: gammas_chunk(:, :)

    call get_proc_info(proc_id, n_procs)
    call get_proc_elem_range(params % num_states, proc_first_state, proc_states, recv_counts, recv_shifts)

    ! Calculate chunks
    allocate(gammas_chunk(proc_states, 2))
    do proc_k = 1, proc_states
      gammas_chunk(proc_k, :) = calculate_gamma_channels(params, As, Bs, proc_Cs, cap, proc_k)
    end do

    ! Gather results
    if (n_procs == 1) then
      gammas = gammas_chunk
      return
    end if
    
    allocate(gammas(params % num_states, 2))
    do ind = 1, size(gammas, 2)
      call MPI_Gatherv(gammas_chunk(1, ind), proc_states, MPI_DOUBLE_PRECISION, gammas(1, ind), recv_counts, recv_shifts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Writes results to file
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_results(params, energies, sym_probs, K_dists, gammas)
    class(input_params), intent(in) :: params
    real*8, intent(in) :: energies(:), sym_probs(:)
    real*8, intent(in) :: K_dists(:, :), gammas(:, :)
    integer :: file_unit, k, K_ind, col_width
    character(:), allocatable :: sym_path, properties_result_path
    ! Parallel
    integer :: proc_id, n_procs

    call get_proc_info(proc_id, n_procs)
    if (proc_id /= 0) then
      return
    end if

    call assert(size(sym_probs) == size(K_dists, 1), 'Sizes of input arrays have to be the same')
    sym_path = get_sym_path_params(params)
    properties_result_path = get_properties_result_path(sym_path)

    open(newunit = file_unit, file = properties_result_path)
    col_width = 25
    write(file_unit, '(2A' // num2str(col_width) // ')', advance = 'no') align_center('Energy (cm-1)', col_width), align_center('Symmetric probability', col_width)
    do K_ind = 1, params % J + 1
      write(file_unit, '(A' // num2str(col_width) // ')', advance = 'no') align_center('K' // num2str(K_ind - 1) // ' probability', col_width)
    end do
    write(file_unit, '(2A' // num2str(col_width) // ')', advance = 'no') align_center('Gamma A (cm-1)', col_width), align_center('Gamma B (cm-1)', col_width)
    write(file_unit, *)

    do k = 1, size(sym_probs)
      write(file_unit, '(2G' // num2str(col_width) // '.16)', advance = 'no') energies(k), sym_probs(k)
      write(file_unit, '(' // num2str(params % J + 1) // 'G' // num2str(col_width) // '.16)', advance = 'no') K_dists(k, :)
      ! do K_ind = 1, params % J + 1
      !   write(file_unit, '(G' // num2str(col_width) // '.16)', advance = 'no') K_dists(k, K_ind)
      ! end do
      write(file_unit, '(2G' // num2str(col_width) // '.16)', advance = 'no') gammas(k, :) * autown
      write(file_unit, *)
    end do
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks orthonormality of 1D solutions
!-------------------------------------------------------------------------------------------------------------------------------------------
  function check_orthonormality_1D(As) result(max_deviation)
    type(array_2d_real), intent(in) :: As(:, :, :)
    real*8 :: max_deviation ! from orthonormality
    integer :: K_ind, n, l, i1, i2
    real*8 :: overlap
    real*8, allocatable :: solution_set(:, :)

    max_deviation = 0
    do K_ind = 1, size(As, 1)
      do n = 1, size(As, 2)
        do l = 1, size(As, 3)
          solution_set = As(K_ind, n, l) % p
          do i1 = 1, size(solution_set, 2)
            do i2 = i1, size(solution_set, 2)           
              overlap = sum(solution_set(:, i1) * solution_set(:, i2))
              max_deviation = max(max_deviation, abs(0 ** (i2 - i1) - overlap))
            end do
          end do
        end do
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks orthonormality of 2D solutions
!-------------------------------------------------------------------------------------------------------------------------------------------
  function check_orthonormality_2D(As, Bs) result(max_deviation)
    type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :)
    real*8 :: max_deviation ! from orthonormality
    integer :: K_ind, n, l, m_ind, j1, j2
    real*8 :: i1_sum, i2_sum, overlap

    max_deviation = 0
    do K_ind = 1, size(As, 1)
      do n = 1, size(As, 2)
        print *, 'Doing n: ', n
        do j1 = 1, size(Bs(K_ind, n, 1) % p, 2)
          do j2 = j1, size(Bs(K_ind, n, 1) % p, 2)

            overlap = 0
            do l = 1, size(As, 3)
              do m_ind = 1, size(As(K_ind, n, l) % p, 1)
                i1_sum = sum(Bs(K_ind, n, l) % p(:, j1) * As(K_ind, n, l) % p(m_ind, :))
                i2_sum = sum(Bs(K_ind, n, l) % p(:, j2) * As(K_ind, n, l) % p(m_ind, :))
                overlap = overlap + i1_sum * i2_sum
              end do
            end do
            max_deviation = max(max_deviation, abs(0 ** (j2 - j1) - overlap))
          end do
        end do
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks orthonormality of 3D solutions. Returns deviation from expected value (0 for different ks, 1 for same ks).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function check_orthonormality_3D(As, Bs, Cs, params, k1_first, k2_first, k1_last, k2_last) result(max_deviation)
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B.
    type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    type(array_2d_complex), intent(in) :: Cs(:, :) ! 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S.
    class(input_params), intent(in) :: params
    integer, intent(in) :: k1_first, k2_first
    integer, optional, intent(in) :: k1_last, k2_last
    real*8 :: max_deviation
    integer :: K_val, K_sym, K_ind, K_ind_comp, n, l, m_ind, j, k1, k2, k1_last_act, k2_last_act
    real*8 :: i_sum
    complex*16 :: j1_sum, j2_sum
    complex*16 :: overlap

    k1_last_act = arg_or_default(k1_last, k1_first)
    k2_last_act = arg_or_default(k2_last, k2_first)
    max_deviation = 0
    do k1 = k1_first, k1_last_act ! size(Cs(1, 1) % p, 2)
      do k2 = k2_first, k2_last_act ! k1, size(Cs(1, 1) % p, 2)
        if (k2 < k1) then
          print *, 'Skipping ovelap between', k1, k2, 'since k1 has to be >= than k2'
          cycle
        end if

        print *, 'Calculating overlap between', k1, k2
        overlap = 0
        do K_val = params % K(1), params % K(2)
          call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_comp)
          do n = 1, size(As, 2)
            ! call track_progress(n * 1d0 / size(As, 2))
            do l = 1, size(As, 3)
              do m_ind = 1, size(As(K_ind_comp, n, l) % p, 1)
                j1_sum = 0
                j2_sum = 0
                do j = 1, size(Bs(K_ind_comp, n, l) % p, 2)
                  i_sum = sum(Bs(K_ind_comp, n, l) % p(:, j) * As(K_ind_comp, n, l) % p(m_ind, :))
                  j1_sum = j1_sum + conjg(Cs(K_ind, n) % p(j, k1)) * i_sum
                  j2_sum = j2_sum + Cs(K_ind, n) % p(j, k2) * i_sum
                end do
                overlap = overlap + j1_sum * j2_sum
              end do
            end do
          end do
        end do

        max_deviation = max(max_deviation, abs(0 ** (k2 - k1) - overlap))
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! The main procedure for calculation of states properties
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_state_properties(params, N, L, cap)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N, L ! Num of points along rho and theta
    real*8, intent(in) :: cap(:)
    integer, allocatable :: num_solutions_2d(:, :)
    integer, allocatable :: num_solutions_1d(:, :, :)
    real*8, allocatable :: energies_3d(:), sym_probs(:)
    real*8, allocatable :: K_dists(:, :) ! n_state x (J + 1)
    real*8, allocatable :: gammas(:, :) ! n_state x 2
    character(:), allocatable :: sym_path, spectrum_path
    ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B
    type(array_2d_real), allocatable :: As(:, :, :), Bs(:, :, :)
    type(array_2d_complex), allocatable :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S

    ! Load expansion coefficients
    if (params % fix_basis_jk == 1) then
      call load_1D_expansion_coefficients_fixed_basis(params, N, L, As, num_solutions_1d)
      call load_2D_expansion_coefficients_fixed_basis(params, N, L, num_solutions_1d, Bs, num_solutions_2d)
    else
      call load_1D_expansion_coefficients(params, N, L, As, num_solutions_1d)
      call load_2D_expansion_coefficients(params, N, L, num_solutions_1d, Bs, num_solutions_2d)
    end if
    call print_parallel('1D and 2D expansion coefficients are loaded. Loading 3D.')
    call load_3D_expansion_coefficients(params, N, num_solutions_2d, proc_Cs)
    call print_parallel('Done loading expansion coefficients')

    ! Calculate properties
    sym_probs = calculate_symmetric_molecule_probability_all(params, As, Bs, proc_Cs)
    call print_parallel('Done calculating symmetric isotopomer probability')
    K_dists = calculate_K_dist_all(params, As, Bs, proc_Cs)
    call print_parallel('Done calculating K-distributions')
    gammas = calculate_gamma_channels_all(params, As, Bs, proc_Cs, cap)
    call print_parallel('Done calculating channel-specific gammas')

    ! Load energies
    sym_path = get_sym_path_params(params)
    spectrum_path = get_spectrum_path(sym_path)
    energies_3d = load_energies_3D(spectrum_path)

    call print_parallel('Writing results...')
    call write_results(params, energies_3d, sym_probs, K_dists, gammas)
  end subroutine
end module
