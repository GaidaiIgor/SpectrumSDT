module state_properties_mod
  use algorithms
  use array_1d_mod
  use array_2d_mod
  use constants
  use general_utils
  use input_params_mod
  use io_utils
  use mpi
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
! Selects the values of *J_low* and *J_high* used to interpolate a given *J*
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine select_interpolating_Js(J, Js_interp)
    integer, intent(in) :: J
    integer, intent(out) :: Js_interp(2)

    if (J < 4) then
      Js_interp(1) = 4
      Js_interp(2) = 8
    else if (J > 56) then
      Js_interp(1) = 52
      Js_interp(2) = 56
    else
      Js_interp(1) = J / 4 * 4 ! integer division
      Js_interp(2) = Js_interp(1) + 4
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Finds n-indices corresponding to the border (inclusive) of well region in channels B, A and S as well as common border of vdw region
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine find_barriers_n_indices(params, rho_grid, vdw_max, cov_n_inds, vdw_n_ind)
    class(input_params), intent(in) :: params
    real*8, intent(in) :: rho_grid(:)
    real*8, intent(in) :: vdw_max
    integer, intent(out) :: cov_n_inds(3)
    integer, intent(out) :: vdw_n_ind
    integer :: i
    integer :: Js_interp(2)
    real*8 :: positions_interp(3, 2), positions_estim(3) ! Barrier positions corresponding to low and high Js and estimation for current J
    character(:), allocatable :: channels_file_path

    call select_interpolating_Js(params % J, Js_interp)
    ! Load interpolating Js barrier positions
    do i = 1, 2
      channels_file_path = get_channels_file_path(params % channels_root, Js_interp(i), 0, params % symmetry)
      call load_lowest_barrier_info(channels_file_path, positions_interp(:, i))
    end do
    ! Interpolate for the current J
    do i = 1, 3
      call linear_estimation(Js_interp * 1d0, positions_interp(i, :), params % J * 1d0, positions_estim(i))
      cov_n_inds(i) = minloc(abs(rho_grid - positions_estim(i)), 1)
    end do
    vdw_n_ind = minloc(abs(rho_grid - vdw_max), 1)
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
    if (compare_reals(phi_range(1), 0d0) == 0 .and. compare_reals(phi_range(2), 2*pi) == 0) then
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
! Calculates <Psi|Psi> expressions in all regions, for all proc states
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_wf_integral_all_regions(params, As, Bs, proc_Cs) result(p_dist)
    class(input_params), intent(in) :: params
    ! Arrays of size 2 x N x L (or K x N x L if compression is disabled). Each element is a matrix with expansion coefficients - M x S_Knl for A (1D), S_Knl x S_Kn for B (2D)
    class(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    class(array_2d_complex), intent(in) :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    real*8, allocatable :: p_dist(:, :, :, :) ! S x K x N x 3
    integer :: proc_first_state, proc_states, proc_state_ind, K_val, K_ind, K_ind_comp, K_sym, n, l, m_ind, j, phi_group_ind
    real*8 :: i_sum
    real*8 :: phi_range(2), phi_borders(4)
    complex*16 :: j_sums(params % basis_size_phi) ! j-sums for different ms

    if (params % molecule == '686') then
      phi_borders = [0d0, 60d0, 117.65d0, 180d0] / 180 * pi
    else if (params % molecule == '868') then
      phi_borders = [0d0, 60d0, 122.35d0, 180d0] / 180 * pi
    else
      stop 'Error: Unknown molecule'
    end if

    call get_proc_elem_range(params % num_states, proc_first_state, proc_states)
    allocate(p_dist(proc_states, params % K(2) - params % K(1) + 1, size(As, 2), 3)) ! S x K x N x 3
    p_dist = 0
    do proc_state_ind = 1, size(p_dist, 1)
      do K_val = params % K(1), params % K(2)
        call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_comp)
        do n = 1, size(p_dist, 3)
          ! Skip empty n-slices
          if (.not. allocated(Bs(K_ind_comp, n, 1) % p)) then
            cycle
          end if

          do phi_group_ind = 1, size(p_dist, 4)
            phi_range = phi_borders(phi_group_ind : phi_group_ind+1)
            do l = 1, size(As, 3)
              j_sums = 0
              do m_ind = 1, size(As(K_ind_comp, n, l) % p, 1)
                do j = 1, size(Bs(K_ind_comp, n, l) % p, 2)
                  i_sum = sum(Bs(K_ind_comp, n, l) % p(:, j) * As(K_ind_comp, n, l) % p(m_ind, :))
                  j_sums(m_ind) = j_sums(m_ind) + proc_Cs(K_ind, n) % p(j, proc_state_ind) * i_sum
                end do
              end do
              p_dist(proc_state_ind, K_ind, n, phi_group_ind) = p_dist(proc_state_ind, K_ind, n, phi_group_ind) + calc_m_sum(params, As, K_val, n, l, phi_range, j_sums)
            end do
          end do
        end do
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates channel-specific values of gamma using probability distributions of a given state
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_gamma_channels(p_dist_k, cap) result(gammas)
    real*8, intent(in) :: p_dist_k(:, :, :) ! K x n x 3
    real*8, intent(in) :: cap(:)
    real*8 :: gammas(2) ! 1 - channel A, 2 - channel B
    integer :: K_ind

    call assert(size(cap) == size(p_dist_k, 2), 'Size of cap array does not match p_dist')
    gammas = 0
    do K_ind = 1, size(p_dist_k, 1)
      gammas(1) = gammas(1) + sum(p_dist_k(K_ind, :, 2) * cap)
      gammas(1) = gammas(1) + sum(p_dist_k(K_ind, :, 3) * cap)
      gammas(2) = gammas(2) + sum(p_dist_k(K_ind, :, 1) * cap)
    end do
    ! One factor of 2 is due to phi symmetry, another is due to gamma definition
    gammas = gammas * 2 * 2
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates channel-specific gammas for all states
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_gamma_channels_all(params, p_dist, cap) result(gammas)
    class(input_params), intent(in) :: params
    real*8, intent(in) :: p_dist(:, :, :, :)
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
      gammas_chunk(proc_k, :) = calculate_gamma_channels(p_dist(proc_k, :, :, :), cap)
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
! Calculates total probability of k-th state in different parts of the PES.
! *region_probs* contains probabilities in the following regions:
! 1 - symmetric molecule, covalent well
! 2 - asymmetric molecule, covalent well
! 3 - vdw A (sym half)
! 4 - vdw A (asym half)
! 5 - vdw B
! 6 - infinity
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_pes_region_probabilities(p_dist_k, cov_n_inds, vdw_n_ind) result(region_probs)
    real*8, intent(in) :: p_dist_k(:, :, :) ! K x n x 3
    integer, intent(in) :: cov_n_inds(3) ! indexes of the rho points that mark the border of well region in B, A and S regions respectively
    integer, intent(in) :: vdw_n_ind ! rho-index that marks the border of vdw area (common for all regions)
    real*8 :: region_probs(6)

    region_probs(1) = sum(p_dist_k(:, 1:cov_n_inds(3), 3))
    region_probs(2) = sum(p_dist_k(:, 1:cov_n_inds(2), 2)) + sum(p_dist_k(:, 1:cov_n_inds(1), 1))
    region_probs(3) = sum(p_dist_k(:, cov_n_inds(3)+1:vdw_n_ind, 3))
    region_probs(4) = sum(p_dist_k(:, cov_n_inds(2)+1:vdw_n_ind, 2))
    region_probs(5) = sum(p_dist_k(:, cov_n_inds(1)+1:vdw_n_ind, 1))
    region_probs(6) = sum(p_dist_k(:, vdw_n_ind+1:size(p_dist_k, 2), :)) ! Calculated explicitly to make sure probabilities add up properly
    region_probs = region_probs * 2 ! Due to phi symmetry
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calls calculate_pes_region_probabilities for all states and gathers results
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_pes_region_probabilities_all(params, p_dist, cov_n_inds, vdw_n_ind) result(region_probs)
    class(input_params), intent(in) :: params
    real*8, intent(in) :: p_dist(:, :, :, :) ! state x K x n x 3
    integer, intent(in) :: cov_n_inds(3)
    integer, intent(in) :: vdw_n_ind
    real*8, allocatable :: region_probs(:, :)
    integer :: proc_id, n_procs, proc_first_state, proc_states, proc_k, ind, ierr
    integer, allocatable :: recv_counts(:), recv_shifts(:)
    real*8, allocatable :: region_probs_chunk(:, :)

    call get_proc_info(proc_id, n_procs)
    call get_proc_elem_range(params % num_states, proc_first_state, proc_states, recv_counts, recv_shifts)

    ! Calculate chunks
    allocate(region_probs_chunk(proc_states, 6))
    do proc_k = 1, proc_states
      region_probs_chunk(proc_k, :) = calculate_pes_region_probabilities(p_dist(proc_k, :, :, :), cov_n_inds, vdw_n_ind)
    end do

    ! Gather results
    if (n_procs == 1) then
      region_probs = region_probs_chunk
      return
    end if
    
    allocate(region_probs(params % num_states, size(region_probs_chunk, 2)))
    do ind = 1, size(region_probs, 2)
      call MPI_Gatherv(region_probs_chunk(1, ind), proc_states, MPI_DOUBLE_PRECISION, region_probs(1, ind), recv_counts, recv_shifts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates K probability distribution for a given state
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_K_dist(params, p_dist_k) result(K_dist)
    class(input_params), intent(in) :: params
    real*8, intent(in) :: p_dist_k(:, :, :) ! K x n x 3
    real*8 :: K_dist(params % J + 1) ! Distibution over all Ks, always includes all Ks
    integer :: K_val, K_sym, K_ind, K_ind_comp

    K_dist = 0
    if (size(p_dist_k, 1) == 1) then
      K_dist(params % K(1) + 1) = 1
      return
    end if

    do K_val = params % K(1), params % K(2)
      call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_comp)
      K_dist(K_val + 1) = sum(p_dist_k(K_ind, :, :))
    end do
    K_dist = K_dist * 2 ! due to phi-symmetry
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates K probability distribution for all states
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_K_dist_all(params, p_dist) result(K_dists)
    class(input_params), intent(in) :: params
    real*8, intent(in) :: p_dist(:, :, :, :) ! states x K x n x 3
    real*8, allocatable :: K_dists(:, :)
    integer :: proc_id, n_procs, proc_first_state, proc_states, proc_k, ind, ierr
    integer, allocatable :: recv_counts(:), recv_shifts(:)
    real*8 :: K_dist_error
    real*8, allocatable :: K_dists_chunk(:, :)

    call get_proc_info(proc_id, n_procs)
    call get_proc_elem_range(params % num_states, proc_first_state, proc_states, recv_counts, recv_shifts)

    ! Calculate chunks
    allocate(K_dists_chunk(proc_states, params % J + 1))
    do proc_k = 1, proc_states
      K_dists_chunk(proc_k, :) = calculate_K_dist(params, p_dist(proc_k, :, :, :))
      K_dist_error = 1 - sum(K_dists_chunk(proc_k, :))
      if (abs(K_dist_error) > 1d-10) print *, 'Warning: K_dist_error is ' // num2str(K_dist_error) // ' on state ' // num2str(proc_first_state + proc_k - 1)
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
! The main procedure for calculation of states properties
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_state_properties(params, rho_grid, L, cap)
    class(input_params), intent(in) :: params
    real*8, intent(in) :: rho_grid(:)
    integer, intent(in) :: L ! Num of points along theta
    real*8, intent(in) :: cap(:)
    integer :: N, vdw_n_ind
    integer :: cov_n_inds(3)
    integer, allocatable :: num_solutions_2d(:, :) ! K x n
    integer, allocatable :: num_solutions_1d(:, :, :)
    real*8 :: vdw_max
    real*8, allocatable :: energies_3d(:)
    real*8, allocatable :: region_probs(:, :) ! nstate x 5
    real*8, allocatable :: K_dists(:, :) ! n_state x (J + 1)
    real*8, allocatable :: gammas(:, :) ! n_state x 2
    real*8, allocatable :: p_dist(:, :, :, :) ! n_state x K x n x 3
    character(:), allocatable :: sym_path, spectrum_path
    ! Arrays of size 2 x N x L (or K x N x L if compressing is disabled). Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B
    type(array_2d_real), allocatable :: As(:, :, :), Bs(:, :, :)
    type(array_2d_complex), allocatable :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S

    ! Load energies
    sym_path = get_sym_path_params(params)
    spectrum_path = get_spectrum_path(sym_path)
    energies_3d = load_energies_3D(spectrum_path)

    vdw_max = 11
    N = size(rho_grid)
    ! Load expansion coefficients
    if (params % fix_basis_jk == 1) then
      call load_1D_expansion_coefficients_fixed_basis(params, N, L, As, num_solutions_1d)
      call load_2D_expansion_coefficients_fixed_basis(params, N, L, num_solutions_1d, Bs, num_solutions_2d)
    else
      call load_1D_expansion_coefficients(params, N, L, As, num_solutions_1d)
      call load_2D_expansion_coefficients(params, N, L, num_solutions_1d, Bs, num_solutions_2d)
    end if

    call load_3D_expansion_coefficients(params, N, num_solutions_2d, proc_Cs)
    call print_parallel('Done loading expansion coefficients')

    p_dist = calculate_wf_integral_all_regions(params, As, Bs, proc_Cs)
    call print_parallel('Done calculating probability distributions')

    ! Calculate properties
    gammas = calculate_gamma_channels_all(params, p_dist, cap)
    call print_parallel('Done calculating channel-specific gammas')

    call find_barriers_n_indices(params, rho_grid, vdw_max, cov_n_inds, vdw_n_ind)
    region_probs = calculate_pes_region_probabilities_all(params, p_dist, cov_n_inds, vdw_n_ind)
    call print_parallel('Done calculating PES region probabilities')

    K_dists = calculate_K_dist_all(params, p_dist)
    call print_parallel('Done calculating K-distributions')

    call print_parallel('Writing results...')
    call write_state_properties(params, energies_3d, gammas * autown, region_probs, K_dists)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! The three functions below can be used for debug purposes to check orthnormality of 1D, 2D and 3D solutions
!-------------------------------------------------------------------------------------------------------------------------------------------

! !-------------------------------------------------------------------------------------------------------------------------------------------
! ! Checks orthonormality of 1D solutions
! !-------------------------------------------------------------------------------------------------------------------------------------------
!   function check_orthonormality_1D(As) result(max_deviation)
!     type(array_2d_real), intent(in) :: As(:, :, :)
!     real*8 :: max_deviation ! from orthonormality
!     integer :: K_ind, n, l, i1, i2
!     real*8 :: overlap
!     real*8, allocatable :: solution_set(:, :)
!
!     max_deviation = 0
!     do K_ind = 1, size(As, 1)
!       do n = 1, size(As, 2)
!         do l = 1, size(As, 3)
!           solution_set = As(K_ind, n, l) % p
!           do i1 = 1, size(solution_set, 2)
!             do i2 = i1, size(solution_set, 2)
!               overlap = sum(solution_set(:, i1) * solution_set(:, i2))
!               max_deviation = max(max_deviation, abs(0 ** (i2 - i1) - overlap))
!             end do
!           end do
!         end do
!       end do
!     end do
!   end function
!
! !-------------------------------------------------------------------------------------------------------------------------------------------
! ! Checks orthonormality of 2D solutions
! !-------------------------------------------------------------------------------------------------------------------------------------------
!   function check_orthonormality_2D(As, Bs) result(max_deviation)
!     type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :)
!     real*8 :: max_deviation ! from orthonormality
!     integer :: K_ind, n, l, m_ind, j1, j2
!     real*8 :: i1_sum, i2_sum, overlap
!
!     max_deviation = 0
!     do K_ind = 1, size(As, 1)
!       do n = 1, size(As, 2)
!         print *, 'Doing n: ', n
!         do j1 = 1, size(Bs(K_ind, n, 1) % p, 2)
!           do j2 = j1, size(Bs(K_ind, n, 1) % p, 2)
!
!             overlap = 0
!             do l = 1, size(As, 3)
!               do m_ind = 1, size(As(K_ind, n, l) % p, 1)
!                 i1_sum = sum(Bs(K_ind, n, l) % p(:, j1) * As(K_ind, n, l) % p(m_ind, :))
!                 i2_sum = sum(Bs(K_ind, n, l) % p(:, j2) * As(K_ind, n, l) % p(m_ind, :))
!                 overlap = overlap + i1_sum * i2_sum
!               end do
!             end do
!             max_deviation = max(max_deviation, abs(0 ** (j2 - j1) - overlap))
!           end do
!         end do
!       end do
!     end do
!   end function
!
! !-------------------------------------------------------------------------------------------------------------------------------------------
! ! Checks orthonormality of 3D solutions. Returns deviation from expected value (0 for different ks, 1 for same ks).
! !-------------------------------------------------------------------------------------------------------------------------------------------
!   function check_orthonormality_3D(As, Bs, Cs, params, k1_first, k2_first, k1_last, k2_last) result(max_deviation)
!     ! Arrays of size K x N x L. Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B.
!     type(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :)
!     type(array_2d_complex), intent(in) :: Cs(:, :) ! 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S.
!     class(input_params), intent(in) :: params
!     integer, intent(in) :: k1_first, k2_first
!     integer, optional, intent(in) :: k1_last, k2_last
!     real*8 :: max_deviation
!     integer :: K_val, K_sym, K_ind, K_ind_comp, n, l, m_ind, j, k1, k2, k1_last_act, k2_last_act
!     real*8 :: i_sum
!     complex*16 :: j1_sum, j2_sum
!     complex*16 :: overlap
!
!     k1_last_act = arg_or_default(k1_last, k1_first)
!     k2_last_act = arg_or_default(k2_last, k2_first)
!     max_deviation = 0
!     do k1 = k1_first, k1_last_act ! size(Cs(1, 1) % p, 2)
!       do k2 = k2_first, k2_last_act ! k1, size(Cs(1, 1) % p, 2)
!         if (k2 < k1) then
!           print *, 'Skipping ovelap between', k1, k2, 'since k1 has to be >= than k2'
!           cycle
!         end if
!
!         print *, 'Calculating overlap between', k1, k2
!         overlap = 0
!         do K_val = params % K(1), params % K(2)
!           call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_comp)
!           do n = 1, size(As, 2)
!             ! call track_progress(n * 1d0 / size(As, 2))
!             do l = 1, size(As, 3)
!               do m_ind = 1, size(As(K_ind_comp, n, l) % p, 1)
!                 j1_sum = 0
!                 j2_sum = 0
!                 do j = 1, size(Bs(K_ind_comp, n, l) % p, 2)
!                   i_sum = sum(Bs(K_ind_comp, n, l) % p(:, j) * As(K_ind_comp, n, l) % p(m_ind, :))
!                   j1_sum = j1_sum + conjg(Cs(K_ind, n) % p(j, k1)) * i_sum
!                   j2_sum = j2_sum + Cs(K_ind, n) % p(j, k2) * i_sum
!                 end do
!                 overlap = overlap + j1_sum * j2_sum
!               end do
!             end do
!           end do
!         end do
!
!         max_deviation = max(max_deviation, abs(0 ** (k2 - k1) - overlap))
!       end do
!     end do
!   end function

end module
