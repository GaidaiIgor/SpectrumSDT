!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures related to post-processing of wave functions obtained during the `eigencalc` stage
!-------------------------------------------------------------------------------------------------------------------------------------------
module state_properties_mod
  use algorithms_mod
  use array_1d_mod
  use array_2d_mod
  use constants
  use general_utils
  use input_params_mod
  use io_utils
  use iso_fortran_env, only: real64
  use mpi
  use parallel_utils
  use path_utils
  use rovib_io_mod
  use rovib_utils_mod
  use vector_mod
  implicit none

  private
  public :: calculate_state_properties

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Packs phi borders of all sections into an array.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function pack_phi_borders(params) result(phi_borders)
    class(input_params), intent(in) :: params
    real(real64), allocatable :: phi_borders(:)
    integer :: i

    allocate(phi_borders(2*size(params % wf_sections)))
    do i = 1, size(params % wf_sections)
      phi_borders(2*i-1 : 2*i) = params % wf_sections(i) % phi
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Maps wf_section *borders* to indices of corresponding points on *grid*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function map_borders_to_point_indices(borders, grid) result(inds)
    real(real64), intent(in) :: borders(2)
    real(real64), intent(in) :: grid(:)
    integer :: inds(2)

    inds(1) = findloc(borders(1) .ale. grid, .true., dim = 1)
    inds(2) = findloc(borders(2) .age. grid, .true., dim = 1, back = .true.)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Maps each wf_section to a set of indices of `p_dist`.
! *wf_sections_dist_ind* - a processed form of wave function sections. It is a 3D array, where the 1st dim - section index,
! 2nd dim - parameter index: 1 - K, 2 - rho, 3 - theta, 4 - phi,
! 3rd dim - start and end indices in `p_dist`, corresponding to the borders of a given parameter of a given section.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_wf_sections_dist_inds(params, rho_grid, theta_grid, phi_borders) result(wf_sections_dist_inds)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: rho_grid(:), theta_grid(:), phi_borders(:)
    integer, allocatable :: wf_sections_dist_inds(:, :, :)
    integer :: i

    allocate(wf_sections_dist_inds(size(params % wf_sections), 4, 2))
    do i = 1, size(params % wf_sections)
      wf_sections_dist_inds(i, 1, :) = params % wf_sections(i) % K - params % K(1) + 1
      wf_sections_dist_inds(i, 2, :) = map_borders_to_point_indices(params % wf_sections(i) % rho, rho_grid)
      wf_sections_dist_inds(i, 3, :) = map_borders_to_point_indices(params % wf_sections(i) % theta, theta_grid)
      wf_sections_dist_inds(i, 4, :) = findloc_all(phi_borders, params % wf_sections(i) % phi)
      ! i-th element of p_dist is integral from i-th to i+1 border, so the second index has to be reduced by 1 to make the range inclusive.
      wf_sections_dist_inds(i, 4, 2) = wf_sections_dist_inds(i, 4, 2) - 1
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Processes wf_sections described in the *params* to generate a border array for phi and map wf_sections to `p_dist` indices.
! *phi_borders* - array of phi values corresponding to integration borders.
! See `get_wf_sections_dist_inds` for details of *wf_sections_dist_inds*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine process_wf_sections(params, rho_grid, theta_grid, phi_borders, wf_sections_dist_inds)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)
    real(real64), allocatable, intent(out) :: phi_borders(:)
    integer, allocatable, intent(out) :: wf_sections_dist_inds(:, :, :)

    phi_borders = pack_phi_borders(params)
    phi_borders = get_unique(phi_borders)
    wf_sections_dist_inds = get_wf_sections_dist_inds(params, rho_grid, theta_grid, phi_borders)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Evaluates the value of phi integral in the range [a, b] (F~_Ks(m,m') = int(F_m1^+-(phi) * F_m2^+-(phi), a, b)).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_phi_integral(m1, m2, a, b, F_sym) result(integral_value)
    integer, intent(in) :: m1, m2 ! Frequency of integrand functions
    real(real64), intent(in) :: a, b ! Range of integration
    integer, intent(in) :: F_sym ! Symmetry of integrand functions: 0 or 1 (same for both)
    real(real64) :: integral_value ! Result
    integer :: sign
    real(real64) :: norm

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
! Calculates m-sum in wf integral.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_m_sum(params, As, K_val, n, l, phi_range, j_sums) result(m_sum)
    class(input_params), intent(in) :: params
    class(array_2d_real), intent(in) :: As(:, :, :)
    integer, intent(in) :: K_val, n, l
    real(real64), intent(in) :: phi_range(2)
    complex(real64), intent(in) :: j_sums(:)
    real(real64) :: m_sum
    integer :: K_ind, K_sym, K_ind_comp, m1, m2, m1_ind, m2_ind
    real(real64) :: product_coeff, phi_integral

    call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_comp)
    m_sum = 0
    if ((phi_range(1) .aeq. 0d0) .and. (phi_range(2) .aeq. 2*pi)) then
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
! Calculates <Psi|Psi> expressions in all regions, for all proc states.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_wf_integral_all_regions(params, As, Bs, proc_Cs, phi_borders) result(p_dist)
    class(input_params), intent(in) :: params
    ! Arrays of size 2 x N x L (or K x N x L if compression is disabled). Each element is a matrix with expansion coefficients - M x S_Knl for A (1D), S_Knl x S_Kn for B (2D)
    class(array_2d_real), intent(in) :: As(:, :, :), Bs(:, :, :) 
    class(array_2d_complex), intent(in) :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    real(real64), intent(in) :: phi_borders(:)
    real(real64), allocatable :: p_dist(:, :, :, :, :)
    integer :: proc_first_state, proc_states, proc_state_ind, K_val, K_ind, K_ind_comp, K_sym, n, l, phi_group_ind, m_ind, j
    real(real64) :: i_sum
    real(real64) :: phi_range(2)
    complex(real64) :: j_sums(params % basis_size_phi) ! j-sums for different ms

    call get_proc_elem_range(params % num_states, proc_first_state, proc_states)
    allocate(p_dist(proc_states, params % K(2) - params % K(1) + 1, size(As, 2), size(As, 3), size(phi_borders) - 1))
    p_dist = 0
    do proc_state_ind = 1, size(p_dist, 1)
      do K_val = params % K(1), params % K(2)
        call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_comp)
        do n = 1, size(As, 2)
          ! Skip empty n-slices
          if (.not. allocated(Bs(K_ind_comp, n, 1) % p)) then
            cycle
          end if

          do l = 1, size(As, 3)
            do phi_group_ind = 1, size(p_dist, 5)
              phi_range = phi_borders(phi_group_ind : phi_group_ind + 1)
              j_sums = 0
              do m_ind = 1, size(As(K_ind_comp, n, l) % p, 1)
                do j = 1, size(Bs(K_ind_comp, n, l) % p, 2)
                  i_sum = sum(Bs(K_ind_comp, n, l) % p(:, j) * As(K_ind_comp, n, l) % p(m_ind, :))
                  j_sums(m_ind) = j_sums(m_ind) + proc_Cs(K_ind, n) % p(j, proc_state_ind) * i_sum
                end do
              end do
              p_dist(proc_state_ind, K_ind, n, l, phi_group_ind) = p_dist(proc_state_ind, K_ind, n, l, phi_group_ind) + calc_m_sum(params, As, K_val, n, l, phi_range, j_sums)
            end do
          end do
        end do
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates total probability of k-th state, for all specified wave function sections in *wf_sections_dist_ind*.
! *p_dist_k* - probability of k-th state in each region defined by user borders.
! *wf_sections_dist_ind* - processed wf_sections, where borders are mapped to p_dist indices.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_wf_section_probabilities(p_dist_k, wf_sections_dist_inds) result(section_probs)
    real(real64), intent(in) :: p_dist_k(:, :, :, :)
    integer, intent(in) :: wf_sections_dist_inds(:, :, :)
    real(real64) :: section_probs(size(wf_sections_dist_inds, 1))
    integer :: i

    associate(p => wf_sections_dist_inds)
      do i = 1, size(section_probs)
        section_probs(i) = sum(p_dist_k(p(i, 1, 1):p(i, 1, 2), p(i, 2, 1):p(i, 2, 2), p(i, 3, 1):p(i, 3, 2), p(i, 4, 1):p(i, 4, 2)))
      end do
    end associate
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calls calculate_wf_section_probabilities for all states in parallel and gathers results.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_wf_section_probabilities_all(params, p_dist, wf_sections_dist_inds, section_probs)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: p_dist(:, :, :, :, :) ! state x sizes of all sections in order K, rho, theta, phi
    integer, intent(in) :: wf_sections_dist_inds(:, :, :)
    real(real64), allocatable, intent(out) :: section_probs(:, :)
    integer :: proc_first_state, proc_states, proc_k
    integer, allocatable :: recv_counts(:), recv_shifts(:)
    real(real64), allocatable :: section_probs_chunk(:, :)

    call get_proc_elem_range(params % num_states, proc_first_state, proc_states, recv_counts, recv_shifts)
    ! Calculate local chunks
    allocate(section_probs_chunk(proc_states, size(params % wf_sections)))
    do proc_k = 1, proc_states
      section_probs_chunk(proc_k, :) = calculate_wf_section_probabilities(p_dist(proc_k, :, :, :, :), wf_sections_dist_inds)
    end do
    call gather_chunks(section_probs_chunk, recv_counts, recv_shifts, section_probs)
  end subroutine

! !-------------------------------------------------------------------------------------------------------------------------------------------
! ! Calculates channel-specific values of gamma using probability distributions of a given state.
! !-------------------------------------------------------------------------------------------------------------------------------------------
!   function calculate_gamma_channels(p_dist_k, cap) result(gammas)
!     real(real64), intent(in) :: p_dist_k(:, :, :) ! K x n x 3
!     real(real64), intent(in) :: cap(:)
!     real(real64) :: gammas(2) ! 1 - channel A, 2 - channel B
!     integer :: K_ind
!
!     call assert(size(cap) == size(p_dist_k, 2), 'Size of cap array does not match p_dist')
!     gammas = 0
!     do K_ind = 1, size(p_dist_k, 1)
!       gammas(1) = gammas(1) + sum(p_dist_k(K_ind, :, 2) * cap)
!       gammas(1) = gammas(1) + sum(p_dist_k(K_ind, :, 3) * cap)
!       gammas(2) = gammas(2) + sum(p_dist_k(K_ind, :, 1) * cap)
!     end do
!     ! One factor of 2 is due to phi symmetry, another is due to gamma definition
!     gammas = gammas * 2 * 2
!   end function
!
! !-------------------------------------------------------------------------------------------------------------------------------------------
! ! Calculates channel-specific gammas for all states.
! !-------------------------------------------------------------------------------------------------------------------------------------------
!   subroutine calculate_gamma_channels_all(params, p_dist, cap, gammas)
!     class(input_params), intent(in) :: params
!     real(real64), intent(in) :: p_dist(:, :, :, :)
!     real(real64), intent(in) :: cap(:)
!     real(real64), allocatable, intent(out) :: gammas(:, :)
!     integer :: proc_first_state, proc_states, proc_k
!     integer, allocatable :: recv_counts(:), recv_shifts(:)
!     real(real64), allocatable :: gammas_chunk(:, :)
!
!     call get_proc_elem_range(params % num_states, proc_first_state, proc_states, recv_counts, recv_shifts)
!     ! Calculate proc chunks
!     allocate(gammas_chunk(proc_states, 2))
!     do proc_k = 1, proc_states
!       gammas_chunk(proc_k, :) = calculate_gamma_channels(p_dist(proc_k, :, :, :), cap)
!     end do
!     call gather_chunks(gammas_chunk, recv_counts, recv_shifts, gammas)
!   end subroutine
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! The main procedure for calculation of states properties.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_state_properties(params, rho_grid, theta_grid, cap)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: rho_grid(:), theta_grid(:)
    real(real64), optional, intent(in) :: cap(:)
    integer :: N, L
    integer, allocatable :: num_solutions_2d(:, :) ! K x n
    integer, allocatable :: num_solutions_1d(:, :, :), wf_sections_dist_inds(:, :, :)
    real(real64), allocatable :: energies_3d(:), phi_borders(:)
    real(real64), allocatable :: section_probs(:, :) ! num_states x size(wf_sections)
    ! real(real64), allocatable :: K_dists(:, :) ! num_states x (J + 1)
    real(real64), allocatable :: gammas(:, :) ! num_states x 2
    real(real64), allocatable :: p_dist(:, :, :, :, :) ! num_states x each section in K, rho, theta, phi
    character(:), allocatable :: sym_path, spectrum_path
    ! Arrays of size 2 x N x L (or K x N x L if compressing is disabled). Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B
    type(array_2d_real), allocatable :: As(:, :, :), Bs(:, :, :)
    type(array_2d_complex), allocatable :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S

    ! Load energies
    sym_path = get_sym_path(params)
    spectrum_path = get_spectrum_path(sym_path)
    energies_3d = load_energies_3D(spectrum_path)

    N = size(rho_grid)
    L = size(theta_grid)
    ! Load expansion coefficients
    if (params % use_fixed_basis_jk == 1) then
      call load_1D_expansion_coefficients_fixed_basis(params, N, L, As, num_solutions_1d)
      call load_2D_expansion_coefficients_fixed_basis(params, N, L, num_solutions_1d, Bs, num_solutions_2d)
    else
      call load_1D_expansion_coefficients(params, N, L, As, num_solutions_1d)
      call load_2D_expansion_coefficients(params, N, L, num_solutions_1d, Bs, num_solutions_2d)
    end if

    call load_3D_expansion_coefficients(params, N, num_solutions_2d, proc_Cs)
    call print_parallel('Done loading expansion coefficients')

    call process_wf_sections(params, rho_grid, theta_grid, phi_borders, wf_sections_dist_inds)
    p_dist = calculate_wf_integral_all_regions(params, As, Bs, proc_Cs, phi_borders)
    call print_parallel('Done calculating probability distributions')

    ! if (present(cap)) then
    !   call calculate_gamma_channels_all(params, p_dist, cap, gammas)
    ! else
      allocate(gammas(params % num_states, 2))
      gammas = 0
    ! end if
    ! call print_parallel('Done calculating channel-specific gammas')

    call calculate_wf_section_probabilities_all(params, p_dist, wf_sections_dist_inds, section_probs)
    call print_parallel('Done calculating wave functions section probabilities')

    call print_parallel('Writing results...')
    call write_state_properties(params, energies_3d, section_probs, gammas * au_to_wn)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! The three functions below can be used for debug purposes to check orthnormality of 1D, 2D and 3D solutions
!-------------------------------------------------------------------------------------------------------------------------------------------

! !-------------------------------------------------------------------------------------------------------------------------------------------
! ! Checks orthonormality of 1D solutions
! !-------------------------------------------------------------------------------------------------------------------------------------------
!   function check_orthonormality_1D(As) result(max_deviation)
!     type(array_2d_real), intent(in) :: As(:, :, :)
!     real(real64) :: max_deviation ! from orthonormality
!     integer :: K_ind, n, l, i1, i2
!     real(real64) :: overlap
!     real(real64), allocatable :: solution_set(:, :)
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
!     real(real64) :: max_deviation ! from orthonormality
!     integer :: K_ind, n, l, m_ind, j1, j2
!     real(real64) :: i1_sum, i2_sum, overlap
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
!     real(real64) :: max_deviation
!     integer :: K_val, K_sym, K_ind, K_ind_comp, n, l, m_ind, j, k1, k2, k1_last_act, k2_last_act
!     real(real64) :: i_sum
!     complex(real64) :: j1_sum, j2_sum
!     complex(real64) :: overlap
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
