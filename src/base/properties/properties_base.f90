!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures related to post-processing of wave functions obtained during the `eigensolve` stage.
!-------------------------------------------------------------------------------------------------------------------------------------------
module properties_base_mod
  use algorithms_mod
  use array_1d_mod
  use array_2d_mod
  use constants
  use general_utils_mod
  use grid_info_mod
  use input_params_mod
  use io_utils_mod
  use iso_fortran_env, only: real64
  use mpi
  use parallel_utils_mod
  use path_utils_mod
  use spectrumsdt_io_mod
  use spectrumsdt_utils_mod
  use rovib_utils_mod
  use vector_mod
  implicit none

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

    inds(1) = minloc(abs(borders(1) - grid), dim = 1)
    inds(2) = minloc(abs(borders(2) - grid), dim = 1)
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
! Returns weight of the border points in `p_dist`.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_border_weights(borders, grid, grid_ends) result(weights)
    real(real64), intent(in) :: borders(2)
    real(real64), intent(in) :: grid(:)
    real(real64), intent(in) :: grid_ends(2)
    real(real64) :: weights(2)
    integer :: i
    integer :: border_inds(2)
    real(real64) :: closest_points(2), edges_left(2), edges_right(2), refs(2), dists(2)
    real(real64), allocatable :: grid_ext(:)

    ! Extend grid with virtual points so that midpoints coincide with grid ends. Makes handling special cases with grid ends less tedious.
    grid_ext = [grid(1) - 2*(grid(1) - grid_ends(1)), grid, grid(size(grid)) + 2*(grid_ends(2) - grid(size(grid)))]
    border_inds = map_borders_to_point_indices(borders, grid) + 1 ! Find on old grid and shift to avoid finding virtual points
    closest_points = grid_ext(border_inds)
    edges_left = (closest_points + grid_ext(border_inds - 1)) / 2
    edges_right = (closest_points + grid_ext(border_inds + 1)) / 2
    refs = [(iff(borders(i) < closest_points(i), edges_left(i), edges_right(i)), i = 1, 2)]

    if (border_inds(1) == border_inds(2)) then
      ! If on borders are on the same side from the grid point
      if ((closest_points(1) - borders(1)) * (closest_points(1) - borders(2)) > 0) then
        ! sqrt since they are going to multiply the same point in p_dist
        weights = sqrt((borders(2) - borders(1)) / abs(closest_points(1) - refs(1)) / 2)
      ! Borders are on different sides from the grid point
      else
        weights(1) = (closest_points(1) - borders(1)) / (closest_points(1) - refs(1)) / 2
        weights(2) = (borders(2) - closest_points(2)) / (refs(2) - closest_points(2)) / 2
        weights = sqrt(sum(weights))
      end if
    else
      dists(1) = iff(borders(1) < closest_points(1), closest_points(1) - borders(1), refs(1) - borders(1))
      dists(2) = iff(borders(2) > closest_points(2), borders(2) - closest_points(2), borders(2) - refs(2))
      weights(1) = iff(borders(1) < closest_points(1), 0.5d0, 0d0)
      weights(2) = iff(borders(2) > closest_points(2), 0.5d0, 0d0)
      weights = weights + dists / abs(closest_points - refs) / 2
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns masks corresponding to wave function sections. The masks define weights for each element of `p_dist` (and has the same size).
! The first dimension of section_masks corresponds to different sections. 
! The remaining 4 have the same meaning as dimensions 2-5 of `p_dist`.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function generate_wf_sections_dist_mask(params, rho_info, theta_info, phi_borders, cap) result(wf_sections_dist_mask)
    class(input_params), intent(in) :: params
    class(grid_info), intent(in) :: rho_info, theta_info
    real(real64), intent(in) :: phi_borders(:)
    real(real64), optional, intent(in) :: cap(:)
    real(real64), allocatable :: wf_sections_dist_mask(:, :, :, :, :)
    integer :: i, j
    integer, allocatable :: section_inds(:, :, :)
    real(real64) :: weights_rho(2), weights_theta(2)

    allocate(wf_sections_dist_mask(size(params % wf_sections), params % K(2) - params % K(1) + 1, size(rho_info % points), size(theta_info % points), size(phi_borders) - 1))
    wf_sections_dist_mask = 0
    section_inds = get_wf_sections_dist_inds(params, rho_info % points, theta_info % points, phi_borders)

    associate(si => section_inds)
      do i = 1, size(wf_sections_dist_mask, 1)
        wf_sections_dist_mask(i, si(i, 1, 1):si(i, 1, 2), si(i, 2, 1):si(i, 2, 2), si(i, 3, 1):si(i, 3, 2), si(i, 4, 1):si(i, 4, 2)) = 1

        weights_rho = get_border_weights(params % wf_sections(i) % rho, rho_info % points, [rho_info % from, rho_info % to])
        wf_sections_dist_mask(i, :, si(i, 2, 1), :, :) = wf_sections_dist_mask(i, :, si(i, 2, 1), :, :) * weights_rho(1)
        wf_sections_dist_mask(i, :, si(i, 2, 2), :, :) = wf_sections_dist_mask(i, :, si(i, 2, 2), :, :) * weights_rho(2)

        weights_theta = get_border_weights(params % wf_sections(i) % theta, theta_info % points, [theta_info % from, theta_info % to])
        wf_sections_dist_mask(i, :, :, si(i, 3, 1), :) = wf_sections_dist_mask(i, :, :, si(i, 3, 1), :) * weights_theta(1)
        wf_sections_dist_mask(i, :, :, si(i, 3, 2), :) = wf_sections_dist_mask(i, :, :, si(i, 3, 2), :) * weights_theta(2)

        if (params % wf_sections(i) % stat == 'gamma') then
          call assert(present(cap), 'Error: cap has to be given for gamma statistics')
          call assert(size(cap) == size(rho_info % points))
          do j = 1, size(wf_sections_dist_mask, 3)
            wf_sections_dist_mask(i, :, j, :, :) = wf_sections_dist_mask(i, :, j, :, :) * 2 * cap(j) ! 2 due to definition of gamma
          end do
        end if
      end do
    end associate
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Processes wf_sections described in the *params* to generate a border array for phi and `p_dist` masks for each requested wf property.
! *phi_borders* - array of phi values corresponding to integration borders.
! See `generate_wf_sections_dist_mask` for details of *wf_sections_dist_mask*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine process_wf_sections(params, rho_info, theta_info, cap, phi_borders, wf_sections_dist_mask)
    class(input_params), intent(in) :: params
    class(grid_info), intent(in) :: rho_info, theta_info
    real(real64), optional, intent(in) :: cap(:)
    real(real64), allocatable, intent(out) :: phi_borders(:)
    real(real64), allocatable, intent(out) :: wf_sections_dist_mask(:, :, :, :, :)

    phi_borders = pack_phi_borders(params)
    phi_borders = get_unique(phi_borders)
    wf_sections_dist_mask = generate_wf_sections_dist_mask(params, rho_info, theta_info, phi_borders, cap)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calls calculate_wf_section_probabilities for all states in parallel and gathers results.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_wf_sections_statistics_all(params, wf_sections_dist_mask, proc_p_dist, section_stats)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: wf_sections_dist_mask(:, :, :, :, :), proc_p_dist(:, :, :, :, :) ! wf_section/state x sizes of all sections in order K, rho, theta, phi
    real(real64), allocatable, intent(out) :: section_stats(:, :)
    integer :: proc_first_state, proc_states, i, j
    real(real64), allocatable :: proc_section_stats(:, :)

    call get_proc_elem_range(params % eigensolve % num_states, proc_first_state, proc_states)
    call assert(size(proc_p_dist, 1) == proc_states, 'Error: wrong size of proc_p_dist')
    ! Calculate local chunks
    allocate(proc_section_stats(proc_states, size(params % wf_sections)))
    do j = 1, size(proc_section_stats, 2)
      do i = 1, size(proc_section_stats, 1)
        proc_section_stats(i, j) = sum(wf_sections_dist_mask(j, :, :, :, :) * proc_p_dist(i, :, :, :, :))
      end do
    end do
    call gather_rows_array_2d_real(proc_section_stats, section_stats)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Evaluates the value of phi integral for functions of given symmetry in the range [a, b].
! m1 and m2 >= 0. m1 or m2 can be 0 only if corresponding F_sym is 0. F_sym is either 0 (cos) or 1 (sin). a <= b. Both a and b are from 0 to 2*pi.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_phi_integral(m1, m2, F1_sym, F2_sym, a, b) result(integral_value)
    integer, intent(in) :: m1, m2 ! Frequency of integrand functions
    real(real64), intent(in) :: a, b ! Range of integration
    integer, intent(in) :: F1_sym, F2_sym ! Symmetry of integrand functions: 0 for cos, 1 for sin
    real(real64) :: integral_value ! Result
    integer :: sign
    real(real64) :: norm

    norm = 1 / (pi * sqrt(1d0 * (delta(m1, 0) + 1) * (delta(m2, 0) + 1)))
    sign = (-1) ** F1_sym
    if (m1 == 0 .and. m2 == 0) then
      integral_value = norm * (b - a)

    else if (m1 == m2) then
      if (F1_sym == F2_sym) then
        integral_value = norm * 0.5d0 * (b - a - sign*(sin(2*a*m1) - sin(2*b*m1)) / (2*m1))
      else
        integral_value = norm * 0.5d0 * (cos(2*a*m1) - cos(2*b*m1)) / (2*m1)
      end if

    else if (m1 /= m2) then
      if (F1_sym == F2_sym) then
        integral_value = norm * 0.5d0 * ((sin(b*(m1 - m2)) - sin(a*(m1 - m2))) / (m1 - m2) + (sign*sin(b*(m1 + m2)) - sign*sin(a*(m1 + m2))) / (m1 + m2))
      else
        integral_value = norm * 0.5d0 * ((-sign*cos(a*(m1 - m2)) + sign*cos(b*(m1 - m2))) / (m1 - m2) + (cos(a*(m1 + m2)) - cos(b*(m1 + m2))) / (m1 + m2))
      end if
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates phi-integral matrix for each phi range and symmetry. Only upper triangle is calculated since matrix is symmetric.
! num_funcs_per_sym - number of phi basis functions per symmetry.
! phi_sym - symmetry. 0 for cos, 1 for sin, 2 for both.
! phi_borders - array of considered phi borders.
! phi_integral_matrix - 3D matrix. 1st dim - m1, 2nd - m2, 3rd - phi_range_ind.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_phi_integral_matrix(num_funcs_phi_per_sym, phi_sym, phi_borders) result(phi_integral_matrix)
    integer, intent(in) :: num_funcs_phi_per_sym, phi_sym
    real(real64), intent(in) :: phi_borders(:)
    real(real64), allocatable :: phi_integral_matrix(:, :, :)
    integer :: num_funcs_phi, range_ind, m1_ind, m2_ind, m1, m2, m1_sym, m2_sym
    real(real64) :: a, b

    num_funcs_phi = get_num_funcs_phi(num_funcs_phi_per_sym, phi_sym)
    allocate(phi_integral_matrix(num_funcs_phi, num_funcs_phi, size(phi_borders) - 1))
    do range_ind = 1, size(phi_integral_matrix, 3)
      a = phi_borders(range_ind)
      b = phi_borders(range_ind + 1)
      do m2_ind = 1, size(phi_integral_matrix, 2)
        call get_m_ind_info(m2_ind, phi_sym, num_funcs_phi_per_sym, m = m2, m_sym = m2_sym)
        do m1_ind = 1, m2_ind
          call get_m_ind_info(m1_ind, phi_sym, num_funcs_phi_per_sym, m = m1, m_sym = m1_sym)
          phi_integral_matrix(m1_ind, m2_ind, range_ind) = calc_phi_integral(m1, m2, m1_sym, m2_sym, a, b)
        end do
      end do
    end do
  end function

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
