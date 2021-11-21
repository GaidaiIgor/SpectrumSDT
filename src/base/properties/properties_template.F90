#include "funcs.macro"
  use properties_base_mod
  implicit none

  interface calculate_wf_integral_all_regions
    module procedure :: CONCAT2(calculate_wf_integral_all_regions_,TEMPLATE_TYPE_NAME)
  end interface

  private
  public :: CONCAT2(calculate_state_properties_,TEMPLATE_TYPE_NAME)

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates <Psi|Psi> expressions in all regions, for all proc states.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function CONCAT2(calculate_wf_integral_all_regions_,TEMPLATE_TYPE_NAME)(params, As, Bs, proc_Cs, phi_borders) result(proc_p_dist)
    class(input_params), intent(in) :: params
    ! Arrays of size (1 or 2) x N x L (or K x N x L if compression is disabled). Each element is a matrix with expansion coefficients - nphi_total x S_Knl for A (1D), S_Knl x S_Kn for B (2D)
    class(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), intent(in) :: As(:, :, :), Bs(:, :, :) 
    class(array_2d_complex), intent(in) :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S
    real(real64), intent(in) :: phi_borders(:)
    real(real64), allocatable :: proc_p_dist(:, :, :, :, :)
    logical :: use_simplified_method
    integer :: nphi_basis_type, nphi_total, proc_first_state, proc_states, K_ind, K1_sym, K_sym, K_ind_smart, proc_state_ind, K_val, n, l, phi_matrix_ind, phi_range_ind, m_ind, j
    real(real64) :: i_sum, m_sum
    real(real64) :: phi_range(2)
    type(array_3d_real) :: phi_integral_matrix(2) ! Phi-integal matrices for each symmetry of phi basis functions. Up to 2 may be needed if phi symmetry changes with K.
    complex(real64), allocatable :: j_sums(:) ! j-sums for different ms
    character(:), allocatable :: molecule_type

    nphi_basis_type = params % get_num_funcs_phi_per_basis_type()
    molecule_type = get_molecule_type(params % mass)
    call get_k_attributes(params % K(1), params, K_sym = K1_sym)
    call get_proc_elem_range(params % eigensolve % num_states, proc_first_state, proc_states)
    allocate(proc_p_dist(proc_states, params % K(2) - params % K(1) + 1, size(As, 2), size(As, 3), size(phi_borders) - 1))
    proc_p_dist = 0

    ! Use simplified method if only the whole range [0; 2*pi] is of interest
    use_simplified_method = size(phi_borders) == 2 .and. (phi_borders(1) .aeq. 0d0) .and. (phi_borders(2) .aeq. 2*pi)
    if (.not. use_simplified_method) then
      ! Calculate necessary phi-integral matrices
      phi_integral_matrix(1) % p = calc_phi_integral_matrix(nphi_basis_type, K1_sym, molecule_type, params % use_geometric_phase, params % basis % use_half_integers, phi_borders)
      if (params % use_geometric_phase == 0 .and. params % K(2) > params % K(1)) then
        call get_k_attributes(params % K(1) + 1, params, K_sym = K_sym)
        phi_integral_matrix(2) % p = calc_phi_integral_matrix(nphi_basis_type, K_sym, molecule_type, params % use_geometric_phase, params % basis % use_half_integers, phi_borders)
      end if
    end if

    nphi_total = params % get_num_funcs_phi_total()
    allocate(j_sums(nphi_total))
    do proc_state_ind = 1, size(proc_p_dist, 1)
      do K_val = params % K(1), params % K(2)
        call get_k_attributes(K_val, params, K_ind, K_sym, K_ind_smart)
        phi_matrix_ind = iff(K_sym == K1_sym, 1, 2)
        do n = 1, size(proc_p_dist, 3)
          ! Skip empty n-slices
          if (.not. allocated(Bs(K_ind_smart, n, 1) % p)) then
            cycle
          end if

          do l = 1, size(proc_p_dist, 4)
            j_sums = 0
            do m_ind = 1, size(As(K_ind_smart, n, l) % p, 1)
              do j = 1, size(Bs(K_ind_smart, n, l) % p, 2)
                i_sum = dot_product(Bs(K_ind_smart, n, l) % p(:, j), As(K_ind_smart, n, l) % p(m_ind, :))
                j_sums(m_ind) = j_sums(m_ind) + proc_Cs(K_ind, n) % p(j, proc_state_ind) * i_sum
              end do
            end do

            do phi_range_ind = 1, size(proc_p_dist, 5)
              if (use_simplified_method) then
                m_sum = real(dot_product(j_sums, j_sums))
              else
                m_sum = calculate_vT_A_v_symmetric_real(phi_integral_matrix(phi_matrix_ind) % p(:, :, phi_range_ind), j_sums)
              end if
              proc_p_dist(proc_state_ind, K_ind, n, l, phi_range_ind) = proc_p_dist(proc_state_ind, K_ind, n, l, phi_range_ind) + m_sum
            end do
          end do
        end do
      end do

      if (get_proc_id() == 0) then
        call track_progress(proc_state_ind * 1d0 / size(proc_p_dist, 1), 1d-2)
      end if
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! The main procedure for calculation of states properties.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine CONCAT2(calculate_state_properties_,TEMPLATE_TYPE_NAME)(params, rho_info, theta_info, cap)
    class(input_params), intent(in) :: params
    class(grid_info), intent(in) :: rho_info, theta_info
    real(real64), optional, intent(in) :: cap(:)
    integer :: N, L
    integer, allocatable :: num_solutions_2d(:, :) ! K x n
    integer, allocatable :: num_solutions_1d(:, :, :)
    real(real64), allocatable :: phi_borders(:)
    real(real64), allocatable :: eigenvalues_3d(:, :)
    real(real64), allocatable :: section_stats(:, :) ! rows - states, columns - statistics
    real(real64), allocatable :: wf_sections_dist_mask(:, :, :, :, :) ! For each wf_section (1st dim) stores weights of each point in K x rho x theta x phi grid corresponding to that section
    real(real64), allocatable :: proc_p_dist(:, :, :, :, :) ! For each state (1st dim) stores probabilities in each point of K x rho x theta x phi grid
    character(:), allocatable :: states_path
    ! Arrays of size 2 x N x L (or K x N x L for adiabatic basis). Each element is a matrix with expansion coefficients - M x S_Knl for A, S_Knl x S_Kn for B
    type(CONCAT2(array_2d_,TEMPLATE_TYPE_NAME)), allocatable :: As(:, :, :), Bs(:, :, :)
    type(array_2d_complex), allocatable :: proc_Cs(:, :) ! local 3D expansion coefficients K x n. Inner expansion matrix is S_Kn x S

    ! Load energies and gammas
    states_path = get_states_path(params)
    eigenvalues_3d = read_matrix_real(states_path, skip_lines = 1)
    call assert(size(eigenvalues_3d, 1) >= params % eigensolve % num_states, 'Error: not enough eigenstates computed')

    N = size(rho_info % points)
    L = size(theta_info % points)
    ! Load expansion coefficients
    call print_parallel('Loading expansion coefficients...')
    call load_1D_expansion_coefficients(params, N, L, As, num_solutions_1d)
    call load_2D_expansion_coefficients(params, N, L, num_solutions_1d, Bs, num_solutions_2d)
    call load_3D_expansion_coefficients(params, N, num_solutions_2d, proc_Cs)

    call print_parallel('Calculating probability distributions...')
    call process_wf_sections(params, rho_info, theta_info, cap, phi_borders, wf_sections_dist_mask)
    proc_p_dist = calculate_wf_integral_all_regions(params, As, Bs, proc_Cs, phi_borders)

    call print_parallel('Calculating statistics...')
    call calculate_wf_sections_statistics_all(params, wf_sections_dist_mask, proc_p_dist, section_stats)

    call print_parallel('Writing results...')
    call write_state_properties(params, eigenvalues_3d(1:size(section_stats, 1), :), section_stats)
  end subroutine

