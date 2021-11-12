module basis_base_mod
  use algorithms_mod, only: prefix_sum_exclusive
  use array_1d_mod, only: array_1d_real
  use array_2d_mod
  use constants_mod, only: au_to_wn, pi
  use fourier_transform_mod, only: dft_derivative2_optimized_dvr, dft_derivative2_equidistant_dvr, dft_derivative2_equidistant_dvr_analytical
  use general_utils_mod, only: iff, identity_matrix
  use input_params_mod
  use iso_fortran_env, only: real64
  use lapack_interface_mod
  use mpi
  use parallel_utils_mod
  use spectrumsdt_paths_mod
  use spectrumsdt_utils_ext_mod

  use debug_tools_mod
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates phi basis on the grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_phi_basis_grid(params, grid_phi) result(basis)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: grid_phi(:)
    real(real64), allocatable :: basis(:, :)
    integer :: nphi_total, m_ind, m, m_type

    nphi_total = params % get_num_funcs_phi_total()
    allocate(basis(size(grid_phi), nphi_total))
    do m_ind = 1, size(basis, 2)
      call get_m_ind_info(m_ind, params, m, m_type)
      if (m == 0) then
        basis(:, m_ind) = sqrt(0.5d0)
      else if (m_type == 0) then
        basis(:, m_ind) = cos(m * grid_phi)
      else if (m_type == 1) then
        basis(:, m_ind) = sin(m * grid_phi)
      end if
    end do
    basis = basis / sqrt(pi)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes kinetic energy matrix for a particle with given reduced mass *mu* in DVR basis described by the remaining arguments.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function compute_kinetic_energy_dvr(mu, num_points, period, jac) result(matrix)
    real(real64), intent(in) :: mu
    integer, intent(in) :: num_points
    real(real64), intent(in) :: period
    real(real64), optional, intent(in) :: jac(:)
    complex(real64), allocatable :: matrix(:, :)

    if (present(jac)) then
      matrix = dft_derivative2_optimized_dvr(jac, period)
    else
      if (mod(num_points, 2) == 0) then
        matrix = dft_derivative2_equidistant_dvr_analytical(num_points, period)
      else
        matrix = dft_derivative2_equidistant_dvr(num_points, period)
      end if
    end if
    matrix = -matrix / (2 * mu)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates kinetic energy matrix for theta.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function compute_kinetic_theta(mu, rho_val, period_theta, num_points_theta) result(ham)
    real(real64), intent(in) :: mu, rho_val, period_theta
    integer, intent(in) :: num_points_theta
    real(real64), allocatable :: ham(:, :)
    complex(real64), allocatable :: ham_complex(:, :)

    ham_complex = compute_kinetic_energy_dvr(mu, num_points_theta, num_points_theta * period_theta)
    call assert(maxval(abs(aimag(ham_complex))) < 1d-10, 'Error: unexpected imaginary components of equidistant theta DVR')
    ham = 4 / rho_val**2 * real(ham_complex)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Gathers and writes all 1D basis energies from all values of theta and rho.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine gather_write_energies_1d(proc_val1, file_path)
    class(array_1d_real), intent(in) :: proc_val1(:, :)
    character(*), intent(in) :: file_path
    integer :: file_unit, i
    integer, allocatable :: proc_shifts(:), flat_global_index_map(:)
    integer, allocatable :: proc_index_map(:, :), global_index_map(:, :)
    real(real64), allocatable :: flat_proc_val1(:), flat_global_val1(:), val1_section(:)

    call flatten_2d_array_1d_real(proc_val1, flat_proc_val1, proc_index_map)
    call gather_array_1d_real(flat_proc_val1, flat_global_val1, proc_shifts = proc_shifts)
    call gather_cols_array_2d_integer(proc_index_map + proc_shifts(get_proc_id() + 1), global_index_map)

    if (get_proc_id() == 0) then
      flat_global_index_map = [global_index_map, size(flat_global_val1) + 1]
      open(newunit = file_unit, file = file_path)
      write(file_unit, '(2I5)') size(global_index_map, 2), size(global_index_map, 1) ! Number of points in rho and theta
      do i = 1, size(flat_global_index_map) - 1
        val1_section = flat_global_val1(flat_global_index_map(i) : flat_global_index_map(i + 1) - 1)
        write(file_unit, '(' // num2str(size(val1_section)) // 'G25.15)') val1_section * au_to_wn
      end do
      close(file_unit)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Gathers and writes all 2D basis energies from all values of rho.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine gather_write_energies_2d(proc_val2, file_path)
    class(array_1d_real), intent(in) :: proc_val2(:)
    character(*), intent(in) :: file_path
    integer :: file_unit, i
    integer, allocatable :: proc_index_map(:), proc_shifts(:), global_index_map(:)
    real(real64), allocatable :: flat_proc_val2(:), flat_global_val2(:), val2_section(:)

    call flatten_1d_array_1d_real(proc_val2, flat_proc_val2, proc_index_map)
    call gather_array_1d_real(flat_proc_val2, flat_global_val2, proc_shifts = proc_shifts)
    call gather_array_1d_integer(proc_index_map + proc_shifts(get_proc_id() + 1), global_index_map)

    if (get_proc_id() == 0) then
      global_index_map = [global_index_map, size(flat_global_val2) + 1]
      open(newunit = file_unit, file = file_path)
      write(file_unit, '(I5)') size(global_index_map) - 1 ! Number of points in rho
      do i = 1, size(global_index_map) - 1
        val2_section = flat_global_val2(global_index_map(i) : global_index_map(i + 1) - 1)
        write(file_unit, '(' // num2str(size(val2_section)) // 'G25.15)') val2_section * au_to_wn
      end do
      close(file_unit)
    end if
  end subroutine

end module
