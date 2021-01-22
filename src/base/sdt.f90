!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures related to calculation of basis for Hamiltonian and base overlaps.
!-------------------------------------------------------------------------------------------------------------------------------------------
module sdt
  use algorithms_mod, only: prefix_sum_exclusive
  use array_1d_mod, only: array_1d_real
  use array_2d_mod, only: array_2d_real
  use constants, only: au_to_wn, pi
  use formulas_mod, only: get_reduced_mass
  use fourier_transform_mod, only: dft_derivative2_optimized_dvr, dft_derivative2_equidistant_dvr, dft_derivative2_equidistant_dvr_analytical
  use general_utils, only: identity_matrix
  use input_params_mod
  use iso_fortran_env, only: real64
  use lapack_interface_mod
  use mpi
  use parallel_utils
  use potential_mod, only: pottot
  use rovib_io_mod, only: load_basis_size_2d, load_solutions_1D, load_solutions_2D
  use spectrumsdt_paths_mod
  implicit none
  integer, parameter :: nvec1min = 3 ! Min number of 1D states in thread before truncation

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates phi basis on the grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_phi_basis_grid(params, grid_phi) result(basis)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: grid_phi(:)
    real(real64) :: basis(size(grid_phi), params % basis_size_phi)
    integer :: i, j
    real(real64) :: norm

    norm = sqrt(1 / pi)
    do j = 1, size(basis, 2)
      do i = 1, size(basis, 1)
        select case (params % symmetry)
          case (0)
            if (j == 1) then
              basis(i, j) = norm / sqrt(2d0)
            else
              basis(i, j) = norm * cos((j-1) * grid_phi(i))
            end if

          case (1)
            basis(i, j) = norm * sin(j * grid_phi(i))

          case default
            stop 'Impossible error: wrong symmetry value'
        end select
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates 1D Hamiltonian (for phi).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_hamiltonian_1d(params, rho_ind, theta_ind, rho_val, theta_val, grid_phi) result(ham)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind, theta_ind
    real(real64), intent(in) :: rho_val, theta_val
    real(real64), intent(in) :: grid_phi(:)
    real(real64) :: ham(params % basis_size_phi, params % basis_size_phi)
    integer :: i, j, phi_ind, shift
    real(real64) :: mu, coeff, step_phi, sum
    real(real64), allocatable :: basis(:, :)

    basis = get_phi_basis_grid(params, grid_phi)
    ! Coefficient in Hamiltonian operator
    mu = get_reduced_mass(params % mass)
    coeff = -1/(2*mu) * 4 / rho_val**2 / sin(theta_val)**2
    step_phi = grid_phi(2) - grid_phi(1)

    ! Build potential energy matrix
    do j = 1, size(ham, 2)
      do i = 1, size(ham, 1)
        sum = 0
        do phi_ind = 1, size(pottot, 1)
          sum = sum + basis(phi_ind, i)*pottot(phi_ind, theta_ind, rho_ind)*basis(phi_ind, j)
        end do
        ham(i, j) = sum * step_phi
      end do
    end do

    ! Add kinetic energy matrix
    shift = iff(params % symmetry == 0, 1, 0)
    do i = 1, size(ham, 1)
      ! Extra minus comes from basis function derivative
      ham(i, i) = ham(i, i) + coeff * (-1) * (i-shift)**2
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Solves 1D problems for each thread in slice specified by *rho_ind*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_1d(params, rho_ind, rho_val, grid_theta, grid_phi, nvec1, val1, vec1)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind
    real(real64), intent(in) :: rho_val
    real(real64), intent(in) :: grid_theta(:), grid_phi(:)
    integer, allocatable, intent(out) :: nvec1(:) ! Num of 1D vecs in each thread
    type(array_1d_real), allocatable, intent(out) :: val1(:) ! 1D values in each thread
    type(array_2d_real), allocatable, intent(out) :: vec1(:) ! 1D vectors in each thread
    integer :: theta_ind, file_unit
    real(real64), allocatable :: val1_all(:) ! All eigenvalues
    real(real64), allocatable :: ham1(:, :) ! 1D Hamiltonian

    allocate(nvec1(size(grid_theta)), val1(size(grid_theta)), vec1(size(grid_theta)))
    ! Solve eigenvalue problem for each thread
    do theta_ind = 1, size(grid_theta)
      ham1 = get_hamiltonian_1d(params, rho_ind, theta_ind, rho_val, grid_theta(theta_ind), grid_phi)
      call lapack_eigensolver(ham1, val1_all)

      ! Save results
      nvec1(theta_ind) = max(nvec1min, findloc(val1_all < params % cutoff_energy, .true., dim = 1, back = .true.))
      val1(theta_ind) % p = val1_all(:nvec1(theta_ind))
      vec1(theta_ind) % p = ham1(:, :nvec1(theta_ind))
    end do

    ! Save results in binary file
    open(newunit = file_unit, file = get_solutions_1d_path(get_sym_path(params), rho_ind), form = 'unformatted')
    write(file_unit) nvec1
    do theta_ind = 1, size(grid_theta)
      if (nvec1(theta_ind) > 0) then
        write(file_unit) val1(theta_ind) % p
        write(file_unit) vec1(theta_ind) % p
      end if
    end do
    close(file_unit)
  end subroutine

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
! Builds a 2D Hamiltonian in 1D basis in phi.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function build_hamiltonian_2d(mu, rho_val, period_theta, nvec1, val1, vec1) result(ham2)
    real(real64), intent(in) :: mu, rho_val, period_theta
    integer, intent(in) :: nvec1(:) ! Number of 1D eigenpairs in each thread
    type(array_1d_real), intent(in) :: val1(:) ! 1D eigenvalues for each thread
    type(array_2d_real), intent(in) :: vec1(:) ! 1D eigenvectors for each thread
    real(real64), allocatable :: ham2(:, :)
    integer :: ham2_size, ic, ir, i
    integer, allocatable :: offset(:) ! Block offsets in 2D Hamiltonian matrix
    real(real64), allocatable :: kin(:, :), block(:, :)

    offset = prefix_sum_exclusive(nvec1)
    kin = compute_kinetic_theta(mu, rho_val, period_theta, size(nvec1))
    ham2_size = sum(nvec1)
    allocate(ham2(ham2_size, ham2_size))
    ! Loop over columns
    do ic = 1, size(nvec1)
      if (nvec1(ic) == 0) then
        cycle
      end if

      ! Loop over rows
      do ir = 1, size(nvec1)
        if (nvec1(ir) == 0) then
          cycle
        end if
        
        ! Calculate block overlaps
        if (ic == ir) then
          block = identity_matrix(nvec1(ic))
        else
          block = matmul(transpose(vec1(ir) % p), vec1(ic) % p)
        end if

        ! Factor in kinetic energy
        block = block * kin(ir, ic)
        ! Factor in 1D eigenvalues
        if (ic == ir) then
          do i = 1, nvec1(ic)
            block(i, i) = block(i, i) + val1(ic) % p(i)
          end do
        end if
        ! Write complete block
        ham2(offset(ir)+1 : offset(ir)+nvec1(ir), offset(ic)+1 : offset(ic)+nvec1(ic)) = block
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Solves 2D problem for *rho_ind*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_2d(params, rho_ind, rho_val, period_theta, nvec1, val1, vec1, val2)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind
    real(real64), intent(in) :: rho_val, period_theta
    integer, intent(in) :: nvec1(:) ! Num of 1D vecs in each thread
    type(array_1d_real), intent(in) :: val1(:) ! 1D eigenvalues for each thread
    type(array_2d_real), intent(in) :: vec1(:) ! 1D eigenvectors for each thread
    real(real64), allocatable, intent(out) :: val2(:) ! 2D values
    integer :: nvec2, file_unit
    real(real64) :: mu
    real(real64), allocatable :: val2_all(:) ! Eigenvalues
    real(real64), allocatable :: ham2(:, :), vec2(:, :)

    mu = get_reduced_mass(params % mass)
    ham2 = build_hamiltonian_2d(mu, rho_val, period_theta, nvec1, val1, vec1)
    call lapack_eigensolver(ham2, val2_all)
    nvec2 = findloc(val2_all < params % cutoff_energy, .true., dim = 1, back = .true.)
    val2 = val2_all(:nvec2)
    vec2 = ham2(:, :nvec2)

    ! Save results in binary file for 3D solution
    open(newunit = file_unit, file = get_solutions_2d_path(get_sym_path(params), rho_ind), form = 'unformatted')
    write(file_unit) nvec2, size(ham2, 1)
    if (nvec2 > 0) then
      write(file_unit) val2
      write(file_unit) vec2
    end if
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates 1D and 2D basis.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_basis(params, grid_rho, grid_theta, grid_phi)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: grid_rho(:), grid_theta(:), grid_phi(:)
    integer :: proc_id, proc_rho_ind, ierr, file_unit
    integer, allocatable :: nvec1(:), val2_counts(:)
    real(real64) :: step_theta
    real(real64), allocatable :: val2(:)
    type(array_1d_real), allocatable :: val1(:)
    type(array_2d_real), allocatable :: vec1(:)

    call assert(get_num_procs() == size(grid_rho), 'Error: number of processes has to be equal to number of points in grid_rho.dat (' // num2str(size(grid_rho)) // ')')
    proc_id = get_proc_id()
    proc_rho_ind = proc_id + 1
    step_theta = grid_theta(2) - grid_theta(1)

    call calc_1d(params, proc_rho_ind, grid_rho(proc_rho_ind), grid_theta, grid_phi, nvec1, val1, vec1)
    call calc_2d(params, proc_rho_ind, grid_rho(proc_rho_ind), step_theta, nvec1, val1, vec1, val2)

    if (proc_id == 0) then
      allocate(val2_counts(size(grid_rho)))
    end if
    call MPI_Gather(size(val2), 1, MPI_INTEGER, val2_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Write total number of 2D vectors in each slice
    if (proc_id == 0) then
      open(newunit = file_unit, file = get_block_info_path(get_sym_path(params)))
      write(file_unit, '(I0)') val2_counts
      close(file_unit)
      print *, 'Total number of 2D basis functions: ', sum(val2_counts)
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Takes a 2D solution *vec2* expressed over 1D basis and expresses it over FBR of sin/cos, using information about 1D basis given in *nvec1* and *vec1*.
! nvec1 - Number of 1D solutions in each thread (theta-slice).
! vec1 - All 1D solutions from all threads expressed over sin/cos. i-th element contains 1D solutions from i-th thread.
! vec2 - 2D solution expressed over 1D solutions. Expansion coefficients over solutions in all threads are stacked together in a single vector.
! vec2_fbr - The same 2D solution expressed in the basis of sin/cos.
! vec2_fbr consists of M blocks of length L each. Each M-block contains expansion coefficient over the same sin/cos function in different ls.
! Each element of the vector is sum over i of a_nlm^i * b_nli^j (j is index of vec2 and is fixed within this subroutine).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function transform_basis_1d_to_fbr(nvec1, vec1, vec2) result(vec2_fbr)
    integer, intent(in) :: nvec1(:)
    type(array_2d_real), intent(in) :: vec1(:)
    real(real64), intent(in) :: vec2(:) ! 2D solution expressed over 1D solutions
    real(real64), allocatable :: vec2_fbr(:) ! vec2 expressed over FBR of sin/cos
    integer :: basis_size_phi, i, j, theta_ind

    basis_size_phi = size(vec1(1) % p, 1) ! All elements of vec1 are matrices with the same number of rows (basis size phi)
    allocate(vec2_fbr(size(nvec1) * basis_size_phi))
    i = 0
    j = 0
    do theta_ind = 1, size(nvec1)
      vec2_fbr(j+1 : j+basis_size_phi) = matmul(vec1(theta_ind) % p, vec2(i+1 : i+nvec1(theta_ind)))
      i = i + nvec1(theta_ind)
      j = j + basis_size_phi
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calculate_overlap_block(params, rho_ind_row, rho_ind_col, num_points_theta, sym_path) result(overlap_block)
    class(input_params), intent(in) :: params
    integer, intent(in) :: rho_ind_row, rho_ind_col, num_points_theta
    character(*), intent(in) :: sym_path
    real(real64), allocatable :: overlap_block(:, :)
    integer :: i, j
    integer, allocatable :: num_solutions_1d_col(:), num_solutions_1d_row(:)
    real(real64), allocatable :: energies_2d_col(:), energies_2d_row(:), solution_2d_fbr_col(:), solution_2d_fbr_row(:)
    real(real64), allocatable :: solutions_2d_col(:, :), solutions_2d_row(:, :)
    type(array_1d_real), allocatable :: energies_1d_col(:), energies_1d_row(:)
    type(array_2d_real), allocatable :: solutions_1d_col(:), solutions_1d_row(:)

    ! Load bases
    call load_solutions_1D(get_solutions_1d_path(sym_path, rho_ind_col), num_points_theta, params % basis_size_phi, num_solutions_1d_col, energies_1d_col, solutions_1d_col)
    call load_solutions_2D(get_solutions_2d_path(sym_path, rho_ind_col), energies_2d_col, solutions_2d_col)
    call load_solutions_1D(get_solutions_1d_path(sym_path, rho_ind_row), num_points_theta, params % basis_size_phi, num_solutions_1d_row, energies_1d_row, solutions_1d_row)
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
  subroutine calculate_overlaps(params, num_points_theta)
    class(input_params), intent(in) :: params
    integer, intent(in) :: num_points_theta
    integer :: proc_id, num_procs, rho_ind_col, rho_ind_row, block_ind, file_unit
    integer, allocatable :: num_solutions_2d(:)
    real(real64), allocatable :: overlap_block(:, :)
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
        overlap_block = calculate_overlap_block(params, rho_ind_row, rho_ind_col, num_points_theta, sym_path)
        open(newunit = file_unit, file = get_regular_overlap_file_path(sym_path, rho_ind_row, rho_ind_col), form = 'unformatted')
        write(file_unit) overlap_block
        close(file_unit)
      end do
    end do
  end subroutine

end module
