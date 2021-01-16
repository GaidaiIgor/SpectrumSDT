module rovib_io_mod
  use array_1d_mod
  use array_2d_mod
  use constants
  use io_utils
  use input_params_mod
  use iso_fortran_env, only: real64
  use parallel_utils
  use rovib_utils_mod
  use spectrumsdt_paths_mod
  implicit none
  
contains 

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads number of 2D states kept in each slice.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function load_basis_size_2d(basis_size_info_path) result(basis_size_2d)
    character(*), intent(in) :: basis_size_info_path
    integer, allocatable :: basis_size_2d(:)
    integer, allocatable :: basis_size_matrix(:, :)

    basis_size_matrix = read_matrix_integer(basis_size_info_path)
    basis_size_2d = basis_size_matrix(:, 1)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads 1D solutions.
! solutions_1d_path - path to the file with eigenvalues and eigenvectors of 1D solutions for specific rho index (n) and K.
! theta_size - number of points on theta grid.
! basis_size_1d - number of 1D elemental basis functions (sin/cos) used to express 1D solutions (params % basis_size_phi).
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_solutions_1D(solutions_1d_path, theta_size, basis_size_1d, num_solutions_1d, energies_1d, exp_coeffs_1d)
    character(*), intent(in) :: solutions_1d_path
    integer, intent(in) :: theta_size, basis_size_1d
    integer, allocatable, intent(out) :: num_solutions_1d(:) ! Num of 1D solutions (S_Knl) in each thread (l)
    type(array_1d_real), allocatable, intent(out) :: energies_1d(:) ! 1D solutions eigenvalues for each thread (l)
    type(array_2d_real), allocatable, intent(out) :: exp_coeffs_1d(:) ! 1D solutions expansion coefficients (over sin/cos) for each thread (l)
    integer :: i, file_unit

    allocate(num_solutions_1d(theta_size), energies_1d(theta_size), exp_coeffs_1d(theta_size))
    open(newunit = file_unit, file = solutions_1d_path, form = 'unformatted')
    read(file_unit) num_solutions_1d
    do i = 1, theta_size
      if (num_solutions_1d(i) == 0) then
        cycle
      end if
      allocate(energies_1d(i) % p(num_solutions_1d(i)), exp_coeffs_1d(i) % p(basis_size_1d, num_solutions_1d(i)))
      read(file_unit) energies_1d(i) % p
      read(file_unit) exp_coeffs_1d(i) % p
    end do
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads 2D solutions.
! solutions_2d_path - path to the file with 2D solutions for specific K and n.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_solutions_2D(solutions_2d_path, energies_2d, exp_coeffs_2d)
    character(*), intent(in) :: solutions_2d_path
    real(real64), allocatable, intent(out) :: energies_2d(:)
    real(real64), allocatable, intent(out) :: exp_coeffs_2d(:, :) ! 2D solutions expansion coefficients over 1D solutions
    integer :: file_unit, num_solutions, exp_size ! exp_size - number of coefficients in the basis expansion over 1D solutions

    open(newunit = file_unit, file = solutions_2d_path, form = 'unformatted')
    read(file_unit) num_solutions, exp_size
    if (num_solutions == 0) then
      return ! Exit right away, if empty
    end if
    allocate(energies_2d(num_solutions), exp_coeffs_2d(exp_size, num_solutions))
    read(file_unit) energies_2d
    read(file_unit) exp_coeffs_2d
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads a 3D solution.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_solution_3D(solution_3d_path, exp_coeffs_size, exp_coeffs_3d)
    character(*), intent(in) :: solution_3d_path ! a path to the file with expansion coefficients
    integer, intent(in) :: exp_coeffs_size ! size of the vector (total number of 2D solutions) needs to be provided from the calling code
    complex(real64), allocatable, intent(out) :: exp_coeffs_3d(:) ! 3D solutions expansion coefficients over 2D solutions
    integer :: file_unit

    allocate(exp_coeffs_3d(exp_coeffs_size))
    open(newunit = file_unit, file = solution_3d_path, form = 'unformatted')
    read(file_unit) exp_coeffs_3d
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads 1D expansion coefficients of a specific K and sym stored at *sym_path* into *K_ind* of *num_solutions_1d* and *As*.
! N, L, M - basis size in rho, theta and phi.
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
! Loads and rearranges 1D expansion coefficients.
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
! Loads and rearranges 1D expansion coefficients for fixed_basis is enabled
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
      sym_path = get_sym_path_root(params % fixed_basis % root_path, params % fixed_basis % K, K_sym)
      call load_1D_expansion_coefficients_K(sym_path, K_ind, N, L, params % basis_size_phi, num_solutions_1d, As)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 2D expansion coefficients for a given K at *sym_path* into *K_ind* of corresponding arrays.
! N, L - number of points along rho and theta.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_2D_expansion_coefficients_K(sym_path, K_ind, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: K_ind, N, L
    integer, intent(in) :: num_solutions_1d(:, :, :)
    type(array_2d_real), intent(inout) :: Bs(:, :, :)
    integer, intent(inout) :: num_solutions_2d(:, :)
    type(array_2d_real), allocatable, optional, intent(inout) :: Bs_plain(:, :) ! all Ls are stacked together
    integer :: n_val, l_val, next_start
    real(real64), allocatable :: energies_2d(:)
    real(real64), allocatable :: Bs_raw(:, :) ! Unarranged matrix of expansion coefficients. Single element of Bs_plain.
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
! Loads and rearranges 2D expansion coefficients.
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
    num_solutions_2d = 0
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
! Loads and rearranges 2D expansion coefficients when fixed_basis is enabled.
! Bs: Outer dimensions: 2 x N x L. Inner dimensions: S_Knl x S_Kn. All coeffs of a single 2D solution should be collected over different Ls.
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
    num_solutions_2d = 0
    if (present(Bs_plain)) then
      allocate(Bs_plain(2, N))
    end if

    do K_sym = 0, 1
      K_ind = K_sym + 1
      sym_path = get_sym_path_root(params % fixed_basis % root_path, params % fixed_basis % K, K_sym)
      call load_2D_expansion_coefficients_K(sym_path, K_ind, N, L, num_solutions_1d, Bs, num_solutions_2d, Bs_plain)
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads and rearranges 3D expansion coefficients.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_3D_expansion_coefficients(params, N, num_solutions_2d, Cs, Cs_plain)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N ! Num of points along rho
    integer, intent(in) :: num_solutions_2d(:, :)
    type(array_2d_complex), allocatable, intent(out) :: Cs(:, :) ! K x N. Inner dimensions: S_Kn x S
    complex(real64), allocatable, optional, intent(out) :: Cs_plain(:, :) ! Stacked over Ks and ns (in this order) version of Cs
    integer :: total_solutions_2d, first_elem, proc_elems, sln_ind, total_Ks, start_ind, K, K_ind, K_sym, K_ind_comp, n_val
    complex(real64), allocatable :: Cs_raw_col(:) ! column in Cs_raw
    complex(real64), allocatable :: Cs_raw(:, :) ! unarranged Cs
    character(:), allocatable :: sym_path, exp_coeffs_path

    ! Count total number of 2D solutions
    total_solutions_2d = 0
    do K = params % K(1), params % K(2)
      call get_k_attributes(K, params, K_ind, K_sym, K_ind_comp)
      total_solutions_2d = total_solutions_2d + sum(num_solutions_2d(K_ind_comp, :))
    end do

    ! Load raw matrix
    call get_proc_elem_range(params % eigencalc % num_states, first_elem, proc_elems)
    allocate(Cs_raw(total_solutions_2d, proc_elems))
    sym_path = get_sym_path(params)
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
! Writes state properties to file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_state_properties(params, eigenvalues_3D, section_stats)
    class(input_params), intent(in) :: params
    real(real64), intent(in) :: eigenvalues_3D(:, :), section_stats(:, :)
    integer :: file_unit, i, j, col_width, num_digits
    real(real64) :: next_output
    character(:), allocatable :: sym_path, properties_result_path, format_string

    ! Sequential print
    if (get_proc_id() /= 0) then
      return
    end if

    call assert(size(eigenvalues_3D, 1) == size(section_stats, 1), 'Error: sizes of input arrays have to be the same')
    call assert(size(params % wf_sections) == size(section_stats, 2), 'Error: sizes of wf_sections and section_stats are inconsistent')
    sym_path = get_sym_path(params)
    properties_result_path = get_properties_result_path(sym_path)

    col_width = 25
    num_digits = 15
    format_string = 'G' // num2str(col_width) // '.' // num2str(num_digits)
    open(newunit = file_unit, file = properties_result_path)

    ! Write header
    write(file_unit, '(2A)', advance = 'no') align_center('Energy (cm^-1)', col_width), align_center('Total Gamma (cm^-1)', col_width)
    do j = 1, size(params % wf_sections)
      write(file_unit, '(A)', advance = 'no') align_center(params % wf_sections(j) % name, col_width)
    end do
    write(file_unit, *)

    ! Write values
    do i = 1, size(eigenvalues_3D, 1)
      write(file_unit, '(2' // format_string // ')', advance = 'no') eigenvalues_3D(i, :)
      do j = 1, size(params % wf_sections)
        next_output = section_stats(i, j)
        if (params % wf_sections(j) % stat == 'gamma') then
          next_output = next_output * au_to_wn
        end if
        write(file_unit, '(' // format_string // ')', advance = 'no') next_output
      end do
      write(file_unit, *)
    end do

    close(file_unit)
  end subroutine

end module
