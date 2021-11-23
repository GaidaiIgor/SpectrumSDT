module spectrumsdt_io_mod
  use spectrumsdt_io_real_mod
  use spectrumsdt_io_complex_mod
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
! Loads and rearranges 3D expansion coefficients.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_3D_expansion_coefficients(params, N, num_solutions_2d, Cs, Cs_plain)
    class(input_params), intent(in) :: params
    integer, intent(in) :: N ! Num of points along rho
    integer, intent(in) :: num_solutions_2d(:, :)
    type(array_2d_complex), allocatable, intent(out) :: Cs(:, :) ! K x N. Inner dimensions: S_Kn x S
    complex(real64), allocatable, optional, intent(out) :: Cs_plain(:, :) ! Stacked over Ks and ns (in this order) version of Cs
    integer :: total_solutions_2d, first_elem, proc_elems, sln_ind, total_Ks, start_ind, K, K_ind, K_ind_smart, n_val
    complex(real64), allocatable :: Cs_raw_col(:) ! column in Cs_raw
    complex(real64), allocatable :: Cs_raw(:, :) ! unarranged Cs
    character(:), allocatable :: exp_coeffs_path

    ! Count total number of 2D solutions
    total_solutions_2d = 0
    do K = params % K(1), params % K(2)
      K_ind_smart = get_k_ind_smart(K, params)
      total_solutions_2d = total_solutions_2d + sum(num_solutions_2d(K_ind_smart, :))
    end do

    ! Load raw matrix
    call get_proc_elem_range(params % get_num_states_total(), first_elem, proc_elems)
    allocate(Cs_raw(total_solutions_2d, proc_elems))
    do sln_ind = first_elem, first_elem + proc_elems - 1
      exp_coeffs_path = get_solution_3d_path(params, sln_ind)
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
      call get_k_attributes(K, params, K_ind = K_ind, K_ind_smart = K_ind_smart)
      do n_val = 1, N
        Cs(K_ind, n_val) = array_2d_complex(Cs_raw(start_ind : start_ind + num_solutions_2d(K_ind_smart, n_val) - 1, :))
        start_ind = start_ind + num_solutions_2d(K_ind_smart, n_val)
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
    character(:), allocatable :: state_properties_path, format_string

    ! Sequential print
    if (get_proc_id() /= 0) then
      return
    end if

    call assert(size(eigenvalues_3D, 1) == size(section_stats, 1), 'Error: sizes of input arrays have to be the same')
    call assert(size(params % wf_sections) == size(section_stats, 2), 'Error: sizes of wf_sections and section_stats are inconsistent')
    state_properties_path = get_state_properties_path(params)

    col_width = 25
    num_digits = 15
    format_string = 'G' // num2str(col_width) // '.' // num2str(num_digits)
    open(newunit = file_unit, file = state_properties_path)

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
        if (any(params % wf_sections(j) % stat == [character(100) :: 'gamma', 'A', 'B', 'C'])) then
          next_output = next_output * au_to_wn
        end if
        write(file_unit, '(' // format_string // ')', advance = 'no') next_output
      end do
      write(file_unit, *)
    end do

    close(file_unit)
  end subroutine

end module
