module rovib_io_mod
  use array_1d_mod
  use array_2d_mod
  use io_utils
  
contains 

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads 1D solutions
! solutions_1d_path - path to the file with eigenvalues and eigenvectors of 1D solutions for specific n and K
! theta_size - number of points on theta grid
! basis_size_1d - number of 1D basis functions used to express 1D solutions
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_solutions_1D(solutions_1d_path, theta_size, basis_size_1d, num_solutions_1d, energies_1d, exp_coeffs_1d)
    character(*), intent(in) :: solutions_1d_path
    integer, intent(in) :: theta_size, basis_size_1d
    integer, allocatable, intent(out) :: num_solutions_1d(:) ! Num of 1D solutions (S_Knl) in each thread (l)
    type(array_1d_real), allocatable, intent(out) :: energies_1d(:) ! 1D solutions eigenvalues for each thread (l)
    type(array_2d_real), allocatable, intent(out) :: exp_coeffs_1d(:) ! 1D solutions expansion coefficients (over sin/cos) for each thread (l)
    integer :: i, file_unit
    
    allocate(num_solutions_1d(theta_size), energies_1d(theta_size), exp_coeffs_1d(theta_size))
    open(newunit=file_unit, file=solutions_1d_path, form='unformatted')
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
! Loads 2D solutions
! solutions_2d_path - path to the file with 2D solutions for specific K and n
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_solutions_2D(solutions_2d_path, energies_2d, exp_coeffs_2d)
    character(*), intent(in) :: solutions_2d_path
    real*8, allocatable, intent(out) :: energies_2d(:)
    real*8, allocatable, intent(out) :: exp_coeffs_2d(:, :) ! 2D solutions expansion coefficients over 1D solutions
    integer :: file_unit, num_solutions, exp_size ! exp_size - number of coefficients in the basis expansion over 1D solutions

    open(newunit=file_unit, file=solutions_2d_path, form='unformatted')
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
! Loads 3D energies from spec-file
! spec_path - path to the file with 3D energies
!-------------------------------------------------------------------------------------------------------------------------------------------
  function load_energies_3D(spec_path) result(energies_3d)
    character(*), intent(in) :: spec_path
    real*8, allocatable :: energies_3d(:)
    real*8, allocatable :: file_content(:, :)
    
    file_content = read_matrix_real(spec_path)
    energies_3d = file_content(:, 2)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads a 3D solution
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_solution_3D(solution_3d_path, exp_coeffs_size, exp_coeffs_3d)
    character(*), intent(in) :: solution_3d_path ! a path to the file with expansion coefficients
    integer, intent(in) :: exp_coeffs_size ! size of the vector (total number of 2D solutions) needs to be provided from the calling code
    complex*16, allocatable, intent(out) :: exp_coeffs_3d(:) ! 3D solutions expansion coefficients over 2D solutions
    integer :: file_unit
    
    allocate(exp_coeffs_3d(exp_coeffs_size))
    open(newunit=file_unit, file=solution_3d_path, form='unformatted')
    read(file_unit) exp_coeffs_3d
    close(file_unit)
  end subroutine

end module
