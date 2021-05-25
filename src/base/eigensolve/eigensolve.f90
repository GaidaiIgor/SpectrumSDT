!-------------------------------------------------------------------------------------------------------------------------------------------
! High level procedures directly related to rovibrational spectrum manipulations.
!-------------------------------------------------------------------------------------------------------------------------------------------
module eigensolve_mod
  use basis_mod, only: compute_kinetic_energy_dvr
  use cap_mod, only: calc_complex_cap
  use constants, only: au_to_wn
  use general_utils_mod
  use grid_info_mod
  use input_params_mod
  use io_utils_mod
  use matmul_operator_mod, only: rovib_ham, init_matmul
  use parallel_utils_mod
  use path_utils_mod
  use slepc_solver_mod, only: find_eigenpairs_slepc
  use spectrumsdt_paths_mod
  use spectrumsdt_utils_mod, only: get_reduced_mass
  implicit none

contains
   
!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints spectrum.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_states(params, eivals, eivecs)
    class(input_params), intent(in) :: params
    complex(real64), allocatable, intent(in) :: eivals(:)
    complex(real64), allocatable, intent(in) :: eivecs(:, :)
    integer :: proc_first_state, proc_states, file_unit, i, global_state_ind, col_width
    real(real64) :: energy, gamma
    character(:), allocatable :: file_path

    call get_proc_elem_range(size(eivals), proc_first_state, proc_states)
    ! Write each eigenvector in a separate binary file
    do i = 1, proc_states
      ! Get global state number
      global_state_ind = proc_first_state + i - 1
      file_path = get_solution_3d_path(params, global_state_ind)

      open(newunit = file_unit, file = file_path, form = 'unformatted')
      write(file_unit) eivecs(:, i)
      close(file_unit)
    end do

    ! Final file is written by the 0th proc
    if (get_proc_id() /= 0) then
      return
    end if

    ! Write states
    file_path = get_states_path(params)
    col_width = 25
    open(newunit = file_unit, file = file_path)
    write(file_unit, '(2A)') align_center('Energy (cm^-1)', col_width), align_center('Total Gamma (cm^-1)', col_width)
    do i = 1, size(eivals)
      energy = real(eivals(i)) * au_to_wn
      gamma = aimag(eivals(i)) * au_to_wn * (-2)
      write(file_unit, '(2G' // num2str(col_width) // '.15)') energy, gamma
    end do
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets up rovib_ham and computes (rotational-)vibrational states.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_states(params, rho_info)
    type(input_params), intent(in) :: params
    class(grid_info), intent(in) :: rho_info
    complex(real64), allocatable :: eivals(:), cap(:)
    complex(real64), allocatable :: eivecs(:, :), kinetic(:, :)

    rovib_ham % compression = params % debug % optimized_mult ! Global in matmul_operator_mod
    if (rovib_ham % compression /= 1) then
      call print_parallel('Warning: using uncompressed Hamiltonian matrix')
    end if

    if (any(.not. (rho_info % jac .aeq. 1d0))) then
      kinetic = compute_kinetic_energy_dvr(get_reduced_mass(params % mass), size(rho_info % points), size(rho_info % points) * rho_info % step, rho_info % jac)
    else
      kinetic = compute_kinetic_energy_dvr(get_reduced_mass(params % mass), size(rho_info % points), size(rho_info % points) * rho_info % step)
    end if

    cap = calc_complex_cap(params, rho_info)
    if (params % cap % type /= 'none') then
      if (get_proc_id() == 0) then
        call write_array(-aimag(cap) * au_to_wn, get_cap_path(params))
        print *, 'CAP is printed'
      end if
    end if
    call rovib_ham % build(params, kinetic, cap)
    call print_parallel('Hamiltonian matrix is built. Size: ' // num2str(rovib_ham % global_chunk_info % columns) // ' x ' // num2str(size(rovib_ham % proc_chunk, 2)))

    call init_matmul(params) ! rovib_ham is already set
    call print_parallel('Eigenvalue solver has started')
    call find_eigenpairs_slepc(params % eigensolve % num_states, params % eigensolve % ncv, params % eigensolve % mpd, eivals, eivecs)
    call print_states(params, eivals, eivecs)
  end subroutine

end module
