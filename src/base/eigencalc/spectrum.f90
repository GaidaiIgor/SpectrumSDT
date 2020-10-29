!-------------------------------------------------------------------------------------------------------------------------------------------
! High level procedures directly related to rovibrational spectrum manipulations.
!-------------------------------------------------------------------------------------------------------------------------------------------
module spectrum_mod
  use cap_mod, only: get_complex_cap
  use constants, only: au_to_wn
  use general_utils
  use general_vars, only: mu, n1, alpha1, jac1
  use input_params_mod
  use matmul_operator_mod, only: rovib_ham, init_matmul
  use parallel_utils
  use path_utils
  use sdt, only: compute_kinetic_energy_dvr
  use slepc_solver_mod, only: find_eigenpairs_slepc
  use spectrumsdt_paths_mod
  implicit none

contains
   
!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints spectrum
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_spectrum(params, eivals, eivecs)
    class(input_params), intent(in) :: params
    complex(real64), allocatable, intent(in) :: eivals(:)
    complex(real64), allocatable, intent(in) :: eivecs(:, :)
    integer :: proc_first_state, proc_states, file_unit, i, global_state_ind
    real(real64) :: energy, gamma
    character(:), allocatable :: sym_path, file_path

    call get_proc_elem_range(size(eivals), proc_first_state, proc_states)
    sym_path = get_sym_path(params)

    ! Write each eigenvector in a separate binary file
    do i = 1, proc_states
      ! Get global state number
      global_state_ind = proc_first_state + i - 1
      file_path = get_solution_3d_path(sym_path, global_state_ind)

      open(newunit = file_unit, file = file_path, form = 'unformatted')
      write(file_unit) eivecs(:, i)
      close(file_unit)
    end do

    ! Final file is written by the 0th proc
    if (get_proc_id() /= 0) return
    ! Write spectrum
    file_path = get_spectrum_path(sym_path)

    open(newunit = file_unit, file = file_path)
    do i = 1, size(eivals)
      energy = real(eivals(i)) * au_to_wn
      gamma = aimag(eivals(i)) * au_to_wn * (-2)
      write(file_unit, '(I5,2G25.15)') i, energy, gamma
    end do
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets up rovib_ham and computes (rotational-)vibrational states
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_states(params)
    type(input_params), intent(in) :: params
    complex(real64), allocatable :: eivals(:), cap(:)
    complex(real64), allocatable :: eivecs(:, :), kinetic(:, :)

    call print_parallel('Using rovib coupling')
    rovib_ham % compression = merge(1, 0, params % optimized_mult == 1) ! Global in matmul_operator_mod
    if (rovib_ham % compression == 0) then
      call print_parallel('Warning: using uncompressed Hamiltonian matrix')
    end if

    if (params % use_optimized_grid_rho == 1) then
      kinetic = compute_kinetic_energy_dvr(mu, n1, n1 * alpha1, jac1)
    else
      kinetic = compute_kinetic_energy_dvr(mu, n1, n1 * alpha1)
    end if

    if (params % cap_type /= 'none') then
      cap = get_complex_cap()
      call rovib_ham % build(params, kinetic, cap)
    else
      call rovib_ham % build(params, kinetic)
    end if
    call print_parallel('Hamiltonian matrix is built. Size: ' // num2str(rovib_ham % global_chunk_info % columns) // ' x ' // num2str(size(rovib_ham % proc_chunk, 2)))

    call init_matmul(params) ! rovib_ham is already set
    call print_parallel('Eigenvalue solver has started')
    call find_eigenpairs_slepc(params % num_states, params % ncv, params % mpd, eivals, eivecs)
    call print_spectrum(params, eivals, eivecs)
  end subroutine

end module
