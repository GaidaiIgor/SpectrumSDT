!-----------------------------------------------------------------------
! Calculation of rovibrational states of ozone in APH coordinates.
! Uses Sequential Diagonalization Truncation approach.
! Authors: Alexander Teplukhin, Igor Gayday
!-----------------------------------------------------------------------
program spectrumsdt
  use cap_mod
  use config
  use constants
  use debug_tools
  use distributed_rovib_hamiltonian_mod
  use fourier_transform_mod
  use general_vars
  use index_conversion_mod
  use input_params_mod
  use io_utils
  use iso_fortran_env, only : real64
  use matmul_operator_mod
  use mpi
  use overlaps_extra_mod
  use parallel_utils
  use path_utils
  use pesgeneral
  use sdt
  use slepc_solver_mod
  use state_properties_mod
  implicit none

  type(input_params) :: params
  params = process_user_settings('spectrumsdt.config')
  call init_parameters(params)
  call init_pottot(params)
  call calc_kin
  call calc_sdt(params)

contains

!-----------------------------------------------------------------------
!  Initilialization.
!-----------------------------------------------------------------------
  subroutine init_sdt(params)
    type(input_params), intent(in) :: params
    character(:), allocatable :: sym_path

    nstate = params % num_states
    n3b = params % basis_size_phi
    trecut = params % cutoff_energy

    ! Old params are stored in plain variables defined in sdt.f90
    if (params % mode == 'basis') then
      mode = 1
    else if (params % mode == 'overlaps') then
      mode = 2
    else if (params % mode == 'diagonalization') then
      mode = 3
    else if (params % mode == 'properties') then
      mode = 4
    end if

    if (params % symmetry == 0) then
      sy = 5
    else if (params % symmetry == 1) then
      sy = 6
    end if

    ! Paths
    gpath = params % grid_path
    sym_path = get_sym_path_params(params)
    bpath = get_basis_path(sym_path)
    opath = get_overlaps_path(sym_path)
    dpath = get_diagonalization_path(sym_path)

    ! Set output directory
    if (.not. (mode == 0)) then
      outdir = getdir(mode)
    end if

    ! Load grids
    call load_grids

    ! Setup shortcuts for products
    nn   = n1 * n2 * n3
    n12  = n1 * n2
    n23  = n2 * n3
    n23b = n2 * n3b

    ! Set ranges of basis sizes
    nvec1min = 3
    nvec1max = n3
    nvec2max = 400

    ! Setting global debug params (debug_tools module)
    debug_mode = params % debug_mode
    test_mode = params % test_mode
    debug_int_1 = str2int(params % debug_param_1)
    parity = params % parity
  end subroutine

!-----------------------------------------------------------------------
!  Init parameters.
!-----------------------------------------------------------------------
  subroutine init_parameters(params)
    type(input_params), intent(in) :: params

    call init_pots_general(params)
    call init_sdt(params)
  end subroutine

!-----------------------------------------------------------------------
!  Init total potential.
!-----------------------------------------------------------------------
  subroutine init_pottot(params)
    class(input_params), intent(in) :: params
    integer i1, i2, i3

    if(mode /= MODE_BASIS) return

    ! Load vibrational potential
    allocate(pottot(n3, n2, n1))
    open(1, file = append_path_token(gpath, 'potvib.dat'), form = 'unformatted')
    read(1) pottot
    close(1)

    ! Add rotational and extra potentials
    do i1 = 1, n1
      do i2 = 1, n2
        do i3 = 1, n3
          pottot(i3, i2, i1) = pottot(i3, i2, i1) + calc_potrot(g1(i1), g2(i2), params % J, params % K(1)) + calc_potxtr(g1(i1), g2(i2))
        end do
      end do
    end do
  end subroutine

!-----------------------------------------------------------------------
!  Calculates kinetic energy matrices
!-----------------------------------------------------------------------
  subroutine calc_kin
    integer :: i
    complex(real64) :: dvr_basis(n1, n1)

    allocate(sintet2(n2), grho2(n1))
    do i = 1, n1
      grho2(i) = g1(i)**2
    end do
    do i = 1, n2
      sintet2(i) = sin(g2(i))**2
    end do

    dvr_basis = 0d0
    do i = 1, n1
      dvr_basis(i, i) = 1d0
    end do
    der1z = dft_derivative2_jac(dvr_basis, jac1, n1 * alpha1)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes kinetic energy matrix of given size (number of points in the grid)
!-------------------------------------------------------------------------------------------------------------------------------------------
  function compute_kinetic_energy_matrix(matrix_size) result(matrix)
    integer, intent(in) :: matrix_size
    complex*16, allocatable :: matrix(:, :)
    matrix = -der1z / (2d0 * mu)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints spectrum
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_spectrum(params, eivals, eivecs)
    class(input_params), intent(in) :: params
    complex*16, allocatable, intent(in) :: eivals(:)
    complex*16, allocatable, intent(in) :: eivecs(:, :)
    integer :: proc_id, n_procs, proc_first_state, proc_states, file_unit, i, global_state_ind
    real*8 :: energy, gamma
    character(:), allocatable :: sym_path, file_folder, file_path
    external :: numroc
    integer :: numroc

    call get_proc_info(proc_id, n_procs)
    call get_proc_elem_range(size(eivals), proc_first_state, proc_states)

    sym_path = get_sym_path_params(params)
    file_folder = get_expansion_coefficients_3d_path(sym_path)
    call create_path(file_folder)

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
    if (proc_id /= 0) return
    ! Write spectrum
    file_path = get_spectrum_path(sym_path)

    open(newunit = file_unit, file = file_path)
    do i = 1, size(eivals)
      energy = real(eivals(i)) * autown
      gamma = aimag(eivals(i)) * autown * (-2)
      write(file_unit, '(I5,2G25.15)') i, energy, gamma
    end do
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets up rovib_ham and computes (rotational-)vibrational states
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_states(params)
    type(input_params), intent(in) :: params
    complex*16, allocatable :: eivals(:), cap(:)
    complex*16, allocatable :: eivecs(:, :), kinetic(:, :)

    call print_parallel('Using rovib coupling')
    kinetic = compute_kinetic_energy_matrix(n1)
    rovib_ham % compression = merge(1, 0, params % optimized_mult == 1) ! Global in matmul_operator_mod
    if (rovib_ham % compression == 0) then
      call print_parallel('Warning: using uncompressed Hamiltonian matrix')
    end if

    if (params % cap_type /= 'none') then
      cap = get_complex_cap()
      call rovib_ham % build(params, kinetic, cap)
    else
      call rovib_ham % build(params, kinetic)
    end if
    call print_parallel('Hamiltonian matrix is built. Size: ' // num2str(rovib_ham % global_chunk_info % columns) // ' x ' // num2str(size(rovib_ham % proc_chunk, 2)))

    call set_matmul_variables(params) ! rovib_ham is already set
    call print_parallel('Eigenvalue solver has started')
    if (params % solver == 'slepc') then
      call find_eigenpairs_slepc(params % num_states, params % ncv, params % mpd, eivals, eivecs)
      call print_spectrum(params, eivals, eivecs)
    end if
  end subroutine

!-----------------------------------------------------------------------
!  Calculation.
!-----------------------------------------------------------------------
  subroutine calc_sdt(params)
    type(input_params), intent(in) :: params
    integer :: ready, ierr

    call get_proc_info(myid, nprocs)
    ! Write information
    call print_parallel('Mode: ' // params % mode)
    call print_parallel('Processes: ' // num2str(nprocs))

    ! Setup directories
    if (params % mode == 'properties') then
      ready = 1
    else
      if (myid == 0) then
        ready = maindir()
        call assert(ready == 1, 'Error: cannot create directory structure')
      end if
    end if
    if (params % sequential == 0) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end if

    ! Call appropriate subroutine
    select case(mode)
      case(MODE_BASIS)
        call calc_basis

      case(MODE_OVERLAP)
        call calc_overlap
        call calculate_overlaps_extra(params, mu, g1, g2)

      case(MODE_3DSDT)
        call init_caps(params, 0d0)
        call calculate_states(params)

      case(MODE_3DSDT_POST)
        call init_caps(params, 0d0)
        call calculate_state_properties(params, g1, size(g2), get_real_cap())

      case default
        call print_parallel('Config check has failed. Mode does not exist.')
    end select

    call print_parallel('Done')
    call MPI_Finalize(ierr)
  end subroutine
end program
