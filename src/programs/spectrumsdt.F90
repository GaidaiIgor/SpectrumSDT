!-----------------------------------------------------------------------
! Calculation of rovibrational states of ozone in APH coordinates.
! Uses Sequential Diagonalization Truncation approach.
! Authors: Alexander Teplukhin, Igor Gayday
!-----------------------------------------------------------------------
program spectrumsdt
  use cap_mod
  use config_mod
  use constants
  use debug_tools
  use distributed_rovib_hamiltonian_mod
  use general_vars
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
  integer :: ierr
  type(input_params) :: params

  ! Processes input parameters and enables MPI if sequential mode is not requested
  params = process_user_settings('spectrumsdt.config')
  call init_parameters(params)
  call init_pottot(params)
  call calc_sdt(params)

  if (params % sequential == 0) then
    call MPI_Finalize(ierr)
  end if

contains

!-----------------------------------------------------------------------
!  Loads grids.
!-----------------------------------------------------------------------
  subroutine load_grids()
    integer :: i
    
    open(2,file=gpath//'/grid_rho.dat')
    read(2,*)n1,alpha1
    allocate(g1(n1),jac1(n1))
    do i=1,n1
      read(2,*)g1(i),jac1(i)
    end do
    close(2)
    
    open(2,file=gpath//'/grid_theta.dat')
    read(2,*)n2,alpha2
    allocate(g2(n2),jac2(n2))
    do i=1,n2
      read(2,*)g2(i),jac2(i)
    end do
    close(2)
    
    open(2,file=gpath//'/grid_phi.dat')
    read(2,*)n3,alpha3
    allocate(g3(n3),jac3(n3))
    do i=1,n3
      read(2,*)g3(i),jac3(i)
    end do
    close(2)
    
    ! Treat alpha2(3) as a grid step size
    alpha2 = alpha2 * jac2(1)
    alpha3 = alpha3 * jac3(1)
  end subroutine

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
    if (params % stage == 'basis') then
      mode = 1
    else if (params % stage == 'overlaps') then
      mode = 2
    else if (params % stage == 'eigencalc') then
      mode = 3
    else if (params % stage == 'properties') then
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
    dpath = get_eigencalc_path(sym_path)

    ! Set output directory
    if (.not. (mode == 0)) then
      outdir = getdir(mode)
    end if

    ! Load grids
    call load_grids()

    ! Init grid derivatives
    grho2 = g1**2
    sintet2 = sin(g2)**2

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
    integer :: file_unit, iostat, i1, i2, i3

    if (mode /= MODE_BASIS) return

    ! Load vibrational potential
    allocate(pottot(n3, n2, n1))
    open(newunit = file_unit, file = append_path_token(gpath, 'pes.out'))
    read(file_unit, *, iostat = iostat) pottot
    close(file_unit)
    call assert(.not. is_iostat_end(iostat), 'Error: size of pes.out is not sufficient for the specified grids')

    ! Add rotational and extra potentials
    do i1 = 1, n1
      do i2 = 1, n2
        do i3 = 1, n3
          pottot(i3, i2, i1) = pottot(i3, i2, i1) + calc_potrot(g1(i1), g2(i2), params % J, params % K(1)) + calc_potxtr(g1(i1), g2(i2))
        end do
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints spectrum
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_spectrum(params, eivals, eivecs)
    class(input_params), intent(in) :: params
    complex(real64), allocatable, intent(in) :: eivals(:)
    complex(real64), allocatable, intent(in) :: eivecs(:, :)
    integer :: proc_first_state, proc_states, file_unit, i, global_state_ind
    real(real64) :: energy, gamma
    character(:), allocatable :: sym_path, file_folder, file_path
    external :: numroc
    integer :: numroc

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
    if (get_proc_id() /= 0) return
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
    complex(real64), allocatable :: eivals(:), cap(:)
    complex(real64), allocatable :: eivecs(:, :), kinetic(:, :)

    call print_parallel('Using rovib coupling')
    if (params % use_optimized_grid_rho == 1) then
      kinetic = compute_kinetic_energy_dvr(mu, n1, n1 * alpha1, jac1)
    else
      kinetic = compute_kinetic_energy_dvr(mu, n1, n1 * alpha1)
    end if

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
    call find_eigenpairs_slepc(params % num_states, params % ncv, params % mpd, eivals, eivecs)
    call print_spectrum(params, eivals, eivecs)
  end subroutine

!-----------------------------------------------------------------------
!  Calculation.
!-----------------------------------------------------------------------
  subroutine calc_sdt(params)
    type(input_params), intent(in) :: params
    integer :: ready, ierr

    ! Write information
    call print_parallel('Stage: ' // params % stage)
    call print_parallel('Processes: ' // num2str(get_num_procs()))

    ! Setup directories
    if (params % stage == 'properties') then
      ready = 1
    else
      if (get_proc_id() == 0) then
        ready = maindir()
        call assert(ready == 1, 'Error: cannot create directory structure')
      end if
    end if
    if (params % sequential == 0) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end if

    ! Call appropriate subroutine
    select case(params % stage)
      case('basis')
        call calc_basis()

      case('overlaps')
        call calc_overlap()
        call calculate_overlaps_extra(params, mu, g1, g2)

      case('eigencalc')
        call init_caps(params, 0d0)
        call calculate_states(params)

      case('properties')
        if (params % cap_type /= 'none') then
          call init_caps(params, 0d0)
          call calculate_state_properties(params, g1, size(g2), get_real_cap())
        else
          call calculate_state_properties(params, g1, size(g2))
        end if
    end select

    call print_parallel('Done')
  end subroutine
end program
