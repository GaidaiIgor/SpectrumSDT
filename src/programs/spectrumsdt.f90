!-----------------------------------------------------------------------
!  Calculation of rovibrational states of ozone in APH coordinates.
!  Uses Sequential Diagonalization Truncation approach.
!  Authors: Alexander Teplukhin, Igor Gayday
!  To-do: asym (1, 1), asym move sym
!  Debug:
!-----------------------------------------------------------------------
program spectrumsdt
  use arnoldi_operator_mod
  use config
  use constants
  use debug_tools
  use dictionary
  use dict_utils
  use distributed_rovib_hamiltonian_mod
  use general_vars
  use index_conversion
  use input_params_mod
  use io_utils
  use mpi_f08
  use overlaps_extra_mod
  use parallel_utils
  use path_utils
  use pesgeneral
  use sdt
  use slepc_solver_mod
  use state_properties_mod

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
    trnchan = 160
    nvec2prt = 160
    nvec2sch = 160
    nch = 10
    ncv = params % basis_size_arnoldi
    maxitr = params % max_iterations
    bst1 = 1
    bstn = 5

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

    solver   = 1
    realver  = .false.
    onewell  = .false.
    dvr      = .false.
    if (params % symmetry == 0) then
      sy = 5
    else if (params % symmetry == 1) then
      sy = 6
    end if
    basout   = 0
    trmeth   = 0
    oltype   = 0
    recstart = 0
    recld    = 0
    ham1type = 3 ! complex fft
    ham2type = 1 ! analytical

    ! Set adiabatic control variable
    adiab = oltype == OVERLAP_ALL ! OVERLAP_ALL=0

    ! Paths
    gpath = params % grid_path
    sym_path = get_sym_path_params(params)
    bpath = get_basis_path(sym_path)
    opath = get_overlaps_path(sym_path)
    dpath = get_diagonalization_path(sym_path)
    rpath = get_channels_folder_parent_path(params % channels_root, params % J, params % K(1), params % symmetry)

    ! Set output directory
    if (.not. (mode == 0)) then 
      outdir = getdir(mode)
    end if

    ! Load grids
    call load_grids

    ! Set DVR basis size equal to number of points
    if (dvr) n3b = n3

    ! Setup shortcuts for products
    nn   = n1 * n2 * n3
    n12  = n1 * n2
    n23  = n2 * n3
    n23b = n2 * n3b

    ! Set ranges of basis sizes
    nvec1min = 3
    if (sy == SY_NONE) nvec1min = nvec1min * 2
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
    implicit none
    class(input_params), intent(in) :: params
    integer i1, i2, i3

    if(mode /= MODE_DPROD .and. mode /= MODE_BASIS) return

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
!  Calculates N second derivatives for N "delta" functions.
!  The function #i is equal 1 in point #i and 0 in all other points.
!-----------------------------------------------------------------------
  subroutine calc_der(der,n,freq,jac,gt)
    implicit none
    real*8 der(n,n),freq(n),jac(n)
    integer i,n,gt
    
    der = 0d0
    do i=1,n
      der(i,i) = 1d0
    end do
    if(gt==1)then
      call calc_derivd_jac_2nd(n,n,der,freq,jac)
    else
      call calc_derivd(2,n,n,freq,der)
    end if
  end subroutine
  
!-----------------------------------------------------------------------
!  Calculates kinetic energy matrices
!-----------------------------------------------------------------------
  subroutine calc_kin
    implicit none
    integer i1,i2
    
    allocate(freq1(n1),freq2(n2),freq3(n3), der1(n1,n1),der2(n2,n2),der3(n3,n3), sintet2(n2),grho2(n1))
    call init_derivd(1,n1,n1*alpha1,freq1) ! alpha1 is not multiplied by jac, but alpha2(3) were previously multiplied
    call init_derivd(2,n2,n2*alpha2,freq2)
    call init_derivd(2,n3,n3*alpha3,freq3)
    call calc_der(der1,n1,freq1,jac1,1) ! calculates 2nd derivative, last argument specifies method, 1 is for optimal grid, 2 for equidistant
    call calc_der(der2,n2,freq2,jac2,2)
    call calc_der(der3,n3,freq3,jac3,2)
    do i1=1,n1
      grho2(i1) = g1(i1)**2
    end do
    do i2=1,n2
      sintet2(i2) = sin(g2(i2))**2
    end do

    ! For complex FFT
    allocate(freq1z(n1),der1z(n1,n1))
    call init_derivz(n1,n1*alpha1,freq1z)
    der1z = 0d0
    do i1=1,n1
      der1z(i1,i1) = 1d0
    end do
    call calc_derivz_jac_2nd(n1,n1,der1z,freq1z,jac1)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes kinetic energy matrix of given size (number of points in the grid)
!-------------------------------------------------------------------------------------------------------------------------------------------
  function compute_kinetic_energy_matrix(matrix_size) result(matrix)
    integer, intent(in) :: matrix_size
    complex*16, allocatable :: matrix(:, :)

    allocate(matrix(matrix_size, matrix_size))
    call init_matrix1z(matrix)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns complex absorbing potential specified by global variable capid
!-------------------------------------------------------------------------------------------------------------------------------------------
  function return_cap() result(cap)
    complex*16, allocatable :: cap(:)
    real*8, allocatable :: imag_part(:)
    allocate(cap(n1), imag_part(n1))

    call get_cap(imag_part)
    cap = (0, -1d0) * imag_part
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Sets global variables in arnoldi_operator module
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine set_arnoldi_variables(params)
    class(input_params), intent(in) :: params
    
    rog = rovib_ham % global_chunk_info % subblocks(1, 1) % borders % top - 1
    msize = rovib_ham % global_chunk_info % columns
    if (params % rovib_coupling == 0) then
      hamz = rovib_ham % proc_chunk
    end if
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Prints spectrum
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine print_spectrum(params, eivals, eivecs, global_mapping)
    class(input_params), intent(in) :: params
    complex*16, allocatable, intent(in) :: eivals(:)
    complex*16, allocatable, intent(in) :: eivecs(:, :)
    integer, optional, intent(in) :: global_mapping ! regulates local-global index mapping behavior
    integer :: proc_id, n_procs, proc_first_state, proc_states, file_unit, i, global_state_ind, global_mapping_act
    real*8 :: energy, gamma
    character(:), allocatable :: sym_path, file_folder, file_path
    external :: numroc
    integer :: numroc

    global_mapping_act = arg_or_default(global_mapping, 0)
    call get_proc_info(proc_id, n_procs)
    if (global_mapping_act == 0) then
      proc_states = numroc(size(eivals), 1, proc_id, 0, n_procs)
    else if (global_mapping_act == 1) then
      call get_proc_elem_range(size(eivals), proc_first_state, proc_states)
    end if

    sym_path = get_sym_path_params(params)
    file_folder = get_expansion_coefficients_3d_path(sym_path)
    call create_path(file_folder)

    ! Write each eigenvector in a separate binary file
    do i = 1, proc_states
      ! Get global state number
      if (global_mapping_act == 0) then
        call l2g(i, proc_id, size(eivals), n_procs, 1, global_state_ind)
      else if (global_mapping_act == 1) then
        global_state_ind = proc_first_state + i - 1
      end if
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
      gamma = aimag(eivals(i)) * autown * -2
      write(file_unit, '(I4,2F30.17)') i, energy, gamma
    end do
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes vibrational states with rovib coupling included
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_states_rovib(params, context)
    type(input_params), intent(in) :: params
    integer, intent(in) :: context
    complex*16, allocatable :: eivals(:), cap(:)
    complex*16, allocatable :: eivecs(:, :), kinetic(:, :)

    call print_parallel('Using rovib coupling')
    kinetic = compute_kinetic_energy_matrix(n1)
    cap = return_cap()
    rovib_ham % compression = merge(1, 0, params % optimized_mult == 1)
    if (rovib_ham % compression == 0) then
      call print_parallel('Warning: using uncompressed Hamiltonian matrix')
    end if
    call rovib_ham % build(params, kinetic, cap)
    call set_arnoldi_variables(params)

    if (test_mode == 'slepc') then
      call find_eigenpairs_slepc(params % num_states, params % basis_size_arnoldi, eivals, eivecs)
      call print_spectrum(params, eivals, eivecs, 1)
    else
      if (params % optimized_mult == 1) then
        call rovib_ham % diagonalize(context, ham_mult_compressed, params % num_states, params % basis_size_arnoldi, params % max_iterations, eivals, eivecs)
      else
        call rovib_ham % diagonalize(context, ham_mult_old, params % num_states, params % basis_size_arnoldi, params % max_iterations, eivals, eivecs)
      end if
      call print_spectrum(params, eivals, eivecs)
    end if
  end subroutine

!-----------------------------------------------------------------------
!  Solution of 3D problem using SDT.
!-----------------------------------------------------------------------
  subroutine calc_3dsdt(params, context)
    type(input_params), intent(in) :: params
    integer, optional, intent(in) :: context ! for parallel solvers

    ! Initialize CAP for complex version
    if (.not.realver) call init_caps(0d0)
    
    ! Log calculation type and run
    if (myid == 0) then
      if (adiab) then
        write(*,*) 'Adiabatic calculation'
      else
        write(*,*) 'Pathway S calculation'
      end if
    end if

    ! Set procedure pointer
    select case(solver)
    case(SOLVER_SCAL)
      call calc_3dsdt_scal
      
    case(SOLVER_PARP)
      if (params % mode == 'diagonalization' .and. params % rovib_coupling == 1) then
        call calculate_states_rovib(params, context)
      else
        call calc_3dsdt_parp
      end if
      
    case default
      call calc_3dsdt_lapa
    end select
  end subroutine
  
!-----------------------------------------------------------------------
!  Calculation.
!-----------------------------------------------------------------------
  subroutine calc_sdt(params)
    type(input_params), intent(in) :: params
    character(256) fn
    integer ready

    ! Get BLACS info
    call get_proc_info(myid, nprocs)
    ! By default process grid is a row
    nprow = 1
    npcol = nprocs

    ! For ScaLapack the grid must be rectangular
    if (mode == MODE_3DSDT .and. solver == SOLVER_SCAL) then
      nprow = int(sqrt(real(nprocs)))
      npcol = nprocs / nprow
    end if
    nprocs_rect = nprow * npcol

    ! For basis the grid size must be equal to the number of slices
    if (mode == MODE_BASIS) then
      npcol = n1
    end if

    ! Call BLACS to setup process grid
    if (params % sequential == 1) then
      context = -1
      myrow = 1
      mycol = 1
    else
      call blacs_get(0,0,context)
      call blacs_gridinit(context,'R',nprow,npcol)
      call blacs_gridinfo(context,nprow,npcol,myrow,mycol)
    end if

    ! Write information
    if (myid == 0) then
      write(*,*)'Mode: ', params % mode
      write(*,*)'Processes:',nprocs
      ! write(*,*)'Grid: ',nprow,' x ',npcol, '=', nprocs_rect
    end if

    ! Setup directories
    if (params % mode == 'properties' .and. params % rovib_coupling == 1) then
      ready = 1
    else
      if (params % sequential == 1) then
        ready = maindir()
      else
        if (myid == 0) then
          ready = maindir()
          call igebs2d(context,'A',' ',1,1,ready,1)
        else
          call igebr2d(context,'A',' ',1,1,ready,1,0,0)
        end if
      end if

      ! Open log file and log process coordinates
      if (ready .and. myrow /= -1) then
        write(fn,'(4A,I5.5,A)')outdir,'/',logdir,'/log',myid,'.out'
        open(LG,file=fn)
        write(LG,*)'Process: ',myid,myrow,mycol
      end if
    end if

    ! If ready to start and process is in the grid, then solve
    if (ready .and. myrow /= -1) then
      ! Call appropriate subroutine
      select case(mode)
        case(MODE_DPROD)
          call calc_dprod
        
        case(MODE_BASIS)
          call calc_basis
        
        case(MODE_OVERLAP)
          call calc_overlap
          call calculate_overlaps_extra(params, mu, g1, g2)
        
        case(MODE_3DSDT)
          call calc_3dsdt(params, context)
        
        case(MODE_3DSDT_POST)
          if (params % rovib_coupling == 0) then
            call calc_3dsdt_post
          else
            call calculate_state_properties(params, size(g1), size(g2))
            ! call calculate_state_properties_debug()
          end if
        
        case(MODE_3DSDT_STATES)
          call calc_3dsdt_states

        case(MODE_CHRECOG)
          call calc_chrecog

        case(MODE_CHDIAG)
          call calc_chdiag

        case(MODE_CHDIAG_STATES)
          call calc_chdiag_states

        case(MODE_BASPROBS)
          call calc_basprobs
        
        case default
          if (myid == 0) write(*,*)'Mode does not exist'
      end select
      
      ! Finish
      if (params % sequential == 0) then
        call blacs_gridexit(context)
      end if
      close(LG)
    end if

    ! Log work done
    if (params % sequential == 0) then
      call blacs_get(0,0,context)
      call blacs_gridinit(context,'R',1,nprocs)
      call blacs_barrier(context,'A')
      call blacs_exit(0)
    end if
    if (myid == 0) write(*,*)'Done: ',outdir
  end subroutine
end program
