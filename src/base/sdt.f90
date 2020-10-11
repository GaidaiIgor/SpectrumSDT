!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures related to calculation of basis for Hamiltonian and base overlaps
!-------------------------------------------------------------------------------------------------------------------------------------------
module sdt
  use constants, only: autown, pi
  use fourier_transform_mod, only: dft_derivative2_optimized_dvr, dft_derivative2_equidistant_dvr, dft_derivative2_equidistant_dvr_analytical
  use general_vars
  use input_params_mod
  use iso_fortran_env, only: real64
  use lapack_interface_mod
  use mpi
  use parallel_utils
  use potential_mod, only: pottot
  use spectrumsdt_paths_mod
  implicit none

  ! Data types for arrays of variable size
  type array1d
    real(real64),allocatable :: a(:)
  endtype array1d

  type array2d
    real(real64),allocatable :: a(:,:)
  endtype array2d

  ! Mode constants
  integer,parameter::MODE_BASIS         = 1
  integer,parameter::MODE_OVERLAP       = 2
  integer,parameter::MODE_3DSDT         = 3
  integer,parameter::MODE_3DSDT_POST    = 4

  ! Symmetry constants
  integer,parameter::SY_SY          = 5
  integer,parameter::SY_AS          = 6

  ! Control variables
  integer mode     ! Mode
  integer sy       ! Symmetry/basis

  ! Sizes
  integer nstate   ! Number of states to calculate
  integer nn       ! Size of direct-product Hamiltonian matrix
  integer n12      ! n1 * n2
  integer n23      ! n2 * n3
  integer n3b      ! Size of 1D matrices
  integer n23b     ! n2 * n3b
  integer nvec1min ! Min number of 1D states in thread before trnc
  integer nvec1max ! Max number of 1D states in thread before trnc
  integer nvec2max ! Max number of 2D states in slice  before trnc

  ! Truncation parameters
  real(real64)  trecut      ! Cut off energy

contains

!-----------------------------------------------------------------------
!  Initilialization.
!-----------------------------------------------------------------------
  subroutine init_sdt(params)
    type(input_params), intent(in) :: params

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
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Creates output directories
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_output_directories(params)
    class(input_params), intent(in) :: params
    integer :: ierr

    if (get_proc_id() == 0) then
      if (params % stage == 'basis') then
        call create_path(get_basis_results_path(get_sym_path(params)))
      end if
      if (params % stage == 'overlaps') then
        call create_path(get_overlaps_results_path(get_sym_path(params)))
      end if
      if (params % stage == 'eigencalc') then
        call create_path(get_eigencalc_results_path(get_sym_path(params)))
      end if
    end if

    if (params % sequential == 0) then
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates FBR basis on the grid.
  !     Symmetry       Cos(nx)         Sin(nx)
  !    A1 and A2       0,3,6,9           3,6,9,12
  !    E1 and E2       1,2,4,5           1,2,4,5
  !    SY and AS       0,1,2,3           1,2,3,4
  !-----------------------------------------------------------------------
  subroutine init_fbrbasis(basis)
    real(real64) basis(n3,n3b)
    real(real64) norm
    integer i,j,np

    ! Normalization coefficient
    norm = sqrt(1/pi)

    ! Periods coefficient
    np = 3
    if (sy==SY_SY.or.sy==SY_AS)np = 1

    ! Fill in basis array
    do j=1,n3b
      do i=1,n3
        if (sy==SY_SY) then
          basis(i,j) = cos(np * (j-1) * g3(i)) * norm
          if (j==1) basis(i,j) = norm / sqrt(2d0)
        else
          basis(i,j) = sin(np * j * g3(i)) * norm
        endif
      enddo
    enddo
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates 1D Hamiltonian for coordinate #3 using FBR.
  !-----------------------------------------------------------------------
  subroutine init_matrix3fbr(ham,ic1,ic2)
    real(real64) ham(n3b,n3b)
    real(real64) basis(n3,n3b)
    integer ic1,ic2,ic3
    integer i,j,np
    real(real64) coeff,sum

    ! Initialize basis
    call init_fbrbasis(basis)
    ! Coefficient in Hamiltonian
    coeff = - 1 / (2d0 * mu) * 4 / grho2(ic1) / sintet2(ic2)

    ! Build potential energy matrix
    do i=1,n3b
      do j=1,n3b
        sum = 0
        do ic3=1,n3
          sum = sum + basis(ic3,i)*pottot(ic3,ic2,ic1)*basis(ic3,j)
        enddo
        ham(i,j) = sum * alpha3
      enddo
    enddo

    ! Periods coefficient
    np = 3
    if (sy==SY_SY.or.sy==SY_AS)np = 1

    ! Build kinetic energy matrix
    j = 0
    if (sy==SY_SY)j = 1
    do i=1,n3b
      ham(i,i) = ham(i,i) + coeff * (-1) * np**2 * (i-j)**2
    enddo
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates symmetry of 1D state.
  !-----------------------------------------------------------------------
  subroutine symmetry_1d(psi,s3)
    real(real64) psi(n3),s3
    integer i3,k1,k2
    s3 = 0
    do i3=1,n3
      k1 = i3
      k2 = n3-i3+1
      s3 = s3 + psi(k1)*psi(k2) * alpha3
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Solves 1D problems for each thread in the slice.
  !  Serial solution.
  !-----------------------------------------------------------------------
  subroutine calc_1d(params, val, vec, nvec, nb2, i1)
    class(input_params), intent(in) :: params
    type(array1d),allocatable::val(:) ! 1D values  for each thread
    type(array2d),allocatable::vec(:) ! 1D vectors for each thread
    integer,allocatable::nvec(:)      ! Num of 1D vecs in each thread
    integer,allocatable::nbr(:)       ! Nums of good vecs in vecraw
    integer nb2                       ! Size of 2D vector in basis
    integer i1                        ! Slice number
    integer i2                        ! Thread number in slice
    integer i3                        ! State number in thread
    integer ivec                      ! Good state number
    integer nnode                     ! Number of nodes
    integer :: k, file_unit
    real(real64),allocatable::valraw(:)     ! Eigenvalues
    real(real64),allocatable::vecraw(:,:)   ! Eigenvectors on the grid
    real(real64),allocatable::vecrawb(:,:)  ! Basis expansion of eivectors
    real(real64),allocatable::basis(:,:)    ! FBR basis
    real(real64),allocatable::psi(:)        ! State

    ! Deallocate solution if allocated
    call free_1d(nvec,val,vec)

    ! Allocate ararys
    allocate(vecrawb(n3b,n3b), basis(n3,n3b), psi(n3), val(n2), vec(n2), nbr(n3b), nvec(n2))
    call init_fbrbasis(basis)

    ! Solve eigenvalue problem for each thread
    nb2 = 0
    do i2=1,n2
      ! Initialize matrix
      call init_matrix3fbr(vecrawb,i1,i2)

      ! Solve matrix
      call lapack_eigensolver(vecrawb, valraw)

      ! Get normalized grid functions
      vecraw = matmul(basis, vecrawb)

      ! Find good vectors
      ivec = 0
      do i3=1,n3b
        ! Skip if higher then Ecut and
        ! minimum number of vectors is taken already
        if (ivec>=nvec1min.and.valraw(i3)>trecut)exit

        ! Calculate number of nodes
        psi = vecraw(:,i3)
        nnode = 0
        do k=2,n3
          if (psi(k-1)*psi(k)<0) nnode = nnode + 1
        end do

        ! If it is a weird vector, then skip it
        if (nnode/=n3-1) then
          ivec = ivec + 1
          nbr(ivec) = i3
        end if
      ! Loop over vectors
      end do

      ! Save good vectors
      if (ivec/=0) then
        allocate(val(i2)%a(ivec))
        allocate(vec(i2)%a(n3b,ivec))
        do i3=1,ivec
          val(i2)%a(i3)   = valraw(nbr(i3))
          vec(i2)%a(:,i3) = vecrawb(:,nbr(i3))
        end do
      end if
      nvec(i2) = ivec
      nb2 = nb2 + nvec(i2)
    ! Loop over threads
    end do

    ! Save results in binary file
    open(newunit = file_unit, file = get_solutions_1d_path(get_sym_path(params), i1), form = 'unformatted')
    write(file_unit) nvec
    do i2 = 1, n2
      if (nvec(i2) /= 0) then
        write(file_unit) val(i2) % a
        write(file_unit) vec(i2) % a
      end if
    end do
    close(file_unit)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Computes kinetic energy matrix for a particle with given reduced *mass* in DVR basis described by the remaining arguments
!-------------------------------------------------------------------------------------------------------------------------------------------
  function compute_kinetic_energy_dvr(mass, num_points, period, jac) result(matrix)
    real(real64), intent(in) :: mass
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
    matrix = -matrix / (2 * mass)
  end function

  !-----------------------------------------------------------------------
  ! Calculates 1D Hamiltonian for theta.
  !-----------------------------------------------------------------------
  function compute_hamiltonian_theta(rho_ind) result(ham)
    integer, intent(in) :: rho_ind
    real(real64), allocatable :: ham(:, :)
    complex(real64), allocatable :: ham_complex(:, :)

    ham_complex = compute_kinetic_energy_dvr(mu, n2, n2 * alpha2)
    call assert(maxval(abs(aimag(ham_complex))) < 1d-10, 'Error: unexpected imaginary components of equidistant theta DVR')
    ham = 4 / g1(rho_ind)**2 * real(ham_complex)
  end function

  !-----------------------------------------------------------------------
  !  Calculates symmetry of 2D state.
  !-----------------------------------------------------------------------
  subroutine symmetry_2d(psi,s3)
    real(real64) psi(n23),s3
    integer i2,i3,j1,j2,k1,k2
    s3 = 0
    do i2=1,n2
      do i3=1,n3
        j1 = i3
        j2 = n3-i3+1
        k1 = j1 + (i2-1)*n3
        k2 = j2 + (i2-1)*n3
        s3 = s3 + psi(k1)*psi(k2) * alpha2 * alpha3
      end do
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Transform 2D state from basis to grid. Real state.
  !-----------------------------------------------------------------------
  subroutine trans_bas_grd_d(basis,vec2bas,vec2grd)
    real(real64) basis(:,:)    ! FBR basis
    real(real64) vec2bas(:)    ! 2D vector in basis
    real(real64) vec2grd(:)    ! 2D vector on grid
    integer is,i,j

    do is=1,n2
      i = (is-1) * n3b
      j = (is-1) * n3
      vec2grd(j+1 : j+n3) = matmul(basis, vec2bas(i+1 : i+n3b))
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Solves 2D problem for the slice.
  !  Serial solution.
  !  Returns 2D vectors in basis representation to save memory.
  !-----------------------------------------------------------------------
  subroutine calc_2d(params, val2, vec2, sym2, nvec2, i1, val1, vec1, nvec1)
    class(input_params), intent(in) :: params
    type(array1d),allocatable::val1(:) ! 1D values  for each thread
    type(array2d),allocatable::vec1(:) ! 1D vectors for each thread
    real(real64),allocatable::vec2(:,:)      ! 2D vectors in basis
    real(real64),allocatable::val2(:)        ! 2D values
    real(real64),allocatable::sym2(:)        ! 2D symmetries
    integer nvec1(n2)                  ! Num of 1D vecs in each thrd
    integer nvec2                      ! Num of 2D vecs
    integer i1                         ! Slice number
    integer ic,ir                      ! Block indices
    integer nb2                        ! Size of 2D vector in basis
    integer ivec                       ! Good state number
    integer nvec1c                     ! Number of 1D vectors in ic
    integer nvec1r                     ! Number of 1D vectors in ir
    integer :: i, j, l, is, file_unit
    integer,allocatable::offset(:)     ! Offsets in final matrix
    integer,allocatable::nbr(:)        ! Nums of good vecs in vecraw
    real(real64),allocatable::kin(:,:)       ! KEO matrix
    real(real64),allocatable::ham1(:,:)      ! One block
    real(real64),allocatable::ham2(:,:)      ! Contracted ham matrix
    real(real64),allocatable::valraw(:)      ! Eigenvalues
    real(real64),allocatable::symraw(:)      ! Symmetries
    real(real64),allocatable::vecraw(:)      ! One eigenvector  on grid
    real(real64),allocatable::vecrawb(:)     ! One eigenvector  in basis
    real(real64),allocatable::vecrawe(:,:)   ! All eigenvectors in eibasis
    real(real64),allocatable::basis(:,:)     ! FBR basis
    real(real64),allocatable::w(:)           ! Temporary array of eivalues
    real(real64) s3                          ! Symmetry along third coord
    real(real64) frac                        ! Fraction in equilat config

    nb2 = sum(nvec1)
    ! Deallocate solution if allocated
    call free_2d(nvec2,val2,vec2)
    if (allocated(sym2)) deallocate(sym2)

    ! Allocate arrays
    allocate(ham2(nb2, nb2), offset(n2))

    ! Calculate offsets in final matrix
    offset(1) = 0
    do i = 2, n2
      offset(i) = offset(i-1) + nvec1(i-1)
    end do
    ! Initialize kinetic matrix
    kin = compute_hamiltonian_theta(i1)

    ! Prepare hamiltonian in basis
    ! Loop over columns
    do ic=1,n2
      nvec1c = nvec1(ic)
      if (nvec1c==0) cycle

      ! Loop over rows
      do ir=1,n2
        nvec1r = nvec1(ir)
        if (nvec1r==0) cycle
        
        ! Calculate overlap matrix
        if (ic==ir) then
          allocate(ham1(nvec1r, nvec1c))
          ham1 = 0d0
          do i=1,nvec1c
            ham1(i,i) = 1d0
          end do
        else
          ham1 = matmul(transpose(vec1(ir) % a), vec1(ic) % a)
        end if

        ! Calculate block
        ham1 = ham1 * kin(ir,ic)
        if (ic==ir) then
          do i=1,nvec1c
            ham1(i,i) = ham1(i,i) + val1(ic) % a(i)
          end do
        end if

        ! Save block
        do i=1,nvec1c
          ham2(offset(ir)+1:offset(ir)+nvec1r, offset(ic)+i) = ham1(:,i)
        end do

        ! Deallocate temporary array
        deallocate(ham1)
      ! Loop over rows
      end do
    ! Loop over columns
    end do

    ! Solve eigenvalue problem
    nvec2 = nvec2max
    if (nvec2>nb2) nvec2 = nb2
    allocate(valraw(nvec2), vecrawe(nb2, nvec2))

    call lapack_eigensolver(ham2, w)
    valraw  = w(1:nvec2)
    vecrawe = ham2(:,1:nvec2)

    ! No need in Hamiltonian matrix, so deallocate
    deallocate(ham2)
    allocate(basis(n3,n3b))
    call init_fbrbasis(basis)

    ! Process calculated 2D vectors
    allocate(vecrawb(n23b),vecraw(n23),nbr(nvec2),symraw(nvec2))
    ivec = 0

    do is=1,nvec2
      ! Transform from eigenvector basis
      call trans_eibas_bas(vec1,nvec1,vecrawe(:,is),vecrawb)
      ! Get normalized grid function
      call trans_bas_grd_d(basis,vecrawb,vecraw)
      vecraw = vecraw / sqrt(alpha2)

      ! Truncate basis
      if (valraw(is) > trecut) then
        exit
      end if

      ! Filter out eigenvectors with fraction in equilateral config
      frac = 0
      do i=1,n2/10
        do j=1,n3
          l = (i-1)*n3 + j
          frac = frac + vecraw(l)**2
        end do
      end do
      frac = frac * alpha2 * alpha3
      if (frac>0.25d0) then
        cycle
      end if

      ! Calculate and store symmetry
      call symmetry_2d(vecraw,s3)
      symraw(is) = s3

      ! Remember eigenvector number
      ivec = ivec + 1
      nbr(ivec) = is
    ! Loop over states
    end do
    nvec2 = ivec

    ! Save good states
    if (ivec/=0) then
      allocate(val2(ivec),sym2(ivec),vec2(nb2,ivec))
      do i=1,ivec
        val2(i)   = valraw(nbr(i))
        sym2(i)   = symraw(nbr(i))
        vec2(:,i) = vecrawe(:,nbr(i))
      end do
    end if

    ! Save results in binary file for 3D solution
    open(newunit = file_unit, file = get_solutions_2d_path(get_sym_path(params), i1), form = 'unformatted')
    write(file_unit) ivec, nb2
    if (ivec /= 0) then
      write(file_unit) val2
      write(file_unit) vec2
    end if
    close(file_unit)
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates 1D and 2D basis.
  !-----------------------------------------------------------------------
  subroutine calc_basis(params)
    class(input_params), intent(in) :: params
    ! Variables for one slice
    type(array1d),allocatable::val1(:) ! 1D values
    type(array2d),allocatable::vec1(:) ! 1D vectors
    real(real64),allocatable::val2(:)        ! 2D values
    real(real64),allocatable::vec2(:,:)      ! 2D vectors in basis
    real(real64),allocatable::sym2(:)        ! 2D symmeties
    real(real64),allocatable::buf(:)         ! Temporary buffer
    integer,allocatable::nvec1(:)      ! Number of 1D vectors
    integer nvec2                      ! Number of 2D vectors
    integer nb2                        ! Basis size
    integer mysl                       ! SLice number
    integer :: i, j, ierr, proc_id, file_unit

    ! Collective variables
    real(real64),allocatable::val2all(:,:)   ! 2D values
    real(real64),allocatable::sym2all(:,:)   ! 2D symmetries
    integer,allocatable::nvec1all(:,:) ! Number of 1D vectors
    integer,allocatable::nvec2all(:)   ! Number of 2D vectors

    ! Set slice number
    proc_id = get_proc_id()
    mysl = proc_id + 1

    ! Calculate 1D
    call calc_1d(params, val1, vec1, nvec1, nb2, mysl)

    ! Calculate 2D
    if (nb2 /= 0) then
      call calc_2d(params, val2, vec2, sym2, nvec2, mysl, val1, vec1, nvec1)
    end if

    ! Allocate collective arrays on root
    allocate(val2all(nvec2max,n1), sym2all(nvec2max,n1), nvec1all(n2,n1), nvec2all(n1))
    val2all  = 0
    sym2all  = 0
    nvec1all = 0
    nvec2all = 0

    ! Allocate send buffer
    allocate(buf(nvec2max))
    buf = 0

    call MPI_Gather(nvec2, 1, MPI_INTEGER, nvec2all, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Gather(nvec1, n2, MPI_INTEGER, nvec1all, n2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    buf(1:nvec2) = val2
    call MPI_Gather(buf, nvec2max, MPI_DOUBLE_PRECISION, val2all, nvec2max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    buf(1:nvec2) = sym2
    call MPI_Gather(buf, nvec2max, MPI_DOUBLE_PRECISION, sym2all, nvec2max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Only root continues with printing
    if (proc_id == 0) then
      print *, 'Basis done, writing summary'

      ! Write number of 2D vectors
      open(newunit = file_unit, file = get_block_info_path(get_sym_path(params)))
      do i = 1, n1
        write(file_unit, '(2I10)') i, nvec2all(i)
      end do
      close(file_unit)

      ! Write 2D eivalues and symmetries
      open(newunit = file_unit, file = get_2d_energies_path(get_sym_path(params)))
      write(file_unit, '(10X)', advance = 'no')
      do j = 1, n1
        write(file_unit, '(I25)', advance = 'no') j
      end do
      write(file_unit, *)
      do i = 1, nvec2max
        write(file_unit, '(I10)', advance = 'no') i
        do j = 1, n1
          write(file_unit, '(F25.17)', advance = 'no') val2all(i, j) * autown ! constants
        end do
        write(file_unit, *)
      end do
      close(file_unit)

      ! Write total number of 2D basis functions
      print *, 'nvec2tot: ', sum(nvec2all)
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Frees 1D solution.
  !-----------------------------------------------------------------------
  subroutine free_1d(nvec,val,vec)
    type(array1d),allocatable::val(:) ! 1D values  for each thread
    type(array2d),allocatable::vec(:) ! 1D vectors for each thread
    integer,allocatable::nvec(:)      ! Num of 1D vecs
    integer i
    if (allocated(nvec))deallocate(nvec)
    if (allocated(val)) then
      do i=1,n2
        if (allocated(val(i)%a))deallocate(val(i)%a)
      end do
      deallocate(val)
    end if
    if (allocated(vec)) then
      do i=1,n2
        if (allocated(vec(i)%a))deallocate(vec(i)%a)
      end do
      deallocate(vec)
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Frees 2D solution.
  !-----------------------------------------------------------------------
  subroutine free_2d(nvec,val,vec)
    real(real64),allocatable::vec(:,:)      ! 2D vectors in basis
    real(real64),allocatable::val(:)        ! 2D values
    integer nvec                      ! Num of 2D vecs
    nvec = 0
    if (allocated(val)) deallocate(val)
    if (allocated(vec)) deallocate(vec)
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads number of 2D states only.
  !-----------------------------------------------------------------------
  subroutine load_nvec2(sym_path, nvec2)
    character(*), intent(in) :: sym_path
    integer, allocatable :: nvec2(:)
    integer :: i, file_unit

    allocate(nvec2(n1))
    open(newunit = file_unit, file = get_block_info_path(sym_path))
    do i = 1, n1
      read(file_unit, '(10X,I10)') nvec2(i)
    end do
    close(file_unit)
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads 1D and 2D eigenvector basis for given slice.
  !-----------------------------------------------------------------------
  subroutine load_eibasis(sym_path, isl, nvec1, val1, vec1, val2, vec2)
    character(*), intent(in) :: sym_path
    type(array1d),allocatable::val1(:) ! 1D solutions eigenvalues for each thread
    type(array2d),allocatable::vec1(:) ! 1D solutions expansion coefficients (over sin/cos) for each thread
    real(real64),allocatable::vec2(:,:)      ! 2D solutions expansion coefficients (over 1D solutions)
    real(real64),allocatable::val2(:)        ! 2D eigenvalues
    integer,allocatable::nvec1(:)      ! Num of 1D solutions in each thread
    integer nvec2                      ! Num of 2D vecs
    integer nb2                        ! Size of 2D vector in basis
    integer isl                        ! Slice number
    integer :: i, file_unit
    
    ! Deallocate arrays if allocated
    call free_1d(nvec1,val1,vec1)
    call free_2d(nvec2,val2,vec2)

    ! Load 2D solution
    open(newunit = file_unit, file = get_solutions_2d_path(sym_path, isl), form = 'unformatted')
    read(file_unit) nvec2, nb2
    if (nvec2 == 0) return ! Exit right away, if empty
    allocate(val2(nvec2),vec2(nb2,nvec2))
    read(file_unit) val2
    read(file_unit) vec2
    close(file_unit)

    ! Load 1D solution
    allocate(nvec1(n2),val1(n2),vec1(n2))
    open(newunit = file_unit, file = get_solutions_1d_path(sym_path, isl), form = 'unformatted')
    read(file_unit) nvec1
    do i = 1, n2
      if (nvec1(i) == 0) cycle
      allocate(val1(i)%a(nvec1(i)),vec1(i)%a(n3b,nvec1(i)))
      read(file_unit) val1(i) % a
      read(file_unit) vec1(i) % a
    end do
    close(file_unit)
  end subroutine

  !-----------------------------------------------------------------------
  !  Transforms 2D state from 1D eigenvector basis to DVR or FBR basis.
  !  On exit vec2b contains a vector of expansion coefficients of a given 2D solution over sin/cos
  !  Vector has M blocks of length L each. Each M-block contains expansion coefficient over the same harmonic in different ls.
  !  Each element of the vector is sum over i of a_nlm^i * b_nli^j (j is index of vec2e and is fixed within this subroutine)
  !-----------------------------------------------------------------------
  subroutine trans_eibas_bas(vec1,nvec1,vec2e,vec2b)
    type(array2d) vec1(:)    ! i-th element is a 2D array of 1D solutions in the i-th thread given as expansion coefficients over sin/cos basis
    real(real64)        vec2b(:)   ! 2D vector in basis
    real(real64)        vec2e(:)   ! expansion coefficients over 1D solutions
    integer nvec1(n2)        ! Number of 1D vectors in each thread
    integer is,i,j

    i = 0
    do is=1,n2
      j = (is-1) * n3b
      vec2b(j+1 : j+n3b) = matmul(vec1(is) % a, vec2e(i+1 : i+nvec1(is)))
      i = i + nvec1(is)
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates overlaps and saves them on disk in binary form.
  !  Only blocks in upper triangle.
  !-----------------------------------------------------------------------
  subroutine calc_overlap(params)
    class(input_params), intent(in) :: params
    ! Arrays for 1D and 2D eigenstates
    type(array1d),allocatable::val1c(:)   ! 1D values for ic
    type(array1d),allocatable::val1r(:)   ! 1D values for ir
    type(array2d),allocatable::vec1c(:)   ! 1D vecs for ic
    type(array2d),allocatable::vec1r(:)   ! 1D vecs for ir
    real(real64),allocatable::vec2c(:,:)        ! 2D vecs in basis for ic
    real(real64),allocatable::vec2r(:,:)        ! 2D vecs in basis for ir
    real(real64),allocatable::val2c(:)          ! 2D eivalues for ic slice
    real(real64),allocatable::val2r(:)          ! 2D eivalues for ir slice

    ! 2D state arrays
    real(real64),allocatable::lambdac(:)        ! 2D state on the grid, ic
    real(real64),allocatable::lambdar(:)        ! 2D state on the grid, ir

    ! Arrays for data arrays length
    integer,allocatable::nvec1c(:)        ! Number of 1D vectors, ic
    integer,allocatable::nvec1r(:)        ! Number of 1D vectors, ir
    integer,allocatable::nvec2(:)         ! Number of 2D basis vectors

    ! Overlap
    real(real64),allocatable::olap(:,:)         ! Overlap matrix
    integer ic,ir                         ! Running block indices
    integer icmax,irmax                   ! Running block indices,max
    integer ibl                           ! Global block index

    ! Miscellaneous
    integer :: i, j, file_unit

    ! Allocate arrays for vectors
    allocate(lambdac(n23b),lambdar(n23b))

    ! Load nvec2, calculate offsets and final matrix size
    call load_nvec2(get_sym_path(params), nvec2)

    ! Setup global block index
    ibl = 0
    ! Set maximum ir
    irmax = n1

    ! Loop over n block rows
    do ir=1,irmax
      ! Skip if no basis
      if (nvec2(ir) == 0) then
        cycle
      end if

      icmax = n1
      ! Loop over n block columns
      do ic = ir+1, icmax
        ! Skip if no basis
        if (nvec2(ic) == 0) then
          cycle
        end if

        ! Update block index
        ibl = ibl + 1
        ! Assign process
        if (mod(ibl - 1, get_num_procs()) /= get_proc_id()) cycle

        ! Load bases
        call load_eibasis(get_sym_path(params), ic, nvec1c, val1c, vec1c, val2c, vec2c)
        call load_eibasis(get_sym_path(params), ir, nvec1r, val1r, vec1r, val2r, vec2r)

        ! Calculate overlap matrix, olap = vec2r * vec2c
        allocate(olap(nvec2(ir),nvec2(ic)))
        do j=1,nvec2(ic)
          call trans_eibas_bas(vec1c,nvec1c,vec2c(:,j),lambdac)
          do i=1,nvec2(ir)
            call trans_eibas_bas(vec1r,nvec1r,vec2r(:,i),lambdar)
            olap(i, j) = dot_product(lambdar, lambdac)
          end do
        end do

        ! Save block
        open(newunit = file_unit, file = get_regular_overlap_file_path(get_sym_path(params), ir, ic), form = 'unformatted')
        write(file_unit) olap
        close(file_unit)

        ! Deallocate matrix
        deallocate(olap)
      end do
    end do
  end subroutine

end module
