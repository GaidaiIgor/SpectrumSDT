!  Sequential Diagonalization Truncation, 3D (sdt.f)
!  Uses optimal grid for first coordinate and equidistant for other two.
!  Calculation is controlled by "sdtcode" variable.
!
!  1st digit: mode
!    0. Calculates 3D eigenpairs of direct product matrix
!    1. Calculates 1D and 2D basis using SDT
!    2. Calculates overlap matrix between 2D vectors
!    3. Calculates 3D eigenpairs using SDT
!    4. Post-processing
!    5. Build states, computed in 3
!    6. Channel recognition
!    7. Channel diagonalization
!    8. Build states, computed in 7
!    9. Compute pathway probs for basis functions
!
!  2nd digit: solver
!    0. ScaLapack
!    1. Parpack
!    2. Lapack
!
!  3rd digit: real or complex
!    0. Real
!    1. Complex
!
!  4th digit: number of wells
!    1. One well
!    3. Three wells
!
!  5th digit: representation
!    0. DVR
!    1. FBR
!
!  6th digit: symmetry/basis
!    0. No symmetry restriction
!    1. A1       Cos  3x
!    2. A2       Sin  3x
!    3. E1       Cos  No 3x
!    4. E2       Sin  No 3x
!    5. Sym.     Cos  All
!    6. Antisym. Sin  All
!
!  7th digit: output control for mode 1
!    0. Save basis in binary form
!    1. Print 1D solutions
!    2. Print 2D solutions
!    3. Print both 1D and 2D solutions
!
!  8th digit: Truncation method
!    0. By Ecut value
!    1. Fixed number of channels
!
!  9th digit: overlaps between
!    0. Adiabatic states in all slices
!    1. Adiabatic states in adjacent slices
!    2. States of one  diabatic channel  in all slices
!    3. States of many diabatic channels in all slices in pathway S
!
! 10th digit: recognition start
!    0. The last slice
!    1. The bottom of the well
!
! 11th digit: recognition load directory
!    0. Default
!    1. Dense grid
!
! 12th digit: hamiltonian type for coordinate #1
!    0. Computed using FFT
!    1. Computed using finite differences, accuracy order 6
!    2. Computed using finite differences, accuracy order 8
!    3. Computed using FFT, complex
!
! 13th digit: hamiltonian type for coordinate #2
!    0. Computed using FFT
!    1. Computed analytically
!
!  Author: Alexander Teplukhin
!-----------------------------------------------------------------------
module sdt
  use algorithms
  use arnoldi_operator_mod
  use blas95
  use constants
  use formulas_mod
  use general_vars
  use index_conversion
  use io_utils
  use lapack95
  use mkl
  use parabola
  use parpack
  use pesgeneral
  use scalapack

  use debug_tools

  implicit none

  ! Data types for arrays of variable size
  type array1d
    real*8,allocatable :: a(:)
  endtype array1d

  type array2d
    real*8,allocatable :: a(:,:)
  endtype array2d

  type array2z
    complex*16,allocatable :: a(:,:)
  endtype array2z

  ! Conversion factors
  real*8,parameter::by2mb = 1024 * 1024  ! Bytes to megabytes

  ! Mode constants
  integer,parameter::MODE_DPROD         = 0
  integer,parameter::MODE_BASIS         = 1
  integer,parameter::MODE_OVERLAP       = 2
  integer,parameter::MODE_3DSDT         = 3
  integer,parameter::MODE_3DSDT_POST    = 4
  integer,parameter::MODE_3DSDT_STATES  = 5
  integer,parameter::MODE_CHRECOG       = 6
  integer,parameter::MODE_CHDIAG        = 7
  integer,parameter::MODE_CHDIAG_STATES = 8
  integer,parameter::MODE_BASPROBS      = 9

  ! Solver constants
  integer,parameter::SOLVER_SCAL    = 0
  integer,parameter::SOLVER_PARP    = 1
  integer,parameter::SOLVER_LAPA    = 2

  ! Symmetry constants
  integer,parameter::SY_NONE        = 0
  integer,parameter::SY_A1          = 1
  integer,parameter::SY_A2          = 2
  integer,parameter::SY_E1          = 3
  integer,parameter::SY_E2          = 4
  integer,parameter::SY_SY          = 5
  integer,parameter::SY_AS          = 6

  ! Constants for level of basis output
  integer,parameter::BASOUT_MINIMAL = 0
  integer,parameter::BASOUT_1D      = 1
  integer,parameter::BASOUT_2D      = 2
  integer,parameter::BASOUT_BOTH    = 3

  ! Constants for truncation method
  integer,parameter::TRMETH_ECUT    = 0
  integer,parameter::TRMETH_NCHAN   = 1

  ! Constants for overlaps
  integer,parameter::OVERLAP_ALL    = 0
  integer,parameter::OVERLAP_ADJ    = 1
  integer,parameter::OVERLAP_DIA    = 2
  integer,parameter::OVERLAP_PTW    = 3

  ! Constants for recognition start
  integer,parameter::RECSTART_ASYM  = 0
  integer,parameter::RECSTART_WELL  = 1

  ! Constants for recognition load directory
  integer,parameter::RECLD_DEF      = 0
  integer,parameter::RECLD_DENSE    = 1

  ! Hamiltonian type for coordinate #1
  integer,parameter::HAM1_FFT       = 0
  integer,parameter::HAM1_ANALYTIC6 = 1
  integer,parameter::HAM1_ANALYTIC8 = 2
  integer,parameter::HAM1_FFT_CMPL  = 3

  ! Hamiltonian type for coordinate #2
  integer,parameter::HAM2_FFT       = 0
  integer,parameter::HAM2_ANALYTIC  = 1

  ! Log file unit
  integer,parameter::LG             = 9

  ! Directories
  character(:),allocatable::gpath ! Grid
  character(:),allocatable::bpath ! Basis
  character(:),allocatable::opath ! Overlap
  character(:),allocatable::rpath ! Recognition
  character(:),allocatable::dpath ! Diagonalization
  character(:),allocatable::outdir
  character(7),parameter  ::dirs(*) = (/ 'basis', 'overlap', '3dsdt', '3dsdtp', '3dsdts', 'chrecog', 'chdiag', 'chdiags', 'basprobs' /)
  character(*),parameter::logdir='logs'
  character(*),parameter::capdir='caps'
  character(*),parameter::expdir='exps'

  ! Number of groups (pathways)
  integer,parameter::ngr = 3

  ! Control variables
  integer mode     ! Mode
  integer solver   ! Solver
  logical realver  ! Real or complex version
  logical onewell  ! One well or three wells
  logical dvr      ! DVR or FBR
  integer sy       ! Symmetry/basis
  integer basout   ! Output control for mode 1
  integer trmeth   ! Truncation method
  integer oltype   ! Overlap type
  logical adiab    ! Adiabatic calculation or pathway S
  integer recstart ! Recognition start
  integer recld    ! Recognition load directory
  integer ham1type ! Type of Hamiltonian for coordinate #1
  integer ham2type ! Type of Hamiltonian for coordinate #2

  ! Sizes
  integer nstate   ! Number of states to calculate
  integer n1,n2,n3 ! Problem dimensions
  integer nn       ! Size of direct-product Hamiltonian matrix
  integer n12      ! n1 * n2
  integer n23      ! n2 * n3
  integer n3b      ! Size of 1D matrices
  integer n23b     ! n2 * n3b
  integer nvec1min ! Min number of 1D states in thread before trnc
  integer nvec1max ! Max number of 1D states in thread before trnc
  integer nvec2max ! Max number of 2D states in slice  before trnc
  integer nvec2prt ! Number of 2D states to print
  integer nvec2sch ! Number of 2D states to search through for rcg

  ! Network
  integer myid        ! Process id
  integer nprocs      ! Number of processes
  integer nprocs_rect ! Number of processes forming rectangle
  integer myrow,mycol ! Process coordinates
  integer nprow,npcol ! Process grid sizes

  ! Parpack parameters
  integer ncv         ! Number of Lanczos basis vectors
  integer maxitr      ! Maximum number of iterations

  ! Range of states to build and print
  integer bst1        ! First state
  integer bstn        ! Number of states

  ! Truncation parameters
  integer trnchan     ! Number of channels
  real*8  trecut      ! Cut off energy

  ! Channel grouping after recognition
  integer nch         ! Number of channels within group

  ! Grids
  real*8,allocatable::g1(:),g2(:),g3(:)
  real*8,allocatable::jac1(:),jac2(:),jac3(:)
  real*8 alpha1,alpha2,alpha3

contains

!-----------------------------------------------------------------------
!  Loads grids.
!-----------------------------------------------------------------------
  subroutine load_grids
    integer i
    
    open(2,file=gpath//'/grid1.dat')
    read(2,*)n1,alpha1
    allocate(g1(n1),jac1(n1))
    do i=1,n1
      read(2,*)g1(i),jac1(i)
    end do
    close(2)
    
    open(2,file=gpath//'/grid2.dat')
    read(2,*)n2,alpha2
    allocate(g2(n2),jac2(n2))
    do i=1,n2
      read(2,*)g2(i),jac2(i)
    end do
    close(2)
    
    open(2,file=gpath//'/grid3.dat')
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
  !  Creates main directory.
  !-----------------------------------------------------------------------
  integer function maindir()
    logical direxists
    inquire(directory=outdir,exist=direxists)
    if (direxists) then
      write(*,*)'Output directory already exists'
      maindir = 0
    else
      call subdir(logdir)
      if (mode == MODE_3DSDT .or. mode == MODE_CHDIAG) then
        call subdir(capdir)
        call subdir(expdir)
      elseif (mode == MODE_3DSDT_POST) then
        call subdir(capdir)
      end if
      maindir = 1
    end if
  end function

  !-----------------------------------------------------------------------
  !  Solves 1D problems for each thread in the slice.
  !  Serial solution.
  !-----------------------------------------------------------------------
  subroutine calc_1d(val,vec,nvec,nb2,i1)
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
    integer i,j,k,l
    real*8 s3                         ! Symmetry along third coord
    real*8,allocatable::valraw(:)     ! Eigenvalues
    real*8,allocatable::vecraw(:,:)   ! Eigenvectors on the grid
    real*8,allocatable::vecrawb(:,:)  ! Basis expansion of eivectors
    real*8,allocatable::basis(:,:)    ! FBR basis
    real*8,allocatable::psi(:)        ! State
    logical output1d                  ! Output level
    character(256) fn                 ! File name

    write(LG,*)'Enter calc_1d'
    ! Deallocate solution if allocated
    call free_1d(nvec,val,vec)

    ! Allocate ararys
    allocate(valraw(n3b),vecraw(n3,n3b), vecrawb(n3b,n3b),basis(n3,n3b),psi(n3), val(n2),vec(n2),nbr(n3b),nvec(n2))
    write(LG,'(A10,F10.3,A3)')'valraw:', sizeof(valraw) /by2mb,' MB'
    write(LG,'(A10,F10.3,A3)')'vecraw:', sizeof(vecraw) /by2mb,' MB'
    write(LG,'(A10,F10.3,A3)')'vecrawb:',sizeof(vecrawb)/by2mb,' MB'
    write(LG,'(A10,F10.3,A3)')'basis:' , sizeof(basis)  /by2mb,' MB'
    write(LG,'(A10,F10.3,A3)')'psi:',    sizeof(psi)    /by2mb,' MB'

    ! Load basis for FBR
    if (.not.dvr)call init_fbrbasis(basis)

    ! Solve eigenvalue problem for each thread
    nb2 = 0
    do i2=1,n2
      ! Decide output level
      if ( ((basout==BASOUT_1D) .or. (basout==BASOUT_BOTH)) .and. ((i1==12.and.i2==28) .or. (i1==82.and.i2==54)  ) ) then
        output1d = .true.
      else
        output1d = .false.
      end if

      ! DVR case
      if (dvr) then
        ! Initialize matrix
        call init_matrix3dvr(vecrawb,i1,i2)
        ! Solve matrix
        call syev(vecrawb,valraw,'V')
        ! Get normalized eigenvectors on the grid
        vecraw = vecrawb / sqrt(alpha3)
      ! FBR case
      else
        ! Initialize matrix
        call init_matrix3fbr(vecrawb,i1,i2)
        ! Print matrix
        if (output1d) then
          write(fn,'(2A,I4.4,A,I4.4,A)')outdir, '/mat1.',i1,'.',i2,'.out'
          open(1,file=fn)
          do i3=1,n3b
            do j=1,n3b
              write(1,'(F25.17)',advance='no')vecrawb(i3,j)
            end do
            write(1,*)
          end do
          close(1)
        end if

        ! Solve matrix
        call syev(vecrawb,valraw,'V')
        ! Get normalized grid functions
        call gemm(basis,vecrawb,vecraw)
      ! DVR and FBR
      end if

      ! Print solution
      if (output1d) then
        ! Print eigenvalues
        write(fn,'(2A,I4.4,A,I4.4,A)')outdir, '/val1.',i1,'.',i2,'.out'
        open(1,file=fn)
        do i3=1,n3b
          call symmetry_1d(vecraw(1,i3),s3)
          write(1,'(I10,2F25.17)')i3,valraw(i3)*autown,s3
        end do
        close(1)

        ! Print eigenvectors
        write(fn,'(2A,I4.4,A,I4.4,A)')outdir, '/vec1.',i1,'.',i2,'.out'
        open(1,file=fn)
        do i3=1,n3
          do j=1,n3b
            write(1,'(F25.17)',advance='no')vecraw(i3,j)
          end do
          write(1,*)
        end do
        close(1)

        ! Print expansion
        write(fn,'(2A,I4.4,A,I4.4,A)')outdir, '/exp1.',i1,'.',i2,'.out'
        open(1,file=fn)
        do i3=1,n3b
          do j=1,n3b
            write(1,'(F25.17)',advance='no')vecrawb(i3,j)
          end do
          write(1,*)
        end do
        close(1)
      end if

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
        else
          write(LG,*)'Skipped 1D vector: ',i1,i2,i3
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
      write(LG,*)'Thread ',i1,i2,nvec(i2)
    ! Loop over threads
    end do

    ! Save results in binary file
    write(fn,'(2A,I4.4,A)')outdir,'/bas1.',i1,'.bin.out'
    open(1,file=fn,form='unformatted')
    write(1)nvec
    do i2=1,n2
      if (nvec(i2) /= 0) then
        write(1)val(i2)%a
        write(1)vec(i2)%a
      end if
    end do
    close(1)

    ! Print sizes
    write(LG,'(A10,F10.3,A3)')'vec1:',nb2*n3*8/by2mb,' MB'
    write(LG,*)'nb2:',nb2
    write(LG,*)'Exit calc_1d'
  end subroutine

  !-----------------------------------------------------------------------
  !  Solves 2D problem for the slice.
  !  Serial solution.
  !  Returns 2D vectors in basis representation to save memory.
  !-----------------------------------------------------------------------
  subroutine calc_2d(val2,vec2,sym2,nvec2,i1,val1,vec1,nvec1)
    type(array1d),allocatable::val1(:) ! 1D values  for each thread
    type(array2d),allocatable::vec1(:) ! 1D vectors for each thread
    real*8,allocatable::vec2(:,:)      ! 2D vectors in basis
    real*8,allocatable::val2(:)        ! 2D values
    real*8,allocatable::sym2(:)        ! 2D symmetries
    integer nvec1(n2)                  ! Num of 1D vecs in each thrd
    integer nvec2                      ! Num of 2D vecs
    integer i1                         ! Slice number
    integer ic,ir                      ! Block indices
    integer nb2                        ! Size of 2D vector in basis
    integer ivec                       ! Good state number
    integer nnode                      ! Number of nodes
    integer nvec1c                     ! Number of 1D vectors in ic
    integer nvec1r                     ! Number of 1D vectors in ir
    integer i,j,k,l,is
    integer,allocatable::offset(:)     ! Offsets in final matrix
    integer,allocatable::nbr(:)        ! Nums of good vecs in vecraw
    real*8,allocatable::kin(:,:)       ! KEO matrix
    real*8,allocatable::ham1(:,:)      ! One block
    real*8,allocatable::ham2(:,:)      ! Contracted ham matrix
    real*8,allocatable::valraw(:)      ! Eigenvalues
    real*8,allocatable::symraw(:)      ! Symmetries
    real*8,allocatable::vecraw(:)      ! One eigenvector  on grid
    real*8,allocatable::vecrawb(:)     ! One eigenvector  in basis
    real*8,allocatable::vecrawe(:,:)   ! All eigenvectors in eibasis
    real*8,allocatable::basis(:,:)     ! FBR basis
    real*8,allocatable::w(:)           ! Temporary array of eivalues
    real*8 s3                          ! Symmetry along third coord
    real*8 frac                        ! Fraction in equilat config
    logical output2d                   ! Output level
    character(256) fn                  ! File name

    nb2 = sum(nvec1)
    write(LG,*)'Enter calc_2d'
    write(LG,*)'Slice ',i1
    write(LG,*)'Matrix size:',nb2
    write(LG,*)'nvec1:'
    do i=1,n2
      write(LG,'(2I5)')i,nvec1(i)
    end do

    ! Deallocate solution if allocated
    call free_2d(nvec2,val2,vec2)
    if (allocated(sym2)) deallocate(sym2)

    ! Allocate arrays
    allocate(kin(n2,n2),ham2(nb2,nb2),offset(n2))
    write(LG,'(A10,F10.3,A3)')'kin:',    sizeof(kin)    /by2mb,' MB'
    write(LG,'(A10,F10.3,A3)')'ham2:',   sizeof(ham2)   /by2mb,' MB'
    write(LG,'(A10,F10.3,A3)')'offset:', sizeof(offset) /by2mb,' MB'

    ! Calculate offsets in final matrix
    offset(1) = 0
    do i=1,n2
      if (i/=1)offset(i) = offset(i-1) + nvec1(i-1)
    end do

    ! Initialize kinetic matrix
    call init_matrix2(kin,i1)
    ! Prepare hamiltonian in basis
    ! Loop over columns
    do ic=1,n2
      nvec1c = nvec1(ic)
      if (nvec1c==0)cycle

      ! Loop over rows
      do ir=1,n2
        nvec1r = nvec1(ir)
        if (nvec1r==0)cycle
        
        ! Calculate overlap matrix
        allocate(ham1(nvec1r,nvec1c))
        if (ic==ir) then
          ham1 = 0d0
          do i=1,nvec1c
            ham1(i,i) = 1d0
          end do
        else
          call gemm(vec1(ir)%a,vec1(ic)%a,ham1,'T')
        end if

        ! Calculate block
        ham1 = ham1 * kin(ir,ic)
        if (ic==ir) then
          do i=1,nvec1c
            ham1(i,i) = ham1(i,i) + val1(ic)%a(i)
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
    if (nvec2>nb2)nvec2 = nb2
    allocate(w(nb2),valraw(nvec2),vecrawe(nb2,nvec2))
    write(LG,'(A10,F10.3,A3)')'valraw:', sizeof(valraw) /by2mb,' MB'
    write(LG,'(A10,F10.3,A3)')'vecrawe:',sizeof(vecrawe)/by2mb,' MB'
    call syev(ham2,w,'V')
    valraw  = w(1:nvec2)
    vecrawe = ham2(:,1:nvec2)

    ! No need in Hamiltonian matrix, so deallocate
    deallocate(ham2)
    ! Load basis for FBR
    if (.not.dvr) then
      allocate(basis(n3,n3b))
      call init_fbrbasis(basis)
    end if

    ! Process calculated 2D vectors
    allocate(vecrawb(n23b),vecraw(n23),nbr(nvec2),symraw(nvec2))
    write(LG,'(A10,F10.3,A3)')'vecrawb:',sizeof(vecrawb)/by2mb,' MB'
    write(LG,'(A10,F10.3,A3)')'vecraw:', sizeof(vecraw) /by2mb,' MB'
    ivec = 0
    do is=1,nvec2
      ! Decide output level
      if ( (basout==BASOUT_2D .or. basout==BASOUT_BOTH) .and. (ivec < nvec2prt) ) then
        output2d = .true.
      else
        output2d = .false.
      end if

      ! Transform from eigenvector basis
      call trans_eibas_bas(vec1,nvec1,vecrawe(:,is),vecrawb)
      ! Get normalized grid function
      call trans_bas_grd_d(basis,vecrawb,vecraw)
      vecraw = vecraw / sqrt(alpha2)

      ! Print eigenvectors
      if (output2d) then
        call writ_vec2(vecraw,i1,is,1)
      end if

      ! Truncate basis
      if (trmeth == TRMETH_ECUT) then
        if (valraw(is) > trecut) then
          exit
        end if
      else
        if (ivec == trnchan)exit
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
        write(LG,*)'Skipped 2D vector: ',i1,is
        cycle
      end if

      ! Calculate and store symmetry
      call symmetry_2d(vecraw,s3)
      symraw(is) = s3

      ! Filter by symmetry in DVR case
      if (dvr.and.sy==SY_A1.and.s3<0)cycle
      if (dvr.and.sy==SY_A2.and.s3>0)cycle

      ! Remember eigenvector number
      ivec = ivec + 1
      nbr(ivec) = is
    ! Loop over states
    end do

    ! Save good states
    if (ivec/=0) then
      allocate(val2(ivec),sym2(ivec),vec2(nb2,ivec))
      do i=1,ivec
        val2(i)   = valraw(  nbr(i))
        sym2(i)   = symraw(  nbr(i))
        vec2(:,i) = vecrawe(:,nbr(i))
      end do
    end if

    ! Save results in binary file for 3D solution
    write(fn,'(2A,I4.4,A)')outdir,'/bas2.',i1,'.bin.out'
    open(1,file=fn,form='unformatted')
    write(1)ivec,nb2
    if (ivec /=0 ) then
      write(1) val2
      write(1) vec2
    end if
    close(1)

    ! Print sizes
    write(LG,'(A10,F10.3,A3)')'vec2:', sizeof(vec2) /by2mb,' MB'
    write(LG,*)'valraw:'
    write(LG,*)valraw*autown
    nvec2 = ivec
    write(LG,*)'nvec2: ',nvec2
    write(LG,*)'Exit calc_2d'
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates 1D and 2D basis.
  !-----------------------------------------------------------------------
  subroutine calc_basis
    ! Variables for one slice
    type(array1d),allocatable::val1(:) ! 1D values
    type(array2d),allocatable::vec1(:) ! 1D vectors
    real*8,allocatable::val2(:)        ! 2D values
    real*8,allocatable::vec2(:,:)      ! 2D vectors in basis
    real*8,allocatable::sym2(:)        ! 2D symmeties
    real*8,allocatable::buf(:)         ! Temporary buffer
    integer,allocatable::nvec1(:)      ! Number of 1D vectors
    integer nvec2                      ! Number of 2D vectors
    integer nb2                        ! Basis size
    integer mysl                       ! SLice number
    integer i,j

    ! Collective variables
    real*8,allocatable::val2all(:,:)   ! 2D values
    real*8,allocatable::sym2all(:,:)   ! 2D symmetries
    integer,allocatable::nvec1all(:,:) ! Number of 1D vectors
    integer,allocatable::nvec2all(:)   ! Number of 2D vectors

    integer :: ierr
    integer, allocatable :: nvec1_test(:)
    integer, allocatable :: nvec1all_test(:, :)
    integer :: my_rank

    ! Set slice number
    mysl = myid + 1

    ! Calculate 1D
    call calc_1d(val1,vec1,nvec1,nb2,mysl)

    ! Calculate 2D
    if (nb2 /= 0) then
      call calc_2d(val2,vec2,sym2,nvec2,mysl,val1,vec1,nvec1)
    end if

    ! Allocate collective arrays on root
    ! if (myid == 0) then
    allocate(val2all(nvec2max,n1), sym2all(nvec2max,n1), nvec1all(n2,n1), nvec2all(n1))
    val2all  = 0
    sym2all  = 0
    nvec1all = 0
    nvec2all = 0
  ! end if

    ! Allocate send buffer
    allocate(buf(nvec2max))
    buf = 0

    call MPI_Gather(nvec2, 1, MPI_INTEGER, nvec2all, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)

    ! allocate(nvec1_test(130), nvec1all_test(130, 90))
    ! nvec1_test = nvec1
    ! call MPI_Gather(nvec1, n2, MPI_INTEGER, nvec1all, n2, MPI_INTEGER, 0, MPI_COMM_WORLD)

    ! print *, myid, size(nvec1)
    ! if (myid == 0) then
    !   print *, 'nvec1all size', size(nvec1all, 1), size(nvec1all, 2)
    ! end if
    call MPI_Gather(nvec1, n2, MPI_INTEGER, nvec1all, n2, MPI_INTEGER, 0, MPI_COMM_WORLD)

    buf(1:nvec2) = val2
    call MPI_Gather(buf, nvec2max, MPI_DOUBLE_PRECISION, val2all, nvec2max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)
    buf(1:nvec2) = sym2
    call MPI_Gather(buf, nvec2max, MPI_DOUBLE_PRECISION, sym2all, nvec2max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)

    ! ! Colect results old
    ! if (myid == 0) then
    !   nvec2all(  1) = nvec2
    !   nvec1all(:,1) = nvec1
    !   val2all (1:nvec2,1) = val2
    !   sym2all (1:nvec2,1) = sym2
    !   do i=2,n1
    !     call igerv2d(context,1,1,nvec2all(i),1,0,i-1)
    !     call igerv2d(context,n2,1,nvec1all(1,i),n2,0,i-1)
    !     call dgerv2d(context,nvec2max,1,val2all(1,i),nvec2max,0,i-1)
    !     call dgerv2d(context,nvec2max,1,sym2all(1,i),nvec2max,0,i-1)
    !   end do
    ! else
    !   call igesd2d(context,1,1,nvec2,1,0,0)
    !   call igesd2d(context,n2,1,nvec1,n2,0,0)
    !   if (nvec2 /= 0)buf(1:nvec2) = val2
    !   call dgesd2d(context,nvec2max,1,buf,nvec2max,0,0)
    !   if (nvec2 /= 0)buf(1:nvec2) = sym2
    !   call dgesd2d(context,nvec2max,1,buf,nvec2max,0,0)
    ! end if

    ! Only root continues with printing
    if (myid == 0) then
      print *, 'Basis done, writing summary'

      ! Write number of 1D vectors
      open(1,file=outdir//'/nvec1.dat')
      write(1,'(10X)',advance='no')
      do j=1,n1
        write(1,'(I10)',advance='no')j
      end do
      write(1,*)
      do i=1,n2
        write(1,'(I10)',advance='no')i
        do j=1,n1
          write(1,'(I10)',advance='no')nvec1all(i,j)
        end do
        write(1,*)
      end do
      close(1)

      ! Write number of 2D vectors
      open(1,file=outdir//'/nvec2.dat')
      do i=1,n1
        write(1,'(2I10)')i,nvec2all(i)
      end do
      close(1)

      ! Write 2D eivalues and symmetries
      open(1,file=outdir//'/val2.out')
      open(2,file=outdir//'/sym2.out')
      write(1,'(10X)',advance='no')
      write(2,'(10X)',advance='no')
      do j=1,n1
        write(1,'(I25)',advance='no')j
        write(2,'(I25)',advance='no')j
      end do
      write(1,*)
      write(2,*)
      do i=1,nvec2max
        write(1,'(I10)',advance='no')i
        write(2,'(I10)',advance='no')i
        do j=1,n1
          write(1,'(F25.17)',advance='no')val2all(i,j) * autown ! constants
          write(2,'(F25.17)',advance='no')sym2all(i,j)
        end do
        write(1,*)
        write(2,*)
      end do
      close(1)
      close(2)

      ! Write total number of 2D basis functions
      write(*,*)'nvec2tot: ',sum(nvec2all)
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates overlaps and saves them on disk in binary form.
  !  Only blocks in upper triangle.
  !-----------------------------------------------------------------------
  subroutine calc_overlap
    ! Arrays for 1D and 2D eigenstates
    type(array1d),allocatable::val1c(:)   ! 1D values for ic
    type(array1d),allocatable::val1r(:)   ! 1D values for ir
    type(array2d),allocatable::vec1c(:)   ! 1D vecs for ic
    type(array2d),allocatable::vec1r(:)   ! 1D vecs for ir
    real*8,allocatable::vec2c(:,:)        ! 2D vecs in basis for ic
    real*8,allocatable::vec2r(:,:)        ! 2D vecs in basis for ir
    real*8,allocatable::val2c(:)          ! 2D eivalues for ic slice
    real*8,allocatable::val2r(:)          ! 2D eivalues for ir slice

    ! 2D state arrays
    real*8,allocatable::lambdac(:)        ! 2D state on the grid, ic
    real*8,allocatable::lambdar(:)        ! 2D state on the grid, ir

    ! Arrays for data arrays length
    integer,allocatable::nvec1c(:)        ! Number of 1D vectors, ic
    integer,allocatable::nvec1r(:)        ! Number of 1D vectors, ir
    integer,allocatable::nvec2(:)         ! Number of 2D basis vectors

    ! Channel data and variables
    integer,allocatable::chs(:)           ! Channel numbers for S
    integer,allocatable::ind(:,:)         ! Channel indices
    integer,allocatable::phs(:,:)         ! Channel phases
    integer ichan                         ! Channel number
    integer nchan                         ! Number of channels
    integer ichc,ichr                     ! Channel indices
    integer nchmax                        ! Maximum nch for pathway S

    ! Overlap
    real*8,allocatable::olap(:,:)         ! Overlap matrix
    integer ic,ir                         ! Running block indices
    integer icmax,irmax                   ! Running block indices,max
    integer ibl                           ! Global block index

    ! Miscellaneous
    integer numroc                        ! Calculates num of blocks
    character(256) fn
    integer i,j

    ! Allocate arrays for vectors
    allocate(lambdac(n23b),lambdar(n23b))

    ! One-channel diabatic overlaps
    if (oltype == OVERLAP_DIA) then
      ! Load recognition
      call load_recognition(ind,phs,nchan)
      ! Get my channel number
      ichan = myid + 1
      if (ichan > nchan) then
        write(LG,*)'Nothing to do, exiting'
        return
      else
        write(LG,*)'Channel number ',ichan
      end if

      ! Calculate overlap matrix
      allocate(olap(n1,n1))
      do ic=1,n1
        call load_vec2(lambdac,ic,ind(ic,ichan),.false.)
        do ir=1,n1
          call load_vec2(lambdar,ir,ind(ir,ichan),.false.)
          olap(ir,ic) = dot(lambdar,lambdac) * phs(ic,ichan) * phs(ir,ichan)
        end do
      end do

      ! Save block
      write(fn,'(2A,I4.4,A)')outdir,'/overlap.',ichan,'.bin.out'
      open(1,file=fn,form='unformatted')
      write(1)olap
      close(1)
    ! Many-channel pathway overlaps
    else if (oltype == OVERLAP_PTW) then
      ! Load recognition
      call load_recognition(ind,phs,nchan)
      ! Load channel numbers for pathway S
      call load_pathway_s(chs,nchmax)
      if (myid == 0)write(*,*)'nchmax: ',nchmax
      ! Allocate overlap matrix
      allocate(olap(n1,n1))
      ! Setup global block index
      ibl = 0
      ! Loop over rows
      do ichr=1,nchmax
        ! Loop over columns
        do ichc=ichr,nchmax
          ! Update block index
          ibl = ibl + 1
          ! Assign process
          if (mod(ibl-1,nprocs) /= myid) cycle
          ! Log block calculation
          write(LG,'(A,3I5)')'Overlap block ',ibl,ichr,ichc
          ! Get channel numbers
          i = chs(ichr)
          j = chs(ichc)

          ! Calculate overlap matrix
          do ic=1,n1
            call load_vec2(lambdac,ic,ind(ic,i),.false.)
            do ir=1,n1
              call load_vec2(lambdar,ir,ind(ir,j),.false.)
              olap(ir,ic) = dot(lambdar,lambdac) * phs(ic,i) * phs(ir,j)
            end do
          end do

          ! Save block
          write(fn,'(2A,I4.4,A,I4.4,A)') outdir,'/overlap.',ichr,'.',ichc,'.bin.out'
          open(1,file=fn,form='unformatted')
          write(1)olap
          close(1)
        end do
      end do
    ! Any other overlap type
    else
      ! Load nvec2, calculate offsets and final matrix size
      call load_nvec2(nvec2)

      ! Setup global block index
      ibl = 0
      ! Set maximum ir
      if (oltype == OVERLAP_ALL) then
        irmax = n1
      else
        irmax = n1 - 1
      end if

      ! Loop over n block rows
      do ir=1,irmax
        ! Skip if no basis
        if (nvec2(ir) == 0) then
          if (myid == 0) write(*,*)'Skipped row:',ir
          cycle
        end if

        ! Set maximum ic
        if (oltype == OVERLAP_ALL) then
          icmax = n1
        else
          icmax = ir + 1
        end if

        ! Loop over n block columns
        do ic = ir+1, icmax
          ! Skip if no basis
          if (nvec2(ic) == 0) then
            if (myid == 0) write(*,*)'Skipped col:',ic
            cycle
          end if

          ! Update block index
          ibl = ibl + 1
          ! Assign process
          if (mod(ibl-1,nprocs) /= myid) cycle
          ! Log block calculation
          write(LG,'(A,3I5,A)')'Overlap block ',ibl,ir,ic

          ! Load bases
          call load_eibasis(ic,nvec1c,val1c,vec1c,val2c,vec2c)
          call load_eibasis(ir,nvec1r,val1r,vec1r,val2r,vec2r)

          ! Calculate overlap matrix, olap = vec2r * vec2c
          allocate(olap(nvec2(ir),nvec2(ic)))
          do j=1,nvec2(ic)
            call trans_eibas_bas(vec1c,nvec1c,vec2c(:,j),lambdac)
            do i=1,nvec2(ir)
              call trans_eibas_bas(vec1r,nvec1r,vec2r(:,i),lambdar)
              olap(i, j) = dot(lambdar, lambdac)
            end do
          end do

          ! Save block
          write(fn,'(2A,I4.4,A,I4.4,A)')outdir,'/overlap.', ir,'.',ic,'.bin.out'
          open(1,file=fn,form='unformatted')
          write(1)olap
          close(1)
          ! Deallocate matrix
          deallocate(olap)
        end do
      end do
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads overlap block from the disk. Adiabatic version.
  !-----------------------------------------------------------------------
  subroutine load_overlap(gblr,gblc,blsize,mtx,unity)
    real*8,allocatable::mtx(:,:)
    real*8,allocatable::tmp(:,:)
    integer,allocatable::blsize(:)
    integer gblr,gblc
    integer gblrt,gblct
    integer i
    character(:),allocatable::ldir
    character(256) fn
    logical unity
    logical ex

    ! Exit if no basis
    if (blsize(gblr) == 0) return
    if (blsize(gblc) == 0) return

    ! Allocate array
    if (allocated(mtx))deallocate(mtx)
    allocate(mtx(blsize(gblr),blsize(gblc)))

    ! Write unity matrix or read from the disk
    if (unity .and. gblr == gblc) then
      mtx = 0d0
      do i=1,blsize(gblc)
        mtx(i,i) = 1d0
      end do
    else
      ldir = opath // '/' // getdir(MODE_OVERLAP)
      ! only upper triangle blocks are stored on disk since the matrix is hermitian 
      ! if a block from lower triangle is requested, swap indices to load corresponding block from the upper triangle
      if (gblr < gblc) then
        gblrt = gblr
        gblct = gblc
      else
        gblrt = gblc
        gblct = gblr
      end if

      ! Get filename
      write(fn,'(2A,I4.4,A,I4.4,A)')ldir, '/overlap.',gblrt,'.',gblct,'.bin.out'
      ! Check if file exists
      inquire(file=fn,exist=ex)
      if (.not.ex) then
        write(*,*)'No overlap: ',fn
        stop
      end if
      
      ! Allocate temporary array
      allocate(tmp(blsize(gblrt),blsize(gblct)))
      ! Read data
      open(1,file=fn,form='unformatted')
      read(1)tmp
      close(1)

      ! Assign data
      if (gblr <= gblc) then
        mtx = tmp
      else
        mtx = transpose(tmp)
      end if
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads asym block for K = 1
  !-----------------------------------------------------------------------
  ! function load_asym(block_num, block_sizes) result(matrix)
  !   integer, intent(in) :: block_num
  !   integer, intent(in) :: block_sizes(:)
  !   real*8, allocatable :: matrix(:, :)
  !   character(:), allocatable :: file_path
  !
  !   ! Exit if no basis
  !   if (block_sizes(block_num) == 0) return
  !   allocate(matrix(block_sizes(block_num), block_sizes(block_num)))
  !   file_path = opath // '/' // getdir(MODE_OVERLAP) // '/asym.' // num2str(-block_num) // '.bin.out'
  !   open(1, file=file_path, form='unformatted')
  !   read(1) matrix
  !   close(1)
  ! end function

  !-----------------------------------------------------------------------
  !  Loads asym block from the disk. Adiabatic version.
  !-----------------------------------------------------------------------
  ! function load_cor(block_num, size1, size2) result(matrix)
  !   integer, intent(in) :: block_num
  !   integer, intent(in) :: size1, size2
  !   real*8, allocatable :: matrix(:, :)
  !   character(:), allocatable :: file_path
  !
  !   allocate(matrix(size1, size2))
  !   file_path = opath // '/' // getdir(MODE_OVERLAP) // '/coriolis.' // num2str(block_num) // '.bin.out'
  !   open(1, file=file_path, form='unformatted')
  !   read(1) matrix
  !   close(1)
  ! end function

  !-----------------------------------------------------------------------
  !  Loads 2D energies.
  !-----------------------------------------------------------------------
  subroutine load_val2(isl,val2)
    real*8,allocatable::val2(:)        ! 2D values, output
    integer nvec2                      ! Num of 2D vecs
    integer isl                        ! Slice number
    character(:),allocatable::ldir
    character(256) fn

    ! Allocate array
    if (allocated(val2)) deallocate(val2)
    
    ! Get load directory
    ldir = bpath // '/' // getdir(MODE_BASIS)
    ! Get filename
    write(fn,'(2A,I4.4,A)')ldir,'/bas2.',isl,'.bin.out'
    ! Read data
    open(1,file=fn,form='unformatted')
    read(1)nvec2
    if (nvec2 == 0) return ! Exit right away, if empty
    
    allocate(val2(nvec2))
    read(1)val2
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads 2D energies.
  !-----------------------------------------------------------------------
  subroutine load_val2_grouped(ich,val2)
    real*8,allocatable::val2(:)        ! 2D values
    real*8,allocatable::val2sl(:)      ! 2D values in slice
    integer isl                        ! Slice number
    integer ich                        ! Channel number
    integer nchmax                     ! Maximum nch for pathway S

    ! Grouping
    integer,allocatable :: chs(:)      ! Channel numbers for pthway S
    integer,allocatable :: ind(:)      ! Channel indices
    integer,allocatable :: phs(:)      ! Channel phases

    ! Allocate array
    if (allocated(val2))deallocate(val2)
    allocate(val2(n1))
    ! Load channel grouping
    call load_pathway_s(chs,nchmax)

    ! Load channel recognition
    if (.not.load_recognition_chan(ind,phs,chs(ich))) then
      write(*,*)'Wrong channel: ',myid,ich
      stop
    end if

    ! Fill in val2
    do isl=1,n1
      call load_val2(isl,val2sl)
      val2(isl) = val2sl(ind(isl))
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Solution of 3D problem using SDT and Lapack
  !  Variables are real (d) and/or complex (z).
  !-----------------------------------------------------------------------
  subroutine calc_3dsdt_lapa
    ! Arrays for 2D states
    real*8,allocatable::val2(:)           ! 2D eivalues
    ! Arrays for real 3D states
    real*8,allocatable::val3d(:)          ! 3D eivalues
    real*8,allocatable::vec3d(:,:)        ! 3D vecs in basis
    ! Arrays for complex 3D states
    complex*16,allocatable::val3z(:)      ! 3D eivalues
    complex*16,allocatable::vec3z(:,:)    ! 3D vecs in basis
    ! Real hamiltonian arrays
    real*8,allocatable::kin(:,:)          ! KEO matrix
    real*8,allocatable::ham1d(:,:)        ! One local block
    ! Complex hamiltonian arrays
    complex*16,allocatable::kinz(:,:)     ! KEO matrix
    complex*16,allocatable::ham1z(:,:)    ! One local block
    ! Matrix partitioning
    integer,allocatable::blsize(:)        ! Block sizes
    integer,allocatable::offset(:)        ! Block offsets
    integer nbl                           ! Number of blocks in r/c
    integer gblr,gblc                     ! Global block indices
    integer blsizemax                     ! Maximum block size
    ! Miscellaneous
    integer    kc,kr                      ! Elem indices in fnl mtrx
    real*8, allocatable   ::s(:)          ! Values to sort
    integer,allocatable   ::ind(:)        ! Sorted indices
    real*8, allocatable   ::wd(:)         ! 3D eivalues, raw
    complex*16,allocatable::wz(:)         ! 3D eivalues, raw
    complex*16,allocatable::vr(:,:)       ! 3D vecs, raw
    integer    i
    
    ! Load matrix partitioning
    call load_partitioning(nbl,blsize,blsizemax,offset,msize)
    
    ! Stop if number of states exceeds matrix size
    if (nstate > msize) then
      write(*,*)'Requested number of states > matrix size, stop'
      return
    end if

    ! Allocate real arrays
    if (realver) then
      allocate(kin(n1,n1))
      write(LG,'(A10,F10.3,A3)')'kin:',   sizeof(kin)   /by2mb,' MB'
    ! Allocate complex arrays
    else
      allocate(kinz(n1,n1))
      write(LG,'(A10,F10.3,A3)')'kinz:',  sizeof(kinz)  /by2mb,' MB'
    end if

    ! Initialize KEO matrix
    ! In complex version it also includes CAP
    if (realver) then
      call init_matrix1d(kin)
    else
      call init_matrix1z(kinz)
    end if

    ! Allocate final matrix
    if (myid == 0) write(*,*) 'Final matrix size: ',msize
    write(LG,*)'Final matrix size: ',msize
    if (realver) then
      if (allocated(hamd))deallocate(hamd)
      allocate(hamd(msize,msize))
      write(LG,'(A10,F10.3,A3)')'hamd:',sizeof(hamd)/by2mb,' MB'
    else
      if (allocated(hamz))deallocate(hamz)
      allocate(hamz(msize,msize))
      write(LG,'(A10,F10.3,A3)')'hamz:',sizeof(hamz)/by2mb,' MB'
    end if

    ! Calculate blocks of Hamiltonian matrix
    ! Loop over columns
    do gblc=1,nbl
      ! Check if column is not empty and log
      if (adiab .and. blsize(gblc) == 0)cycle
      write(LG,*)'Column: ',gblc
      ! Loop over rows
      do gblr=1,nbl
        ! Check if row is not empty and log
        if (adiab .and. blsize(gblr) == 0)cycle
        write(LG,*)'Row: ',gblr
        ! Load block
        call load_overlap(gblr,gblc,blsize,ham1d,adiab)
        ! Finish block construction and save, real version
        if (realver) then
          ! Multiply by KEO element
          if (adiab) then
            ham1d = ham1d * kin(gblr,gblc)
          else
            ham1d = ham1d * kin
          end if

          ! Add 2D energies to the diagonal of diagonal block
          if (gblc == gblr) then
            if (adiab) then
              call load_val2(gblc,val2)
            else
              call load_val2_grouped(gblc,val2)
            end if
            do i=1,blsize(gblc)
              ham1d(i,i) = ham1d(i,i) + val2(i)
            end do
          end if

          ! Save block
          kr = offset(gblr)
          kc = offset(gblc)
          hamd(kr+1:kr+blsize(gblr),kc+1:kc+blsize(gblc)) = ham1d
        ! Finish block construction and save, complex version
        else
          ! Allocate complex array
          allocate(ham1z(blsize(gblr),blsize(gblc)))
          ! Multiply by KEO element
          if (adiab) then
            ham1z = ham1d * kinz(gblr,gblc)
          else
            ham1z = ham1d * kinz
          end if

          ! Add 2D energies to the diagonal of diagonal block
          if (gblc == gblr) then
            if (adiab) then
              call load_val2(gblc,val2)
            else
              call load_val2_grouped(gblc,val2)
            end if
            do i=1,blsize(gblc)
              ham1z(i,i) = ham1z(i,i) + val2(i)
            end do
          end if

          ! Save block
          kr = offset(gblr)
          kc = offset(gblc)
          do i=1,blsize(gblc)
            hamz(kr+1:kr+blsize(gblr),kc+i) = ham1z(:,i)
          end do
          ! Deallocate temporary array
          deallocate(ham1z)
        end if
        ! Deallocate temporary array
        deallocate(ham1d)
      ! Loop over rows
      end do
    ! Loop over columns
    end do

    ! Allocate memory for values, vectors and symmetries
    if (realver) then
      allocate(val3d(nstate),wd(msize),vec3d(msize,nstate))
      write(LG,'(A10,F10.3,A3)')'val3d:', sizeof(val3d) /by2mb,' MB'
      write(LG,'(A10,F10.3,A3)')'vec3d:', sizeof(vec3d) /by2mb,' MB'
    else
      allocate(val3z(nstate),wz(msize),vec3z(msize,nstate))
      write(LG,'(A10,F10.3,A3)')'val3z:', sizeof(val3z) /by2mb,' MB'
      write(LG,'(A10,F10.3,A3)')'vec3z:', sizeof(vec3z) /by2mb,' MB'
      allocate(s(msize),ind(msize),vr(msize,msize))
    end if

    ! Solve matrix
    if (realver) then
      call syev(hamd,wd,'V')
      val3d = wd(1:nstate)
      vec3z = hamd(:,1:nstate)
    else
      call geev(hamz,wz,vr,vr)
      s = real(wz)
      s = bubble_sort(s, ind) ! algorithms
      do i=1,nstate
        val3z(i)   = wz(ind(i))
        vec3z(:,i) = vr(:,ind(i))
      end do
    end if

    ! Print spectrum
    call prnt_3dsdt(val3d,vec3d,val3z,vec3z,nstate,nprocs)
  end subroutine

  !-----------------------------------------------------------------------
  !  Solution of 3D problem using SDT and ScaLapack.
  !  Variables are real (d) and/or complex (z).
  !-----------------------------------------------------------------------
  subroutine calc_3dsdt_scal
    ! Arrays for 2D states
    real*8,allocatable::val2(:)           ! 2D eivalues
    ! Arrays for real 3D states
    real*8,allocatable::val3d(:)          ! 3D eivalues
    real*8,allocatable::vec3d(:,:)        ! 3D vecs in basis
    ! Arrays for complex 3D states
    complex*16,allocatable::val3z(:)      ! 3D eivalues
    complex*16,allocatable::vec3z(:,:)    ! 3D vecs in basis
    ! Real hamiltonian arrays
    real*8,allocatable::kin(:,:)          ! KEO matrix
    real*8,allocatable::ham1d(:,:)        ! One local block
    real*8,allocatable::ham2d(:,:,:,:)    ! All local blocks
    real*8,allocatable::hamtd(:,:,:,:)    ! Stores bcasted blockes
    ! Complex hamiltonian arrays
    complex*16,allocatable::kinz(:,:)     ! KEO matrix
    complex*16,allocatable::ham1z(:,:)    ! One local block
    complex*16,allocatable::ham2z(:,:,:,:)! All local blocks
    complex*16,allocatable::hamtz(:,:,:,:)! Stores bcasted blockes
    ! Matrix partitioning
    integer,allocatable::blsize(:)        ! Block sizes
    integer,allocatable::offset(:)        ! Block offsets
    integer nbl                           ! Number of blocks in r/c
    integer nblr,nblc                     ! Local number of blocks
    integer iblr,iblc                     ! Local  block indices
    integer gblr,gblc                     ! Global block indices
    integer blsizemax                     ! Maximum block size
    integer nblcmax                       ! nblc, max
    integer nblrmax                       ! nblr, max
    integer nblrt,nblct                   ! Local num of blocks,temp
    ! Matrix distribution for ScaLapack
    integer mblr,mblc                     ! Local num of blocks
    integer fblc,fblr                     ! Local indices in fnl mtrx
    integer pcol,prow                     ! Process coords
    integer fpcol,fprow                   ! Process coords, fnl mtrx
    integer jc,jr                         ! Elem indices in block
    integer kc,kr                         ! Elem indices in fnl mtrx
    integer nb                            ! ScaLapack block size
    integer numroc                        ! Calculates num of blocks
    integer nstloc                        ! Local num of states
    ! Miscellaneous
    character(256) fn
    integer i

    ! Load matrix partitioning
    call load_partitioning(nbl,blsize,blsizemax,offset,msize)
    
    ! Stop if number of states exceeds matrix size
    if (nstate > msize) then
      if (myid == 0) then
        write(*,*)'Requested number of states > matrix size, stop'
      end if
      return
    end if

    ! Get number of blocks stored locally
    nblr = numroc(nbl,1,myrow,0,nprow)
    nblc = numroc(nbl,1,mycol,0,npcol)
    write(LG,*)'Number of blocks: ',nblr,'x',nblc

    ! Allocate real arrays
    if (realver) then
      allocate(kin(n1,n1),ham2d(blsizemax,blsizemax,nblr,nblc))
      write(LG,'(A10,F10.3,A3)')'kin:',   sizeof(kin)   /by2mb,' MB'
      write(LG,'(A10,F10.3,A3)')'ham2d:', sizeof(ham2d) /by2mb,' MB'
      ham2d = 0
    ! Allocate complex arrays
    else
      allocate(kinz(n1,n1),ham2z(blsizemax,blsizemax,nblr,nblc))
      write(LG,'(A10,F10.3,A3)')'kinz:',  sizeof(kinz)  /by2mb,' MB'
      write(LG,'(A10,F10.3,A3)')'ham2z:', sizeof(ham2z) /by2mb,' MB'
      ham2z = 0
    end if

    ! Initialize KEO matrix
    ! In complex version it also includes CAP
    if (realver) then
      call init_matrix1d(kin)
    else
      call init_matrix1z(kinz)
    end if

    ! Calculate blocks of Hamiltonian matrix
    ! Loop over columns
    do iblc=1,nblc
      ! Get global column number
      call l2g(iblc,mycol,nbl,npcol,1,gblc) ! index_conversion
      ! Check if column is not empty and log
      if (adiab .and. blsize(gblc) == 0)cycle
      write(LG,*)'Column: ',iblc,gblc

      ! Loop over rows
      do iblr=1,nblr
        ! Get global row number
        call l2g(iblr,myrow,nbl,nprow,1,gblr)
        ! Check if row is not empty and log
        if (adiab .and. blsize(gblr) == 0)cycle
        write(LG,*)'Row: ',iblr,gblr
        ! Load block
        call load_overlap(gblr,gblc,blsize,ham1d,adiab)
        ! Finish block construction and save, real version
        if (realver) then
          ! Multiply by KEO element
          if (adiab) then
            ham1d = ham1d * kin(gblr,gblc)
          else
            ham1d = ham1d * kin
          end if

          ! Add 2D energies to the diagonal of diagonal block
          if (gblc == gblr) then
            if (adiab) then
              call load_val2(gblc,val2)
            else
              call load_val2_grouped(gblc,val2)
            end if
            do i=1,blsize(gblc)
              ham1d(i,i) = ham1d(i,i) + val2(i)
            end do
          end if

          ! Save block
          do i=1,blsize(gblc)
            ham2d(1:blsize(gblr),i,iblr,iblc) = ham1d(:,i)
          end do
        ! Finish block construction and save, complex version
        else
          ! Allocate complex array
          allocate(ham1z(blsize(gblr),blsize(gblc)))
          ! Multiply by KEO element
          if (adiab) then
            ham1z = ham1d * kinz(gblr,gblc)
          else
            ham1z = ham1d * kinz
          end if

          ! Add 2D energies to the diagonal of diagonal block
          if (gblc == gblr) then
            if (adiab) then
              call load_val2(gblc,val2)
            else
              call load_val2_grouped(gblc,val2)
            end if
            do i=1,blsize(gblc)
              ham1z(i,i) = ham1z(i,i) + val2(i)
            end do
          end if

          ! Save block
          do i=1,blsize(gblc)
            ham2z(1:blsize(gblr),i,iblr,iblc) = ham1z(:,i)
          end do
          ! Deallocate temporary array
          deallocate(ham1z)
        end if
        ! Deallocate temporary array
        deallocate(ham1d)
      ! Loop over rows
      end do
    ! Loop over columns
    end do

    ! Define ScaLapack block size for matrix
    ! For non-hermitian solver, minimum block size is 6
    if (realver) then
      nb = 1
    else
      nb = 6
    end if

    ! Allocate final matrix
    mblr = numroc(msize,nb,myrow,0,nprow)
    mblc = numroc(msize,nb,mycol,0,npcol)
    if (myid == 0)write(*,*)'Final matrix size: ',msize
    write(LG,*)'Final matrix size: ',msize
    write(LG,*)'Local chunk size:',mblr,'x',mblc
    if (realver) then
      if (allocated(hamd))deallocate(hamd)
      allocate(hamd(mblr,mblc))
      write(LG,'(A10,F10.3,A3)')'hamd:',sizeof(hamd)/by2mb,' MB'
    else
      if (allocated(hamz))deallocate(hamz)
      allocate(hamz(mblr,mblc))
      write(LG,'(A10,F10.3,A3)')'hamz:',sizeof(hamz)/by2mb,' MB'
    end if

    ! Allocate space for broadcasted blocks
    nblcmax = nblc
    nblrmax = nblr
    call igamx2d(context,'A',' ',1,1,nblcmax,1,i,i,-1,-1,-1)
    call igamx2d(context,'A',' ',1,1,nblrmax,1,i,i,-1,-1,-1)
    write(LG,*)'Maximum number of blocks: ',nblrmax,'x',nblcmax
    if (realver) then
      allocate(hamtd(blsizemax,blsizemax,nblrmax,nblcmax))
      write(LG,'(A10,F10.3,A3)')'hamtd:',sizeof(hamtd)/by2mb,' MB'
    else
      allocate(hamtz(blsizemax,blsizemax,nblrmax,nblcmax))
      write(LG,'(A10,F10.3,A3)')'hamtz:',sizeof(hamtz)/by2mb,' MB'
    end if
    i = blsizemax**2 * nblrmax

    ! Redistribute elements of matrix in 2d cyclic way
    ! Spread elements of broadcasted blocks
    ! Loop over processes
    do pcol=0,npcol-1
      do prow=0,nprow-1
        ! Broadcast all local blocks of one process
        if (pcol==mycol.and.prow==myrow) then
          if (realver) then
            hamtd = 0d0
          else
            hamtz = 0d0
          end if
          do iblc=1,nblc
            do iblr=1,nblr
              if (realver) then
                hamtd(:,:,iblr,iblc) = ham2d(:,:,iblr,iblc)
              else
                hamtz(:,:,iblr,iblc) = ham2z(:,:,iblr,iblc)
              end if
            end do
          end do
          nblct = nblc
          nblrt = nblr
          if (realver) then
            call dgebs2d(context,'A',' ',i,nblcmax,hamtd,i)
          else
            call zgebs2d(context,'A',' ',i,nblcmax,hamtz,i)
          end if
          call igebs2d(context,'A',' ',1,1,nblct,1)
          call igebs2d(context,'A',' ',1,1,nblrt,1)
        else
          if (realver) then
            call dgebr2d(context,'A',' ',i,nblcmax,hamtd,i,prow,pcol)
          else
            call zgebr2d(context,'A',' ',i,nblcmax,hamtz,i,prow,pcol)
          end if
          call igebr2d(context,'A',' ',1,1,nblct,1,prow,pcol)
          call igebr2d(context,'A',' ',1,1,nblrt,1,prow,pcol)
        end if

        ! Fill in corresponding local elements of final matrix
        do iblc=1,nblct
          call l2g(iblc,pcol,nbl,npcol,1,gblc)
          do iblr=1,nblrt
            call l2g(iblr,prow,nbl,nprow,1,gblr)
            do jc=1,blsize(gblc)
              kc = offset(gblc) + jc
              call g2l(kc,msize,npcol,nb,fpcol,fblc)
              do jr=1,blsize(gblr)
                kr = offset(gblr) + jr
                call g2l(kr,msize,nprow,nb,fprow,fblr)
                if (fpcol==mycol.and.fprow==myrow) then
                  if (realver) then
                    hamd(fblr,fblc) = hamtd(jr,jc,iblr,iblc)
                  else
                    hamz(fblr,fblc) = hamtz(jr,jc,iblr,iblc)
                  end if
                end if
              end do
            end do
          end do
        end do
      ! Loop over processes
      end do
    end do

    ! Allocate memory for values, vectors and symmetries
    nstloc = numroc(nstate,1,myid,0,nprocs_rect)
    if (realver) then
      allocate(val3d(nstate),vec3d(msize,nstloc))
      write(LG,'(A10,F10.3,A3)')'val3d:', sizeof(val3d) /by2mb,' MB'
      write(LG,'(A10,F10.3,A3)')'vec3d:', sizeof(vec3d) /by2mb,' MB'
    else
      allocate(val3z(nstate),vec3z(msize,nstloc))
      write(LG,'(A10,F10.3,A3)')'val3z:', sizeof(val3z) /by2mb,' MB'
      write(LG,'(A10,F10.3,A3)')'vec3z:', sizeof(vec3z) /by2mb,' MB'
    end if

    ! Solve matrix
    if (realver) then
      call scald(context,msize,nstate,vec3d,val3d,hamd,mblr,mblc,nstloc,nb) ! tools/eicalc/scalapack.f90
    else
      call scalz(context,msize,nstate,vec3z,val3z,hamz,mblr,mblc,nstloc,nb)
    end if

    ! Print spectrum
    call prnt_3dsdt(val3d,vec3d,val3z,vec3z,nstloc,nprocs_rect)
  end subroutine

  !-----------------------------------------------------------------------
  !  Solution of 3D problem using SDT and PARPACK for a specific value of K
  !  Variables are real (d) and/or complex (z).
  !-----------------------------------------------------------------------
  subroutine calc_3dsdt_parp
    ! Arrays for 2D states
    real*8,allocatable::val2(:)           ! 2D eivalues
    ! Arrays for real 3D states
    real*8,allocatable::val3d(:)          ! 3D eivalues
    real*8,allocatable::vec3d(:,:)        ! 3D vecs in basis
    ! Arrays for complex 3D states
    complex*16,allocatable::val3z(:)      ! 3D eivalues
    complex*16,allocatable::vec3z(:,:)    ! 3D vecs in basis
    ! Saved matrix
    real*8,allocatable::ham(:,:)          ! One full block of overlaps
    ! Real hamiltonian arrays
    real*8,allocatable::kin(:,:)          ! KEO matrix
    real*8,allocatable::ham1d(:,:)        ! One chunk
    ! Complex hamiltonian arrays
    complex*16,allocatable::kinz(:,:)     ! KEO matrix
    complex*16,allocatable::ham1z(:,:)    ! Part of overlap block relevant for this process
    ! Matrix partitioning
    integer,allocatable::blsize(:)        ! Block sizes (in rows)
    integer,allocatable::offset(:)        ! starting row number of each block (block offsets)
    integer nbl                           ! Number of blocks along rows or columns
    integer blsizemax                     ! Maximum block size (in rows)
    integer gblr,gblc                     ! Global block indices
    ! Variables for Parpack matrix distrubution
    integer mloc    ! Local number of matrix rows
    integer mrem    ! Number of remaining rows if not divisible
    integer rng     ! Global row number
    integer rns     ! Starting row number in starting block
    integer rne     ! Ending   row number in ending   block
    integer bns     ! Block number, starting (for this proc)
    integer bne     ! Block number, ending (for this proc)
    integer rn1     ! Starting row number in current block
    integer rn2     ! Ending   row number in current block
    integer rc      ! Row count in [rn1,rn2]
    integer rct     ! Row count, total
    integer numroc  ! Calculates num of blocks
    integer nstloc  ! Local num of states
    ! Miscellaneous
    integer i,j

    ! test adding asym on overlaps stage
    integer :: sym, msize_1, J_my, K1, K2, par
    integer :: first_row, last_row, first_col, last_col
    integer, allocatable :: blsize_1(:)        ! Block sizes (in rows)
    integer, allocatable :: offset_1(:)        ! starting row number of each block (block offsets)
    real*8, allocatable :: asym(:, :) ! asym block
    real*8, allocatable :: cor(:, :) ! coriolis block
    complex*16, allocatable :: hamz_total(:, :), hamz_1(:, :), hamz_01(:, :)
    character(:), allocatable :: sym_label, sym_label_op

! Original code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Load matrix partitioning
    call load_partitioning(nbl,blsize,blsizemax,offset,msize)
    
    ! Stop if number of states exceeds matrix size
    if (nstate > msize) then
      if (myid == 0) then
        write(*,*)'Requested number of states > matrix size, stop'
      end if
      return
    end if

    ! Find local number of rows and offset
    mloc = msize / nprocs
    mrem = mod(msize, nprocs)
    if (myid < mrem) mloc = mloc + 1 ! distribute remaining rows
    rog = msize / nprocs * myid + min(myid,mrem) ! starting row for current proc (starting from 0)

    ! Check that ncv <= nloc, as required by pzneupd
    if (ncv > mloc) then
      write(*,*)'NCV exceeds MLOC',myid,ncv,mloc
      stop
    end if

    ! Allocate KEO matrix
    if (realver) then
      allocate(kin(n1,n1))
      write(LG,'(A10,F10.3,A3)')'kin:', sizeof(kin) /by2mb,' MB'
    else
      allocate(kinz(n1,n1))
      write(LG,'(A10,F10.3,A3)')'kinz:',sizeof(kinz)/by2mb,' MB'
    end if

    ! Initialize KEO matrix
    ! In complex version it also includes CAP
    if (realver) then
      call init_matrix1d(kin)
    else
      call init_matrix1z(kinz)
    end if

    ! Allocate final matrix
    if (myid == 0) write(*,*)'Final matrix size: ',msize
    write(LG,*)'Final matrix size: ',msize
    write(LG,*)'Local chunk size:',mloc,'x',msize
    if (realver) then
      if (allocated(hamd))deallocate(hamd)
      allocate(hamd(mloc,msize))
      write(LG,'(A10,F10.3,A3)')'hamd:',sizeof(hamd)/by2mb,' MB'
    else
      if (allocated(hamz))deallocate(hamz)
      allocate(hamz(mloc,msize))
      write(LG,'(A10,F10.3,A3)')'hamz:',sizeof(hamz)/by2mb,' MB'
    end if

    ! Find out starting block and starting row in it
    j = rog + 1
    bns = 0
    rns = 0
    do i=1,nbl
      if (offset(i)+1 <= j .and. j <= offset(i)+blsize(i)) then
        bns = i
        rns = j - offset(i)
        exit
      end if
    end do

    ! Stop if no block/row found
    if (bns == 0) then
      write(LG,*)'No starting block/row found'
      stop
    else
      write(LG,*)'Starting block:', bns
      write(LG,*)'Starting row:',   rns
    end if

    ! Find out ending block and ending row in it
    j = rog + mloc
    bne = 0
    rne = 0
    do i=bns,nbl
      if (offset(i)+1 <= j .and. j <= offset(i)+blsize(i)) then
        bne = i
        rne = j - offset(i)
        exit
      end if
    end do

    ! Stop if no block/row found
    if (bne == 0) then
      write(LG,*)'No ending block/row found'
      stop
    else
      write(LG,*)'Ending block:', bne
      write(LG,*)'Ending row:',   rne
    end if

    ! Setup number of processed rows
    rct = 0
    ! Calculate blocks of Hamiltonian matrix
    ! Loop over needed rows
    do gblr=bns,bne
      ! Check if row is not empty
      if (adiab .and. blsize(gblr) == 0) cycle
      
      ! Find out starting row number in current block
      if (gblr == bns) then
        rn1 = rns
      else
        rn1 = 1
      end if

      ! Find out ending row number in current block
      if (gblr == bne) then
        rn2 = rne
      else
        rn2 = blsize(gblr)
      end if

      ! Get row count in current block
      rc = rn2 - rn1 + 1
      ! Loop over all block columns
      do gblc=1,nbl
        ! Check if column is not empty
        if (adiab .and. blsize(gblc) == 0) cycle
        call load_overlap(gblr,gblc,blsize,ham,adiab)

        ! Finish block construction and save, real version
        if (realver) then
          ! Allocate real array
          allocate(ham1d(rc,blsize(gblc)))
          ! Multiply by KEO element
          if (adiab) then
            ham1d = ham(rn1:rn2,:) * kin(gblr,gblc)
          else
            ham1d = ham(rn1:rn2,:) * kin(rn1:rn2,:)
          end if

          ! Add 2D energies to the diagonal of diagonal block
          if (gblc == gblr) then
            if (adiab) then
              call load_val2(gblc,val2)
            else
              call load_val2_grouped(gblc,val2)
            end if
            do i=1,rc
              j = rn1 + i - 1
              ham1d(i,j) = ham1d(i,j) + val2(j)
            end do
          end if

          ! Save block
          do i=1,blsize(gblc)
            hamd( rct+1 : rct+rc, offset(gblc)+i ) = ham1d(:,i)
          end do
          ! Deallocate temporary array
          deallocate(ham1d)
        ! Finish block construction and save, complex version
        else
          ! Allocate complex array
          allocate(ham1z(rc,blsize(gblc)))
          ! Multiply by KEO element
          if (adiab) then
            ham1z = ham(rn1:rn2,:) * kinz(gblr,gblc)
          else
            ham1z = ham(rn1:rn2,:) * kinz(rn1:rn2,:)
          end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! if (test_mode == 'asym') then
          !   if (gblc == gblr) then
          !     asym = load_asym(gblr, blsize)
          !     ham1z = ham1z + asym(rn1:rn2, :) * calculate_U(jlarge, klarge, klarge, parity)
          !   end if
          ! end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! Add 2D energies to the diagonal of diagonal block
          if (gblc == gblr) then
            if (adiab) then
              call load_val2(gblc,val2)
            else
              call load_val2_grouped(gblc,val2)
            end if
            do i=1,rc
              j = rn1 + i - 1
              ham1z(i,j) = ham1z(i,j) + val2(j)
            end do
          end if

          ! Save block
          do i=1,blsize(gblc)
            hamz( rct+1 : rct+rc, offset(gblc)+i ) = ham1z(:,i)
          end do
          ! Deallocate temporary array
          deallocate(ham1z)
        end if
      ! Loop over columns
      end do
      ! Update number of processed rows
      rct = rct + rc
    ! Loop over rows
    end do


! Coriolis hacking
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if (test_mode == 'cor') then
    !   ! Loading second block (K=1) + asym
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   print *, 'Loading K=1 block'
    !   sym_label = iff('even', 'odd', sy == 5)
    !   sym_label_op = iff('odd', 'even', sy == 5)
    !   par = 1
    !   bpath = '/global/cscratch1/sd/gaidai/research/spectrumsdt/calcs/686/rmax_5/alpha_0.5/ecut_0/J_' // num2str(jlarge) // '/K_1/' // sym_label_op // '/basis'
    !   opath = '/global/cscratch1/sd/gaidai/research/spectrumsdt/calcs/686/rmax_5/alpha_0.5/ecut_0/J_' // num2str(jlarge) // '/K_1/' // sym_label_op // '/overlaps'
    !   call load_partitioning(nbl, blsize_1, blsizemax, offset_1, msize_1)
    !
    !   allocate(hamz_1(msize_1, msize_1))
    !   ! Loop over needed row blocks
    !   do gblr = 1, n1
    !     if (blsize_1(gblr) == 0) cycle
    !     ! Loop over column blocks
    !     do gblc = 1, n1
    !       if (blsize_1(gblc) == 0) cycle
    !       call load_overlap(gblr, gblc, blsize_1, ham, adiab)
    !       ! Multiply by KEO element
    !       ham = ham * kinz(gblr,gblc)
    !
    !       ! Add asym term
    !       if (gblc == gblr) then
    !         asym = load_asym(gblr, blsize_1)
    !         asym = asym * calculate_U(jlarge, 1, 1, par)
    !         ham = ham + asym
    !       end if
    !
    !       ! Add 2D energies to the diagonal of diagonal block
    !       if (gblc == gblr) then
    !         call load_val2(gblc,val2)
    !         do i=1, blsize_1(gblr)
    !           ham(i, i) = ham(i, i) + val2(i)
    !         end do
    !       end if
    !
    !       ! Save block
    !       hamz_1(offset_1(gblr) + 1 : offset_1(gblr) + blsize_1(gblr), offset_1(gblc) + 1 : offset_1(gblc) + blsize_1(gblc)) = ham
    !     end do
    !   end do
    !
    !   ! Loading offdiagonal block
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   print *, 'Loading coriolis block'
    !   opath = '/global/cscratch1/sd/gaidai/research/spectrumsdt/calcs/686/rmax_5/alpha_0.5/ecut_0/J_' // num2str(jlarge) // '/K_' // num2str(jlarge - 1) // '/' // sym_label_op // '/overlaps'
    !
    !   allocate(hamz_01(msize, msize_1))
    !   ! Loop over needed row blocks
    !   do gblr = 1, n1
    !     if (blsize(gblr) == 0 .or. blsize_1(gblr) == 0) cycle
    !     ham = load_cor(gblr, blsize(gblr), blsize_1(gblr))
    !
    !     if (jlarge == 1) then
    !       ham = ham * calculate_W(jlarge, 0, 1, par)
    !     end if
    !     if (jlarge == 2) then
    !       ham = ham * calculate_W(jlarge, 1, 2, par)
    !     end if
    !
    !     ! test factor
    !     ! ham = ham / sqrt(2d0)
    !
    !     ! Save block
    !     hamz_01(offset(gblr) + 1 : offset(gblr) + blsize(gblr), offset_1(gblr) + 1 : offset_1(gblr) + blsize_1(gblr)) = ham
    !   end do
    !
    !   ! Putting all the blocks together
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   print *, 'Putting the blocks together'
    !   allocate(hamz_total(msize + msize_1, msize + msize_1))
    !   if (jlarge == 1) then
    !     hamz_total(: msize, : msize) = hamz
    !     hamz_total(msize + 1 :, msize + 1 :) = hamz_1
    !     hamz_total(: msize, msize + 1 :) = hamz_01
    !     hamz_total(msize + 1 :, : msize) = transpose(hamz_01)
    !   else if (jlarge == 2) then
    !     hamz_total(: msize_1, : msize_1) = hamz_1
    !     hamz_total(msize_1 + 1 :, msize_1 + 1 :) = hamz
    !     hamz_total(: msize_1, msize_1 + 1 :) = transpose(hamz_01)
    !     hamz_total(msize_1 + 1 :, : msize_1) = hamz_01
    !   end if
    !
    !   ! Override the variables
    !   deallocate(hamz)
    !   allocate(hamz, source = hamz_total)
    !   msize = msize + msize_1
    !   mloc = msize
    ! end if

! Original code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate memory for values, vectors and symmetries
    nstloc = numroc(nstate,1,myid,0,nprocs)
    if (realver) then
      allocate(val3d(nstate),vec3d(msize,nstloc))
      write(LG,'(A10,F10.3,A3)')'val3d:', sizeof(val3d) /by2mb,' MB'
      write(LG,'(A10,F10.3,A3)')'vec3d:', sizeof(vec3d) /by2mb,' MB'
    else
      allocate(val3z(nstate),vec3z(msize,nstloc))
      write(LG,'(A10,F10.3,A3)')'val3z:', sizeof(val3z) /by2mb,' MB'
      write(LG,'(A10,F10.3,A3)')'vec3z:', sizeof(vec3z) /by2mb,' MB'
    end if

    ! Solve matrix
    if (realver) then
      call pard(context,val3d,vec3d,msize,mloc,nstate,nstloc,ncv,maxitr,opd)
    else
      call parz(context,val3z,vec3z,msize,mloc,nstate,nstloc,ncv,maxitr,opz) ! tools/eicalc/parpack.f90
    end if
    
    ! Print spectrum
    call prnt_3dsdt(val3d,vec3d,val3z,vec3z,nstloc,nprocs)
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates states.
  !-----------------------------------------------------------------------
  subroutine calc_3dsdt_states
    ! State arrays
    real*8,allocatable     :: stated(:)
    complex*16,allocatable :: statez(:)
    ! Arrays for data arrays length
    integer,allocatable::nvec1(:)         ! Number of 1D vectors
    integer,allocatable::nvec2(:)         ! Number of 2D vectors
    integer,allocatable::offset(:)        ! Offsets in final matrix
    ! Miscellaneous
    integer numroc                        ! Calculates num of blocks
    integer nstloc                        ! Local num of states
    character(256) fn1
    character(256) fn2
    character(256) fn3
    integer bstnloc
    integer is,myis
    integer i
    real*8 a,b

    ! Allocate arrays
    if (realver) then
      allocate(stated(nn))
      write(LG,'(A10,F10.3,A3)')'stated:', sizeof(stated)/by2mb,' MB'
    else
      allocate(statez(nn))
      write(LG,'(A10,F10.3,A3)')'statez:', sizeof(statez)/by2mb,' MB'
    end if

    ! Local number of states to build
    bstnloc = numroc(bstn,1,myid,0,nprocs)
    write(LG,*)'Local number of states to build: ',bstnloc
    ! Loop over states
    do myis=1,bstnloc
      ! Get global state number
      call l2g(myis,myid,bstn,nprocs,1,is)
      is = is + bst1 - 1
      ! Build state
      call build_state(is,stated,statez)
      
      ! Build file names
      write(fn1,'(2A,I4.4,A)')outdir,'/state.',is,'.out'
      write(fn2,'(2A,I4.4,A)')outdir,'/state.',is,'.p.out'
      write(fn3,'(2A,I4.4,A)')outdir,'/state.',is,'.n.out'

      ! Open files
      open(1,file=fn1,buffered='yes')
      open(2,file=fn2,buffered='yes')
      open(3,file=fn3,buffered='yes')

      ! Write data
      do i=1,nn
        if (realver) then
          a = abs(stated(i)**2)
          b = stated(i)
        else
          a = abs(statez(i)**2)
          b = dble(statez(i))
        end if
        write(1,'(F25.17)')a
        write(2,'(F25.17)')max(0d0,b)
        write(3,'(F25.17)')min(0d0,b)
      end do

      ! Close files
      close(1)
      close(2)
      close(3)

      ! Log state written
      write(LG,*)'Wrote state ',is
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Post-processing.
  !-----------------------------------------------------------------------
  subroutine calc_3dsdt_post
    ! Barrier
    real*8  baren(ngr)
    real*8  barps(ngr)
    real*8  barpsd                    ! For double well
    real*8  barpsv                    ! For    vdw well
    ! Probabilities
    integer,parameter      :: npb = 4
    real*8,allocatable     :: pb(:,:) ! Probabilities
    real*8,allocatable     :: pd(:,:) ! Probability for each (n, igr) 
    integer i1,i2,i3
    integer i3min,i3max
    real*8  w
    ! Gammas
    integer,parameter      :: ngm = 2
    real*8,allocatable     :: gm(:,:)
    real*8,allocatable     :: cap(:)
    ! State arrays
    real*8,allocatable     :: stated(:)
    complex*16,allocatable :: statez(:) ! 2-levels blocks vector. 1st - n, inside n - l, inside l values at m~. Total size N*L*M~
    ! Miscellaneous
    real*8, parameter      :: vdwmax = 11d0
    integer ist,myist
    integer numroc
    integer nstloc
    integer isl
    integer igr
    integer k
    real*8 en,ga
    character(256) fn

    ! Real version and feshbash are not supported
    if (realver .or. .not. adiab) return

    ! Local number of states to process
    nstloc = numroc(nstate,1,myid,0,nprocs)
    write(LG,*)'Local number of states to process: ',nstloc
    
    ! Load lowest channel
    call load_lowest_barrier(barps,baren)
    ! barps = 5d0
    ! baren = 0d0

    ! Setup arrays
    allocate(pd(n1,ngr), pb(nstate,npb), gm(nstate,ngm), statez(nn), cap(n1))
    pb = 0
    gm = 0
    write(LG,'(A,F10.3,A)')'statez:', sizeof(statez)/by2mb,' MB'

    ! Load CAP
    call init_caps(0d0)
    call get_cap(cap)

    ! Process my states
    do myist=1,nstloc
      ! Get global state number
      call l2g(myist,myid,nstate,nprocs,1,ist)
      write(LG,'(A,I5)')'Processing state ',ist

      ! Build state
      call build_state(ist,stated,statez)

      ! Calculate probability distribution in all phi ranges
      ! Phi is symmetric with respect to pi, so only half of the total range needs to be calculated
      do igr=1,ngr
        i3min = (igr-1) * n3 / 6 + 1
        i3max = (igr)   * n3 / 6
        do i1=1,n1
          w = 0
          do i2=1,n2
            do i3=i3min,i3max
              k = n23*(i1-1) + n3*(i2-1) + i3 ! Index in statez corresponding to current (n, l, m~)
              w = w + abs( statez(k)**2 )
            end do
          end do
          pd(i1,igr) = w * alpha2 * alpha3
        end do
      end do

      ! Double probability distribution due to symmetry
      pd = pd * 2
      ! Calculate gammas, integral of cap * (-2)
      do isl=1,n1
        ! In single channel (B)
        gm(ist,1) = gm(ist,1) + cap(isl) * pd(isl,1) * jac1(isl) * alpha1 * 2
        ! In double channel (A)
        gm(ist,2) = gm(ist,2) + cap(isl) * (pd(isl,2) + pd(isl,3)) * jac1(isl) * alpha1 * 2
      end do

      ! Decide what barrier to use for double well probabilities
      w = abs(gm(ist,1)-gm(ist,2)) / sum(gm(ist,:)) * 2
      if (w < 0.5d0) then
        barpsd = (barps(1) + barps(2)) / 2
      elseif (gm(ist,1) > gm(ist,2)) then
        barpsd = barps(1)
      else
        barpsd = barps(2)
      end if

      ! Decide what barrier to use for vdw well probabilities
      barpsv = (barps(2) + barps(3)) / 2
      ! Calculate probabilities in wells
      do isl=1,n1
        ! Single well
        if (g1(isl) < barps(3)) then
          pb(ist,1) = pb(ist,1) + pd(isl,3) * jac1(isl) * alpha1
        end if

        ! Double well
        if (g1(isl) < barpsd) then
          pb(ist,2) = pb(ist,2) + (pd(isl,1) + pd(isl,2)) * jac1(isl) * alpha1
        end if

        ! Vdw wells associated with channel B
        if (barps(1) < g1(isl) .and. g1(isl) < vdwmax) then
          pb(ist,3) = pb(ist,3) + pd(isl,1) * jac1(isl) * alpha1
        end if

        ! Vdw wells associated with channel A
        if (barpsv < g1(isl) .and. g1(isl) < vdwmax) then
          pb(ist,4) = pb(ist,4) + (pd(isl,2) + pd(isl,3)) * jac1(isl) * alpha1
        end if
      end do
    end do

    ! Gather probabilities
    if (solver /= SOLVER_LAPA) then
      call dgsum2d(context,'A',' ',nstate,npb,pb,nstate,0,0)
      call dgsum2d(context,'A',' ',nstate,ngm,gm,nstate,0,0)
    end if

    ! Exit if not a root
    if (myid /= 0)return

    ! Write detailed spectrum file
    write(fn, '(4A)') dpath, '/', getdir(MODE_3DSDT), '/spec.out'
    open(1,file=fn)
    write(fn,'(2A)')outdir,'/spec.out'
    open(2,file=fn)
    do ist=1,nstate
      read (1,'(4X, 2F30.17)')en,ga
      write(2,'(I4,11F30.17)')ist, en, en - baren(1), en - baren(2), en - baren(3), ga, gm(ist,:) * autown, pb(ist,:)
    end do
    close(1)
    close(2)
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates probability distribution (pd) as a function of 1st coord.
  !-----------------------------------------------------------------------
  subroutine calc_pd(pd,stated,statez,ir)
    real*8,allocatable     :: pd(:)
    real*8,allocatable     :: stated(:)
    complex*16,allocatable :: statez(:)
    integer ir,i1,i2,i3,l
    integer i3min,i3max
    real*8 w
    i3min = (ir-1) * n3 / 6 + 1
    i3max = (ir)   * n3 / 6
    do i1=1,n1
      w = 0
      do i2=1,n2
        do i3=i3min,i3max
          l = n2*n3*(i1-1) + n3*(i2-1) + i3
          if (realver) then
            w = w + abs( stated(l)**2 )
          else
            w = w + abs( statez(l)**2 )
          end if
        end do
      end do
      pd(i1) = w * alpha2 * alpha3
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Recognizes channels.
  !-----------------------------------------------------------------------
  subroutine calc_chrecog
    integer,parameter :: nu = 9         ! Number of output files
    integer,allocatable :: nvec2(:)     ! Number of 2D vectors
    real*8, allocatable :: val2all(:,:) ! 2D eigenvalues
    real*8, allocatable :: vec2(:)      ! 2D eigenvector
    real*8, allocatable :: sim(:,:)     ! Similarities
    real*8, allocatable :: prof(:)      ! Recognized energy profile
    integer,allocatable :: ind(:)       ! Index of recog'd state
    integer,allocatable :: phs(:)       ! Phase of recog'd state
    integer nrec                        ! Number of states to recog
    integer irec                        ! State number to recognize
    integer ist                         ! State number
    integer isl                         ! Slice number
    integer islmin                      ! Slice number of global min
    integer istmin                      ! State number of global min
    integer iu                          ! Output file number
    integer u                           ! Output file unit
    logical recognized                  ! State was recognized
    integer ncase,i,j
    character(256) fn
    ! Barrier properties
    integer gr                          ! Group number
    real*8  bar(3)                      ! Barrier props: pos,nrg,coe
    ! Variables for rexnext
    integer islp,istp,phsp              ! Prev slice
    integer isln,istn,phsn              ! Next slice
    ! Cumulative arrays on root
    integer,allocatable :: indall(:,:)  ! Indices of all states
    integer,allocatable :: phsall(:,:)  ! Phases  of all states
    integer,allocatable :: grall(:)     ! Group numbers of all states
    real*8, allocatable :: barall(:,:)  ! Barrier props of all states

    ! Load number of 2D eigenvalues
    call load_nvec2(nvec2)
    
    ! Load 2D eigenvalues
    call load_val2all(val2all)

    ! Get number of states to recognize
    nrec = nprocs
    if (nrec > nvec2(n1)) then
      write(LG,*)'No states to recognize, exiting'
      return
    end if

    ! Get my state number to recognize
    irec = myid + 1
    write(LG,*)'State number to recognize ',irec

    ! Allocate arrays
    allocate(sim(nvec2sch,n1),ind(n1),phs(n1),prof(n1))
    ! Start from asymptote
    if (recstart == RECSTART_ASYM) then
      ! Set up previous slice
      ind(n1) = irec
      phs(n1) = 1
      istp = irec
      phsp = 1

      ! Loop over slices from right to left
      do islp=n1,2,-1
        isln = islp - 1
        call recnext(nvec2,sim,islp,istp,phsp,isln,istn,phsn)
        ind(isln) = istn
        phs(isln) = phsn
        istp = istn
        phsp = phsn
      end do
    ! Start from the bottom of the well
    else
      ! Find global minimum approximately
      islmin = 1
      do while(g1(islmin) < 4d0)
        islmin = islmin + 1
      end do

      ! Set up previous slice
      ind(islmin) = irec
      phs(islmin) = 1
      istp = irec
      phsp = 1

      ! Loop over slices to the left
      do islp=islmin,2,-1
        isln = islp - 1
        call recnext(nvec2,sim,islp,istp,phsp,isln,istn,phsn)
        ind(isln) = istn
        phs(isln) = phsn
        istp = istn
        phsp = phsn
      end do

      ! Set up previous slice
      istp = irec
      phsp = 1

      ! Loop over slices to the right
      do islp=islmin,n1-1
        isln = islp + 1
        call recnext(nvec2,sim,islp,istp,phsp,isln,istn,phsn)
        ind(isln) = istn
        phs(isln) = phsn
        istp = istn
        phsp = phsn
      end do
    end if

    ! Write recognized states
    if ( basout == BASOUT_2D .or. basout == BASOUT_BOTH ) then
      do isl=1,n1
        call load_vec2(vec2,isl,ind(isl),.true.)
        call writ_vec2(vec2,irec,isl,phs(isl))
      end do
    end if

    ! Write similarities
    write(fn,'(2A,I4.4,A)')outdir,'/sim.',irec,'.out'
    open(1,file=fn)
    do ist=nvec2sch,1,-1
      do isl=1,n1
        write(1,'(F25.17)',advance='no')sim(ist,isl)
      end do
      write(1,*)
    end do
    close(1)

    ! Build energy profile
    do isl=1,n1
      prof(isl) = val2all(ind(isl),isl)
    end do

    ! Find global minimum
    islmin = 1
    do isl=1,n1
      if (prof(isl) < prof(islmin)) then
        islmin = isl
      end if
    end do

    ! Log found minimum
    write(LG,*)'Minimum: ',islmin,prof(islmin)*autown
    ! Go right from minimum until barrier
    do isl=islmin,n1-1
      if (prof(isl) > prof(isl+1))exit
    end do
    if (isl == n1) then
      write(*,*)'No barrier found: ',irec
      do isl=1,n1
        if (g1(isl) > 5.65d0)exit
      end do
      write(*,*)'Approximate barrier slice: ',isl
    end if

    ! Find parabola
    call find_parabola(g1(isl-1),prof(isl-1), g1(isl  ),prof(isl  ), g1(isl+1),prof(isl+1), bar(1),bar(2),bar(3)) ! tools/parabola.f90
    ! Log found barrier
    write(LG,*)'Barrier: ',isl,bar(1),bar(2)*autown,bar(3)*autown

    ! Find a group number the state belongs to
    call load_vec2(vec2,isl,ind(isl),.true.)
    gr = find_group(vec2)

    ! Allocate collective arrays on root
    if (myid == 0) then
      allocate(indall(n1,nrec), phsall(n1,nrec), grall(nrec), barall(3,nrec))
      indall = 0
      phsall = 0
      grall  = 0
      barall = 0
    end if

    ! Colect indices
    if (myid == 0) then
      indall(:,1) = ind
      phsall(:,1) = phs
      barall(:,1) = bar
      grall   (1) = gr
      do irec=2,nrec
        call igerv2d(context,n1,1,indall(1,irec),n1,0,irec-1)
        call igerv2d(context,n1,1,phsall(1,irec),n1,0,irec-1)
        call dgerv2d(context, 3,1,barall(1,irec), 3,0,irec-1)
        call igerv2d(context, 1,1,grall   (irec), 1,0,irec-1)
      end do
    else
      call igesd2d(context,n1,1,ind,n1,0,0)
      call igesd2d(context,n1,1,phs,n1,0,0)
      call dgesd2d(context, 3,1,bar, 3,0,0)
      call igesd2d(context, 1,1,gr,  1,0,0)
    end if

    ! Only root continues
    if (myid /= 0)return

    ! Open files for results of recognition
    open(11,file=outdir//'/val2.rec.out')
    open(12,file=outdir//'/val2.rec.1.out')
    open(13,file=outdir//'/val2.rec.2.out')
    open(14,file=outdir//'/val2.rec.3.out')
    open(15,file=outdir//'/recind.dat')
    open(16,file=outdir//'/recind.1.out')
    open(17,file=outdir//'/recind.2.out')
    open(18,file=outdir//'/recind.3.out')
    open(19,file=outdir//'/recphs.dat')

    ! Write header
    do iu=1,nu
      u = 10 + iu
      write(u,'(10X)',advance='no')
      do isl=1,n1
        write(u,'(I25)',advance='no')isl
      end do
      write(u,*)
    end do

    ! Write content
    do irec=1,nrec
      gr = grall(irec)
      ! Write state number
      write(11,   '(I10)',advance='no')irec
      write(11+gr,'(I10)',advance='no')irec
      write(15,   '(I10)',advance='no')irec
      write(15+gr,'(I10)',advance='no')irec
      write(19,   '(I10)',advance='no')irec

      ! Write results
      do isl=1,n1
        ist = indall(isl,irec)
        write(11,   '(F25.17)',advance='no')val2all(ist,isl) * autown
        write(11+gr,'(F25.17)',advance='no')val2all(ist,isl) * autown
        write(15,   '(I25)',   advance='no')indall (isl,irec)
        write(15+gr,'(I25)',   advance='no')indall (isl,irec)
        write(19,   '(I25)',   advance='no')phsall (isl,irec)
      end do

      ! Write new line
      write(11,   *)
      write(11+gr,*)
      write(15,   *)
      write(15+gr,*)
      write(19,   *)
    end do

    ! Close files
    do iu=1,nu
      u = 10 + iu
      close(u)
    end do

    ! Write channels properties
    ! 1. Channel number
    ! 2. Group number
    ! 3. Barrier position
    ! 4. Barrier energy
    ! 5. Threshold energy
    ! 6. Barrier energy wrt threshold
    ! 7. Parabola coefficient for well
    open(1,file=outdir//'/channels.dat')
    do irec=1,nrec
      write(1,'(2I35,5F35.17)') irec, grall(irec), barall(1,irec), barall(2,irec) * autown, val2all(irec,n1) * autown, barall(2,irec) * autown - val2all(irec,n1) * autown, barall(3,irec) * autown
    end do
    close(1)

    ! Check results for ambiguity
    open(1,file=outdir//'/ambiguities.out')
    ncase = 0
    do isl=1,n1-1
      do i=1,nrec-1
        do j=i+1,nrec
          if (indall(isl,i) == indall(isl,j)) then
            ncase = ncase + 1
            write(1,*)isl,i,j
          end if
        end do
      end do
    end do
    close(1)
    write(*,*)'Number of ambiguous cases: ',ncase

    ! Print initial energies, excluding those recognized
    open(1,file=outdir//'/val2.rest.out')
    write(1,'(10X)',advance='no')
    do j=1,n1
      write(1,'(I25)',advance='no')j
    end do
    write(1,*)
    do i=1,nvec2max
      write(1,'(I10)',advance='no')i
      do j=1,n1
        ! Check if the state number was taken
        recognized = .false.
        do ist=1,nrec
          if (indall(j,ist) == i) then
            recognized = .true.
          end if
        end do

        ! Write space or energy
        if (recognized) then
          write(1,'(25X)',advance='no')
        else
          write(1,'(F25.17)',advance='no')val2all(i,j) * autown
        end if
      end do
      write(1,*)
    end do
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  !  Recognizes the state in the next slice.
  !-----------------------------------------------------------------------
  subroutine recnext(nvec2,sim,islp,istp,phsp,isln,istn,phsn)
    integer,allocatable :: nvec2(:)     ! Number of 2D vectors
    real*8, allocatable :: sim(:,:)     ! Similarities
    real*8, allocatable :: omtx(:,:)    ! Overlap matrix
    integer islp,istp,phsp              ! Prev slice
    integer isln,istn,phsn              ! Next slice
    integer ist                         ! State number
    real*8  olap                        ! Overlap
    real*8  maxolap                     ! Maximum overlap

    ! Load overlap
    call load_overlap(islp,isln,nvec2,omtx,.true.)
    ! Set up maximum overlap
    maxolap = 0

    ! Loop over states in slice
    do ist=1,min(nvec2(islp),nvec2sch)
      ! Get overlap
      olap = omtx(istp,ist)
      ! Save similarity
      if (islp < isln) then
        sim(ist,islp) = olap
      else
        sim(ist,isln) = olap
      end if

      ! Update next state
      if (abs(olap) > abs(maxolap)) then
        maxolap = olap
        istn = ist
        if (olap > 0) then
          phsn = + phsp
        else
          phsn = - phsp
        end if
      end if
    end do

    ! Log processed slice
    write(LG,*)'Slice processed: ',isln,istn,phsn,maxolap
  end subroutine

  !-----------------------------------------------------------------------
  !  Performs diagonalization of 1D channel Hamiltonians.
  !-----------------------------------------------------------------------
  subroutine calc_chdiag
    ! Digonalization
    real*8, allocatable :: olap(:,:)    ! Overlap matrix
    integer,allocatable :: ind(:)       ! Channel indices
    integer,allocatable :: phs(:)       ! Channel phases
    real*8, allocatable :: val2all(:,:) ! 2D eigenvalues
    integer ic,ir                       ! Running block numbers
    ! Real branch arrays
    real*8, allocatable :: hamd(:,:)    ! Hamiltonian matrix
    real*8, allocatable :: chvald(:)    ! Channel eigenvalues
    real*8, allocatable :: chvalalld(:,:)! Cumulative eivals
    ! Complex branch arrays
    complex*16, allocatable :: hamz(:,:)! Hamiltonian matrix
    complex*16, allocatable :: chvalz(:)! Channel eigenvalues
    complex*16, allocatable :: w(:)     ! Channel eivals, unsorted
    complex*16, allocatable :: vr(:,:)  ! Channel eivecs, unsorted
    real*8,     allocatable :: prob(:,:)! Probabilities
    ! Sorting
    real*8, allocatable :: sortval(:)   ! Sorted values
    integer,allocatable :: sortind(:)   ! Sorted indices
    ! Channel properties
    integer nchan                       ! Number of channels
    integer ichan                       ! Channel number
    real*8  barps                       ! Barrier position
    real*8  baren                       ! Barirer energy
    real*8  capebar                     ! Barrier energy for CAP
    ! Miscellaneous
    character(:),allocatable::ldir      ! Load directory
    integer ist                         ! 2D state number
    integer isl                         ! Slice number
    real*8  en                          ! Energy
    character(256) fn
    integer i,j

    ! Get number of channels
    nchan = nprocs
    ! Get my channel number and load recognition
    ichan = myid + 1
    if (load_recognition_chan(ind,phs,ichan)) then
      write(LG,*)'Channel number ',ichan
    else
      write(LG,*)'Nothing to do, exiting'
      return
    end if

    ! Load 2D eigenvalues
    call load_val2all(val2all)
    ! Complex version preparation
    if (.not.realver) then
      ! Get channel properties
      ldir = getrecldpath()
      open(1,file=ldir//'/channels.dat')
      do i=1,ichan
        read(1,'(70X,2F35.17,35X,F35.17)')barps,baren,capebar
      end do
      close(1)
      write(LG,*)'Barrier Position: ',barps
      write(LG,*)'CAP Ebar: ',capebar
      ! Initialize CAPs
      call init_caps(capebar / autown)
    end if

    ! Load overlap matrix
    allocate(olap(n1,n1))
    ldir = opath // '/' // getdir(MODE_OVERLAP)
    write(fn,'(2A,I4.4,A)')ldir,'/overlap.',ichan,'.bin.out'
    open(1,file=fn,form='unformatted')
    read(1)olap
    close(1)

    ! Real branch
    if (realver) then
      ! Initialize real KEO matrix
      allocate(hamd(n1,n1))
      write(LG,'(A10,F10.3,A3)')'hamd:',  sizeof(hamd)  /by2mb,' MB'
      call init_matrix1d(hamd)

      ! Multiply KEO by overlap matrix element-wise
      hamd = olap * hamd
      ! Add energy to diagonal
      do isl=1,n1
        hamd(isl,isl) = hamd(isl,isl) + val2all(ind(isl),isl)
      end do

      ! Solve matrix
      allocate(chvald(n1))
      call syev(hamd,chvald,'V')

      ! Log eigenvalues
      write(LG,*)'Spectrum:'
      do isl=1,n1
        write(LG,'(I5,F25.17)')isl,chvald(isl) * autown
      end do

      ! Normalize states
      do ist=1,n1
        hamd(:,ist) = hamd(:,ist) / sqrt(jac1 * alpha1)
      end do

      ! Write states in binary form
      write(fn,'(2A,I4.4,A)')outdir,'/exp.',ichan,'.bin.out'
      open(1,file=fn,form='unformatted')
      write(1)hamd
      close(1)

      ! Write states in text form
      write(fn,'(2A,I4.4,A)')outdir,'/exp.',ichan,'.out'
      open(1,file=fn)
      do i=1,n1
        do j=1,n1
          write(1,'(F25.17)',advance='no')hamd(i,j)
        end do
        write(1,*)
      end do
      close(1)

      ! Allocate collective arrays on root
      if (myid == 0) then
        allocate(chvalalld(n1,nchan))
        chvalalld = 0
      end if

      ! Colect eigenvalues
      if (myid == 0) then
        chvalalld(:,1) = chvald
        do ichan=2,nchan
          call dgerv2d(context,n1,1,chvalalld(1,ichan),n1,0,ichan-1)
        end do
      else
        call dgesd2d(context,n1,1,chvald,n1,0,0)
      end if

      ! Only root continues
      if (myid /= 0)return

      ! Write channel eigenvalues
      open(1,file='en.out')
      write(1,'(10X)',advance='no')
      do ichan=1,nchan
        write(1,'(I25)',advance='no')ichan
      end do
      write(1,*)
      do ist=1,n1
        write(1,'(I10)',advance='no')ist
        do ichan=1,nchan
          write(1,'(F25.17)',advance='no')chvalalld(ist,ichan)*autown
        end do
        write(1,*)
      end do
      close(1)
    ! Complex branch
    else
      ! Initialize complex KEO matrix
      allocate(hamz(n1,n1))
      write(LG,'(A10,F10.3,A3)')'hamz:',  sizeof(hamz)  /by2mb,' MB'
      call init_matrix1z(hamz)

      ! Multiply KEO by overlap matrix element-wise
      hamz = olap * hamz
      ! Add energy to diagonal
      do isl=1,n1
        hamz(isl,isl) = hamz(isl,isl) + val2all(ind(isl),isl)
      end do

      ! Solve matrix
      allocate(chvalz(n1),w(n1),vr(n1,n1),prob(n1,5))
      call geev(hamz,w,vr,vr)

      ! Sort eigenpairs
      allocate(sortval(n1),sortind(n1))
      sortval = real(w)
      sortval = bubble_sort(sortval, sortind)
      do i=1,n1
        chvalz(i) = w (  sortind(i))
        hamz(:,i) = vr(:,sortind(i))
      end do

      ! Log eigenvalues
      write(LG,*)'Spectrum:'
      do isl=1,n1
        write(LG,'(I5,2F25.17)')isl,chvalz(isl) * autown
      end do

      ! Calculate probabilities
      prob = 0
      do ist=1,n1
        do isl=1,n1
          ! Real part only
          prob(ist,2) = prob(ist,2) + real(hamz(isl,ist))**2
          ! Imaginary part only
          prob(ist,3) = prob(ist,3) + aimag(hamz(isl,ist))**2

          if (g1(isl) < barps) then
            ! Probability in the well
            prob(ist,1) = prob(ist,1) + conjg(hamz(isl,ist)) * hamz(isl,ist)
            ! Real part in the well
            prob(ist,5) = prob(ist,5) + real(hamz(isl,ist))**2
          end if
        end do

        ! Difference between real and imaginary parts
        prob(ist,4) = prob(ist,2) - prob(ist,3)
        ! Divide real part in the well by real part
        prob(ist,5) = prob(ist,5) / prob(ist,2)
      end do

      ! Save spectrum
      write(fn,'(2A,I4.4,A)')outdir,'/spec.',ichan,'.out'
      open(1,file=fn)
        do ist=1,n1
          write(1,'(8F30.17)') real(chvalz(ist))  * autown, real(chvalz(ist))  * autown - baren, aimag(chvalz(ist)) * autown * (-2), (prob(ist,i),i=1,5)
        end do
      close(1)

      ! Normalize states
      do ist=1,n1
        hamz(:,ist) = hamz(:,ist) / sqrt(jac1 * alpha1)
      end do

      ! Write states in binary form
      write(fn,'(4A,I4.4,A)') outdir,'/',expdir,'/exp.',ichan,'.bin.out'
      open(1,file=fn,form='unformatted')
      write(1)hamz
      close(1)

      ! Write states in text form
      write(fn,'(4A,I4.4,A)') outdir,'/',expdir,'/exp.',ichan,'.re.out'
      open(1,file=fn)
      write(fn,'(4A,I4.4,A)') outdir,'/',expdir,'/exp.',ichan,'.im.out'
      open(2,file=fn)
      do i=1,n1
        do j=1,n1
          write(1,'(F25.17)',advance='no')real(hamz(i,j))
          write(2,'(F25.17)',advance='no')aimag(hamz(i,j))
        end do
        write(1,*)
        write(2,*)
      end do
      close(1)
      close(2)
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Returns a group number 2D state belongs to.
  !-----------------------------------------------------------------------
  integer function find_group(vec2)
    real*8, allocatable::vec2(:)
    integer isl     ! Slice number
    integer ist     ! State number
    integer ir      ! Range number for phi coordinate
    integer i2      ! Second coordinate index
    integer i3      ! Third  coordinate index
    integer i3min   ! Third  coordinate index, min
    integer i3max   ! Third  coordinate index, max
    integer k       ! Point index
    real*8  frac    ! Fraction
    integer irmax   ! Range number of maximum
    real*8  fracmax ! Maximum fraction

    ! Find where state is localized
    ! Calculates probability in a half of 2D space due to symmetry
    fracmax = 0
    do ir=1,3
      ! Figure out range
      i3min = (ir-1) * n3 / 6 + 1
      i3max = (ir)   * n3 / 6

      ! Compute fraction in the range
      frac = 0
      do i2=1,n2
        do i3=i3min,i3max
          k = (i2-1) * n3 + i3
          frac = frac + vec2(k)**2
        end do
      end do
      frac = frac * alpha2 * alpha3

      ! Log fraction
      write(LG,*)'Fraction: ',ir,frac
      ! Update range number of maximum fraction
      if (frac > fracmax) then
        irmax = ir
        fracmax = frac
      end if
    end do

    ! Log group number
    write(LG,*)'Group number: ',irmax
    ! Return group number
    find_group = irmax
  end function

  !-----------------------------------------------------------------------
  !  Calculates states after channel diagonalization.
  !-----------------------------------------------------------------------
  subroutine calc_chdiag_states
    ! Real branch arrays
    real*8,allocatable :: vec3d(:)     ! Vector in grid
    real*8,allocatable :: vec3ed(:,:)  ! Vector in eibasis
    real*8,allocatable :: sliced(:)    ! Slice in basis
    ! Miscellaneous
    integer,allocatable :: ind(:)      ! Channel indices
    integer,allocatable :: phs(:)      ! Channel phases
    integer nchan                      ! Number of channels
    integer ichan                      ! Channel number
    integer is                         ! 3D state number
    integer ist                        ! 2D state number
    integer isl                        ! Slice number
    real*8  vol12                      ! Volume element for 1 and 2
    real*8  a,b                        ! Wave function value
    real*8  sym                        ! Symmetry
    character(256)fn                   ! File name
    integer i,l

    ! Get number of channels
    nchan = nprocs
    ! Get my channel number and load recognition
    ichan = myid + 1
    if (load_recognition_chan(ind,phs,ichan)) then
      write(LG,*)'Channel number ',ichan
    else
      write(LG,*)'Nothing to do, exiting'
      return
    end if

    ! Allocatge arrays
    allocate(vec3d(nn),vec3ed(n1,n1),sliced(n23b))
    ! Load expansion
    write(fn,'(4A,I4.4,A)') getdir(MODE_CHDIAG),'/',expdir,'/exp.',ichan,'.bin.out'
    open(1,file=fn,form='unformatted')
    read(1)vec3ed
    close(1)

    ! Loop over states
    do is=bst1,bst1+bstn-1
      ! Loop over slices
      do isl=1,n1
        ! Load 2d state
        call load_vec2(sliced,isl,ind(isl),.true.)

        ! Multiply by expansion coefficient and phase
        l = (isl-1) * n23
        vec3d(l+1:l+n23) = sliced * vec3ed(isl,is) * phs(ind(isl))
      end do

      ! Open files
      write(fn,'(2A,I4.4,A,I4.4,A)') outdir,'/state',ichan,'.',is,'.out'
      open(1,file=fn,buffered='yes')
      write(fn,'(2A,I4.4,A,I4.4,A)') outdir,'/state',ichan,'.',is,'.p.out'
      open(2,file=fn,buffered='yes')
      write(fn,'(2A,I4.4,A,I4.4,A)') outdir,'/state',ichan,'.',is,'.n.out'
      open(3,file=fn,buffered='yes')

      ! Write data
      do i=1,nn
        a = abs(vec3d(i)**2)
        b = vec3d(i)
        write(1,'(F25.17)')a
        write(2,'(F25.17)')max(0d0,b)
        write(3,'(F25.17)')min(0d0,b)
      end do

      ! Close files
      close(1)
      close(2)
      close(3)

      ! Check symmetry
      call symmetryd(vec3d,sym)
      ! Log state written and symmetry
      write(LG,'(A,I4.4,F25.17)')'Wrote state ',is,sym
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates FBR basis on the grid.
  !     Symmetry       Cos(nx)         Sin(nx)
  !    A1 and A2       0,3,6,9           3,6,9,12
  !    E1 and E2       1,2,4,5           1,2,4,5
  !    SY and AS       0,1,2,3           1,2,3,4
  !-----------------------------------------------------------------------
  subroutine init_fbrbasis(basis)
    implicit none
    real*8 basis(n3,n3b)
    real*8 norm
    integer i,j,k,np

    ! Normalization coefficient
    norm = sqrt(1/pi)
    if (onewell)norm = norm * sqrt(3d0)

    ! Symmetries E1 and E2: exclude 3-fold basis functions
    if (sy==SY_E1.or.sy==SY_E2) then
      k = 0
      do j=1,n3b
        k = k + 1
        if (mod(k,3)==0)k = k + 1
        do i=1,n3
          if (sy==SY_E1) then
            basis(i,j) = cos(k * g3(i)) * norm
          else
            basis(i,j) = sin(k * g3(i)) * norm
          endif
        enddo
      enddo
    ! Other symmetries: use 1-fold or 3-fold basis functions
    else
      ! Periods coefficient
      np = 3
      if (sy==SY_SY.or.sy==SY_AS)np = 1

      ! Fill in basis array
      do j=1,n3b
        do i=1,n3
          if (sy==SY_A1.or.sy==SY_SY) then
            basis(i,j) = cos(np * (j-1) * g3(i)) * norm
            if (j==1) basis(i,j) = norm / sqrt(2d0)
          else
            basis(i,j) = sin(np * j * g3(i)) * norm
          endif
        enddo
      enddo
    endif
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads 1D and 2D eigenvector basis for given slice.
  !-----------------------------------------------------------------------
  subroutine load_eibasis(isl,nvec1,val1,vec1,val2,vec2)
    type(array1d),allocatable::val1(:) ! 1D solutions eigenvalues for each thread
    type(array2d),allocatable::vec1(:) ! 1D solutions expansion coefficients (over sin/cos) for each thread
    real*8,allocatable::vec2(:,:)      ! 2D solutions expansion coefficients (over 1D solutions)
    real*8,allocatable::val2(:)        ! 2D eigenvalues
    integer,allocatable::nvec1(:)      ! Num of 1D solutions in each thread
    integer nvec2                      ! Num of 2D vecs
    integer nb2                        ! Size of 2D vector in basis
    integer isl                        ! Slice number
    integer i
    character(:),allocatable::ldir
    character(256) fn
    
    ! Deallocate arrays if allocated
    call free_1d(nvec1,val1,vec1)
    call free_2d(nvec2,val2,vec2)

    ! Get load directory
    ldir = bpath // '/' // getdir(MODE_BASIS)
    ! Get file name
    write(fn,'(2A,I4.4,A)')ldir,'/bas2.',isl,'.bin.out'
    ! Load 2D solution
    open(1,file=fn,form='unformatted')
    read(1)nvec2,nb2
    if (nvec2 == 0)return ! Exit right away, if empty
    allocate(val2(nvec2),vec2(nb2,nvec2))
    read(1)val2
    read(1)vec2
    close(1)

    ! Get file name
    write(fn,'(2A,I4.4,A)')ldir,'/bas1.',isl,'.bin.out'
    ! Load 1D solution
    allocate(nvec1(n2),val1(n2),vec1(n2))
    open(1,file=fn,form='unformatted')
    read(1)nvec1
    do i=1,n2
      if (nvec1(i) == 0)cycle
      allocate(val1(i)%a(nvec1(i)),vec1(i)%a(n3b,nvec1(i)))
      read(1)val1(i)%a
      read(1)vec1(i)%a
    end do
    close(1)
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
    real*8,allocatable::vec(:,:)      ! 2D vectors in basis
    real*8,allocatable::val(:)        ! 2D values
    integer nvec                      ! Num of 2D vecs
    nvec = 0
    if (allocated(val)) deallocate(val)
    if (allocated(vec)) deallocate(vec)
  end subroutine

  !-----------------------------------------------------------------------
  !  Transforms 2D state from 1D eigenvector basis to DVR or FBR basis.
  !  On exit vec2b contains a vector of expansion coefficients of a given 2D solution over sin/cos
  !  Vector has M blocks of length L each. Each M-block contains expansion coefficient over the same harmonic in different ls.
  !  Each element of the vector is sum over i of a_nlm^i * b_nli^j (j is index of vec2e and is fixed within this subroutine)
  !-----------------------------------------------------------------------
  subroutine trans_eibas_bas(vec1,nvec1,vec2e,vec2b)
    type(array2d) vec1(:)    ! i-th element is a 2D array of 1D solutions in the i-th thread given as expansion coefficients over sin/cos basis
    real*8        vec2b(:)   ! 2D vector in basis
    real*8        vec2e(:)   ! expansion coefficients over 1D solutions
    integer nvec1(n2)        ! Number of 1D vectors in each thread
    integer is,i,j
    i = 0
    do is=1,n2
      j = (is-1) * n3b
      call gemv(vec1(is)%a,vec2e(i+1:i+nvec1(is)),vec2b(j+1:j+n3b))
      i = i + nvec1(is)
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Transform 2D state from basis to grid. Real state.
  !-----------------------------------------------------------------------
  subroutine trans_bas_grd_d(basis,vec2bas,vec2grd)
    real*8 basis(:,:)    ! FBR basis
    real*8 vec2bas(:)    ! 2D vector in basis
    real*8 vec2grd(:)    ! 2D vector on grid
    integer is,i,j
    if (dvr) then
      vec2grd = vec2bas / sqrt(alpha3)
    else
      do is=1,n2
        i = (is-1) * n3b
        j = (is-1) * n3
        call gemv(basis,vec2bas(i+1:i+n3b),vec2grd(j+1:j+n3))
      end do
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Transform 2D state from basis to grid. Complex state.
  !-----------------------------------------------------------------------
  subroutine trans_bas_grd_z(basis,vec2bas,vec2grd)
    real*8     basis(:,:)    ! FBR basis
    complex*16 vec2bas(:)    ! 2D vector in basis
    complex*16 vec2grd(:)    ! 2D vector on grid
    integer is,i,j
    if (dvr) then
      vec2grd = vec2bas / sqrt(alpha3)
    else
      do is=1,n2
        i = (is-1) * n3b
        j = (is-1) * n3
        call gemv(basis,vec2bas(i+1:i+n3b),vec2grd(j+1:j+n3))
      end do
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Builds 3D state on grid from basis representation.
  !-----------------------------------------------------------------------
  subroutine build_state(ist,vec3d,vec3z)
    ! 3D vectors, real
    real*8,allocatable::vec3d(:)       ! Vector in grid
    real*8,allocatable::vec3ed(:)      ! Vector in eibasis
    real*8,allocatable::vd(:)          ! 1 n-block of 3D expansion
    real*8,allocatable::sliced(:)      ! Slice in basis
    ! 3D vectors, complex
    complex*16,allocatable::vec3z(:)   ! Vector in grid (output)
    complex*16,allocatable::vec3ez(:)  ! Vector in eibasis
    complex*16,allocatable::vz(:)      ! 1 n-block of 3D expansion
    complex*16,allocatable::slicez(:)  ! Slice in basis
    ! 1D and 2D basis
    type(array1d),allocatable::val1(:) ! 1D values  for each thread
    type(array2d),allocatable::vec1(:) ! 1D vectors for each thread
    real*8,allocatable::vec2e(:,:)     ! 2D vectors in eibasis repr
    real*8,allocatable::vec2b(:,:)     ! 2D vectors in   basis repr
    real*8,allocatable::val2(:)        ! 2D values
    real*8,allocatable::basis(:,:)     ! FBR basis
    integer,allocatable::nvec1(:)      ! Num of 1D vectors for each l within given n
    ! Matrix partitioning
    integer,allocatable::blsize(:)     ! n-block sizes
    integer,allocatable::offset(:)     ! n-block offsets
    integer msize                      ! Matrix size
    integer nbl                        ! Number of n-blocks along rows or columns
    integer blsizemax                  ! Maximum n-block size
    ! Channel data and variables
    integer,allocatable::chi(:)        ! Channel numbers for pthway S
    integer,allocatable::ind(:,:)      ! Channel indices
    integer,allocatable::phs(:,:)      ! Channel phases
    integer ichan                      ! Channel number
    integer nchan                      ! Number of channels
    integer nchmax                     ! Minimum nch over pathways
    ! Miscellaneous
    integer ist                        ! 3D state number
    integer isl                        ! Slice number
    integer bassize                    ! Basis size for one n-slice
    real*8,allocatable::vec2t(:)       ! Temporary 2D state
    real*8  vol12                      ! Volume element for 1 and 2
    integer k,l

    ! Load matrix partitioning
    call load_partitioning(nbl,blsize,blsizemax,offset,msize)
    
    ! Load recognition and channel grouping if needed
    if (.not.adiab) then
      call load_recognition(ind,phs,nchan)
      call load_pathway_s(chi,nchmax)
    end if

    ! Allocate arrays for expansion coefficients of 3D state
    if (realver) then
      allocate(vec3ed(msize),sliced(n23b))
    else
      allocate(vec3ez(msize),slicez(n23b))
    end if

    ! Load 3D expansion. Just reads the file with coefficients, nothing else.
    call load_expansion(vec3ed,vec3ez,ist)
    ! Load basis for FBR. DVR means DVR along phi, which is always false now.
    if (.not.dvr) then
      allocate(basis(n3,n3b))
      call init_fbrbasis(basis) ! Evaluates FBR (normalized) on phi grid
    end if

    ! Loop over slice (n-block)
    do isl=1,n1
      ! Set basis size in slice
      if (adiab) then
        bassize = blsize(isl) ! Number of 2D functions in current n-block
        if (bassize == 0) cycle
      else
        bassize = nch
      end if

      ! Load basis - all 1D and 2D energies and wave functions
      call load_eibasis(isl,nvec1,val1,vec1,val2,vec2e)
      ! Allocate array for 2D functions in DVR/FBR basis
      allocate(vec2b(n23b,bassize))
      write(LG,'(A10,I3,A,F10.3,A3)')'vec2b:', isl, ' ', sizeof(vec2b) / by2mb, ' MB'

      ! Allocate array for temporary 2D function
      allocate(vec2t(size(vec2e,1)))
      ! Loop over 2D functions in current n-block
      do k=1,bassize
        ! Get wave function of recognized state
        if (adiab) then
          vec2t = vec2e(:,k)
        else
          ichan = chi(k)
          vec2t = vec2e(:,ind(isl,ichan)) * phs(isl,ichan)
        end if
        ! Transform current 2D solution to DVR/FBR basis
        call trans_eibas_bas(vec1,nvec1,vec2t,vec2b(:,k))
      end do

      ! Compute slice
      l = (isl-1) * n23
      vol12 = jac1(isl) * alpha1 * alpha2
      if (realver) then
        ! Get vec3e chunk
        allocate(vd(bassize))
        do k=1,bassize
          if (adiab) then
            vd(k) = vec3ed(offset(isl) + k)
          else
            vd(k) = vec3ed(offset(k) + isl)
          end if
        end do

        ! Convert from eibasis to DVR/FBR basis
        call gemv(vec2b, vd, sliced)
        ! Get normalized grid function
        sliced = sliced / sqrt(vol12)
        call trans_bas_grd_d(basis,sliced,vec3d(l+1:l+n23))
        ! Deallocate chuck
        deallocate(vd)
      else
        ! Get vec3e chunk
        allocate(vz(bassize))
        do k=1,bassize
          if (adiab) then
            ! Copy expansion of current n-block
            vz(k) = vec3ez(offset(isl) + k)
          else
            vz(k) = vec3ez(offset(k) + isl)
          end if
        end do

        ! Convert from eibasis to DVR/FBR basis
        call gemv(vec2b, vz, slicez)

        ! Get normalized grid function
        slicez = slicez / sqrt(vol12)
        ! Transforms expansion over FBR to expansion over grid (DVR)
        call trans_bas_grd_z(basis,slicez,vec3z(l+1:l+n23))
        ! Deallocate chuck
        deallocate(vz)
      end if

      ! Deallocate arrays
      deallocate(vec2b,vec2t)
    end do
    ! Log that state was built
    write(LG,*)'Built state ',ist
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates basis probabilities Pa, Pb and Ps.
  !-----------------------------------------------------------------------
  subroutine calc_basprobs
    ! Probabilities
    real*8,allocatable::pall(:,:,:)    ! All probabilities
    real*8,allocatable::myp(:,:)       ! Probabilities for one slice
    real*8 p                           ! Probability
    ! Basis
    type(array1d),allocatable::val1(:) ! 1D values  for each thread
    type(array2d),allocatable::vec1(:) ! 1D vectors for each thread
    real*8,allocatable::vec2e(:,:)     ! 2D vectors in eibasis
    real*8,allocatable::vec2b(:)       ! One eigenvector in basis
    real*8,allocatable::vec2(:)        ! One eigenvector on grid
    real*8,allocatable::val2(:)        ! 2D values
    real*8,allocatable::val2all(:,:)   ! 2D values, all
    integer,allocatable::nvec1(:)      ! Num of 1D vecs in each thrd
    integer nvec2                      ! Num of 2D vecs
    real*8,allocatable::basis(:,:)     ! FBR basis
    ! Miscellaneous
    integer igr                        ! Group
    integer isl                        ! Slice number
    integer ist                        ! State number
    integer ip,ip1,ip2,pln             ! Phi loop indices
    integer it                         ! Tet loop indices
    integer i,j,k
    ! Get my slice number
    isl = myid + 1
    if (isl > n1) then
      write(LG,*)'Nothing to do, exiting'
      return
    else
      write(LG,*)'Slice number ',isl
    end if

    ! Load basis
    call load_eibasis(isl,nvec1,val1,vec1,val2,vec2e)
    nvec2 = size(vec2e,2)

    ! Load basis for FBR
    if (.not.dvr) then
      allocate(basis(n3,n3b))
      call init_fbrbasis(basis)
    end if

    ! Allocate arrays
    allocate(myp(nvec2max,ngr), vec2b(n23b), vec2(n23))
    myp = 0
    ! Set phi length, or number of points in one phi range
    pln = n3 / 6
    ! Process all states in the slice
    do ist=1,nvec2
      ! Exit loop if the state is higher than Ecut
      if (val2(ist) > trecut)exit

      ! Transform from eigenvector basis
      call trans_eibas_bas(vec1,nvec1,vec2e(:,ist),vec2b)
      ! Get normalized grid function
      call trans_bas_grd_d(basis,vec2b,vec2)
      vec2 = vec2 / sqrt(alpha2)

      ! Loop over groups
      do igr=1,ngr
        ! Get phi range
        ip1 = (igr-1) * pln + 1
        ip2 = igr * pln
        ! Calculate probability in phi range
        p = 0
        do it=1,n2
          do ip=ip1,ip2
            ! Add contribution from left half
            k = (it-1) * n3 + ip
            p = p + vec2(k)**2
            ! Add contribution from right half
            k = (it-1) * n3 + n3 - ip + 1
            p = p + vec2(k)**2
          end do
        end do
        p = p * alpha2 * alpha3
        ! Save probability
        myp(ist,igr) = p
      end do
      ! Log state processed
      write(LG,*)'Processed state ',ist,'/',nvec2
    end do

    ! Allocate collective array on root
    if (myid == 0) then
      allocate(pall(nvec2max,ngr,n1))
      pall = 0
    end if

    ! Colect results
    if (myid == 0) then
      pall(:,:,1) = myp
      do i=2,n1
        call dgerv2d(context,nvec2max,ngr,pall(1,1,i),nvec2max,0,i-1)
      end do
    else
      call dgesd2d(context,nvec2max,ngr,myp,nvec2max,0,0)
    end if

    ! Only root continues with printing
    if (myid /= 0)return

    ! Load 2D eigenvalues
    call load_val2all(val2all)
    ! Open files
    open(1,file=outdir//'/pb.out')
    open(2,file=outdir//'/pa.out')
    open(3,file=outdir//'/ps.out')
    do igr=1,ngr
      ! Write header
      write(igr,'(10X)',advance='no')
      do j=1,n1
        write(igr,'(I25)',advance='no')j
      end do
      write(igr,*)

      ! Write probabilities
      do i=1,nvec2max
        write(igr,'(I10)',advance='no')i
        do j=1,n1
          if (val2all(i,j) > trecut) then
            write(igr,'(A25)',advance='no')'N'
          else
            write(igr,'(F25.17)',advance='no')pall(i,igr,j)
          end if
        end do
        write(igr,*)
      end do
    ! Close files
    end do
    close(1)
    close(2)
    close(3)
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates spectrum using direct-product approach.
  !  Since exact direct-product is very expensive, the configuration space
  !  is reduced by excluding points with a large potential energy.
  !-----------------------------------------------------------------------
  subroutine calc_dprod
    real*8,allocatable::a(:,:)      ! Chunk of reduced matrix
    real*8,allocatable::psi(:)      ! Initial basis function
    real*8,allocatable::hpsi(:)     ! Hamilt * psi
    real*8,allocatable::val(:)      ! Resulted eigenvalues
    real*8,allocatable::vec(:,:)    ! Resulted eigenvectors
    real*8,allocatable::vecred(:,:) ! Eigenvectors in reduced form
    real*8,allocatable::idx(:,:)    ! Indices of included points
    integer np,ip                   ! Number of included points
    integer i1,i2,i3,k,istate
    character(256) fn
    real*8 stsy,s(n1)
    integer nploc,numroc,j,id,myj

    ! Calculate indices of included points
    call calc_reduction(np,idx)
    ! Calculate reduced matrix
    nploc = numroc( np, 1, myid, 0, nprocs)
    allocate(a(np,nploc),psi(nn),hpsi(nn), val(nstate),vecred(np,nstate),vec(nn,nstate))

    ! Calculate matrix blocks
    do j=1,np
      call g2l(j,np,nprocs,1,id,myj)
      if (myid==id) then
       psi = 0d0
       k = (idx(j,1)-1)*n23 + (idx(j,2)-1)*n3 + idx(j,3)
       psi(k) = 1d0
       call hamilt3D(psi,hpsi)
       do ip=1,np
         k = (idx(ip,1)-1)*n23 + (idx(ip,2)-1)*n3 + idx(ip,3)
         a(ip,myj) = hpsi(k)
       end do
      end if
    end do
    if (myid==0)write(*,*)'Reduced matrix size: ',np

    ! Calculate matrix eigen decomposition
    call scald_cd(context,np,nstate,vecred,val,a,nploc)
    vec = 0
    do istate=1,nstate
      do ip=1,np
        k = (idx(ip,1)-1)*n23 + (idx(ip,2)-1)*n3 + idx(ip,3)
        vec(k,istate) = vecred(ip,istate)
      end do
    end do

    ! Only root continues with printing
    if (myid /= 0)return

    ! Convert states
    s = sqrt(jac1 * alpha1)
    do istate=1,nstate
      do k=1,nn
        i1 = (k-1) / n23 + 1
        vec(k,istate) = vec(k,istate) / s(i1)
      end do
    end do

    ! Print results
    open(1,file='eigenvalues.out')
    do k=1,nstate
      call symmetryd(vec(:,k),stsy)
      write(1,'(I4,2F25.17)')k,val(k)*autown,stsy
    end do
    close(1)
    do istate=1,nstate
      write(fn,'(A5,I4.4,A4)')'state',istate,'.out'
      open(1,file=fn)
      do i1=1,n1
      do i2=1,n2
      do i3=1,n3
        k = i3 + (i2-1)*n3 + (i1-1)*n2*n3
        write(1,*)vec(k,istate)**2
      end do
      end do
      end do
      close(1)
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates symmetry of 3D state. Real version.
  !-----------------------------------------------------------------------
  subroutine symmetryd(psi,s3)
    real*8 psi(nn),s3
    integer i1,i2,i3,j1,j2,k1,k2
    s3 = 0
    do i1=1,n1
      do i2=1,n2
        do i3=1,n3
          j1 = i3
          j2 = n3-i3+1
          k1 = j1 + (i2-1)*n3 + (i1-1)*n2*n3
          k2 = j2 + (i2-1)*n3 + (i1-1)*n2*n3
          s3 = s3 + psi(k1)*psi(k2) * jac1(i1) * alpha1 * alpha2 * alpha3
        end do
      end do
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates symmetry of 3D state. Complex version.
  !-----------------------------------------------------------------------
  subroutine symmetryz(psi,s3)
    complex*16 psi(nn),s3
    integer i1,i2,i3,j1,j2,k1,k2
    real*8 p2
    s3 = 0
    do i1=1,n1
      do i2=1,n2
        do i3=1,n3
          j1 = i3
          j2 = n3-i3+1
          k1 = j1 + (i2-1)*n3 + (i1-1)*n2*n3
          k2 = j2 + (i2-1)*n3 + (i1-1)*n2*n3
          p2 = conjg(psi(k1)) * psi(k2)
          s3 = s3 + p2 * jac1(i1) * alpha1 * alpha2 * alpha3
        end do
      end do
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates symmetry of 2D state.
  !-----------------------------------------------------------------------
  subroutine symmetry_2d(psi,s3)
    real*8 psi(n23),s3
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
  !  Calculates symmetry of 1D state.
  !-----------------------------------------------------------------------
  subroutine symmetry_1d(psi,s3)
    real*8 psi(n3),s3
    integer i3,k1,k2
    s3 = 0
    do i3=1,n3
      k1 = i3
      k2 = n3-i3+1
      s3 = s3 + psi(k1)*psi(k2) * alpha3
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Prints results of 3SDT.
  !-----------------------------------------------------------------------
  subroutine prnt_3dsdt(val3d,vec3d,val3z,vec3z,nstloc,np)
    real*8,allocatable     :: val3d(:)
    real*8,allocatable     :: vec3d(:,:)
    complex*16,allocatable :: val3z(:)
    complex*16,allocatable :: vec3z(:,:)
    real*8,allocatable     :: pb(:,:)
    integer, parameter     :: npb = 3
    real*8  baren(ngr)
    real*8  barps(ngr)
    integer ist,myist
    integer nstloc ! num of states in this process
    integer np ! num of processes
    integer isl
    integer ich
    integer k
    character(256) fn

    ! Log spectrum
    write(LG,*)'Calculated spectrum:'
    if (realver) then
      do ist=1,nstate
        write(LG,'(I4, F25.17)')ist,val3d(ist) * autown
      end do
    else
      do ist=1,nstate
        write(LG,'(I4,2F25.17)')ist,val3z(ist) * autown
      end do
    end if

    ! Load lowest channel in pathway S and setup arrays
    if (.not.adiab) then
      allocate(pb(nstate,npb))
      pb = 0
      call load_lowest_barrier(barps,baren)
    end if

    ! Write each state in a separate binary file
    do myist=1,nstloc
      ! Get global state number
      call l2g(myist,myid,nstate,np,1,ist)
      ! Build file name
      write(fn,'(4A,I4.4,A)')outdir,'/',expdir,'/exp.',ist,'.bin.out'

      ! Write data
      open(1,file=fn,form='unformatted')
      if (realver) then
        write(1)vec3d(:,myist)
      else
        write(1)vec3z(:,myist)
      end if
      close(1)

      ! Do not compute probabilities in adiab or real mode
      if (adiab .or. realver)cycle

      ! Calculate probabilities
      do isl=1,n1
        do ich=1,nch
          k = (ich-1) * n1 + isl
          ! Real part only
          pb(ist,2) = pb(ist,2) + real(vec3z(k,ist))**2
          ! Imaginary part only
          pb(ist,3) = pb(ist,3) + aimag(vec3z(k,ist))**2
          ! Probability in the well
          if (g1(isl) < barps(3)) then
            pb(ist,1) = pb(ist,1) + conjg(vec3z(k,myist)) * vec3z(k,myist)
          end if
        end do
      end do
    end do

    ! Gather probabilities if not adiab
    if (solver /= SOLVER_LAPA .and. .not.adiab) then
      call dgsum2d(context,'A',' ',nstate,npb,pb,nstate,0,0)
    end if

    ! Exit if not a root
    if (myid /= 0)return

    ! Write to file
    write(fn,'(2A)')outdir,'/spec.out'
    open(1,file=fn)
    if (realver) then
      do ist=1,nstate
        write(1,'(I4,F30.17)')ist,val3d(ist) * autown
      end do
    else
      do ist=1,nstate
        if (adiab) then
          write(1,'(I4,2F30.17)')ist, real( val3z(ist)) * autown, aimag(val3z(ist)) * autown * (-2) ! Imaginary part is multiplied by -2 to convert it to gamma (resonance width)
        else
          write(1,'(I4,6F30.17)')ist, real( val3z(ist)) * autown, real( val3z(ist)) * autown - baren(3), aimag(val3z(ist)) * autown * (-2), pb(ist,:)
        end if
      end do
    end if
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads expansion.
  !-----------------------------------------------------------------------
  subroutine load_expansion(vec3d,vec3z,is)
    real*8,allocatable     :: vec3d(:)
    complex*16,allocatable :: vec3z(:)
    integer is
    character(256) fn

    ! Build file name
    write(fn,'(6A,I4.4,A)') dpath, '/', getdir(MODE_3DSDT),'/',expdir,'/exp.',is,'.bin.out'
    ! Read data
    open(1,file=fn,form='unformatted')
    if (realver) then
      read(1)vec3d
    else
      read(1)vec3z
    end if
    close(1)
  end subroutine
  
  !-----------------------------------------------------------------------
  !  Loads matrix partitioning.
  !-----------------------------------------------------------------------
  subroutine load_partitioning(nbl,blsize,blsizemax,offset,msize)
    integer,allocatable::blsize(:)
    integer,allocatable::offset(:)
    integer nbl
    integer blsizemax
    integer msize
    integer i

    ! Setup block structure
    if(adiab)then
      nbl = n1
      call load_nvec2(blsize)
      blsizemax = nvec2max
    else
      nbl = nch
      allocate(blsize(nch))
      blsize = n1
      blsizemax = n1
    endif

    ! Calculate offsets
    allocate(offset(nbl))
    offset(1) = 0
    do i=2,nbl
      offset(i) = offset(i-1) + blsize(i-1)
    enddo

    ! Calculate matrix size
    msize = sum(blsize)
  end subroutine
       
  !-----------------------------------------------------------------------
  !  Loads number of 2D states only.
  !-----------------------------------------------------------------------
  subroutine load_nvec2(nvec2)
    integer,allocatable::nvec2(:)
    character(:),allocatable::ldir
    integer i
    allocate(nvec2(n1))
    ldir = bpath // '/' // getdir(MODE_BASIS)
    open(1,file=ldir//'/nvec2.dat')
    do i=1,n1
      read(1,'(10X,I10)')nvec2(i)
    enddo
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads all 2D eigenvalues.
  !-----------------------------------------------------------------------
  subroutine load_val2all(val2)
    real*8, allocatable::val2(:,:)
    real*8, allocatable::val2s(:)
    integer isl
    allocate(val2(nvec2max,n1))
    val2 = 0
    do isl=1,n1
      call load_val2(isl,val2s)
      val2(1:size(val2s),isl) = val2s
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads results of recognition.
  !-----------------------------------------------------------------------
  subroutine load_recognition(ind,phs,nchan)
    integer,allocatable::ind(:,:)
    integer,allocatable::phs(:,:)
    character(:),allocatable::ldir
    integer nchan
    integer ichan
    integer i,isl,ios

    ! Get channel recognition directory
    ldir = getrecldpath()
    ! Get number of channels (number of lines in file)
    open(1,file=ldir//'/channels.dat')
    nchan = 0
    do while(.true.)
      read(1,*,iostat=ios)
      if (ios /= 0) then
        exit
      else
        nchan = nchan + 1
      end if
    end do
    close(1)

    ! Allocate arrays
    allocate(ind(n1,nchan),phs(n1,nchan))
    ! Open indices and phase files and skip headers
    open(1,file=ldir//'/recind.dat')
    open(2,file=ldir//'/recphs.dat')
    read(1,*)
    read(2,*)

    ! Read data
    do ichan=1,nchan
      read(1,'(10X)',advance='no')
      read(2,'(10X)',advance='no')
      do isl=1,n1
        read(1,'(I25)',advance='no')ind(isl,ichan)
        read(2,'(I25)',advance='no')phs(isl,ichan)
      end do
      read(1,*)
      read(2,*)
    end do
    close(1)
    close(2)
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads results of recognition for a specific channel.
  !-----------------------------------------------------------------------
  function load_recognition_chan(ind,phs,ichan) result(res)
    logical res
    integer,allocatable::ind(:)
    integer,allocatable::indall(:,:)
    integer,allocatable::phs(:)
    integer,allocatable::phsall(:,:)
    integer ichan
    integer nchan
    call load_recognition(indall,phsall,nchan)
    if (ichan > nchan) then
      res = .false.
    else
      allocate(ind(n1),phs(n1))
      ind = indall(:,ichan)
      phs = phsall(:,ichan)
      res = .true.
    end if
  end function

  !-----------------------------------------------------------------------
  !  Loads channel numbers for pathway S.
  !-----------------------------------------------------------------------
  subroutine load_pathway_s(chs,nchmax)
    integer,allocatable::chs(:)
    character(:),allocatable::ldir
    integer ichan
    integer nchmax
    integer igr
    integer ios

    ! Get channel recognition directory
    ldir = getrecldpath()
    ! Allocate arrays
    allocate(chs(nch))
    chs = 0
    nchmax = 0

    ! Read channel grouping
    open(1,file=ldir//'/channels.dat')
    do while(.true.)
      read(1,*,iostat=ios)ichan,igr
      if (ios /= 0) then
        exit
      else if (igr /= 3) then
        cycle
      else
        nchmax = nchmax + 1
        chs(nchmax) = ichan
        if (nchmax == nch)exit
      end if
    end do
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  !  Writes 2D vector.
  !-----------------------------------------------------------------------
  subroutine writ_vec2(vec,i1,is,sgn)
    real*8, allocatable::vec(:)
    integer i1,is,sgn,i,j,l
    character(256) fn
    write(fn,'(2A,I4.4,A,I4.4,A)')outdir,'/vec2.',i1,'.',is,'.out'
    open(1,file=fn)
    do i=1,n2
      do j=1,n3
        l = (i-1)*n3 + j
        write(1,'(F25.17)',advance='no')vec(l) * sgn
      end do
      write(1,*)
    end do
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads 2D vector in basis repr or on grid.
  !-----------------------------------------------------------------------
  subroutine load_vec2(vec2,isl,ist,ongrid)
    integer isl,ist
    logical ongrid
    ! 1D basis
    type(array1d),allocatable::val1(:) ! 1D values  for each thread
    type(array2d),allocatable::vec1(:) ! 1D vectors for each thread
    integer,allocatable::nvec1(:)      ! Num of 1D vectors
    ! 2D basis
    real*8,allocatable::val2(:)        ! 2D values
    real*8,allocatable::vec2e(:,:)     ! 2D vectors in eibasis repr
    real*8,allocatable::vec2b(:)       ! 2D vector  in   basis repr
    real*8,allocatable::vec2(:)        ! 2D vector  for output
    real*8,allocatable::basis(:,:)     ! FBR basis

    ! Deallocate if already allocated
    if (allocated(vec2)) deallocate(vec2)

    ! Allocate basis repr
    allocate(vec2b(n23b))
    ! Load eibasis
    call load_eibasis(isl,nvec1,val1,vec1,val2,vec2e)
    ! Transform from eigenvector basis
    call trans_eibas_bas(vec1,nvec1,vec2e(:,ist),vec2b)
    ! Transform to grid if requested
    if (ongrid) then
      ! Load basis for FBR
      if (.not.dvr) then
        allocate(basis(n3,n3b))
        call init_fbrbasis(basis)
      end if

      ! Get normalized grid function
      allocate(vec2(n23))
      call trans_bas_grd_d(basis,vec2b,vec2)
      vec2 = vec2 / sqrt(alpha2)
    ! Or return basis repr
    else
      allocate(vec2(n23b))
      vec2 = vec2b
    end if
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads 2D vector in basis repr or on grid.
  !-----------------------------------------------------------------------
  function getdir(m) result(res)
    character(:),allocatable :: res
    integer   m ! Mode
    character c ! Character
    res = trim(dirs(m))
    if (m == MODE_OVERLAP) then
      write(c,'(I1)')oltype
      res = res // c
    end if
  end function

  !-----------------------------------------------------------------------
  !  Creates a sub directory.
  !-----------------------------------------------------------------------
  subroutine subdir(dir)
    character(*) dir
    call execute_command_line('mkdir -p ' // outdir // '/' // dir)
  end subroutine

  !-----------------------------------------------------------------------
  !  Loads properties of the lowest channel.
  !-----------------------------------------------------------------------
  subroutine load_lowest_barrier(barps,baren)
    character(:),allocatable::ldir
    real*8  barps(ngr)
    real*8  baren(ngr)
    integer found(ngr)
    integer igr
    integer isl
    real*8  ps,en

    ! Get channel recognition directory
    ldir = getrecldpath()
    ! Find lowest channel in channels file
    found = 0
    open(1,file=ldir//'/channels.dat')
    do
      read(1,*)igr,igr,ps,en
      if (found(igr) == 0) then
        found(igr) = 1
        barps(igr) = ps
        baren(igr) = en
      end if
      if (sum(found) == 3)exit
    end do
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  !  Returns recognition load directory path.
  !-----------------------------------------------------------------------
  function getrecldpath() result(res)
    character(:),allocatable::res
    if (recld == RECLD_DEF) then
      res = rpath // '/' // getdir(MODE_CHRECOG)
    else
      res = rpath // '/reccontr'
    end if
  end function
  
  subroutine init_caps(capebarin)
    implicit none
    real*8 capebarin
    character(256)fn
    capebar = capebarin
    call calc_cap(n1,g1,LG)
    write(fn,'(4A,I5.5,A)')outdir,'/',capdir,'/cap',myid+1,'.out'
    call prnt_cap(fn)
  end subroutine

  !-----------------------------------------------------------------------
  !  Provides CAP to sdt for integration.
  !-----------------------------------------------------------------------
  subroutine get_cap(p)
    implicit none
    real*8 p(n1)
    p = cap(:,capid)
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates 1D Hamiltonian for coordinate #1. Real version.
  !-----------------------------------------------------------------------
  subroutine init_matrix1d(ham)
    implicit none
    real*8 ham(n1,n1)
    real*8 L
    integer i,j

    ! FFT
    if (ham1type == HAM1_FFT) then
      ham = - der1 / (2d0 * mu) ! /tools/pesgeneral.f90
    ! Finite differences
    else
      L = n1 * alpha1
      do i=1,n1
        do j=1,n1
          if (i == j) then 
            ham(i,j) = pi**2 / (mu * L**2) * (n1**2 + 2) / 6 / jac1(i)**2
          else
            ham(i,j) = pi**2 / (mu * L**2 * 2) * (-1)**(i-j) / (sin((i-j)*pi/n1))**2 * (1/jac1(i)**2 + 1/jac1(j)**2)
          endif
        enddo
      enddo
    endif
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates 1D Hamiltonian for coordinate #1. Complex version.
  !-----------------------------------------------------------------------
  subroutine init_matrix1z(ham)
    implicit none
    complex*16 ham(n1,n1)
    real*8 tmp(n1,n1)
    integer i

    ! Variables for jac potential
    real*8 fd
    real*8 sd
    real*8 t

    ! Complex fft
    if (ham1type == HAM1_FFT_CMPL) then
      ham = - der1z / (2d0 * mu)
    else
      call init_matrix1d(tmp)
      ham = tmp
    endif

    ! Add complex potential
    if (capid /= 0) then
      do i=1,n1
        ham(i,i) = ham(i,i) + (0,-1) * cap(i,capid)
      enddo
    endif

    ! Exit if FFT
    if (ham1type == HAM1_FFT .or. ham1type == HAM1_FFT_CMPL) then
      return
    end if

    ! If not FFT, then add potential term due to grid jacobian
    do i=1,n1
      ! Calculate 1st derivative
      if (i == 1) then
        fd = jac1(i)   * (-1     ) + jac1(i+1) * ( 1     )
      elseif (i == n1) then
        fd = jac1(i-1) * (-1     ) + jac1(i)   * ( 1     )
      elseif (i == 2 .or. i == n1-1) then
        fd = jac1(i-1) * (-1d0/2 ) + jac1(i)   * (     0 ) + jac1(i+1) * ( 1d0/2 )
      elseif (i == 3 .or. i == n1-2) then
        fd = jac1(i-2) * ( 1d0/12) + jac1(i-1) * (-2d0/3 ) + jac1(i)   * (     0 ) + jac1(i+1) * ( 2d0/3 ) + jac1(i+2) * (-1d0/12)
      elseif (ham1type == HAM1_ANALYTIC6 .or. (ham1type == HAM1_ANALYTIC8 .and. (i == 4 .or. i == n1-3))) then
        fd = jac1(i-3) * (-1d0/60) + jac1(i-2) * ( 3d0/20) + jac1(i-1) * (-3d0/4 ) + jac1(i)   * (     0 ) + jac1(i+1) * ( 3d0/4 ) + jac1(i+2) * (-3d0/20) + jac1(i+3) * ( 1d0/60)
      elseif (ham1type == HAM1_ANALYTIC8) then
        fd = jac1(i-4) * ( 1d0/280) + jac1(i-3) * (-4d0/105) + jac1(i-2) * (1d0/5) + jac1(i-1) * (-4d0/5) + jac1(i) * (0) + jac1(i+1) * (4d0/5) + jac1(i+2) * (-1d0/5) + jac1(i+3) * (4d0/105) + jac1(i+4) * (-1d0/280)
      endif
      fd = fd / alpha1

      ! Calculate 2nd derivative
      if (i == 1) then
        sd = jac1(i)   * ( 1     ) + jac1(i+1) * (-2     ) + jac1(i+2) * ( 1     )
      elseif (i == n1) then
        sd = jac1(i-2) * ( 1     ) + jac1(i-1) * (-2     ) + jac1(i)   * ( 1     )
      elseif (i == 2 .or. i == n1-1) then
        sd = jac1(i-1) * ( 1     ) + jac1(i)   * (-2     ) + jac1(i+1) * ( 1     )
      elseif (i == 3 .or. i == n1-2) then
        sd = jac1(i-2) * (-1d0/12) + jac1(i-1) * ( 4d0/3 ) + jac1(i)   * (-5d0/2 ) + jac1(i+1) * ( 4d0/3 ) + jac1(i+2) * (-1d0/12)
      elseif (ham1type == HAM1_ANALYTIC6 .or. (ham1type == HAM1_ANALYTIC8 .and. (i == 4 .or. i == n1-3))) then
        sd = jac1(i-3) * ( 1d0/90) + jac1(i-2) * (-3d0/20) + jac1(i-1) * ( 3d0/2 ) + jac1(i)   * (-49d0/18 ) + jac1(i+1) * ( 3d0/2 ) + jac1(i+2) * (-3d0/20) + jac1(i+3) * ( 1d0/90)
      elseif (ham1type == HAM1_ANALYTIC8) then
        sd = jac1(i-4) * (-1d0/560) + jac1(i-3) * (8d0/315) + jac1(i-2) * (-1d0/5) + jac1(i-1) * (8d0/5) + jac1(i) * (-205d0/72) + jac1(i+1) * (8d0/5) + jac1(i+2) * (-1d0/5) + jac1(i+3) * (8d0/315) + jac1(i+4) * (-1d0/560)
      endif
      sd = sd / alpha1**2
      t = 7d0/4 * fd**2 / jac1(i)**4 - sd / jac1(i)**3 / 2
      ham(i,i) = ham(i,i) + t / (2d0 * mu)
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates 1D Hamiltonian for coordinate #2.
  !-----------------------------------------------------------------------
  subroutine init_matrix2(ham,ic1)
    implicit none
    real*8 ham(n2,n2)
    integer ic1
    real*8 coeff
    real*8 L
    integer i,j

    ! FFT
    if (ham2type == HAM2_FFT) then
      coeff = - 1 / (2d0 * mu) * 4 / grho2(ic1)
      ham = der2 * coeff
    ! Analytic
    else
      L = n2 * alpha2
      coeff = pi**2 / (mu * L**2) * 4 / grho2(ic1)
      do i=1,n2
        ham(i,i) = coeff * (n2**2 + 2) / 6d0
        do j=i+1,n2
          ham(i,j) = (-1)**(i-j) * coeff / sin((i-j) * pi / n2)**2
          ham(j,i) = ham(i,j)
        end do
      end do
    endif
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates 1D Hamiltonian for coordinate #3 using DVR.
  !-----------------------------------------------------------------------
  subroutine init_matrix3dvr(ham,ic1,ic2)
    implicit none
    real*8 ham(n3,n3)
    integer ic1,ic2,ic3
    real*8 coeff
    coeff = - 1 / (2d0 * mu) * 4 / grho2(ic1) / sintet2(ic2)
    ham = der3 * coeff
    do ic3=1,n3
      ham(ic3,ic3) = ham(ic3,ic3) + pottot(ic3,ic2,ic1)
    enddo
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates 1D Hamiltonian for coordinate #3 using FBR.
  !-----------------------------------------------------------------------
  subroutine init_matrix3fbr(ham,ic1,ic2)
    implicit none
    real*8 ham(n3b,n3b)
    real*8 basis(n3,n3b)
    integer ic1,ic2,ic3
    integer i,j,k,np
    real*8 coeff,sum

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

    ! Symmetries E1 and E2: exclude 3-fold basis functions
    if (sy==SY_E1.or.sy==SY_E2) then
      k = 0
      do i=1,n3b
        k = k + 1
        if (mod(k,3)==0)k = k + 1
        ham(i,i) = ham(i,i) + coeff * (-1) * k**2
      enddo
    ! Other symmetries: use 1-fold or 3-fold basis functions
    else
      ! Periods coefficient
      np = 3
      if (sy==SY_SY.or.sy==SY_AS)np = 1

      ! Build kinetic energy matrix
      j = 0
      if (sy==SY_A1.or.sy==SY_SY)j = 1
      do i=1,n3b
        ham(i,i) = ham(i,i) + coeff * (-1) * np**2 * (i-j)**2
      enddo
    endif
  end subroutine

  !-----------------------------------------------------------------------
  !  Applies 3D hamiltonian operator to 3D wave function.
  !-----------------------------------------------------------------------
  subroutine hamilt3D(psi,hpsi)
    implicit none
    real*8 psi(nn),hpsi(nn)
    real*8 hpsi1(n1,n3,n2)
    real*8 hpsi2(n2,n3,n1)
    real*8 hpsi3(n3,n2,n1)
    integer i1,i2,i3,k
    hpsi1 = 0
    hpsi2 = 0
    hpsi3 = 0
    do i1=1,n1
      do i2=1,n2
        do i3=1,n3
          k = i3 + (i2-1)*n3 + (i1-1)*n2*n3
          hpsi1(i1,i3,i2) = psi(k)
          hpsi2(i2,i3,i1) = psi(k)
          hpsi3(i3,i2,i1) = psi(k)
        enddo
      enddo
    enddo
    call calc_derivd_jac_2nd(n1,n2*n3,hpsi1,freq1,jac1) ! general_vars
    call calc_derivd(2,n1*n3,n2,freq2,hpsi2)
    call calc_derivd(2,n1*n2,n3,freq3,hpsi3)
    do i1=1,n1
      do i2=1,n2
        do i3=1,n3
          k = i3 + (i2-1)*n3 + (i1-1)*n2*n3
          hpsi(k) =  hpsi1(i1,i3,i2) + 4/grho2(i1) * ( hpsi2(i2,i3,i1) + hpsi3(i3,i2,i1) / sintet2(i2) )
          hpsi(k) = -hpsi(k)/(2.0d0*mu) + psi(k)*pottot(i3,i2,i1)
          if (capid /= 0)hpsi(k)= hpsi(k) + psi(k)*(0,-1)*cap(i1,capid)
        enddo
      enddo
    enddo
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates second derivative on optimal grid. Real
  !-----------------------------------------------------------------------
  subroutine calc_derivd_jac_2nd(ni,nj,psi,freq,jac)
    implicit none
    integer ni,nj,i,j
    real*8 psi(ni,nj),freq(ni),jac(ni),sqrtjac(ni)
    
    sqrtjac = sqrt(jac)
    do j=1,nj
      psi(:,j) = psi(:,j) / sqrtjac
    enddo
    call calc_derivd(1,nj,ni,freq,psi)
    do j=1,nj
      psi(:,j) = psi(:,j) / jac
    enddo
    call calc_derivd(1,nj,ni,freq,psi)
    do j=1,nj
      psi(:,j) = psi(:,j) / sqrtjac
    enddo
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates second derivative on optimal grid. Complex
  !-----------------------------------------------------------------------
  subroutine calc_derivz_jac_2nd(ni,nj,psi,freq,jac)
    implicit none
    integer    ni,nj,i,j
    complex*16 psi(ni,nj),freq(ni)
    real*8     jac(ni),sqrtjac(ni)
    
    sqrtjac = sqrt(jac)
    do j=1,nj
      psi(:,j) = psi(:,j) / sqrtjac
    enddo
    call calc_derivz(nj,ni,freq,psi)
    do j=1,nj
      psi(:,j) = psi(:,j) / jac
    enddo
    call calc_derivz(nj,ni,freq,psi)
    do j=1,nj
      psi(:,j) = psi(:,j) / sqrtjac
    enddo
  end subroutine

  !-----------------------------------------------------------------------
  !  Finds indices of reduced configuration space for direct product.
  !-----------------------------------------------------------------------
  subroutine calc_reduction(np,idx)
    implicit none
    integer np,ip,i1,i2,i3
    real*8,allocatable::idx(:,:)        ! Indices of included points
    
    np = 0
    do i1=1,n1
      do i2=1,n2
        do i3=1,n3
          if (pottot(i3,i2,i1)<trecut) np = np + 1
        enddo
      enddo
    enddo
    allocate(idx(np,3))
    ip = 0
    do i1=1,n1
      do i2=1,n2
        do i3=1,n3
          if (pottot(i3,i2,i1)<trecut) then
            ip = ip + 1
            idx(ip,1) = i1
            idx(ip,2) = i2
            idx(ip,3) = i3
          endif
        enddo
      enddo
    enddo
  end subroutine
end module
