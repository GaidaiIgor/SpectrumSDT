module parpack
  use algorithms
  use index_conversion_mod
  use matmul_operator_mod
  use parallel_utils

contains

  !-----------------------------------------------------------------------
  !  This file contains drivers to PARPACK (parpack.f)
  !  Calculates eigenpairs, adopted from sample code.
  !  Author: Alexander Teplukhin
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !  PARPACK eigenproblem solver, real verison.
  !-----------------------------------------------------------------------
  subroutine pard(comm,eivals,eivecs,n,nloc,nev,nevloc,ncv,maxitr,op)
    implicit none
    include 'debug.h'
    ! include 'stat.h'

    ! Input variables
    integer   n, nloc, nev, nevloc, ncv, maxitr
    procedure(matmul_operator_real) :: op
    ! Output arrays
    real*8    eivals(nev)
    real*8    eivecs(n,nevloc)
    ! BLACS variables
    integer   comm, iam, nprocs
    integer   nprow, npcol, myrow, mycol
    ! Parpack variables
    logical   rvec
    character bmat*1, which*2
    integer   ido, lworkl, info, ierr, j, nconv, ldv
    integer   iparam(11), ipntr(11)
    real*8    tol, sigma, pdnorm2
    ! Parpack arrays
    real*8,allocatable::resid(:)
    real*8,allocatable::v(:,:)
    real*8,allocatable::workd(:)
    real*8,allocatable::workl(:)
    real*8,allocatable::d(:,:)
    real*8,allocatable::ax(:)
    logical,allocatable::select(:)
    ! Variables for eigenvectors reorganization
    real*8,allocatable::bv(:,:)
    integer bproc
    integer bnloc
    integer offset
    integer myj

    ! Get BLACS information
    call blacs_gridinfo(comm,nprow,npcol,myrow,mycol)
    ! Exit if not in the grid
    if (myrow == -1)return

    ! Check that grid is a row of processes
    if (nprow /= 1) then
      if (mycol == 0)write(*,*)'Process grid is not a row'
      return
    end if

    ! Set process number and number of processes
    iam    = mycol
    nprocs = npcol

    ! Allocate arrays
    ldv    = nloc
    lworkl = ncv*(ncv+8)
    allocate(resid(nloc), v(ldv,ncv), workd(3*nloc), workl(lworkl), d(nev,2), ax(nloc), select(ncv))

    ! Set Parpack variables
    ndigit = -3
    logfil = 6
    msaupd = 1
    msaup2 = 1
    bmat = 'I'
    which = 'SA'
    tol = 0
    info = 0
    ido = 0
    iparam(1) = 1
    iparam(3) = maxitr 
    iparam(7) = 1

    ! Reverse communication loop
    do
      call pdsaupd ( comm, ido, bmat, nloc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
      if (ido /= -1 .and. ido /= 1) exit
      call op(nloc, workd(ipntr(1)), workd(ipntr(2)))
    end do

    ! Check for errors
    if (info < 0) then
      if (iam == 0) then
        print *, ' '
        print *, ' Error with _saupd, info = ', info
        print *, ' Check documentation in _saupd '
        print *, ' '
      end if
      return
    end if

    ! Post-processing
    rvec = .true.
    call pdseupd ( comm, rvec, 'All', select, d, v, ldv, sigma, bmat, nloc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr )

    ! Check for errors
    if (ierr /= 0) then
      if (iam == 0) then
        print *, ' '
        print *, ' Error with _seupd, info = ', ierr
        print *, ' Check the documentation of _seupd. '
        print *, ' '
      end if
      return
    end if

    ! Print converged eigenvalues and residuals
    nconv = iparam(5)
    do j=1,nconv
      call op(nloc, v(1,j), ax)
      call daxpy(nloc, -d(j,1), v(1,j), 1, ax, 1)
      d(j,2) = pdnorm2( comm, nloc, ax, 1 )
    end do
    call pdmout(comm, 6, nconv, 2, d, nev, -6, 'Ritz values and direct residuals')

    ! Print additional convergence information
    if (iam == 0) then
      if (info == 1) then
        print *, ' '
        print *, ' Maximum number of iterations reached.'
        print *, ' '
      elseif (info == 3) then
        print *, ' ' 
        print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
        print *, ' '
      end if      
      print *, ' '
      print *, '_SDRV1 '
      print *, '====== '
      print *, ' '
      print *, ' Size of the matrix is ', n
      print *, ' The number of processors is ', nprocs
      print *, ' The number of Ritz values requested is ', nev
      print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
      print *, ' What portion of the spectrum: ', which
      print *, ' The number of converged Ritz values is ', nconv 
      print *, ' The number of Implicit Arnoldi update', ' iterations taken is ', iparam(3)
      print *, ' The number of OP*x is ', iparam(9)
      print *, ' The convergence criterion is ', tol
      print *, ' '
    end if

    ! Save eigenvalues
    eivals = d(:,1)
    ! Collect eigenvectors
    ! Each process broadcasts his chunk of eigenvectors.
    ! Then all processes save columns they only are responsible for.
    offset = 0
    do bproc=0,nprocs-1
      ! Bcast nloc
      if (bproc == iam) then
        bnloc = nloc
        call igebs2d(comm, 'A', ' ', 1, 1, bnloc, 1)
      else
        call igebr2d(comm, 'A', ' ', 1, 1, bnloc, 1, 0, bproc)
      end if

      ! Allocate space for chunk
      allocate(bv(bnloc,nev))
      ! Bcast chunk
      if (bproc == iam) then
        bv = v(:,1:nev)
        call dgebs2d(comm, 'A', ' ', bnloc, nev, bv, bnloc)
      else
        call dgebr2d(comm, 'A', ' ', bnloc, nev, bv, bnloc, 0, bproc)
      end if

      ! Save needed columns from chunk
      do myj=1,nevloc
        call l2g(myj,iam,nev,npcol,1,j)
        eivecs(offset+1:offset+bnloc,myj) = bv(:,j)
      end do

      ! Deallocate chunk
      deallocate(bv)
      ! Update offset
      offset = offset + bnloc
    end do
  end subroutine

  !-----------------------------------------------------------------------
  ! PARPACK eigenproblem solver, complex verison.
  !-----------------------------------------------------------------------
  subroutine parz(comm,eivals,eivecs,n,nloc,nev,nevloc,ncv,maxitr)
    implicit none

    ! Input variables. 
    ! n - full matrix size
    ! nloc - size of hamiltonian chunk belonging to this processor (in rows)
    ! nev - number of states requested
    ! nevloc - number of states belonging to this processor
    ! ncv - arnoldi basis size
    ! maxitr - maximum number of iterations
    integer   n, nloc, nev, nevloc, ncv, maxitr
    ! Output arrays
    complex*16 eivals(nev)
    complex*16 eivecs(n,nevloc) ! each process has only some columns of the overall eigenvectors matrix, distributed in a round-robin fashion (cyclically)
    ! BLACS variables
    integer    comm, iam, nprocs
    integer    nprow, npcol, myrow, mycol ! total number of processors in the grid along rows and columns; this processor's row and column
    ! Parpack variables
    logical    rvec
    character  bmat*1, which*2
    integer    ido, lworkl, info, ierr, j, nconv, ldv
    integer    iparam(11), ipntr(14) ! iparam contains algorithm parameters; ipntr contains starting indexes of different data parts in workd and workl
    real*8     tol, pdznorm2
    complex*16 sigma
    ! Parpack arrays
    complex*16,allocatable::resid(:) ! residual vector
    complex*16,allocatable::v(:,:) ! arnoldi basis matrix / Ritz vectors
    complex*16,allocatable::workd(:) ! communication array
    complex*16,allocatable::workl(:) ! some internal storage for parpack ?
    complex*16,allocatable::d(:) ! Ritz values (eigenvalues)
    complex*16,allocatable::ax(:)
    complex*16,allocatable::workev(:) ! some internal storage for parpack ?
    real*8,    allocatable::rwork(:) ! some internal storage for parpack ?
    real*8,    allocatable::rd(:,:) ! used to store results and residual norms
    logical,   allocatable::select(:)
    ! Variables for eigenvectors reorganization
    complex*16,allocatable::bv(:,:)
    real*8  s(nev)
    integer, allocatable :: ind(:)
    integer bproc
    integer bnloc ! nloc of currently broadcasting processor
    integer offset
    integer myj
    ! Variables to fix eigenvector phase
    integer i
    real*8 maxv
    complex*16 c
    
    call print_parallel('Starting diagonalization')

    ! Get BLACS information
    call get_grid_info(comm, nprow, npcol, myrow, mycol)
    ! Exit if not in the grid
    if (myrow == -1) then
      print *, 'Not in grid, exiting'
      return
    end if

    ! Check that grid is a row of processes
    if (nprow /= 1) then
      if (mycol == 0) then
        print *, 'Process grid is not a row'
      end if
      return
    end if

    ! Set process number and number of processes
    iam    = mycol
    nprocs = npcol

    ! Allocate arrays
    ldv    = nloc ! leading dimension of v
    lworkl = ncv*(3*ncv+5) ! length of workl
    allocate(resid(nloc), v(ldv,ncv), workd(3*nloc), workl(lworkl), d(nev), ax(nloc), workev(3*ncv), rwork(ncv), rd(nev,3), select(ncv))
    
    bmat = 'I' ! indicates standard eigenvalue problem A*x=lambda*x
    which = 'SR' ! smallest real part of spectrum
    tol = 0 ! convergence precision up to machine epsilon
    info = 0 ! use random initial residual vector. on output contains status of operation
    ido = 0 ! type of operation to be performed
    iparam(1) = 1 ! automatic shifts
    iparam(3) = maxitr ! max number of iterations
    iparam(7) = 1 ! specifies standard eigenvalue problem Ax=ax

    ! Reverse communication loop
    do
      call pznaupd(comm, ido, bmat, nloc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)
      if (ido /= -1 .and. ido /= 1) then
        exit
      end if

      call active_matmul_operator(nloc, workd(ipntr(1)), workd(ipntr(2)))
    end do

    ! Check for errors
    if (info < 0) then
      if (iam == 0) then
        print *, ' '
        print *, ' Error with _naupd, info = ', info
        print *, ' Check documentation in _naupd '
        print *, ' '
      end if
      return
    end if

    ! Post-processing
    rvec = .true. ! whether to compute eigenvectors in addition to eigenvalues
    ! select is not used since we chose to compute All vectors
    ! sigma is used only when iparam(7) = 3
    call pzneupd ( comm, rvec, 'All', select, d, v, ldv, sigma, workev, bmat, nloc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, ierr )

    ! Check for errors
    if (ierr /= 0) then
      if (iam == 0) then
        print *, ' '
        print *, ' Error with _neupd, info = ', ierr
        print *, ' Check the documentation of _neupd. '
        print *, ' '
      end if
      return
    end if

    ! Print converged eigenvalues and residuals
    nconv = iparam(5)
    do j=1,nconv
      call active_matmul_operator(nloc, v(1,j), ax)
      ! computes Y =     A   *  X   +    Y (residual vector)
      call zaxpy(nloc, -d(j), v(1,j), 1, ax, 1) ! 1 and 1 are increments in X and Y respectively
      rd(j,1) = dble(d(j)) ! real part of a complex number
      rd(j,2) = dimag(d(j)) ! imaginary part
      rd(j,3) = pdznorm2(comm, nloc, ax, 1) ! norm of residual vector
    end do
    ! prints rd to file, first 6 - out file unit, second 6 - num of digits to print
    call pdmout(comm, 6, nconv, size(rd, 2), rd, nev, 6, 'Ritz values (Real, Imag) and direct residuals')

    ! Fix phase: rotate the functions in complex plane to align largest component of each function with real axis
    do j=1,nconv
      maxv = 0 ! max abs value of j-th vector (among available rows)
      do i=1,nloc
        if (abs(v(i,j)) > maxv) then ! abs works in the sense of 2-norm of complex vector
          maxv = abs(v(i,j))
          c = conjg(v(i,j)) / maxv ! normalized complex-conjugate of the largest component of vector. c * v(i, j) makes v(i, j) real while preserving modulus of v(i, j)
        end if
      end do
      ! Compares maxv between processes. The maximum value is stored in maxv
      ! i and bproc will contain row and column index of the process that had the maximum
      call dgamx2d(comm,'A',' ',1,1,maxv,1,i,bproc,1,-1,-1)
      if (bproc == iam) then
        ! sends c to other processes
        call zgebs2d(comm, 'A', ' ', 1, 1, c, 1)
      else
        ! receives c from other processes
        call zgebr2d(comm, 'A', ' ', 1, 1, c, 1, 0, bproc)
      end if
      v(:,j) = c * v(:,j)
    end do

    ! Print additional convergence information
    if (iam == 0) then
      if (info == 1) then
        print *, ' '
        print *, ' Maximum number of iterations reached.'
        print *, ' '
      elseif (info == 3) then
        print *, ' ' 
        print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
        print *, ' '
      end if      
      print *, ' '
      print *, '_NDRV1 '
      print *, '====== '
      print *, ' '
      print *, ' Size of the matrix is ', n
      print *, ' The number of processors is ', nprocs
      print *, ' The number of Ritz values requested is ', nev
      print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
      print *, ' What portion of the spectrum: ', which
      print *, ' The number of converged Ritz values is ', nconv 
      print *, ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
      print *, ' The number of OP*x is ', iparam(9)
      print *, ' The convergence criterion is ', tol
      print *, ' '
    end if
    
    ! Sort and save eigenvalues
    s = real(d)
    s = bubble_sort(s, ind)
    do i=1,nev
      eivals(i) = d(ind(i))
    end do

    ! Now each process has a horizontal chunk of v (eigenvector matrix)
    ! We want them to have vertical chunks instead
    ! Each process broadcasts his chunk of eigenvectors
    ! Then all processes only save columns they are responsible for
    offset = 0 ! controls what part of eivecs is currently being filled out
    do bproc=0,nprocs-1 ! bproc = broadcasting processor id
      ! broadcast nloc
      if (bproc == iam) then
        bnloc = nloc
        call igebs2d(comm, 'A', ' ', 1, 1, bnloc, 1)
      else
        call igebr2d(comm, 'A', ' ', 1, 1, bnloc, 1, 0, bproc)
      end if
      ! Allocate space for chunk
      allocate(bv(bnloc,nev))
      ! broadcast horizontal chunk of v current process has
      if (bproc == iam) then
        bv = v(:,1:nev)
        call zgebs2d(comm, 'A', ' ', bnloc, nev, bv, bnloc)
      else
        call zgebr2d(comm, 'A', ' ', bnloc, nev, bv, bnloc, 0, bproc)
      end if

      ! Each process stores only some columns of eigenvector matrix (nevloc)
      ! Copies columns relevant for current process from the shared chunk of v into eivecs
      do myj=1,nevloc
        call l2g(myj,iam,nev,npcol,1,j) ! calculates global column number based on local column number and process number (columns are distributed cyclically) and stores result in j
        eivecs(offset+1:offset+bnloc,myj) = bv(:,ind(j)) ! saves relevant part of the shared chunk
      end do

      ! prepare to fill out next portion of eivecs
      deallocate(bv)
      offset = offset + bnloc 
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Find eigenvalues and eigenvectors of global variable rovib_ham (has to be built)
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine find_eigenpairs_parpack(context, n_states, ncv, max_iterations, eigenvalues, eigenvectors)
    integer, intent(in) :: context, n_states, ncv, max_iterations
    complex*16, allocatable, intent(out) :: eigenvalues(:), eigenvectors(:, :)
    integer :: proc_id, n_procs, proc_states
    external :: numroc
    integer :: numroc
  
    call get_proc_info(proc_id, n_procs)
    ! Computes number of states for this process. Workload is distributed in a round-robin fashion.
    proc_states = numroc(n_states, 1, proc_id, 0, n_procs)
    allocate(eigenvalues(n_states), eigenvectors(rovib_ham % global_chunk_info % columns, proc_states))
    call parz(context, eigenvalues, eigenvectors, rovib_ham % global_chunk_info % columns, size(rovib_ham % proc_chunk, 1), n_states, proc_states, ncv, max_iterations)
  end subroutine
end module
