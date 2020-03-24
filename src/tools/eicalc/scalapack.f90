module scalapack
!-----------------------------------------------------------------------
!  This file contains drivers to ScaLapack
!  Calculates eigenpairs, adopted from sample code.
!  Author: Alexander Teplukhin
!-----------------------------------------------------------------------
  use index_conversion
  use algorithms
contains
  !-----------------------------------------------------------------------
  !  ScaLapack eigenproblem solver, real version.
  !  One-dimensional block column distribution.
  !  How matrix should be built in the calling program:
  !      do j=1,n
  !        call g2l(j,n,nprocs,1,id,myj)
  !        if (myid==id) then
  !          psi = 0d0
  !          psi(j) = 1d0
  !          call mult(psi,hpsi)
  !          A(:,myj) = hpsi
  !        endif
  !      enddo
  !  Number of columns in A should be calculated using
  !      lda = numroc( n, 1, iam, 0, nprocs)
  !-----------------------------------------------------------------------
  subroutine scald_cd(context,n,ns,eivecs,eivals,a,nc)
    implicit none
    integer              n,ns,nc
    integer              context,iam,info,lwork
    integer              npcol,nprow,mycol,myrow,nprocs
    integer              desc(50)
    real*8,allocatable:: w(:),work(:),z(:,:)
    real*8               eivals(ns),eivecs(n,ns),a(n,nc)
    real*8               lsize
    integer              j,myj,jproc,numroc
    external             blacs_pinfo,blacs_get,blacs_gridinfo
    external             descinit,pdsyev

    ! Get BLACS information
    call blacs_pinfo(iam,nprocs)
    call blacs_gridinfo(context,nprow,npcol,myrow,mycol)
    if(myrow==0.and.mycol==0)then
      print '(A,I0)','Matrix size is ',n
      print '(A,I0,A,I0)','Procs grid is ',nprow,'x',npcol
    endif

    ! Bail out if this process is not a part of this context.
    if(myrow==-1)return
    nprocs = nprow * npcol
    allocate(w(n),z(n,nc))

    ! Initialize descriptors
    call descinit( desc, n, n, 1, 1, 0, 0, context, n, info)

    ! Ask pdsyev to compute the entire eigendecomposition
    lwork = -1
    call pdsyev( 'V', 'U', n, a, 1, 1, desc, w, z, 1, 1, desc, lsize, lwork, info)
    lwork = INT(lsize)
    allocate(work(lwork))
    call pdsyev( 'V', 'U', n, a, 1, 1, desc, w, z, 1, 1, desc, work, lwork, info)

    ! Collect vectors
    do j=1,ns
      call g2l(j,n,npcol,1,jproc,myj)
      if(mycol==0)then
        if(jproc==0)then
          eivecs(:,j) = z(:,myj)
        else
          call dgerv2d(context,n,1,eivecs(1,j),n,0,jproc)
        endif
      elseif(mycol==jproc)then
        call dgesd2d(context,n,1,z(1,myj),n,0,0)
      endif
    enddo
    eivals = w(1:ns)
  end subroutine


  !-----------------------------------------------------------------------
  !  ScaLapack eigenproblem solver, real version.
  !  Two-dimensional block cyclic distribution.
  !  BLACS initialization and finish must be done in calling program.
  !-----------------------------------------------------------------------
  subroutine scald(context,n,ns,eivecs,eivals,a,nr,nc,nstloc,nb)
    implicit none
    integer              n,ns,nr,nc,nstloc
    integer              context,iam,info,lwork
    integer              npcol,nprow,mycol,myrow,nprocs
    integer              desc(50)
    real*8,allocatable:: w(:),work(:),z(:,:)
    real*8               eivals(ns),eivecs(n,nstloc),a(nr,nc)
    real*8               lsize
    integer              sproc,sproci,sprocj,myis
    integer              i,j,myi,myj,iproc,jproc,numroc,nb
    external             blacs_pinfo,blacs_get,blacs_gridinfo
    external             descinit,pdsyev

    ! Get BLACS information
    call blacs_pinfo(iam,nprocs)
    call blacs_gridinfo(context,nprow,npcol,myrow,mycol)
    if(myrow==0.and.mycol==0)then
      print '(A,I0)','Matrix size is ',n
      print '(A,I0)','Block size is ', nb
      print '(A,I0,A,I0)','Procs grid is ',nprow,'x',npcol
    endif

    ! Bail out if this process is not a part of this context.
    if(myrow==-1)return
    nprocs = nprow * npcol
    allocate(w(n),z(nr,nc))

    ! Initialize descriptors
    call descinit( desc, n, n, nb, nb, 0, 0, context, nr, info)

    ! Ask pdsyev to compute the entire eigendecomposition
    lwork = -1
    call pdsyev( 'V', 'U', n, a, 1, 1, desc, w, z, 1, 1, desc, lsize, lwork, info)
    lwork = INT(lsize)
    allocate(work(lwork))
    call pdsyev( 'V', 'U', n, a, 1, 1, desc, w, z, 1, 1, desc, work, lwork, info)

    ! Collect vectors
    eivecs = 0
    do myj=1,nstloc
      do myi=1,nr
        call l2g(myi,myrow,n,nprow,nb,i)
        eivecs(i,myj) = z(myi,myj)
      enddo
    enddo
    call dgsum2d(context,'C',' ',n,nstloc,eivecs,n,-1,-1)
    eivals = w(1:ns)
  end subroutine

  !-----------------------------------------------------------------------
  !  ScaLapack eigenproblem solver, complex hermitian version.
  !  Two-dimensional block cyclic distribution.
  !  BLACS initialization and finish must be done in calling program.
  !-----------------------------------------------------------------------
  subroutine scalh(context,n,ns,eivecs,eivals,a,nr,nc,nstloc,nb)
    implicit none
    integer                 n,ns,nr,nc,nstloc
    integer                 context,iam,info,lwork,lrwork
    integer                 npcol,nprow,mycol,myrow,nprocs
    integer                 desc(50)
    real*8,allocatable::    w(:)
    complex*16,allocatable::work(:),rwork(:),z(:,:)
    complex*16              eivals(ns),eivecs(n,nstloc),a(nr,nc)
    complex*16              work1,rwork1
    integer                 i,j,myi,myj,jproc,numroc,nb

    ! Get BLACS information
    call blacs_pinfo(iam,nprocs)
    call blacs_gridinfo(context,nprow,npcol,myrow,mycol)
    if(myrow==0.and.mycol==0)then
      print '(A,I0)','Matrix size is ',n
      print '(A,I0)','Block size is ', nb
      print '(A,I0,A,I0)','Procs grid is ',nprow,'x',npcol
    endif

    ! Bail out if this process is not a part of this context.
    if(myrow==-1)return
    nprocs = nprow * npcol
    allocate(w(n),z(nr,nc))

    ! Initialize descriptors
    call descinit( desc, n, n, nb, nb, 0, 0, context, nr, info)

    ! Ask pdsyev to compute the entire eigendecomposition
    lwork  = -1
    lrwork = -1
    call pzheev( 'V', 'U', n, a, 1, 1, desc, w, z, 1, 1, desc, work1, lwork, rwork1, lrwork, info)
    lwork  = INT(work1)
    lrwork = INT(rwork1)
    allocate(work(lwork))
    allocate(rwork(lrwork))
    call pzheev( 'V', 'U', n, a, 1, 1, desc, w, z, 1, 1, desc, work, lwork, rwork, lrwork, info)

    ! Collect vectors
    eivecs = 0
    do myj=1,nstloc
      do myi=1,nr
        call l2g(myi,myrow,n,nprow,nb,i)
        eivecs(i,myj) = z(myi,myj)
      enddo
    enddo
    call zgsum2d(context,'C',' ',n,nstloc,eivecs,n,-1,-1)
    eivals = w(1:ns)
  end subroutine

  !-----------------------------------------------------------------------
  !  ScaLapack eigenproblem solver, complex nonhermitian version.
  !  Two-dimensional block cyclic distribution.
  !  BLACS initialization and finish must be done in calling program.
  !-----------------------------------------------------------------------
  subroutine scalz(context,n,ns,eivecs,eivals,a,nr,nc,nstloc,nb)
    implicit none
    integer                 n,ns,nr,nc,nstloc
    integer                 context,iam,info,lwork,lrwork,ltau
    integer                 npcol,nprow,mycol,myrow,nprocs
    integer                 desca(50),descv(50)
    complex*16,allocatable::q(:,:),z(:,:)
    complex*16,allocatable::w(:),tau(:),work(:)
    real*8, allocatable::   rwork(:)
    complex*16              eivals(ns),eivecs(n,nstloc),a(nr,nc)
    complex*16              work1
    real*8                  s(n)
    integer, allocatable :: ind(:)
    integer                 i,j,myi,myj,numroc,nb
    integer                 procs,procd,myjs,myjd

    ! One-dimensional cyclic column distribution for eigenvectors
    complex*16,allocatable::v(:,:)
    integer                 ncl
    integer                 lcontext

    ! Get BLACS information
    ! Block size must be square and at least six
    call blacs_pinfo(iam,nprocs)
    call blacs_gridinfo(context,nprow,npcol,myrow,mycol)
    if(myrow==0.and.mycol==0)then
      print '(A,I0)','Matrix size is ',n
      print '(A,I0)','Block size is ', nb
      print '(A,I0,A,I0)','Procs grid is ',nprow,'x',npcol
    endif

    ! Bail out if this process is not a part of this context.
    if(myrow==-1)return
    nprocs = nprow * npcol
    ltau = numroc(n-1, nb, mycol, 0, npcol)
    allocate(w(n),q(nr,nc),z(nr,nc),tau(ltau),rwork(n))

    ! Initialize descriptors
    call descinit( desca, n, n, nb, nb, 0, 0, context, nr, info)

    ! Hessenberg reduction
    lwork = -1
    call pzgehrd( n, 1, n, a, 1, 1, desca, tau, work1, lwork, info )
    lwork = int(work1)
    allocate(work(lwork))
    call pzgehrd( n, 1, n, a, 1, 1, desca, tau, work, lwork, info )
    Q = A

    ! Add zeros to hessenberg matrix and set Z to identity matrix
    z = 0
    do myj=1,nc
      call l2g(myj,mycol,n,npcol,nb,j)
      do myi=1,nr
        call l2g(myi,myrow,n,nprow,nb,i)
        if(i>j+1)a(myi,myj) = 0
        if(i==j)z(myi,myj) = 1
      enddo
    enddo

    ! Schur decomposition
    lwork = -1
    call pzlahqr(.true.,.true., n, 1, n, a, desca, w, 1, n, z, desca, work1, lwork, 0, 0, info )
    lwork = int(work1)
    deallocate(work)
    allocate(work(lwork))
    call pzlahqr(.true.,.true., n, 1, n, a, desca, w, 1, n, z, desca, work, lwork, 0, 0, info )

    ! Extract eigenvectors
    deallocate(work)
    allocate(work(2*nr))
    call pztrevc('R', 'B', 0, n, a, desca, 0, 0, z, desca, n, n, work, rwork, info )
    lwork = -1
    call pzunmhr('l', 'n', n, n, 1, n, q, 1, 1, desca, tau, z, 1, 1, desca, work1, lwork, info )
    lwork = int(work1)
    deallocate(work)
    allocate(work(lwork))
    call pzunmhr('l', 'n', n, n, 1, n, q, 1, 1, desca, tau, z, 1, 1, desca, work, lwork, info )

    ! Setup linear process grid
    call blacs_get(0,0,lcontext)
    call blacs_gridinit(lcontext,'R',1,nprocs)

    ! Initialize descriptor
    call descinit(descv, n, n, n, 1, 0, 0, lcontext, n, info)

    ! Allocate storage for new distribution
    ncl = numroc(n, 1, iam, 0, nprocs)
    allocate(v(n,ncl))

    ! Convert to one-dimensional cyclic column distribution
    call pzgemr2d(n, n, z, 1, 1, desca, v, 1, 1, descv, context)

    ! Sorting
    s = real(w)
    s = bubble_sort(s, ind)
    eivecs = 0
    do j=1,ns
      eivals(j) = w(ind(j))
      call g2l(j     ,ns,nprocs,1,procd,myjd)
      call g2l(ind(j),ns,nprocs,1,procs,myjs)
      if(iam == procd)then
        if(iam == procs)then
          eivecs(:,myjd) = v(:,myjs)
        else
          call zgerv2d(lcontext,n,1,eivecs(1,myjd),n,0,procs)
        endif
      elseif(iam == procs)then
        call zgesd2d(lcontext,n,1,v(1,myjs),n,0,procd)
      endif
    enddo

    ! Release liner process grid
    call blacs_gridexit(lcontext)
  end subroutine
end module