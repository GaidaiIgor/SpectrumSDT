!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Description: Test Fortran interface of spectrum-slicing Krylov-Schur.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepceps.h>
      use slepceps
      implicit none

#define MAXSUB 16
#define MAXSHI 16

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Mat            A,B,As,Bs,Au
      EPS            eps
      ST             st
      KSP            ksp
      PC             pc
      Vec            v
      PetscScalar    value
      PetscInt       n,m,i,j,k,Istart,Iend
      PetscInt       nev,ncv,mpd,nval
      PetscInt       row,col,nloc,nlocs,mlocs
      PetscInt       II,npart,inertias(MAXSHI)
      PetscBool      flg,lock
      PetscMPIInt    size,rank
      PetscReal      int0,int1,keep,subint(MAXSUB)
      PetscReal      shifts(MAXSHI)
      PetscScalar    eval,one,mone,zero
      PetscErrorCode ierr
      MPI_Comm       comm

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'SlepcInitialize failed'
        stop
      endif
      call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRA(ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRA(ierr)
      n = 35
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr);CHKERRA(ierr)
      m = n*n
      if (rank .eq. 0) then
        write(*,100) n
      endif
 100  format (/'Spectrum-slicing test, n =',I3,' (Fortran)'/)

      call MatCreate(PETSC_COMM_WORLD,A,ierr);CHKERRA(ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m,ierr);CHKERRA(ierr)
      call MatSetFromOptions(A,ierr);CHKERRA(ierr)
      call MatSetUp(A,ierr);CHKERRA(ierr)
      call MatCreate(PETSC_COMM_WORLD,B,ierr);CHKERRA(ierr)
      call MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,m,m,ierr);CHKERRA(ierr)
      call MatSetFromOptions(B,ierr);CHKERRA(ierr)
      call MatSetUp(B,ierr);CHKERRA(ierr)
      call MatGetOwnershipRange(A,Istart,Iend,ierr);CHKERRA(ierr)
      do II=Istart,Iend-1
        i = II/n
        j = II-i*n
        value = -1.0
        row = II
        if (i>0) then
          col = II-n
          call MatSetValue(A,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        if (i<n-1) then
          col = II+n
          call MatSetValue(A,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        if (j>0) then
          col = II-1
          call MatSetValue(A,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        if (j<n-1) then
          col = II+1
          call MatSetValue(A,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        col = II
        value = 4.0
        call MatSetValue(A,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        value = 2.0
        call MatSetValue(B,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
      enddo
      if (Istart .eq. 0) then
        row = 0
        col = 0
        value = 6.0
        call MatSetValue(B,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        row = 0
        col = 1
        value = -1.0
        call MatSetValue(B,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        row = 1
        col = 0
        value = -1.0
        call MatSetValue(B,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        row = 1
        col = 1
        value = 1.0
        call MatSetValue(B,row,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
      endif
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call EPSCreate(PETSC_COMM_WORLD,eps,ierr);CHKERRA(ierr)
      call EPSSetOperators(eps,A,B,ierr);CHKERRA(ierr)
      call EPSSetProblemType(eps,EPS_GHEP,ierr);CHKERRA(ierr)
      call EPSSetType(eps,EPSKRYLOVSCHUR,ierr);CHKERRA(ierr)

!     Set interval and other settings for spectrum slicing

      call EPSSetWhichEigenpairs(eps,EPS_ALL,ierr);CHKERRA(ierr)
      int0 = 1.1
      int1 = 1.3
      call EPSSetInterval(eps,int0,int1,ierr);CHKERRA(ierr)
      call EPSGetST(eps,st,ierr);CHKERRA(ierr)
      call STSetType(st,STSINVERT,ierr);CHKERRA(ierr)
      call STGetKSP(st,ksp,ierr);CHKERRA(ierr)
      call KSPGetPC(ksp,pc,ierr);CHKERRA(ierr)
      call KSPSetType(ksp,KSPPREONLY,ierr);CHKERRA(ierr)
      call PCSetType(pc,PCCHOLESKY,ierr);CHKERRA(ierr)

!     Test interface functions of Krylov-Schur solver

      call EPSKrylovSchurGetRestart(eps,keep,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,110) keep
      endif
 110  format (' Restart parameter before changing = ',f7.4)
      keep = 0.4
      call EPSKrylovSchurSetRestart(eps,keep,ierr);CHKERRA(ierr)
      call EPSKrylovSchurGetRestart(eps,keep,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,120) keep
      endif
 120  format (' ... changed to ',f7.4)

      call EPSKrylovSchurGetLocking(eps,lock,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,130) lock
      endif
 130  format (' Locking flag before changing = ',L)
      call EPSKrylovSchurSetLocking(eps,PETSC_FALSE,ierr);CHKERRA(ierr)
      call EPSKrylovSchurGetLocking(eps,lock,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,140) lock
      endif
 140  format (' ... changed to ',L)

      call EPSKrylovSchurGetDimensions(eps,nev,ncv,mpd,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,150) nev,ncv,mpd
      endif
 150  format (' Sub-solve dimensions before changing: nev=',I2,', ncv=',I2,', mpd=',I2)
      nev = 30
      ncv = 60
      mpd = 60
      call EPSKrylovSchurSetDimensions(eps,nev,ncv,mpd,ierr);CHKERRA(ierr)
      call EPSKrylovSchurGetDimensions(eps,nev,ncv,mpd,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,160) nev,ncv,mpd
      endif
 160  format (' ... changed to: nev=',I2,', ncv=',I2,', mpd=',I2)

      if (size>0) then
        npart = size
        call EPSKrylovSchurSetPartitions(eps,npart,ierr);CHKERRA(ierr)
        call EPSKrylovSchurGetPartitions(eps,npart,ierr);CHKERRA(ierr)
        if (rank .eq. 0) then
          write(*,170) npart
        endif
 170    format (' Using ',I2,' partitions')
        if (npart>MAXSUB) then; SETERRA(PETSC_COMM_SELF,1,'Too many subintervals'); endif

        subint(1) = int0
        subint(npart+1) = int1
        do i=2,npart
          subint(i) = int0+(i-1)*(int1-int0)/npart
        enddo
        call EPSKrylovSchurSetSubintervals(eps,subint,ierr);CHKERRA(ierr)
        call EPSKrylovSchurGetSubintervals(eps,subint,ierr);CHKERRA(ierr)
        if (rank .eq. 0) then
          write(*,*) 'Using sub-interval separations ='
          do i=2,npart
            write(*,180) subint(i)
          enddo
        endif
 180    format (f7.4)
      endif

      call EPSSetFromOptions(eps,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute all eigenvalues in interval and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call EPSSetUp(eps,ierr);CHKERRA(ierr)
      call EPSKrylovSchurGetInertias(eps,k,PETSC_NULL_REAL,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      if (k>MAXSHI) then; SETERRA(PETSC_COMM_SELF,1,'Too many shifts'); endif
      call EPSKrylovSchurGetInertias(eps,k,shifts,inertias,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,*) 'Inertias after EPSSetUp:'
        do i=1,k
          write(*,185) shifts(i),inertias(i)
        enddo
      endif
 185  format (' .. ',f4.1,' (',I3,')')

      call EPSSolve(eps,ierr);CHKERRA(ierr)
      call EPSGetDimensions(eps,nev,ncv,mpd,ierr);CHKERRA(ierr)
      call EPSGetInterval(eps,int0,int1,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,190) nev,int0,int1
      endif
 190  format (' Found ',I2,' eigenvalues in interval [',f7.4,',',f7.4,']')

      if (size>0) then
        call EPSKrylovSchurGetSubcommInfo(eps,k,nval,v,ierr);CHKERRA(ierr)
        if (rank .eq. 0) then
          write(*,200) rank,k,nval
          do i=0,nval-1
            call EPSKrylovSchurGetSubcommPairs(eps,i,eval,v,ierr);CHKERRA(ierr)
            write(*,210) PetscRealPart(eval)
          enddo
        endif
 200    format (' Process ',I2,' has worked in sub-interval ',I2,', containing ',I2,' eigenvalues')
 210    format (f7.4)
        call VecDestroy(v,ierr);CHKERRA(ierr)

        call EPSKrylovSchurGetSubcommMats(eps,As,Bs,ierr);CHKERRA(ierr)
        call MatGetLocalSize(A,nloc,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
        call MatGetLocalSize(As,nlocs,mlocs,ierr);CHKERRA(ierr)
        if (rank .eq. 0) then
          write(*,220) rank,nloc,nlocs
        endif
 220    format (' Process ',I2,' owns ',I5,', rows of the global',' matrices, and ',I5,' rows in the subcommunicator')


!       modify A on subcommunicators
        call PetscObjectGetComm(As,comm,ierr);CHKERRA(ierr)
        call MatCreate(comm,Au,ierr);CHKERRA(ierr)
        call MatSetSizes(Au,nlocs,mlocs,m,m,ierr);CHKERRA(ierr)
        call MatSetFromOptions(Au,ierr);CHKERRA(ierr)
        call MatSetUp(Au,ierr);CHKERRA(ierr)
        call MatGetOwnershipRange(Au,Istart,Iend,ierr);CHKERRA(ierr)
        do II=Istart,Iend-1
          value = 0.5
          call MatSetValue(Au,II,II,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        end do
        call MatAssemblyBegin(Au,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
        call MatAssemblyEnd(Au,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
        one = 1.0
        mone = -1.0
        zero = 0.0
        call EPSKrylovSchurUpdateSubcommMats(eps,one,mone,Au,zero,zero, &
     & PETSC_NULL_MAT,DIFFERENT_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRA(ierr)
        call MatDestroy(Au,ierr);CHKERRA(ierr)
      endif

      call EPSDestroy(eps,ierr);CHKERRA(ierr)
      call MatDestroy(A,ierr);CHKERRA(ierr)
      call MatDestroy(B,ierr);CHKERRA(ierr)

      call SlepcFinalize(ierr)
      end

!/*TEST
!
!   test:
!      suffix: 1
!      nsize: 2
!
!TEST*/
