!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex16f90 [-help] [-n <n>] [-m <m>] [SLEPc opts]
!
!  Description: Simple example that solves a quadratic eigensystem with the
!  PEP object. This is the Fortran90 equivalent to ex16.c
!
!  The command line options are:
!    -n <n>, where <n> = number of grid subdivisions in x dimension
!    -m <m>, where <m> = number of grid subdivisions in y dimension
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcpep.h>
      use slepcpep
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     M,C,K  problem matrices
!     pep    polynomial eigenproblem solver context

      Mat            M, C, K, A(3)
      PEP            pep
      PEPType        tname
      PetscInt       N, nx, ny, i, j, Istart, Iend, II
      PetscInt       nev, ithree
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg, terse
      PetscScalar    mone, two, four, val

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'SlepcInitialize failed'
        stop
      endif
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRA(ierr)
      nx = 10
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',nx,flg,ierr);CHKERRA(ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-m',ny,flg,ierr);CHKERRA(ierr)
      if (.not. flg) then
        ny = nx
      endif
      N = nx*ny
      if (rank .eq. 0) then
        write(*,100) N, nx, ny
      endif
 100  format (/'Quadratic Eigenproblem, N=',I6,' (',I4,'x',I4,' grid)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** K is the 2-D Laplacian
      call MatCreate(PETSC_COMM_WORLD,K,ierr);CHKERRA(ierr)
      call MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr);CHKERRA(ierr)
      call MatSetFromOptions(K,ierr);CHKERRA(ierr)
      call MatSetUp(K,ierr);CHKERRA(ierr)
      call MatGetOwnershipRange(K,Istart,Iend,ierr);CHKERRA(ierr)
      mone = -1.0
      four = 4.0
      do II=Istart,Iend-1
        i = II/nx
        j = II-i*nx
        if (i .gt. 0) then
          call MatSetValue(K,II,II-nx,mone,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        if (i .lt. ny-1) then
          call MatSetValue(K,II,II+nx,mone,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        if (j .gt. 0) then
          call MatSetValue(K,II,II-1,mone,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        if (j .lt. nx-1) then
          call MatSetValue(K,II,II+1,mone,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        call MatSetValue(K,II,II,four,INSERT_VALUES,ierr);CHKERRA(ierr)
      end do
      call MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

!     ** C is the 1-D Laplacian on horizontal lines
      call MatCreate(PETSC_COMM_WORLD,C,ierr);CHKERRA(ierr)
      call MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr);CHKERRA(ierr)
      call MatSetFromOptions(C,ierr);CHKERRA(ierr)
      call MatSetUp(C,ierr);CHKERRA(ierr)
      call MatGetOwnershipRange(C,Istart,Iend,ierr);CHKERRA(ierr)
      two = 2.0
      do II=Istart,Iend-1
        i = II/nx
        j = II-i*nx
        if (j .gt. 0) then
          call MatSetValue(C,II,II-1,mone,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        if (j .lt. nx-1) then
          call MatSetValue(C,II,II+1,mone,INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        call MatSetValue(C,II,II,two,INSERT_VALUES,ierr);CHKERRA(ierr)
      end do
      call MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

!     ** M is a diagonal matrix
      call MatCreate(PETSC_COMM_WORLD,M,ierr);CHKERRA(ierr)
      call MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr);CHKERRA(ierr)
      call MatSetFromOptions(M,ierr);CHKERRA(ierr)
      call MatSetUp(M,ierr);CHKERRA(ierr)
      call MatGetOwnershipRange(M,Istart,Iend,ierr);CHKERRA(ierr)
      do II=Istart,Iend-1
        val = II+1
        call MatSetValue(M,II,II,val,INSERT_VALUES,ierr);CHKERRA(ierr)
      end do
      call MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      call PEPCreate(PETSC_COMM_WORLD,pep,ierr);CHKERRA(ierr)

!     ** Set matrices and problem type
      A(1) = K
      A(2) = C
      A(3) = M
      ithree = 3
      call PEPSetOperators(pep,ithree,A,ierr);CHKERRA(ierr)
      call PEPSetProblemType(pep,PEP_GENERAL,ierr);CHKERRA(ierr)

!     ** Set solver parameters at runtime
      call PEPSetFromOptions(pep,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call PEPSolve(pep,ierr);CHKERRA(ierr)

!     ** Optional: Get some information from the solver and display it
      call PEPGetType(pep,tname,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,120) tname
      endif
 120  format (' Solution method: ',A)
      call PEPGetDimensions(pep,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,130) nev
      endif
 130  format (' Number of requested eigenvalues:',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** show detailed info unless -terse option is given by user
      call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-terse',terse,ierr);CHKERRA(ierr)
      if (terse) then
        call PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_NULL_VIEWER,ierr);CHKERRA(ierr)
      else
        call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr);CHKERRA(ierr)
        call PEPReasonView(pep,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
        call PEPErrorView(pep,PEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
        call PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
      endif
      call PEPDestroy(pep,ierr);CHKERRA(ierr)
      call MatDestroy(K,ierr);CHKERRA(ierr)
      call MatDestroy(C,ierr);CHKERRA(ierr)
      call MatDestroy(M,ierr);CHKERRA(ierr)
      call SlepcFinalize(ierr)
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -pep_nev 4 -pep_ncv 19 -terse
!      requires: !complex
!
!TEST*/
