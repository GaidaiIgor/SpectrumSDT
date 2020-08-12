!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test14f [-help] [-n <n>] [all SLEPc options]
!
!  Description: Simple example that tests solving a DSNHEP problem.
!
!  The command line options are:
!    -n <n>, where <n> = matrix size
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcds.h>
      use slepcds
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     A     problem matrix
!     ds    dense solver context

      Mat            A
      DS             ds
      PetscInt       n, i, ld, zero
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg
      PetscScalar    aa(1), wr(100), wi(100)
      PetscReal      re, im
      PetscOffset    ia

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      zero = 0
      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'SlepcInitialize failed'
        stop
      endif
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRA(ierr)
      n = 10
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr);CHKERRA(ierr)
      if (n .gt. 100) then; SETERRA(PETSC_COMM_SELF,1,'Program currently limited to n=100'); endif

      if (rank .eq. 0) then
        write(*,110) n
      endif
 110  format (/'Solve a Dense System of type NHEP, n =',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create DS object
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call DSCreate(PETSC_COMM_WORLD,ds,ierr);CHKERRA(ierr)
      call DSSetType(ds,DSNHEP,ierr);CHKERRA(ierr)
      call DSSetFromOptions(ds,ierr);CHKERRA(ierr)
      ld = n
      call DSAllocate(ds,ld,ierr);CHKERRA(ierr)
      call DSSetDimensions(ds,n,zero,zero,zero,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Fill with Grcar matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call DSGetMat(ds,DS_MAT_A,A,ierr);CHKERRA(ierr)
      call MatDenseGetArray(A,aa,ia,ierr);CHKERRA(ierr)
      call FillUpMatrix(n,aa(ia+1))
      call MatDenseRestoreArray(A,aa,ia,ierr);CHKERRA(ierr)
      call DSRestoreMat(ds,DS_MAT_A,A,ierr);CHKERRA(ierr)
      call DSSetState(ds,DS_STATE_INTERMEDIATE,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the problem and show eigenvalues
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call DSSolve(ds,wr,wi,ierr);CHKERRA(ierr)
!     call DSSort(ds,wr,wi,PETSC_NULL_SCALAR,PETSC_NULL_SCALAR,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)

      if (rank .eq. 0) then
        write(*,*) 'Computed eigenvalues ='
        do i=1,n
#if defined(PETSC_USE_COMPLEX)
          re = PetscRealPart(wr(i))
          im = PetscImaginaryPart(wr(i))
#else
          re = wr(i)
          im = wi(i)
#endif
          if (abs(im).lt.1.d-10) then
            write(*,120) re
          else
            write(*,130) re, im
          endif
        end do
      endif
 120  format ('  ',F8.5)
 130  format ('  ',F8.5,SP,F8.5,'i')

!     *** Clean up
      call DSDestroy(ds,ierr);CHKERRA(ierr)
      call SlepcFinalize(ierr)
      end

! -----------------------------------------------------------------

      subroutine FillUpMatrix(n,X)
      PetscInt    n,i,j
      PetscScalar X(n,n)

      do i=2,n
        X(i,i-1) = -1.d0
      end do
      do j=0,3
        do i=1,n-j
          X(i,i+j) = 1.d0
        end do
      end do
      return
      end

!/*TEST
!
!   test:
!      suffix: 1
!      requires: !complex
!
!TEST*/
