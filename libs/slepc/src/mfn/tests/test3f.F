!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test3f [-help] [-n <n>] [all SLEPc options]
!
!  Description: square root of the 2-D Laplacian.
!
!  The command line options are:
!    -n <n>, where <n> = matrix rows and columns
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcmfn.h>
      use slepcmfn
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      Mat                A, B
      MFN                mfn
      FN                 f
      MFNConvergedReason reason;
      Vec                v, y
      PetscInt           Nt, n, i, j, II
      PetscInt           Istart, maxit, ncv
      PetscInt           col, its, Iend
      PetscScalar        val
      PetscReal          tol, norm
      PetscMPIInt        rank
      PetscErrorCode     ierr
      PetscBool          flg

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      n = 4
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,  &
     &                        '-n',n,flg,ierr)
      Nt = n*n

      if (rank .eq. 0) then
        write(*,100) n
      endif
 100  format (/'nSquare root of Laplacian, n=',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute the discrete 2-D Laplacian
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call MatCreate(PETSC_COMM_WORLD,A,ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Nt,Nt,ierr)
      call MatSetFromOptions(A,ierr)
      call MatSetUp(A,ierr)

      call MatGetOwnershipRange(A,Istart,Iend,ierr)
      do II=Istart,Iend-1
        i = II/n
        j = II-i*n
        val = -1.0
        if (i .gt. 0) then
          col = II-n
          call MatSetValue(A,II,col,val,INSERT_VALUES,ierr)
        end if
        if (i .lt. n-1) then
          col = II+n
          call MatSetValue(A,II,col,val,INSERT_VALUES,ierr)
        end if
        if (j .gt. 0) then
          col = II-1
          call MatSetValue(A,II,col,val,INSERT_VALUES,ierr)
        end if
        if (j .lt. n-1) then
          col = II+1
          call MatSetValue(A,II,col,val,INSERT_VALUES,ierr)
        end if
        val = 4.0
        call MatSetValue(A,II,II,val,INSERT_VALUES,ierr)
      enddo

      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

      call MatCreateVecs(A,PETSC_NULL_VEC,v,ierr)
      i = 0
      val = 1.0
      call VecSetValue(v,i,val,INSERT_VALUES,ierr)
      call VecAssemblyBegin(v,ierr)
      call VecAssemblyEnd(v,ierr)
      call VecDuplicate(v,y,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute singular values
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call MFNCreate(PETSC_COMM_WORLD,mfn,ierr)
      call MFNSetOperator(mfn,A,ierr)
      call MFNGetFN(mfn,f,ierr)
      call FNSetType(f,FNSQRT,ierr)

!     ** test some interface functions
      call MFNGetOperator(mfn,B,ierr)
      call MatView(B,PETSC_VIEWER_STDOUT_WORLD,ierr)
      call MFNSetOptionsPrefix(mfn,'myprefix_',ierr)
      tol = 1e-4
      maxit = 500
      call MFNSetTolerances(mfn,tol,maxit,ierr)
      ncv = 6
      call MFNSetDimensions(mfn,ncv,ierr)
      call MFNSetErrorIfNotConverged(mfn,PETSC_TRUE,ierr)
      call MFNSetFromOptions(mfn,ierr)

!     ** query properties and print them
      call MFNGetTolerances(mfn,tol,maxit,ierr)
      if (rank .eq. 0) then
        write(*,110) tol,maxit
      endif
 110  format (/' Tolerance: ',F7.4,', maxit: ',I4)
      call MFNGetDimensions(mfn,ncv,ierr)
      if (rank .eq. 0) then
        write(*,120) ncv
      endif
 120  format (' Subspace dimension: ',I3)
      call MFNGetErrorIfNotConverged(mfn,flg,ierr)
      if (rank .eq. 0 .and. flg) then
        write(*,*) 'Erroring out if convergence fails'
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Call the solver
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call MFNSolve(mfn,v,y,ierr)
      call MFNGetConvergedReason(mfn,reason,ierr)
      if (rank .eq. 0) then
        write(*,130) reason
      endif
 130  format (' Converged reason:',I2)
      call MFNGetIterationNumber(mfn,its,ierr)
!     if (rank .eq. 0) then
!       write(*,140) its
!     endif
!140  format (' Number of iterations of the method:',I4)

      call VecNorm(y,NORM_2,norm,ierr)
      if (rank .eq. 0) then
        write(*,150) norm
      endif
 150  format (' sqrt(A)*v has norm ',F7.4)

      call MFNDestroy(mfn,ierr)
      call MatDestroy(A,ierr)
      call VecDestroy(v,ierr)
      call VecDestroy(y,ierr)

      call SlepcFinalize(ierr)
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -log_exclude mfn
!
!TEST*/
