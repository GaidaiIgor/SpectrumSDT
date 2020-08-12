!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test7f [-help] [-n <n>] [-verbose] [-inplace]
!
!  Description: Simple example that tests the matrix square root.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcfn.h>
      use slepcfn
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Mat            A,S,R
      FN             fn
      PetscInt       n
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg,verbose,inplace
      PetscReal      re,im,nrm
      PetscScalar    aa(1),tau,eta,alpha,x,y,yp
      PetscOffset    ia

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      n = 10
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,  &
     &                        '-n',n,flg,ierr)
      call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
     &                        '-verbose',verbose,ierr)
      call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
     &                        '-inplace',inplace,ierr)

      if (rank .eq. 0) then
        write(*,100) n
      endif
 100  format (/'Matrix square root, n =',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create FN object and matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call FNCreate(PETSC_COMM_WORLD,fn,ierr)
      call FNSetType(fn,FNSQRT,ierr)
      tau = 0.15
      eta = 1.0
      call FNSetScale(fn,tau,eta,ierr)
      call FNSetFromOptions(fn,ierr)
      call FNGetScale(fn,tau,eta,ierr)
      call FNView(fn,PETSC_NULL_VIEWER,ierr)

      call MatCreateSeqDense(PETSC_COMM_SELF,n,n,PETSC_NULL_SCALAR,A,   &
     &     ierr)
      call PetscObjectSetName(A,'A',ierr)
      call MatDenseGetArray(A,aa,ia,ierr)
      call FillUpMatrix(n,aa(ia+1))
      call MatDenseRestoreArray(A,aa,ia,ierr)
      call MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Scalar evaluation
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      x = 2.2
      call FNEvaluateFunction(fn,x,y,ierr)
      call FNEvaluateDerivative(fn,x,yp,ierr)

      if (rank .eq. 0) then
        re = PetscRealPart(y)
        im = PetscImaginaryPart(y)
        if (abs(im).lt.1.d-10) then
          write(*,110) 'f', PetscRealPart(x), re
        else
          write(*,120) 'f', PetscRealPart(x), re, im
        endif
        re = PetscRealPart(yp)
        im = PetscImaginaryPart(yp)
        if (abs(im).lt.1.d-10) then
          write(*,110) 'f''', PetscRealPart(x), re
        else
          write(*,120) 'f''', PetscRealPart(x), re, im
        endif
      endif
 110  format (A2,'(',F4.1,') = ',F8.5)
 120  format (A2,'(',F4.1,') = ',F8.5,SP,F8.5,'i')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute matrix square root
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call MatCreateSeqDense(PETSC_COMM_SELF,n,n,PETSC_NULL_SCALAR,S,   &
     &     ierr)
      call PetscObjectSetName(S,'S',ierr)
      if (inplace) then
        call MatCopy(A,S,SAME_NONZERO_PATTERN,ierr)
        call MatSetOption(S,MAT_HERMITIAN,PETSC_TRUE,ierr)
        call FNEvaluateFunctionMat(fn,S,PETSC_NULL_MAT,ierr)
      else
        call FNEvaluateFunctionMat(fn,A,S,ierr)
      endif
      if (verbose) then
        if (rank .eq. 0) write (*,*) 'Matrix A - - - - - - - -'
        call MatView(A,PETSC_NULL_VIEWER,ierr)
        if (rank .eq. 0) write (*,*) 'Computed sqrtm(A) - - - - - - - -'
        call MatView(S,PETSC_NULL_VIEWER,ierr)
      endif

!     *** check error ||S*S-A||_F
      call MatMatMult(S,S,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,R,ierr)
      if (eta .ne. 1.0) then
        alpha = 1.0/(eta*eta)
        call MatScale(R,alpha,ierr)
      endif
      alpha = -tau
      call MatAXPY(R,alpha,A,SAME_NONZERO_PATTERN,ierr)
      call MatNorm(R,NORM_FROBENIUS,nrm,ierr)
      if (nrm<100*PETSC_MACHINE_EPSILON) then
        write (*,*) '||S*S-A||_F < 100*eps'
      else
        write (*,130) nrm
      endif
 130  format ('||S*S-A||_F = ',F8.5)

!     *** Clean up
      call MatDestroy(S,ierr)
      call MatDestroy(R,ierr)
      call MatDestroy(A,ierr)
      call FNDestroy(fn,ierr)
      call SlepcFinalize(ierr)
      end

! -----------------------------------------------------------------

      subroutine FillUpMatrix(n,X)
      PetscInt    n,i,j
      PetscScalar X(n,n)

      do i=1,n
        X(i,i) = 2.5
      end do
      do j=1,2
        do i=1,n-j
          X(i,i+j) = 1.d0
          X(i+j,i) = 1.d0
        end do
      end do
      return
      end

!/*TEST
!
!   test:
!      suffix: 1
!      nsize: 1
!      args: -fn_scale .13,2 -n 19 -fn_method {{0 1 2 3}shared output}
!      filter: grep -v "computing matrix functions"
!      output_file: output/test7f_1.out
!
!TEST*/
