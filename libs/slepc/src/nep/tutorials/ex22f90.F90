!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex22f90 [-n <n>] [-tau <tau>] [SLEPc opts]
!
!  Description: Delay differential equation. Fortran90 equivalent of ex22.c
!
!  The command line options are:
!    -n <n>, where <n> = number of grid subdivisions
!    -tau <tau>, where <tau> = delay parameter
!
! ----------------------------------------------------------------------
!  Solve parabolic partial differential equation with time delay tau
!
!           u_t = u_xx + aa*u(t) + bb*u(t-tau)
!           u(0,t) = u(pi,t) = 0
!
!  with aa = 20 and bb(x) = -4.1+x*(1-exp(x-pi)).
!
!  Discretization leads to a DDE of dimension n
!
!           -u' = A*u(t) + B*u(t-tau)
!
!  which results in the nonlinear eigenproblem
!
!           (-lambda*I + A + exp(-tau*lambda)*B)*u = 0
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcnep.h>
      use slepcnep
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     nep       nonlinear eigensolver context
!     Id,A,B    problem matrices
!     f1,f2,f3  functions to define the nonlinear operator

      Mat            Id, A, B, mats(3)
      FN             f1, f2, f3, funs(3)
      NEP            nep
      NEPType        tname
      PetscScalar    one, bb, coeffs(2), scal
      PetscReal      tau, h, aa, xi, tol
      PetscInt       n, i, k, nev, Istart, Iend
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg, terse

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'SlepcInitialize failed'
        stop
      endif
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRA(ierr)
      n = 128
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr);CHKERRA(ierr)
      tau = 0.001
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-tau',tau,flg,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,100) n, tau
      endif
 100  format (/'Delay Eigenproblem, n =',I4,', tau =',F6.3)

      one = 1.0
      aa = 20.0
      h = PETSC_PI/real(n+1)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create problem matrices
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Id is the identity matrix
      call MatCreate(PETSC_COMM_WORLD,Id,ierr);CHKERRA(ierr)
      call MatSetSizes(Id,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
      call MatSetFromOptions(Id,ierr);CHKERRA(ierr)
      call MatSetUp(Id,ierr);CHKERRA(ierr)
      call MatGetOwnershipRange(Id,Istart,Iend,ierr);CHKERRA(ierr)
      do i=Istart,Iend-1
        call MatSetValue(Id,i,i,one,INSERT_VALUES,ierr);CHKERRA(ierr)
      end do
      call MatAssemblyBegin(Id,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(Id,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatSetOption(Id,MAT_HERMITIAN,PETSC_TRUE,ierr);CHKERRA(ierr)

!     ** A = 1/h^2*tridiag(1,-2,1) + aa*I
      call MatCreate(PETSC_COMM_WORLD,A,ierr);CHKERRA(ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
      call MatSetFromOptions(A,ierr);CHKERRA(ierr)
      call MatSetUp(A,ierr);CHKERRA(ierr)
      call MatGetOwnershipRange(A,Istart,Iend,ierr);CHKERRA(ierr)
      coeffs(1) = 1.0/(h*h)
      coeffs(2) = -2.0/(h*h)+aa
      do i=Istart,Iend-1
        if (i .gt. 0) then
          call MatSetValue(A,i,i-1,coeffs(1),INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        if (i .lt. n-1) then
          call MatSetValue(A,i,i+1,coeffs(1),INSERT_VALUES,ierr);CHKERRA(ierr)
        endif
        call MatSetValue(A,i,i,coeffs(2),INSERT_VALUES,ierr);CHKERRA(ierr)
      end do
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatSetOption(A,MAT_HERMITIAN,PETSC_TRUE,ierr);CHKERRA(ierr)

!     ** B = diag(bb(xi))
      call MatCreate(PETSC_COMM_WORLD,B,ierr);CHKERRA(ierr)
      call MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
      call MatSetFromOptions(B,ierr);CHKERRA(ierr)
      call MatSetUp(B,ierr);CHKERRA(ierr)
      call MatGetOwnershipRange(B,Istart,Iend,ierr);CHKERRA(ierr)
      do i=Istart,Iend-1
        xi = (i+1)*h
        bb  = -4.1+xi*(1.0-exp(xi-PETSC_PI))
        call MatSetValue(B,i,i,bb,INSERT_VALUES,ierr);CHKERRA(ierr)
      end do
      call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatSetOption(B,MAT_HERMITIAN,PETSC_TRUE,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create problem functions, f1=-lambda, f2=1.0, f3=exp(-tau*lambda)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call FNCreate(PETSC_COMM_WORLD,f1,ierr);CHKERRA(ierr)
      call FNSetType(f1,FNRATIONAL,ierr);CHKERRA(ierr)
      k = 2
      coeffs(1) = -1.0
      coeffs(2) = 0.0
      call FNRationalSetNumerator(f1,k,coeffs,ierr);CHKERRA(ierr)

      call FNCreate(PETSC_COMM_WORLD,f2,ierr);CHKERRA(ierr)
      call FNSetType(f2,FNRATIONAL,ierr);CHKERRA(ierr)
      k = 1
      coeffs(1) = 1.0
      call FNRationalSetNumerator(f2,k,coeffs,ierr);CHKERRA(ierr)

      call FNCreate(PETSC_COMM_WORLD,f3,ierr);CHKERRA(ierr)
      call FNSetType(f3,FNEXP,ierr);CHKERRA(ierr)
      scal = -tau
      call FNSetScale(f3,scal,one,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      call NEPCreate(PETSC_COMM_WORLD,nep,ierr);CHKERRA(ierr)

!     ** Set the split operator. Note that A is passed first so that
!     ** SUBSET_NONZERO_PATTERN can be used
      k = 3
      mats(1) = A
      mats(2) = Id
      mats(3) = B
      funs(1) = f2
      funs(2) = f1
      funs(3) = f3
      call NEPSetSplitOperator(nep,k,mats,funs,SUBSET_NONZERO_PATTERN,ierr);CHKERRA(ierr)
      call NEPSetProblemType(nep,NEP_GENERAL,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Customize nonlinear solver; set runtime options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      tol = 1e-9
      call NEPSetTolerances(nep,tol,PETSC_DEFAULT_INTEGER,ierr);CHKERRA(ierr)
      k = 1
      call NEPSetDimensions(nep,k,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr);CHKERRA(ierr)
      k = 0
      call NEPRIISetLagPreconditioner(nep,k,ierr);CHKERRA(ierr)

!     ** Set solver parameters at runtime
      call NEPSetFromOptions(nep,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call NEPSolve(nep,ierr);CHKERRA(ierr)

!     ** Optional: Get some information from the solver and display it
      call NEPGetType(nep,tname,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,120) tname
      endif
 120  format (' Solution method: ',A)
      call NEPGetDimensions(nep,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
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
        call NEPErrorView(nep,PEP_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr);CHKERRA(ierr)
      else
        call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr);CHKERRA(ierr)
        call NEPReasonView(nep,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
        call NEPErrorView(nep,PEP_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
        call PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
      endif
      call NEPDestroy(nep,ierr);CHKERRA(ierr)
      call MatDestroy(Id,ierr);CHKERRA(ierr)
      call MatDestroy(A,ierr);CHKERRA(ierr)
      call MatDestroy(B,ierr);CHKERRA(ierr)
      call FNDestroy(f1,ierr);CHKERRA(ierr)
      call FNDestroy(f2,ierr);CHKERRA(ierr)
      call FNDestroy(f3,ierr);CHKERRA(ierr)
      call SlepcFinalize(ierr)
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -terse
!      requires: !single
!
!TEST*/
