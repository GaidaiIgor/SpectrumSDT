!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex27f90 [-help] [-n <n>] [all SLEPc options]
!
!  Description: Simple NLEIGS example. Fortran90 equivalent of ex27.c
!
!  The command line options are:
!    -n <n>, where <n> = matrix dimension
!
! ----------------------------------------------------------------------
!   Solve T(lambda)x=0 using NLEIGS solver
!      with T(lambda) = -D+sqrt(lambda)*I
!      where D is the Laplacian operator in 1 dimension
!      and with the interpolation interval [.01,16]
! ----------------------------------------------------------------------
!
PROGRAM main
#include <slepc/finclude/slepcnep.h>
  USE slepcnep
  implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  NEP                :: nep
  Mat                :: A(2),F,J
  NEPType            :: ntype
  PetscInt           :: n=100,nev,Istart,Iend,i,col,one,two,three
  PetscErrorCode     :: ierr
  PetscBool          :: terse,flg,split=PETSC_TRUE
  PetscReal          :: ia,ib,ic,id
  RG                 :: rg
  FN                 :: fn(2)
  PetscScalar        :: coeffs,sigma,done
  CHARACTER(LEN=128) :: string

  ! NOTE: Any user-defined Fortran routines (such as ComputeSingularities)
  !       MUST be declared as external.
  external ComputeSingularities, FormFunction, FormJacobian

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    print*,'SlepcInitialize failed'
    stop
  end if
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-n",n,flg,ierr);CHKERRA(ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-split",split,flg,ierr);CHKERRA(ierr)
  if (split) then
     write(string,*) 'Square root eigenproblem, n=',n,' (split-form)\n'
  else
     write(string,*) 'Square root eigenproblem, n=',n,'\n'
  end if
  call PetscPrintf(PETSC_COMM_WORLD,trim(string),ierr);CHKERRA(ierr)
  done  = 1.0
  one   = 1
  two   = 2
  three = 3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create nonlinear eigensolver context and set options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  call NEPCreate(PETSC_COMM_WORLD,nep,ierr);CHKERRA(ierr)
  call NEPSetType(nep,NEPNLEIGS,ierr);CHKERRA(ierr)
  call NEPNLEIGSSetSingularitiesFunction(nep,ComputeSingularities,0,ierr);CHKERRA(ierr)
  call NEPGetRG(nep,rg,ierr);CHKERRA(ierr)
  call RGSetType(rg,RGINTERVAL,ierr);CHKERRA(ierr)
  ia = 0.01
  ib = 16.0
#if defined(PETSC_USE_COMPLEX)
  ic = -0.001
  id = 0.001
#else
  ic = 0.0
  id = 0.0
#endif
  call RGIntervalSetEndpoints(rg,ia,ib,ic,id,ierr);CHKERRA(ierr)
  sigma = 1.1
  call NEPSetTarget(nep,sigma,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Define the nonlinear problem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (split) then
     ! ** Create matrices for the split form
     call MatCreate(PETSC_COMM_WORLD,A(1),ierr);CHKERRA(ierr)
     call MatSetSizes(A(1),PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
     call MatSetFromOptions(A(1),ierr);CHKERRA(ierr)
     call MatSetUp(A(1),ierr);CHKERRA(ierr)
     call MatGetOwnershipRange(A(1),Istart,Iend,ierr);CHKERRA(ierr)
     coeffs = -2.0
     do i=Istart,Iend-1
        if (i.gt.0) then
           col = i-1
           call MatSetValue(A(1),i,col,done,INSERT_VALUES,ierr);CHKERRA(ierr)
        end if
        if (i.lt.n-1) then
           col = i+1
           call MatSetValue(A(1),i,col,done,INSERT_VALUES,ierr);CHKERRA(ierr)
        end if
        call MatSetValue(A(1),i,i,coeffs,INSERT_VALUES,ierr);CHKERRA(ierr)
     end do
     call MatAssemblyBegin(A(1),MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
     call MatAssemblyEnd(A(1),MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

     call MatCreate(PETSC_COMM_WORLD,A(2),ierr);CHKERRA(ierr)
     call MatSetSizes(A(2),PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
     call MatSetFromOptions(A(2),ierr);CHKERRA(ierr)
     call MatSetUp(A(2),ierr);CHKERRA(ierr)
     call MatAssemblyBegin(A(2),MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
     call MatAssemblyEnd(A(2),MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
     call MatShift(A(2),done,ierr);CHKERRA(ierr)

     ! ** Define functions for the split form
     call FNCreate(PETSC_COMM_WORLD,fn(1),ierr);CHKERRA(ierr)
     call FNSetType(fn(1),FNRATIONAL,ierr);CHKERRA(ierr)
     call FNRationalSetNumerator(fn(1),one,done,ierr);CHKERRA(ierr)
     call FNCreate(PETSC_COMM_WORLD,fn(2),ierr);CHKERRA(ierr)
     call FNSetType(fn(2),FNSQRT,ierr);CHKERRA(ierr)
     call NEPSetSplitOperator(nep,two,A,fn,SUBSET_NONZERO_PATTERN,ierr);CHKERRA(ierr)
  else
    ! ** Callback form: create matrix and set Function evaluation routine
    call MatCreate(PETSC_COMM_WORLD,F,ierr);CHKERRA(ierr)
    call MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
    call MatSetFromOptions(F,ierr);CHKERRA(ierr)
    call MatSeqAIJSetPreallocation(F,three,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
    call MatMPIAIJSetPreallocation(F,three,PETSC_NULL_INTEGER,one,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
    Call MatSetUp(F,ierr);CHKERRA(ierr)
    call NEPSetFunction(nep,F,F,FormFunction,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)

    call MatCreate(PETSC_COMM_WORLD,J,ierr);CHKERRA(ierr)
    call MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
    call MatSetFromOptions(J,ierr);CHKERRA(ierr)
    call MatSeqAIJSetPreallocation(J,one,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
    call MatMPIAIJSetPreallocation(J,one,PETSC_NULL_INTEGER,one,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
    call MatSetUp(J,ierr);CHKERRA(ierr)
    call NEPSetJacobian(nep,J,FormJacobian,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
  end if

  call NEPSetFromOptions(nep,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  call NEPSolve(nep,ierr);CHKERRA(ierr)
  call NEPGetType(nep,ntype,ierr);CHKERRA(ierr)
  write(string,*) 'Solution method: ',ntype,'\n'
  call PetscPrintf(PETSC_COMM_WORLD,trim(string),ierr);CHKERRA(ierr)
  call NEPGetDimensions(nep,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
  write(string,*) 'Number of requested eigenvalues:',nev,'\n'
  call PetscPrintf(PETSC_COMM_WORLD,trim(string),ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! ** show detailed info unless -terse option is given by user
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-terse',terse,ierr);CHKERRA(ierr)
  if (terse) then
    call NEPErrorView(nep,NEP_ERROR_BACKWARD,PETSC_NULL_VIEWER,ierr);CHKERRA(ierr)
  else
    call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr);CHKERRA(ierr)
    call NEPReasonView(nep,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
    call NEPErrorView(nep,NEP_ERROR_BACKWARD,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
    call PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRA(ierr)
  end if

  if (split) then
    call MatDestroy(A(1),ierr);CHKERRA(ierr)
    call MatDestroy(A(2),ierr);CHKERRA(ierr)
    call FNDestroy(fn(1),ierr);CHKERRA(ierr)
    call FNDestroy(fn(2),ierr);CHKERRA(ierr)
  else
    call MatDestroy(F,ierr);CHKERRA(ierr)
    call MatDestroy(J,ierr);CHKERRA(ierr)
  end if
  call NEPDestroy(nep,ierr)
  call SlepcFinalize(ierr)

END PROGRAM main

! --------------------------------------------------------------
!
!   FormFunction - Computes Function matrix  T(lambda)
!
SUBROUTINE FormFunction(nep,lambda,fun,B,ctx,ierr)
#include <slepc/finclude/slepcnep.h>
  use slepcnep
  implicit none

  NEP            :: nep
  PetscScalar    :: lambda,val(0:2),t
  Mat            :: fun,B
  PetscInt       :: ctx,i,n,col(0:2),Istart,Iend,Istart0,Iend0,one,two,three
  PetscErrorCode :: ierr
  PetscBool      :: FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE

  one   = 1
  two   = 2
  three = 3

  ! ** Compute Function entries and insert into matrix
  t = sqrt(lambda)
  call MatGetSize(fun,n,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
  call MatGetOwnershipRange(fun,Istart,Iend,ierr);CHKERRA(ierr)
  if (Istart.eq.0) FirstBlock=PETSC_TRUE;
  if (Iend.eq.n) LastBlock=PETSC_TRUE;
  val(0)=1.0; val(1)=t-2.0; val(2)=1.0;

  Istart0 = Istart
  if (FirstBlock) Istart0 = Istart+1
  Iend0 = Iend
  if (LastBlock) Iend0 = Iend-1

  do i=Istart0,Iend0-1
     col(0) = i-1
     col(1) = i
     col(2) = i+1
     call MatSetValues(fun,one,i,three,col,val,INSERT_VALUES,ierr);CHKERRA(ierr)
  end do

  if (LastBlock) then
     i = n-1
     col(0) = n-2
     col(1) = n-1
     val(0) = 1.0
     val(1) = t-2.0
     call MatSetValues(fun,one,i,two,col,val,INSERT_VALUES,ierr);CHKERRA(ierr)
  end if

  if (FirstBlock) then
     i = 0
     col(0) = 0
     col(1) = 1
     val(0) = t-2.0
     val(1) = 1.0
     call MatSetValues(fun,one,i,two,col,val,INSERT_VALUES,ierr);CHKERRA(ierr)
  end if

  ! ** Assemble matrix
  call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(fun,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(fun,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

END SUBROUTINE FormFunction

! --------------------------------------------------------------
!
!   FormJacobian - Computes Jacobian matrix  T'(lambda)
!
SUBROUTINE FormJacobian(nep,lambda,jac,ctx,ierr)
#include <slepc/finclude/slepcnep.h>
  USE slepcnep
  implicit none

  NEP            :: nep
  PetscScalar    :: lambda,t
  Mat            :: jac
  PetscInt       :: ctx
  PetscErrorCode :: ierr
  Vec            :: d

  call MatCreateVecs(jac,d,PETSC_NULL_VEC,ierr);CHKERRA(ierr)
  t = 0.5/sqrt(lambda)
  call VecSet(d,t,ierr);CHKERRA(ierr)
  call MatDiagonalSet(jac,d,INSERT_VALUES,ierr);CHKERRA(ierr)
  calL VecDestroy(d,ierr);CHKERRA(ierr)

END SUBROUTINE FormJacobian

! --------------------------------------------------------------
!
!  ComputeSingularities - This is a user-defined routine to compute maxnp
!  points (at most) in the complex plane where the function T(.) is not analytic.
!
!  In this case, we discretize the singularity region (-inf,0)~(-10e+6,-10e-6)
!
!  Input Parameters:
!    nep   - nonlinear eigensolver context
!    maxnp - on input number of requested points in the discretization (can be set)
!    xi    - computed values of the discretization
!    dummy - optional user-defined monitor context (unused here)
!
SUBROUTINE ComputeSingularities(nep,maxnp,xi,dummy,ierr)
#include <slepc/finclude/slepcnep.h>
  use slepcnep
  implicit none

  NEP            :: nep
  PetscInt       :: maxnp, dummy
  PetscScalar    :: xi(0:maxnp-1)
  PetscErrorCode :: ierr
  PetscReal      :: h
  PetscInt       :: i

  h = 11.0/real(maxnp-1)
  xi(0) = -1e-5
  xi(maxnp-1) = -1e+6
  do i=1,maxnp-2
     xi(i) = -10**(-5+h*i)
  end do

END SUBROUTINE ComputeSingularities

!/*TEST
!
!   test:
!      suffix: 1
!      args: -nep_nev 3 -nep_nleigs_interpolation_degree 90 -terse
!      requires: !single
!
!   test:
!      suffix: 2
!      args: -split 0 -nep_nev 3 -nep_nleigs_interpolation_degree 90 -terse
!      requires: !single
!
!TEST*/
