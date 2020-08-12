!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Description: Simple example to test the NEP Fortran interface.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcnep.h>
      use slepcnep
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Mat                A(3),B
      FN                 f(3),g
      NEP                nep
      DS                 ds
      RG                 rg
      PetscReal          tol
      PetscScalar        coeffs(2),tget,val
      PetscInt           n,i,its,Istart,Iend
      PetscInt           nev,ncv,mpd,nterm
      PetscInt           nc,np
      NEPWhich           which
      NEPConvergedReason reason
      NEPType            tname
      NEPRefine          refine
      NEPRefineScheme    rscheme
      NEPConv            conv
      NEPStop            stp
      NEPProblemType     ptype
      MatStructure       mstr
      PetscMPIInt        rank
      PetscErrorCode     ierr
      SlepcConvMonitor   ctx
      PetscViewerAndFormat vf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      n = 20
      if (rank .eq. 0) then
        write(*,100) n
      endif
 100  format (/'Diagonal Nonlinear Eigenproblem, n =',I3,' (Fortran)')

!     Matrices
      call MatCreate(PETSC_COMM_WORLD,A(1),ierr)
      call MatSetSizes(A(1),PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(A(1),ierr)
      call MatSetUp(A(1),ierr)
      call MatGetOwnershipRange(A(1),Istart,Iend,ierr)
      do i=Istart,Iend-1
        val = i+1
        call MatSetValue(A(1),i,i,val,INSERT_VALUES,ierr)
      enddo
      call MatAssemblyBegin(A(1),MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A(1),MAT_FINAL_ASSEMBLY,ierr)

      call MatCreate(PETSC_COMM_WORLD,A(2),ierr)
      call MatSetSizes(A(2),PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(A(2),ierr)
      call MatSetUp(A(2),ierr)
      call MatGetOwnershipRange(A(2),Istart,Iend,ierr)
      do i=Istart,Iend-1
        val = 1
        call MatSetValue(A(2),i,i,val,INSERT_VALUES,ierr)
      enddo
      call MatAssemblyBegin(A(2),MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A(2),MAT_FINAL_ASSEMBLY,ierr)

      call MatCreate(PETSC_COMM_WORLD,A(3),ierr)
      call MatSetSizes(A(3),PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(A(3),ierr)
      call MatSetUp(A(3),ierr)
      call MatGetOwnershipRange(A(3),Istart,Iend,ierr)
      do i=Istart,Iend-1
        val = real(n)/real(i+1)
        call MatSetValue(A(3),i,i,val,INSERT_VALUES,ierr)
      enddo
      call MatAssemblyBegin(A(3),MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A(3),MAT_FINAL_ASSEMBLY,ierr)

!     Functions: f0=-lambda, f1=1.0, f2=sqrt(lambda)
      call FNCreate(PETSC_COMM_WORLD,f(1),ierr)
      call FNSetType(f(1),FNRATIONAL,ierr)
      nc = 2
      coeffs(1) = -1.0
      coeffs(2) = 0.0
      call FNRationalSetNumerator(f(1),nc,coeffs,ierr)

      call FNCreate(PETSC_COMM_WORLD,f(2),ierr)
      call FNSetType(f(2),FNRATIONAL,ierr)
      nc = 1
      coeffs(1) = 1.0
      call FNRationalSetNumerator(f(2),nc,coeffs,ierr)

      call FNCreate(PETSC_COMM_WORLD,f(3),ierr)
      call FNSetType(f(3),FNSQRT,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create eigensolver and test interface functions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call NEPCreate(PETSC_COMM_WORLD,nep,ierr)
      nterm = 3
      mstr = SAME_NONZERO_PATTERN
      call NEPSetSplitOperator(nep,nterm,A,f,mstr,ierr)
      call NEPGetSplitOperatorInfo(nep,nterm,mstr,ierr)
      if (rank .eq. 0) then
        write(*,110) nterm
      endif
 110  format (' Nonlinear function with ',I2,' terms')
      i = 0
      call NEPGetSplitOperatorTerm(nep,i,B,g,ierr)
      call MatView(B,PETSC_NULL_VIEWER,ierr)
      call FNView(g,PETSC_NULL_VIEWER,ierr)

      call NEPSetType(nep,NEPRII,ierr)
      call NEPGetType(nep,tname,ierr)
      if (rank .eq. 0) then
        write(*,120) tname
      endif
 120  format (' Type set to ',A)

      call NEPGetProblemType(nep,ptype,ierr)
      if (rank .eq. 0) then
        write(*,130) ptype
      endif
 130  format (' Problem type before changing = ',I2)
      call NEPSetProblemType(nep,NEP_RATIONAL,ierr)
      call NEPGetProblemType(nep,ptype,ierr)
      if (rank .eq. 0) then
        write(*,140) ptype
      endif
 140  format (' ... changed to ',I2)

      np = 1
      tol = 1e-9
      its = 2
      refine = NEP_REFINE_SIMPLE
      rscheme = NEP_REFINE_SCHEME_EXPLICIT
      call NEPSetRefine(nep,refine,np,tol,its,rscheme,ierr)
      call NEPGetRefine(nep,refine,np,tol,its,rscheme,ierr)
      if (rank .eq. 0) then
        write(*,190) refine,tol,its,rscheme
      endif
 190  format (' Refinement: ',I2,', tol=',F12.9,', its=',I2', scheme=', &
     &   I2)

      tget = 1.1
      call NEPSetTarget(nep,tget,ierr)
      call NEPGetTarget(nep,tget,ierr)
      call NEPSetWhichEigenpairs(nep,NEP_TARGET_MAGNITUDE,ierr)
      call NEPGetWhichEigenpairs(nep,which,ierr)
      if (rank .eq. 0) then
        write(*,200) which,PetscRealPart(tget)
      endif
 200  format (' Which = ',I2,', target = ',F4.1)

      nev = 1
      ncv = 12
      call NEPSetDimensions(nep,nev,ncv,PETSC_DEFAULT_INTEGER,ierr)
      call NEPGetDimensions(nep,nev,ncv,mpd,ierr)
      if (rank .eq. 0) then
        write(*,210) nev,ncv,mpd
      endif
 210  format (' Dimensions: nev=',I2,', ncv=',I2,', mpd=',I2)

      tol = 1.0e-6
      its = 200
      call NEPSetTolerances(nep,tol,its,ierr)
      call NEPGetTolerances(nep,tol,its,ierr)
      if (rank .eq. 0) then
        write(*,220) tol,its
      endif
 220  format (' Tolerance =',F9.6,', max_its =',I4)

      call NEPSetConvergenceTest(nep,NEP_CONV_ABS,ierr)
      call NEPGetConvergenceTest(nep,conv,ierr)
      call NEPSetStoppingTest(nep,NEP_STOP_BASIC,ierr)
      call NEPGetStoppingTest(nep,stp,ierr)
      if (rank .eq. 0) then
        write(*,230) conv,stp
      endif
 230  format (' Convergence test =',I2,', stopping test =',I2)

      call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,        &
     &                   PETSC_VIEWER_DEFAULT,vf,ierr)
      call NEPMonitorSet(nep,NEPMONITORFIRST,vf,                        &
     &                   PetscViewerAndFormatDestroy,ierr)
      call SlepcConvMonitorCreate(PETSC_VIEWER_STDOUT_WORLD,            &
     &                   PETSC_VIEWER_DEFAULT,ctx,ierr)
      call NEPMonitorSet(nep,NEPMONITORCONVERGED,ctx,                   &
     &                   SlepcConvMonitorDestroy,ierr)
      call NEPMonitorCancel(nep,ierr)

      call NEPGetDS(nep,ds,ierr)
      call DSView(ds,PETSC_NULL_VIEWER,ierr)
      call NEPSetFromOptions(nep,ierr)

      call NEPGetRG(nep,rg,ierr)
      call RGView(rg,PETSC_NULL_VIEWER,ierr)

      call NEPSolve(nep,ierr)
      call NEPGetConvergedReason(nep,reason,ierr)
      if (rank .eq. 0) then
        write(*,240) reason
      endif
 240  format (' Finished - converged reason =',I2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call NEPErrorView(nep,NEP_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr)
      call NEPDestroy(nep,ierr)
      call MatDestroy(A(1),ierr)
      call MatDestroy(A(2),ierr)
      call MatDestroy(A(3),ierr)
      call FNDestroy(f(1),ierr)
      call FNDestroy(f(2),ierr)
      call FNDestroy(f(3),ierr)

      call SlepcFinalize(ierr)
      end

!/*TEST
!
!   test:
!      suffix: 1
!      requires: !single
!
!TEST*/
