!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Description: Simple example to test the EPS Fortran interface.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepceps.h>
      use slepceps
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Mat                A,B
      EPS                eps
      ST                 st
      KSP                ksp
      DS                 ds
      PetscReal          cut,tol,tolabs
      PetscScalar        tget,value
      PetscInt           n,i,its,Istart,Iend
      PetscInt           nev,ncv,mpd
      PetscBool          flg
      EPSConvergedReason reason
      EPSType            tname
      EPSExtraction      extr
      EPSBalance         bal
      EPSWhich           which
      EPSConv            conv
      EPSStop            stp
      EPSProblemType     ptype
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
 100  format (/'Diagonal Eigenproblem, n =',I3,' (Fortran)')

      call MatCreate(PETSC_COMM_WORLD,A,ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(A,ierr)
      call MatSetUp(A,ierr)
      call MatGetOwnershipRange(A,Istart,Iend,ierr)
      do i=Istart,Iend-1
        value = i+1
        call MatSetValue(A,i,i,value,INSERT_VALUES,ierr)
      enddo
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create eigensolver and test interface functions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call EPSCreate(PETSC_COMM_WORLD,eps,ierr)
      call EPSSetOperators(eps,A,PETSC_NULL_MAT,ierr)
      call EPSGetOperators(eps,B,PETSC_NULL_MAT,ierr)
      call MatView(B,PETSC_NULL_VIEWER,ierr)

      call EPSSetType(eps,EPSKRYLOVSCHUR,ierr)
      call EPSGetType(eps,tname,ierr)
      if (rank .eq. 0) then
        write(*,110) tname
      endif
 110  format (' Type set to ',A)

      call EPSGetProblemType(eps,ptype,ierr)
      if (rank .eq. 0) then
        write(*,120) ptype
      endif
 120  format (' Problem type before changing = ',I2)
      call EPSSetProblemType(eps,EPS_HEP,ierr)
      call EPSGetProblemType(eps,ptype,ierr)
      if (rank .eq. 0) then
        write(*,130) ptype
      endif
 130  format (' ... changed to ',I2)
      call EPSIsGeneralized(eps,flg,ierr)
      if (flg .and. rank .eq. 0) then
        write(*,*) 'generalized'
      endif
      call EPSIsHermitian(eps,flg,ierr)
      if (flg .and. rank .eq. 0) then
        write(*,*) 'hermitian'
      endif
      call EPSIsPositive(eps,flg,ierr)
      if (flg .and. rank .eq. 0) then
        write(*,*) 'positive'
      endif

      call EPSGetExtraction(eps,extr,ierr)
      if (rank .eq. 0) then
        write(*,140) extr
      endif
 140  format (' Extraction before changing = ',I2)
      call EPSSetExtraction(eps,EPS_HARMONIC,ierr)
      call EPSGetExtraction(eps,extr,ierr)
      if (rank .eq. 0) then
        write(*,150) extr
      endif
 150  format (' ... changed to ',I2)

      its = 8
      cut = 1.0e-6
      bal = EPS_BALANCE_ONESIDE
      call EPSSetBalance(eps,bal,its,cut,ierr)
      call EPSGetBalance(eps,bal,its,cut,ierr)
      if (rank .eq. 0) then
        write(*,160) bal,its,cut
      endif
 160  format (' Balance: ',I2,', its=',I2,', cutoff=',F9.6)

      tget = 4.8
      call EPSSetTarget(eps,tget,ierr)
      call EPSGetTarget(eps,tget,ierr)
      call EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE,ierr)
      call EPSGetWhichEigenpairs(eps,which,ierr)
      if (rank .eq. 0) then
        write(*,170) which,PetscRealPart(tget)
      endif
 170  format (' Which = ',I2,', target = ',F4.1)

      nev = 4
      call EPSSetDimensions(eps,nev,PETSC_DEFAULT_INTEGER,              &
     &                      PETSC_DEFAULT_INTEGER,ierr)
      call EPSGetDimensions(eps,nev,ncv,mpd,ierr)
      if (rank .eq. 0) then
        write(*,180) nev,ncv,mpd
      endif
 180  format (' Dimensions: nev=',I2,', ncv=',I2,', mpd=',I2)

      tol = 2.2e-4
      its = 200
      call EPSSetTolerances(eps,tol,its,ierr)
      call EPSGetTolerances(eps,tol,its,ierr)
      if (rank .eq. 0) then
        write(*,190) tol,its
      endif
 190  format (' Tolerance =',F8.5,', max_its =',I4)

      call EPSSetConvergenceTest(eps,EPS_CONV_ABS,ierr)
      call EPSGetConvergenceTest(eps,conv,ierr)
      call EPSSetStoppingTest(eps,EPS_STOP_BASIC,ierr)
      call EPSGetStoppingTest(eps,stp,ierr)
      if (rank .eq. 0) then
        write(*,200) conv,stp
      endif
 200  format (' Convergence test =',I2,', stopping test =',I2)

      call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,        &
     &                   PETSC_VIEWER_DEFAULT,vf,ierr)
      call EPSMonitorSet(eps,EPSMONITORFIRST,vf,                        &
     &                   PetscViewerAndFormatDestroy,ierr)
      call SlepcConvMonitorCreate(PETSC_VIEWER_STDOUT_WORLD,            &
     &                   PETSC_VIEWER_DEFAULT,ctx,ierr)
      call EPSMonitorSet(eps,EPSMONITORCONVERGED,ctx,                   &
     &                   SlepcConvMonitorDestroy,ierr)
      call EPSMonitorCancel(eps,ierr)

      call EPSGetST(eps,st,ierr)
      call STGetKSP(st,ksp,ierr)
      tol = 1.e-8
      tolabs = 1.e-35
      call KSPSetTolerances(ksp,tol,tolabs,PETSC_DEFAULT_REAL,          &
     &                      PETSC_DEFAULT_INTEGER,ierr)
      call STView(st,PETSC_NULL_VIEWER,ierr)
      call EPSGetDS(eps,ds,ierr)
      call DSView(ds,PETSC_NULL_VIEWER,ierr)

      call EPSSetFromOptions(eps,ierr)
      call EPSSolve(eps,ierr)
      call EPSGetConvergedReason(eps,reason,ierr)
      call EPSGetIterationNumber(eps,its,ierr)
      if (rank .eq. 0) then
        write(*,210) reason,its
      endif
 210  format (' Finished - converged reason =',I2,', its=',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_NULL_VIEWER,ierr)
      call EPSDestroy(eps,ierr)
      call MatDestroy(A,ierr)

      call SlepcFinalize(ierr)
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -eps_ncv 14
!
!TEST*/
