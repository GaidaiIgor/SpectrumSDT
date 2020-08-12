!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex23f90 [-help] [-t <t>] [-m <m>] [SLEPc opts]
!
!  Description: Computes exp(t*A)*v for a matrix associated with a
!  Markov model. This is the Fortran90 equivalent to ex23.c
!
!  The command line options are:
!    -t <t>, where <t> = time parameter (multiplies the matrix)
!    -m <m>, where <m> = number of grid subdivisions in each dimension
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
!  Variables:
!     A      problem matrix
!     mfn    matrix function solver context

      Mat            A
      MFN            mfn
      FN             f
      PetscReal      tol, norm, cst, pd, pu
      PetscScalar    t, z
      Vec            v, y
      PetscInt       N, m, ncv, maxit, its, ii, jj
      PetscInt       i, j, jmax, ix, Istart, Iend
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg
      MFNConvergedReason reason

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'SlepcInitialize failed'
        stop
      endif
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRA(ierr)
      m = 15
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-m',m,flg,ierr);CHKERRA(ierr)
      t = 2.0
      call PetscOptionsGetScalar(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-t',t,flg,ierr);CHKERRA(ierr)
      N = m*(m+1)/2
      if (rank .eq. 0) then
        write(*,100) N, m
      endif
 100  format (/'Markov y=exp(t*A)*e_1, N=',I6,' (m=',I4,')')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Compute the transition probability matrix, A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call MatCreate(PETSC_COMM_WORLD,A,ierr);CHKERRA(ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr);CHKERRA(ierr)
      call MatSetFromOptions(A,ierr);CHKERRA(ierr)
      call MatSetUp(A,ierr);CHKERRA(ierr)
      call MatGetOwnershipRange(A,Istart,Iend,ierr);CHKERRA(ierr)
      ix = 0
      cst = 0.5/real(m-1)
      do i=1,m
        jmax = m-i+1
        do j=1,jmax
          ix = ix + 1
          ii = ix - 1
          if (ix-1.ge.Istart .and. ix.le.Iend) then
            if (j.ne.jmax) then
              pd = cst*(i+j-1)
              !** north
              if (i.eq.1) then
                z = 2.0*pd
              else
                z = pd
              end if
              call MatSetValue(A,ii,ix,z,INSERT_VALUES,ierr);CHKERRA(ierr)
              !** east
              if (j.eq.1) then
                z = 2.0*pd
              else
                z = pd
              end if
              jj = ix+jmax-1
              call MatSetValue(A,ii,jj,z,INSERT_VALUES,ierr);CHKERRA(ierr)
            end if

            !** south
            pu = 0.5 - cst*(i+j-3)
            z = pu
            if (j.gt.1) then
              jj = ix-2
              call MatSetValue(A,ii,jj,z,INSERT_VALUES,ierr);CHKERRA(ierr)
            end if
            !** west
            if (i.gt.1) then
              jj = ix-jmax-2
              call MatSetValue(A,ii,jj,z,INSERT_VALUES,ierr);CHKERRA(ierr)
            end if
          end if
        end do
      end do
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

!     ** Set v = e_1
      call MatCreateVecs(A,y,v,ierr);CHKERRA(ierr)
      ii = 0
      z = 1.0
      call VecSetValue(v,ii,z,INSERT_VALUES,ierr);CHKERRA(ierr)
      call VecAssemblyBegin(v,ierr);CHKERRA(ierr)
      call VecAssemblyEnd(v,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create matrix function solver context
      call MFNCreate(PETSC_COMM_WORLD,mfn,ierr);CHKERRA(ierr)

!     ** Set operator matrix, the function to compute, and other options
      call MFNSetOperator(mfn,A,ierr);CHKERRA(ierr)
      call MFNGetFN(mfn,f,ierr);CHKERRA(ierr)
      call FNSetType(f,FNEXP,ierr);CHKERRA(ierr)
      z = 1.0
      call FNSetScale(f,t,z,ierr);CHKERRA(ierr)
      tol = 0.0000001
      call MFNSetTolerances(mfn,tol,PETSC_DEFAULT_INTEGER,ierr);CHKERRA(ierr)

!     ** Set solver parameters at runtime
      call MFNSetFromOptions(mfn,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the problem, y=exp(t*A)*v
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call MFNSolve(mfn,v,y,ierr);CHKERRA(ierr)
      call MFNGetConvergedReason(mfn,reason,ierr);CHKERRA(ierr)
      if (reason.lt.0) then; SETERRA(PETSC_COMM_WORLD,1,'Solver did not converge'); endif
      call VecNorm(y,NORM_2,norm,ierr);CHKERRA(ierr)

!     ** Optional: Get some information from the solver and display it
      call MFNGetIterationNumber(mfn,its,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,120) its
      endif
 120  format (' Number of iterations of the method: ',I4)
      call MFNGetDimensions(mfn,ncv,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,130) ncv
      endif
 130  format (' Subspace dimension:',I4)
      call MFNGetTolerances(mfn,tol,maxit,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,140) tol,maxit
      endif
 140  format (' Stopping condition: tol=',f10.7,' maxit=',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (rank .eq. 0) then
        write(*,150) PetscRealPart(t),norm
      endif
 150  format (' Computed vector at time t=',f4.1,' has norm ',f8.5)

      call MFNDestroy(mfn,ierr);CHKERRA(ierr)
      call MatDestroy(A,ierr);CHKERRA(ierr)
      call VecDestroy(v,ierr);CHKERRA(ierr)
      call VecDestroy(y,ierr);CHKERRA(ierr)
      call SlepcFinalize(ierr)
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -mfn_ncv 6
!
!TEST*/
