!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./test1f [-help]
!
!  Description: Simple example that tests BV interface functions.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcbv.h>
      use slepcbv
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define KMAX 35

      Vec            t,v
      Mat            Q,M
      BV             X,Y;
      PetscMPIInt    rank
      PetscInt       i,j,n,k,l,izero,ione
      PetscScalar    qq(1),z(KMAX),val
      PetscScalar    one,mone,two,zero
      PetscOffset    iq
      PetscReal      nrm
      PetscBool      flg
      PetscErrorCode ierr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      n = 10
      k = 5
      l = 3
      one = 1.0
      mone = -1.0
      two = 2.0
      zero = 0.0
      izero = 0
      ione = 1
      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'SlepcInitialize failed'
        stop
      endif
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRA(ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr);CHKERRA(ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-k',k,flg,ierr);CHKERRA(ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-l',l,flg,ierr);CHKERRA(ierr)
      if (k .gt. KMAX) then; SETERRA(PETSC_COMM_SELF,1,'Program currently limited to k=35'); endif
      if (rank .eq. 0) then
        write(*,110) k,n
      endif
 110  format (/'Test BV with',I3,' columns of length',I3,' (Fortran)')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initialize data
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create template vector
      call VecCreate(PETSC_COMM_WORLD,t,ierr);CHKERRA(ierr)
      call VecSetSizes(t,PETSC_DECIDE,n,ierr);CHKERRA(ierr)
      call VecSetFromOptions(t,ierr);CHKERRA(ierr)

!     ** Create BV object X
      call BVCreate(PETSC_COMM_WORLD,X,ierr);CHKERRA(ierr)
      call PetscObjectSetName(X,'X',ierr);CHKERRA(ierr)
      call BVSetSizesFromVec(X,t,k,ierr);CHKERRA(ierr)
      call BVSetFromOptions(X,ierr);CHKERRA(ierr)

!     ** Fill X entries
      do j=0,k-1
        call BVGetColumn(X,j,v,ierr);CHKERRA(ierr)
        call VecSet(v,zero,ierr);CHKERRA(ierr)
        do i=0,3
          if (i+j<n) then
            val = 3*i+j-2
            call VecSetValue(v,i+j,val,INSERT_VALUES,ierr);CHKERRA(ierr)
          end if
        end do
        call VecAssemblyBegin(v,ierr);CHKERRA(ierr)
        call VecAssemblyEnd(v,ierr);CHKERRA(ierr)
        call BVRestoreColumn(X,j,v,ierr);CHKERRA(ierr)
      end do

!     ** Create BV object Y
      call BVCreate(PETSC_COMM_WORLD,Y,ierr);CHKERRA(ierr)
      call PetscObjectSetName(Y,'Y',ierr);CHKERRA(ierr)
      call BVSetSizesFromVec(Y,t,l,ierr);CHKERRA(ierr)
      call BVSetFromOptions(Y,ierr);CHKERRA(ierr)

!     ** Fill Y entries
      do j=0,l-1
        call BVGetColumn(Y,j,v,ierr);CHKERRA(ierr)
        val = real(j+1)/4.0
        call VecSet(v,val,ierr);CHKERRA(ierr)
        call BVRestoreColumn(Y,j,v,ierr);CHKERRA(ierr)
      end do

!     ** Create Mat
      call MatCreateSeqDense(PETSC_COMM_SELF,k,l,PETSC_NULL_SCALAR,Q,ierr);CHKERRA(ierr)
      call PetscObjectSetName(Q,'Q',ierr);CHKERRA(ierr)
      call MatDenseGetArray(Q,qq,iq,ierr);CHKERRA(ierr)
      do i=0,k-1
        do j=0,l-1
          if (i<j) then
            qq(iq+1+i+j*k) = 2.0
          else
            qq(iq+1+i+j*k) = -0.5
          end if
        end do
      end do
      call MatDenseRestoreArray(Q,qq,iq,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Test several operations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Test BVMult
      call BVMult(Y,two,one,X,Q,ierr);CHKERRA(ierr)

!     ** Test BVMultVec
      call BVGetColumn(Y,izero,v,ierr);CHKERRA(ierr)
      z(1) = 2.0
      do i=2,k
        z(i) = -0.5*z(i-1)
      end do
      call BVMultVec(X,mone,one,v,z,ierr);CHKERRA(ierr)
      call BVRestoreColumn(Y,izero,v,ierr);CHKERRA(ierr)

!     ** Test BVDot
      call MatCreateSeqDense(PETSC_COMM_SELF,l,k,PETSC_NULL_SCALAR,M,ierr);CHKERRA(ierr)
      call PetscObjectSetName(M,'M',ierr);CHKERRA(ierr)
      call BVDot(X,Y,M,ierr);CHKERRA(ierr)

!     ** Test BVDotVec
      call BVGetColumn(Y,izero,v,ierr);CHKERRA(ierr)
      call BVDotVec(X,v,z,ierr);CHKERRA(ierr)
      call BVRestoreColumn(Y,izero,v,ierr);CHKERRA(ierr)

!     ** Test BVMultInPlace and BVScale
      call BVMultInPlace(X,Q,ione,l,ierr);CHKERRA(ierr)
      call BVScale(X,two,ierr);CHKERRA(ierr)

!     ** Test BVNorm
      call BVNormColumn(X,izero,NORM_2,nrm,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,120) nrm
      endif
 120  format ('2-Norm of X[0] = ',f8.4)
      call BVNorm(X,NORM_FROBENIUS,nrm,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,130) nrm
      endif
 130  format ('Frobenius Norm of X = ',f8.4)

!     *** Clean up
      call BVDestroy(X,ierr);CHKERRA(ierr)
      call BVDestroy(Y,ierr);CHKERRA(ierr)
      call VecDestroy(t,ierr);CHKERRA(ierr)
      call MatDestroy(Q,ierr);CHKERRA(ierr)
      call MatDestroy(M,ierr);CHKERRA(ierr)
      call SlepcFinalize(ierr)
      end

!/*TEST
!
!   test:
!      suffix: 1
!      nsize: 1
!      args: -bv_type {{vecs contiguous svec mat}separate output}
!      output_file: output/test1f_1.out
!
!   test:
!      suffix: 2
!      nsize: 2
!      args: -bv_type {{vecs contiguous svec mat}separate output}
!      output_file: output/test1f_1.out
!
!TEST*/
