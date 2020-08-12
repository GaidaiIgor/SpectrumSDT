!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Program usage: mpiexec -n <np> ./ex20f90 [-n <n>] [SLEPc opts]
!
!  Description: Simple 1-D nonlinear eigenproblem. Fortran90 equivalent of ex20.c
!
!  The command line options are:
!    -n <n>, where <n> = number of grid subdivisions
!
! ----------------------------------------------------------------------
!  Solve 1-D PDE
!           -u'' = lambda*u
!  on [0,1] subject to
!           u(0)=0, u'(1)=u(1)*lambda*kappa/(kappa-lambda)
! ----------------------------------------------------------------------
!

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     User-defined application context
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      module UserModule
#include <slepc/finclude/slepcnep.h>
      use slepcnep
      type User
        PetscScalar kappa
        PetscReal   h
      end type User
      end module

      program main
#include <slepc/finclude/slepcnep.h>
      use UserModule
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     nep       nonlinear eigensolver context
!     x         eigenvector
!     lambda    eigenvalue
!     F,J       Function and Jacobian matrices
!     ctx       user-defined context

      NEP            nep
      Vec            x, v(1)
      PetscScalar    lambda
      Mat            F, J
      type(User)     ctx
      NEPType        tname
      PetscInt       n, i, k, nev, its, maxit, nconv, three, one
      PetscReal      tol, norm
      PetscScalar    alpha
      PetscMPIInt    rank
      PetscBool      flg
      PetscErrorCode ierr
!  Note: Any user-defined Fortran routines (such as FormJacobian)
!  MUST be declared as external.
      external       FormFunction, FormJacobian

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
      if (rank .eq. 0) then
        write(*,'(/A,I4)') 'Nonlinear Eigenproblem, n =',n
      endif

      ctx%h = 1.0/real(n)
      ctx%kappa = 1.0

      three = 3
      one = 1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create matrix data structure to hold the Function and the Jacobian
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call MatCreate(PETSC_COMM_WORLD,F,ierr);CHKERRA(ierr)
      call MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
      call MatSetFromOptions(F,ierr);CHKERRA(ierr)
      call MatSeqAIJSetPreallocation(F,three,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      call MatMPIAIJSetPreallocation(F,three,PETSC_NULL_INTEGER,one,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      call MatSetUp(F,ierr);CHKERRA(ierr)

      call MatCreate(PETSC_COMM_WORLD,J,ierr);CHKERRA(ierr)
      call MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr);CHKERRA(ierr)
      call MatSetFromOptions(J,ierr);CHKERRA(ierr)
      call MatSeqAIJSetPreallocation(J,three,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      call MatMPIAIJSetPreallocation(J,three,PETSC_NULL_INTEGER,one,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      call MatSetUp(J,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Create the eigensolver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Create eigensolver context
      call NEPCreate(PETSC_COMM_WORLD,nep,ierr);CHKERRA(ierr)

!     ** Set routines for evaluation of Function and Jacobian
      call NEPSetFunction(nep,F,F,FormFunction,ctx,ierr);CHKERRA(ierr)
      call NEPSetJacobian(nep,J,FormJacobian,ctx,ierr);CHKERRA(ierr)

!     ** Customize nonlinear solver
      tol = 1e-9
      call NEPSetTolerances(nep,tol,PETSC_DEFAULT_INTEGER,ierr);CHKERRA(ierr)
      k = 1
      call NEPSetDimensions(nep,k,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr);CHKERRA(ierr)

!     ** Set solver parameters at runtime
      call NEPSetFromOptions(nep,ierr);CHKERRA(ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     ** Evaluate initial guess
      call MatCreateVecs(F,x,PETSC_NULL_VEC,ierr);CHKERRA(ierr)
      call VecDuplicate(x,v(1),ierr);CHKERRA(ierr)
      alpha = 1.0
      call VecSet(v(1),alpha,ierr);CHKERRA(ierr)
      k = 1
      call NEPSetInitialSpace(nep,k,v,ierr);CHKERRA(ierr)
      call VecDestroy(v(1),ierr);CHKERRA(ierr)

!     ** Call the solver
      call NEPSolve(nep,ierr);CHKERRA(ierr)
      call NEPGetIterationNumber(nep,its,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,'(A,I3)') ' Number of NEP iterations =',its
      endif

!     ** Optional: Get some information from the solver and display it
      call NEPGetType(nep,tname,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,'(A,A10)') ' Solution method: ',tname
      endif
      call NEPGetDimensions(nep,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,'(A,I4)') ' Number of requested eigenvalues:',nev
      endif
      call NEPGetTolerances(nep,tol,maxit,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,'(A,F12.9,A,I5)') ' Stopping condition: tol=',tol,', maxit=',maxit
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call NEPGetConverged(nep,nconv,ierr);CHKERRA(ierr)
      if (rank .eq. 0) then
        write(*,'(A,I2/)') ' Number of converged approximate eigenpairs:',nconv
      endif

!     ** Display eigenvalues and relative errors
      if (nconv .gt. 0) then
        if (rank .eq. 0) then
          write(*,*) '        k              ||T(k)x||'
          write(*,*) '----------------- ------------------'
        endif
        do i=0,nconv-1
!         ** Get converged eigenpairs: (in this example they are always real)
          call NEPGetEigenpair(nep,i,lambda,PETSC_NULL_SCALAR,x,PETSC_NULL_VEC,ierr);CHKERRA(ierr)

!         ** Compute residual norm and error
          call NEPComputeError(nep,i,NEP_ERROR_RELATIVE,norm,ierr);CHKERRA(ierr)
          if (rank .eq. 0) then
            write(*,'(1P,E15.4,E18.4)') PetscRealPart(lambda), norm
          endif
        enddo
        if (rank .eq. 0) then
          write(*,*)
        endif
      endif

      call NEPDestroy(nep,ierr);CHKERRA(ierr)
      call MatDestroy(F,ierr);CHKERRA(ierr)
      call MatDestroy(J,ierr);CHKERRA(ierr)
      call VecDestroy(x,ierr);CHKERRA(ierr)
      call SlepcFinalize(ierr)
      end

! ---------------  Evaluate Function matrix  T(lambda)  ----------------

      subroutine FormFunction(nep,lambda,fun,B,ctx,ierr)
      use UserModule
      implicit none
      NEP            nep
      PetscScalar    lambda, A(3), c, d
      Mat            fun,B
      type(User)     ctx
      PetscReal      h
      PetscInt       i, n, j(3), Istart, Iend, one, two, three
      PetscErrorCode ierr

!     ** Compute Function entries and insert into matrix
      call MatGetSize(fun,n,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
      call MatGetOwnershipRange(fun,Istart,Iend,ierr);CHKERRQ(ierr)
      h = ctx%h
      c = ctx%kappa/(lambda-ctx%kappa)
      d = n
      one = 1
      two = 2
      three = 3

!     ** Boundary points
      if (Istart .eq. 0) then
        i = 0
        j(1) = 0
        j(2) = 1
        A(1) = 2.0*(d-lambda*h/3.0)
        A(2) = -d-lambda*h/6.0
        call MatSetValues(fun,one,i,two,j,A,INSERT_VALUES,ierr);CHKERRQ(ierr)
        Istart = Istart + 1
      endif

      if (Iend .eq. n) then
        i = n-1
        j(1) = n-2
        j(2) = n-1
        A(1) = -d-lambda*h/6.0
        A(2) = d-lambda*h/3.0+lambda*c
        call MatSetValues(fun,one,i,two,j,A,INSERT_VALUES,ierr);CHKERRQ(ierr)
        Iend = Iend - 1
      endif

!     ** Interior grid points
      do i=Istart,Iend-1
        j(1) = i-1
        j(2) = i
        j(3) = i+1
        A(1) = -d-lambda*h/6.0
        A(2) = 2.0*(d-lambda*h/3.0)
        A(3) = -d-lambda*h/6.0
        call MatSetValues(fun,one,i,three,j,A,INSERT_VALUES,ierr);CHKERRQ(ierr)
      enddo

!     ** Assemble matrix
      call MatAssemblyBegin(fun,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      call MatAssemblyEnd(fun,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      return
      end

! ---------------  Evaluate Jacobian matrix  T'(lambda)  ---------------

      subroutine FormJacobian(nep,lambda,jac,ctx,ierr)
      use UserModule
      implicit none
      NEP            nep
      PetscScalar    lambda, A(3), c
      Mat            jac
      type(User)     ctx
      PetscReal      h
      PetscInt       i, n, j(3), Istart, Iend, one, two, three
      PetscErrorCode ierr

!     ** Compute Jacobian entries and insert into matrix
      call MatGetSize(jac,n,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
      call MatGetOwnershipRange(jac,Istart,Iend,ierr);CHKERRQ(ierr)
      h = ctx%h
      c = ctx%kappa/(lambda-ctx%kappa)
      one = 1
      two = 2
      three = 3

!     ** Boundary points
      if (Istart .eq. 0) then
        i = 0
        j(1) = 0
        j(2) = 1
        A(1) = -2.0*h/3.0
        A(2) = -h/6.0
        call MatSetValues(jac,one,i,two,j,A,INSERT_VALUES,ierr);CHKERRQ(ierr)
        Istart = Istart + 1
      endif

      if (Iend .eq. n) then
        i = n-1
        j(1) = n-2
        j(2) = n-1
        A(1) = -h/6.0
        A(2) = -h/3.0-c*c
        call MatSetValues(jac,one,i,two,j,A,INSERT_VALUES,ierr);CHKERRQ(ierr)
        Iend = Iend - 1
      endif

!     ** Interior grid points
      do i=Istart,Iend-1
        j(1) = i-1
        j(2) = i
        j(3) = i+1
        A(1) = -h/6.0
        A(2) = -2.0*h/3.0
        A(3) = -h/6.0
        call MatSetValues(jac,one,i,three,j,A,INSERT_VALUES,ierr);CHKERRQ(ierr)
      enddo

!     ** Assemble matrix
      call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      return
      end

!/*TEST
!
!   test:
!      suffix: 1
!      args: -nep_target 4
!      filter: sed -e "s/[0-9]\.[0-9]*E-[0-9]*/removed/g" -e "s/ Number of NEP iterations = [ 0-9]*/ Number of NEP iterations = /"
!      requires: !single
!
!TEST*/
