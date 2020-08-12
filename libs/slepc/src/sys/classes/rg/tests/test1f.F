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
!  Description: Simple example that tests RG interface functions.
!
! ----------------------------------------------------------------------
!
      program main
#include <slepc/finclude/slepcrg.h>
      use slepcrg
      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      RG             rg
      PetscInt       i,n,inside,one
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscReal      re,im
      PetscScalar    ar,ai,cr(10),ci(10)
      PetscScalar    vr(7),vi(7)
      PetscScalar    center
      PetscReal      radius,vscale,a,b,c,d

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      one = 1
      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      call RGCreate(PETSC_COMM_WORLD,rg,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Ellipse
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call RGSetType(rg,RGELLIPSE,ierr)
      center = 1.1
      radius = 2
      vscale = 0.1
      call RGEllipseSetParameters(rg,center,radius,vscale,ierr)
      call RGSetFromOptions(rg,ierr)
      call RGView(rg,PETSC_NULL_VIEWER,ierr)
      re = 0.1
      im = 0.3
#if defined(PETSC_USE_COMPLEX)
      ar = re+im*PETSC_i
      ai = 0.0
#else
      ar = re
      ai = im
#endif
      call RGCheckInside(rg,one,ar,ai,inside,ierr)
      if (rank .eq. 0) then
        if (inside >= 0) then
          write(*,110) re, im, 'inside'
        else
          write(*,110) re, im, 'outside'
        endif
      endif
 110  format ('Point (',F4.1,',',F4.1,') is ',A7,' the region')

      call RGComputeBoundingBox(rg,a,b,c,d,ierr)
      if (rank .eq. 0) then
        write(*,115) a, b, c, d
      endif
 115  format ('Bounding box: [',F4.1,',',F4.1,']x[',F4.1,',',F4.1,']')

      if (rank .eq. 0) then
        write (*,*) 'Contour points:'
      endif
      n = 10
      call RGComputeContour(rg,n,cr,ci,ierr)
      do i=1,n
#if defined(PETSC_USE_COMPLEX)
        re = PetscRealPart(cr(i))
        im = PetscImaginaryPart(cr(i))
#else
        re = cr(i)
        im = ci(i)
#endif
        if (rank .eq. 0) then
          write(*,120) re, im
        endif
      enddo
 120  format ('(',F7.4,',',F7.4,')')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Interval
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      call RGSetType(rg,RGINTERVAL,ierr)
      a = -1
      b =  1
      c = -0.1
      d =  0.1
      call RGIntervalSetEndpoints(rg,a,b,c,d,ierr)
      call RGSetFromOptions(rg,ierr)
      call RGView(rg,PETSC_NULL_VIEWER,ierr)
      re = 0.2
      im = 0
#if defined(PETSC_USE_COMPLEX)
      ar = re+im*PETSC_i
      ai = 0.0
#else
      ar = re
      ai = im
#endif
      call RGCheckInside(rg,one,ar,ai,inside,ierr)
      if (rank .eq. 0) then
        if (inside >= 0) then
          write(*,110) re, im, 'inside'
        else
          write(*,110) re, im, 'outside'
        endif
      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Polygon
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if defined(PETSC_USE_COMPLEX)
      vr(1) = (0.0,2.0)
      vr(2) = (1.0,4.0)
      vr(3) = (2.0,5.0)
      vr(4) = (4.0,3.0)
      vr(5) = (5.0,4.0)
      vr(6) = (6.0,1.0)
      vr(7) = (2.0,0.0)
#else
      vr(1) = 0.0
      vi(1) = 1.0
      vr(2) = 0.0
      vi(2) = -1.0
      vr(3) = 0.6
      vi(3) = -0.8
      vr(4) = 1.0
      vi(4) = -1.0
      vr(5) = 2.0
      vi(5) = 0.0
      vr(6) = 1.0
      vi(6) = 1.0
      vr(7) = 0.6
      vi(7) = 0.8
#endif
      call RGSetType(rg,RGPOLYGON,ierr)
      n = 7
      call RGPolygonSetVertices(rg,n,vr,vi,ierr)
      call RGSetFromOptions(rg,ierr)
      call RGView(rg,PETSC_NULL_VIEWER,ierr)
      re = 5
      im = 0.9
#if defined(PETSC_USE_COMPLEX)
      ar = re+im*PETSC_i
      ai = 0.0
#else
      ar = re
      ai = im
#endif
      call RGCheckInside(rg,one,ar,ai,inside,ierr)
      if (rank .eq. 0) then
        if (inside >= 0) then
          write(*,110) re, im, 'inside'
        else
          write(*,110) re, im, 'outside'
        endif
      endif

!     *** Clean up
      call RGDestroy(rg,ierr)
      call SlepcFinalize(ierr)
      end

!/*TEST
!
!   test:
!      suffix: 1
!      requires: !complex
!
!   test:
!      suffix: 1_complex
!      requires: complex
!
!TEST*/
