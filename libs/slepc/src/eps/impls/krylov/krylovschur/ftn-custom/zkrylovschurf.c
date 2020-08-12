/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepc/private/epsimpl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define epskrylovschursetsubintervals_    EPSKRYLOVSCHURSETSUBINTERVALS
#define epskrylovschurgetsubintervals_    EPSKRYLOVSCHURGETSUBINTERVALS
#define epskrylovschurgetinertias_        EPSKRYLOVSCHURGETINERTIAS
#define epskrylovschurgetsubcomminfo_     EPSKRYLOVSCHURGETSUBCOMMINFO
#define epskrylovschurgetsubcommmats_     EPSKRYLOVSCHURGETSUBCOMMMATS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define epskrylovschursetsubintervals_    epskrylovschursetsubintervals
#define epskrylovschurgetsubintervals_    epskrylovschurgetsubintervals
#define epskrylovschurgetinertias_        epskrylovschurgetinertias
#define epskrylovschurgetsubcomminfo_     epskrylovschurgetsubcomminfo
#define epskrylovschurgetsubcommmats_     epskrylovschurgetsubcommmats
#endif

SLEPC_EXTERN void epskrylovschursetsubintervals_(EPS *eps,PetscReal *subint,PetscErrorCode *ierr)
{
  CHKFORTRANNULLREAL(subint);
  *ierr = EPSKrylovSchurSetSubintervals(*eps,subint);
}

SLEPC_EXTERN void epskrylovschurgetsubintervals_(EPS *eps,PetscReal *subint,PetscErrorCode *ierr)
{
  PetscReal *osubint;
  PetscInt  npart;

  CHKFORTRANNULLREAL(subint);
  *ierr = EPSKrylovSchurGetSubintervals(*eps,&osubint); if (*ierr) return;
  *ierr = EPSKrylovSchurGetPartitions(*eps,&npart); if (*ierr) return;
  *ierr = PetscArraycpy(subint,osubint,npart+1); if (*ierr) return;
  *ierr = PetscFree(osubint);
}

SLEPC_EXTERN void epskrylovschurgetinertias_(EPS *eps,PetscInt *nshift,PetscReal *shifts,PetscInt *inertias,PetscErrorCode *ierr)
{
  PetscReal *oshifts;
  PetscInt  *oinertias;
  PetscInt  n;

  CHKFORTRANNULLREAL(shifts);
  CHKFORTRANNULLINTEGER(inertias);
  *ierr = EPSKrylovSchurGetInertias(*eps,&n,&oshifts,&oinertias); if (*ierr) return;
  if (shifts) { *ierr = PetscArraycpy(shifts,oshifts,n); if (*ierr) return; }
  if (inertias) { *ierr = PetscArraycpy(inertias,oinertias,n); if (*ierr) return; }
  *nshift = n;
  *ierr = PetscFree(oshifts);if (*ierr) return;
  *ierr = PetscFree(oinertias);
}

SLEPC_EXTERN void epskrylovschurgetsubcomminfo_(EPS *eps,PetscInt *k,PetscInt *n,Vec *v,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(v);
  *ierr = EPSKrylovSchurGetSubcommInfo(*eps,k,n,v);
}

SLEPC_EXTERN void epskrylovschurgetsubcommmats_(EPS *eps,Mat *A,Mat *B,PetscErrorCode *ierr)
{
  CHKFORTRANNULLOBJECT(A);
  CHKFORTRANNULLOBJECT(B);
  *ierr = EPSKrylovSchurGetSubcommMats(*eps,A,B);
}

