/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepc/private/fnimpl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define fnrationalsetnumerator0_    FNRATIONALSETNUMERATOR0
#define fnrationalsetnumerator1_    FNRATIONALSETNUMERATOR1
#define fnrationalsetdenominator0_  FNRATIONALSETDENOMINATOR0
#define fnrationalsetdenominator1_  FNRATIONALSETDENOMINATOR1
#define fnrationalgetnumerator0_    FNRATIONALGETNUMERATOR0
#define fnrationalgetnumerator1_    FNRATIONALGETNUMERATOR1
#define fnrationalgetdenominator0_  FNRATIONALGETDENOMINATOR0
#define fnrationalgetdenominator1_  FNRATIONALGETDENOMINATOR1
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define fnrationalsetnumerator0_    fnrationalsetnumerator0
#define fnrationalsetnumerator1_    fnrationalsetnumerator1
#define fnrationalsetdenominator0_  fnrationalsetdenominator0
#define fnrationalsetdenominator1_  fnrationalsetdenominator1
#define fnrationalgetnumerator0_    fnrationalgetnumerator0
#define fnrationalgetnumerator1_    fnrationalgetnumerator1
#define fnrationalgetdenominator0_  fnrationalgetdenominator0
#define fnrationalgetdenominator1_  fnrationalgetdenominator1
#endif

SLEPC_EXTERN void fnrationalsetnumerator_(FN *fn,PetscInt *np,PetscScalar *pcoeff,PetscErrorCode *ierr)
{
  CHKFORTRANNULLSCALAR(pcoeff);
  *ierr = FNRationalSetNumerator(*fn,*np,pcoeff);
}

SLEPC_EXTERN void fnrationalsetnumerator0_(FN *fn,PetscInt *np,PetscScalar *pcoeff,PetscErrorCode *ierr)
{
  fnrationalsetnumerator_(fn,np,pcoeff,ierr);
}

SLEPC_EXTERN void fnrationalsetnumerator1_(FN *fn,PetscInt *np,PetscScalar *pcoeff,PetscErrorCode *ierr)
{
  fnrationalsetnumerator_(fn,np,pcoeff,ierr);
}

SLEPC_EXTERN void fnrationalsetdenominator_(FN *fn,PetscInt *nq,PetscScalar *qcoeff,PetscErrorCode *ierr)
{
  CHKFORTRANNULLSCALAR(qcoeff);
  *ierr = FNRationalSetDenominator(*fn,*nq,qcoeff);
}

SLEPC_EXTERN void fnrationalsetdenominator0_(FN *fn,PetscInt *nq,PetscScalar *qcoeff,PetscErrorCode *ierr)
{
  fnrationalsetdenominator_(fn,nq,qcoeff,ierr);
}

SLEPC_EXTERN void fnrationalsetdenominator1_(FN *fn,PetscInt *nq,PetscScalar *qcoeff,PetscErrorCode *ierr)
{
  fnrationalsetdenominator_(fn,nq,qcoeff,ierr);
}

SLEPC_EXTERN void fnrationalgetnumerator_(FN *fn,PetscInt *np,PetscScalar *pcoeff,PetscErrorCode *ierr)
{
  PetscScalar *ocoeff;
  PetscInt    n;

  CHKFORTRANNULLINTEGER(np);
  CHKFORTRANNULLSCALAR(pcoeff);
  *ierr = FNRationalGetNumerator(*fn,&n,&ocoeff); if (*ierr) return;
  if (pcoeff && ocoeff) { *ierr = PetscArraycpy(pcoeff,ocoeff,n); if (*ierr) return; }
  if (np) *np = n;
  *ierr = PetscFree(ocoeff);
}

SLEPC_EXTERN void fnrationalgetnumerator0_(FN *fn,PetscInt *np,PetscScalar *qcoeff,PetscErrorCode *ierr)
{
  fnrationalgetnumerator_(fn,np,qcoeff,ierr);
}

SLEPC_EXTERN void fnrationalgetnumerator1_(FN *fn,PetscInt *np,PetscScalar *qcoeff,PetscErrorCode *ierr)
{
  fnrationalgetnumerator_(fn,np,qcoeff,ierr);
}

SLEPC_EXTERN void fnrationalgetdenominator_(FN *fn,PetscInt *nq,PetscScalar *qcoeff,PetscErrorCode *ierr)
{
  PetscScalar *ocoeff;
  PetscInt    n;

  CHKFORTRANNULLINTEGER(nq);
  CHKFORTRANNULLSCALAR(qcoeff);
  *ierr = FNRationalGetDenominator(*fn,&n,&ocoeff); if (*ierr) return;
  if (qcoeff && ocoeff) { *ierr = PetscArraycpy(qcoeff,ocoeff,n); if (*ierr) return; }
  if (nq) *nq = n;
  *ierr = PetscFree(ocoeff);
}

SLEPC_EXTERN void fnrationalgetdenominator0_(FN *fn,PetscInt *nq,PetscScalar *qcoeff,PetscErrorCode *ierr)
{
  fnrationalgetdenominator_(fn,nq,qcoeff,ierr);
}

SLEPC_EXTERN void fnrationalgetdenominator1_(FN *fn,PetscInt *nq,PetscScalar *qcoeff,PetscErrorCode *ierr)
{
  fnrationalgetdenominator_(fn,nq,qcoeff,ierr);
}

