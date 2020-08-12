#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* slepcsc.c */
/* Fortran interface file */

/*
* This file was generated automatically by bfort from the C source
* file.  
 */

#ifdef PETSC_USE_POINTER_CONVERSION
#if defined(__cplusplus)
extern "C" { 
#endif 
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *);
extern void PetscRmPointer(void*);
#if defined(__cplusplus)
} 
#endif 

#else

#define PetscToPointer(a) (*(PetscFortranAddr *)(a))
#define PetscFromPointer(a) (PetscFortranAddr)(a)
#define PetscRmPointer(a)
#endif

#include "slepcsys.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define slepcsccompare_ SLEPCSCCOMPARE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define slepcsccompare_ slepcsccompare
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define slepcsorteigenvalues_ SLEPCSORTEIGENVALUES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define slepcsorteigenvalues_ slepcsorteigenvalues
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  slepcsccompare_(SlepcSC *sc,PetscScalar *ar,PetscScalar *ai,PetscScalar *br,PetscScalar *bi,PetscInt *res, int *__ierr){
*__ierr = SlepcSCCompare(*sc,*ar,*ai,*br,*bi,res);
}
SLEPC_EXTERN void  slepcsorteigenvalues_(SlepcSC *sc,PetscInt *n,PetscScalar *eigr,PetscScalar *eigi,PetscInt *perm, int *__ierr){
*__ierr = SlepcSortEigenvalues(*sc,*n,eigr,eigi,perm);
}
#if defined(__cplusplus)
}
#endif
