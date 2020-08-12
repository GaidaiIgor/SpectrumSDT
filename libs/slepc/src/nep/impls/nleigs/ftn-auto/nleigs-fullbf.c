#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* nleigs-fullb.c */
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

#include "slepcnep.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigsseteps_ NEPNLEIGSSETEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigsseteps_ nepnleigsseteps
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigsgeteps_ NEPNLEIGSGETEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigsgeteps_ nepnleigsgeteps
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  nepnleigsseteps_(NEP nep,EPS eps, int *__ierr){
*__ierr = NEPNLEIGSSetEPS(
	(NEP)PetscToPointer((nep) ),
	(EPS)PetscToPointer((eps) ));
}
SLEPC_EXTERN void  nepnleigsgeteps_(NEP nep,EPS *eps, int *__ierr){
*__ierr = NEPNLEIGSGetEPS(
	(NEP)PetscToPointer((nep) ),eps);
}
#if defined(__cplusplus)
}
#endif
