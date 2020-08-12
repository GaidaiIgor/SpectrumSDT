#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* qslice.c */
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

#include "slepcpep.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepcheckdefiniteqep_ PEPCHECKDEFINITEQEP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepcheckdefiniteqep_ pepcheckdefiniteqep
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  pepcheckdefiniteqep_(PEP pep,PetscReal *xi,PetscReal *mu,PetscInt *definite,PetscInt *hyperbolic, int *__ierr){
*__ierr = PEPCheckDefiniteQEP(
	(PEP)PetscToPointer((pep) ),xi,mu,definite,hyperbolic);
}
#if defined(__cplusplus)
}
#endif
