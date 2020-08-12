#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* nepsetup.c */
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
#define nepsetup_ NEPSETUP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetup_ nepsetup
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepallocatesolution_ NEPALLOCATESOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepallocatesolution_ nepallocatesolution
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  nepsetup_(NEP nep, int *__ierr){
*__ierr = NEPSetUp(
	(NEP)PetscToPointer((nep) ));
}
SLEPC_EXTERN void  nepallocatesolution_(NEP nep,PetscInt *extra, int *__ierr){
*__ierr = NEPAllocateSolution(
	(NEP)PetscToPointer((nep) ),*extra);
}
#if defined(__cplusplus)
}
#endif
