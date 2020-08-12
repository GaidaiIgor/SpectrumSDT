#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* rqcg.c */
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

#include "slepceps.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsrqcgsetreset_ EPSRQCGSETRESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsrqcgsetreset_ epsrqcgsetreset
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsrqcggetreset_ EPSRQCGGETRESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsrqcggetreset_ epsrqcggetreset
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epsrqcgsetreset_(EPS eps,PetscInt *nrest, int *__ierr){
*__ierr = EPSRQCGSetReset(
	(EPS)PetscToPointer((eps) ),*nrest);
}
SLEPC_EXTERN void  epsrqcggetreset_(EPS eps,PetscInt *nrest, int *__ierr){
*__ierr = EPSRQCGGetReset(
	(EPS)PetscToPointer((eps) ),nrest);
}
#if defined(__cplusplus)
}
#endif
