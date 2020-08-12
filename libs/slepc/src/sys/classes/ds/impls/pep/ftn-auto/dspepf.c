#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* dspep.c */
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

#include "slepcds.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dspepsetdegree_ DSPEPSETDEGREE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dspepsetdegree_ dspepsetdegree
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dspepgetdegree_ DSPEPGETDEGREE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dspepgetdegree_ dspepgetdegree
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  dspepsetdegree_(DS ds,PetscInt *d, int *__ierr){
*__ierr = DSPEPSetDegree(
	(DS)PetscToPointer((ds) ),*d);
}
SLEPC_EXTERN void  dspepgetdegree_(DS ds,PetscInt *d, int *__ierr){
*__ierr = DSPEPGetDegree(
	(DS)PetscToPointer((ds) ),d);
}
#if defined(__cplusplus)
}
#endif
