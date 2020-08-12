#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* lanczos.c */
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
#define epslanczossetreorthog_ EPSLANCZOSSETREORTHOG
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epslanczossetreorthog_ epslanczossetreorthog
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epslanczosgetreorthog_ EPSLANCZOSGETREORTHOG
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epslanczosgetreorthog_ epslanczosgetreorthog
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epslanczossetreorthog_(EPS eps,EPSLanczosReorthogType *reorthog, int *__ierr){
*__ierr = EPSLanczosSetReorthog(
	(EPS)PetscToPointer((eps) ),*reorthog);
}
SLEPC_EXTERN void  epslanczosgetreorthog_(EPS eps,EPSLanczosReorthogType *reorthog, int *__ierr){
*__ierr = EPSLanczosGetReorthog(
	(EPS)PetscToPointer((eps) ),reorthog);
}
#if defined(__cplusplus)
}
#endif
