#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* arnoldi.c */
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
#define epsarnoldisetdelayed_ EPSARNOLDISETDELAYED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsarnoldisetdelayed_ epsarnoldisetdelayed
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsarnoldigetdelayed_ EPSARNOLDIGETDELAYED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsarnoldigetdelayed_ epsarnoldigetdelayed
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epsarnoldisetdelayed_(EPS eps,PetscBool *delayed, int *__ierr){
*__ierr = EPSArnoldiSetDelayed(
	(EPS)PetscToPointer((eps) ),*delayed);
}
SLEPC_EXTERN void  epsarnoldigetdelayed_(EPS eps,PetscBool *delayed, int *__ierr){
*__ierr = EPSArnoldiGetDelayed(
	(EPS)PetscToPointer((eps) ),delayed);
}
#if defined(__cplusplus)
}
#endif
