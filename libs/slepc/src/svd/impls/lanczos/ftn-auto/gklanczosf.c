#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* gklanczos.c */
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

#include "slepcsvd.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdlanczossetoneside_ SVDLANCZOSSETONESIDE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdlanczossetoneside_ svdlanczossetoneside
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdlanczosgetoneside_ SVDLANCZOSGETONESIDE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdlanczosgetoneside_ svdlanczosgetoneside
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  svdlanczossetoneside_(SVD svd,PetscBool *oneside, int *__ierr){
*__ierr = SVDLanczosSetOneSide(
	(SVD)PetscToPointer((svd) ),*oneside);
}
SLEPC_EXTERN void  svdlanczosgetoneside_(SVD svd,PetscBool *oneside, int *__ierr){
*__ierr = SVDLanczosGetOneSide(
	(SVD)PetscToPointer((svd) ),oneside);
}
#if defined(__cplusplus)
}
#endif
