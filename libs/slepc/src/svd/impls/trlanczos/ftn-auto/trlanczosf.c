#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* trlanczos.c */
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
#define svdtrlanczossetoneside_ SVDTRLANCZOSSETONESIDE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdtrlanczossetoneside_ svdtrlanczossetoneside
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdtrlanczosgetoneside_ SVDTRLANCZOSGETONESIDE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdtrlanczosgetoneside_ svdtrlanczosgetoneside
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  svdtrlanczossetoneside_(SVD svd,PetscBool *oneside, int *__ierr){
*__ierr = SVDTRLanczosSetOneSide(
	(SVD)PetscToPointer((svd) ),*oneside);
}
SLEPC_EXTERN void  svdtrlanczosgetoneside_(SVD svd,PetscBool *oneside, int *__ierr){
*__ierr = SVDTRLanczosGetOneSide(
	(SVD)PetscToPointer((svd) ),oneside);
}
#if defined(__cplusplus)
}
#endif
