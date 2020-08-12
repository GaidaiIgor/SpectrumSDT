#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* fncombine.c */
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

#include "slepcfn.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fncombinesetchildren_ FNCOMBINESETCHILDREN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fncombinesetchildren_ fncombinesetchildren
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fncombinegetchildren_ FNCOMBINEGETCHILDREN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fncombinegetchildren_ fncombinegetchildren
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  fncombinesetchildren_(FN fn,FNCombineType *comb,FN f1,FN f2, int *__ierr){
*__ierr = FNCombineSetChildren(
	(FN)PetscToPointer((fn) ),*comb,
	(FN)PetscToPointer((f1) ),
	(FN)PetscToPointer((f2) ));
}
SLEPC_EXTERN void  fncombinegetchildren_(FN fn,FNCombineType *comb,FN *f1,FN *f2, int *__ierr){
*__ierr = FNCombineGetChildren(
	(FN)PetscToPointer((fn) ),comb,f1,f2);
}
#if defined(__cplusplus)
}
#endif
