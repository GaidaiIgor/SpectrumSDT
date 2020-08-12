#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* matutil.c */
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

#include "slepcsys.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define matcreatetile_ MATCREATETILE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define matcreatetile_ matcreatetile
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  matcreatetile_(PetscScalar *a,Mat A,PetscScalar *b,Mat B,PetscScalar *c,Mat C,PetscScalar *d,Mat D,Mat *G, int *__ierr){
*__ierr = MatCreateTile(*a,
	(Mat)PetscToPointer((A) ),*b,
	(Mat)PetscToPointer((B) ),*c,
	(Mat)PetscToPointer((C) ),*d,
	(Mat)PetscToPointer((D) ),G);
}
#if defined(__cplusplus)
}
#endif
