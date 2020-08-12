#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* bvbiorthog.c */
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

#include "slepcbv.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvbiorthogonalizecolumn_ BVBIORTHOGONALIZECOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvbiorthogonalizecolumn_ bvbiorthogonalizecolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvbiorthonormalizecolumn_ BVBIORTHONORMALIZECOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvbiorthonormalizecolumn_ bvbiorthonormalizecolumn
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  bvbiorthogonalizecolumn_(BV V,BV W,PetscInt *j, int *__ierr){
*__ierr = BVBiorthogonalizeColumn(
	(BV)PetscToPointer((V) ),
	(BV)PetscToPointer((W) ),*j);
}
SLEPC_EXTERN void  bvbiorthonormalizecolumn_(BV V,BV W,PetscInt *j,PetscReal *delta, int *__ierr){
*__ierr = BVBiorthonormalizeColumn(
	(BV)PetscToPointer((V) ),
	(BV)PetscToPointer((W) ),*j,delta);
}
#if defined(__cplusplus)
}
#endif
