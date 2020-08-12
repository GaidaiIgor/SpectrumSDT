#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* rgpolygon.c */
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

#include "slepcrg.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgpolygonsetvertices_ RGPOLYGONSETVERTICES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgpolygonsetvertices_ rgpolygonsetvertices
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  rgpolygonsetvertices_(RG rg,PetscInt *n,PetscScalar vr[],PetscScalar vi[], int *__ierr){
*__ierr = RGPolygonSetVertices(
	(RG)PetscToPointer((rg) ),*n,vr,vi);
}
#if defined(__cplusplus)
}
#endif
