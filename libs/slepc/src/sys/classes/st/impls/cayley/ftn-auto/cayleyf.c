#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* cayley.c */
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

#include "slepcst.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stcayleysetantishift_ STCAYLEYSETANTISHIFT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stcayleysetantishift_ stcayleysetantishift
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stcayleygetantishift_ STCAYLEYGETANTISHIFT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stcayleygetantishift_ stcayleygetantishift
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  stcayleysetantishift_(ST st,PetscScalar *nu, int *__ierr){
*__ierr = STCayleySetAntishift(
	(ST)PetscToPointer((st) ),*nu);
}
SLEPC_EXTERN void  stcayleygetantishift_(ST st,PetscScalar *nu, int *__ierr){
*__ierr = STCayleyGetAntishift(
	(ST)PetscToPointer((st) ),nu);
}
#if defined(__cplusplus)
}
#endif
