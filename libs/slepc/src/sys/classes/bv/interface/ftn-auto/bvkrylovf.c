#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* bvkrylov.c */
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
#define bvmatarnoldi_ BVMATARNOLDI
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmatarnoldi_ bvmatarnoldi
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  bvmatarnoldi_(BV V,Mat A,PetscScalar *H,PetscInt *ldh,PetscInt *k,PetscInt *m,PetscReal *beta,PetscBool *breakdown, int *__ierr){
*__ierr = BVMatArnoldi(
	(BV)PetscToPointer((V) ),
	(Mat)PetscToPointer((A) ),H,*ldh,*k,m,beta,breakdown);
}
#if defined(__cplusplus)
}
#endif
