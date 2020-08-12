#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* epssetup.c */
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
#define epssetup_ EPSSETUP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetup_ epssetup
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetoperators_ EPSSETOPERATORS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetoperators_ epssetoperators
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetoperators_ EPSGETOPERATORS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetoperators_ epsgetoperators
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsallocatesolution_ EPSALLOCATESOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsallocatesolution_ epsallocatesolution
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epssetup_(EPS eps, int *__ierr){
*__ierr = EPSSetUp(
	(EPS)PetscToPointer((eps) ));
}
SLEPC_EXTERN void  epssetoperators_(EPS eps,Mat A,Mat B, int *__ierr){
*__ierr = EPSSetOperators(
	(EPS)PetscToPointer((eps) ),
	(Mat)PetscToPointer((A) ),
	(Mat)PetscToPointer((B) ));
}
SLEPC_EXTERN void  epsgetoperators_(EPS eps,Mat *A,Mat *B, int *__ierr){
*__ierr = EPSGetOperators(
	(EPS)PetscToPointer((eps) ),A,B);
}
SLEPC_EXTERN void  epsallocatesolution_(EPS eps,PetscInt *extra, int *__ierr){
*__ierr = EPSAllocateSolution(
	(EPS)PetscToPointer((eps) ),*extra);
}
#if defined(__cplusplus)
}
#endif
