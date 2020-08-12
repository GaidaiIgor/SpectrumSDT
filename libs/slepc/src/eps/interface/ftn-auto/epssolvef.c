#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* epssolve.c */
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
#include "slepcbv.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssolve_ EPSSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssolve_ epssolve
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetiterationnumber_ EPSGETITERATIONNUMBER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetiterationnumber_ epsgetiterationnumber
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetconverged_ EPSGETCONVERGED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetconverged_ epsgetconverged
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetconvergedreason_ EPSGETCONVERGEDREASON
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetconvergedreason_ epsgetconvergedreason
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetinvariantsubspace_ EPSGETINVARIANTSUBSPACE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetinvariantsubspace_ epsgetinvariantsubspace
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgeteigenvector_ EPSGETEIGENVECTOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgeteigenvector_ epsgeteigenvector
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetlefteigenvector_ EPSGETLEFTEIGENVECTOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetlefteigenvector_ epsgetlefteigenvector
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgeterrorestimate_ EPSGETERRORESTIMATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgeterrorestimate_ epsgeterrorestimate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscomputeerror_ EPSCOMPUTEERROR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscomputeerror_ epscomputeerror
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epssolve_(EPS eps, int *__ierr){
*__ierr = EPSSolve(
	(EPS)PetscToPointer((eps) ));
}
SLEPC_EXTERN void  epsgetiterationnumber_(EPS eps,PetscInt *its, int *__ierr){
*__ierr = EPSGetIterationNumber(
	(EPS)PetscToPointer((eps) ),its);
}
SLEPC_EXTERN void  epsgetconverged_(EPS eps,PetscInt *nconv, int *__ierr){
*__ierr = EPSGetConverged(
	(EPS)PetscToPointer((eps) ),nconv);
}
SLEPC_EXTERN void  epsgetconvergedreason_(EPS eps,EPSConvergedReason *reason, int *__ierr){
*__ierr = EPSGetConvergedReason(
	(EPS)PetscToPointer((eps) ),reason);
}
SLEPC_EXTERN void  epsgetinvariantsubspace_(EPS eps,Vec v[], int *__ierr){
*__ierr = EPSGetInvariantSubspace(
	(EPS)PetscToPointer((eps) ),v);
}
SLEPC_EXTERN void  epsgeteigenvector_(EPS eps,PetscInt *i,Vec Vr,Vec Vi, int *__ierr){
*__ierr = EPSGetEigenvector(
	(EPS)PetscToPointer((eps) ),*i,
	(Vec)PetscToPointer((Vr) ),
	(Vec)PetscToPointer((Vi) ));
}
SLEPC_EXTERN void  epsgetlefteigenvector_(EPS eps,PetscInt *i,Vec Wr,Vec Wi, int *__ierr){
*__ierr = EPSGetLeftEigenvector(
	(EPS)PetscToPointer((eps) ),*i,
	(Vec)PetscToPointer((Wr) ),
	(Vec)PetscToPointer((Wi) ));
}
SLEPC_EXTERN void  epsgeterrorestimate_(EPS eps,PetscInt *i,PetscReal *errest, int *__ierr){
*__ierr = EPSGetErrorEstimate(
	(EPS)PetscToPointer((eps) ),*i,errest);
}
SLEPC_EXTERN void  epscomputeerror_(EPS eps,PetscInt *i,EPSErrorType *type,PetscReal *error, int *__ierr){
*__ierr = EPSComputeError(
	(EPS)PetscToPointer((eps) ),*i,*type,error);
}
#if defined(__cplusplus)
}
#endif
