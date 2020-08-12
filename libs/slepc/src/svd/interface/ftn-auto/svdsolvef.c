#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* svdsolve.c */
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
#define svdsolve_ SVDSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsolve_ svdsolve
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdgetiterationnumber_ SVDGETITERATIONNUMBER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdgetiterationnumber_ svdgetiterationnumber
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdgetconvergedreason_ SVDGETCONVERGEDREASON
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdgetconvergedreason_ svdgetconvergedreason
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdgetconverged_ SVDGETCONVERGED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdgetconverged_ svdgetconverged
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdcomputeerror_ SVDCOMPUTEERROR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdcomputeerror_ svdcomputeerror
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  svdsolve_(SVD svd, int *__ierr){
*__ierr = SVDSolve(
	(SVD)PetscToPointer((svd) ));
}
SLEPC_EXTERN void  svdgetiterationnumber_(SVD svd,PetscInt *its, int *__ierr){
*__ierr = SVDGetIterationNumber(
	(SVD)PetscToPointer((svd) ),its);
}
SLEPC_EXTERN void  svdgetconvergedreason_(SVD svd,SVDConvergedReason *reason, int *__ierr){
*__ierr = SVDGetConvergedReason(
	(SVD)PetscToPointer((svd) ),reason);
}
SLEPC_EXTERN void  svdgetconverged_(SVD svd,PetscInt *nconv, int *__ierr){
*__ierr = SVDGetConverged(
	(SVD)PetscToPointer((svd) ),nconv);
}
SLEPC_EXTERN void  svdcomputeerror_(SVD svd,PetscInt *i,SVDErrorType *type,PetscReal *error, int *__ierr){
*__ierr = SVDComputeError(
	(SVD)PetscToPointer((svd) ),*i,*type,error);
}
#if defined(__cplusplus)
}
#endif
