#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* mfnsolve.c */
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

#include "slepcmfn.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfnsolve_ MFNSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnsolve_ mfnsolve
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfnsolvetranspose_ MFNSOLVETRANSPOSE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnsolvetranspose_ mfnsolvetranspose
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfngetiterationnumber_ MFNGETITERATIONNUMBER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfngetiterationnumber_ mfngetiterationnumber
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfngetconvergedreason_ MFNGETCONVERGEDREASON
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfngetconvergedreason_ mfngetconvergedreason
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  mfnsolve_(MFN mfn,Vec b,Vec x, int *__ierr){
*__ierr = MFNSolve(
	(MFN)PetscToPointer((mfn) ),
	(Vec)PetscToPointer((b) ),
	(Vec)PetscToPointer((x) ));
}
SLEPC_EXTERN void  mfnsolvetranspose_(MFN mfn,Vec b,Vec x, int *__ierr){
*__ierr = MFNSolveTranspose(
	(MFN)PetscToPointer((mfn) ),
	(Vec)PetscToPointer((b) ),
	(Vec)PetscToPointer((x) ));
}
SLEPC_EXTERN void  mfngetiterationnumber_(MFN mfn,PetscInt *its, int *__ierr){
*__ierr = MFNGetIterationNumber(
	(MFN)PetscToPointer((mfn) ),its);
}
SLEPC_EXTERN void  mfngetconvergedreason_(MFN mfn,MFNConvergedReason *reason, int *__ierr){
*__ierr = MFNGetConvergedReason(
	(MFN)PetscToPointer((mfn) ),reason);
}
#if defined(__cplusplus)
}
#endif
