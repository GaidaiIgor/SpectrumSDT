#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* stsles.c */
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
#define stmatmult_ STMATMULT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stmatmult_ stmatmult
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stmatmulttranspose_ STMATMULTTRANSPOSE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stmatmulttranspose_ stmatmulttranspose
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stmatsolve_ STMATSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stmatsolve_ stmatsolve
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stmatsolvetranspose_ STMATSOLVETRANSPOSE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stmatsolvetranspose_ stmatsolvetranspose
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stsetksp_ STSETKSP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stsetksp_ stsetksp
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stgetksp_ STGETKSP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stgetksp_ stgetksp
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stchecknullspace_ STCHECKNULLSPACE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stchecknullspace_ stchecknullspace
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  stmatmult_(ST st,PetscInt *k,Vec x,Vec y, int *__ierr){
*__ierr = STMatMult(
	(ST)PetscToPointer((st) ),*k,
	(Vec)PetscToPointer((x) ),
	(Vec)PetscToPointer((y) ));
}
SLEPC_EXTERN void  stmatmulttranspose_(ST st,PetscInt *k,Vec x,Vec y, int *__ierr){
*__ierr = STMatMultTranspose(
	(ST)PetscToPointer((st) ),*k,
	(Vec)PetscToPointer((x) ),
	(Vec)PetscToPointer((y) ));
}
SLEPC_EXTERN void  stmatsolve_(ST st,Vec b,Vec x, int *__ierr){
*__ierr = STMatSolve(
	(ST)PetscToPointer((st) ),
	(Vec)PetscToPointer((b) ),
	(Vec)PetscToPointer((x) ));
}
SLEPC_EXTERN void  stmatsolvetranspose_(ST st,Vec b,Vec x, int *__ierr){
*__ierr = STMatSolveTranspose(
	(ST)PetscToPointer((st) ),
	(Vec)PetscToPointer((b) ),
	(Vec)PetscToPointer((x) ));
}
SLEPC_EXTERN void  stsetksp_(ST st,KSP ksp, int *__ierr){
*__ierr = STSetKSP(
	(ST)PetscToPointer((st) ),
	(KSP)PetscToPointer((ksp) ));
}
SLEPC_EXTERN void  stgetksp_(ST st,KSP* ksp, int *__ierr){
*__ierr = STGetKSP(
	(ST)PetscToPointer((st) ),ksp);
}
SLEPC_EXTERN void  stchecknullspace_(ST st,BV V, int *__ierr){
*__ierr = STCheckNullSpace(
	(ST)PetscToPointer((st) ),
	(BV)PetscToPointer((V) ));
}
#if defined(__cplusplus)
}
#endif
