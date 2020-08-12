#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* lmesetup.c */
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

#include "slepclme.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmesetup_ LMESETUP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmesetup_ lmesetup
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmesetcoefficients_ LMESETCOEFFICIENTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmesetcoefficients_ lmesetcoefficients
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmegetcoefficients_ LMEGETCOEFFICIENTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmegetcoefficients_ lmegetcoefficients
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmesetrhs_ LMESETRHS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmesetrhs_ lmesetrhs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmegetrhs_ LMEGETRHS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmegetrhs_ lmegetrhs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmesetsolution_ LMESETSOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmesetsolution_ lmesetsolution
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmegetsolution_ LMEGETSOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmegetsolution_ lmegetsolution
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmeallocatesolution_ LMEALLOCATESOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmeallocatesolution_ lmeallocatesolution
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  lmesetup_(LME lme, int *__ierr){
*__ierr = LMESetUp(
	(LME)PetscToPointer((lme) ));
}
SLEPC_EXTERN void  lmesetcoefficients_(LME lme,Mat A,Mat B,Mat D,Mat E, int *__ierr){
*__ierr = LMESetCoefficients(
	(LME)PetscToPointer((lme) ),
	(Mat)PetscToPointer((A) ),
	(Mat)PetscToPointer((B) ),
	(Mat)PetscToPointer((D) ),
	(Mat)PetscToPointer((E) ));
}
SLEPC_EXTERN void  lmegetcoefficients_(LME lme,Mat *A,Mat *B,Mat *D,Mat *E, int *__ierr){
*__ierr = LMEGetCoefficients(
	(LME)PetscToPointer((lme) ),A,B,D,E);
}
SLEPC_EXTERN void  lmesetrhs_(LME lme,Mat C, int *__ierr){
*__ierr = LMESetRHS(
	(LME)PetscToPointer((lme) ),
	(Mat)PetscToPointer((C) ));
}
SLEPC_EXTERN void  lmegetrhs_(LME lme,Mat *C, int *__ierr){
*__ierr = LMEGetRHS(
	(LME)PetscToPointer((lme) ),C);
}
SLEPC_EXTERN void  lmesetsolution_(LME lme,Mat X, int *__ierr){
*__ierr = LMESetSolution(
	(LME)PetscToPointer((lme) ),
	(Mat)PetscToPointer((X) ));
}
SLEPC_EXTERN void  lmegetsolution_(LME lme,Mat *X, int *__ierr){
*__ierr = LMEGetSolution(
	(LME)PetscToPointer((lme) ),X);
}
SLEPC_EXTERN void  lmeallocatesolution_(LME lme,PetscInt *extra, int *__ierr){
*__ierr = LMEAllocateSolution(
	(LME)PetscToPointer((lme) ),*extra);
}
#if defined(__cplusplus)
}
#endif
