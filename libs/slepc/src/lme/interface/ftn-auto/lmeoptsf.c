#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* lmeopts.c */
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
#define lmesetfromoptions_ LMESETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmesetfromoptions_ lmesetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmesetproblemtype_ LMESETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmesetproblemtype_ lmesetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmegetproblemtype_ LMEGETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmegetproblemtype_ lmegetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmesettolerances_ LMESETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmesettolerances_ lmesettolerances
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmegetdimensions_ LMEGETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmegetdimensions_ lmegetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmesetdimensions_ LMESETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmesetdimensions_ lmesetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmeseterrorifnotconverged_ LMESETERRORIFNOTCONVERGED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmeseterrorifnotconverged_ lmeseterrorifnotconverged
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmegeterrorifnotconverged_ LMEGETERRORIFNOTCONVERGED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmegeterrorifnotconverged_ lmegeterrorifnotconverged
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  lmesetfromoptions_(LME lme, int *__ierr){
*__ierr = LMESetFromOptions(
	(LME)PetscToPointer((lme) ));
}
SLEPC_EXTERN void  lmesetproblemtype_(LME lme,LMEProblemType *type, int *__ierr){
*__ierr = LMESetProblemType(
	(LME)PetscToPointer((lme) ),*type);
}
SLEPC_EXTERN void  lmegetproblemtype_(LME lme,LMEProblemType *type, int *__ierr){
*__ierr = LMEGetProblemType(
	(LME)PetscToPointer((lme) ),type);
}
SLEPC_EXTERN void  lmesettolerances_(LME lme,PetscReal *tol,PetscInt *maxits, int *__ierr){
*__ierr = LMESetTolerances(
	(LME)PetscToPointer((lme) ),*tol,*maxits);
}
SLEPC_EXTERN void  lmegetdimensions_(LME lme,PetscInt *ncv, int *__ierr){
*__ierr = LMEGetDimensions(
	(LME)PetscToPointer((lme) ),ncv);
}
SLEPC_EXTERN void  lmesetdimensions_(LME lme,PetscInt *ncv, int *__ierr){
*__ierr = LMESetDimensions(
	(LME)PetscToPointer((lme) ),*ncv);
}
SLEPC_EXTERN void  lmeseterrorifnotconverged_(LME lme,PetscBool *flg, int *__ierr){
*__ierr = LMESetErrorIfNotConverged(
	(LME)PetscToPointer((lme) ),*flg);
}
SLEPC_EXTERN void  lmegeterrorifnotconverged_(LME lme,PetscBool *flag, int *__ierr){
*__ierr = LMEGetErrorIfNotConverged(
	(LME)PetscToPointer((lme) ),flag);
}
#if defined(__cplusplus)
}
#endif
