#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* dsnep.c */
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

#include "slepcds.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsnepsetfn_ DSNEPSETFN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsnepsetfn_ dsnepsetfn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsnepgetfn_ DSNEPGETFN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsnepgetfn_ dsnepgetfn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsnepgetnumfn_ DSNEPGETNUMFN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsnepgetnumfn_ dsnepgetnumfn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsnepsetcomputematrixfunction_ DSNEPSETCOMPUTEMATRIXFUNCTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsnepsetcomputematrixfunction_ dsnepsetcomputematrixfunction
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  dsnepsetfn_(DS ds,PetscInt *n,FN fn[], int *__ierr){
*__ierr = DSNEPSetFN(
	(DS)PetscToPointer((ds) ),*n,fn);
}
SLEPC_EXTERN void  dsnepgetfn_(DS ds,PetscInt *k,FN *fn, int *__ierr){
*__ierr = DSNEPGetFN(
	(DS)PetscToPointer((ds) ),*k,fn);
}
SLEPC_EXTERN void  dsnepgetnumfn_(DS ds,PetscInt *n, int *__ierr){
*__ierr = DSNEPGetNumFN(
	(DS)PetscToPointer((ds) ),n);
}
SLEPC_EXTERN void  dsnepsetcomputematrixfunction_(DS ds,PetscErrorCode (*fun)(DS,PetscScalar,PetscBool,DSMatType,void*),void*ctx, int *__ierr){
*__ierr = DSNEPSetComputeMatrixFunction(
	(DS)PetscToPointer((ds) ),fun,ctx);
}
#if defined(__cplusplus)
}
#endif
