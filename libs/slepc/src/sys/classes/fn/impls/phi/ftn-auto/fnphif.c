#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* fnphi.c */
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

#include "slepcfn.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnphisetindex_ FNPHISETINDEX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnphisetindex_ fnphisetindex
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnphigetindex_ FNPHIGETINDEX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnphigetindex_ fnphigetindex
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  fnphisetindex_(FN fn,PetscInt *k, int *__ierr){
*__ierr = FNPhiSetIndex(
	(FN)PetscToPointer((fn) ),*k);
}
SLEPC_EXTERN void  fnphigetindex_(FN fn,PetscInt *k, int *__ierr){
*__ierr = FNPhiGetIndex(
	(FN)PetscToPointer((fn) ),k);
}
#if defined(__cplusplus)
}
#endif
