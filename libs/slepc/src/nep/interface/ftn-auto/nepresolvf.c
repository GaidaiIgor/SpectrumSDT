#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* nepresolv.c */
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

#include "slepcnep.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepapplyresolvent_ NEPAPPLYRESOLVENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepapplyresolvent_ nepapplyresolvent
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  nepapplyresolvent_(NEP nep,RG rg,PetscScalar *omega,Vec v,Vec r, int *__ierr){
*__ierr = NEPApplyResolvent(
	(NEP)PetscToPointer((nep) ),
	(RG)PetscToPointer((rg) ),*omega,
	(Vec)PetscToPointer((v) ),
	(Vec)PetscToPointer((r) ));
}
#if defined(__cplusplus)
}
#endif
