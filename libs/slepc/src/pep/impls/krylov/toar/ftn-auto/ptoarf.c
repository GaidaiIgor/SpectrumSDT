#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* ptoar.c */
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

#include "slepcpep.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peptoarsetrestart_ PEPTOARSETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peptoarsetrestart_ peptoarsetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peptoargetrestart_ PEPTOARGETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peptoargetrestart_ peptoargetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peptoarsetlocking_ PEPTOARSETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peptoarsetlocking_ peptoarsetlocking
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peptoargetlocking_ PEPTOARGETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peptoargetlocking_ peptoargetlocking
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  peptoarsetrestart_(PEP pep,PetscReal *keep, int *__ierr){
*__ierr = PEPTOARSetRestart(
	(PEP)PetscToPointer((pep) ),*keep);
}
SLEPC_EXTERN void  peptoargetrestart_(PEP pep,PetscReal *keep, int *__ierr){
*__ierr = PEPTOARGetRestart(
	(PEP)PetscToPointer((pep) ),keep);
}
SLEPC_EXTERN void  peptoarsetlocking_(PEP pep,PetscBool *lock, int *__ierr){
*__ierr = PEPTOARSetLocking(
	(PEP)PetscToPointer((pep) ),*lock);
}
SLEPC_EXTERN void  peptoargetlocking_(PEP pep,PetscBool *lock, int *__ierr){
*__ierr = PEPTOARGetLocking(
	(PEP)PetscToPointer((pep) ),lock);
}
#if defined(__cplusplus)
}
#endif
