#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* qarnoldi.c */
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
#define pepqarnoldisetrestart_ PEPQARNOLDISETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepqarnoldisetrestart_ pepqarnoldisetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepqarnoldigetrestart_ PEPQARNOLDIGETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepqarnoldigetrestart_ pepqarnoldigetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepqarnoldisetlocking_ PEPQARNOLDISETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepqarnoldisetlocking_ pepqarnoldisetlocking
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepqarnoldigetlocking_ PEPQARNOLDIGETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepqarnoldigetlocking_ pepqarnoldigetlocking
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  pepqarnoldisetrestart_(PEP pep,PetscReal *keep, int *__ierr){
*__ierr = PEPQArnoldiSetRestart(
	(PEP)PetscToPointer((pep) ),*keep);
}
SLEPC_EXTERN void  pepqarnoldigetrestart_(PEP pep,PetscReal *keep, int *__ierr){
*__ierr = PEPQArnoldiGetRestart(
	(PEP)PetscToPointer((pep) ),keep);
}
SLEPC_EXTERN void  pepqarnoldisetlocking_(PEP pep,PetscBool *lock, int *__ierr){
*__ierr = PEPQArnoldiSetLocking(
	(PEP)PetscToPointer((pep) ),*lock);
}
SLEPC_EXTERN void  pepqarnoldigetlocking_(PEP pep,PetscBool *lock, int *__ierr){
*__ierr = PEPQArnoldiGetLocking(
	(PEP)PetscToPointer((pep) ),lock);
}
#if defined(__cplusplus)
}
#endif
