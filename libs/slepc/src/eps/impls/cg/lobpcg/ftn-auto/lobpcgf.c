#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* lobpcg.c */
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

#include "slepceps.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epslobpcgsetblocksize_ EPSLOBPCGSETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epslobpcgsetblocksize_ epslobpcgsetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epslobpcggetblocksize_ EPSLOBPCGGETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epslobpcggetblocksize_ epslobpcggetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epslobpcgsetrestart_ EPSLOBPCGSETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epslobpcgsetrestart_ epslobpcgsetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epslobpcggetrestart_ EPSLOBPCGGETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epslobpcggetrestart_ epslobpcggetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epslobpcgsetlocking_ EPSLOBPCGSETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epslobpcgsetlocking_ epslobpcgsetlocking
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epslobpcggetlocking_ EPSLOBPCGGETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epslobpcggetlocking_ epslobpcggetlocking
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epslobpcgsetblocksize_(EPS eps,PetscInt *bs, int *__ierr){
*__ierr = EPSLOBPCGSetBlockSize(
	(EPS)PetscToPointer((eps) ),*bs);
}
SLEPC_EXTERN void  epslobpcggetblocksize_(EPS eps,PetscInt *bs, int *__ierr){
*__ierr = EPSLOBPCGGetBlockSize(
	(EPS)PetscToPointer((eps) ),bs);
}
SLEPC_EXTERN void  epslobpcgsetrestart_(EPS eps,PetscReal *restart, int *__ierr){
*__ierr = EPSLOBPCGSetRestart(
	(EPS)PetscToPointer((eps) ),*restart);
}
SLEPC_EXTERN void  epslobpcggetrestart_(EPS eps,PetscReal *restart, int *__ierr){
*__ierr = EPSLOBPCGGetRestart(
	(EPS)PetscToPointer((eps) ),restart);
}
SLEPC_EXTERN void  epslobpcgsetlocking_(EPS eps,PetscBool *lock, int *__ierr){
*__ierr = EPSLOBPCGSetLocking(
	(EPS)PetscToPointer((eps) ),*lock);
}
SLEPC_EXTERN void  epslobpcggetlocking_(EPS eps,PetscBool *lock, int *__ierr){
*__ierr = EPSLOBPCGGetLocking(
	(EPS)PetscToPointer((eps) ),lock);
}
#if defined(__cplusplus)
}
#endif
