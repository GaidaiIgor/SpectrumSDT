#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* gd.c */
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
#define epsgdsetkrylovstart_ EPSGDSETKRYLOVSTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdsetkrylovstart_ epsgdsetkrylovstart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdgetkrylovstart_ EPSGDGETKRYLOVSTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdgetkrylovstart_ epsgdgetkrylovstart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdsetblocksize_ EPSGDSETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdsetblocksize_ epsgdsetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdgetblocksize_ EPSGDGETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdgetblocksize_ epsgdgetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdsetrestart_ EPSGDSETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdsetrestart_ epsgdsetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdgetrestart_ EPSGDGETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdgetrestart_ epsgdgetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdsetinitialsize_ EPSGDSETINITIALSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdsetinitialsize_ epsgdsetinitialsize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdgetinitialsize_ EPSGDGETINITIALSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdgetinitialsize_ epsgdgetinitialsize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdsetborth_ EPSGDSETBORTH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdsetborth_ epsgdsetborth
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdgetborth_ EPSGDGETBORTH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdgetborth_ epsgdgetborth
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdsetdoubleexpansion_ EPSGDSETDOUBLEEXPANSION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdsetdoubleexpansion_ epsgdsetdoubleexpansion
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgdgetdoubleexpansion_ EPSGDGETDOUBLEEXPANSION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgdgetdoubleexpansion_ epsgdgetdoubleexpansion
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epsgdsetkrylovstart_(EPS eps,PetscBool *krylovstart, int *__ierr){
*__ierr = EPSGDSetKrylovStart(
	(EPS)PetscToPointer((eps) ),*krylovstart);
}
SLEPC_EXTERN void  epsgdgetkrylovstart_(EPS eps,PetscBool *krylovstart, int *__ierr){
*__ierr = EPSGDGetKrylovStart(
	(EPS)PetscToPointer((eps) ),krylovstart);
}
SLEPC_EXTERN void  epsgdsetblocksize_(EPS eps,PetscInt *blocksize, int *__ierr){
*__ierr = EPSGDSetBlockSize(
	(EPS)PetscToPointer((eps) ),*blocksize);
}
SLEPC_EXTERN void  epsgdgetblocksize_(EPS eps,PetscInt *blocksize, int *__ierr){
*__ierr = EPSGDGetBlockSize(
	(EPS)PetscToPointer((eps) ),blocksize);
}
SLEPC_EXTERN void  epsgdsetrestart_(EPS eps,PetscInt *minv,PetscInt *plusk, int *__ierr){
*__ierr = EPSGDSetRestart(
	(EPS)PetscToPointer((eps) ),*minv,*plusk);
}
SLEPC_EXTERN void  epsgdgetrestart_(EPS eps,PetscInt *minv,PetscInt *plusk, int *__ierr){
*__ierr = EPSGDGetRestart(
	(EPS)PetscToPointer((eps) ),minv,plusk);
}
SLEPC_EXTERN void  epsgdsetinitialsize_(EPS eps,PetscInt *initialsize, int *__ierr){
*__ierr = EPSGDSetInitialSize(
	(EPS)PetscToPointer((eps) ),*initialsize);
}
SLEPC_EXTERN void  epsgdgetinitialsize_(EPS eps,PetscInt *initialsize, int *__ierr){
*__ierr = EPSGDGetInitialSize(
	(EPS)PetscToPointer((eps) ),initialsize);
}
SLEPC_EXTERN void  epsgdsetborth_(EPS eps,PetscBool *borth, int *__ierr){
*__ierr = EPSGDSetBOrth(
	(EPS)PetscToPointer((eps) ),*borth);
}
SLEPC_EXTERN void  epsgdgetborth_(EPS eps,PetscBool *borth, int *__ierr){
*__ierr = EPSGDGetBOrth(
	(EPS)PetscToPointer((eps) ),borth);
}
SLEPC_EXTERN void  epsgdsetdoubleexpansion_(EPS eps,PetscBool *doubleexp, int *__ierr){
*__ierr = EPSGDSetDoubleExpansion(
	(EPS)PetscToPointer((eps) ),*doubleexp);
}
SLEPC_EXTERN void  epsgdgetdoubleexpansion_(EPS eps,PetscBool *doubleexp, int *__ierr){
*__ierr = EPSGDGetDoubleExpansion(
	(EPS)PetscToPointer((eps) ),doubleexp);
}
#if defined(__cplusplus)
}
#endif
