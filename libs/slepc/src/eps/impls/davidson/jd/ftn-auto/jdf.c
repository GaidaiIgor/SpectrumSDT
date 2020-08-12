#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* jd.c */
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
#define epsjdsetkrylovstart_ EPSJDSETKRYLOVSTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdsetkrylovstart_ epsjdsetkrylovstart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdgetkrylovstart_ EPSJDGETKRYLOVSTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdgetkrylovstart_ epsjdgetkrylovstart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdsetblocksize_ EPSJDSETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdsetblocksize_ epsjdsetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdgetblocksize_ EPSJDGETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdgetblocksize_ epsjdgetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdsetrestart_ EPSJDSETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdsetrestart_ epsjdsetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdgetrestart_ EPSJDGETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdgetrestart_ epsjdgetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdsetinitialsize_ EPSJDSETINITIALSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdsetinitialsize_ epsjdsetinitialsize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdgetinitialsize_ EPSJDGETINITIALSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdgetinitialsize_ epsjdgetinitialsize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdsetfix_ EPSJDSETFIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdsetfix_ epsjdsetfix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdgetfix_ EPSJDGETFIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdgetfix_ epsjdgetfix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdsetconstcorrectiontol_ EPSJDSETCONSTCORRECTIONTOL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdsetconstcorrectiontol_ epsjdsetconstcorrectiontol
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdgetconstcorrectiontol_ EPSJDGETCONSTCORRECTIONTOL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdgetconstcorrectiontol_ epsjdgetconstcorrectiontol
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdsetborth_ EPSJDSETBORTH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdsetborth_ epsjdsetborth
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsjdgetborth_ EPSJDGETBORTH
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsjdgetborth_ epsjdgetborth
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epsjdsetkrylovstart_(EPS eps,PetscBool *krylovstart, int *__ierr){
*__ierr = EPSJDSetKrylovStart(
	(EPS)PetscToPointer((eps) ),*krylovstart);
}
SLEPC_EXTERN void  epsjdgetkrylovstart_(EPS eps,PetscBool *krylovstart, int *__ierr){
*__ierr = EPSJDGetKrylovStart(
	(EPS)PetscToPointer((eps) ),krylovstart);
}
SLEPC_EXTERN void  epsjdsetblocksize_(EPS eps,PetscInt *blocksize, int *__ierr){
*__ierr = EPSJDSetBlockSize(
	(EPS)PetscToPointer((eps) ),*blocksize);
}
SLEPC_EXTERN void  epsjdgetblocksize_(EPS eps,PetscInt *blocksize, int *__ierr){
*__ierr = EPSJDGetBlockSize(
	(EPS)PetscToPointer((eps) ),blocksize);
}
SLEPC_EXTERN void  epsjdsetrestart_(EPS eps,PetscInt *minv,PetscInt *plusk, int *__ierr){
*__ierr = EPSJDSetRestart(
	(EPS)PetscToPointer((eps) ),*minv,*plusk);
}
SLEPC_EXTERN void  epsjdgetrestart_(EPS eps,PetscInt *minv,PetscInt *plusk, int *__ierr){
*__ierr = EPSJDGetRestart(
	(EPS)PetscToPointer((eps) ),minv,plusk);
}
SLEPC_EXTERN void  epsjdsetinitialsize_(EPS eps,PetscInt *initialsize, int *__ierr){
*__ierr = EPSJDSetInitialSize(
	(EPS)PetscToPointer((eps) ),*initialsize);
}
SLEPC_EXTERN void  epsjdgetinitialsize_(EPS eps,PetscInt *initialsize, int *__ierr){
*__ierr = EPSJDGetInitialSize(
	(EPS)PetscToPointer((eps) ),initialsize);
}
SLEPC_EXTERN void  epsjdsetfix_(EPS eps,PetscReal *fix, int *__ierr){
*__ierr = EPSJDSetFix(
	(EPS)PetscToPointer((eps) ),*fix);
}
SLEPC_EXTERN void  epsjdgetfix_(EPS eps,PetscReal *fix, int *__ierr){
*__ierr = EPSJDGetFix(
	(EPS)PetscToPointer((eps) ),fix);
}
SLEPC_EXTERN void  epsjdsetconstcorrectiontol_(EPS eps,PetscBool *constant, int *__ierr){
*__ierr = EPSJDSetConstCorrectionTol(
	(EPS)PetscToPointer((eps) ),*constant);
}
SLEPC_EXTERN void  epsjdgetconstcorrectiontol_(EPS eps,PetscBool *constant, int *__ierr){
*__ierr = EPSJDGetConstCorrectionTol(
	(EPS)PetscToPointer((eps) ),constant);
}
SLEPC_EXTERN void  epsjdsetborth_(EPS eps,PetscBool *borth, int *__ierr){
*__ierr = EPSJDSetBOrth(
	(EPS)PetscToPointer((eps) ),*borth);
}
SLEPC_EXTERN void  epsjdgetborth_(EPS eps,PetscBool *borth, int *__ierr){
*__ierr = EPSJDGetBOrth(
	(EPS)PetscToPointer((eps) ),borth);
}
#if defined(__cplusplus)
}
#endif
