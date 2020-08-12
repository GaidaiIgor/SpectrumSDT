#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* blzpack.c */
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
#define epsblzpacksetblocksize_ EPSBLZPACKSETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsblzpacksetblocksize_ epsblzpacksetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsblzpackgetblocksize_ EPSBLZPACKGETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsblzpackgetblocksize_ epsblzpackgetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsblzpacksetnsteps_ EPSBLZPACKSETNSTEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsblzpacksetnsteps_ epsblzpacksetnsteps
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsblzpackgetnsteps_ EPSBLZPACKGETNSTEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsblzpackgetnsteps_ epsblzpackgetnsteps
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epsblzpacksetblocksize_(EPS eps,PetscInt *bs, int *__ierr){
*__ierr = EPSBlzpackSetBlockSize(
	(EPS)PetscToPointer((eps) ),*bs);
}
SLEPC_EXTERN void  epsblzpackgetblocksize_(EPS eps,PetscInt *bs, int *__ierr){
*__ierr = EPSBlzpackGetBlockSize(
	(EPS)PetscToPointer((eps) ),bs);
}
SLEPC_EXTERN void  epsblzpacksetnsteps_(EPS eps,PetscInt *nsteps, int *__ierr){
*__ierr = EPSBlzpackSetNSteps(
	(EPS)PetscToPointer((eps) ),*nsteps);
}
SLEPC_EXTERN void  epsblzpackgetnsteps_(EPS eps,PetscInt *nsteps, int *__ierr){
*__ierr = EPSBlzpackGetNSteps(
	(EPS)PetscToPointer((eps) ),nsteps);
}
#if defined(__cplusplus)
}
#endif
