#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* nepopts.c */
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
#define nepsetfromoptions_ NEPSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetfromoptions_ nepsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsettolerances_ NEPSETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsettolerances_ nepsettolerances
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetdimensions_ NEPSETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetdimensions_ nepsetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetwhicheigenpairs_ NEPSETWHICHEIGENPAIRS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetwhicheigenpairs_ nepsetwhicheigenpairs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgetwhicheigenpairs_ NEPGETWHICHEIGENPAIRS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgetwhicheigenpairs_ nepgetwhicheigenpairs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetproblemtype_ NEPSETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetproblemtype_ nepsetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgetproblemtype_ NEPGETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgetproblemtype_ nepgetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsettwosided_ NEPSETTWOSIDED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsettwosided_ nepsettwosided
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgettwosided_ NEPGETTWOSIDED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgettwosided_ nepgettwosided
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetconvergencetest_ NEPSETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetconvergencetest_ nepsetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgetconvergencetest_ NEPGETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgetconvergencetest_ nepgetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetstoppingtest_ NEPSETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetstoppingtest_ nepsetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgetstoppingtest_ NEPGETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgetstoppingtest_ nepgetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsettrackall_ NEPSETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsettrackall_ nepsettrackall
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgettrackall_ NEPGETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgettrackall_ nepgettrackall
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetrefine_ NEPSETREFINE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetrefine_ nepsetrefine
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  nepsetfromoptions_(NEP nep, int *__ierr){
*__ierr = NEPSetFromOptions(
	(NEP)PetscToPointer((nep) ));
}
SLEPC_EXTERN void  nepsettolerances_(NEP nep,PetscReal *tol,PetscInt *maxits, int *__ierr){
*__ierr = NEPSetTolerances(
	(NEP)PetscToPointer((nep) ),*tol,*maxits);
}
SLEPC_EXTERN void  nepsetdimensions_(NEP nep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd, int *__ierr){
*__ierr = NEPSetDimensions(
	(NEP)PetscToPointer((nep) ),*nev,*ncv,*mpd);
}
SLEPC_EXTERN void  nepsetwhicheigenpairs_(NEP nep,NEPWhich *which, int *__ierr){
*__ierr = NEPSetWhichEigenpairs(
	(NEP)PetscToPointer((nep) ),*which);
}
SLEPC_EXTERN void  nepgetwhicheigenpairs_(NEP nep,NEPWhich *which, int *__ierr){
*__ierr = NEPGetWhichEigenpairs(
	(NEP)PetscToPointer((nep) ),which);
}
SLEPC_EXTERN void  nepsetproblemtype_(NEP nep,NEPProblemType *type, int *__ierr){
*__ierr = NEPSetProblemType(
	(NEP)PetscToPointer((nep) ),*type);
}
SLEPC_EXTERN void  nepgetproblemtype_(NEP nep,NEPProblemType *type, int *__ierr){
*__ierr = NEPGetProblemType(
	(NEP)PetscToPointer((nep) ),type);
}
SLEPC_EXTERN void  nepsettwosided_(NEP nep,PetscBool *twosided, int *__ierr){
*__ierr = NEPSetTwoSided(
	(NEP)PetscToPointer((nep) ),*twosided);
}
SLEPC_EXTERN void  nepgettwosided_(NEP nep,PetscBool *twosided, int *__ierr){
*__ierr = NEPGetTwoSided(
	(NEP)PetscToPointer((nep) ),twosided);
}
SLEPC_EXTERN void  nepsetconvergencetest_(NEP nep,NEPConv *conv, int *__ierr){
*__ierr = NEPSetConvergenceTest(
	(NEP)PetscToPointer((nep) ),*conv);
}
SLEPC_EXTERN void  nepgetconvergencetest_(NEP nep,NEPConv *conv, int *__ierr){
*__ierr = NEPGetConvergenceTest(
	(NEP)PetscToPointer((nep) ),conv);
}
SLEPC_EXTERN void  nepsetstoppingtest_(NEP nep,NEPStop *stop, int *__ierr){
*__ierr = NEPSetStoppingTest(
	(NEP)PetscToPointer((nep) ),*stop);
}
SLEPC_EXTERN void  nepgetstoppingtest_(NEP nep,NEPStop *stop, int *__ierr){
*__ierr = NEPGetStoppingTest(
	(NEP)PetscToPointer((nep) ),stop);
}
SLEPC_EXTERN void  nepsettrackall_(NEP nep,PetscBool *trackall, int *__ierr){
*__ierr = NEPSetTrackAll(
	(NEP)PetscToPointer((nep) ),*trackall);
}
SLEPC_EXTERN void  nepgettrackall_(NEP nep,PetscBool *trackall, int *__ierr){
*__ierr = NEPGetTrackAll(
	(NEP)PetscToPointer((nep) ),trackall);
}
SLEPC_EXTERN void  nepsetrefine_(NEP nep,NEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,NEPRefineScheme *scheme, int *__ierr){
*__ierr = NEPSetRefine(
	(NEP)PetscToPointer((nep) ),*refine,*npart,*tol,*its,*scheme);
}
#if defined(__cplusplus)
}
#endif
