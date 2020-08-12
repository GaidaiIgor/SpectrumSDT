#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* epsopts.c */
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
#define epssetfromoptions_ EPSSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetfromoptions_ epssetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssettolerances_ EPSSETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssettolerances_ epssettolerances
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetdimensions_ EPSSETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetdimensions_ epssetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetwhicheigenpairs_ EPSSETWHICHEIGENPAIRS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetwhicheigenpairs_ epssetwhicheigenpairs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetwhicheigenpairs_ EPSGETWHICHEIGENPAIRS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetwhicheigenpairs_ epsgetwhicheigenpairs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetconvergencetest_ EPSSETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetconvergencetest_ epssetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetconvergencetest_ EPSGETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetconvergencetest_ epsgetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetstoppingtest_ EPSSETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetstoppingtest_ epssetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetstoppingtest_ EPSGETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetstoppingtest_ epsgetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetproblemtype_ EPSSETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetproblemtype_ epssetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetproblemtype_ EPSGETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetproblemtype_ epsgetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetextraction_ EPSSETEXTRACTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetextraction_ epssetextraction
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetextraction_ EPSGETEXTRACTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetextraction_ epsgetextraction
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetbalance_ EPSSETBALANCE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetbalance_ epssetbalance
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssettwosided_ EPSSETTWOSIDED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssettwosided_ epssettwosided
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgettwosided_ EPSGETTWOSIDED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgettwosided_ epsgettwosided
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssettrueresidual_ EPSSETTRUERESIDUAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssettrueresidual_ epssettrueresidual
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgettrueresidual_ EPSGETTRUERESIDUAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgettrueresidual_ epsgettrueresidual
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssettrackall_ EPSSETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssettrackall_ epssettrackall
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgettrackall_ EPSGETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgettrackall_ epsgettrackall
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetpurify_ EPSSETPURIFY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetpurify_ epssetpurify
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetpurify_ EPSGETPURIFY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetpurify_ epsgetpurify
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epssetfromoptions_(EPS eps, int *__ierr){
*__ierr = EPSSetFromOptions(
	(EPS)PetscToPointer((eps) ));
}
SLEPC_EXTERN void  epssettolerances_(EPS eps,PetscReal *tol,PetscInt *maxits, int *__ierr){
*__ierr = EPSSetTolerances(
	(EPS)PetscToPointer((eps) ),*tol,*maxits);
}
SLEPC_EXTERN void  epssetdimensions_(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd, int *__ierr){
*__ierr = EPSSetDimensions(
	(EPS)PetscToPointer((eps) ),*nev,*ncv,*mpd);
}
SLEPC_EXTERN void  epssetwhicheigenpairs_(EPS eps,EPSWhich *which, int *__ierr){
*__ierr = EPSSetWhichEigenpairs(
	(EPS)PetscToPointer((eps) ),*which);
}
SLEPC_EXTERN void  epsgetwhicheigenpairs_(EPS eps,EPSWhich *which, int *__ierr){
*__ierr = EPSGetWhichEigenpairs(
	(EPS)PetscToPointer((eps) ),which);
}
SLEPC_EXTERN void  epssetconvergencetest_(EPS eps,EPSConv *conv, int *__ierr){
*__ierr = EPSSetConvergenceTest(
	(EPS)PetscToPointer((eps) ),*conv);
}
SLEPC_EXTERN void  epsgetconvergencetest_(EPS eps,EPSConv *conv, int *__ierr){
*__ierr = EPSGetConvergenceTest(
	(EPS)PetscToPointer((eps) ),conv);
}
SLEPC_EXTERN void  epssetstoppingtest_(EPS eps,EPSStop *stop, int *__ierr){
*__ierr = EPSSetStoppingTest(
	(EPS)PetscToPointer((eps) ),*stop);
}
SLEPC_EXTERN void  epsgetstoppingtest_(EPS eps,EPSStop *stop, int *__ierr){
*__ierr = EPSGetStoppingTest(
	(EPS)PetscToPointer((eps) ),stop);
}
SLEPC_EXTERN void  epssetproblemtype_(EPS eps,EPSProblemType *type, int *__ierr){
*__ierr = EPSSetProblemType(
	(EPS)PetscToPointer((eps) ),*type);
}
SLEPC_EXTERN void  epsgetproblemtype_(EPS eps,EPSProblemType *type, int *__ierr){
*__ierr = EPSGetProblemType(
	(EPS)PetscToPointer((eps) ),type);
}
SLEPC_EXTERN void  epssetextraction_(EPS eps,EPSExtraction *extr, int *__ierr){
*__ierr = EPSSetExtraction(
	(EPS)PetscToPointer((eps) ),*extr);
}
SLEPC_EXTERN void  epsgetextraction_(EPS eps,EPSExtraction *extr, int *__ierr){
*__ierr = EPSGetExtraction(
	(EPS)PetscToPointer((eps) ),extr);
}
SLEPC_EXTERN void  epssetbalance_(EPS eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff, int *__ierr){
*__ierr = EPSSetBalance(
	(EPS)PetscToPointer((eps) ),*bal,*its,*cutoff);
}
SLEPC_EXTERN void  epssettwosided_(EPS eps,PetscBool *twosided, int *__ierr){
*__ierr = EPSSetTwoSided(
	(EPS)PetscToPointer((eps) ),*twosided);
}
SLEPC_EXTERN void  epsgettwosided_(EPS eps,PetscBool *twosided, int *__ierr){
*__ierr = EPSGetTwoSided(
	(EPS)PetscToPointer((eps) ),twosided);
}
SLEPC_EXTERN void  epssettrueresidual_(EPS eps,PetscBool *trueres, int *__ierr){
*__ierr = EPSSetTrueResidual(
	(EPS)PetscToPointer((eps) ),*trueres);
}
SLEPC_EXTERN void  epsgettrueresidual_(EPS eps,PetscBool *trueres, int *__ierr){
*__ierr = EPSGetTrueResidual(
	(EPS)PetscToPointer((eps) ),trueres);
}
SLEPC_EXTERN void  epssettrackall_(EPS eps,PetscBool *trackall, int *__ierr){
*__ierr = EPSSetTrackAll(
	(EPS)PetscToPointer((eps) ),*trackall);
}
SLEPC_EXTERN void  epsgettrackall_(EPS eps,PetscBool *trackall, int *__ierr){
*__ierr = EPSGetTrackAll(
	(EPS)PetscToPointer((eps) ),trackall);
}
SLEPC_EXTERN void  epssetpurify_(EPS eps,PetscBool *purify, int *__ierr){
*__ierr = EPSSetPurify(
	(EPS)PetscToPointer((eps) ),*purify);
}
SLEPC_EXTERN void  epsgetpurify_(EPS eps,PetscBool *purify, int *__ierr){
*__ierr = EPSGetPurify(
	(EPS)PetscToPointer((eps) ),purify);
}
#if defined(__cplusplus)
}
#endif
