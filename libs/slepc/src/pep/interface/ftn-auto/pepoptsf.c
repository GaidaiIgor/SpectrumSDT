#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* pepopts.c */
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
#define pepsetfromoptions_ PEPSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetfromoptions_ pepsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsettolerances_ PEPSETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsettolerances_ pepsettolerances
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetdimensions_ PEPSETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetdimensions_ pepsetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetwhicheigenpairs_ PEPSETWHICHEIGENPAIRS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetwhicheigenpairs_ pepsetwhicheigenpairs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetwhicheigenpairs_ PEPGETWHICHEIGENPAIRS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetwhicheigenpairs_ pepgetwhicheigenpairs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetproblemtype_ PEPSETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetproblemtype_ pepsetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetproblemtype_ PEPGETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetproblemtype_ pepgetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetbasis_ PEPSETBASIS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetbasis_ pepsetbasis
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetbasis_ PEPGETBASIS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetbasis_ pepgetbasis
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsettrackall_ PEPSETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsettrackall_ pepsettrackall
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgettrackall_ PEPGETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgettrackall_ pepgettrackall
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetconvergencetest_ PEPSETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetconvergencetest_ pepsetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetconvergencetest_ PEPGETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetconvergencetest_ pepgetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetstoppingtest_ PEPSETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetstoppingtest_ pepsetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetstoppingtest_ PEPGETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetstoppingtest_ pepgetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetscale_ PEPSETSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetscale_ pepsetscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetextract_ PEPSETEXTRACT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetextract_ pepsetextract
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetextract_ PEPGETEXTRACT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetextract_ pepgetextract
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetrefine_ PEPSETREFINE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetrefine_ pepsetrefine
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  pepsetfromoptions_(PEP pep, int *__ierr){
*__ierr = PEPSetFromOptions(
	(PEP)PetscToPointer((pep) ));
}
SLEPC_EXTERN void  pepsettolerances_(PEP pep,PetscReal *tol,PetscInt *maxits, int *__ierr){
*__ierr = PEPSetTolerances(
	(PEP)PetscToPointer((pep) ),*tol,*maxits);
}
SLEPC_EXTERN void  pepsetdimensions_(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd, int *__ierr){
*__ierr = PEPSetDimensions(
	(PEP)PetscToPointer((pep) ),*nev,*ncv,*mpd);
}
SLEPC_EXTERN void  pepsetwhicheigenpairs_(PEP pep,PEPWhich *which, int *__ierr){
*__ierr = PEPSetWhichEigenpairs(
	(PEP)PetscToPointer((pep) ),*which);
}
SLEPC_EXTERN void  pepgetwhicheigenpairs_(PEP pep,PEPWhich *which, int *__ierr){
*__ierr = PEPGetWhichEigenpairs(
	(PEP)PetscToPointer((pep) ),which);
}
SLEPC_EXTERN void  pepsetproblemtype_(PEP pep,PEPProblemType *type, int *__ierr){
*__ierr = PEPSetProblemType(
	(PEP)PetscToPointer((pep) ),*type);
}
SLEPC_EXTERN void  pepgetproblemtype_(PEP pep,PEPProblemType *type, int *__ierr){
*__ierr = PEPGetProblemType(
	(PEP)PetscToPointer((pep) ),type);
}
SLEPC_EXTERN void  pepsetbasis_(PEP pep,PEPBasis *basis, int *__ierr){
*__ierr = PEPSetBasis(
	(PEP)PetscToPointer((pep) ),*basis);
}
SLEPC_EXTERN void  pepgetbasis_(PEP pep,PEPBasis *basis, int *__ierr){
*__ierr = PEPGetBasis(
	(PEP)PetscToPointer((pep) ),basis);
}
SLEPC_EXTERN void  pepsettrackall_(PEP pep,PetscBool *trackall, int *__ierr){
*__ierr = PEPSetTrackAll(
	(PEP)PetscToPointer((pep) ),*trackall);
}
SLEPC_EXTERN void  pepgettrackall_(PEP pep,PetscBool *trackall, int *__ierr){
*__ierr = PEPGetTrackAll(
	(PEP)PetscToPointer((pep) ),trackall);
}
SLEPC_EXTERN void  pepsetconvergencetest_(PEP pep,PEPConv *conv, int *__ierr){
*__ierr = PEPSetConvergenceTest(
	(PEP)PetscToPointer((pep) ),*conv);
}
SLEPC_EXTERN void  pepgetconvergencetest_(PEP pep,PEPConv *conv, int *__ierr){
*__ierr = PEPGetConvergenceTest(
	(PEP)PetscToPointer((pep) ),conv);
}
SLEPC_EXTERN void  pepsetstoppingtest_(PEP pep,PEPStop *stop, int *__ierr){
*__ierr = PEPSetStoppingTest(
	(PEP)PetscToPointer((pep) ),*stop);
}
SLEPC_EXTERN void  pepgetstoppingtest_(PEP pep,PEPStop *stop, int *__ierr){
*__ierr = PEPGetStoppingTest(
	(PEP)PetscToPointer((pep) ),stop);
}
SLEPC_EXTERN void  pepsetscale_(PEP pep,PEPScale *scale,PetscReal *alpha,Vec Dl,Vec Dr,PetscInt *its,PetscReal *lambda, int *__ierr){
*__ierr = PEPSetScale(
	(PEP)PetscToPointer((pep) ),*scale,*alpha,
	(Vec)PetscToPointer((Dl) ),
	(Vec)PetscToPointer((Dr) ),*its,*lambda);
}
SLEPC_EXTERN void  pepsetextract_(PEP pep,PEPExtract *extract, int *__ierr){
*__ierr = PEPSetExtract(
	(PEP)PetscToPointer((pep) ),*extract);
}
SLEPC_EXTERN void  pepgetextract_(PEP pep,PEPExtract *extract, int *__ierr){
*__ierr = PEPGetExtract(
	(PEP)PetscToPointer((pep) ),extract);
}
SLEPC_EXTERN void  pepsetrefine_(PEP pep,PEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,PEPRefineScheme *scheme, int *__ierr){
*__ierr = PEPSetRefine(
	(PEP)PetscToPointer((pep) ),*refine,*npart,*tol,*its,*scheme);
}
#if defined(__cplusplus)
}
#endif
