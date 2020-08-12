#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* svdopts.c */
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

#include "slepcsvd.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsetimplicittranspose_ SVDSETIMPLICITTRANSPOSE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsetimplicittranspose_ svdsetimplicittranspose
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdgetimplicittranspose_ SVDGETIMPLICITTRANSPOSE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdgetimplicittranspose_ svdgetimplicittranspose
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsettolerances_ SVDSETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsettolerances_ svdsettolerances
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsetdimensions_ SVDSETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsetdimensions_ svdsetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsetwhichsingulartriplets_ SVDSETWHICHSINGULARTRIPLETS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsetwhichsingulartriplets_ svdsetwhichsingulartriplets
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdgetwhichsingulartriplets_ SVDGETWHICHSINGULARTRIPLETS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdgetwhichsingulartriplets_ svdgetwhichsingulartriplets
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsetconvergencetest_ SVDSETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsetconvergencetest_ svdsetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdgetconvergencetest_ SVDGETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdgetconvergencetest_ svdgetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsetstoppingtest_ SVDSETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsetstoppingtest_ svdsetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdgetstoppingtest_ SVDGETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdgetstoppingtest_ svdgetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsetfromoptions_ SVDSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsetfromoptions_ svdsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsettrackall_ SVDSETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsettrackall_ svdsettrackall
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdgettrackall_ SVDGETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdgettrackall_ svdgettrackall
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  svdsetimplicittranspose_(SVD svd,PetscBool *impl, int *__ierr){
*__ierr = SVDSetImplicitTranspose(
	(SVD)PetscToPointer((svd) ),*impl);
}
SLEPC_EXTERN void  svdgetimplicittranspose_(SVD svd,PetscBool *impl, int *__ierr){
*__ierr = SVDGetImplicitTranspose(
	(SVD)PetscToPointer((svd) ),impl);
}
SLEPC_EXTERN void  svdsettolerances_(SVD svd,PetscReal *tol,PetscInt *maxits, int *__ierr){
*__ierr = SVDSetTolerances(
	(SVD)PetscToPointer((svd) ),*tol,*maxits);
}
SLEPC_EXTERN void  svdsetdimensions_(SVD svd,PetscInt *nsv,PetscInt *ncv,PetscInt *mpd, int *__ierr){
*__ierr = SVDSetDimensions(
	(SVD)PetscToPointer((svd) ),*nsv,*ncv,*mpd);
}
SLEPC_EXTERN void  svdsetwhichsingulartriplets_(SVD svd,SVDWhich *which, int *__ierr){
*__ierr = SVDSetWhichSingularTriplets(
	(SVD)PetscToPointer((svd) ),*which);
}
SLEPC_EXTERN void  svdgetwhichsingulartriplets_(SVD svd,SVDWhich *which, int *__ierr){
*__ierr = SVDGetWhichSingularTriplets(
	(SVD)PetscToPointer((svd) ),which);
}
SLEPC_EXTERN void  svdsetconvergencetest_(SVD svd,SVDConv *conv, int *__ierr){
*__ierr = SVDSetConvergenceTest(
	(SVD)PetscToPointer((svd) ),*conv);
}
SLEPC_EXTERN void  svdgetconvergencetest_(SVD svd,SVDConv *conv, int *__ierr){
*__ierr = SVDGetConvergenceTest(
	(SVD)PetscToPointer((svd) ),conv);
}
SLEPC_EXTERN void  svdsetstoppingtest_(SVD svd,SVDStop *stop, int *__ierr){
*__ierr = SVDSetStoppingTest(
	(SVD)PetscToPointer((svd) ),*stop);
}
SLEPC_EXTERN void  svdgetstoppingtest_(SVD svd,SVDStop *stop, int *__ierr){
*__ierr = SVDGetStoppingTest(
	(SVD)PetscToPointer((svd) ),stop);
}
SLEPC_EXTERN void  svdsetfromoptions_(SVD svd, int *__ierr){
*__ierr = SVDSetFromOptions(
	(SVD)PetscToPointer((svd) ));
}
SLEPC_EXTERN void  svdsettrackall_(SVD svd,PetscBool *trackall, int *__ierr){
*__ierr = SVDSetTrackAll(
	(SVD)PetscToPointer((svd) ),*trackall);
}
SLEPC_EXTERN void  svdgettrackall_(SVD svd,PetscBool *trackall, int *__ierr){
*__ierr = SVDGetTrackAll(
	(SVD)PetscToPointer((svd) ),trackall);
}
#if defined(__cplusplus)
}
#endif
