#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* stoar.c */
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
#define pepstoarsetlocking_ PEPSTOARSETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoarsetlocking_ pepstoarsetlocking
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepstoargetlocking_ PEPSTOARGETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoargetlocking_ pepstoargetlocking
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepstoarsetdetectzeros_ PEPSTOARSETDETECTZEROS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoarsetdetectzeros_ pepstoarsetdetectzeros
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepstoargetdetectzeros_ PEPSTOARGETDETECTZEROS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoargetdetectzeros_ pepstoargetdetectzeros
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepstoarsetlinearization_ PEPSTOARSETLINEARIZATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoarsetlinearization_ pepstoarsetlinearization
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepstoargetlinearization_ PEPSTOARGETLINEARIZATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoargetlinearization_ pepstoargetlinearization
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepstoarsetdimensions_ PEPSTOARSETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoarsetdimensions_ pepstoarsetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepstoargetdimensions_ PEPSTOARGETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoargetdimensions_ pepstoargetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepstoarsetcheckeigenvaluetype_ PEPSTOARSETCHECKEIGENVALUETYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoarsetcheckeigenvaluetype_ pepstoarsetcheckeigenvaluetype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepstoargetcheckeigenvaluetype_ PEPSTOARGETCHECKEIGENVALUETYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepstoargetcheckeigenvaluetype_ pepstoargetcheckeigenvaluetype
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  pepstoarsetlocking_(PEP pep,PetscBool *lock, int *__ierr){
*__ierr = PEPSTOARSetLocking(
	(PEP)PetscToPointer((pep) ),*lock);
}
SLEPC_EXTERN void  pepstoargetlocking_(PEP pep,PetscBool *lock, int *__ierr){
*__ierr = PEPSTOARGetLocking(
	(PEP)PetscToPointer((pep) ),lock);
}
SLEPC_EXTERN void  pepstoarsetdetectzeros_(PEP pep,PetscBool *detect, int *__ierr){
*__ierr = PEPSTOARSetDetectZeros(
	(PEP)PetscToPointer((pep) ),*detect);
}
SLEPC_EXTERN void  pepstoargetdetectzeros_(PEP pep,PetscBool *detect, int *__ierr){
*__ierr = PEPSTOARGetDetectZeros(
	(PEP)PetscToPointer((pep) ),detect);
}
SLEPC_EXTERN void  pepstoarsetlinearization_(PEP pep,PetscReal *alpha,PetscReal *beta, int *__ierr){
*__ierr = PEPSTOARSetLinearization(
	(PEP)PetscToPointer((pep) ),*alpha,*beta);
}
SLEPC_EXTERN void  pepstoargetlinearization_(PEP pep,PetscReal *alpha,PetscReal *beta, int *__ierr){
*__ierr = PEPSTOARGetLinearization(
	(PEP)PetscToPointer((pep) ),alpha,beta);
}
SLEPC_EXTERN void  pepstoarsetdimensions_(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd, int *__ierr){
*__ierr = PEPSTOARSetDimensions(
	(PEP)PetscToPointer((pep) ),*nev,*ncv,*mpd);
}
SLEPC_EXTERN void  pepstoargetdimensions_(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd, int *__ierr){
*__ierr = PEPSTOARGetDimensions(
	(PEP)PetscToPointer((pep) ),nev,ncv,mpd);
}
SLEPC_EXTERN void  pepstoarsetcheckeigenvaluetype_(PEP pep,PetscBool *checket, int *__ierr){
*__ierr = PEPSTOARSetCheckEigenvalueType(
	(PEP)PetscToPointer((pep) ),*checket);
}
SLEPC_EXTERN void  pepstoargetcheckeigenvaluetype_(PEP pep,PetscBool *checket, int *__ierr){
*__ierr = PEPSTOARGetCheckEigenvalueType(
	(PEP)PetscToPointer((pep) ),checket);
}
#if defined(__cplusplus)
}
#endif
