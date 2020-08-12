#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* krylovschur.c */
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
#define epskrylovschursetrestart_ EPSKRYLOVSCHURSETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschursetrestart_ epskrylovschursetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschurgetrestart_ EPSKRYLOVSCHURGETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschurgetrestart_ epskrylovschurgetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschursetlocking_ EPSKRYLOVSCHURSETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschursetlocking_ epskrylovschursetlocking
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschurgetlocking_ EPSKRYLOVSCHURGETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschurgetlocking_ epskrylovschurgetlocking
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschursetpartitions_ EPSKRYLOVSCHURSETPARTITIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschursetpartitions_ epskrylovschursetpartitions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschurgetpartitions_ EPSKRYLOVSCHURGETPARTITIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschurgetpartitions_ epskrylovschurgetpartitions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschursetdetectzeros_ EPSKRYLOVSCHURSETDETECTZEROS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschursetdetectzeros_ epskrylovschursetdetectzeros
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschurgetdetectzeros_ EPSKRYLOVSCHURGETDETECTZEROS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschurgetdetectzeros_ epskrylovschurgetdetectzeros
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschursetdimensions_ EPSKRYLOVSCHURSETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschursetdimensions_ epskrylovschursetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschurgetdimensions_ EPSKRYLOVSCHURGETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschurgetdimensions_ epskrylovschurgetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschurgetsubcommpairs_ EPSKRYLOVSCHURGETSUBCOMMPAIRS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschurgetsubcommpairs_ epskrylovschurgetsubcommpairs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epskrylovschurupdatesubcommmats_ EPSKRYLOVSCHURUPDATESUBCOMMMATS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epskrylovschurupdatesubcommmats_ epskrylovschurupdatesubcommmats
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epskrylovschursetrestart_(EPS eps,PetscReal *keep, int *__ierr){
*__ierr = EPSKrylovSchurSetRestart(
	(EPS)PetscToPointer((eps) ),*keep);
}
SLEPC_EXTERN void  epskrylovschurgetrestart_(EPS eps,PetscReal *keep, int *__ierr){
*__ierr = EPSKrylovSchurGetRestart(
	(EPS)PetscToPointer((eps) ),keep);
}
SLEPC_EXTERN void  epskrylovschursetlocking_(EPS eps,PetscBool *lock, int *__ierr){
*__ierr = EPSKrylovSchurSetLocking(
	(EPS)PetscToPointer((eps) ),*lock);
}
SLEPC_EXTERN void  epskrylovschurgetlocking_(EPS eps,PetscBool *lock, int *__ierr){
*__ierr = EPSKrylovSchurGetLocking(
	(EPS)PetscToPointer((eps) ),lock);
}
SLEPC_EXTERN void  epskrylovschursetpartitions_(EPS eps,PetscInt *npart, int *__ierr){
*__ierr = EPSKrylovSchurSetPartitions(
	(EPS)PetscToPointer((eps) ),*npart);
}
SLEPC_EXTERN void  epskrylovschurgetpartitions_(EPS eps,PetscInt *npart, int *__ierr){
*__ierr = EPSKrylovSchurGetPartitions(
	(EPS)PetscToPointer((eps) ),npart);
}
SLEPC_EXTERN void  epskrylovschursetdetectzeros_(EPS eps,PetscBool *detect, int *__ierr){
*__ierr = EPSKrylovSchurSetDetectZeros(
	(EPS)PetscToPointer((eps) ),*detect);
}
SLEPC_EXTERN void  epskrylovschurgetdetectzeros_(EPS eps,PetscBool *detect, int *__ierr){
*__ierr = EPSKrylovSchurGetDetectZeros(
	(EPS)PetscToPointer((eps) ),detect);
}
SLEPC_EXTERN void  epskrylovschursetdimensions_(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd, int *__ierr){
*__ierr = EPSKrylovSchurSetDimensions(
	(EPS)PetscToPointer((eps) ),*nev,*ncv,*mpd);
}
SLEPC_EXTERN void  epskrylovschurgetdimensions_(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd, int *__ierr){
*__ierr = EPSKrylovSchurGetDimensions(
	(EPS)PetscToPointer((eps) ),nev,ncv,mpd);
}
SLEPC_EXTERN void  epskrylovschurgetsubcommpairs_(EPS eps,PetscInt *i,PetscScalar *eig,Vec v, int *__ierr){
*__ierr = EPSKrylovSchurGetSubcommPairs(
	(EPS)PetscToPointer((eps) ),*i,eig,
	(Vec)PetscToPointer((v) ));
}
SLEPC_EXTERN void  epskrylovschurupdatesubcommmats_(EPS eps,PetscScalar *s,PetscScalar *a,Mat Au,PetscScalar *t,PetscScalar *b,Mat Bu,MatStructure *str,PetscBool *globalup, int *__ierr){
*__ierr = EPSKrylovSchurUpdateSubcommMats(
	(EPS)PetscToPointer((eps) ),*s,*a,
	(Mat)PetscToPointer((Au) ),*t,*b,
	(Mat)PetscToPointer((Bu) ),*str,*globalup);
}
#if defined(__cplusplus)
}
#endif
