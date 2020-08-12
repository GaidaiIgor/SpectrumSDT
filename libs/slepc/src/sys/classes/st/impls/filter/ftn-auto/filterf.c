#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* filter.c */
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

#include "slepcst.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stfiltersetinterval_ STFILTERSETINTERVAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stfiltersetinterval_ stfiltersetinterval
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stfiltergetinterval_ STFILTERGETINTERVAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stfiltergetinterval_ stfiltergetinterval
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stfiltersetrange_ STFILTERSETRANGE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stfiltersetrange_ stfiltersetrange
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stfiltergetrange_ STFILTERGETRANGE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stfiltergetrange_ stfiltergetrange
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stfiltersetdegree_ STFILTERSETDEGREE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stfiltersetdegree_ stfiltersetdegree
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stfiltergetdegree_ STFILTERGETDEGREE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stfiltergetdegree_ stfiltergetdegree
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stfiltergetthreshold_ STFILTERGETTHRESHOLD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stfiltergetthreshold_ stfiltergetthreshold
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  stfiltersetinterval_(ST st,PetscReal *inta,PetscReal *intb, int *__ierr){
*__ierr = STFilterSetInterval(
	(ST)PetscToPointer((st) ),*inta,*intb);
}
SLEPC_EXTERN void  stfiltergetinterval_(ST st,PetscReal *inta,PetscReal *intb, int *__ierr){
*__ierr = STFilterGetInterval(
	(ST)PetscToPointer((st) ),inta,intb);
}
SLEPC_EXTERN void  stfiltersetrange_(ST st,PetscReal *left,PetscReal *right, int *__ierr){
*__ierr = STFilterSetRange(
	(ST)PetscToPointer((st) ),*left,*right);
}
SLEPC_EXTERN void  stfiltergetrange_(ST st,PetscReal *left,PetscReal *right, int *__ierr){
*__ierr = STFilterGetRange(
	(ST)PetscToPointer((st) ),left,right);
}
SLEPC_EXTERN void  stfiltersetdegree_(ST st,PetscInt *deg, int *__ierr){
*__ierr = STFilterSetDegree(
	(ST)PetscToPointer((st) ),*deg);
}
SLEPC_EXTERN void  stfiltergetdegree_(ST st,PetscInt *deg, int *__ierr){
*__ierr = STFilterGetDegree(
	(ST)PetscToPointer((st) ),deg);
}
SLEPC_EXTERN void  stfiltergetthreshold_(ST st,PetscReal *gamma, int *__ierr){
*__ierr = STFilterGetThreshold(
	(ST)PetscToPointer((st) ),gamma);
}
#if defined(__cplusplus)
}
#endif
