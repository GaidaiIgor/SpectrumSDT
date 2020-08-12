#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* rgbasic.c */
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

#include "slepcrg.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgcreate_ RGCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgcreate_ rgcreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgsetfromoptions_ RGSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgsetfromoptions_ rgsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgistrivial_ RGISTRIVIAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgistrivial_ rgistrivial
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgcheckinside_ RGCHECKINSIDE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgcheckinside_ rgcheckinside
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgcomputecontour_ RGCOMPUTECONTOUR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgcomputecontour_ rgcomputecontour
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgcomputeboundingbox_ RGCOMPUTEBOUNDINGBOX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgcomputeboundingbox_ rgcomputeboundingbox
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgsetcomplement_ RGSETCOMPLEMENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgsetcomplement_ rgsetcomplement
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rggetcomplement_ RGGETCOMPLEMENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rggetcomplement_ rggetcomplement
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgsetscale_ RGSETSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgsetscale_ rgsetscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rggetscale_ RGGETSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rggetscale_ rggetscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgpushscale_ RGPUSHSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgpushscale_ rgpushscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgpopscale_ RGPOPSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgpopscale_ rgpopscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgdestroy_ RGDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgdestroy_ rgdestroy
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  rgcreate_(MPI_Fint * comm,RG *newrg, int *__ierr){
*__ierr = RGCreate(
	MPI_Comm_f2c(*(comm)),newrg);
}
SLEPC_EXTERN void  rgsetfromoptions_(RG rg, int *__ierr){
*__ierr = RGSetFromOptions(
	(RG)PetscToPointer((rg) ));
}
SLEPC_EXTERN void  rgistrivial_(RG rg,PetscBool *trivial, int *__ierr){
*__ierr = RGIsTrivial(
	(RG)PetscToPointer((rg) ),trivial);
}
SLEPC_EXTERN void  rgcheckinside_(RG rg,PetscInt *n,PetscScalar *ar,PetscScalar *ai,PetscInt *inside, int *__ierr){
*__ierr = RGCheckInside(
	(RG)PetscToPointer((rg) ),*n,ar,ai,inside);
}
SLEPC_EXTERN void  rgcomputecontour_(RG rg,PetscInt *n,PetscScalar cr[],PetscScalar ci[], int *__ierr){
*__ierr = RGComputeContour(
	(RG)PetscToPointer((rg) ),*n,cr,ci);
}
SLEPC_EXTERN void  rgcomputeboundingbox_(RG rg,PetscReal *a,PetscReal *b,PetscReal *c,PetscReal *d, int *__ierr){
*__ierr = RGComputeBoundingBox(
	(RG)PetscToPointer((rg) ),a,b,c,d);
}
SLEPC_EXTERN void  rgsetcomplement_(RG rg,PetscBool *flg, int *__ierr){
*__ierr = RGSetComplement(
	(RG)PetscToPointer((rg) ),*flg);
}
SLEPC_EXTERN void  rggetcomplement_(RG rg,PetscBool *flg, int *__ierr){
*__ierr = RGGetComplement(
	(RG)PetscToPointer((rg) ),flg);
}
SLEPC_EXTERN void  rgsetscale_(RG rg,PetscReal *sfactor, int *__ierr){
*__ierr = RGSetScale(
	(RG)PetscToPointer((rg) ),*sfactor);
}
SLEPC_EXTERN void  rggetscale_(RG rg,PetscReal *sfactor, int *__ierr){
*__ierr = RGGetScale(
	(RG)PetscToPointer((rg) ),sfactor);
}
SLEPC_EXTERN void  rgpushscale_(RG rg,PetscReal *sfactor, int *__ierr){
*__ierr = RGPushScale(
	(RG)PetscToPointer((rg) ),*sfactor);
}
SLEPC_EXTERN void  rgpopscale_(RG rg, int *__ierr){
*__ierr = RGPopScale(
	(RG)PetscToPointer((rg) ));
}
SLEPC_EXTERN void  rgdestroy_(RG *rg, int *__ierr){
*__ierr = RGDestroy(rg);
}
#if defined(__cplusplus)
}
#endif
