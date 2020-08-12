#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* ciss.c */
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
#define epscisssetsizes_ EPSCISSSETSIZES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscisssetsizes_ epscisssetsizes
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscissgetsizes_ EPSCISSGETSIZES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscissgetsizes_ epscissgetsizes
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscisssetthreshold_ EPSCISSSETTHRESHOLD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscisssetthreshold_ epscisssetthreshold
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscissgetthreshold_ EPSCISSGETTHRESHOLD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscissgetthreshold_ epscissgetthreshold
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscisssetrefinement_ EPSCISSSETREFINEMENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscisssetrefinement_ epscisssetrefinement
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscissgetrefinement_ EPSCISSGETREFINEMENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscissgetrefinement_ epscissgetrefinement
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscisssetusest_ EPSCISSSETUSEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscisssetusest_ epscisssetusest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscissgetusest_ EPSCISSGETUSEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscissgetusest_ epscissgetusest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscisssetquadrule_ EPSCISSSETQUADRULE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscisssetquadrule_ epscisssetquadrule
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscissgetquadrule_ EPSCISSGETQUADRULE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscissgetquadrule_ epscissgetquadrule
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscisssetextraction_ EPSCISSSETEXTRACTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscisssetextraction_ epscisssetextraction
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscissgetextraction_ EPSCISSGETEXTRACTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscissgetextraction_ epscissgetextraction
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epscisssetsizes_(EPS eps,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats, int *__ierr){
*__ierr = EPSCISSSetSizes(
	(EPS)PetscToPointer((eps) ),*ip,*bs,*ms,*npart,*bsmax,*realmats);
}
SLEPC_EXTERN void  epscissgetsizes_(EPS eps,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats, int *__ierr){
*__ierr = EPSCISSGetSizes(
	(EPS)PetscToPointer((eps) ),ip,bs,ms,npart,bsmax,realmats);
}
SLEPC_EXTERN void  epscisssetthreshold_(EPS eps,PetscReal *delta,PetscReal *spur, int *__ierr){
*__ierr = EPSCISSSetThreshold(
	(EPS)PetscToPointer((eps) ),*delta,*spur);
}
SLEPC_EXTERN void  epscissgetthreshold_(EPS eps,PetscReal *delta,PetscReal *spur, int *__ierr){
*__ierr = EPSCISSGetThreshold(
	(EPS)PetscToPointer((eps) ),delta,spur);
}
SLEPC_EXTERN void  epscisssetrefinement_(EPS eps,PetscInt *inner,PetscInt *blsize, int *__ierr){
*__ierr = EPSCISSSetRefinement(
	(EPS)PetscToPointer((eps) ),*inner,*blsize);
}
SLEPC_EXTERN void  epscissgetrefinement_(EPS eps,PetscInt *inner,PetscInt *blsize, int *__ierr){
*__ierr = EPSCISSGetRefinement(
	(EPS)PetscToPointer((eps) ),inner,blsize);
}
SLEPC_EXTERN void  epscisssetusest_(EPS eps,PetscBool *usest, int *__ierr){
*__ierr = EPSCISSSetUseST(
	(EPS)PetscToPointer((eps) ),*usest);
}
SLEPC_EXTERN void  epscissgetusest_(EPS eps,PetscBool *usest, int *__ierr){
*__ierr = EPSCISSGetUseST(
	(EPS)PetscToPointer((eps) ),usest);
}
SLEPC_EXTERN void  epscisssetquadrule_(EPS eps,EPSCISSQuadRule *quad, int *__ierr){
*__ierr = EPSCISSSetQuadRule(
	(EPS)PetscToPointer((eps) ),*quad);
}
SLEPC_EXTERN void  epscissgetquadrule_(EPS eps,EPSCISSQuadRule *quad, int *__ierr){
*__ierr = EPSCISSGetQuadRule(
	(EPS)PetscToPointer((eps) ),quad);
}
SLEPC_EXTERN void  epscisssetextraction_(EPS eps,EPSCISSExtraction *extraction, int *__ierr){
*__ierr = EPSCISSSetExtraction(
	(EPS)PetscToPointer((eps) ),*extraction);
}
SLEPC_EXTERN void  epscissgetextraction_(EPS eps,EPSCISSExtraction *extraction, int *__ierr){
*__ierr = EPSCISSGetExtraction(
	(EPS)PetscToPointer((eps) ),extraction);
}
#if defined(__cplusplus)
}
#endif
