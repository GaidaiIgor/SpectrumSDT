#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* slp.c */
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
#define nepslpsetdeflationthreshold_ NEPSLPSETDEFLATIONTHRESHOLD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepslpsetdeflationthreshold_ nepslpsetdeflationthreshold
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepslpgetdeflationthreshold_ NEPSLPGETDEFLATIONTHRESHOLD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepslpgetdeflationthreshold_ nepslpgetdeflationthreshold
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepslpseteps_ NEPSLPSETEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepslpseteps_ nepslpseteps
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepslpgeteps_ NEPSLPGETEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepslpgeteps_ nepslpgeteps
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepslpsetepsleft_ NEPSLPSETEPSLEFT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepslpsetepsleft_ nepslpsetepsleft
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepslpgetepsleft_ NEPSLPGETEPSLEFT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepslpgetepsleft_ nepslpgetepsleft
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepslpsetksp_ NEPSLPSETKSP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepslpsetksp_ nepslpsetksp
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepslpgetksp_ NEPSLPGETKSP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepslpgetksp_ nepslpgetksp
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  nepslpsetdeflationthreshold_(NEP nep,PetscReal *deftol, int *__ierr){
*__ierr = NEPSLPSetDeflationThreshold(
	(NEP)PetscToPointer((nep) ),*deftol);
}
SLEPC_EXTERN void  nepslpgetdeflationthreshold_(NEP nep,PetscReal *deftol, int *__ierr){
*__ierr = NEPSLPGetDeflationThreshold(
	(NEP)PetscToPointer((nep) ),deftol);
}
SLEPC_EXTERN void  nepslpseteps_(NEP nep,EPS eps, int *__ierr){
*__ierr = NEPSLPSetEPS(
	(NEP)PetscToPointer((nep) ),
	(EPS)PetscToPointer((eps) ));
}
SLEPC_EXTERN void  nepslpgeteps_(NEP nep,EPS *eps, int *__ierr){
*__ierr = NEPSLPGetEPS(
	(NEP)PetscToPointer((nep) ),eps);
}
SLEPC_EXTERN void  nepslpsetepsleft_(NEP nep,EPS eps, int *__ierr){
*__ierr = NEPSLPSetEPSLeft(
	(NEP)PetscToPointer((nep) ),
	(EPS)PetscToPointer((eps) ));
}
SLEPC_EXTERN void  nepslpgetepsleft_(NEP nep,EPS *eps, int *__ierr){
*__ierr = NEPSLPGetEPSLeft(
	(NEP)PetscToPointer((nep) ),eps);
}
SLEPC_EXTERN void  nepslpsetksp_(NEP nep,KSP ksp, int *__ierr){
*__ierr = NEPSLPSetKSP(
	(NEP)PetscToPointer((nep) ),
	(KSP)PetscToPointer((ksp) ));
}
SLEPC_EXTERN void  nepslpgetksp_(NEP nep,KSP *ksp, int *__ierr){
*__ierr = NEPSLPGetKSP(
	(NEP)PetscToPointer((nep) ),ksp);
}
#if defined(__cplusplus)
}
#endif
