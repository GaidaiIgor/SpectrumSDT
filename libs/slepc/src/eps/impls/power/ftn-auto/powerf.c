#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* power.c */
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
#define epspowersetshifttype_ EPSPOWERSETSHIFTTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epspowersetshifttype_ epspowersetshifttype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epspowergetshifttype_ EPSPOWERGETSHIFTTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epspowergetshifttype_ epspowergetshifttype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epspowersetnonlinear_ EPSPOWERSETNONLINEAR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epspowersetnonlinear_ epspowersetnonlinear
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epspowergetnonlinear_ EPSPOWERGETNONLINEAR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epspowergetnonlinear_ epspowergetnonlinear
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epspowersetupdate_ EPSPOWERSETUPDATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epspowersetupdate_ epspowersetupdate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epspowergetupdate_ EPSPOWERGETUPDATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epspowergetupdate_ epspowergetupdate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epspowersetsnes_ EPSPOWERSETSNES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epspowersetsnes_ epspowersetsnes
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epspowergetsnes_ EPSPOWERGETSNES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epspowergetsnes_ epspowergetsnes
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epspowersetshifttype_(EPS eps,EPSPowerShiftType *shift, int *__ierr){
*__ierr = EPSPowerSetShiftType(
	(EPS)PetscToPointer((eps) ),*shift);
}
SLEPC_EXTERN void  epspowergetshifttype_(EPS eps,EPSPowerShiftType *shift, int *__ierr){
*__ierr = EPSPowerGetShiftType(
	(EPS)PetscToPointer((eps) ),shift);
}
SLEPC_EXTERN void  epspowersetnonlinear_(EPS eps,PetscBool *nonlinear, int *__ierr){
*__ierr = EPSPowerSetNonlinear(
	(EPS)PetscToPointer((eps) ),*nonlinear);
}
SLEPC_EXTERN void  epspowergetnonlinear_(EPS eps,PetscBool *nonlinear, int *__ierr){
*__ierr = EPSPowerGetNonlinear(
	(EPS)PetscToPointer((eps) ),nonlinear);
}
SLEPC_EXTERN void  epspowersetupdate_(EPS eps,PetscBool *update, int *__ierr){
*__ierr = EPSPowerSetUpdate(
	(EPS)PetscToPointer((eps) ),*update);
}
SLEPC_EXTERN void  epspowergetupdate_(EPS eps,PetscBool *update, int *__ierr){
*__ierr = EPSPowerGetUpdate(
	(EPS)PetscToPointer((eps) ),update);
}
SLEPC_EXTERN void  epspowersetsnes_(EPS eps,SNES snes, int *__ierr){
*__ierr = EPSPowerSetSNES(
	(EPS)PetscToPointer((eps) ),
	(SNES)PetscToPointer((snes) ));
}
SLEPC_EXTERN void  epspowergetsnes_(EPS eps,SNES *snes, int *__ierr){
*__ierr = EPSPowerGetSNES(
	(EPS)PetscToPointer((eps) ),snes);
}
#if defined(__cplusplus)
}
#endif
