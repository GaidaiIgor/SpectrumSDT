#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* nciss.c */
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
#define nepcisssetsizes_ NEPCISSSETSIZES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepcisssetsizes_ nepcisssetsizes
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepcissgetsizes_ NEPCISSGETSIZES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepcissgetsizes_ nepcissgetsizes
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepcisssetthreshold_ NEPCISSSETTHRESHOLD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepcisssetthreshold_ nepcisssetthreshold
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepcissgetthreshold_ NEPCISSGETTHRESHOLD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepcissgetthreshold_ nepcissgetthreshold
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepcisssetrefinement_ NEPCISSSETREFINEMENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepcisssetrefinement_ nepcisssetrefinement
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepcissgetrefinement_ NEPCISSGETREFINEMENT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepcissgetrefinement_ nepcissgetrefinement
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  nepcisssetsizes_(NEP nep,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats, int *__ierr){
*__ierr = NEPCISSSetSizes(
	(NEP)PetscToPointer((nep) ),*ip,*bs,*ms,*npart,*bsmax,*realmats);
}
SLEPC_EXTERN void  nepcissgetsizes_(NEP nep,PetscInt *ip,PetscInt *bs,PetscInt *ms,PetscInt *npart,PetscInt *bsmax,PetscBool *realmats, int *__ierr){
*__ierr = NEPCISSGetSizes(
	(NEP)PetscToPointer((nep) ),ip,bs,ms,npart,bsmax,realmats);
}
SLEPC_EXTERN void  nepcisssetthreshold_(NEP nep,PetscReal *delta,PetscReal *spur, int *__ierr){
*__ierr = NEPCISSSetThreshold(
	(NEP)PetscToPointer((nep) ),*delta,*spur);
}
SLEPC_EXTERN void  nepcissgetthreshold_(NEP nep,PetscReal *delta,PetscReal *spur, int *__ierr){
*__ierr = NEPCISSGetThreshold(
	(NEP)PetscToPointer((nep) ),delta,spur);
}
SLEPC_EXTERN void  nepcisssetrefinement_(NEP nep,PetscInt *inner,PetscInt *blsize, int *__ierr){
*__ierr = NEPCISSSetRefinement(
	(NEP)PetscToPointer((nep) ),*inner,*blsize);
}
SLEPC_EXTERN void  nepcissgetrefinement_(NEP nep,PetscInt *inner,PetscInt *blsize, int *__ierr){
*__ierr = NEPCISSGetRefinement(
	(NEP)PetscToPointer((nep) ),inner,blsize);
}
#if defined(__cplusplus)
}
#endif
