#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* rgring.c */
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
#define rgringsetparameters_ RGRINGSETPARAMETERS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgringsetparameters_ rgringsetparameters
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define rgringgetparameters_ RGRINGGETPARAMETERS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define rgringgetparameters_ rgringgetparameters
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  rgringsetparameters_(RG rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscReal *start_ang,PetscReal *end_ang,PetscReal *width, int *__ierr){
*__ierr = RGRingSetParameters(
	(RG)PetscToPointer((rg) ),*center,*radius,*vscale,*start_ang,*end_ang,*width);
}
SLEPC_EXTERN void  rgringgetparameters_(RG rg,PetscScalar *center,PetscReal *radius,PetscReal *vscale,PetscReal *start_ang,PetscReal *end_ang,PetscReal *width, int *__ierr){
*__ierr = RGRingGetParameters(
	(RG)PetscToPointer((rg) ),center,radius,vscale,start_ang,end_ang,width);
}
#if defined(__cplusplus)
}
#endif
