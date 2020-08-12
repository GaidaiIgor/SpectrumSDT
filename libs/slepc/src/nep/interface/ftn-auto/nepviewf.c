#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* nepview.c */
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
#define nepreasonviewfromoptions_ NEPREASONVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepreasonviewfromoptions_ nepreasonviewfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define neperrorviewfromoptions_ NEPERRORVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define neperrorviewfromoptions_ neperrorviewfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepvaluesviewfromoptions_ NEPVALUESVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepvaluesviewfromoptions_ nepvaluesviewfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepvectorsviewfromoptions_ NEPVECTORSVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepvectorsviewfromoptions_ nepvectorsviewfromoptions
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  nepreasonviewfromoptions_(NEP nep, int *__ierr){
*__ierr = NEPReasonViewFromOptions(
	(NEP)PetscToPointer((nep) ));
}
SLEPC_EXTERN void  neperrorviewfromoptions_(NEP nep, int *__ierr){
*__ierr = NEPErrorViewFromOptions(
	(NEP)PetscToPointer((nep) ));
}
SLEPC_EXTERN void  nepvaluesviewfromoptions_(NEP nep, int *__ierr){
*__ierr = NEPValuesViewFromOptions(
	(NEP)PetscToPointer((nep) ));
}
SLEPC_EXTERN void  nepvectorsviewfromoptions_(NEP nep, int *__ierr){
*__ierr = NEPVectorsViewFromOptions(
	(NEP)PetscToPointer((nep) ));
}
#if defined(__cplusplus)
}
#endif
