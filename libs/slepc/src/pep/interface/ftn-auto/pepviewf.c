#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* pepview.c */
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
#define pepreasonviewfromoptions_ PEPREASONVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepreasonviewfromoptions_ pepreasonviewfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peperrorviewfromoptions_ PEPERRORVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peperrorviewfromoptions_ peperrorviewfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepvaluesviewfromoptions_ PEPVALUESVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepvaluesviewfromoptions_ pepvaluesviewfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepvectorsviewfromoptions_ PEPVECTORSVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepvectorsviewfromoptions_ pepvectorsviewfromoptions
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  pepreasonviewfromoptions_(PEP pep, int *__ierr){
*__ierr = PEPReasonViewFromOptions(
	(PEP)PetscToPointer((pep) ));
}
SLEPC_EXTERN void  peperrorviewfromoptions_(PEP pep, int *__ierr){
*__ierr = PEPErrorViewFromOptions(
	(PEP)PetscToPointer((pep) ));
}
SLEPC_EXTERN void  pepvaluesviewfromoptions_(PEP pep, int *__ierr){
*__ierr = PEPValuesViewFromOptions(
	(PEP)PetscToPointer((pep) ));
}
SLEPC_EXTERN void  pepvectorsviewfromoptions_(PEP pep, int *__ierr){
*__ierr = PEPVectorsViewFromOptions(
	(PEP)PetscToPointer((pep) ));
}
#if defined(__cplusplus)
}
#endif
