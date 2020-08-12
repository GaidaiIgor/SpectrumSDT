#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* pepsetup.c */
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
#define pepsetup_ PEPSETUP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetup_ pepsetup
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetoperators_ PEPSETOPERATORS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetoperators_ pepsetoperators
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetoperators_ PEPGETOPERATORS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetoperators_ pepgetoperators
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetnummatrices_ PEPGETNUMMATRICES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetnummatrices_ pepgetnummatrices
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepallocatesolution_ PEPALLOCATESOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepallocatesolution_ pepallocatesolution
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  pepsetup_(PEP pep, int *__ierr){
*__ierr = PEPSetUp(
	(PEP)PetscToPointer((pep) ));
}
SLEPC_EXTERN void  pepsetoperators_(PEP pep,PetscInt *nmat,Mat A[], int *__ierr){
*__ierr = PEPSetOperators(
	(PEP)PetscToPointer((pep) ),*nmat,A);
}
SLEPC_EXTERN void  pepgetoperators_(PEP pep,PetscInt *k,Mat *A, int *__ierr){
*__ierr = PEPGetOperators(
	(PEP)PetscToPointer((pep) ),*k,A);
}
SLEPC_EXTERN void  pepgetnummatrices_(PEP pep,PetscInt *nmat, int *__ierr){
*__ierr = PEPGetNumMatrices(
	(PEP)PetscToPointer((pep) ),nmat);
}
SLEPC_EXTERN void  pepallocatesolution_(PEP pep,PetscInt *extra, int *__ierr){
*__ierr = PEPAllocateSolution(
	(PEP)PetscToPointer((pep) ),*extra);
}
#if defined(__cplusplus)
}
#endif
