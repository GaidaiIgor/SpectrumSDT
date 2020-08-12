#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* linear.c */
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
#define peplinearsetlinearization_ PEPLINEARSETLINEARIZATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peplinearsetlinearization_ peplinearsetlinearization
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peplineargetlinearization_ PEPLINEARGETLINEARIZATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peplineargetlinearization_ peplineargetlinearization
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peplinearsetexplicitmatrix_ PEPLINEARSETEXPLICITMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peplinearsetexplicitmatrix_ peplinearsetexplicitmatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peplineargetexplicitmatrix_ PEPLINEARGETEXPLICITMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peplineargetexplicitmatrix_ peplineargetexplicitmatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peplinearseteps_ PEPLINEARSETEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peplinearseteps_ peplinearseteps
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peplineargeteps_ PEPLINEARGETEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peplineargeteps_ peplineargeteps
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  peplinearsetlinearization_(PEP pep,PetscReal *alpha,PetscReal *beta, int *__ierr){
*__ierr = PEPLinearSetLinearization(
	(PEP)PetscToPointer((pep) ),*alpha,*beta);
}
SLEPC_EXTERN void  peplineargetlinearization_(PEP pep,PetscReal *alpha,PetscReal *beta, int *__ierr){
*__ierr = PEPLinearGetLinearization(
	(PEP)PetscToPointer((pep) ),alpha,beta);
}
SLEPC_EXTERN void  peplinearsetexplicitmatrix_(PEP pep,PetscBool *explicitmatrix, int *__ierr){
*__ierr = PEPLinearSetExplicitMatrix(
	(PEP)PetscToPointer((pep) ),*explicitmatrix);
}
SLEPC_EXTERN void  peplineargetexplicitmatrix_(PEP pep,PetscBool *explicitmatrix, int *__ierr){
*__ierr = PEPLinearGetExplicitMatrix(
	(PEP)PetscToPointer((pep) ),explicitmatrix);
}
SLEPC_EXTERN void  peplinearseteps_(PEP pep,EPS eps, int *__ierr){
*__ierr = PEPLinearSetEPS(
	(PEP)PetscToPointer((pep) ),
	(EPS)PetscToPointer((eps) ));
}
SLEPC_EXTERN void  peplineargeteps_(PEP pep,EPS *eps, int *__ierr){
*__ierr = PEPLinearGetEPS(
	(PEP)PetscToPointer((pep) ),eps);
}
#if defined(__cplusplus)
}
#endif
