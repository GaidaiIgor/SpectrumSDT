#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* cross.c */
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

#include "slepcsvd.h"
#include "slepceps.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdcrosssetexplicitmatrix_ SVDCROSSSETEXPLICITMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdcrosssetexplicitmatrix_ svdcrosssetexplicitmatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdcrossgetexplicitmatrix_ SVDCROSSGETEXPLICITMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdcrossgetexplicitmatrix_ svdcrossgetexplicitmatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdcrossseteps_ SVDCROSSSETEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdcrossseteps_ svdcrossseteps
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdcrossgeteps_ SVDCROSSGETEPS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdcrossgeteps_ svdcrossgeteps
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  svdcrosssetexplicitmatrix_(SVD svd,PetscBool *explicitmatrix, int *__ierr){
*__ierr = SVDCrossSetExplicitMatrix(
	(SVD)PetscToPointer((svd) ),*explicitmatrix);
}
SLEPC_EXTERN void  svdcrossgetexplicitmatrix_(SVD svd,PetscBool *explicitmatrix, int *__ierr){
*__ierr = SVDCrossGetExplicitMatrix(
	(SVD)PetscToPointer((svd) ),explicitmatrix);
}
SLEPC_EXTERN void  svdcrossseteps_(SVD svd,EPS eps, int *__ierr){
*__ierr = SVDCrossSetEPS(
	(SVD)PetscToPointer((svd) ),
	(EPS)PetscToPointer((eps) ));
}
SLEPC_EXTERN void  svdcrossgeteps_(SVD svd,EPS *eps, int *__ierr){
*__ierr = SVDCrossGetEPS(
	(SVD)PetscToPointer((svd) ),eps);
}
#if defined(__cplusplus)
}
#endif
