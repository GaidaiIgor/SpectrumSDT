#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* svdsetup.c */
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
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsetoperator_ SVDSETOPERATOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsetoperator_ svdsetoperator
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdgetoperator_ SVDGETOPERATOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdgetoperator_ svdgetoperator
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdsetup_ SVDSETUP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdsetup_ svdsetup
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdallocatesolution_ SVDALLOCATESOLUTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdallocatesolution_ svdallocatesolution
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  svdsetoperator_(SVD svd,Mat mat, int *__ierr){
*__ierr = SVDSetOperator(
	(SVD)PetscToPointer((svd) ),
	(Mat)PetscToPointer((mat) ));
}
SLEPC_EXTERN void  svdgetoperator_(SVD svd,Mat *A, int *__ierr){
*__ierr = SVDGetOperator(
	(SVD)PetscToPointer((svd) ),A);
}
SLEPC_EXTERN void  svdsetup_(SVD svd, int *__ierr){
*__ierr = SVDSetUp(
	(SVD)PetscToPointer((svd) ));
}
SLEPC_EXTERN void  svdallocatesolution_(SVD svd,PetscInt *extra, int *__ierr){
*__ierr = SVDAllocateSolution(
	(SVD)PetscToPointer((svd) ),*extra);
}
#if defined(__cplusplus)
}
#endif
