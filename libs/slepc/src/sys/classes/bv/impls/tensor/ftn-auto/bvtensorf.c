#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* bvtensor.c */
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

#include "slepcbv.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvtensorbuildfirstcolumn_ BVTENSORBUILDFIRSTCOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvtensorbuildfirstcolumn_ bvtensorbuildfirstcolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvtensorcompress_ BVTENSORCOMPRESS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvtensorcompress_ bvtensorcompress
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvtensorgetdegree_ BVTENSORGETDEGREE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvtensorgetdegree_ bvtensorgetdegree
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvcreatetensor_ BVCREATETENSOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvcreatetensor_ bvcreatetensor
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  bvtensorbuildfirstcolumn_(BV V,PetscInt *k, int *__ierr){
*__ierr = BVTensorBuildFirstColumn(
	(BV)PetscToPointer((V) ),*k);
}
SLEPC_EXTERN void  bvtensorcompress_(BV V,PetscInt *newc, int *__ierr){
*__ierr = BVTensorCompress(
	(BV)PetscToPointer((V) ),*newc);
}
SLEPC_EXTERN void  bvtensorgetdegree_(BV bv,PetscInt *d, int *__ierr){
*__ierr = BVTensorGetDegree(
	(BV)PetscToPointer((bv) ),d);
}
SLEPC_EXTERN void  bvcreatetensor_(BV U,PetscInt *d,BV *V, int *__ierr){
*__ierr = BVCreateTensor(
	(BV)PetscToPointer((U) ),*d,V);
}
#if defined(__cplusplus)
}
#endif
