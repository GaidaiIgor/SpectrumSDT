#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* svdprimme.c */
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
#define svdprimmesetblocksize_ SVDPRIMMESETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdprimmesetblocksize_ svdprimmesetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdprimmegetblocksize_ SVDPRIMMEGETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdprimmegetblocksize_ svdprimmegetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdprimmesetmethod_ SVDPRIMMESETMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdprimmesetmethod_ svdprimmesetmethod
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define svdprimmegetmethod_ SVDPRIMMEGETMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define svdprimmegetmethod_ svdprimmegetmethod
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  svdprimmesetblocksize_(SVD svd,PetscInt *bs, int *__ierr){
*__ierr = SVDPRIMMESetBlockSize(
	(SVD)PetscToPointer((svd) ),*bs);
}
SLEPC_EXTERN void  svdprimmegetblocksize_(SVD svd,PetscInt *bs, int *__ierr){
*__ierr = SVDPRIMMEGetBlockSize(
	(SVD)PetscToPointer((svd) ),bs);
}
SLEPC_EXTERN void  svdprimmesetmethod_(SVD svd,SVDPRIMMEMethod *method, int *__ierr){
*__ierr = SVDPRIMMESetMethod(
	(SVD)PetscToPointer((svd) ),*method);
}
SLEPC_EXTERN void  svdprimmegetmethod_(SVD svd,SVDPRIMMEMethod *method, int *__ierr){
*__ierr = SVDPRIMMEGetMethod(
	(SVD)PetscToPointer((svd) ),method);
}
#if defined(__cplusplus)
}
#endif
