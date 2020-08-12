#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* mfnopts.c */
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

#include "slepcmfn.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfnsetfromoptions_ MFNSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnsetfromoptions_ mfnsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfnsettolerances_ MFNSETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnsettolerances_ mfnsettolerances
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfngetdimensions_ MFNGETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfngetdimensions_ mfngetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfnsetdimensions_ MFNSETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnsetdimensions_ mfnsetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfnseterrorifnotconverged_ MFNSETERRORIFNOTCONVERGED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnseterrorifnotconverged_ mfnseterrorifnotconverged
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfngeterrorifnotconverged_ MFNGETERRORIFNOTCONVERGED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfngeterrorifnotconverged_ mfngeterrorifnotconverged
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  mfnsetfromoptions_(MFN mfn, int *__ierr){
*__ierr = MFNSetFromOptions(
	(MFN)PetscToPointer((mfn) ));
}
SLEPC_EXTERN void  mfnsettolerances_(MFN mfn,PetscReal *tol,PetscInt *maxits, int *__ierr){
*__ierr = MFNSetTolerances(
	(MFN)PetscToPointer((mfn) ),*tol,*maxits);
}
SLEPC_EXTERN void  mfngetdimensions_(MFN mfn,PetscInt *ncv, int *__ierr){
*__ierr = MFNGetDimensions(
	(MFN)PetscToPointer((mfn) ),ncv);
}
SLEPC_EXTERN void  mfnsetdimensions_(MFN mfn,PetscInt *ncv, int *__ierr){
*__ierr = MFNSetDimensions(
	(MFN)PetscToPointer((mfn) ),*ncv);
}
SLEPC_EXTERN void  mfnseterrorifnotconverged_(MFN mfn,PetscBool *flg, int *__ierr){
*__ierr = MFNSetErrorIfNotConverged(
	(MFN)PetscToPointer((mfn) ),*flg);
}
SLEPC_EXTERN void  mfngeterrorifnotconverged_(MFN mfn,PetscBool *flag, int *__ierr){
*__ierr = MFNGetErrorIfNotConverged(
	(MFN)PetscToPointer((mfn) ),flag);
}
#if defined(__cplusplus)
}
#endif
