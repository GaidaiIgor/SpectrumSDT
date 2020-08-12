#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* vecutil.c */
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

#include "slepcvec.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define vecnormalizecomplex_ VECNORMALIZECOMPLEX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vecnormalizecomplex_ vecnormalizecomplex
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define veccheckorthogonality_ VECCHECKORTHOGONALITY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define veccheckorthogonality_ veccheckorthogonality
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define vecduplicateempty_ VECDUPLICATEEMPTY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define vecduplicateempty_ vecduplicateempty
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  vecnormalizecomplex_(Vec xr,Vec xi,PetscBool *iscomplex,PetscReal *norm, int *__ierr){
*__ierr = VecNormalizeComplex(
	(Vec)PetscToPointer((xr) ),
	(Vec)PetscToPointer((xi) ),*iscomplex,norm);
}
SLEPC_EXTERN void  veccheckorthogonality_(Vec V[],PetscInt *nv,Vec W[],PetscInt *nw,Mat B,PetscViewer viewer,PetscReal *lev, int *__ierr){
*__ierr = VecCheckOrthogonality(V,*nv,W,*nw,
	(Mat)PetscToPointer((B) ),
	(PetscViewer)PetscToPointer((viewer) ),lev);
}
SLEPC_EXTERN void  vecduplicateempty_(Vec v,Vec *newv, int *__ierr){
*__ierr = VecDuplicateEmpty(
	(Vec)PetscToPointer((v) ),newv);
}
#if defined(__cplusplus)
}
#endif
