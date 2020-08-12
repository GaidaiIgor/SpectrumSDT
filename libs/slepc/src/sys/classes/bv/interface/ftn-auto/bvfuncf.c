#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* bvfunc.c */
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
#define bvdestroy_ BVDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvdestroy_ bvdestroy
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvcreate_ BVCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvcreate_ bvcreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvcreatefrommat_ BVCREATEFROMMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvcreatefrommat_ bvcreatefrommat
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvinsertvec_ BVINSERTVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvinsertvec_ bvinsertvec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvinsertvecs_ BVINSERTVECS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvinsertvecs_ bvinsertvecs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvinsertconstraints_ BVINSERTCONSTRAINTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvinsertconstraints_ bvinsertconstraints
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  bvdestroy_(BV *bv, int *__ierr){
*__ierr = BVDestroy(bv);
}
SLEPC_EXTERN void  bvcreate_(MPI_Fint * comm,BV *newbv, int *__ierr){
*__ierr = BVCreate(
	MPI_Comm_f2c(*(comm)),newbv);
}
SLEPC_EXTERN void  bvcreatefrommat_(Mat A,BV *bv, int *__ierr){
*__ierr = BVCreateFromMat(
	(Mat)PetscToPointer((A) ),bv);
}
SLEPC_EXTERN void  bvinsertvec_(BV V,PetscInt *j,Vec w, int *__ierr){
*__ierr = BVInsertVec(
	(BV)PetscToPointer((V) ),*j,
	(Vec)PetscToPointer((w) ));
}
SLEPC_EXTERN void  bvinsertvecs_(BV V,PetscInt *s,PetscInt *m,Vec *W,PetscBool *orth, int *__ierr){
*__ierr = BVInsertVecs(
	(BV)PetscToPointer((V) ),*s,m,W,*orth);
}
SLEPC_EXTERN void  bvinsertconstraints_(BV V,PetscInt *nc,Vec *C, int *__ierr){
*__ierr = BVInsertConstraints(
	(BV)PetscToPointer((V) ),nc,C);
}
#if defined(__cplusplus)
}
#endif
