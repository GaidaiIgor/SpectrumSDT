#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* bvglobal.c */
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
#define bvdot_ BVDOT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvdot_ bvdot
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvdotvec_ BVDOTVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvdotvec_ bvdotvec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvdotvecbegin_ BVDOTVECBEGIN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvdotvecbegin_ bvdotvecbegin
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvdotvecend_ BVDOTVECEND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvdotvecend_ bvdotvecend
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvdotcolumn_ BVDOTCOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvdotcolumn_ bvdotcolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvdotcolumnbegin_ BVDOTCOLUMNBEGIN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvdotcolumnbegin_ bvdotcolumnbegin
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvdotcolumnend_ BVDOTCOLUMNEND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvdotcolumnend_ bvdotcolumnend
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvnorm_ BVNORM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvnorm_ bvnorm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvnormvec_ BVNORMVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvnormvec_ bvnormvec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvnormvecbegin_ BVNORMVECBEGIN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvnormvecbegin_ bvnormvecbegin
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvnormvecend_ BVNORMVECEND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvnormvecend_ bvnormvecend
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvnormcolumn_ BVNORMCOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvnormcolumn_ bvnormcolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvnormcolumnbegin_ BVNORMCOLUMNBEGIN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvnormcolumnbegin_ bvnormcolumnbegin
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvnormcolumnend_ BVNORMCOLUMNEND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvnormcolumnend_ bvnormcolumnend
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmatproject_ BVMATPROJECT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmatproject_ bvmatproject
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  bvdot_(BV X,BV Y,Mat M, int *__ierr){
*__ierr = BVDot(
	(BV)PetscToPointer((X) ),
	(BV)PetscToPointer((Y) ),
	(Mat)PetscToPointer((M) ));
}
SLEPC_EXTERN void  bvdotvec_(BV X,Vec y,PetscScalar m[], int *__ierr){
*__ierr = BVDotVec(
	(BV)PetscToPointer((X) ),
	(Vec)PetscToPointer((y) ),m);
}
SLEPC_EXTERN void  bvdotvecbegin_(BV X,Vec y,PetscScalar *m, int *__ierr){
*__ierr = BVDotVecBegin(
	(BV)PetscToPointer((X) ),
	(Vec)PetscToPointer((y) ),m);
}
SLEPC_EXTERN void  bvdotvecend_(BV X,Vec y,PetscScalar *m, int *__ierr){
*__ierr = BVDotVecEnd(
	(BV)PetscToPointer((X) ),
	(Vec)PetscToPointer((y) ),m);
}
SLEPC_EXTERN void  bvdotcolumn_(BV X,PetscInt *j,PetscScalar *q, int *__ierr){
*__ierr = BVDotColumn(
	(BV)PetscToPointer((X) ),*j,q);
}
SLEPC_EXTERN void  bvdotcolumnbegin_(BV X,PetscInt *j,PetscScalar *m, int *__ierr){
*__ierr = BVDotColumnBegin(
	(BV)PetscToPointer((X) ),*j,m);
}
SLEPC_EXTERN void  bvdotcolumnend_(BV X,PetscInt *j,PetscScalar *m, int *__ierr){
*__ierr = BVDotColumnEnd(
	(BV)PetscToPointer((X) ),*j,m);
}
SLEPC_EXTERN void  bvnorm_(BV bv,NormType *type,PetscReal *val, int *__ierr){
*__ierr = BVNorm(
	(BV)PetscToPointer((bv) ),*type,val);
}
SLEPC_EXTERN void  bvnormvec_(BV bv,Vec v,NormType *type,PetscReal *val, int *__ierr){
*__ierr = BVNormVec(
	(BV)PetscToPointer((bv) ),
	(Vec)PetscToPointer((v) ),*type,val);
}
SLEPC_EXTERN void  bvnormvecbegin_(BV bv,Vec v,NormType *type,PetscReal *val, int *__ierr){
*__ierr = BVNormVecBegin(
	(BV)PetscToPointer((bv) ),
	(Vec)PetscToPointer((v) ),*type,val);
}
SLEPC_EXTERN void  bvnormvecend_(BV bv,Vec v,NormType *type,PetscReal *val, int *__ierr){
*__ierr = BVNormVecEnd(
	(BV)PetscToPointer((bv) ),
	(Vec)PetscToPointer((v) ),*type,val);
}
SLEPC_EXTERN void  bvnormcolumn_(BV bv,PetscInt *j,NormType *type,PetscReal *val, int *__ierr){
*__ierr = BVNormColumn(
	(BV)PetscToPointer((bv) ),*j,*type,val);
}
SLEPC_EXTERN void  bvnormcolumnbegin_(BV bv,PetscInt *j,NormType *type,PetscReal *val, int *__ierr){
*__ierr = BVNormColumnBegin(
	(BV)PetscToPointer((bv) ),*j,*type,val);
}
SLEPC_EXTERN void  bvnormcolumnend_(BV bv,PetscInt *j,NormType *type,PetscReal *val, int *__ierr){
*__ierr = BVNormColumnEnd(
	(BV)PetscToPointer((bv) ),*j,*type,val);
}
SLEPC_EXTERN void  bvmatproject_(BV X,Mat A,BV Y,Mat M, int *__ierr){
*__ierr = BVMatProject(
	(BV)PetscToPointer((X) ),
	(Mat)PetscToPointer((A) ),
	(BV)PetscToPointer((Y) ),
	(Mat)PetscToPointer((M) ));
}
#if defined(__cplusplus)
}
#endif
