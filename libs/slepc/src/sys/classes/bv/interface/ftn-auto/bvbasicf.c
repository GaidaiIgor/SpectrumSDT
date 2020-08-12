#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* bvbasic.c */
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
#define bvsetsizes_ BVSETSIZES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetsizes_ bvsetsizes
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetsizesfromvec_ BVSETSIZESFROMVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetsizesfromvec_ bvsetsizesfromvec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetsizes_ BVGETSIZES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetsizes_ bvgetsizes
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetnumconstraints_ BVSETNUMCONSTRAINTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetnumconstraints_ bvsetnumconstraints
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetnumconstraints_ BVGETNUMCONSTRAINTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetnumconstraints_ bvgetnumconstraints
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvresize_ BVRESIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvresize_ bvresize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetactivecolumns_ BVSETACTIVECOLUMNS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetactivecolumns_ bvsetactivecolumns
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetactivecolumns_ BVGETACTIVECOLUMNS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetactivecolumns_ bvgetactivecolumns
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetmatrix_ BVSETMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetmatrix_ bvsetmatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetmatrix_ BVGETMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetmatrix_ bvgetmatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvapplymatrix_ BVAPPLYMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvapplymatrix_ bvapplymatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvapplymatrixbv_ BVAPPLYMATRIXBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvapplymatrixbv_ bvapplymatrixbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetsignature_ BVSETSIGNATURE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetsignature_ bvsetsignature
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetsignature_ BVGETSIGNATURE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetsignature_ bvgetsignature
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetbuffervec_ BVSETBUFFERVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetbuffervec_ bvsetbuffervec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetbuffervec_ BVGETBUFFERVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetbuffervec_ bvgetbuffervec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetrandomcontext_ BVSETRANDOMCONTEXT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetrandomcontext_ bvsetrandomcontext
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetrandomcontext_ BVGETRANDOMCONTEXT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetrandomcontext_ bvgetrandomcontext
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetfromoptions_ BVSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetfromoptions_ bvsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetorthogonalization_ BVSETORTHOGONALIZATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetorthogonalization_ bvsetorthogonalization
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetorthogonalization_ BVGETORTHOGONALIZATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetorthogonalization_ bvgetorthogonalization
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetmatmultmethod_ BVSETMATMULTMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetmatmultmethod_ bvsetmatmultmethod
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetmatmultmethod_ BVGETMATMULTMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetmatmultmethod_ bvgetmatmultmethod
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetcolumn_ BVGETCOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetcolumn_ bvgetcolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvrestorecolumn_ BVRESTORECOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvrestorecolumn_ bvrestorecolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvcreatevec_ BVCREATEVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvcreatevec_ bvcreatevec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvcreatemat_ BVCREATEMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvcreatemat_ bvcreatemat
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetmat_ BVGETMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetmat_ bvgetmat
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvrestoremat_ BVRESTOREMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvrestoremat_ bvrestoremat
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvduplicate_ BVDUPLICATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvduplicate_ bvduplicate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvduplicateresize_ BVDUPLICATERESIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvduplicateresize_ bvduplicateresize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetcachedbv_ BVGETCACHEDBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetcachedbv_ bvgetcachedbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvcopy_ BVCOPY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvcopy_ bvcopy
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvcopyvec_ BVCOPYVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvcopyvec_ bvcopyvec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvcopycolumn_ BVCOPYCOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvcopycolumn_ bvcopycolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvgetsplit_ BVGETSPLIT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvgetsplit_ bvgetsplit
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvrestoresplit_ BVRESTORESPLIT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvrestoresplit_ bvrestoresplit
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  bvsetsizes_(BV bv,PetscInt *n,PetscInt *N,PetscInt *m, int *__ierr){
*__ierr = BVSetSizes(
	(BV)PetscToPointer((bv) ),*n,*N,*m);
}
SLEPC_EXTERN void  bvsetsizesfromvec_(BV bv,Vec t,PetscInt *m, int *__ierr){
*__ierr = BVSetSizesFromVec(
	(BV)PetscToPointer((bv) ),
	(Vec)PetscToPointer((t) ),*m);
}
SLEPC_EXTERN void  bvgetsizes_(BV bv,PetscInt *n,PetscInt *N,PetscInt *m, int *__ierr){
*__ierr = BVGetSizes(
	(BV)PetscToPointer((bv) ),n,N,m);
}
SLEPC_EXTERN void  bvsetnumconstraints_(BV V,PetscInt *nc, int *__ierr){
*__ierr = BVSetNumConstraints(
	(BV)PetscToPointer((V) ),*nc);
}
SLEPC_EXTERN void  bvgetnumconstraints_(BV bv,PetscInt *nc, int *__ierr){
*__ierr = BVGetNumConstraints(
	(BV)PetscToPointer((bv) ),nc);
}
SLEPC_EXTERN void  bvresize_(BV bv,PetscInt *m,PetscBool *copy, int *__ierr){
*__ierr = BVResize(
	(BV)PetscToPointer((bv) ),*m,*copy);
}
SLEPC_EXTERN void  bvsetactivecolumns_(BV bv,PetscInt *l,PetscInt *k, int *__ierr){
*__ierr = BVSetActiveColumns(
	(BV)PetscToPointer((bv) ),*l,*k);
}
SLEPC_EXTERN void  bvgetactivecolumns_(BV bv,PetscInt *l,PetscInt *k, int *__ierr){
*__ierr = BVGetActiveColumns(
	(BV)PetscToPointer((bv) ),l,k);
}
SLEPC_EXTERN void  bvsetmatrix_(BV bv,Mat B,PetscBool *indef, int *__ierr){
*__ierr = BVSetMatrix(
	(BV)PetscToPointer((bv) ),
	(Mat)PetscToPointer((B) ),*indef);
}
SLEPC_EXTERN void  bvgetmatrix_(BV bv,Mat *B,PetscBool *indef, int *__ierr){
*__ierr = BVGetMatrix(
	(BV)PetscToPointer((bv) ),B,indef);
}
SLEPC_EXTERN void  bvapplymatrix_(BV bv,Vec x,Vec y, int *__ierr){
*__ierr = BVApplyMatrix(
	(BV)PetscToPointer((bv) ),
	(Vec)PetscToPointer((x) ),
	(Vec)PetscToPointer((y) ));
}
SLEPC_EXTERN void  bvapplymatrixbv_(BV X,BV Y, int *__ierr){
*__ierr = BVApplyMatrixBV(
	(BV)PetscToPointer((X) ),
	(BV)PetscToPointer((Y) ));
}
SLEPC_EXTERN void  bvsetsignature_(BV bv,Vec omega, int *__ierr){
*__ierr = BVSetSignature(
	(BV)PetscToPointer((bv) ),
	(Vec)PetscToPointer((omega) ));
}
SLEPC_EXTERN void  bvgetsignature_(BV bv,Vec omega, int *__ierr){
*__ierr = BVGetSignature(
	(BV)PetscToPointer((bv) ),
	(Vec)PetscToPointer((omega) ));
}
SLEPC_EXTERN void  bvsetbuffervec_(BV bv,Vec buffer, int *__ierr){
*__ierr = BVSetBufferVec(
	(BV)PetscToPointer((bv) ),
	(Vec)PetscToPointer((buffer) ));
}
SLEPC_EXTERN void  bvgetbuffervec_(BV bv,Vec *buffer, int *__ierr){
*__ierr = BVGetBufferVec(
	(BV)PetscToPointer((bv) ),buffer);
}
SLEPC_EXTERN void  bvsetrandomcontext_(BV bv,PetscRandom rand, int *__ierr){
*__ierr = BVSetRandomContext(
	(BV)PetscToPointer((bv) ),
	(PetscRandom)PetscToPointer((rand) ));
}
SLEPC_EXTERN void  bvgetrandomcontext_(BV bv,PetscRandom* rand, int *__ierr){
*__ierr = BVGetRandomContext(
	(BV)PetscToPointer((bv) ),rand);
}
SLEPC_EXTERN void  bvsetfromoptions_(BV bv, int *__ierr){
*__ierr = BVSetFromOptions(
	(BV)PetscToPointer((bv) ));
}
SLEPC_EXTERN void  bvsetorthogonalization_(BV bv,BVOrthogType *type,BVOrthogRefineType *refine,PetscReal *eta,BVOrthogBlockType *block, int *__ierr){
*__ierr = BVSetOrthogonalization(
	(BV)PetscToPointer((bv) ),*type,*refine,*eta,*block);
}
SLEPC_EXTERN void  bvgetorthogonalization_(BV bv,BVOrthogType *type,BVOrthogRefineType *refine,PetscReal *eta,BVOrthogBlockType *block, int *__ierr){
*__ierr = BVGetOrthogonalization(
	(BV)PetscToPointer((bv) ),type,refine,eta,block);
}
SLEPC_EXTERN void  bvsetmatmultmethod_(BV bv,BVMatMultType *method, int *__ierr){
*__ierr = BVSetMatMultMethod(
	(BV)PetscToPointer((bv) ),*method);
}
SLEPC_EXTERN void  bvgetmatmultmethod_(BV bv,BVMatMultType *method, int *__ierr){
*__ierr = BVGetMatMultMethod(
	(BV)PetscToPointer((bv) ),method);
}
SLEPC_EXTERN void  bvgetcolumn_(BV bv,PetscInt *j,Vec *v, int *__ierr){
*__ierr = BVGetColumn(
	(BV)PetscToPointer((bv) ),*j,v);
}
SLEPC_EXTERN void  bvrestorecolumn_(BV bv,PetscInt *j,Vec *v, int *__ierr){
*__ierr = BVRestoreColumn(
	(BV)PetscToPointer((bv) ),*j,v);
}
SLEPC_EXTERN void  bvcreatevec_(BV bv,Vec *v, int *__ierr){
*__ierr = BVCreateVec(
	(BV)PetscToPointer((bv) ),v);
}
SLEPC_EXTERN void  bvcreatemat_(BV bv,Mat *A, int *__ierr){
*__ierr = BVCreateMat(
	(BV)PetscToPointer((bv) ),A);
}
SLEPC_EXTERN void  bvgetmat_(BV bv,Mat *A, int *__ierr){
*__ierr = BVGetMat(
	(BV)PetscToPointer((bv) ),A);
}
SLEPC_EXTERN void  bvrestoremat_(BV bv,Mat *A, int *__ierr){
*__ierr = BVRestoreMat(
	(BV)PetscToPointer((bv) ),A);
}
SLEPC_EXTERN void  bvduplicate_(BV V,BV *W, int *__ierr){
*__ierr = BVDuplicate(
	(BV)PetscToPointer((V) ),W);
}
SLEPC_EXTERN void  bvduplicateresize_(BV V,PetscInt *m,BV *W, int *__ierr){
*__ierr = BVDuplicateResize(
	(BV)PetscToPointer((V) ),*m,W);
}
SLEPC_EXTERN void  bvgetcachedbv_(BV bv,BV *cached, int *__ierr){
*__ierr = BVGetCachedBV(
	(BV)PetscToPointer((bv) ),cached);
}
SLEPC_EXTERN void  bvcopy_(BV V,BV W, int *__ierr){
*__ierr = BVCopy(
	(BV)PetscToPointer((V) ),
	(BV)PetscToPointer((W) ));
}
SLEPC_EXTERN void  bvcopyvec_(BV V,PetscInt *j,Vec w, int *__ierr){
*__ierr = BVCopyVec(
	(BV)PetscToPointer((V) ),*j,
	(Vec)PetscToPointer((w) ));
}
SLEPC_EXTERN void  bvcopycolumn_(BV V,PetscInt *j,PetscInt *i, int *__ierr){
*__ierr = BVCopyColumn(
	(BV)PetscToPointer((V) ),*j,*i);
}
SLEPC_EXTERN void  bvgetsplit_(BV bv,BV *L,BV *R, int *__ierr){
*__ierr = BVGetSplit(
	(BV)PetscToPointer((bv) ),L,R);
}
SLEPC_EXTERN void  bvrestoresplit_(BV bv,BV *L,BV *R, int *__ierr){
*__ierr = BVRestoreSplit(
	(BV)PetscToPointer((bv) ),L,R);
}
#if defined(__cplusplus)
}
#endif
