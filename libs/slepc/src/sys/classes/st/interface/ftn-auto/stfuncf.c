#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* stfunc.c */
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

#include "slepcst.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define streset_ STRESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define streset_ streset
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stdestroy_ STDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stdestroy_ stdestroy
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stcreate_ STCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stcreate_ stcreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stsetmatrices_ STSETMATRICES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stsetmatrices_ stsetmatrices
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stgetmatrix_ STGETMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stgetmatrix_ stgetmatrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stgetmatrixtransformed_ STGETMATRIXTRANSFORMED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stgetmatrixtransformed_ stgetmatrixtransformed
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stgetnummatrices_ STGETNUMMATRICES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stgetnummatrices_ stgetnummatrices
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stresetmatrixstate_ STRESETMATRIXSTATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stresetmatrixstate_ stresetmatrixstate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stsetshift_ STSETSHIFT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stsetshift_ stsetshift
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stgetshift_ STGETSHIFT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stgetshift_ stgetshift
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stsetdefaultshift_ STSETDEFAULTSHIFT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stsetdefaultshift_ stsetdefaultshift
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stscaleshift_ STSCALESHIFT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stscaleshift_ stscaleshift
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stsetbalancematrix_ STSETBALANCEMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stsetbalancematrix_ stsetbalancematrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stgetbalancematrix_ STGETBALANCEMATRIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stgetbalancematrix_ stgetbalancematrix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stmatgetsize_ STMATGETSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stmatgetsize_ stmatgetsize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stmatgetlocalsize_ STMATGETLOCALSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stmatgetlocalsize_ stmatgetlocalsize
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  streset_(ST st, int *__ierr){
*__ierr = STReset(
	(ST)PetscToPointer((st) ));
}
SLEPC_EXTERN void  stdestroy_(ST *st, int *__ierr){
*__ierr = STDestroy(st);
}
SLEPC_EXTERN void  stcreate_(MPI_Fint * comm,ST *newst, int *__ierr){
*__ierr = STCreate(
	MPI_Comm_f2c(*(comm)),newst);
}
SLEPC_EXTERN void  stsetmatrices_(ST st,PetscInt *n,Mat A[], int *__ierr){
*__ierr = STSetMatrices(
	(ST)PetscToPointer((st) ),*n,A);
}
SLEPC_EXTERN void  stgetmatrix_(ST st,PetscInt *k,Mat *A, int *__ierr){
*__ierr = STGetMatrix(
	(ST)PetscToPointer((st) ),*k,A);
}
SLEPC_EXTERN void  stgetmatrixtransformed_(ST st,PetscInt *k,Mat *T, int *__ierr){
*__ierr = STGetMatrixTransformed(
	(ST)PetscToPointer((st) ),*k,T);
}
SLEPC_EXTERN void  stgetnummatrices_(ST st,PetscInt *n, int *__ierr){
*__ierr = STGetNumMatrices(
	(ST)PetscToPointer((st) ),n);
}
SLEPC_EXTERN void  stresetmatrixstate_(ST st, int *__ierr){
*__ierr = STResetMatrixState(
	(ST)PetscToPointer((st) ));
}
SLEPC_EXTERN void  stsetshift_(ST st,PetscScalar *shift, int *__ierr){
*__ierr = STSetShift(
	(ST)PetscToPointer((st) ),*shift);
}
SLEPC_EXTERN void  stgetshift_(ST st,PetscScalar* shift, int *__ierr){
*__ierr = STGetShift(
	(ST)PetscToPointer((st) ),shift);
}
SLEPC_EXTERN void  stsetdefaultshift_(ST st,PetscScalar *defaultshift, int *__ierr){
*__ierr = STSetDefaultShift(
	(ST)PetscToPointer((st) ),*defaultshift);
}
SLEPC_EXTERN void  stscaleshift_(ST st,PetscScalar *factor, int *__ierr){
*__ierr = STScaleShift(
	(ST)PetscToPointer((st) ),*factor);
}
SLEPC_EXTERN void  stsetbalancematrix_(ST st,Vec D, int *__ierr){
*__ierr = STSetBalanceMatrix(
	(ST)PetscToPointer((st) ),
	(Vec)PetscToPointer((D) ));
}
SLEPC_EXTERN void  stgetbalancematrix_(ST st,Vec *D, int *__ierr){
*__ierr = STGetBalanceMatrix(
	(ST)PetscToPointer((st) ),D);
}
SLEPC_EXTERN void  stmatgetsize_(ST st,PetscInt *m,PetscInt *n, int *__ierr){
*__ierr = STMatGetSize(
	(ST)PetscToPointer((st) ),m,n);
}
SLEPC_EXTERN void  stmatgetlocalsize_(ST st,PetscInt *m,PetscInt *n, int *__ierr){
*__ierr = STMatGetLocalSize(
	(ST)PetscToPointer((st) ),m,n);
}
#if defined(__cplusplus)
}
#endif
