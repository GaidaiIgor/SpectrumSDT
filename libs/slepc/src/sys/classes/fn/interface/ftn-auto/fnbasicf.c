#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* fnbasic.c */
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

#include "slepcfn.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fncreate_ FNCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fncreate_ fncreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnsetscale_ FNSETSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnsetscale_ fnsetscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fngetscale_ FNGETSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fngetscale_ fngetscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnsetmethod_ FNSETMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnsetmethod_ fnsetmethod
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fngetmethod_ FNGETMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fngetmethod_ fngetmethod
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnsetparallel_ FNSETPARALLEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnsetparallel_ fnsetparallel
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fngetparallel_ FNGETPARALLEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fngetparallel_ fngetparallel
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnevaluatefunction_ FNEVALUATEFUNCTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnevaluatefunction_ fnevaluatefunction
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnevaluatederivative_ FNEVALUATEDERIVATIVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnevaluatederivative_ fnevaluatederivative
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnevaluatefunctionmat_ FNEVALUATEFUNCTIONMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnevaluatefunctionmat_ fnevaluatefunctionmat
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnevaluatefunctionmatvec_ FNEVALUATEFUNCTIONMATVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnevaluatefunctionmatvec_ fnevaluatefunctionmatvec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnsetfromoptions_ FNSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnsetfromoptions_ fnsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fnduplicate_ FNDUPLICATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fnduplicate_ fnduplicate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define fndestroy_ FNDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define fndestroy_ fndestroy
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  fncreate_(MPI_Fint * comm,FN *newfn, int *__ierr){
*__ierr = FNCreate(
	MPI_Comm_f2c(*(comm)),newfn);
}
SLEPC_EXTERN void  fnsetscale_(FN fn,PetscScalar *alpha,PetscScalar *beta, int *__ierr){
*__ierr = FNSetScale(
	(FN)PetscToPointer((fn) ),*alpha,*beta);
}
SLEPC_EXTERN void  fngetscale_(FN fn,PetscScalar *alpha,PetscScalar *beta, int *__ierr){
*__ierr = FNGetScale(
	(FN)PetscToPointer((fn) ),alpha,beta);
}
SLEPC_EXTERN void  fnsetmethod_(FN fn,PetscInt *meth, int *__ierr){
*__ierr = FNSetMethod(
	(FN)PetscToPointer((fn) ),*meth);
}
SLEPC_EXTERN void  fngetmethod_(FN fn,PetscInt *meth, int *__ierr){
*__ierr = FNGetMethod(
	(FN)PetscToPointer((fn) ),meth);
}
SLEPC_EXTERN void  fnsetparallel_(FN fn,FNParallelType *pmode, int *__ierr){
*__ierr = FNSetParallel(
	(FN)PetscToPointer((fn) ),*pmode);
}
SLEPC_EXTERN void  fngetparallel_(FN fn,FNParallelType *pmode, int *__ierr){
*__ierr = FNGetParallel(
	(FN)PetscToPointer((fn) ),pmode);
}
SLEPC_EXTERN void  fnevaluatefunction_(FN fn,PetscScalar *x,PetscScalar *y, int *__ierr){
*__ierr = FNEvaluateFunction(
	(FN)PetscToPointer((fn) ),*x,y);
}
SLEPC_EXTERN void  fnevaluatederivative_(FN fn,PetscScalar *x,PetscScalar *y, int *__ierr){
*__ierr = FNEvaluateDerivative(
	(FN)PetscToPointer((fn) ),*x,y);
}
SLEPC_EXTERN void  fnevaluatefunctionmat_(FN fn,Mat A,Mat B, int *__ierr){
*__ierr = FNEvaluateFunctionMat(
	(FN)PetscToPointer((fn) ),
	(Mat)PetscToPointer((A) ),
	(Mat)PetscToPointer((B) ));
}
SLEPC_EXTERN void  fnevaluatefunctionmatvec_(FN fn,Mat A,Vec v, int *__ierr){
*__ierr = FNEvaluateFunctionMatVec(
	(FN)PetscToPointer((fn) ),
	(Mat)PetscToPointer((A) ),
	(Vec)PetscToPointer((v) ));
}
SLEPC_EXTERN void  fnsetfromoptions_(FN fn, int *__ierr){
*__ierr = FNSetFromOptions(
	(FN)PetscToPointer((fn) ));
}
SLEPC_EXTERN void  fnduplicate_(FN fn,MPI_Fint * comm,FN *newfn, int *__ierr){
*__ierr = FNDuplicate(
	(FN)PetscToPointer((fn) ),
	MPI_Comm_f2c(*(comm)),newfn);
}
SLEPC_EXTERN void  fndestroy_(FN *fn, int *__ierr){
*__ierr = FNDestroy(fn);
}
#if defined(__cplusplus)
}
#endif
