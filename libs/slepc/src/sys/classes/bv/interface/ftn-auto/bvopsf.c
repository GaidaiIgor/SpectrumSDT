#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* bvops.c */
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
#define bvmult_ BVMULT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmult_ bvmult
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmultvec_ BVMULTVEC
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmultvec_ bvmultvec
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmultcolumn_ BVMULTCOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmultcolumn_ bvmultcolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmultinplace_ BVMULTINPLACE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmultinplace_ bvmultinplace
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmultinplacetranspose_ BVMULTINPLACETRANSPOSE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmultinplacetranspose_ bvmultinplacetranspose
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvscale_ BVSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvscale_ bvscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvscalecolumn_ BVSCALECOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvscalecolumn_ bvscalecolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetrandom_ BVSETRANDOM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetrandom_ bvsetrandom
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetrandomcolumn_ BVSETRANDOMCOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetrandomcolumn_ bvsetrandomcolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvsetrandomcond_ BVSETRANDOMCOND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvsetrandomcond_ bvsetrandomcond
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmatmult_ BVMATMULT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmatmult_ bvmatmult
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmatmulttranspose_ BVMATMULTTRANSPOSE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmatmulttranspose_ bvmatmulttranspose
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmatmulthermitiantranspose_ BVMATMULTHERMITIANTRANSPOSE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmatmulthermitiantranspose_ bvmatmulthermitiantranspose
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmatmultcolumn_ BVMATMULTCOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmatmultcolumn_ bvmatmultcolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmatmulttransposecolumn_ BVMATMULTTRANSPOSECOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmatmulttransposecolumn_ bvmatmulttransposecolumn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define bvmatmulthermitiantransposecolumn_ BVMATMULTHERMITIANTRANSPOSECOLUMN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define bvmatmulthermitiantransposecolumn_ bvmatmulthermitiantransposecolumn
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  bvmult_(BV Y,PetscScalar *alpha,PetscScalar *beta,BV X,Mat Q, int *__ierr){
*__ierr = BVMult(
	(BV)PetscToPointer((Y) ),*alpha,*beta,
	(BV)PetscToPointer((X) ),
	(Mat)PetscToPointer((Q) ));
}
SLEPC_EXTERN void  bvmultvec_(BV X,PetscScalar *alpha,PetscScalar *beta,Vec y,PetscScalar q[], int *__ierr){
*__ierr = BVMultVec(
	(BV)PetscToPointer((X) ),*alpha,*beta,
	(Vec)PetscToPointer((y) ),q);
}
SLEPC_EXTERN void  bvmultcolumn_(BV X,PetscScalar *alpha,PetscScalar *beta,PetscInt *j,PetscScalar *q, int *__ierr){
*__ierr = BVMultColumn(
	(BV)PetscToPointer((X) ),*alpha,*beta,*j,q);
}
SLEPC_EXTERN void  bvmultinplace_(BV V,Mat Q,PetscInt *s,PetscInt *e, int *__ierr){
*__ierr = BVMultInPlace(
	(BV)PetscToPointer((V) ),
	(Mat)PetscToPointer((Q) ),*s,*e);
}
SLEPC_EXTERN void  bvmultinplacetranspose_(BV V,Mat Q,PetscInt *s,PetscInt *e, int *__ierr){
*__ierr = BVMultInPlaceTranspose(
	(BV)PetscToPointer((V) ),
	(Mat)PetscToPointer((Q) ),*s,*e);
}
SLEPC_EXTERN void  bvscale_(BV bv,PetscScalar *alpha, int *__ierr){
*__ierr = BVScale(
	(BV)PetscToPointer((bv) ),*alpha);
}
SLEPC_EXTERN void  bvscalecolumn_(BV bv,PetscInt *j,PetscScalar *alpha, int *__ierr){
*__ierr = BVScaleColumn(
	(BV)PetscToPointer((bv) ),*j,*alpha);
}
SLEPC_EXTERN void  bvsetrandom_(BV bv, int *__ierr){
*__ierr = BVSetRandom(
	(BV)PetscToPointer((bv) ));
}
SLEPC_EXTERN void  bvsetrandomcolumn_(BV bv,PetscInt *j, int *__ierr){
*__ierr = BVSetRandomColumn(
	(BV)PetscToPointer((bv) ),*j);
}
SLEPC_EXTERN void  bvsetrandomcond_(BV bv,PetscReal *condn, int *__ierr){
*__ierr = BVSetRandomCond(
	(BV)PetscToPointer((bv) ),*condn);
}
SLEPC_EXTERN void  bvmatmult_(BV V,Mat A,BV Y, int *__ierr){
*__ierr = BVMatMult(
	(BV)PetscToPointer((V) ),
	(Mat)PetscToPointer((A) ),
	(BV)PetscToPointer((Y) ));
}
SLEPC_EXTERN void  bvmatmulttranspose_(BV V,Mat A,BV Y, int *__ierr){
*__ierr = BVMatMultTranspose(
	(BV)PetscToPointer((V) ),
	(Mat)PetscToPointer((A) ),
	(BV)PetscToPointer((Y) ));
}
SLEPC_EXTERN void  bvmatmulthermitiantranspose_(BV V,Mat A,BV Y, int *__ierr){
*__ierr = BVMatMultHermitianTranspose(
	(BV)PetscToPointer((V) ),
	(Mat)PetscToPointer((A) ),
	(BV)PetscToPointer((Y) ));
}
SLEPC_EXTERN void  bvmatmultcolumn_(BV V,Mat A,PetscInt *j, int *__ierr){
*__ierr = BVMatMultColumn(
	(BV)PetscToPointer((V) ),
	(Mat)PetscToPointer((A) ),*j);
}
SLEPC_EXTERN void  bvmatmulttransposecolumn_(BV V,Mat A,PetscInt *j, int *__ierr){
*__ierr = BVMatMultTransposeColumn(
	(BV)PetscToPointer((V) ),
	(Mat)PetscToPointer((A) ),*j);
}
SLEPC_EXTERN void  bvmatmulthermitiantransposecolumn_(BV V,Mat A,PetscInt *j, int *__ierr){
*__ierr = BVMatMultHermitianTransposeColumn(
	(BV)PetscToPointer((V) ),
	(Mat)PetscToPointer((A) ),*j);
}
#if defined(__cplusplus)
}
#endif
