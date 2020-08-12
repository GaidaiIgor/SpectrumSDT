#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* dsops.c */
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

#include "slepcds.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetleadingdimension_ DSGETLEADINGDIMENSION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetleadingdimension_ dsgetleadingdimension
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssetstate_ DSSETSTATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssetstate_ dssetstate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetstate_ DSGETSTATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetstate_ dsgetstate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssetdimensions_ DSSETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssetdimensions_ dssetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetdimensions_ DSGETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetdimensions_ dsgetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dstruncate_ DSTRUNCATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dstruncate_ dstruncate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsmatgetsize_ DSMATGETSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsmatgetsize_ dsmatgetsize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsmatishermitian_ DSMATISHERMITIAN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsmatishermitian_ dsmatishermitian
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetmat_ DSGETMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetmat_ dsgetmat
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsrestoremat_ DSRESTOREMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsrestoremat_ dsrestoremat
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssolve_ DSSOLVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssolve_ dssolve
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssynchronize_ DSSYNCHRONIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssynchronize_ dssynchronize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsupdateextrarow_ DSUPDATEEXTRAROW
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsupdateextrarow_ dsupdateextrarow
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dscond_ DSCOND
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dscond_ dscond
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dstranslaterks_ DSTRANSLATERKS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dstranslaterks_ dstranslaterks
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dscopymat_ DSCOPYMAT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dscopymat_ dscopymat
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  dsgetleadingdimension_(DS ds,PetscInt *ld, int *__ierr){
*__ierr = DSGetLeadingDimension(
	(DS)PetscToPointer((ds) ),ld);
}
SLEPC_EXTERN void  dssetstate_(DS ds,DSStateType *state, int *__ierr){
*__ierr = DSSetState(
	(DS)PetscToPointer((ds) ),*state);
}
SLEPC_EXTERN void  dsgetstate_(DS ds,DSStateType *state, int *__ierr){
*__ierr = DSGetState(
	(DS)PetscToPointer((ds) ),state);
}
SLEPC_EXTERN void  dssetdimensions_(DS ds,PetscInt *n,PetscInt *m,PetscInt *l,PetscInt *k, int *__ierr){
*__ierr = DSSetDimensions(
	(DS)PetscToPointer((ds) ),*n,*m,*l,*k);
}
SLEPC_EXTERN void  dsgetdimensions_(DS ds,PetscInt *n,PetscInt *m,PetscInt *l,PetscInt *k,PetscInt *t, int *__ierr){
*__ierr = DSGetDimensions(
	(DS)PetscToPointer((ds) ),n,m,l,k,t);
}
SLEPC_EXTERN void  dstruncate_(DS ds,PetscInt *n, int *__ierr){
*__ierr = DSTruncate(
	(DS)PetscToPointer((ds) ),*n);
}
SLEPC_EXTERN void  dsmatgetsize_(DS ds,DSMatType *t,PetscInt *m,PetscInt *n, int *__ierr){
*__ierr = DSMatGetSize(
	(DS)PetscToPointer((ds) ),*t,m,n);
}
SLEPC_EXTERN void  dsmatishermitian_(DS ds,DSMatType *t,PetscBool *flg, int *__ierr){
*__ierr = DSMatIsHermitian(
	(DS)PetscToPointer((ds) ),*t,flg);
}
SLEPC_EXTERN void  dsgetmat_(DS ds,DSMatType *m,Mat *A, int *__ierr){
*__ierr = DSGetMat(
	(DS)PetscToPointer((ds) ),*m,A);
}
SLEPC_EXTERN void  dsrestoremat_(DS ds,DSMatType *m,Mat *A, int *__ierr){
*__ierr = DSRestoreMat(
	(DS)PetscToPointer((ds) ),*m,A);
}
SLEPC_EXTERN void  dssolve_(DS ds,PetscScalar eigr[],PetscScalar eigi[], int *__ierr){
*__ierr = DSSolve(
	(DS)PetscToPointer((ds) ),eigr,eigi);
}
SLEPC_EXTERN void  dssynchronize_(DS ds,PetscScalar eigr[],PetscScalar eigi[], int *__ierr){
*__ierr = DSSynchronize(
	(DS)PetscToPointer((ds) ),eigr,eigi);
}
SLEPC_EXTERN void  dsupdateextrarow_(DS ds, int *__ierr){
*__ierr = DSUpdateExtraRow(
	(DS)PetscToPointer((ds) ));
}
SLEPC_EXTERN void  dscond_(DS ds,PetscReal *cond, int *__ierr){
*__ierr = DSCond(
	(DS)PetscToPointer((ds) ),cond);
}
SLEPC_EXTERN void  dstranslaterks_(DS ds,PetscScalar *alpha, int *__ierr){
*__ierr = DSTranslateRKS(
	(DS)PetscToPointer((ds) ),*alpha);
}
SLEPC_EXTERN void  dscopymat_(DS ds,DSMatType *m,PetscInt *mr,PetscInt *mc,Mat A,PetscInt *Ar,PetscInt *Ac,PetscInt *rows,PetscInt *cols,PetscBool *out, int *__ierr){
*__ierr = DSCopyMat(
	(DS)PetscToPointer((ds) ),*m,*mr,*mc,
	(Mat)PetscToPointer((A) ),*Ar,*Ac,*rows,*cols,*out);
}
#if defined(__cplusplus)
}
#endif
