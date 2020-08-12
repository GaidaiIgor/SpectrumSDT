#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* mfnbasic.c */
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
#define mfnreasonviewfromoptions_ MFNREASONVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnreasonviewfromoptions_ mfnreasonviewfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfncreate_ MFNCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfncreate_ mfncreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfnreset_ MFNRESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnreset_ mfnreset
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfndestroy_ MFNDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfndestroy_ mfndestroy
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfnsetbv_ MFNSETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnsetbv_ mfnsetbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfngetbv_ MFNGETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfngetbv_ mfngetbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfnsetfn_ MFNSETFN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfnsetfn_ mfnsetfn
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mfngetfn_ MFNGETFN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mfngetfn_ mfngetfn
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  mfnreasonviewfromoptions_(MFN mfn, int *__ierr){
*__ierr = MFNReasonViewFromOptions(
	(MFN)PetscToPointer((mfn) ));
}
SLEPC_EXTERN void  mfncreate_(MPI_Fint * comm,MFN *outmfn, int *__ierr){
*__ierr = MFNCreate(
	MPI_Comm_f2c(*(comm)),outmfn);
}
SLEPC_EXTERN void  mfnreset_(MFN mfn, int *__ierr){
*__ierr = MFNReset(
	(MFN)PetscToPointer((mfn) ));
}
SLEPC_EXTERN void  mfndestroy_(MFN *mfn, int *__ierr){
*__ierr = MFNDestroy(mfn);
}
SLEPC_EXTERN void  mfnsetbv_(MFN mfn,BV bv, int *__ierr){
*__ierr = MFNSetBV(
	(MFN)PetscToPointer((mfn) ),
	(BV)PetscToPointer((bv) ));
}
SLEPC_EXTERN void  mfngetbv_(MFN mfn,BV *bv, int *__ierr){
*__ierr = MFNGetBV(
	(MFN)PetscToPointer((mfn) ),bv);
}
SLEPC_EXTERN void  mfnsetfn_(MFN mfn,FN fn, int *__ierr){
*__ierr = MFNSetFN(
	(MFN)PetscToPointer((mfn) ),
	(FN)PetscToPointer((fn) ));
}
SLEPC_EXTERN void  mfngetfn_(MFN mfn,FN *fn, int *__ierr){
*__ierr = MFNGetFN(
	(MFN)PetscToPointer((mfn) ),fn);
}
#if defined(__cplusplus)
}
#endif
