#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* epsbasic.c */
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

#include "slepceps.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epscreate_ EPSCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epscreate_ epscreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsreset_ EPSRESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsreset_ epsreset
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsdestroy_ EPSDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsdestroy_ epsdestroy
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssettarget_ EPSSETTARGET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssettarget_ epssettarget
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgettarget_ EPSGETTARGET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgettarget_ epsgettarget
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetinterval_ EPSSETINTERVAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetinterval_ epssetinterval
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetinterval_ EPSGETINTERVAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetinterval_ epsgetinterval
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetst_ EPSSETST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetst_ epssetst
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetst_ EPSGETST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetst_ epsgetst
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetbv_ EPSSETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetbv_ epssetbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetbv_ EPSGETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetbv_ epsgetbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetrg_ EPSSETRG
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetrg_ epssetrg
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetrg_ EPSGETRG
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetrg_ epsgetrg
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epssetds_ EPSSETDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epssetds_ epssetds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsgetds_ EPSGETDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsgetds_ epsgetds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsisgeneralized_ EPSISGENERALIZED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsisgeneralized_ epsisgeneralized
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsishermitian_ EPSISHERMITIAN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsishermitian_ epsishermitian
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define epsispositive_ EPSISPOSITIVE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define epsispositive_ epsispositive
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  epscreate_(MPI_Fint * comm,EPS *outeps, int *__ierr){
*__ierr = EPSCreate(
	MPI_Comm_f2c(*(comm)),outeps);
}
SLEPC_EXTERN void  epsreset_(EPS eps, int *__ierr){
*__ierr = EPSReset(
	(EPS)PetscToPointer((eps) ));
}
SLEPC_EXTERN void  epsdestroy_(EPS *eps, int *__ierr){
*__ierr = EPSDestroy(eps);
}
SLEPC_EXTERN void  epssettarget_(EPS eps,PetscScalar *target, int *__ierr){
*__ierr = EPSSetTarget(
	(EPS)PetscToPointer((eps) ),*target);
}
SLEPC_EXTERN void  epsgettarget_(EPS eps,PetscScalar* target, int *__ierr){
*__ierr = EPSGetTarget(
	(EPS)PetscToPointer((eps) ),target);
}
SLEPC_EXTERN void  epssetinterval_(EPS eps,PetscReal *inta,PetscReal *intb, int *__ierr){
*__ierr = EPSSetInterval(
	(EPS)PetscToPointer((eps) ),*inta,*intb);
}
SLEPC_EXTERN void  epsgetinterval_(EPS eps,PetscReal* inta,PetscReal* intb, int *__ierr){
*__ierr = EPSGetInterval(
	(EPS)PetscToPointer((eps) ),inta,intb);
}
SLEPC_EXTERN void  epssetst_(EPS eps,ST st, int *__ierr){
*__ierr = EPSSetST(
	(EPS)PetscToPointer((eps) ),
	(ST)PetscToPointer((st) ));
}
SLEPC_EXTERN void  epsgetst_(EPS eps,ST *st, int *__ierr){
*__ierr = EPSGetST(
	(EPS)PetscToPointer((eps) ),st);
}
SLEPC_EXTERN void  epssetbv_(EPS eps,BV V, int *__ierr){
*__ierr = EPSSetBV(
	(EPS)PetscToPointer((eps) ),
	(BV)PetscToPointer((V) ));
}
SLEPC_EXTERN void  epsgetbv_(EPS eps,BV *V, int *__ierr){
*__ierr = EPSGetBV(
	(EPS)PetscToPointer((eps) ),V);
}
SLEPC_EXTERN void  epssetrg_(EPS eps,RG rg, int *__ierr){
*__ierr = EPSSetRG(
	(EPS)PetscToPointer((eps) ),
	(RG)PetscToPointer((rg) ));
}
SLEPC_EXTERN void  epsgetrg_(EPS eps,RG *rg, int *__ierr){
*__ierr = EPSGetRG(
	(EPS)PetscToPointer((eps) ),rg);
}
SLEPC_EXTERN void  epssetds_(EPS eps,DS ds, int *__ierr){
*__ierr = EPSSetDS(
	(EPS)PetscToPointer((eps) ),
	(DS)PetscToPointer((ds) ));
}
SLEPC_EXTERN void  epsgetds_(EPS eps,DS *ds, int *__ierr){
*__ierr = EPSGetDS(
	(EPS)PetscToPointer((eps) ),ds);
}
SLEPC_EXTERN void  epsisgeneralized_(EPS eps,PetscBool* is, int *__ierr){
*__ierr = EPSIsGeneralized(
	(EPS)PetscToPointer((eps) ),is);
}
SLEPC_EXTERN void  epsishermitian_(EPS eps,PetscBool* is, int *__ierr){
*__ierr = EPSIsHermitian(
	(EPS)PetscToPointer((eps) ),is);
}
SLEPC_EXTERN void  epsispositive_(EPS eps,PetscBool* is, int *__ierr){
*__ierr = EPSIsPositive(
	(EPS)PetscToPointer((eps) ),is);
}
#if defined(__cplusplus)
}
#endif
