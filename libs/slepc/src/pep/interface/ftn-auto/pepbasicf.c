#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* pepbasic.c */
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

#include "slepcpep.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepcreate_ PEPCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepcreate_ pepcreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepreset_ PEPRESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepreset_ pepreset
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepdestroy_ PEPDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepdestroy_ pepdestroy
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetbv_ PEPSETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetbv_ pepsetbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetbv_ PEPGETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetbv_ pepgetbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetrg_ PEPSETRG
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetrg_ pepsetrg
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetrg_ PEPGETRG
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetrg_ pepgetrg
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetds_ PEPSETDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetds_ pepsetds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetds_ PEPGETDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetds_ pepgetds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetst_ PEPSETST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetst_ pepsetst
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetst_ PEPGETST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetst_ pepgetst
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define peprefinegetksp_ PEPREFINEGETKSP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define peprefinegetksp_ peprefinegetksp
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsettarget_ PEPSETTARGET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsettarget_ pepsettarget
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgettarget_ PEPGETTARGET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgettarget_ pepgettarget
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetinterval_ PEPSETINTERVAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetinterval_ pepsetinterval
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetinterval_ PEPGETINTERVAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetinterval_ pepgetinterval
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  pepcreate_(MPI_Fint * comm,PEP *outpep, int *__ierr){
*__ierr = PEPCreate(
	MPI_Comm_f2c(*(comm)),outpep);
}
SLEPC_EXTERN void  pepreset_(PEP pep, int *__ierr){
*__ierr = PEPReset(
	(PEP)PetscToPointer((pep) ));
}
SLEPC_EXTERN void  pepdestroy_(PEP *pep, int *__ierr){
*__ierr = PEPDestroy(pep);
}
SLEPC_EXTERN void  pepsetbv_(PEP pep,BV bv, int *__ierr){
*__ierr = PEPSetBV(
	(PEP)PetscToPointer((pep) ),
	(BV)PetscToPointer((bv) ));
}
SLEPC_EXTERN void  pepgetbv_(PEP pep,BV *bv, int *__ierr){
*__ierr = PEPGetBV(
	(PEP)PetscToPointer((pep) ),bv);
}
SLEPC_EXTERN void  pepsetrg_(PEP pep,RG rg, int *__ierr){
*__ierr = PEPSetRG(
	(PEP)PetscToPointer((pep) ),
	(RG)PetscToPointer((rg) ));
}
SLEPC_EXTERN void  pepgetrg_(PEP pep,RG *rg, int *__ierr){
*__ierr = PEPGetRG(
	(PEP)PetscToPointer((pep) ),rg);
}
SLEPC_EXTERN void  pepsetds_(PEP pep,DS ds, int *__ierr){
*__ierr = PEPSetDS(
	(PEP)PetscToPointer((pep) ),
	(DS)PetscToPointer((ds) ));
}
SLEPC_EXTERN void  pepgetds_(PEP pep,DS *ds, int *__ierr){
*__ierr = PEPGetDS(
	(PEP)PetscToPointer((pep) ),ds);
}
SLEPC_EXTERN void  pepsetst_(PEP pep,ST st, int *__ierr){
*__ierr = PEPSetST(
	(PEP)PetscToPointer((pep) ),
	(ST)PetscToPointer((st) ));
}
SLEPC_EXTERN void  pepgetst_(PEP pep,ST *st, int *__ierr){
*__ierr = PEPGetST(
	(PEP)PetscToPointer((pep) ),st);
}
SLEPC_EXTERN void  peprefinegetksp_(PEP pep,KSP *ksp, int *__ierr){
*__ierr = PEPRefineGetKSP(
	(PEP)PetscToPointer((pep) ),ksp);
}
SLEPC_EXTERN void  pepsettarget_(PEP pep,PetscScalar *target, int *__ierr){
*__ierr = PEPSetTarget(
	(PEP)PetscToPointer((pep) ),*target);
}
SLEPC_EXTERN void  pepgettarget_(PEP pep,PetscScalar* target, int *__ierr){
*__ierr = PEPGetTarget(
	(PEP)PetscToPointer((pep) ),target);
}
SLEPC_EXTERN void  pepsetinterval_(PEP pep,PetscReal *inta,PetscReal *intb, int *__ierr){
*__ierr = PEPSetInterval(
	(PEP)PetscToPointer((pep) ),*inta,*intb);
}
SLEPC_EXTERN void  pepgetinterval_(PEP pep,PetscReal* inta,PetscReal* intb, int *__ierr){
*__ierr = PEPGetInterval(
	(PEP)PetscToPointer((pep) ),inta,intb);
}
#if defined(__cplusplus)
}
#endif
