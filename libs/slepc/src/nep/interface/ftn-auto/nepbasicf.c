#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* nepbasic.c */
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

#include "slepcnep.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepcreate_ NEPCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepcreate_ nepcreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepreset_ NEPRESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepreset_ nepreset
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepdestroy_ NEPDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepdestroy_ nepdestroy
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetbv_ NEPSETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetbv_ nepsetbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgetbv_ NEPGETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgetbv_ nepgetbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetrg_ NEPSETRG
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetrg_ nepsetrg
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgetrg_ NEPGETRG
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgetrg_ nepgetrg
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetds_ NEPSETDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetds_ nepsetds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgetds_ NEPGETDS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgetds_ nepgetds
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define neprefinegetksp_ NEPREFINEGETKSP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define neprefinegetksp_ neprefinegetksp
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsettarget_ NEPSETTARGET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsettarget_ nepsettarget
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgettarget_ NEPGETTARGET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgettarget_ nepgettarget
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepsetsplitoperator_ NEPSETSPLITOPERATOR
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepsetsplitoperator_ nepsetsplitoperator
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgetsplitoperatorterm_ NEPGETSPLITOPERATORTERM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgetsplitoperatorterm_ nepgetsplitoperatorterm
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepgetsplitoperatorinfo_ NEPGETSPLITOPERATORINFO
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepgetsplitoperatorinfo_ nepgetsplitoperatorinfo
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  nepcreate_(MPI_Fint * comm,NEP *outnep, int *__ierr){
*__ierr = NEPCreate(
	MPI_Comm_f2c(*(comm)),outnep);
}
SLEPC_EXTERN void  nepreset_(NEP nep, int *__ierr){
*__ierr = NEPReset(
	(NEP)PetscToPointer((nep) ));
}
SLEPC_EXTERN void  nepdestroy_(NEP *nep, int *__ierr){
*__ierr = NEPDestroy(nep);
}
SLEPC_EXTERN void  nepsetbv_(NEP nep,BV bv, int *__ierr){
*__ierr = NEPSetBV(
	(NEP)PetscToPointer((nep) ),
	(BV)PetscToPointer((bv) ));
}
SLEPC_EXTERN void  nepgetbv_(NEP nep,BV *bv, int *__ierr){
*__ierr = NEPGetBV(
	(NEP)PetscToPointer((nep) ),bv);
}
SLEPC_EXTERN void  nepsetrg_(NEP nep,RG rg, int *__ierr){
*__ierr = NEPSetRG(
	(NEP)PetscToPointer((nep) ),
	(RG)PetscToPointer((rg) ));
}
SLEPC_EXTERN void  nepgetrg_(NEP nep,RG *rg, int *__ierr){
*__ierr = NEPGetRG(
	(NEP)PetscToPointer((nep) ),rg);
}
SLEPC_EXTERN void  nepsetds_(NEP nep,DS ds, int *__ierr){
*__ierr = NEPSetDS(
	(NEP)PetscToPointer((nep) ),
	(DS)PetscToPointer((ds) ));
}
SLEPC_EXTERN void  nepgetds_(NEP nep,DS *ds, int *__ierr){
*__ierr = NEPGetDS(
	(NEP)PetscToPointer((nep) ),ds);
}
SLEPC_EXTERN void  neprefinegetksp_(NEP nep,KSP *ksp, int *__ierr){
*__ierr = NEPRefineGetKSP(
	(NEP)PetscToPointer((nep) ),ksp);
}
SLEPC_EXTERN void  nepsettarget_(NEP nep,PetscScalar *target, int *__ierr){
*__ierr = NEPSetTarget(
	(NEP)PetscToPointer((nep) ),*target);
}
SLEPC_EXTERN void  nepgettarget_(NEP nep,PetscScalar* target, int *__ierr){
*__ierr = NEPGetTarget(
	(NEP)PetscToPointer((nep) ),target);
}
SLEPC_EXTERN void  nepsetsplitoperator_(NEP nep,PetscInt *n,Mat A[],FN f[],MatStructure *str, int *__ierr){
*__ierr = NEPSetSplitOperator(
	(NEP)PetscToPointer((nep) ),*n,A,f,*str);
}
SLEPC_EXTERN void  nepgetsplitoperatorterm_(NEP nep,PetscInt *k,Mat *A,FN *f, int *__ierr){
*__ierr = NEPGetSplitOperatorTerm(
	(NEP)PetscToPointer((nep) ),*k,A,f);
}
SLEPC_EXTERN void  nepgetsplitoperatorinfo_(NEP nep,PetscInt *n,MatStructure *str, int *__ierr){
*__ierr = NEPGetSplitOperatorInfo(
	(NEP)PetscToPointer((nep) ),n,str);
}
#if defined(__cplusplus)
}
#endif
