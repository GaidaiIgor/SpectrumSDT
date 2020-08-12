#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* lmebasic.c */
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

#include "slepclme.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmereasonviewfromoptions_ LMEREASONVIEWFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmereasonviewfromoptions_ lmereasonviewfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmecreate_ LMECREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmecreate_ lmecreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmereset_ LMERESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmereset_ lmereset
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmedestroy_ LMEDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmedestroy_ lmedestroy
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmesetbv_ LMESETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmesetbv_ lmesetbv
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define lmegetbv_ LMEGETBV
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define lmegetbv_ lmegetbv
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  lmereasonviewfromoptions_(LME lme, int *__ierr){
*__ierr = LMEReasonViewFromOptions(
	(LME)PetscToPointer((lme) ));
}
SLEPC_EXTERN void  lmecreate_(MPI_Fint * comm,LME *outlme, int *__ierr){
*__ierr = LMECreate(
	MPI_Comm_f2c(*(comm)),outlme);
}
SLEPC_EXTERN void  lmereset_(LME lme, int *__ierr){
*__ierr = LMEReset(
	(LME)PetscToPointer((lme) ));
}
SLEPC_EXTERN void  lmedestroy_(LME *lme, int *__ierr){
*__ierr = LMEDestroy(lme);
}
SLEPC_EXTERN void  lmesetbv_(LME lme,BV bv, int *__ierr){
*__ierr = LMESetBV(
	(LME)PetscToPointer((lme) ),
	(BV)PetscToPointer((bv) ));
}
SLEPC_EXTERN void  lmegetbv_(LME lme,BV *bv, int *__ierr){
*__ierr = LMEGetBV(
	(LME)PetscToPointer((lme) ),bv);
}
#if defined(__cplusplus)
}
#endif
