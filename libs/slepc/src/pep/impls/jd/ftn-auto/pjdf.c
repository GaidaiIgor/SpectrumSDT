#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* pjd.c */
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
#define pepjdsetrestart_ PEPJDSETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdsetrestart_ pepjdsetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepjdgetrestart_ PEPJDGETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdgetrestart_ pepjdgetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepjdsetfix_ PEPJDSETFIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdsetfix_ pepjdsetfix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepjdgetfix_ PEPJDGETFIX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdgetfix_ pepjdgetfix
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepjdsetreusepreconditioner_ PEPJDSETREUSEPRECONDITIONER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdsetreusepreconditioner_ pepjdsetreusepreconditioner
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepjdgetreusepreconditioner_ PEPJDGETREUSEPRECONDITIONER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdgetreusepreconditioner_ pepjdgetreusepreconditioner
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepjdsetminimalityindex_ PEPJDSETMINIMALITYINDEX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdsetminimalityindex_ pepjdsetminimalityindex
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepjdgetminimalityindex_ PEPJDGETMINIMALITYINDEX
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdgetminimalityindex_ pepjdgetminimalityindex
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepjdsetprojection_ PEPJDSETPROJECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdsetprojection_ pepjdsetprojection
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepjdgetprojection_ PEPJDGETPROJECTION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepjdgetprojection_ pepjdgetprojection
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  pepjdsetrestart_(PEP pep,PetscReal *keep, int *__ierr){
*__ierr = PEPJDSetRestart(
	(PEP)PetscToPointer((pep) ),*keep);
}
SLEPC_EXTERN void  pepjdgetrestart_(PEP pep,PetscReal *keep, int *__ierr){
*__ierr = PEPJDGetRestart(
	(PEP)PetscToPointer((pep) ),keep);
}
SLEPC_EXTERN void  pepjdsetfix_(PEP pep,PetscReal *fix, int *__ierr){
*__ierr = PEPJDSetFix(
	(PEP)PetscToPointer((pep) ),*fix);
}
SLEPC_EXTERN void  pepjdgetfix_(PEP pep,PetscReal *fix, int *__ierr){
*__ierr = PEPJDGetFix(
	(PEP)PetscToPointer((pep) ),fix);
}
SLEPC_EXTERN void  pepjdsetreusepreconditioner_(PEP pep,PetscBool *reusepc, int *__ierr){
*__ierr = PEPJDSetReusePreconditioner(
	(PEP)PetscToPointer((pep) ),*reusepc);
}
SLEPC_EXTERN void  pepjdgetreusepreconditioner_(PEP pep,PetscBool *reusepc, int *__ierr){
*__ierr = PEPJDGetReusePreconditioner(
	(PEP)PetscToPointer((pep) ),reusepc);
}
SLEPC_EXTERN void  pepjdsetminimalityindex_(PEP pep,PetscInt *mmidx, int *__ierr){
*__ierr = PEPJDSetMinimalityIndex(
	(PEP)PetscToPointer((pep) ),*mmidx);
}
SLEPC_EXTERN void  pepjdgetminimalityindex_(PEP pep,PetscInt *mmidx, int *__ierr){
*__ierr = PEPJDGetMinimalityIndex(
	(PEP)PetscToPointer((pep) ),mmidx);
}
SLEPC_EXTERN void  pepjdsetprojection_(PEP pep,PEPJDProjection *proj, int *__ierr){
*__ierr = PEPJDSetProjection(
	(PEP)PetscToPointer((pep) ),*proj);
}
SLEPC_EXTERN void  pepjdgetprojection_(PEP pep,PEPJDProjection *proj, int *__ierr){
*__ierr = PEPJDGetProjection(
	(PEP)PetscToPointer((pep) ),proj);
}
#if defined(__cplusplus)
}
#endif
