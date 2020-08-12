#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* stset.c */
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
#define stsetfromoptions_ STSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stsetfromoptions_ stsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stsetmatstructure_ STSETMATSTRUCTURE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stsetmatstructure_ stsetmatstructure
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stgetmatstructure_ STGETMATSTRUCTURE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stgetmatstructure_ stgetmatstructure
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stsetmatmode_ STSETMATMODE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stsetmatmode_ stsetmatmode
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stgetmatmode_ STGETMATMODE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stgetmatmode_ stgetmatmode
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stsettransform_ STSETTRANSFORM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stsettransform_ stsettransform
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define stgettransform_ STGETTRANSFORM
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define stgettransform_ stgettransform
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  stsetfromoptions_(ST st, int *__ierr){
*__ierr = STSetFromOptions(
	(ST)PetscToPointer((st) ));
}
SLEPC_EXTERN void  stsetmatstructure_(ST st,MatStructure *str, int *__ierr){
*__ierr = STSetMatStructure(
	(ST)PetscToPointer((st) ),*str);
}
SLEPC_EXTERN void  stgetmatstructure_(ST st,MatStructure *str, int *__ierr){
*__ierr = STGetMatStructure(
	(ST)PetscToPointer((st) ),str);
}
SLEPC_EXTERN void  stsetmatmode_(ST st,STMatMode *mode, int *__ierr){
*__ierr = STSetMatMode(
	(ST)PetscToPointer((st) ),*mode);
}
SLEPC_EXTERN void  stgetmatmode_(ST st,STMatMode *mode, int *__ierr){
*__ierr = STGetMatMode(
	(ST)PetscToPointer((st) ),mode);
}
SLEPC_EXTERN void  stsettransform_(ST st,PetscBool *flg, int *__ierr){
*__ierr = STSetTransform(
	(ST)PetscToPointer((st) ),*flg);
}
SLEPC_EXTERN void  stgettransform_(ST st,PetscBool *flg, int *__ierr){
*__ierr = STGetTransform(
	(ST)PetscToPointer((st) ),flg);
}
#if defined(__cplusplus)
}
#endif
