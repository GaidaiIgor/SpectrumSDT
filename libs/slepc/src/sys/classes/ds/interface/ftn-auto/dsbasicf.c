#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* dsbasic.c */
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
#define dscreate_ DSCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dscreate_ dscreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsduplicate_ DSDUPLICATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsduplicate_ dsduplicate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssetmethod_ DSSETMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssetmethod_ dssetmethod
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetmethod_ DSGETMETHOD
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetmethod_ dsgetmethod
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssetparallel_ DSSETPARALLEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssetparallel_ dssetparallel
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetparallel_ DSGETPARALLEL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetparallel_ dsgetparallel
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssetcompact_ DSSETCOMPACT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssetcompact_ dssetcompact
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetcompact_ DSGETCOMPACT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetcompact_ dsgetcompact
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssetextrarow_ DSSETEXTRAROW
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssetextrarow_ dssetextrarow
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetextrarow_ DSGETEXTRAROW
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetextrarow_ dsgetextrarow
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssetrefined_ DSSETREFINED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssetrefined_ dssetrefined
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetrefined_ DSGETREFINED
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetrefined_ dsgetrefined
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssetblocksize_ DSSETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssetblocksize_ dssetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsgetblocksize_ DSGETBLOCKSIZE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsgetblocksize_ dsgetblocksize
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dssetfromoptions_ DSSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dssetfromoptions_ dssetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsallocate_ DSALLOCATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsallocate_ dsallocate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsreset_ DSRESET
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsreset_ dsreset
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define dsdestroy_ DSDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define dsdestroy_ dsdestroy
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
SLEPC_EXTERN void  dscreate_(MPI_Fint * comm,DS *newds, int *__ierr){
*__ierr = DSCreate(
	MPI_Comm_f2c(*(comm)),newds);
}
SLEPC_EXTERN void  dsduplicate_(DS ds,DS *dsnew, int *__ierr){
*__ierr = DSDuplicate(
	(DS)PetscToPointer((ds) ),dsnew);
}
SLEPC_EXTERN void  dssetmethod_(DS ds,PetscInt *meth, int *__ierr){
*__ierr = DSSetMethod(
	(DS)PetscToPointer((ds) ),*meth);
}
SLEPC_EXTERN void  dsgetmethod_(DS ds,PetscInt *meth, int *__ierr){
*__ierr = DSGetMethod(
	(DS)PetscToPointer((ds) ),meth);
}
SLEPC_EXTERN void  dssetparallel_(DS ds,DSParallelType *pmode, int *__ierr){
*__ierr = DSSetParallel(
	(DS)PetscToPointer((ds) ),*pmode);
}
SLEPC_EXTERN void  dsgetparallel_(DS ds,DSParallelType *pmode, int *__ierr){
*__ierr = DSGetParallel(
	(DS)PetscToPointer((ds) ),pmode);
}
SLEPC_EXTERN void  dssetcompact_(DS ds,PetscBool *comp, int *__ierr){
*__ierr = DSSetCompact(
	(DS)PetscToPointer((ds) ),*comp);
}
SLEPC_EXTERN void  dsgetcompact_(DS ds,PetscBool *comp, int *__ierr){
*__ierr = DSGetCompact(
	(DS)PetscToPointer((ds) ),comp);
}
SLEPC_EXTERN void  dssetextrarow_(DS ds,PetscBool *ext, int *__ierr){
*__ierr = DSSetExtraRow(
	(DS)PetscToPointer((ds) ),*ext);
}
SLEPC_EXTERN void  dsgetextrarow_(DS ds,PetscBool *ext, int *__ierr){
*__ierr = DSGetExtraRow(
	(DS)PetscToPointer((ds) ),ext);
}
SLEPC_EXTERN void  dssetrefined_(DS ds,PetscBool *ref, int *__ierr){
*__ierr = DSSetRefined(
	(DS)PetscToPointer((ds) ),*ref);
}
SLEPC_EXTERN void  dsgetrefined_(DS ds,PetscBool *ref, int *__ierr){
*__ierr = DSGetRefined(
	(DS)PetscToPointer((ds) ),ref);
}
SLEPC_EXTERN void  dssetblocksize_(DS ds,PetscInt *bs, int *__ierr){
*__ierr = DSSetBlockSize(
	(DS)PetscToPointer((ds) ),*bs);
}
SLEPC_EXTERN void  dsgetblocksize_(DS ds,PetscInt *bs, int *__ierr){
*__ierr = DSGetBlockSize(
	(DS)PetscToPointer((ds) ),bs);
}
SLEPC_EXTERN void  dssetfromoptions_(DS ds, int *__ierr){
*__ierr = DSSetFromOptions(
	(DS)PetscToPointer((ds) ));
}
SLEPC_EXTERN void  dsallocate_(DS ds,PetscInt *ld, int *__ierr){
*__ierr = DSAllocate(
	(DS)PetscToPointer((ds) ),*ld);
}
SLEPC_EXTERN void  dsreset_(DS ds, int *__ierr){
*__ierr = DSReset(
	(DS)PetscToPointer((ds) ));
}
SLEPC_EXTERN void  dsdestroy_(DS *ds, int *__ierr){
*__ierr = DSDestroy(ds);
}
#if defined(__cplusplus)
}
#endif
