/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepcst.h>
#include <slepcds.h>
#include <slepcfn.h>
#include <slepcbv.h>
#include <slepcrg.h>

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)

#if defined(PETSC_USE_SINGLE_LIBRARY)
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepceps(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcnep(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcpep(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcsvd(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcmfn(void);
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepclme(void);
#endif

/*
  PetscDLLibraryRegister - This function is called when the dynamic library
  it is in is opened.

  This one registers all the basic objects ST, FN, DS, BV, RG.
 */
#if defined(PETSC_USE_SINGLE_LIBRARY)
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepc(void)
#else
SLEPC_EXTERN PetscErrorCode PetscDLLibraryRegister_slepcsys(void)
#endif
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = STInitializePackage();CHKERRQ(ierr);
  ierr = DSInitializePackage();CHKERRQ(ierr);
  ierr = FNInitializePackage();CHKERRQ(ierr);
  ierr = BVInitializePackage();CHKERRQ(ierr);
  ierr = RGInitializePackage();CHKERRQ(ierr);

#if defined(PETSC_USE_SINGLE_LIBRARY)
  ierr = PetscDLLibraryRegister_slepceps();CHKERRQ(ierr);
  ierr = PetscDLLibraryRegister_slepcnep();CHKERRQ(ierr);
  ierr = PetscDLLibraryRegister_slepcpep();CHKERRQ(ierr);
  ierr = PetscDLLibraryRegister_slepcsvd();CHKERRQ(ierr);
  ierr = PetscDLLibraryRegister_slepcmfn();CHKERRQ(ierr);
  ierr = PetscDLLibraryRegister_slepclme();CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}
#endif /* PETSC_HAVE_DYNAMIC_LIBRARIES */

