/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Modification of the *temp* implementation of the BLOPEX multivector in order
   to wrap created PETSc vectors as multivectors.
*/

#include <slepc/private/bvimpl.h>
#include <stdlib.h>
#include <interpreter.h>
#include <temp_multivector.h>
#include "blopex.h"

static void* mv_TempMultiVectorCreateFromBV(void* ii_,BlopexInt n,void* sample)
{
  PetscErrorCode          ierr;
  BV                      bv = (BV)sample;
  Vec                     v;
  PetscInt                i,l,k,nc,useconstr=PETSC_FALSE,flg;
  mv_TempMultiVector      *x;
  mv_InterfaceInterpreter *ii = (mv_InterfaceInterpreter*)ii_;

  x = (mv_TempMultiVector*)malloc(sizeof(mv_TempMultiVector));
  if (!x) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Allocation for x failed");

  x->interpreter = ii;
  x->numVectors  = n;

  x->vector = (void**)calloc(n,sizeof(void*));
  if (!x->vector) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Allocation for x->vector failed");

  x->ownsVectors = 1;
  x->mask = NULL;
  x->ownsMask = 0;

  ierr = BVGetActiveColumns(bv,&l,&k);CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = PetscObjectComposedDataGetInt((PetscObject)bv,slepc_blopex_useconstr,useconstr,flg);CHKERRABORT(PETSC_COMM_SELF,ierr);
  if (!l && useconstr) {
    ierr = BVGetNumConstraints(bv,&nc);CHKERRABORT(PETSC_COMM_SELF,ierr);
    l = -nc;
  }
  if (n != k-l) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_PLIB,"BV active columns plus constraints do not match argument n");
  for (i=0;i<n;i++) {
    ierr = BVGetColumn(bv,l+i,&v);CHKERRABORT(PETSC_COMM_SELF,ierr);
    ierr = PetscObjectReference((PetscObject)v);CHKERRABORT(PETSC_COMM_SELF,ierr);
    x->vector[i] = (void*)v;
    ierr = BVRestoreColumn(bv,l+i,&v);CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  return x;
}

/*
    Create an InterfaceInterpreter using the PETSc implementation
    but overloading CreateMultiVector that doesn't create any
    new vector.
*/
int SLEPCSetupInterpreter(mv_InterfaceInterpreter *i)
{
  PETSCSetupInterpreter(i);
  i->CreateMultiVector = mv_TempMultiVectorCreateFromBV;

  return 0;
}

