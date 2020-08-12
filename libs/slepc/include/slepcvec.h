/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   User interface for various vector operations added in SLEPc
*/

#if !defined(SLEPCVEC_H)
#define SLEPCVEC_H
#include <slepcsys.h>

/* VecComp: Vec composed of several smaller Vecs */
#define VECCOMP  "comp"

SLEPC_EXTERN PetscErrorCode VecCreateComp(MPI_Comm,PetscInt*,PetscInt,VecType,Vec,Vec*);
SLEPC_EXTERN PetscErrorCode VecCreateCompWithVecs(Vec*,PetscInt,Vec,Vec*);
SLEPC_EXTERN PetscErrorCode VecCompGetSubVecs(Vec,PetscInt*,const Vec**);
SLEPC_EXTERN PetscErrorCode VecCompSetSubVecs(Vec,PetscInt,Vec*);

/* Some auxiliary functions */
SLEPC_EXTERN PetscErrorCode VecNormalizeComplex(Vec,Vec,PetscBool,PetscReal*);
SLEPC_EXTERN PetscErrorCode VecCheckOrthogonality(Vec[],PetscInt,Vec[],PetscInt,Mat,PetscViewer,PetscReal*);
SLEPC_EXTERN PetscErrorCode VecDuplicateEmpty(Vec,Vec*);

/* Deprecated functions */
PETSC_DEPRECATED_FUNCTION("Use VecNormalizeComplex()") PETSC_STATIC_INLINE PetscErrorCode SlepcVecNormalize(Vec xr,Vec xi,PetscBool c,PetscReal *nrm) {return VecNormalizeComplex(xr,xi,c,nrm);}
PETSC_DEPRECATED_FUNCTION("Use VecCheckOrthogonality()") PETSC_STATIC_INLINE PetscErrorCode SlepcCheckOrthogonality(Vec *V,PetscInt nv,Vec *W,PetscInt nw,Mat B,PetscViewer viewer,PetscReal *lev) {return VecCheckOrthogonality(V,nv,W,nw,B,viewer,lev);}

#endif

