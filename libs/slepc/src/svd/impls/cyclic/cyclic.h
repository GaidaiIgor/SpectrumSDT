/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if !defined(SLEPC_CYCLIC_H)
#define SLEPC_CYCLIC_H

typedef struct {
  PetscBool explicitmatrix;
  EPS       eps;
  PetscBool usereps;
  Mat       mat;
  Vec       x1,x2,y1,y2;
} SVD_CYCLIC;

SLEPC_INTERN PetscErrorCode MatMult_Cyclic_CUDA(Mat,Vec,Vec);
SLEPC_INTERN PetscErrorCode MatCreateVecs_Cyclic_CUDA(Mat,Vec*,Vec*);

#endif
