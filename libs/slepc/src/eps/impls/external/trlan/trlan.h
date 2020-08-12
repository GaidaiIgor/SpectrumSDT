/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Private data structure used by the TRLAN interface
*/

#if !defined(SLEPC_TRLAN_H)
#define SLEPC_TRLAN_H

typedef struct {
  PetscBLASInt       maxlan;
  PetscBLASInt       restart;
  PetscReal          *work;
  PetscBLASInt       lwork;
} EPS_TRLAN;

/*
   Definition of routines from the TRLAN package
   These are real case. TRLAN currently only has DOUBLE PRECISION version
*/

#if defined(SLEPC_TRLAN_HAVE_UNDERSCORE)
#define TRLan_ trlan77_
#elif defined(SLEPC_TRLAN_HAVE_CAPS)
#define TRLan_ TRLAN77
#else
#define TRLan_ trlan77
#endif

SLEPC_EXTERN void TRLan_(PetscBLASInt(*op)(PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*),PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);

#endif

