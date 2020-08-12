/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Private data structure used by the BLZPACK interface
*/

#if !defined(SLEPC_BLZPACK_H)
#define SLEPC_BLZPACK_H

typedef struct {
  PetscBLASInt         block_size;      /* block size */
  PetscBLASInt         slice;           /* use spectrum slicing */
  PetscBLASInt         nsteps;          /* maximum number of steps per run */
  PetscBLASInt         *istor;
  PetscReal            *rstor;
  PetscScalar          *u;
  PetscScalar          *v;
  PetscScalar          *eig;
} EPS_BLZPACK;

/*
   Definition of routines from the BLZPACK package
*/

#if defined(SLEPC_BLZPACK_HAVE_UNDERSCORE)
#define SLEPC_BLZPACK(lcase,ucase) lcase##_
#elif defined(SLEPC_BLZPACK_HAVE_CAPS)
#define SLEPC_BLZPACK(lcase,ucase) ucase
#else
#define SLEPC_BLZPACK(lcase,ucase) lcase
#endif

/*
    These are real case, current version of BLZPACK only supports real
    matrices
*/

#if defined(PETSC_USE_SINGLE)
/*
   For these machines we must call the single precision Fortran version
*/
#define BLZpack_ SLEPC_BLZPACK(blzdrs,BLZDRS)
#else
#define BLZpack_ SLEPC_BLZPACK(blzdrd,BLZDRD)
#endif

#define BLZistorr_ SLEPC_BLZPACK(istorr,ISTORR)
#define BLZrstorr_ SLEPC_BLZPACK(rstorr,RSTORR)

SLEPC_EXTERN void BLZpack_(PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*);

SLEPC_EXTERN PetscBLASInt BLZistorr_(PetscBLASInt*,const char*,int);
SLEPC_EXTERN PetscReal BLZrstorr_(PetscReal*,char*,int);

#endif

