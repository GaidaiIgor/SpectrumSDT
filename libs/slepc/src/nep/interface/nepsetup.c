/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to problem setup
*/

#include <slepc/private/nepimpl.h>       /*I "slepcnep.h" I*/

/*@
   NEPSetUp - Sets up all the internal data structures necessary for the
   execution of the NEP solver.

   Collective on nep

   Input Parameter:
.  nep   - solver context

   Notes:
   This function need not be called explicitly in most cases, since NEPSolve()
   calls it. It can be useful when one wants to measure the set-up time
   separately from the solve time.

   Level: developer

.seealso: NEPCreate(), NEPSolve(), NEPDestroy()
@*/
PetscErrorCode NEPSetUp(NEP nep)
{
  PetscErrorCode ierr;
  PetscInt       k;
  SlepcSC        sc;
  Mat            T;
  PetscBool      flg;
  KSP            ksp;
  PC             pc;
  PetscMPIInt    size;
  MatSolverType  stype;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  NEPCheckProblem(nep,1);
  if (nep->state) PetscFunctionReturn(0);
  ierr = PetscLogEventBegin(NEP_SetUp,nep,0,0,0);CHKERRQ(ierr);

  /* reset the convergence flag from the previous solves */
  nep->reason = NEP_CONVERGED_ITERATING;

  /* set default solver type (NEPSetFromOptions was not called) */
  if (!((PetscObject)nep)->type_name) {
    ierr = NEPSetType(nep,NEPRII);CHKERRQ(ierr);
  }
  if (nep->useds && !nep->ds) { ierr = NEPGetDS(nep,&nep->ds);CHKERRQ(ierr); }
  if (!nep->rg) { ierr = NEPGetRG(nep,&nep->rg);CHKERRQ(ierr); }
  if (!((PetscObject)nep->rg)->type_name) {
    ierr = RGSetType(nep->rg,RGINTERVAL);CHKERRQ(ierr);
  }
  if (nep->twosided && !nep->hasts) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver does not support computing left eigenvectors (no two-sided variant)");

  /* set problem dimensions */
  switch (nep->fui) {
  case NEP_USER_INTERFACE_CALLBACK:
    ierr = NEPGetFunction(nep,&T,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = MatGetSize(T,&nep->n,NULL);CHKERRQ(ierr);
    ierr = MatGetLocalSize(T,&nep->nloc,NULL);CHKERRQ(ierr);
    break;
  case NEP_USER_INTERFACE_SPLIT:
    ierr = MatDuplicate(nep->A[0],MAT_DO_NOT_COPY_VALUES,&nep->function);CHKERRQ(ierr);
    ierr = MatDuplicate(nep->A[0],MAT_DO_NOT_COPY_VALUES,&nep->jacobian);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)nep->function);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)nep->jacobian);CHKERRQ(ierr);
    ierr = MatGetSize(nep->A[0],&nep->n,NULL);CHKERRQ(ierr);
    ierr = MatGetLocalSize(nep->A[0],&nep->nloc,NULL);CHKERRQ(ierr);
    break;
  }

  /* set default problem type */
  if (!nep->problem_type) {
    ierr = NEPSetProblemType(nep,NEP_GENERAL);CHKERRQ(ierr);
  }

  /* check consistency of refinement options */
  if (nep->refine) {
    if (nep->fui!=NEP_USER_INTERFACE_SPLIT) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Iterative refinement only implemented in split form");
    if (!nep->scheme) {  /* set default scheme */
      ierr = NEPRefineGetKSP(nep,&ksp);CHKERRQ(ierr);
      ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
      ierr = PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&flg);CHKERRQ(ierr);
      if (flg) {
        ierr = PetscObjectTypeCompareAny((PetscObject)pc,&flg,PCLU,PCCHOLESKY,"");CHKERRQ(ierr);
      }
      nep->scheme = flg? NEP_REFINE_SCHEME_MBE: NEP_REFINE_SCHEME_SCHUR;
    }
    if (nep->scheme==NEP_REFINE_SCHEME_MBE) {
      ierr = NEPRefineGetKSP(nep,&ksp);CHKERRQ(ierr);
      ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
      ierr = PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&flg);CHKERRQ(ierr);
      if (flg) {
        ierr = PetscObjectTypeCompareAny((PetscObject)pc,&flg,PCLU,PCCHOLESKY,"");CHKERRQ(ierr);
      }
      if (!flg) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The MBE scheme for refinement requires a direct solver in KSP");
      ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pc),&size);CHKERRQ(ierr);
      if (size>1) {   /* currently selected PC is a factorization */
        ierr = PCFactorGetMatSolverType(pc,&stype);CHKERRQ(ierr);
        ierr = PetscStrcmp(stype,MATSOLVERPETSC,&flg);CHKERRQ(ierr);
        if (flg) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"For Newton refinement, you chose to solve linear systems with a factorization, but in parallel runs you need to select an external package");
      }
    }
    if (nep->scheme==NEP_REFINE_SCHEME_SCHUR) {
      if (nep->npart>1) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The Schur scheme for refinement does not support subcommunicators");
    }
  }
  /* call specific solver setup */
  ierr = (*nep->ops->setup)(nep);CHKERRQ(ierr);

  /* by default, compute eigenvalues close to target */
  /* nep->target should contain the initial guess for the eigenvalue */
  if (!nep->which) nep->which = NEP_TARGET_MAGNITUDE;

  /* set tolerance if not yet set */
  if (nep->tol==PETSC_DEFAULT) nep->tol = SLEPC_DEFAULT_TOL;
  if (nep->refine) {
    if (nep->rtol==PETSC_DEFAULT) nep->rtol = PetscMax(nep->tol/1000,PETSC_MACHINE_EPSILON);
    if (nep->rits==PETSC_DEFAULT) nep->rits = (nep->refine==NEP_REFINE_SIMPLE)? 10: 1;
  }

  /* fill sorting criterion context */
  switch (nep->which) {
    case NEP_LARGEST_MAGNITUDE:
      nep->sc->comparison    = SlepcCompareLargestMagnitude;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_SMALLEST_MAGNITUDE:
      nep->sc->comparison    = SlepcCompareSmallestMagnitude;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_LARGEST_REAL:
      nep->sc->comparison    = SlepcCompareLargestReal;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_SMALLEST_REAL:
      nep->sc->comparison    = SlepcCompareSmallestReal;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_LARGEST_IMAGINARY:
      nep->sc->comparison    = SlepcCompareLargestImaginary;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_SMALLEST_IMAGINARY:
      nep->sc->comparison    = SlepcCompareSmallestImaginary;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_TARGET_MAGNITUDE:
      nep->sc->comparison    = SlepcCompareTargetMagnitude;
      nep->sc->comparisonctx = &nep->target;
      break;
    case NEP_TARGET_REAL:
      nep->sc->comparison    = SlepcCompareTargetReal;
      nep->sc->comparisonctx = &nep->target;
      break;
    case NEP_TARGET_IMAGINARY:
#if defined(PETSC_USE_COMPLEX)
      nep->sc->comparison    = SlepcCompareTargetImaginary;
      nep->sc->comparisonctx = &nep->target;
#endif
      break;
    case NEP_ALL:
      nep->sc->comparison    = SlepcCompareSmallestReal;
      nep->sc->comparisonctx = NULL;
      break;
    case NEP_WHICH_USER:
      break;
  }

  nep->sc->map    = NULL;
  nep->sc->mapobj = NULL;

  /* fill sorting criterion for DS */
  if (nep->useds) {
    ierr = DSGetSlepcSC(nep->ds,&sc);CHKERRQ(ierr);
    sc->comparison    = nep->sc->comparison;
    sc->comparisonctx = nep->sc->comparisonctx;
    ierr = PetscObjectTypeCompare((PetscObject)nep,NEPNLEIGS,&flg);CHKERRQ(ierr);
    if (!flg) {
      sc->map    = NULL;
      sc->mapobj = NULL;
    }
  }
  if (nep->nev > nep->ncv) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"nev bigger than ncv");

  /* process initial vectors */
  if (nep->nini<0) {
    k = -nep->nini;
    if (k>nep->ncv) SETERRQ(PetscObjectComm((PetscObject)nep),1,"The number of initial vectors is larger than ncv");
    ierr = BVInsertVecs(nep->V,0,&k,nep->IS,PETSC_TRUE);CHKERRQ(ierr);
    ierr = SlepcBasisDestroy_Private(&nep->nini,&nep->IS);CHKERRQ(ierr);
    nep->nini = k;
  }
  ierr = PetscLogEventEnd(NEP_SetUp,nep,0,0,0);CHKERRQ(ierr);
  nep->state = NEP_STATE_SETUP;
  PetscFunctionReturn(0);
}

/*@C
   NEPSetInitialSpace - Specify a basis of vectors that constitute the initial
   space, that is, the subspace from which the solver starts to iterate.

   Collective on nep

   Input Parameter:
+  nep   - the nonlinear eigensolver context
.  n     - number of vectors
-  is    - set of basis vectors of the initial space

   Notes:
   Some solvers start to iterate on a single vector (initial vector). In that case,
   the other vectors are ignored.

   These vectors do not persist from one NEPSolve() call to the other, so the
   initial space should be set every time.

   The vectors do not need to be mutually orthonormal, since they are explicitly
   orthonormalized internally.

   Common usage of this function is when the user can provide a rough approximation
   of the wanted eigenspace. Then, convergence may be faster.

   Level: intermediate
@*/
PetscErrorCode NEPSetInitialSpace(NEP nep,PetscInt n,Vec is[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(nep,n,2);
  if (n<0) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative");
  if (n>0) {
    PetscValidPointer(is,3);
    PetscValidHeaderSpecific(*is,VEC_CLASSID,3);
  }
  ierr = SlepcBasisReference_Private(n,is,&nep->nini,&nep->IS);CHKERRQ(ierr);
  if (n>0) nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*
  NEPSetDimensions_Default - Set reasonable values for ncv, mpd if not set
  by the user. This is called at setup.
 */
PetscErrorCode NEPSetDimensions_Default(NEP nep,PetscInt nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  if (*ncv!=PETSC_DEFAULT) { /* ncv set */
    if (*ncv<nev) SETERRQ(PetscObjectComm((PetscObject)nep),1,"The value of ncv must be at least nev");
  } else if (*mpd!=PETSC_DEFAULT) { /* mpd set */
    *ncv = PetscMin(nep->n,nev+(*mpd));
  } else { /* neither set: defaults depend on nev being small or large */
    if (nev<500) *ncv = PetscMin(nep->n,PetscMax(2*nev,nev+15));
    else {
      *mpd = 500;
      *ncv = PetscMin(nep->n,nev+(*mpd));
    }
  }
  if (*mpd==PETSC_DEFAULT) *mpd = *ncv;
  PetscFunctionReturn(0);
}

/*@
   NEPAllocateSolution - Allocate memory storage for common variables such
   as eigenvalues and eigenvectors.

   Collective on nep

   Input Parameters:
+  nep   - eigensolver context
-  extra - number of additional positions, used for methods that require a
           working basis slightly larger than ncv

   Developers Note:
   This is SLEPC_EXTERN because it may be required by user plugin NEP
   implementations.

   Level: developer
@*/
PetscErrorCode NEPAllocateSolution(NEP nep,PetscInt extra)
{
  PetscErrorCode ierr;
  PetscInt       oldsize,newc,requested;
  PetscLogDouble cnt;
  Mat            T;
  Vec            t;

  PetscFunctionBegin;
  requested = nep->ncv + extra;

  /* oldsize is zero if this is the first time setup is called */
  ierr = BVGetSizes(nep->V,NULL,NULL,&oldsize);CHKERRQ(ierr);
  newc = PetscMax(0,requested-oldsize);

  /* allocate space for eigenvalues and friends */
  if (requested != oldsize || !nep->eigr) {
    ierr = PetscFree4(nep->eigr,nep->eigi,nep->errest,nep->perm);CHKERRQ(ierr);
    ierr = PetscMalloc4(requested,&nep->eigr,requested,&nep->eigi,requested,&nep->errest,requested,&nep->perm);CHKERRQ(ierr);
    cnt = newc*sizeof(PetscScalar) + newc*sizeof(PetscReal) + newc*sizeof(PetscInt);
    ierr = PetscLogObjectMemory((PetscObject)nep,cnt);CHKERRQ(ierr);
  }

  /* allocate V */
  if (!nep->V) { ierr = NEPGetBV(nep,&nep->V);CHKERRQ(ierr); }
  if (!oldsize) {
    if (!((PetscObject)(nep->V))->type_name) {
      ierr = BVSetType(nep->V,BVSVEC);CHKERRQ(ierr);
    }
    if (nep->fui==NEP_USER_INTERFACE_SPLIT) T = nep->A[0];
    else {
      ierr = NEPGetFunction(nep,&T,NULL,NULL,NULL);CHKERRQ(ierr);
    }
    ierr = MatCreateVecsEmpty(T,&t,NULL);CHKERRQ(ierr);
    ierr = BVSetSizesFromVec(nep->V,t,requested);CHKERRQ(ierr);
    ierr = VecDestroy(&t);CHKERRQ(ierr);
  } else {
    ierr = BVResize(nep->V,requested,PETSC_FALSE);CHKERRQ(ierr);
  }

  /* allocate W */
  if (nep->twosided) {
    ierr = BVDestroy(&nep->W);CHKERRQ(ierr);
    ierr = BVDuplicate(nep->V,&nep->W);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

