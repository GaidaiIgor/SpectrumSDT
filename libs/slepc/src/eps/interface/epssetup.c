/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   EPS routines related to problem setup
*/

#include <slepc/private/epsimpl.h>       /*I "slepceps.h" I*/

/*
   Let the solver choose the ST type that should be used by default,
   otherwise set it to SHIFT.
   This is called at EPSSetFromOptions (before STSetFromOptions)
   and also at EPSSetUp (in case EPSSetFromOptions was not called).
*/
PetscErrorCode EPSSetDefaultST(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (eps->ops->setdefaultst) { ierr = (*eps->ops->setdefaultst)(eps);CHKERRQ(ierr); }
  if (!((PetscObject)eps->st)->type_name) {
    ierr = STSetType(eps->st,STSHIFT);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   This is done by preconditioned eigensolvers that use the PC only.
   It sets STPRECOND with KSPPREONLY.
*/
PetscErrorCode EPSSetDefaultST_Precond(EPS eps)
{
  PetscErrorCode ierr;
  KSP            ksp;

  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) {
    ierr = STSetType(eps->st,STPRECOND);CHKERRQ(ierr);
  }
  ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr);
  if (!((PetscObject)ksp)->type_name) {
    ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   This is done by preconditioned eigensolvers that can also use the KSP.
   It sets STPRECOND with the default KSP (GMRES) and maxit=5.
*/
PetscErrorCode EPSSetDefaultST_GMRES(EPS eps)
{
  PetscErrorCode ierr;
  KSP            ksp;

  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) {
    ierr = STSetType(eps->st,STPRECOND);CHKERRQ(ierr);
  }
  ierr = STPrecondSetKSPHasMat(eps->st,PETSC_TRUE);CHKERRQ(ierr);
  ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr);
  if (!((PetscObject)ksp)->type_name) {
    ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,5);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   Check that the ST selected by the user is compatible with the EPS solver and options
*/
PetscErrorCode EPSCheckCompatibleST(EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      precond,shift,sinvert,cayley,lyapii;
#if defined(PETSC_USE_COMPLEX)
  PetscScalar    sigma;
#endif

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STPRECOND,&precond);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STSHIFT,&shift);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STSINVERT,&sinvert);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)eps->st,STCAYLEY,&cayley);CHKERRQ(ierr);

  /* preconditioned eigensolvers */
  if (eps->categ==EPS_CATEGORY_PRECOND && !precond) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver requires ST=PRECOND");
  if (eps->categ!=EPS_CATEGORY_PRECOND && precond) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"STPRECOND is intended for preconditioned eigensolvers only");

  /* harmonic extraction */
  if (!(precond || shift) && eps->extraction && eps->extraction!=EPS_RITZ) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Cannot use a spectral transformation combined with harmonic extraction");

  /* real shifts in Hermitian problems */
#if defined(PETSC_USE_COMPLEX)
  ierr = STGetShift(eps->st,&sigma);CHKERRQ(ierr);
  if (eps->ishermitian && PetscImaginaryPart(sigma) != 0.0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Hermitian problems are not compatible with complex shifts");
#endif

  /* Cayley with PGNHEP */
  if (cayley && eps->problem_type == EPS_PGNHEP) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Cayley spectral transformation is not compatible with PGNHEP");

  /* make sure that the user does not specify smallest magnitude with shift-and-invert */
  if ((cayley || sinvert) && (eps->categ==EPS_CATEGORY_KRYLOV || eps->categ==EPS_CATEGORY_OTHER)) {
    ierr = PetscObjectTypeCompare((PetscObject)eps,EPSLYAPII,&lyapii);CHKERRQ(ierr);
    if (!lyapii && eps->which!=EPS_TARGET_MAGNITUDE && eps->which!=EPS_TARGET_REAL && eps->which!=EPS_TARGET_IMAGINARY && eps->which!=EPS_ALL && eps->which!=EPS_WHICH_USER) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Shift-and-invert requires a target 'which' (see EPSSetWhichEigenpairs), for instance -st_type sinvert -eps_target 0 -eps_target_magnitude");
  }
  PetscFunctionReturn(0);
}

/*@
   EPSSetUp - Sets up all the internal data structures necessary for the
   execution of the eigensolver. Then calls STSetUp() for any set-up
   operations associated to the ST object.

   Collective on eps

   Input Parameter:
.  eps   - eigenproblem solver context

   Notes:
   This function need not be called explicitly in most cases, since EPSSolve()
   calls it. It can be useful when one wants to measure the set-up time
   separately from the solve time.

   Level: developer

.seealso: EPSCreate(), EPSSolve(), EPSDestroy(), STSetUp(), EPSSetInitialSpace()
@*/
PetscErrorCode EPSSetUp(EPS eps)
{
  PetscErrorCode ierr;
  Mat            A,B;
  SlepcSC        sc;
  PetscInt       k,nmat;
  PetscBool      flg,istrivial;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (eps->state) PetscFunctionReturn(0);
  ierr = PetscLogEventBegin(EPS_SetUp,eps,0,0,0);CHKERRQ(ierr);

  /* reset the convergence flag from the previous solves */
  eps->reason = EPS_CONVERGED_ITERATING;

  /* Set default solver type (EPSSetFromOptions was not called) */
  if (!((PetscObject)eps)->type_name) {
    ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRQ(ierr);
  }
  if (!eps->st) { ierr = EPSGetST(eps,&eps->st);CHKERRQ(ierr); }
  ierr = EPSSetDefaultST(eps);CHKERRQ(ierr);

  ierr = STSetTransform(eps->st,PETSC_TRUE);CHKERRQ(ierr);
  if (eps->useds && !eps->ds) { ierr = EPSGetDS(eps,&eps->ds);CHKERRQ(ierr); }
  if (eps->twosided) {
    if (!eps->hasts) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support computing left eigenvectors (no two-sided variant)");
    if (eps->ishermitian) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Two-sided methods are not intended for symmetric problems");
    if (!eps->dsts) { ierr = DSDuplicate(eps->ds,&eps->dsts);CHKERRQ(ierr); }
  }
  if (!eps->rg) { ierr = EPSGetRG(eps,&eps->rg);CHKERRQ(ierr); }
  if (!((PetscObject)eps->rg)->type_name) {
    ierr = RGSetType(eps->rg,RGINTERVAL);CHKERRQ(ierr);
  }

  /* Set problem dimensions */
  ierr = STGetNumMatrices(eps->st,&nmat);CHKERRQ(ierr);
  if (!nmat) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONGSTATE,"EPSSetOperators must be called first");
  ierr = STMatGetSize(eps->st,&eps->n,NULL);CHKERRQ(ierr);
  ierr = STMatGetLocalSize(eps->st,&eps->nloc,NULL);CHKERRQ(ierr);

  /* Set default problem type */
  if (!eps->problem_type) {
    if (nmat==1) {
      ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr);
    } else {
      ierr = EPSSetProblemType(eps,EPS_GNHEP);CHKERRQ(ierr);
    }
  } else if (nmat==1 && eps->isgeneralized) {
    ierr = PetscInfo(eps,"Eigenproblem set as generalized but no matrix B was provided; reverting to a standard eigenproblem\n");CHKERRQ(ierr);
    eps->isgeneralized = PETSC_FALSE;
    eps->problem_type = eps->ishermitian? EPS_HEP: EPS_NHEP;
  } else if (nmat>1 && !eps->isgeneralized) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_INCOMP,"Inconsistent EPS state");

  if (eps->nev > eps->n) eps->nev = eps->n;
  if (eps->ncv > eps->n) eps->ncv = eps->n;

  /* initialization of matrix norms */
  if (eps->conv==EPS_CONV_NORM) {
    if (!eps->nrma) {
      ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
      ierr = MatNorm(A,NORM_INFINITY,&eps->nrma);CHKERRQ(ierr);
    }
    if (nmat>1 && !eps->nrmb) {
      ierr = STGetMatrix(eps->st,1,&B);CHKERRQ(ierr);
      ierr = MatNorm(B,NORM_INFINITY,&eps->nrmb);CHKERRQ(ierr);
    }
  }

  /* call specific solver setup */
  ierr = (*eps->ops->setup)(eps);CHKERRQ(ierr);

  /* if purification is set, check that it really makes sense */
  if (eps->purify) {
    if (eps->categ==EPS_CATEGORY_PRECOND || eps->categ==EPS_CATEGORY_CONTOUR) eps->purify = PETSC_FALSE;
    else {
      if (!eps->isgeneralized) eps->purify = PETSC_FALSE;
      else if (!eps->ishermitian && !eps->ispositive) eps->purify = PETSC_FALSE;
      else {
        ierr = PetscObjectTypeCompare((PetscObject)eps->st,STCAYLEY,&flg);CHKERRQ(ierr);
        if (flg) eps->purify = PETSC_FALSE;
      }
    }
  }

  /* set tolerance if not yet set */
  if (eps->tol==PETSC_DEFAULT) eps->tol = SLEPC_DEFAULT_TOL;

  /* fill sorting criterion context */
  switch (eps->which) {
    case EPS_LARGEST_MAGNITUDE:
      eps->sc->comparison    = SlepcCompareLargestMagnitude;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_SMALLEST_MAGNITUDE:
      eps->sc->comparison    = SlepcCompareSmallestMagnitude;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_LARGEST_REAL:
      eps->sc->comparison    = SlepcCompareLargestReal;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_SMALLEST_REAL:
      eps->sc->comparison    = SlepcCompareSmallestReal;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_LARGEST_IMAGINARY:
      eps->sc->comparison    = SlepcCompareLargestImaginary;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_SMALLEST_IMAGINARY:
      eps->sc->comparison    = SlepcCompareSmallestImaginary;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_TARGET_MAGNITUDE:
      eps->sc->comparison    = SlepcCompareTargetMagnitude;
      eps->sc->comparisonctx = &eps->target;
      break;
    case EPS_TARGET_REAL:
      eps->sc->comparison    = SlepcCompareTargetReal;
      eps->sc->comparisonctx = &eps->target;
      break;
    case EPS_TARGET_IMAGINARY:
#if defined(PETSC_USE_COMPLEX)
      eps->sc->comparison    = SlepcCompareTargetImaginary;
      eps->sc->comparisonctx = &eps->target;
#endif
      break;
    case EPS_ALL:
      eps->sc->comparison    = SlepcCompareSmallestReal;
      eps->sc->comparisonctx = NULL;
      break;
    case EPS_WHICH_USER:
      break;
  }
  eps->sc->map    = NULL;
  eps->sc->mapobj = NULL;

  /* fill sorting criterion for DS */
  if (eps->useds && eps->which!=EPS_ALL) {
    ierr = DSGetSlepcSC(eps->ds,&sc);CHKERRQ(ierr);
    ierr = RGIsTrivial(eps->rg,&istrivial);CHKERRQ(ierr);
    sc->rg            = istrivial? NULL: eps->rg;
    sc->comparison    = eps->sc->comparison;
    sc->comparisonctx = eps->sc->comparisonctx;
    sc->map           = SlepcMap_ST;
    sc->mapobj        = (PetscObject)eps->st;
    if (eps->twosided) {
      ierr = DSSetSlepcSC(eps->dsts,sc);CHKERRQ(ierr);
    }
  }

  /* Build balancing matrix if required */
  if (eps->balance!=EPS_BALANCE_USER) {
    ierr = STSetBalanceMatrix(eps->st,NULL);CHKERRQ(ierr);
    if (!eps->ishermitian && (eps->balance==EPS_BALANCE_ONESIDE || eps->balance==EPS_BALANCE_TWOSIDE)) {
      if (!eps->D) {
        ierr = BVCreateVec(eps->V,&eps->D);CHKERRQ(ierr);
        ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)eps->D);CHKERRQ(ierr);
      }
      ierr = EPSBuildBalance_Krylov(eps);CHKERRQ(ierr);
      ierr = STSetBalanceMatrix(eps->st,eps->D);CHKERRQ(ierr);
    }
  }

  /* Setup ST */
  ierr = STSetUp(eps->st);CHKERRQ(ierr);
  ierr = EPSCheckCompatibleST(eps);CHKERRQ(ierr);

  /* process deflation and initial vectors */
  if (eps->nds<0) {
    k = -eps->nds;
    ierr = BVInsertConstraints(eps->V,&k,eps->defl);CHKERRQ(ierr);
    ierr = SlepcBasisDestroy_Private(&eps->nds,&eps->defl);CHKERRQ(ierr);
    eps->nds = k;
    ierr = STCheckNullSpace(eps->st,eps->V);CHKERRQ(ierr);
  }
  if (eps->nini<0) {
    k = -eps->nini;
    if (k>eps->ncv) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The number of initial vectors is larger than ncv");
    ierr = BVInsertVecs(eps->V,0,&k,eps->IS,PETSC_TRUE);CHKERRQ(ierr);
    ierr = SlepcBasisDestroy_Private(&eps->nini,&eps->IS);CHKERRQ(ierr);
    eps->nini = k;
  }
  if (eps->twosided && eps->ninil<0) {
    k = -eps->ninil;
    if (k>eps->ncv) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The number of left initial vectors is larger than ncv");
    ierr = BVInsertVecs(eps->W,0,&k,eps->ISL,PETSC_TRUE);CHKERRQ(ierr);
    ierr = SlepcBasisDestroy_Private(&eps->ninil,&eps->ISL);CHKERRQ(ierr);
    eps->ninil = k;
  }

  ierr = PetscLogEventEnd(EPS_SetUp,eps,0,0,0);CHKERRQ(ierr);
  eps->state = EPS_STATE_SETUP;
  PetscFunctionReturn(0);
}

/*@
   EPSSetOperators - Sets the matrices associated with the eigenvalue problem.

   Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
.  A  - the matrix associated with the eigensystem
-  B  - the second matrix in the case of generalized eigenproblems

   Notes:
   To specify a standard eigenproblem, use NULL for parameter B.

   It must be called before EPSSetUp(). If it is called again after EPSSetUp() then
   the EPS object is reset.

   Level: beginner

.seealso: EPSSolve(), EPSSetUp(), EPSReset(), EPSGetST(), STGetMatrix()
@*/
PetscErrorCode EPSSetOperators(EPS eps,Mat A,Mat B)
{
  PetscErrorCode ierr;
  PetscInt       m,n,m0,nmat;
  Mat            mat[2];

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  if (B) PetscValidHeaderSpecific(B,MAT_CLASSID,3);
  PetscCheckSameComm(eps,1,A,2);
  if (B) PetscCheckSameComm(eps,1,B,3);

  /* Check for square matrices */
  ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
  if (m!=n) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"A is a non-square matrix");
  if (B) {
    ierr = MatGetSize(B,&m0,&n);CHKERRQ(ierr);
    if (m0!=n) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"B is a non-square matrix");
    if (m!=m0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_INCOMP,"Dimensions of A and B do not match");
  }
  if (eps->state && n!=eps->n) { ierr = EPSReset(eps);CHKERRQ(ierr); }
  eps->nrma = 0.0;
  eps->nrmb = 0.0;
  if (!eps->st) { ierr = EPSGetST(eps,&eps->st);CHKERRQ(ierr); }
  mat[0] = A;
  if (B) {
    mat[1] = B;
    nmat = 2;
  } else nmat = 1;
  ierr = STSetMatrices(eps->st,nmat,mat);CHKERRQ(ierr);
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSGetOperators - Gets the matrices associated with the eigensystem.

   Collective on eps

   Input Parameter:
.  eps - the EPS context

   Output Parameters:
+  A  - the matrix associated with the eigensystem
-  B  - the second matrix in the case of generalized eigenproblems

   Level: intermediate

.seealso: EPSSolve(), EPSGetST(), STGetMatrix(), STSetMatrices()
@*/
PetscErrorCode EPSGetOperators(EPS eps,Mat *A,Mat *B)
{
  PetscErrorCode ierr;
  ST             st;
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  if (A) { ierr = STGetMatrix(st,0,A);CHKERRQ(ierr); }
  if (B) {
    ierr = STGetNumMatrices(st,&k);CHKERRQ(ierr);
    if (k==1) B = NULL;
    else {
      ierr = STGetMatrix(st,1,B);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@C
   EPSSetDeflationSpace - Specify a basis of vectors that constitute the deflation
   space.

   Collective on eps

   Input Parameter:
+  eps - the eigenproblem solver context
.  n   - number of vectors
-  v   - set of basis vectors of the deflation space

   Notes:
   When a deflation space is given, the eigensolver seeks the eigensolution
   in the restriction of the problem to the orthogonal complement of this
   space. This can be used for instance in the case that an invariant
   subspace is known beforehand (such as the nullspace of the matrix).

   These vectors do not persist from one EPSSolve() call to the other, so the
   deflation space should be set every time.

   The vectors do not need to be mutually orthonormal, since they are explicitly
   orthonormalized internally.

   Level: intermediate

.seealso: EPSSetInitialSpace()
@*/
PetscErrorCode EPSSetDeflationSpace(EPS eps,PetscInt n,Vec v[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,n,2);
  if (n<0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument n out of range");
  if (n>0) {
    PetscValidPointer(v,3);
    PetscValidHeaderSpecific(*v,VEC_CLASSID,3);
  }
  ierr = SlepcBasisReference_Private(n,v,&eps->nds,&eps->defl);CHKERRQ(ierr);
  if (n>0) eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@C
   EPSSetInitialSpace - Specify a basis of vectors that constitute the initial
   space, that is, the subspace from which the solver starts to iterate.

   Collective on eps

   Input Parameter:
+  eps - the eigenproblem solver context
.  n   - number of vectors
-  is  - set of basis vectors of the initial space

   Notes:
   Some solvers start to iterate on a single vector (initial vector). In that case,
   the other vectors are ignored.

   These vectors do not persist from one EPSSolve() call to the other, so the
   initial space should be set every time.

   The vectors do not need to be mutually orthonormal, since they are explicitly
   orthonormalized internally.

   Common usage of this function is when the user can provide a rough approximation
   of the wanted eigenspace. Then, convergence may be faster.

   Level: intermediate

.seealso: EPSSetLeftInitialSpace(), EPSSetDeflationSpace()
@*/
PetscErrorCode EPSSetInitialSpace(EPS eps,PetscInt n,Vec is[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,n,2);
  if (n<0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative");
  if (n>0) {
    PetscValidPointer(is,3);
    PetscValidHeaderSpecific(*is,VEC_CLASSID,3);
  }
  ierr = SlepcBasisReference_Private(n,is,&eps->nini,&eps->IS);CHKERRQ(ierr);
  if (n>0) eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@C
   EPSSetLeftInitialSpace - Specify a basis of vectors that constitute the left
   initial space, used by two-sided solvers to start the left subspace.

   Collective on eps

   Input Parameter:
+  eps - the eigenproblem solver context
.  n   - number of vectors
-  isl - set of basis vectors of the left initial space

   Notes:
   Left initial vectors are used to initiate the left search space in two-sided
   eigensolvers. Users should pass here an approximation of the left eigenspace,
   if available.

   The same comments in EPSSetInitialSpace() are applicable here.

   Level: intermediate

.seealso: EPSSetInitialSpace(), EPSSetTwoSided()
@*/
PetscErrorCode EPSSetLeftInitialSpace(EPS eps,PetscInt n,Vec isl[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,n,2);
  if (n<0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Argument n cannot be negative");
  if (n>0) {
    PetscValidPointer(isl,3);
    PetscValidHeaderSpecific(*isl,VEC_CLASSID,3);
  }
  ierr = SlepcBasisReference_Private(n,isl,&eps->ninil,&eps->ISL);CHKERRQ(ierr);
  if (n>0) eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*
  EPSSetDimensions_Default - Set reasonable values for ncv, mpd if not set
  by the user. This is called at setup.
 */
PetscErrorCode EPSSetDimensions_Default(EPS eps,PetscInt nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscErrorCode ierr;
  PetscBool      krylov;

  PetscFunctionBegin;
  if (*ncv!=PETSC_DEFAULT) { /* ncv set */
    ierr = PetscObjectTypeCompareAny((PetscObject)eps,&krylov,EPSKRYLOVSCHUR,EPSARNOLDI,EPSLANCZOS,"");CHKERRQ(ierr);
    if (krylov) {
      if (*ncv<nev+1 && !(*ncv==nev && *ncv==eps->n)) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv must be at least nev+1");
    } else {
      if (*ncv<nev) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv must be at least nev");
    }
  } else if (*mpd!=PETSC_DEFAULT) { /* mpd set */
    *ncv = PetscMin(eps->n,nev+(*mpd));
  } else { /* neither set: defaults depend on nev being small or large */
    if (nev<500) *ncv = PetscMin(eps->n,PetscMax(2*nev,nev+15));
    else {
      *mpd = 500;
      *ncv = PetscMin(eps->n,nev+(*mpd));
    }
  }
  if (*mpd==PETSC_DEFAULT) *mpd = *ncv;
  PetscFunctionReturn(0);
}

/*@
   EPSAllocateSolution - Allocate memory storage for common variables such
   as eigenvalues and eigenvectors.

   Collective on eps

   Input Parameters:
+  eps   - eigensolver context
-  extra - number of additional positions, used for methods that require a
           working basis slightly larger than ncv

   Developers Note:
   This is SLEPC_EXTERN because it may be required by user plugin EPS
   implementations.

   Level: developer
@*/
PetscErrorCode EPSAllocateSolution(EPS eps,PetscInt extra)
{
  PetscErrorCode ierr;
  PetscInt       oldsize,newc,requested;
  PetscLogDouble cnt;
  Vec            t;

  PetscFunctionBegin;
  requested = eps->ncv + extra;

  /* oldsize is zero if this is the first time setup is called */
  ierr = BVGetSizes(eps->V,NULL,NULL,&oldsize);CHKERRQ(ierr);
  newc = PetscMax(0,requested-oldsize);

  /* allocate space for eigenvalues and friends */
  if (requested != oldsize || !eps->eigr) {
    ierr = PetscFree4(eps->eigr,eps->eigi,eps->errest,eps->perm);CHKERRQ(ierr);
    ierr = PetscMalloc4(requested,&eps->eigr,requested,&eps->eigi,requested,&eps->errest,requested,&eps->perm);CHKERRQ(ierr);
    cnt = 2*newc*sizeof(PetscScalar) + 2*newc*sizeof(PetscReal) + newc*sizeof(PetscInt);
    ierr = PetscLogObjectMemory((PetscObject)eps,cnt);CHKERRQ(ierr);
  }

  /* workspace for the case of arbitrary selection */
  if (eps->arbitrary) {
    if (eps->rr) {
      ierr = PetscFree2(eps->rr,eps->ri);CHKERRQ(ierr);
    }
    ierr = PetscMalloc2(requested,&eps->rr,requested,&eps->ri);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)eps,2*newc*sizeof(PetscScalar));CHKERRQ(ierr);
  }

  /* allocate V */
  if (!eps->V) { ierr = EPSGetBV(eps,&eps->V);CHKERRQ(ierr); }
  if (!oldsize) {
    if (!((PetscObject)(eps->V))->type_name) {
      ierr = BVSetType(eps->V,BVSVEC);CHKERRQ(ierr);
    }
    ierr = STMatCreateVecsEmpty(eps->st,&t,NULL);CHKERRQ(ierr);
    ierr = BVSetSizesFromVec(eps->V,t,requested);CHKERRQ(ierr);
    ierr = VecDestroy(&t);CHKERRQ(ierr);
  } else {
    ierr = BVResize(eps->V,requested,PETSC_FALSE);CHKERRQ(ierr);
  }

  /* allocate W */
  if (eps->twosided) {
    ierr = BVDestroy(&eps->W);CHKERRQ(ierr);
    ierr = BVDuplicate(eps->V,&eps->W);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

