/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic EPS routines
*/

#include <slepc/private/epsimpl.h>      /*I "slepceps.h" I*/

PetscFunctionList EPSList = 0;
PetscBool         EPSRegisterAllCalled = PETSC_FALSE;
PetscClassId      EPS_CLASSID = 0;
PetscLogEvent     EPS_SetUp = 0,EPS_Solve = 0;

/*@
   EPSCreate - Creates the default EPS context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  eps - location to put the EPS context

   Note:
   The default EPS type is EPSKRYLOVSCHUR

   Level: beginner

.seealso: EPSSetUp(), EPSSolve(), EPSDestroy(), EPS
@*/
PetscErrorCode EPSCreate(MPI_Comm comm,EPS *outeps)
{
  PetscErrorCode ierr;
  EPS            eps;

  PetscFunctionBegin;
  PetscValidPointer(outeps,2);
  *outeps = 0;
  ierr = EPSInitializePackage();CHKERRQ(ierr);
  ierr = SlepcHeaderCreate(eps,EPS_CLASSID,"EPS","Eigenvalue Problem Solver","EPS",comm,EPSDestroy,EPSView);CHKERRQ(ierr);

  eps->max_it          = PETSC_DEFAULT;
  eps->nev             = 1;
  eps->ncv             = PETSC_DEFAULT;
  eps->mpd             = PETSC_DEFAULT;
  eps->nini            = 0;
  eps->nds             = 0;
  eps->target          = 0.0;
  eps->tol             = PETSC_DEFAULT;
  eps->conv            = EPS_CONV_REL;
  eps->stop            = EPS_STOP_BASIC;
  eps->which           = (EPSWhich)0;
  eps->inta            = 0.0;
  eps->intb            = 0.0;
  eps->problem_type    = (EPSProblemType)0;
  eps->extraction      = EPS_RITZ;
  eps->balance         = EPS_BALANCE_NONE;
  eps->balance_its     = 5;
  eps->balance_cutoff  = 1e-8;
  eps->trueres         = PETSC_FALSE;
  eps->trackall        = PETSC_FALSE;
  eps->purify          = PETSC_TRUE;
  eps->twosided        = PETSC_FALSE;

  eps->converged       = EPSConvergedRelative;
  eps->convergeduser   = NULL;
  eps->convergeddestroy= NULL;
  eps->stopping        = EPSStoppingBasic;
  eps->stoppinguser    = NULL;
  eps->stoppingdestroy = NULL;
  eps->arbitrary       = NULL;
  eps->convergedctx    = NULL;
  eps->stoppingctx     = NULL;
  eps->arbitraryctx    = NULL;
  eps->numbermonitors  = 0;

  eps->st              = NULL;
  eps->ds              = NULL;
  eps->dsts            = NULL;
  eps->V               = NULL;
  eps->W               = NULL;
  eps->rg              = NULL;
  eps->D               = NULL;
  eps->IS              = NULL;
  eps->ISL             = NULL;
  eps->defl            = NULL;
  eps->eigr            = NULL;
  eps->eigi            = NULL;
  eps->errest          = NULL;
  eps->rr              = NULL;
  eps->ri              = NULL;
  eps->perm            = NULL;
  eps->nwork           = 0;
  eps->work            = NULL;
  eps->data            = NULL;

  eps->state           = EPS_STATE_INITIAL;
  eps->categ           = EPS_CATEGORY_KRYLOV;
  eps->nconv           = 0;
  eps->its             = 0;
  eps->nloc            = 0;
  eps->nrma            = 0.0;
  eps->nrmb            = 0.0;
  eps->useds           = PETSC_FALSE;
  eps->hasts           = PETSC_FALSE;
  eps->isgeneralized   = PETSC_FALSE;
  eps->ispositive      = PETSC_FALSE;
  eps->ishermitian     = PETSC_FALSE;
  eps->reason          = EPS_CONVERGED_ITERATING;

  ierr = PetscNewLog(eps,&eps->sc);CHKERRQ(ierr);
  *outeps = eps;
  PetscFunctionReturn(0);
}

/*@C
   EPSSetType - Selects the particular solver to be used in the EPS object.

   Logically Collective on eps

   Input Parameters:
+  eps  - the eigensolver context
-  type - a known method

   Options Database Key:
.  -eps_type <method> - Sets the method; use -help for a list
    of available methods

   Notes:
   See "slepc/include/slepceps.h" for available methods. The default
   is EPSKRYLOVSCHUR.

   Normally, it is best to use the EPSSetFromOptions() command and
   then set the EPS type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different available methods.
   The EPSSetType() routine is provided for those situations where it
   is necessary to set the iterative solver independently of the command
   line or options database.

   Level: intermediate

.seealso: STSetType(), EPSType
@*/
PetscErrorCode EPSSetType(EPS eps,EPSType type)
{
  PetscErrorCode ierr,(*r)(EPS);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidCharPointer(type,2);

  ierr = PetscObjectTypeCompare((PetscObject)eps,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);

  ierr = PetscFunctionListFind(EPSList,type,&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown EPS type given: %s",type);

  if (eps->ops->destroy) { ierr = (*eps->ops->destroy)(eps);CHKERRQ(ierr); }
  ierr = PetscMemzero(eps->ops,sizeof(struct _EPSOps));CHKERRQ(ierr);

  eps->state = EPS_STATE_INITIAL;
  ierr = PetscObjectChangeTypeName((PetscObject)eps,type);CHKERRQ(ierr);
  ierr = (*r)(eps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   EPSGetType - Gets the EPS type as a string from the EPS object.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  name - name of EPS method

   Level: intermediate

.seealso: EPSSetType()
@*/
PetscErrorCode EPSGetType(EPS eps,EPSType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(type,2);
  *type = ((PetscObject)eps)->type_name;
  PetscFunctionReturn(0);
}

/*@C
   EPSRegister - Adds a method to the eigenproblem solver package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined solver
-  function - routine to create the solver context

   Notes:
   EPSRegister() may be called multiple times to add several user-defined solvers.

   Sample usage:
.vb
    EPSRegister("my_solver",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     EPSSetType(eps,"my_solver")
   or at runtime via the option
$     -eps_type my_solver

   Level: advanced

.seealso: EPSRegisterAll()
@*/
PetscErrorCode EPSRegister(const char *name,PetscErrorCode (*function)(EPS))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = EPSInitializePackage();CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&EPSList,name,function);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSReset - Resets the EPS context to the initial state (prior to setup)
   and destroys any allocated Vecs and Mats.

   Collective on eps

   Input Parameter:
.  eps - eigensolver context obtained from EPSCreate()

   Note:
   This can be used when a problem of different matrix size wants to be solved.
   All options that have previously been set are preserved, so in a next use
   the solver configuration is the same, but new sizes for matrices and vectors
   are allowed.

   Level: advanced

.seealso: EPSDestroy()
@*/
PetscErrorCode EPSReset(EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (eps) PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!eps) PetscFunctionReturn(0);
  if (eps->ops->reset) { ierr = (eps->ops->reset)(eps);CHKERRQ(ierr); }
  if (eps->st) { ierr = STReset(eps->st);CHKERRQ(ierr); }
  ierr = VecDestroy(&eps->D);CHKERRQ(ierr);
  ierr = BVDestroy(&eps->V);CHKERRQ(ierr);
  ierr = BVDestroy(&eps->W);CHKERRQ(ierr);
  ierr = VecDestroyVecs(eps->nwork,&eps->work);CHKERRQ(ierr);
  eps->nwork = 0;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSDestroy - Destroys the EPS context.

   Collective on eps

   Input Parameter:
.  eps - eigensolver context obtained from EPSCreate()

   Level: beginner

.seealso: EPSCreate(), EPSSetUp(), EPSSolve()
@*/
PetscErrorCode EPSDestroy(EPS *eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*eps) PetscFunctionReturn(0);
  PetscValidHeaderSpecific(*eps,EPS_CLASSID,1);
  if (--((PetscObject)(*eps))->refct > 0) { *eps = 0; PetscFunctionReturn(0); }
  ierr = EPSReset(*eps);CHKERRQ(ierr);
  if ((*eps)->ops->destroy) { ierr = (*(*eps)->ops->destroy)(*eps);CHKERRQ(ierr); }
  if ((*eps)->eigr) {
    ierr = PetscFree4((*eps)->eigr,(*eps)->eigi,(*eps)->errest,(*eps)->perm);CHKERRQ(ierr);
  }
  if ((*eps)->rr) {
    ierr = PetscFree2((*eps)->rr,(*eps)->ri);CHKERRQ(ierr);
  }
  ierr = STDestroy(&(*eps)->st);CHKERRQ(ierr);
  ierr = RGDestroy(&(*eps)->rg);CHKERRQ(ierr);
  ierr = DSDestroy(&(*eps)->ds);CHKERRQ(ierr);
  ierr = DSDestroy(&(*eps)->dsts);CHKERRQ(ierr);
  ierr = PetscFree((*eps)->sc);CHKERRQ(ierr);
  /* just in case the initial vectors have not been used */
  ierr = SlepcBasisDestroy_Private(&(*eps)->nds,&(*eps)->defl);CHKERRQ(ierr);
  ierr = SlepcBasisDestroy_Private(&(*eps)->nini,&(*eps)->IS);CHKERRQ(ierr);
  ierr = SlepcBasisDestroy_Private(&(*eps)->ninil,&(*eps)->ISL);CHKERRQ(ierr);
  if ((*eps)->convergeddestroy) {
    ierr = (*(*eps)->convergeddestroy)((*eps)->convergedctx);CHKERRQ(ierr);
  }
  ierr = EPSMonitorCancel(*eps);CHKERRQ(ierr);
  ierr = PetscHeaderDestroy(eps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSSetTarget - Sets the value of the target.

   Logically Collective on eps

   Input Parameters:
+  eps    - eigensolver context
-  target - the value of the target

   Options Database Key:
.  -eps_target <scalar> - the value of the target

   Notes:
   The target is a scalar value used to determine the portion of the spectrum
   of interest. It is used in combination with EPSSetWhichEigenpairs().

   In the case of complex scalars, a complex value can be provided in the
   command line with [+/-][realnumber][+/-]realnumberi with no spaces, e.g.
   -eps_target 1.0+2.0i

   Level: intermediate

.seealso: EPSGetTarget(), EPSSetWhichEigenpairs()
@*/
PetscErrorCode EPSSetTarget(EPS eps,PetscScalar target)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveScalar(eps,target,2);
  eps->target = target;
  if (!eps->st) { ierr = EPSGetST(eps,&eps->st);CHKERRQ(ierr); }
  ierr = STSetDefaultShift(eps->st,target);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGetTarget - Gets the value of the target.

   Not Collective

   Input Parameter:
.  eps - eigensolver context

   Output Parameter:
.  target - the value of the target

   Note:
   If the target was not set by the user, then zero is returned.

   Level: intermediate

.seealso: EPSSetTarget()
@*/
PetscErrorCode EPSGetTarget(EPS eps,PetscScalar* target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidScalarPointer(target,2);
  *target = eps->target;
  PetscFunctionReturn(0);
}

/*@
   EPSSetInterval - Defines the computational interval for spectrum slicing.

   Logically Collective on eps

   Input Parameters:
+  eps  - eigensolver context
.  inta - left end of the interval
-  intb - right end of the interval

   Options Database Key:
.  -eps_interval <a,b> - set [a,b] as the interval of interest

   Notes:
   Spectrum slicing is a technique employed for computing all eigenvalues of
   symmetric eigenproblems in a given interval. This function provides the
   interval to be considered. It must be used in combination with EPS_ALL, see
   EPSSetWhichEigenpairs().

   In the command-line option, two values must be provided. For an open interval,
   one can give an infinite, e.g., -eps_interval 1.0,inf or -eps_interval -inf,1.0.
   An open interval in the programmatic interface can be specified with
   PETSC_MAX_REAL and -PETSC_MAX_REAL.

   Level: intermediate

.seealso: EPSGetInterval(), EPSSetWhichEigenpairs()
@*/
PetscErrorCode EPSSetInterval(EPS eps,PetscReal inta,PetscReal intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,inta,2);
  PetscValidLogicalCollectiveReal(eps,intb,3);
  if (inta>intb) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be inta<intb");
  if (eps->inta != inta || eps->intb != intb) {
    eps->inta = inta;
    eps->intb = intb;
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSGetInterval - Gets the computational interval for spectrum slicing.

   Not Collective

   Input Parameter:
.  eps - eigensolver context

   Output Parameters:
+  inta - left end of the interval
-  intb - right end of the interval

   Level: intermediate

   Note:
   If the interval was not set by the user, then zeros are returned.

.seealso: EPSSetInterval()
@*/
PetscErrorCode EPSGetInterval(EPS eps,PetscReal* inta,PetscReal* intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (inta) *inta = eps->inta;
  if (intb) *intb = eps->intb;
  PetscFunctionReturn(0);
}

/*@
   EPSSetST - Associates a spectral transformation object to the eigensolver.

   Collective on eps

   Input Parameters:
+  eps - eigensolver context obtained from EPSCreate()
-  st   - the spectral transformation object

   Note:
   Use EPSGetST() to retrieve the spectral transformation context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: EPSGetST()
@*/
PetscErrorCode EPSSetST(EPS eps,ST st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(st,ST_CLASSID,2);
  PetscCheckSameComm(eps,1,st,2);
  ierr = PetscObjectReference((PetscObject)st);CHKERRQ(ierr);
  ierr = STDestroy(&eps->st);CHKERRQ(ierr);
  eps->st = st;
  ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)eps->st);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGetST - Obtain the spectral transformation (ST) object associated
   to the eigensolver object.

   Not Collective

   Input Parameters:
.  eps - eigensolver context obtained from EPSCreate()

   Output Parameter:
.  st - spectral transformation context

   Level: intermediate

.seealso: EPSSetST()
@*/
PetscErrorCode EPSGetST(EPS eps,ST *st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(st,2);
  if (!eps->st) {
    ierr = STCreate(PetscObjectComm((PetscObject)eps),&eps->st);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)eps->st,(PetscObject)eps,0);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)eps->st);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)eps->st,((PetscObject)eps)->options);CHKERRQ(ierr);
  }
  *st = eps->st;
  PetscFunctionReturn(0);
}

/*@
   EPSSetBV - Associates a basis vectors object to the eigensolver.

   Collective on eps

   Input Parameters:
+  eps - eigensolver context obtained from EPSCreate()
-  V   - the basis vectors object

   Level: advanced

.seealso: EPSGetBV()
@*/
PetscErrorCode EPSSetBV(EPS eps,BV V)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(V,BV_CLASSID,2);
  PetscCheckSameComm(eps,1,V,2);
  ierr = PetscObjectReference((PetscObject)V);CHKERRQ(ierr);
  ierr = BVDestroy(&eps->V);CHKERRQ(ierr);
  eps->V = V;
  ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)eps->V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGetBV - Obtain the basis vectors object associated to the eigensolver object.

   Not Collective

   Input Parameters:
.  eps - eigensolver context obtained from EPSCreate()

   Output Parameter:
.  V - basis vectors context

   Level: advanced

.seealso: EPSSetBV()
@*/
PetscErrorCode EPSGetBV(EPS eps,BV *V)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(V,2);
  if (!eps->V) {
    ierr = BVCreate(PetscObjectComm((PetscObject)eps),&eps->V);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)eps->V,(PetscObject)eps,0);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)eps->V);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)eps->V,((PetscObject)eps)->options);CHKERRQ(ierr);
  }
  *V = eps->V;
  PetscFunctionReturn(0);
}

/*@
   EPSSetRG - Associates a region object to the eigensolver.

   Collective on eps

   Input Parameters:
+  eps - eigensolver context obtained from EPSCreate()
-  rg  - the region object

   Note:
   Use EPSGetRG() to retrieve the region context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: EPSGetRG()
@*/
PetscErrorCode EPSSetRG(EPS eps,RG rg)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(rg,RG_CLASSID,2);
  PetscCheckSameComm(eps,1,rg,2);
  ierr = PetscObjectReference((PetscObject)rg);CHKERRQ(ierr);
  ierr = RGDestroy(&eps->rg);CHKERRQ(ierr);
  eps->rg = rg;
  ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)eps->rg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGetRG - Obtain the region object associated to the eigensolver.

   Not Collective

   Input Parameters:
.  eps - eigensolver context obtained from EPSCreate()

   Output Parameter:
.  rg - region context

   Level: advanced

.seealso: EPSSetRG()
@*/
PetscErrorCode EPSGetRG(EPS eps,RG *rg)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(rg,2);
  if (!eps->rg) {
    ierr = RGCreate(PetscObjectComm((PetscObject)eps),&eps->rg);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)eps->rg,(PetscObject)eps,0);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)eps->rg);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)eps->rg,((PetscObject)eps)->options);CHKERRQ(ierr);
  }
  *rg = eps->rg;
  PetscFunctionReturn(0);
}

/*@
   EPSSetDS - Associates a direct solver object to the eigensolver.

   Collective on eps

   Input Parameters:
+  eps - eigensolver context obtained from EPSCreate()
-  ds  - the direct solver object

   Note:
   Use EPSGetDS() to retrieve the direct solver context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: EPSGetDS()
@*/
PetscErrorCode EPSSetDS(EPS eps,DS ds)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(ds,DS_CLASSID,2);
  PetscCheckSameComm(eps,1,ds,2);
  ierr = PetscObjectReference((PetscObject)ds);CHKERRQ(ierr);
  ierr = DSDestroy(&eps->ds);CHKERRQ(ierr);
  eps->ds = ds;
  ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)eps->ds);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGetDS - Obtain the direct solver object associated to the eigensolver object.

   Not Collective

   Input Parameters:
.  eps - eigensolver context obtained from EPSCreate()

   Output Parameter:
.  ds - direct solver context

   Level: advanced

.seealso: EPSSetDS()
@*/
PetscErrorCode EPSGetDS(EPS eps,DS *ds)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(ds,2);
  if (!eps->ds) {
    ierr = DSCreate(PetscObjectComm((PetscObject)eps),&eps->ds);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)eps->ds,(PetscObject)eps,0);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)eps->ds);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)eps->ds,((PetscObject)eps)->options);CHKERRQ(ierr);
  }
  *ds = eps->ds;
  PetscFunctionReturn(0);
}

/*@
   EPSIsGeneralized - Ask if the EPS object corresponds to a generalized
   eigenvalue problem.

   Not collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  is - the answer

   Level: intermediate

.seealso: EPSIsHermitian(), EPSIsPositive()
@*/
PetscErrorCode EPSIsGeneralized(EPS eps,PetscBool* is)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(is,2);
  *is = eps->isgeneralized;
  PetscFunctionReturn(0);
}

/*@
   EPSIsHermitian - Ask if the EPS object corresponds to a Hermitian
   eigenvalue problem.

   Not collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  is - the answer

   Level: intermediate

.seealso: EPSIsGeneralized(), EPSIsPositive()
@*/
PetscErrorCode EPSIsHermitian(EPS eps,PetscBool* is)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(is,2);
  *is = eps->ishermitian;
  PetscFunctionReturn(0);
}

/*@
   EPSIsPositive - Ask if the EPS object corresponds to an eigenvalue
   problem type that requires a positive (semi-) definite matrix B.

   Not collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  is - the answer

   Level: intermediate

.seealso: EPSIsGeneralized(), EPSIsHermitian()
@*/
PetscErrorCode EPSIsPositive(EPS eps,PetscBool* is)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(is,2);
  *is = eps->ispositive;
  PetscFunctionReturn(0);
}

