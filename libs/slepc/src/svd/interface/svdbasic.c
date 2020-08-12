/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic SVD routines
*/

#include <slepc/private/svdimpl.h>      /*I "slepcsvd.h" I*/

PetscFunctionList SVDList = 0;
PetscBool         SVDRegisterAllCalled = PETSC_FALSE;
PetscClassId      SVD_CLASSID = 0;
PetscLogEvent     SVD_SetUp = 0,SVD_Solve = 0;

/*@
   SVDCreate - Creates the default SVD context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  svd - location to put the SVD context

   Note:
   The default SVD type is SVDCROSS

   Level: beginner

.seealso: SVDSetUp(), SVDSolve(), SVDDestroy(), SVD
@*/
PetscErrorCode SVDCreate(MPI_Comm comm,SVD *outsvd)
{
  PetscErrorCode ierr;
  SVD            svd;

  PetscFunctionBegin;
  PetscValidPointer(outsvd,2);
  *outsvd = 0;
  ierr = SVDInitializePackage();CHKERRQ(ierr);
  ierr = SlepcHeaderCreate(svd,SVD_CLASSID,"SVD","Singular Value Decomposition","SVD",comm,SVDDestroy,SVDView);CHKERRQ(ierr);

  svd->OP               = NULL;
  svd->max_it           = PETSC_DEFAULT;
  svd->nsv              = 1;
  svd->ncv              = PETSC_DEFAULT;
  svd->mpd              = PETSC_DEFAULT;
  svd->nini             = 0;
  svd->ninil            = 0;
  svd->tol              = PETSC_DEFAULT;
  svd->conv             = SVD_CONV_REL;
  svd->stop             = SVD_STOP_BASIC;
  svd->which            = SVD_LARGEST;
  svd->impltrans        = PETSC_FALSE;
  svd->trackall         = PETSC_FALSE;

  svd->converged        = SVDConvergedRelative;
  svd->convergeduser    = NULL;
  svd->convergeddestroy = NULL;
  svd->stopping         = SVDStoppingBasic;
  svd->stoppinguser     = NULL;
  svd->stoppingdestroy  = NULL;
  svd->convergedctx     = NULL;
  svd->stoppingctx      = NULL;
  svd->numbermonitors   = 0;

  svd->ds               = NULL;
  svd->U                = NULL;
  svd->V                = NULL;
  svd->A                = NULL;
  svd->AT               = NULL;
  svd->IS               = NULL;
  svd->ISL              = NULL;
  svd->sigma            = NULL;
  svd->perm             = NULL;
  svd->errest           = NULL;
  svd->data             = NULL;

  svd->state            = SVD_STATE_INITIAL;
  svd->nconv            = 0;
  svd->its              = 0;
  svd->leftbasis        = PETSC_FALSE;
  svd->reason           = SVD_CONVERGED_ITERATING;

  ierr = PetscNewLog(svd,&svd->sc);CHKERRQ(ierr);
  *outsvd = svd;
  PetscFunctionReturn(0);
}

/*@
   SVDReset - Resets the SVD context to the initial state (prior to setup)
   and destroys any allocated Vecs and Mats.

   Collective on svd

   Input Parameter:
.  svd - singular value solver context obtained from SVDCreate()

   Level: advanced

.seealso: SVDDestroy()
@*/
PetscErrorCode SVDReset(SVD svd)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (svd) PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!svd) PetscFunctionReturn(0);
  if (svd->ops->reset) { ierr = (svd->ops->reset)(svd);CHKERRQ(ierr); }
  ierr = MatDestroy(&svd->OP);CHKERRQ(ierr);
  ierr = MatDestroy(&svd->A);CHKERRQ(ierr);
  ierr = MatDestroy(&svd->AT);CHKERRQ(ierr);
  ierr = BVDestroy(&svd->U);CHKERRQ(ierr);
  ierr = BVDestroy(&svd->V);CHKERRQ(ierr);
  svd->state = SVD_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   SVDDestroy - Destroys the SVD context.

   Collective on svd

   Input Parameter:
.  svd - singular value solver context obtained from SVDCreate()

   Level: beginner

.seealso: SVDCreate(), SVDSetUp(), SVDSolve()
@*/
PetscErrorCode SVDDestroy(SVD *svd)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*svd) PetscFunctionReturn(0);
  PetscValidHeaderSpecific(*svd,SVD_CLASSID,1);
  if (--((PetscObject)(*svd))->refct > 0) { *svd = 0; PetscFunctionReturn(0); }
  ierr = SVDReset(*svd);CHKERRQ(ierr);
  if ((*svd)->ops->destroy) { ierr = (*(*svd)->ops->destroy)(*svd);CHKERRQ(ierr); }
  if ((*svd)->sigma) {
    ierr = PetscFree3((*svd)->sigma,(*svd)->perm,(*svd)->errest);CHKERRQ(ierr);
  }
  ierr = DSDestroy(&(*svd)->ds);CHKERRQ(ierr);
  ierr = PetscFree((*svd)->sc);CHKERRQ(ierr);
  /* just in case the initial vectors have not been used */
  ierr = SlepcBasisDestroy_Private(&(*svd)->nini,&(*svd)->IS);CHKERRQ(ierr);
  ierr = SlepcBasisDestroy_Private(&(*svd)->ninil,&(*svd)->ISL);CHKERRQ(ierr);
  ierr = SVDMonitorCancel(*svd);CHKERRQ(ierr);
  ierr = PetscHeaderDestroy(svd);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   SVDSetType - Selects the particular solver to be used in the SVD object.

   Logically Collective on svd

   Input Parameters:
+  svd      - the singular value solver context
-  type     - a known method

   Options Database Key:
.  -svd_type <method> - Sets the method; use -help for a list
    of available methods

   Notes:
   See "slepc/include/slepcsvd.h" for available methods. The default
   is SVDCROSS.

   Normally, it is best to use the SVDSetFromOptions() command and
   then set the SVD type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different available methods.
   The SVDSetType() routine is provided for those situations where it
   is necessary to set the iterative solver independently of the command
   line or options database.

   Level: intermediate

.seealso: SVDType
@*/
PetscErrorCode SVDSetType(SVD svd,SVDType type)
{
  PetscErrorCode ierr,(*r)(SVD);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidCharPointer(type,2);

  ierr = PetscObjectTypeCompare((PetscObject)svd,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);

  ierr = PetscFunctionListFind(SVDList,type,&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown SVD type given: %s",type);

  if (svd->ops->destroy) { ierr = (*svd->ops->destroy)(svd);CHKERRQ(ierr); }
  ierr = PetscMemzero(svd->ops,sizeof(struct _SVDOps));CHKERRQ(ierr);

  svd->state = SVD_STATE_INITIAL;
  ierr = PetscObjectChangeTypeName((PetscObject)svd,type);CHKERRQ(ierr);
  ierr = (*r)(svd);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   SVDGetType - Gets the SVD type as a string from the SVD object.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  name - name of SVD method

   Level: intermediate

.seealso: SVDSetType()
@*/
PetscErrorCode SVDGetType(SVD svd,SVDType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(type,2);
  *type = ((PetscObject)svd)->type_name;
  PetscFunctionReturn(0);
}

/*@C
   SVDRegister - Adds a method to the singular value solver package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined solver
-  function - routine to create the solver context

   Notes:
   SVDRegister() may be called multiple times to add several user-defined solvers.

   Sample usage:
.vb
    SVDRegister("my_solver",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     SVDSetType(svd,"my_solver")
   or at runtime via the option
$     -svd_type my_solver

   Level: advanced

.seealso: SVDRegisterAll()
@*/
PetscErrorCode SVDRegister(const char *name,PetscErrorCode (*function)(SVD))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SVDInitializePackage();CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&SVDList,name,function);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   SVDSetBV - Associates basis vectors objects to the singular value solver.

   Collective on svd

   Input Parameters:
+  svd - singular value solver context obtained from SVDCreate()
.  V   - the basis vectors object for right singular vectors
-  U   - the basis vectors object for left singular vectors

   Note:
   Use SVDGetBV() to retrieve the basis vectors contexts (for example,
   to free them at the end of the computations).

   Level: advanced

.seealso: SVDGetBV()
@*/
PetscErrorCode SVDSetBV(SVD svd,BV V,BV U)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (V) {
    PetscValidHeaderSpecific(V,BV_CLASSID,2);
    PetscCheckSameComm(svd,1,V,2);
    ierr = PetscObjectReference((PetscObject)V);CHKERRQ(ierr);
    ierr = BVDestroy(&svd->V);CHKERRQ(ierr);
    svd->V = V;
    ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)svd->V);CHKERRQ(ierr);
  }
  if (U) {
    PetscValidHeaderSpecific(U,BV_CLASSID,3);
    PetscCheckSameComm(svd,1,U,3);
    ierr = PetscObjectReference((PetscObject)U);CHKERRQ(ierr);
    ierr = BVDestroy(&svd->U);CHKERRQ(ierr);
    svd->U = U;
    ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)svd->U);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   SVDGetBV - Obtain the basis vectors objects associated to the singular
   value solver object.

   Not Collective

   Input Parameters:
.  svd - singular value solver context obtained from SVDCreate()

   Output Parameter:
+  V - basis vectors context for right singular vectors
-  U - basis vectors context for left singular vectors

   Level: advanced

.seealso: SVDSetBV()
@*/
PetscErrorCode SVDGetBV(SVD svd,BV *V,BV *U)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (V) {
    if (!svd->V) {
      ierr = BVCreate(PetscObjectComm((PetscObject)svd),&svd->V);CHKERRQ(ierr);
      ierr = PetscObjectIncrementTabLevel((PetscObject)svd->V,(PetscObject)svd,0);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)svd->V);CHKERRQ(ierr);
      ierr = PetscObjectSetOptions((PetscObject)svd->V,((PetscObject)svd)->options);CHKERRQ(ierr);
    }
    *V = svd->V;
  }
  if (U) {
    if (!svd->U) {
      ierr = BVCreate(PetscObjectComm((PetscObject)svd),&svd->U);CHKERRQ(ierr);
      ierr = PetscObjectIncrementTabLevel((PetscObject)svd->U,(PetscObject)svd,0);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)svd->U);CHKERRQ(ierr);
      ierr = PetscObjectSetOptions((PetscObject)svd->U,((PetscObject)svd)->options);CHKERRQ(ierr);
    }
    *U = svd->U;
  }
  PetscFunctionReturn(0);
}

/*@
   SVDSetDS - Associates a direct solver object to the singular value solver.

   Collective on svd

   Input Parameters:
+  svd - singular value solver context obtained from SVDCreate()
-  ds  - the direct solver object

   Note:
   Use SVDGetDS() to retrieve the direct solver context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: SVDGetDS()
@*/
PetscErrorCode SVDSetDS(SVD svd,DS ds)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(ds,DS_CLASSID,2);
  PetscCheckSameComm(svd,1,ds,2);
  ierr = PetscObjectReference((PetscObject)ds);CHKERRQ(ierr);
  ierr = DSDestroy(&svd->ds);CHKERRQ(ierr);
  svd->ds = ds;
  ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)svd->ds);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   SVDGetDS - Obtain the direct solver object associated to the singular value
   solver object.

   Not Collective

   Input Parameters:
.  svd - singular value solver context obtained from SVDCreate()

   Output Parameter:
.  ds - direct solver context

   Level: advanced

.seealso: SVDSetDS()
@*/
PetscErrorCode SVDGetDS(SVD svd,DS *ds)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(ds,2);
  if (!svd->ds) {
    ierr = DSCreate(PetscObjectComm((PetscObject)svd),&svd->ds);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)svd->ds,(PetscObject)svd,0);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)svd->ds);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)svd->ds,((PetscObject)svd)->options);CHKERRQ(ierr);
  }
  *ds = svd->ds;
  PetscFunctionReturn(0);
}

