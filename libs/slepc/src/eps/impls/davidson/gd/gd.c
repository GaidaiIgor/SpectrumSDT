/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "gd"

   Method: Generalized Davidson

   Algorithm:

       Generalized Davidson with various subspace extraction and
       restart techniques.

   References:

       [1] E.R. Davidson, "Super-matrix methods", Comput. Phys. Commun.
           53(2):49-60, 1989.

       [2] E. Romero and J.E. Roman, "A parallel implementation of
           Davidson methods for large-scale eigenvalue problems in
           SLEPc", ACM Trans. Math. Software 40(2), Article 13, 2014.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include <../src/eps/impls/davidson/davidson.h>

PetscErrorCode EPSSetFromOptions_GD(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      flg,flg2,op,orth;
  PetscInt       opi,opi0;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS Generalized Davidson (GD) Options");CHKERRQ(ierr);

    ierr = EPSGDGetKrylovStart(eps,&op);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-eps_gd_krylov_start","Start the search subspace with a Krylov basis","EPSGDSetKrylovStart",op,&op,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSGDSetKrylovStart(eps,op);CHKERRQ(ierr); }

    ierr = EPSGDGetBOrth(eps,&orth);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-eps_gd_borth","Use B-orthogonalization in the search subspace","EPSGDSetBOrth",op,&op,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSGDSetBOrth(eps,op);CHKERRQ(ierr); }

    ierr = EPSGDGetBlockSize(eps,&opi);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_gd_blocksize","Number of vectors to add to the search subspace","EPSGDSetBlockSize",opi,&opi,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSGDSetBlockSize(eps,opi);CHKERRQ(ierr); }

    ierr = EPSGDGetRestart(eps,&opi,&opi0);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_gd_minv","Size of the search subspace after restarting","EPSGDSetRestart",opi,&opi,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_gd_plusk","Number of eigenvectors saved from the previous iteration when restarting","EPSGDSetRestart",opi0,&opi0,&flg2);CHKERRQ(ierr);
    if (flg || flg2) { ierr = EPSGDSetRestart(eps,opi,opi0);CHKERRQ(ierr); }

    ierr = EPSGDGetInitialSize(eps,&opi);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_gd_initial_size","Initial size of the search subspace","EPSGDSetInitialSize",opi,&opi,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSGDSetInitialSize(eps,opi);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-eps_gd_double_expansion","Use the doble-expansion variant of GD","EPSGDSetDoubleExpansion",PETSC_FALSE,&op,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSGDSetDoubleExpansion(eps,op);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetUp_GD(EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      t;
  KSP            ksp;

  PetscFunctionBegin;
  /* Setup common for all davidson solvers */
  ierr = EPSSetUp_XD(eps);CHKERRQ(ierr);

  /* Check some constraints */
  ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&t);CHKERRQ(ierr);
  if (!t) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"EPSGD only works with KSPPREONLY");
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_GD(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isascii,opb;
  PetscInt       opi,opi0;
  PetscBool      borth;
  EPS_DAVIDSON   *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    if (data->doubleexp) {
      ierr = PetscViewerASCIIPrintf(viewer,"  using double expansion variant (GD2)\n");CHKERRQ(ierr);
    }
    ierr = EPSXDGetBOrth_XD(eps,&borth);CHKERRQ(ierr);
    if (borth) {
      ierr = PetscViewerASCIIPrintf(viewer,"  search subspace is B-orthogonalized\n");CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"  search subspace is orthogonalized\n");CHKERRQ(ierr);
    }
    ierr = EPSXDGetBlockSize_XD(eps,&opi);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  block size=%D\n",opi);CHKERRQ(ierr);
    ierr = EPSXDGetKrylovStart_XD(eps,&opb);CHKERRQ(ierr);
    if (!opb) {
      ierr = PetscViewerASCIIPrintf(viewer,"  type of the initial subspace: non-Krylov\n");CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"  type of the initial subspace: Krylov\n");CHKERRQ(ierr);
    }
    ierr = EPSXDGetRestart_XD(eps,&opi,&opi0);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  size of the subspace after restarting: %D\n",opi);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  number of vectors after restarting from the previous iteration: %D\n",opi0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_GD(EPS eps)
{
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetKrylovStart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetKrylovStart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetBOrth_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetBOrth_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetInitialSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetInitialSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetDoubleExpansion_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetDoubleExpansion_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDSetKrylovStart - Activates or deactivates starting the searching
   subspace with a Krylov basis.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  krylovstart - boolean flag

   Options Database Key:
.  -eps_gd_krylov_start - Activates starting the searching subspace with a
    Krylov basis

   Level: advanced

.seealso: EPSGDGetKrylovStart()
@*/
PetscErrorCode EPSGDSetKrylovStart(EPS eps,PetscBool krylovstart)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,krylovstart,2);
  ierr = PetscTryMethod(eps,"EPSGDSetKrylovStart_C",(EPS,PetscBool),(eps,krylovstart));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDGetKrylovStart - Returns a flag indicating if the search subspace is started with a
   Krylov basis.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  krylovstart - boolean flag indicating if the search subspace is started
   with a Krylov basis

   Level: advanced

.seealso: EPSGDGetKrylovStart()
@*/
PetscErrorCode EPSGDGetKrylovStart(EPS eps,PetscBool *krylovstart)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(krylovstart,2);
  ierr = PetscUseMethod(eps,"EPSGDGetKrylovStart_C",(EPS,PetscBool*),(eps,krylovstart));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDSetBlockSize - Sets the number of vectors to be added to the searching space
   in every iteration.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  blocksize - number of vectors added to the search space in every iteration

   Options Database Key:
.  -eps_gd_blocksize - number of vectors added to the search space in every iteration

   Level: advanced

.seealso: EPSGDSetKrylovStart()
@*/
PetscErrorCode EPSGDSetBlockSize(EPS eps,PetscInt blocksize)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,blocksize,2);
  ierr = PetscTryMethod(eps,"EPSGDSetBlockSize_C",(EPS,PetscInt),(eps,blocksize));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDGetBlockSize - Returns the number of vectors to be added to the searching space
   in every iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  blocksize - number of vectors added to the search space in every iteration

   Level: advanced

.seealso: EPSGDSetBlockSize()
@*/
PetscErrorCode EPSGDGetBlockSize(EPS eps,PetscInt *blocksize)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(blocksize,2);
  ierr = PetscUseMethod(eps,"EPSGDGetBlockSize_C",(EPS,PetscInt*),(eps,blocksize));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDSetRestart - Sets the number of vectors of the searching space after
   restarting and the number of vectors saved from the previous iteration.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
.  minv - number of vectors of the searching subspace after restarting
-  plusk - number of vectors saved from the previous iteration

   Options Database Keys:
+  -eps_gd_minv - number of vectors of the searching subspace after restarting
-  -eps_gd_plusk - number of vectors saved from the previous iteration

   Level: advanced

.seealso: EPSGDSetRestart()
@*/
PetscErrorCode EPSGDSetRestart(EPS eps,PetscInt minv,PetscInt plusk)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,minv,2);
  PetscValidLogicalCollectiveInt(eps,plusk,2);
  ierr = PetscTryMethod(eps,"EPSGDSetRestart_C",(EPS,PetscInt,PetscInt),(eps,minv,plusk));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDGetRestart - Gets the number of vectors of the searching space after
   restarting and the number of vectors saved from the previous iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
+  minv - number of vectors of the searching subspace after restarting
-  plusk - number of vectors saved from the previous iteration

   Level: advanced

.seealso: EPSGDSetRestart()
@*/
PetscErrorCode EPSGDGetRestart(EPS eps,PetscInt *minv,PetscInt *plusk)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = PetscUseMethod(eps,"EPSGDGetRestart_C",(EPS,PetscInt*,PetscInt*),(eps,minv,plusk));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDSetInitialSize - Sets the initial size of the searching space.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  initialsize - number of vectors of the initial searching subspace

   Options Database Key:
.  -eps_gd_initial_size - number of vectors of the initial searching subspace

   Notes:
   If EPSGDGetKrylovStart() is PETSC_FALSE and the user provides vectors with
   EPSSetInitialSpace(), up to initialsize vectors will be used; and if the
   provided vectors are not enough, the solver completes the subspace with
   random vectors. In the case of EPSGDGetKrylovStart() being PETSC_TRUE, the solver
   gets the first vector provided by the user or, if not available, a random vector,
   and expands the Krylov basis up to initialsize vectors.

   Level: advanced

.seealso: EPSGDGetInitialSize(), EPSGDGetKrylovStart()
@*/
PetscErrorCode EPSGDSetInitialSize(EPS eps,PetscInt initialsize)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,initialsize,2);
  ierr = PetscTryMethod(eps,"EPSGDSetInitialSize_C",(EPS,PetscInt),(eps,initialsize));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDGetInitialSize - Returns the initial size of the searching space.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  initialsize - number of vectors of the initial searching subspace

   Notes:
   If EPSGDGetKrylovStart() is PETSC_FALSE and the user provides vectors with
   EPSSetInitialSpace(), up to initialsize vectors will be used; and if the
   provided vectors are not enough, the solver completes the subspace with
   random vectors. In the case of EPSGDGetKrylovStart() being PETSC_TRUE, the solver
   gets the first vector provided by the user or, if not available, a random vector,
   and expands the Krylov basis up to initialsize vectors.

   Level: advanced

.seealso: EPSGDSetInitialSize(), EPSGDGetKrylovStart()
@*/
PetscErrorCode EPSGDGetInitialSize(EPS eps,PetscInt *initialsize)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(initialsize,2);
  ierr = PetscUseMethod(eps,"EPSGDGetInitialSize_C",(EPS,PetscInt*),(eps,initialsize));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDSetBOrth - Selects the orthogonalization that will be used in the search
   subspace in case of generalized Hermitian problems.

   Logically Collective on eps

   Input Parameters:
+  eps   - the eigenproblem solver context
-  borth - whether to B-orthogonalize the search subspace

   Options Database Key:
.  -eps_gd_borth - Set the orthogonalization used in the search subspace

   Level: advanced

.seealso: EPSGDGetBOrth()
@*/
PetscErrorCode EPSGDSetBOrth(EPS eps,PetscBool borth)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,borth,2);
  ierr = PetscTryMethod(eps,"EPSGDSetBOrth_C",(EPS,PetscBool),(eps,borth));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSGDGetBOrth - Returns the orthogonalization used in the search
   subspace in case of generalized Hermitian problems.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  borth - whether to B-orthogonalize the search subspace

   Level: advanced

.seealso: EPSGDSetBOrth()
@*/
PetscErrorCode EPSGDGetBOrth(EPS eps,PetscBool *borth)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(borth,2);
  ierr = PetscUseMethod(eps,"EPSGDGetBOrth_C",(EPS,PetscBool*),(eps,borth));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSGDSetDoubleExpansion_GD(EPS eps,PetscBool doubleexp)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  data->doubleexp = doubleexp;
  PetscFunctionReturn(0);
}

/*@
   EPSGDSetDoubleExpansion - Activate a variant where the search subspace is
   expanded with K*[A*x B*x] (double expansion) instead of the classic K*r,
   where K is the preconditioner, x the selected approximate eigenvector and
   r its associated residual vector.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  doubleexp - the boolean flag

   Options Database Keys:
.  -eps_gd_double_expansion - activate the double-expansion variant of GD

   Level: advanced
@*/
PetscErrorCode EPSGDSetDoubleExpansion(EPS eps,PetscBool doubleexp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,doubleexp,2);
  ierr = PetscTryMethod(eps,"EPSGDSetDoubleExpansion_C",(EPS,PetscBool),(eps,doubleexp));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSGDGetDoubleExpansion_GD(EPS eps,PetscBool *doubleexp)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *doubleexp = data->doubleexp;
  PetscFunctionReturn(0);
}

/*@
   EPSGDGetDoubleExpansion - Gets a flag indicating whether the double
   expansion variant has been activated or not.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  doubleexp - the flag

   Level: advanced

.seealso: EPSGDSetDoubleExpansion()
@*/
PetscErrorCode EPSGDGetDoubleExpansion(EPS eps,PetscBool *doubleexp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(doubleexp,2);
  ierr = PetscUseMethod(eps,"EPSGDGetDoubleExpansion_C",(EPS,PetscBool*),(eps,doubleexp));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_GD(EPS eps)
{
  PetscErrorCode  ierr;
  EPS_DAVIDSON    *data;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&data);CHKERRQ(ierr);
  eps->data = (void*)data;

  data->blocksize   = 1;
  data->initialsize = 6;
  data->minv        = 6;
  data->plusk       = PETSC_DEFAULT;
  data->ipB         = PETSC_TRUE;
  data->fix         = 0.0;
  data->krylovstart = PETSC_FALSE;
  data->dynamic     = PETSC_FALSE;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_XD;
  eps->ops->setup          = EPSSetUp_GD;
  eps->ops->setfromoptions = EPSSetFromOptions_GD;
  eps->ops->destroy        = EPSDestroy_GD;
  eps->ops->reset          = EPSReset_XD;
  eps->ops->view           = EPSView_GD;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_XD;
  eps->ops->setdefaultst   = EPSSetDefaultST_Precond;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetKrylovStart_C",EPSXDSetKrylovStart_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetKrylovStart_C",EPSXDGetKrylovStart_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetBOrth_C",EPSXDSetBOrth_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetBOrth_C",EPSXDGetBOrth_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetBlockSize_C",EPSXDSetBlockSize_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetBlockSize_C",EPSXDGetBlockSize_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetRestart_C",EPSXDSetRestart_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetRestart_C",EPSXDGetRestart_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetInitialSize_C",EPSXDSetInitialSize_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetInitialSize_C",EPSXDGetInitialSize_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDSetDoubleExpansion_C",EPSGDSetDoubleExpansion_GD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSGDGetDoubleExpansion_C",EPSGDGetDoubleExpansion_GD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

