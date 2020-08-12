/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "jd"

   Method: Jacobi-Davidson

   Algorithm:

       Jacobi-Davidson with various subspace extraction and
       restart techniques.

   References:

       [1] G.L.G. Sleijpen and H.A. van der Vorst, "A Jacobi-Davidson
           iteration method for linear eigenvalue problems", SIAM J.
           Matrix Anal. Appl. 17(2):401-425, 1996.

       [2] E. Romero and J.E. Roman, "A parallel implementation of
           Davidson methods for large-scale eigenvalue problems in
           SLEPc", ACM Trans. Math. Software 40(2), Article 13, 2014.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include <../src/eps/impls/davidson/davidson.h>

PetscErrorCode EPSSetFromOptions_JD(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode ierr;
  PetscBool      flg,flg2,op,orth;
  PetscInt       opi,opi0;
  PetscReal      opf;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS Jacobi-Davidson (JD) Options");CHKERRQ(ierr);

    ierr = EPSJDGetKrylovStart(eps,&op);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-eps_jd_krylov_start","Start the search subspace with a Krylov basis","EPSJDSetKrylovStart",op,&op,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSJDSetKrylovStart(eps,op);CHKERRQ(ierr); }

    ierr = EPSJDGetBOrth(eps,&orth);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-eps_jd_borth","Use B-orthogonalization in the search subspace","EPSJDSetBOrth",op,&op,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSJDSetBOrth(eps,op);CHKERRQ(ierr); }

    ierr = EPSJDGetBlockSize(eps,&opi);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_jd_blocksize","Number of vectors to add to the search subspace","EPSJDSetBlockSize",opi,&opi,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSJDSetBlockSize(eps,opi);CHKERRQ(ierr); }

    ierr = EPSJDGetRestart(eps,&opi,&opi0);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_jd_minv","Size of the search subspace after restarting","EPSJDSetRestart",opi,&opi,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_jd_plusk","Number of eigenvectors saved from the previous iteration when restarting","EPSJDSetRestart",opi0,&opi0,&flg2);CHKERRQ(ierr);
    if (flg || flg2) { ierr = EPSJDSetRestart(eps,opi,opi0);CHKERRQ(ierr); }

    ierr = EPSJDGetInitialSize(eps,&opi);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-eps_jd_initial_size","Initial size of the search subspace","EPSJDSetInitialSize",opi,&opi,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSJDSetInitialSize(eps,opi);CHKERRQ(ierr); }

    ierr = EPSJDGetFix(eps,&opf);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-eps_jd_fix","Tolerance for changing the target in the correction equation","EPSJDSetFix",opf,&opf,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSJDSetFix(eps,opf);CHKERRQ(ierr); }

    ierr = EPSJDGetConstCorrectionTol(eps,&op);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-eps_jd_const_correction_tol","Disable the dynamic stopping criterion when solving the correction equation","EPSJDSetConstCorrectionTol",op,&op,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSJDSetConstCorrectionTol(eps,op);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetDefaultST_JD(EPS eps)
{
  PetscErrorCode ierr;
  KSP            ksp;

  PetscFunctionBegin;
  if (!((PetscObject)eps->st)->type_name) {
    ierr = STSetType(eps->st,STPRECOND);CHKERRQ(ierr);
    ierr = STPrecondSetKSPHasMat(eps->st,PETSC_TRUE);CHKERRQ(ierr);
  }
  ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr);
  if (!((PetscObject)ksp)->type_name) {
    ierr = KSPSetType(ksp,KSPBCGSL);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1e-4,PETSC_DEFAULT,PETSC_DEFAULT,90);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetUp_JD(EPS eps)
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
  if (t) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"EPSJD does not work with KSPPREONLY");
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_JD(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isascii,opb;
  PetscReal      opf;
  PetscInt       opi,opi0;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = EPSXDGetBOrth_XD(eps,&opb);CHKERRQ(ierr);
    if (opb) {
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

    ierr = EPSJDGetFix_JD(eps,&opf);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  threshold for changing the target in the correction equation (fix): %g\n",(double)opf);CHKERRQ(ierr);

    ierr = EPSJDGetConstCorrectionTol_JD(eps,&opb);CHKERRQ(ierr);
    if (!opb) { ierr = PetscViewerASCIIPrintf(viewer,"  using dynamic tolerance for the correction equation\n");CHKERRQ(ierr); }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_JD(EPS eps)
{
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetKrylovStart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetKrylovStart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetBlockSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetRestart_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetInitialSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetInitialSize_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetFix_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetFix_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetConstCorrectionTol_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetConstCorrectionTol_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetBOrth_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetBOrth_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDSetKrylovStart - Activates or deactivates starting the searching
   subspace with a Krylov basis.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  krylovstart - boolean flag

   Options Database Key:
.  -eps_jd_krylov_start - Activates starting the searching subspace with a
    Krylov basis

   Level: advanced

.seealso: EPSJDGetKrylovStart()
@*/
PetscErrorCode EPSJDSetKrylovStart(EPS eps,PetscBool krylovstart)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,krylovstart,2);
  ierr = PetscTryMethod(eps,"EPSJDSetKrylovStart_C",(EPS,PetscBool),(eps,krylovstart));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDGetKrylovStart - Returns a flag indicating if the searching subspace is started with a
   Krylov basis.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  krylovstart - boolean flag indicating if the searching subspace is started
   with a Krylov basis

   Level: advanced

.seealso: EPSJDGetKrylovStart()
@*/
PetscErrorCode EPSJDGetKrylovStart(EPS eps,PetscBool *krylovstart)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(krylovstart,2);
  ierr = PetscUseMethod(eps,"EPSJDGetKrylovStart_C",(EPS,PetscBool*),(eps,krylovstart));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDSetBlockSize - Sets the number of vectors to be added to the searching space
   in every iteration.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  blocksize - number of vectors added to the search space in every iteration

   Options Database Key:
.  -eps_jd_blocksize - number of vectors added to the searching space every iteration

   Level: advanced

.seealso: EPSJDSetKrylovStart()
@*/
PetscErrorCode EPSJDSetBlockSize(EPS eps,PetscInt blocksize)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,blocksize,2);
  ierr = PetscTryMethod(eps,"EPSJDSetBlockSize_C",(EPS,PetscInt),(eps,blocksize));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDGetBlockSize - Returns the number of vectors to be added to the searching space
   in every iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  blocksize - number of vectors added to the search space in every iteration

   Level: advanced

.seealso: EPSJDSetBlockSize()
@*/
PetscErrorCode EPSJDGetBlockSize(EPS eps,PetscInt *blocksize)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(blocksize,2);
  ierr = PetscUseMethod(eps,"EPSJDGetBlockSize_C",(EPS,PetscInt*),(eps,blocksize));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDSetRestart - Sets the number of vectors of the searching space after
   restarting and the number of vectors saved from the previous iteration.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
.  minv - number of vectors of the searching subspace after restarting
-  plusk - number of vectors saved from the previous iteration

   Options Database Keys:
+  -eps_jd_minv - number of vectors of the searching subspace after restarting
-  -eps_jd_plusk - number of vectors saved from the previous iteration

   Level: advanced

.seealso: EPSJDGetRestart()
@*/
PetscErrorCode EPSJDSetRestart(EPS eps,PetscInt minv,PetscInt plusk)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,minv,2);
  PetscValidLogicalCollectiveInt(eps,plusk,3);
  ierr = PetscTryMethod(eps,"EPSJDSetRestart_C",(EPS,PetscInt,PetscInt),(eps,minv,plusk));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDGetRestart - Gets the number of vectors of the searching space after
   restarting and the number of vectors saved from the previous iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
+  minv - number of vectors of the searching subspace after restarting
-  plusk - number of vectors saved from the previous iteration

   Level: advanced

.seealso: EPSJDSetRestart()
@*/
PetscErrorCode EPSJDGetRestart(EPS eps,PetscInt *minv,PetscInt *plusk)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = PetscUseMethod(eps,"EPSJDGetRestart_C",(EPS,PetscInt*,PetscInt*),(eps,minv,plusk));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDSetInitialSize - Sets the initial size of the searching space.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  initialsize - number of vectors of the initial searching subspace

   Options Database Key:
.  -eps_jd_initial_size - number of vectors of the initial searching subspace

   Notes:
   If EPSJDGetKrylovStart() is PETSC_FALSE and the user provides vectors with
   EPSSetInitialSpace(), up to initialsize vectors will be used; and if the
   provided vectors are not enough, the solver completes the subspace with
   random vectors. In the case of EPSJDGetKrylovStart() being PETSC_TRUE, the solver
   gets the first vector provided by the user or, if not available, a random vector,
   and expands the Krylov basis up to initialsize vectors.

   Level: advanced

.seealso: EPSJDGetInitialSize(), EPSJDGetKrylovStart()
@*/
PetscErrorCode EPSJDSetInitialSize(EPS eps,PetscInt initialsize)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,initialsize,2);
  ierr = PetscTryMethod(eps,"EPSJDSetInitialSize_C",(EPS,PetscInt),(eps,initialsize));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDGetInitialSize - Returns the initial size of the searching space.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  initialsize - number of vectors of the initial searching subspace

   Notes:
   If EPSJDGetKrylovStart() is PETSC_FALSE and the user provides vectors with
   EPSSetInitialSpace(), up to initialsize vectors will be used; and if the
   provided vectors are not enough, the solver completes the subspace with
   random vectors. In the case of EPSJDGetKrylovStart() being PETSC_TRUE, the solver
   gets the first vector provided by the user or, if not available, a random vector,
   and expands the Krylov basis up to initialsize vectors.

   Level: advanced

.seealso: EPSJDSetInitialSize(), EPSJDGetKrylovStart()
@*/
PetscErrorCode EPSJDGetInitialSize(EPS eps,PetscInt *initialsize)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidIntPointer(initialsize,2);
  ierr = PetscUseMethod(eps,"EPSJDGetInitialSize_C",(EPS,PetscInt*),(eps,initialsize));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSJDSetFix_JD(EPS eps,PetscReal fix)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  if (fix == PETSC_DEFAULT || fix == PETSC_DECIDE) fix = 0.01;
  if (fix < 0.0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid fix value");
  data->fix = fix;
  PetscFunctionReturn(0);
}

/*@
   EPSJDSetFix - Sets the threshold for changing the target in the correction
   equation.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  fix - threshold for changing the target

   Options Database Key:
.  -eps_jd_fix - the fix value

   Note:
   The target in the correction equation is fixed at the first iterations.
   When the norm of the residual vector is lower than the fix value,
   the target is set to the corresponding eigenvalue.

   Level: advanced

.seealso: EPSJDGetFix()
@*/
PetscErrorCode EPSJDSetFix(EPS eps,PetscReal fix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,fix,2);
  ierr = PetscTryMethod(eps,"EPSJDSetFix_C",(EPS,PetscReal),(eps,fix));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSJDGetFix_JD(EPS eps,PetscReal *fix)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *fix = data->fix;
  PetscFunctionReturn(0);
}

/*@
   EPSJDGetFix - Returns the threshold for changing the target in the correction
   equation.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameter:
.  fix - threshold for changing the target

   Note:
   The target in the correction equation is fixed at the first iterations.
   When the norm of the residual vector is lower than the fix value,
   the target is set to the corresponding eigenvalue.

   Level: advanced

.seealso: EPSJDSetFix()
@*/
PetscErrorCode EPSJDGetFix(EPS eps,PetscReal *fix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidRealPointer(fix,2);
  ierr = PetscUseMethod(eps,"EPSJDGetFix_C",(EPS,PetscReal*),(eps,fix));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSJDSetConstCorrectionTol_JD(EPS eps,PetscBool constant)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  data->dynamic = PetscNot(constant);
  PetscFunctionReturn(0);
}

/*@
   EPSJDSetConstCorrectionTol - If true, deactivates the dynamic stopping criterion
   (also called Newton) that sets the KSP relative tolerance
   to 0.5**i, where i is the number of EPS iterations from the last converged value.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  constant - if false, the KSP relative tolerance is set to 0.5**i.

   Options Database Key:
.  -eps_jd_const_correction_tol - Deactivates the dynamic stopping criterion.

   Level: advanced

.seealso: EPSJDGetConstCorrectionTol()
@*/
PetscErrorCode EPSJDSetConstCorrectionTol(EPS eps,PetscBool constant)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,constant,2);
  ierr = PetscTryMethod(eps,"EPSJDSetConstCorrectionTol_C",(EPS,PetscBool),(eps,constant));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSJDGetConstCorrectionTol_JD(EPS eps,PetscBool *constant)
{
  EPS_DAVIDSON *data = (EPS_DAVIDSON*)eps->data;

  PetscFunctionBegin;
  *constant = PetscNot(data->dynamic);
  PetscFunctionReturn(0);
}

/*@
   EPSJDGetConstCorrectionTol - Returns a flag indicating if the dynamic stopping is being used for
   solving the correction equation. If the flag is false the KSP relative tolerance is set
   to 0.5**i, where i is the number of EPS iterations from the last converged value.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  constant - boolean flag indicating if the dynamic stopping criterion is not being used.

   Level: advanced

.seealso: EPSJDGetConstCorrectionTol()
@*/
PetscErrorCode EPSJDGetConstCorrectionTol(EPS eps,PetscBool *constant)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(constant,2);
  ierr = PetscUseMethod(eps,"EPSJDGetConstCorrectionTol_C",(EPS,PetscBool*),(eps,constant));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDSetBOrth - Selects the orthogonalization that will be used in the search
   subspace in case of generalized Hermitian problems.

   Logically Collective on eps

   Input Parameters:
+  eps   - the eigenproblem solver context
-  borth - whether to B-orthogonalize the search subspace

   Options Database Key:
.  -eps_jd_borth - Set the orthogonalization used in the search subspace

   Level: advanced

.seealso: EPSJDGetBOrth()
@*/
PetscErrorCode EPSJDSetBOrth(EPS eps,PetscBool borth)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,borth,2);
  ierr = PetscTryMethod(eps,"EPSJDSetBOrth_C",(EPS,PetscBool),(eps,borth));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   EPSJDGetBOrth - Returns the orthogonalization used in the search
   subspace in case of generalized Hermitian problems.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Output Parameters:
.  borth - whether to B-orthogonalize the search subspace

   Level: advanced

.seealso: EPSJDSetBOrth()
@*/
PetscErrorCode EPSJDGetBOrth(EPS eps,PetscBool *borth)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(borth,2);
  ierr = PetscUseMethod(eps,"EPSJDGetBOrth_C",(EPS,PetscBool*),(eps,borth));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_JD(EPS eps)
{
  PetscErrorCode ierr;
  EPS_DAVIDSON   *data;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&data);CHKERRQ(ierr);
  eps->data = (void*)data;

  data->blocksize   = 1;
  data->initialsize = 6;
  data->minv        = 6;
  data->plusk       = PETSC_DEFAULT;
  data->ipB         = PETSC_TRUE;
  data->fix         = 0.01;
  data->krylovstart = PETSC_FALSE;
  data->dynamic     = PETSC_FALSE;

  eps->useds = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_PRECOND;

  eps->ops->solve          = EPSSolve_XD;
  eps->ops->setup          = EPSSetUp_JD;
  eps->ops->setfromoptions = EPSSetFromOptions_JD;
  eps->ops->destroy        = EPSDestroy_JD;
  eps->ops->reset          = EPSReset_XD;
  eps->ops->view           = EPSView_JD;
  eps->ops->backtransform  = EPSBackTransform_Default;
  eps->ops->computevectors = EPSComputeVectors_XD;
  eps->ops->setdefaultst   = EPSSetDefaultST_JD;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetKrylovStart_C",EPSXDSetKrylovStart_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetKrylovStart_C",EPSXDGetKrylovStart_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetBlockSize_C",EPSXDSetBlockSize_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetBlockSize_C",EPSXDGetBlockSize_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetRestart_C",EPSXDSetRestart_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetRestart_C",EPSXDGetRestart_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetInitialSize_C",EPSXDSetInitialSize_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetInitialSize_C",EPSXDGetInitialSize_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetFix_C",EPSJDSetFix_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetFix_C",EPSJDGetFix_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetConstCorrectionTol_C",EPSJDSetConstCorrectionTol_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetConstCorrectionTol_C",EPSJDGetConstCorrectionTol_JD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDSetBOrth_C",EPSXDSetBOrth_XD);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSJDGetBOrth_C",EPSXDGetBOrth_XD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

