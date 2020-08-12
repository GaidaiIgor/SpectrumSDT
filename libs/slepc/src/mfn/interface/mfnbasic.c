/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic MFN routines
*/

#include <slepc/private/mfnimpl.h>      /*I "slepcmfn.h" I*/

PetscFunctionList MFNList = 0;
PetscBool         MFNRegisterAllCalled = PETSC_FALSE;
PetscClassId      MFN_CLASSID = 0;
PetscLogEvent     MFN_SetUp = 0,MFN_Solve = 0;

/*@C
   MFNView - Prints the MFN data structure.

   Collective on mfn

   Input Parameters:
+  mfn - the matrix function solver context
-  viewer - optional visualization context

   Options Database Key:
.  -mfn_view -  Calls MFNView() at end of MFNSolve()

   Note:
   The available visualization contexts include
+     PETSC_VIEWER_STDOUT_SELF - standard output (default)
-     PETSC_VIEWER_STDOUT_WORLD - synchronized standard
         output where only the first processor opens
         the file.  All other processors send their
         data to the first processor to print.

   The user can open an alternative visualization context with
   PetscViewerASCIIOpen() - output to a specified file.

   Level: beginner
@*/
PetscErrorCode MFNView(MFN mfn,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isascii;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)mfn),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(mfn,1,viewer,2);

  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscObjectPrintClassNamePrefixType((PetscObject)mfn,viewer);CHKERRQ(ierr);
    if (mfn->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*mfn->ops->view)(mfn,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"  number of column vectors (ncv): %D\n",mfn->ncv);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum number of iterations: %D\n",mfn->max_it);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  tolerance: %g\n",(double)mfn->tol);CHKERRQ(ierr);
  } else {
    if (mfn->ops->view) {
      ierr = (*mfn->ops->view)(mfn,viewer);CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
  if (!mfn->V) { ierr = MFNGetFN(mfn,&mfn->fn);CHKERRQ(ierr); }
  ierr = FNView(mfn->fn,viewer);CHKERRQ(ierr);
  if (!mfn->V) { ierr = MFNGetBV(mfn,&mfn->V);CHKERRQ(ierr); }
  ierr = BVView(mfn->V,viewer);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   MFNViewFromOptions - View from options

   Collective on MFN

   Input Parameters:
+  mfn  - the matrix function context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: MFNView(), MFNCreate()
@*/
PetscErrorCode MFNViewFromOptions(MFN mfn,PetscObject obj,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  ierr = PetscObjectViewFromOptions((PetscObject)mfn,obj,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/*@C
   MFNReasonView - Displays the reason an MFN solve converged or diverged.

   Collective on mfn

   Input Parameters:
+  mfn - the matrix function context
-  viewer - the viewer to display the reason

   Options Database Keys:
.  -mfn_converged_reason - print reason for convergence, and number of iterations

   Level: intermediate

.seealso: MFNSetTolerances(), MFNGetIterationNumber()
@*/
PetscErrorCode MFNReasonView(MFN mfn,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isAscii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isAscii);CHKERRQ(ierr);
  if (isAscii) {
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)mfn)->tablevel);CHKERRQ(ierr);
    if (mfn->reason > 0) {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Matrix function solve converged due to %s; iterations %D\n",((PetscObject)mfn)->prefix?((PetscObject)mfn)->prefix:"",MFNConvergedReasons[mfn->reason],mfn->its);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Matrix function solve did not converge due to %s; iterations %D\n",((PetscObject)mfn)->prefix?((PetscObject)mfn)->prefix:"",MFNConvergedReasons[mfn->reason],mfn->its);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)mfn)->tablevel);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   MFNReasonViewFromOptions - Processes command line options to determine if/how
   the MFN converged reason is to be viewed.

   Collective on mfn

   Input Parameter:
.  mfn - the matrix function context

   Level: developer
@*/
PetscErrorCode MFNReasonViewFromOptions(MFN mfn)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)mfn),((PetscObject)mfn)->options,((PetscObject)mfn)->prefix,"-mfn_converged_reason",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = MFNReasonView(mfn,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@
   MFNCreate - Creates the default MFN context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  mfn - location to put the MFN context

   Note:
   The default MFN type is MFNKRYLOV

   Level: beginner

.seealso: MFNSetUp(), MFNSolve(), MFNDestroy(), MFN
@*/
PetscErrorCode MFNCreate(MPI_Comm comm,MFN *outmfn)
{
  PetscErrorCode ierr;
  MFN            mfn;

  PetscFunctionBegin;
  PetscValidPointer(outmfn,2);
  *outmfn = 0;
  ierr = MFNInitializePackage();CHKERRQ(ierr);
  ierr = SlepcHeaderCreate(mfn,MFN_CLASSID,"MFN","Matrix Function","MFN",comm,MFNDestroy,MFNView);CHKERRQ(ierr);

  mfn->A               = NULL;
  mfn->fn              = NULL;
  mfn->max_it          = PETSC_DEFAULT;
  mfn->ncv             = PETSC_DEFAULT;
  mfn->tol             = PETSC_DEFAULT;
  mfn->errorifnotconverged = PETSC_FALSE;

  mfn->numbermonitors  = 0;

  mfn->V               = NULL;
  mfn->nwork           = 0;
  mfn->work            = NULL;
  mfn->data            = NULL;

  mfn->its             = 0;
  mfn->nv              = 0;
  mfn->errest          = 0;
  mfn->setupcalled     = 0;
  mfn->reason          = MFN_CONVERGED_ITERATING;

  *outmfn = mfn;
  PetscFunctionReturn(0);
}

/*@C
   MFNSetType - Selects the particular solver to be used in the MFN object.

   Logically Collective on mfn

   Input Parameters:
+  mfn  - the matrix function context
-  type - a known method

   Options Database Key:
.  -mfn_type <method> - Sets the method; use -help for a list
    of available methods

   Notes:
   See "slepc/include/slepcmfn.h" for available methods. The default
   is MFNKRYLOV

   Normally, it is best to use the MFNSetFromOptions() command and
   then set the MFN type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different available methods.
   The MFNSetType() routine is provided for those situations where it
   is necessary to set the iterative solver independently of the command
   line or options database.

   Level: intermediate

.seealso: MFNType
@*/
PetscErrorCode MFNSetType(MFN mfn,MFNType type)
{
  PetscErrorCode ierr,(*r)(MFN);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidCharPointer(type,2);

  ierr = PetscObjectTypeCompare((PetscObject)mfn,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);

  ierr = PetscFunctionListFind(MFNList,type,&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown MFN type given: %s",type);

  if (mfn->ops->destroy) { ierr = (*mfn->ops->destroy)(mfn);CHKERRQ(ierr); }
  ierr = PetscMemzero(mfn->ops,sizeof(struct _MFNOps));CHKERRQ(ierr);

  mfn->setupcalled = 0;
  ierr = PetscObjectChangeTypeName((PetscObject)mfn,type);CHKERRQ(ierr);
  ierr = (*r)(mfn);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   MFNGetType - Gets the MFN type as a string from the MFN object.

   Not Collective

   Input Parameter:
.  mfn - the matrix function context

   Output Parameter:
.  name - name of MFN method

   Level: intermediate

.seealso: MFNSetType()
@*/
PetscErrorCode MFNGetType(MFN mfn,MFNType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidPointer(type,2);
  *type = ((PetscObject)mfn)->type_name;
  PetscFunctionReturn(0);
}

/*@C
   MFNRegister - Adds a method to the matrix function solver package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined solver
-  function - routine to create the solver context

   Notes:
   MFNRegister() may be called multiple times to add several user-defined solvers.

   Sample usage:
.vb
    MFNRegister("my_solver",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     MFNSetType(mfn,"my_solver")
   or at runtime via the option
$     -mfn_type my_solver

   Level: advanced

.seealso: MFNRegisterAll()
@*/
PetscErrorCode MFNRegister(const char *name,PetscErrorCode (*function)(MFN))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MFNInitializePackage();CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&MFNList,name,function);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   MFNReset - Resets the MFN context to the initial state (prior to setup)
   and destroys any allocated Vecs and Mats.

   Collective on mfn

   Input Parameter:
.  mfn - matrix function context obtained from MFNCreate()

   Level: advanced

.seealso: MFNDestroy()
@*/
PetscErrorCode MFNReset(MFN mfn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (mfn) PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  if (!mfn) PetscFunctionReturn(0);
  if (mfn->ops->reset) { ierr = (mfn->ops->reset)(mfn);CHKERRQ(ierr); }
  ierr = MatDestroy(&mfn->A);CHKERRQ(ierr);
  ierr = BVDestroy(&mfn->V);CHKERRQ(ierr);
  ierr = VecDestroyVecs(mfn->nwork,&mfn->work);CHKERRQ(ierr);
  mfn->nwork = 0;
  mfn->setupcalled = 0;
  PetscFunctionReturn(0);
}

/*@
   MFNDestroy - Destroys the MFN context.

   Collective on mfn

   Input Parameter:
.  mfn - matrix function context obtained from MFNCreate()

   Level: beginner

.seealso: MFNCreate(), MFNSetUp(), MFNSolve()
@*/
PetscErrorCode MFNDestroy(MFN *mfn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*mfn) PetscFunctionReturn(0);
  PetscValidHeaderSpecific(*mfn,MFN_CLASSID,1);
  if (--((PetscObject)(*mfn))->refct > 0) { *mfn = 0; PetscFunctionReturn(0); }
  ierr = MFNReset(*mfn);CHKERRQ(ierr);
  if ((*mfn)->ops->destroy) { ierr = (*(*mfn)->ops->destroy)(*mfn);CHKERRQ(ierr); }
  ierr = FNDestroy(&(*mfn)->fn);CHKERRQ(ierr);
  ierr = MatDestroy(&(*mfn)->AT);CHKERRQ(ierr);
  ierr = MFNMonitorCancel(*mfn);CHKERRQ(ierr);
  ierr = PetscHeaderDestroy(mfn);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   MFNSetBV - Associates a basis vectors object to the matrix function solver.

   Collective on mfn

   Input Parameters:
+  mfn - matrix function context obtained from MFNCreate()
-  bv  - the basis vectors object

   Note:
   Use MFNGetBV() to retrieve the basis vectors context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: MFNGetBV()
@*/
PetscErrorCode MFNSetBV(MFN mfn,BV bv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidHeaderSpecific(bv,BV_CLASSID,2);
  PetscCheckSameComm(mfn,1,bv,2);
  ierr = PetscObjectReference((PetscObject)bv);CHKERRQ(ierr);
  ierr = BVDestroy(&mfn->V);CHKERRQ(ierr);
  mfn->V = bv;
  ierr = PetscLogObjectParent((PetscObject)mfn,(PetscObject)mfn->V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   MFNGetBV - Obtain the basis vectors object associated to the matrix
   function solver.

   Not Collective

   Input Parameters:
.  mfn - matrix function context obtained from MFNCreate()

   Output Parameter:
.  bv - basis vectors context

   Level: advanced

.seealso: MFNSetBV()
@*/
PetscErrorCode MFNGetBV(MFN mfn,BV *bv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidPointer(bv,2);
  if (!mfn->V) {
    ierr = BVCreate(PetscObjectComm((PetscObject)mfn),&mfn->V);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)mfn->V,(PetscObject)mfn,0);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)mfn,(PetscObject)mfn->V);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)mfn->V,((PetscObject)mfn)->options);CHKERRQ(ierr);
  }
  *bv = mfn->V;
  PetscFunctionReturn(0);
}

/*@
   MFNSetFN - Specifies the function to be computed.

   Collective on mfn

   Input Parameters:
+  mfn - matrix function context obtained from MFNCreate()
-  fn  - the math function object

   Note:
   Use MFNGetFN() to retrieve the math function context (for example,
   to free it at the end of the computations).

   Level: beginner

.seealso: MFNGetFN()
@*/
PetscErrorCode MFNSetFN(MFN mfn,FN fn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidHeaderSpecific(fn,FN_CLASSID,2);
  PetscCheckSameComm(mfn,1,fn,2);
  ierr = PetscObjectReference((PetscObject)fn);CHKERRQ(ierr);
  ierr = FNDestroy(&mfn->fn);CHKERRQ(ierr);
  mfn->fn = fn;
  ierr = PetscLogObjectParent((PetscObject)mfn,(PetscObject)mfn->fn);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   MFNGetFN - Obtain the math function object associated to the MFN object.

   Not Collective

   Input Parameters:
.  mfn - matrix function context obtained from MFNCreate()

   Output Parameter:
.  fn - math function context

   Level: beginner

.seealso: MFNSetFN()
@*/
PetscErrorCode MFNGetFN(MFN mfn,FN *fn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidPointer(fn,2);
  if (!mfn->fn) {
    ierr = FNCreate(PetscObjectComm((PetscObject)mfn),&mfn->fn);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)mfn->fn,(PetscObject)mfn,0);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)mfn,(PetscObject)mfn->fn);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)mfn->fn,((PetscObject)mfn)->options);CHKERRQ(ierr);
  }
  *fn = mfn->fn;
  PetscFunctionReturn(0);
}

