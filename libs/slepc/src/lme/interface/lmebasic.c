/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic LME routines
*/

#include <slepc/private/lmeimpl.h>      /*I "slepclme.h" I*/

PetscFunctionList LMEList = 0;
PetscBool         LMERegisterAllCalled = PETSC_FALSE;
PetscClassId      LME_CLASSID = 0;
PetscLogEvent     LME_SetUp = 0,LME_Solve = 0,LME_ComputeError = 0;

/*@C
   LMEView - Prints the LME data structure.

   Collective on lme

   Input Parameters:
+  lme - the linear matrix equation solver context
-  viewer - optional visualization context

   Options Database Key:
.  -lme_view -  Calls LMEView() at end of LMESolve()

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
PetscErrorCode LMEView(LME lme,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isascii;
  const char     *eqname[] = {
                   "continuous-time Lyapunov",
                   "continuous-time Sylvester",
                   "generalized Lyapunov",
                   "generalized Sylvester",
                   "Stein",
                   "discrete-time Lyapunov"
  };

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)lme),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(lme,1,viewer,2);

  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscObjectPrintClassNamePrefixType((PetscObject)lme,viewer);CHKERRQ(ierr);
    if (lme->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*lme->ops->view)(lme,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"  equation type: %s\n",eqname[lme->problem_type]);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  number of column vectors (ncv): %D\n",lme->ncv);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum number of iterations: %D\n",lme->max_it);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  tolerance: %g\n",(double)lme->tol);CHKERRQ(ierr);
  } else {
    if (lme->ops->view) {
      ierr = (*lme->ops->view)(lme,viewer);CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
  if (!lme->V) { ierr = LMEGetBV(lme,&lme->V);CHKERRQ(ierr); }
  ierr = BVView(lme->V,viewer);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   LMEViewFromOptions - View from options

   Collective on LME

   Input Parameters:
+  lme  - the linear matrix equation context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: LMEView(), LMECreate()
@*/
PetscErrorCode LMEViewFromOptions(LME lme,PetscObject obj,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  ierr = PetscObjectViewFromOptions((PetscObject)lme,obj,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/*@C
   LMEReasonView - Displays the reason an LME solve converged or diverged.

   Collective on lme

   Input Parameters:
+  lme - the linear matrix equation context
-  viewer - the viewer to display the reason

   Options Database Keys:
.  -lme_converged_reason - print reason for convergence, and number of iterations

   Level: intermediate

.seealso: LMESetTolerances(), LMEGetIterationNumber()
@*/
PetscErrorCode LMEReasonView(LME lme,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isAscii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isAscii);CHKERRQ(ierr);
  if (isAscii) {
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)lme)->tablevel);CHKERRQ(ierr);
    if (lme->reason > 0) {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Linear matrix equation solve converged due to %s; iterations %D\n",((PetscObject)lme)->prefix?((PetscObject)lme)->prefix:"",LMEConvergedReasons[lme->reason],lme->its);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Linear matrix equation solve did not converge due to %s; iterations %D\n",((PetscObject)lme)->prefix?((PetscObject)lme)->prefix:"",LMEConvergedReasons[lme->reason],lme->its);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)lme)->tablevel);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   LMEReasonViewFromOptions - Processes command line options to determine if/how
   the LME converged reason is to be viewed.

   Collective on lme

   Input Parameter:
.  lme - the linear matrix equation context

   Level: developer
@*/
PetscErrorCode LMEReasonViewFromOptions(LME lme)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)lme),((PetscObject)lme)->options,((PetscObject)lme)->prefix,"-lme_converged_reason",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = LMEReasonView(lme,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@
   LMECreate - Creates the default LME context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  lme - location to put the LME context

   Note:
   The default LME type is LMEKRYLOV

   Level: beginner

.seealso: LMESetUp(), LMESolve(), LMEDestroy(), LME
@*/
PetscErrorCode LMECreate(MPI_Comm comm,LME *outlme)
{
  PetscErrorCode ierr;
  LME            lme;

  PetscFunctionBegin;
  PetscValidPointer(outlme,2);
  *outlme = 0;
  ierr = LMEInitializePackage();CHKERRQ(ierr);
  ierr = SlepcHeaderCreate(lme,LME_CLASSID,"LME","Linear Matrix Equation","LME",comm,LMEDestroy,LMEView);CHKERRQ(ierr);

  lme->A               = NULL;
  lme->B               = NULL;
  lme->D               = NULL;
  lme->E               = NULL;
  lme->C               = NULL;
  lme->X               = NULL;
  lme->problem_type    = LME_LYAPUNOV;
  lme->max_it          = PETSC_DEFAULT;
  lme->ncv             = PETSC_DEFAULT;
  lme->tol             = PETSC_DEFAULT;
  lme->errorifnotconverged = PETSC_FALSE;

  lme->numbermonitors  = 0;

  lme->V               = NULL;
  lme->nwork           = 0;
  lme->work            = NULL;
  lme->data            = NULL;

  lme->its             = 0;
  lme->errest          = 0;
  lme->setupcalled     = 0;
  lme->reason          = LME_CONVERGED_ITERATING;

  *outlme = lme;
  PetscFunctionReturn(0);
}

/*@C
   LMESetType - Selects the particular solver to be used in the LME object.

   Logically Collective on lme

   Input Parameters:
+  lme  - the linear matrix equation context
-  type - a known method

   Options Database Key:
.  -lme_type <method> - Sets the method; use -help for a list
    of available methods

   Notes:
   See "slepc/include/slepclme.h" for available methods. The default
   is LMEKRYLOV

   Normally, it is best to use the LMESetFromOptions() command and
   then set the LME type from the options database rather than by using
   this routine.  Using the options database provides the user with
   maximum flexibility in evaluating the different available methods.
   The LMESetType() routine is provided for those situations where it
   is necessary to set the iterative solver independently of the command
   line or options database.

   Level: intermediate

.seealso: LMEType
@*/
PetscErrorCode LMESetType(LME lme,LMEType type)
{
  PetscErrorCode ierr,(*r)(LME);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidCharPointer(type,2);

  ierr = PetscObjectTypeCompare((PetscObject)lme,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);

  ierr = PetscFunctionListFind(LMEList,type,&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown LME type given: %s",type);

  if (lme->ops->destroy) { ierr = (*lme->ops->destroy)(lme);CHKERRQ(ierr); }
  ierr = PetscMemzero(lme->ops,sizeof(struct _LMEOps));CHKERRQ(ierr);

  lme->setupcalled = 0;
  ierr = PetscObjectChangeTypeName((PetscObject)lme,type);CHKERRQ(ierr);
  ierr = (*r)(lme);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   LMEGetType - Gets the LME type as a string from the LME object.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation context

   Output Parameter:
.  name - name of LME method

   Level: intermediate

.seealso: LMESetType()
@*/
PetscErrorCode LMEGetType(LME lme,LMEType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidPointer(type,2);
  *type = ((PetscObject)lme)->type_name;
  PetscFunctionReturn(0);
}

/*@C
   LMERegister - Adds a method to the linear matrix equation solver package.

   Not Collective

   Input Parameters:
+  name - name of a new user-defined solver
-  function - routine to create the solver context

   Notes:
   LMERegister() may be called multiple times to add several user-defined solvers.

   Sample usage:
.vb
    LMERegister("my_solver",MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     LMESetType(lme,"my_solver")
   or at runtime via the option
$     -lme_type my_solver

   Level: advanced

.seealso: LMERegisterAll()
@*/
PetscErrorCode LMERegister(const char *name,PetscErrorCode (*function)(LME))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = LMEInitializePackage();CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&LMEList,name,function);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   LMEReset - Resets the LME context to the initial state (prior to setup)
   and destroys any allocated Vecs and Mats.

   Collective on lme

   Input Parameter:
.  lme - linear matrix equation context obtained from LMECreate()

   Level: advanced

.seealso: LMEDestroy()
@*/
PetscErrorCode LMEReset(LME lme)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (lme) PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (!lme) PetscFunctionReturn(0);
  if (lme->ops->reset) { ierr = (lme->ops->reset)(lme);CHKERRQ(ierr); }
  ierr = MatDestroy(&lme->A);CHKERRQ(ierr);
  ierr = MatDestroy(&lme->B);CHKERRQ(ierr);
  ierr = MatDestroy(&lme->D);CHKERRQ(ierr);
  ierr = MatDestroy(&lme->E);CHKERRQ(ierr);
  ierr = MatDestroy(&lme->C);CHKERRQ(ierr);
  ierr = MatDestroy(&lme->X);CHKERRQ(ierr);
  ierr = BVDestroy(&lme->V);CHKERRQ(ierr);
  ierr = VecDestroyVecs(lme->nwork,&lme->work);CHKERRQ(ierr);
  lme->nwork = 0;
  lme->setupcalled = 0;
  PetscFunctionReturn(0);
}

/*@
   LMEDestroy - Destroys the LME context.

   Collective on lme

   Input Parameter:
.  lme - linear matrix equation context obtained from LMECreate()

   Level: beginner

.seealso: LMECreate(), LMESetUp(), LMESolve()
@*/
PetscErrorCode LMEDestroy(LME *lme)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*lme) PetscFunctionReturn(0);
  PetscValidHeaderSpecific(*lme,LME_CLASSID,1);
  if (--((PetscObject)(*lme))->refct > 0) { *lme = 0; PetscFunctionReturn(0); }
  ierr = LMEReset(*lme);CHKERRQ(ierr);
  if ((*lme)->ops->destroy) { ierr = (*(*lme)->ops->destroy)(*lme);CHKERRQ(ierr); }
  ierr = LMEMonitorCancel(*lme);CHKERRQ(ierr);
  ierr = PetscHeaderDestroy(lme);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   LMESetBV - Associates a basis vectors object to the linear matrix equation solver.

   Collective on lme

   Input Parameters:
+  lme - linear matrix equation context obtained from LMECreate()
-  bv  - the basis vectors object

   Note:
   Use LMEGetBV() to retrieve the basis vectors context (for example,
   to free it at the end of the computations).

   Level: advanced

.seealso: LMEGetBV()
@*/
PetscErrorCode LMESetBV(LME lme,BV bv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidHeaderSpecific(bv,BV_CLASSID,2);
  PetscCheckSameComm(lme,1,bv,2);
  ierr = PetscObjectReference((PetscObject)bv);CHKERRQ(ierr);
  ierr = BVDestroy(&lme->V);CHKERRQ(ierr);
  lme->V = bv;
  ierr = PetscLogObjectParent((PetscObject)lme,(PetscObject)lme->V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   LMEGetBV - Obtain the basis vectors object associated to the matrix
   function solver.

   Not Collective

   Input Parameters:
.  lme - linear matrix equation context obtained from LMECreate()

   Output Parameter:
.  bv - basis vectors context

   Level: advanced

.seealso: LMESetBV()
@*/
PetscErrorCode LMEGetBV(LME lme,BV *bv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidPointer(bv,2);
  if (!lme->V) {
    ierr = BVCreate(PetscObjectComm((PetscObject)lme),&lme->V);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)lme->V,(PetscObject)lme,0);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)lme,(PetscObject)lme->V);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)lme->V,((PetscObject)lme)->options);CHKERRQ(ierr);
  }
  *bv = lme->V;
  PetscFunctionReturn(0);
}

