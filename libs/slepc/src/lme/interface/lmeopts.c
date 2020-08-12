/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   LME routines related to options that can be set via the command-line
   or procedurally
*/

#include <slepc/private/lmeimpl.h>   /*I "slepclme.h" I*/
#include <petscdraw.h>

/*@C
   LMEMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective on lme

   Input Parameters:
+  lme      - the eigensolver context
.  name     - the monitor option name
.  help     - message indicating what monitoring is done
.  manual   - manual page for the monitor
-  monitor  - the monitor function, whose context is a PetscViewerAndFormat

   Level: developer

.seealso: LMEMonitorSet()
@*/
PetscErrorCode LMEMonitorSetFromOptions(LME lme,const char name[],const char help[],const char manual[],PetscErrorCode (*monitor)(LME,PetscInt,PetscReal,PetscViewerAndFormat*))
{
  PetscErrorCode       ierr;
  PetscBool            flg;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerAndFormat *vf;

  PetscFunctionBegin;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)lme),((PetscObject)lme)->options,((PetscObject)lme)->prefix,name,&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerAndFormatCreate(viewer,format,&vf);CHKERRQ(ierr);
    ierr = PetscObjectDereference((PetscObject)viewer);CHKERRQ(ierr);
    ierr = LMEMonitorSet(lme,(PetscErrorCode (*)(LME,PetscInt,PetscReal,void*))monitor,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   LMESetFromOptions - Sets LME options from the options database.
   This routine must be called before LMESetUp() if the user is to be
   allowed to set the solver type.

   Collective on lme

   Input Parameters:
.  lme - the linear matrix equation solver context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner
@*/
PetscErrorCode LMESetFromOptions(LME lme)
{
  PetscErrorCode ierr;
  char           type[256];
  PetscBool      set,flg,flg1,flg2;
  PetscReal      r;
  PetscInt       i;
  PetscDrawLG    lg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  ierr = LMERegisterAll();CHKERRQ(ierr);
  ierr = PetscObjectOptionsBegin((PetscObject)lme);CHKERRQ(ierr);
    ierr = PetscOptionsFList("-lme_type","Linear matrix equation","LMESetType",LMEList,(char*)(((PetscObject)lme)->type_name?((PetscObject)lme)->type_name:LMEKRYLOV),type,256,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = LMESetType(lme,type);CHKERRQ(ierr);
    } else if (!((PetscObject)lme)->type_name) {
      ierr = LMESetType(lme,LMEKRYLOV);CHKERRQ(ierr);
    }

    ierr = PetscOptionsBoolGroupBegin("-lme_lyapunov","Continuous-time Lyapunov equation","LMESetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = LMESetProblemType(lme,LME_LYAPUNOV);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-lme_sylvester","Continuous-time Sylvester equation","LMESetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = LMESetProblemType(lme,LME_SYLVESTER);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-lme_gen_lyapunov","Generalized Lyapunov equation","LMESetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = LMESetProblemType(lme,LME_GEN_LYAPUNOV);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-lme_gen_sylvester","Generalized Sylvester equation","LMESetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = LMESetProblemType(lme,LME_GEN_SYLVESTER);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-lme_dt_lyapunov","Discrete-time Lyapunov equation","LMESetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = LMESetProblemType(lme,LME_DT_LYAPUNOV);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-lme_stein","Stein equation","LMESetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = LMESetProblemType(lme,LME_STEIN);CHKERRQ(ierr); }

    i = lme->max_it;
    ierr = PetscOptionsInt("-lme_max_it","Maximum number of iterations","LMESetTolerances",lme->max_it,&i,&flg1);CHKERRQ(ierr);
    if (!flg1) i = PETSC_DEFAULT;
    r = lme->tol;
    ierr = PetscOptionsReal("-lme_tol","Tolerance","LMESetTolerances",lme->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL:lme->tol,&r,&flg2);CHKERRQ(ierr);
    if (flg1 || flg2) { ierr = LMESetTolerances(lme,r,i);CHKERRQ(ierr); }

    ierr = PetscOptionsInt("-lme_ncv","Number of basis vectors","LMESetDimensions",lme->ncv,&i,&flg);CHKERRQ(ierr);
    if (flg) { ierr = LMESetDimensions(lme,i);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-lme_error_if_not_converged","Generate error if solver does not converge","LMESetErrorIfNotConverged",lme->errorifnotconverged,&lme->errorifnotconverged,NULL);CHKERRQ(ierr);

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to LMESetFromOptions()
    */
    ierr = PetscOptionsBool("-lme_monitor_cancel","Remove any hardwired monitor routines","LMEMonitorCancel",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = LMEMonitorCancel(lme);CHKERRQ(ierr);
    }
    /*
      Text monitors
    */
    ierr = LMEMonitorSetFromOptions(lme,"-lme_monitor","Monitor error estimate","LMEMonitorDefault",LMEMonitorDefault);CHKERRQ(ierr);
    /*
      Line graph monitors
    */
    ierr = PetscOptionsBool("-lme_monitor_lg","Monitor error estimate graphically","LMEMonitorSet",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = LMEMonitorLGCreate(PetscObjectComm((PetscObject)lme),NULL,"Error estimate",PETSC_DECIDE,PETSC_DECIDE,300,300,&lg);CHKERRQ(ierr);
      ierr = LMEMonitorSet(lme,LMEMonitorLG,lg,(PetscErrorCode (*)(void**))PetscDrawLGDestroy);CHKERRQ(ierr);
    }

    /* -----------------------------------------------------------------------*/
    ierr = PetscOptionsName("-lme_view","Print detailed information on solver used","LMEView",NULL);CHKERRQ(ierr);

    if (lme->ops->setfromoptions) {
      ierr = (*lme->ops->setfromoptions)(PetscOptionsObject,lme);CHKERRQ(ierr);
    }
    ierr = PetscObjectProcessOptionsHandlers(PetscOptionsObject,(PetscObject)lme);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  if (!lme->V) { ierr = LMEGetBV(lme,&lme->V);CHKERRQ(ierr); }
  ierr = BVSetFromOptions(lme->V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   LMESetProblemType - Specifies the type of matrix equation to be solved.

   Logically Collective on lme

   Input Parameters:
+  lme  - the linear matrix equation solver context
-  type - a known type of matrix equation

   Options Database Keys:
+  -lme_lyapunov - continuous-time Lyapunov equation A*X+X*A'=-C
.  -lme_sylvester - continuous-time Sylvester equation A*X+X*B=C
.  -lme_gen_lyapunov - generalized Lyapunov equation A*X*D'+D*X*A'=-C
.  -lme_gen_sylvester - generalized Sylvester equation A*X*E+D*X*B=C
.  -lme_dt_lyapunov - discrete-time Lyapunov equation A*X*A'-X=-C
-  -lme_stein - Stein equation A*X*E+X=C

   Notes:
   The coefficient matrices A, B, D, E must be provided via LMESetCoefficients(),
   but some of them are optional depending on the matrix equation.

.vb
                            equation              A    B    D    E
                          -----------------      ---  ---  ---  ---
       LME_LYAPUNOV        A*X+X*A'=-C           yes (A-t)  -    -
       LME_SYLVESTER       A*X+X*B=C             yes  yes   -    -
       LME_GEN_LYAPUNOV    A*X*D'+D*X*A'=-C      yes (A-t) yes (D-t)
       LME_GEN_SYLVESTER   A*X*E+D*X*B=C         yes  yes  yes  yes
       LME_DT_LYAPUNOV     A*X*A'-X=-C           yes   -    -  (A-t)
       LME_STEIN           A*X*E+X=C             yes   -    -   yes
.ve

   In the above table, the notation (A-t) means that this matrix need
   not be passed, but the user may choose to pass an explicit transpose
   of matrix A (for improved efficiency).

   Also note that some of the equation types impose restrictions on the
   properties of the coefficient matrices and possibly on the right-hand
   side C.

   Level: beginner

.seealso: LMESetCoefficients(), LMESetType(), LMEGetProblemType(), LMEProblemType
@*/
PetscErrorCode LMESetProblemType(LME lme,LMEProblemType type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveEnum(lme,type,2);
  if (type == lme->problem_type) PetscFunctionReturn(0);
  switch (type) {
    case LME_LYAPUNOV:
    case LME_SYLVESTER:
    case LME_GEN_LYAPUNOV:
    case LME_GEN_SYLVESTER:
    case LME_DT_LYAPUNOV:
    case LME_STEIN:
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_WRONG,"Unknown matrix equation type");
  }
  lme->problem_type = type;
  lme->setupcalled  = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@
   LMEGetProblemType - Gets the matrix equation type from the LME object.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  type - name of LME problem type

   Level: intermediate

.seealso: LMESetProblemType(), LMEProblemType
@*/
PetscErrorCode LMEGetProblemType(LME lme,LMEProblemType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidPointer(type,2);
  *type = lme->problem_type;
  PetscFunctionReturn(0);
}

/*@C
   LMEGetTolerances - Gets the tolerance and maximum iteration count used
   by the LME convergence tests.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: LMESetTolerances()
@*/
PetscErrorCode LMEGetTolerances(LME lme,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (tol)    *tol    = lme->tol;
  if (maxits) *maxits = lme->max_it;
  PetscFunctionReturn(0);
}

/*@
   LMESetTolerances - Sets the tolerance and maximum iteration count used
   by the LME convergence tests.

   Logically Collective on lme

   Input Parameters:
+  lme - the linear matrix equation solver context
.  tol - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -lme_tol <tol> - Sets the convergence tolerance
-  -lme_max_it <maxits> - Sets the maximum number of iterations allowed

   Notes:
   Use PETSC_DEFAULT for either argument to assign a reasonably good value.

   Level: intermediate

.seealso: LMEGetTolerances()
@*/
PetscErrorCode LMESetTolerances(LME lme,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveReal(lme,tol,2);
  PetscValidLogicalCollectiveInt(lme,maxits,3);
  if (tol == PETSC_DEFAULT) {
    lme->tol = PETSC_DEFAULT;
    lme->setupcalled = 0;
  } else {
    if (tol <= 0.0) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    lme->tol = tol;
  }
  if (maxits == PETSC_DEFAULT || maxits == PETSC_DECIDE) {
    lme->max_it = PETSC_DEFAULT;
    lme->setupcalled = 0;
  } else {
    if (maxits <= 0) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    lme->max_it = maxits;
  }
  PetscFunctionReturn(0);
}

/*@
   LMEGetDimensions - Gets the dimension of the subspace used by the solver.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  ncv - the maximum dimension of the subspace to be used by the solver

   Level: intermediate

.seealso: LMESetDimensions()
@*/
PetscErrorCode LMEGetDimensions(LME lme,PetscInt *ncv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidIntPointer(ncv,2);
  *ncv = lme->ncv;
  PetscFunctionReturn(0);
}

/*@
   LMESetDimensions - Sets the dimension of the subspace to be used by the solver.

   Logically Collective on lme

   Input Parameters:
+  lme - the linear matrix equation solver context
-  ncv - the maximum dimension of the subspace to be used by the solver

   Options Database Keys:
.  -lme_ncv <ncv> - Sets the dimension of the subspace

   Notes:
   Use PETSC_DEFAULT for ncv to assign a reasonably good value, which is
   dependent on the solution method.

   Level: intermediate

.seealso: LMEGetDimensions()
@*/
PetscErrorCode LMESetDimensions(LME lme,PetscInt ncv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveInt(lme,ncv,2);
  if (ncv == PETSC_DECIDE || ncv == PETSC_DEFAULT) {
    lme->ncv = PETSC_DEFAULT;
  } else {
    if (ncv<1) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    lme->ncv = ncv;
  }
  lme->setupcalled = 0;
  PetscFunctionReturn(0);
}

/*@
   LMESetErrorIfNotConverged - Causes LMESolve() to generate an error if the
   solver has not converged.

   Logically Collective on lme

   Input Parameters:
+  lme - the linear matrix equation solver context
-  flg - PETSC_TRUE indicates you want the error generated

   Options Database Keys:
.  -lme_error_if_not_converged - this takes an optional truth value (0/1/no/yes/true/false)

   Level: intermediate

   Note:
   Normally SLEPc continues if the solver fails to converge, you can call
   LMEGetConvergedReason() after a LMESolve() to determine if it has converged.

.seealso: LMEGetErrorIfNotConverged()
@*/
PetscErrorCode LMESetErrorIfNotConverged(LME lme,PetscBool flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidLogicalCollectiveBool(lme,flg,2);
  lme->errorifnotconverged = flg;
  PetscFunctionReturn(0);
}

/*@
   LMEGetErrorIfNotConverged - Return a flag indicating whether LMESolve() will
   generate an error if the solver does not converge.

   Not Collective

   Input Parameter:
.  lme - the linear matrix equation solver context

   Output Parameter:
.  flag - PETSC_TRUE if it will generate an error, else PETSC_FALSE

   Level: intermediate

.seealso:  LMESetErrorIfNotConverged()
@*/
PetscErrorCode LMEGetErrorIfNotConverged(LME lme,PetscBool *flag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidBoolPointer(flag,2);
  *flag = lme->errorifnotconverged;
  PetscFunctionReturn(0);
}

/*@C
   LMESetOptionsPrefix - Sets the prefix used for searching for all
   LME options in the database.

   Logically Collective on lme

   Input Parameters:
+  lme - the linear matrix equation solver context
-  prefix - the prefix string to prepend to all LME option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different LME contexts, one could call
.vb
      LMESetOptionsPrefix(lme1,"fun1_")
      LMESetOptionsPrefix(lme2,"fun2_")
.ve

   Level: advanced

.seealso: LMEAppendOptionsPrefix(), LMEGetOptionsPrefix()
@*/
PetscErrorCode LMESetOptionsPrefix(LME lme,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (!lme->V) { ierr = LMEGetBV(lme,&lme->V);CHKERRQ(ierr); }
  ierr = BVSetOptionsPrefix(lme->V,prefix);CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)lme,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   LMEAppendOptionsPrefix - Appends to the prefix used for searching for all
   LME options in the database.

   Logically Collective on lme

   Input Parameters:
+  lme - the linear matrix equation solver context
-  prefix - the prefix string to prepend to all LME option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: LMESetOptionsPrefix(), LMEGetOptionsPrefix()
@*/
PetscErrorCode LMEAppendOptionsPrefix(LME lme,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (!lme->V) { ierr = LMEGetBV(lme,&lme->V);CHKERRQ(ierr); }
  ierr = BVAppendOptionsPrefix(lme->V,prefix);CHKERRQ(ierr);
  ierr = PetscObjectAppendOptionsPrefix((PetscObject)lme,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   LMEGetOptionsPrefix - Gets the prefix used for searching for all
   LME options in the database.

   Not Collective

   Input Parameters:
.  lme - the linear matrix equation solver context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: LMESetOptionsPrefix(), LMEAppendOptionsPrefix()
@*/
PetscErrorCode LMEGetOptionsPrefix(LME lme,const char *prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidPointer(prefix,2);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)lme,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
