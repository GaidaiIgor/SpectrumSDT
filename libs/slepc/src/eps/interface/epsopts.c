/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   EPS routines related to options that can be set via the command-line
   or procedurally.
*/

#include <slepc/private/epsimpl.h>   /*I "slepceps.h" I*/
#include <petscdraw.h>

/*@C
   EPSMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective on eps

   Input Parameters:
+  eps      - the eigensolver context
.  name     - the monitor option name
.  help     - message indicating what monitoring is done
.  manual   - manual page for the monitor
.  monitor  - the monitor function, whose context is a PetscViewerAndFormat
-  trackall - whether this monitor tracks all eigenvalues or not

   Level: developer

.seealso: EPSMonitorSet(), EPSSetTrackAll(), EPSConvMonitorSetFromOptions()
@*/
PetscErrorCode EPSMonitorSetFromOptions(EPS eps,const char name[],const char help[],const char manual[],PetscErrorCode (*monitor)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscBool trackall)
{
  PetscErrorCode       ierr;
  PetscBool            flg;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerAndFormat *vf;

  PetscFunctionBegin;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)eps),((PetscObject)eps)->options,((PetscObject)eps)->prefix,name,&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerAndFormatCreate(viewer,format,&vf);CHKERRQ(ierr);
    ierr = PetscObjectDereference((PetscObject)viewer);CHKERRQ(ierr);
    ierr = EPSMonitorSet(eps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))monitor,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);
    if (trackall) {
      ierr = EPSSetTrackAll(eps,PETSC_TRUE);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@C
   EPSConvMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user (for monitors that only show iteration numbers of convergence).

   Collective on eps

   Input Parameters:
+  eps      - the eigensolver context
.  name     - the monitor option name
.  help     - message indicating what monitoring is done
.  manual   - manual page for the monitor
-  monitor  - the monitor function, whose context is a SlepcConvMonitor

   Level: developer

.seealso: EPSMonitorSet(), EPSMonitorSetFromOptions()
@*/
PetscErrorCode EPSConvMonitorSetFromOptions(EPS eps,const char name[],const char help[],const char manual[],PetscErrorCode (*monitor)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,SlepcConvMonitor))
{
  PetscErrorCode    ierr;
  PetscBool         flg;
  PetscViewer       viewer;
  PetscViewerFormat format;
  SlepcConvMonitor  ctx;

  PetscFunctionBegin;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)eps),((PetscObject)eps)->options,((PetscObject)eps)->prefix,name,&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = SlepcConvMonitorCreate(viewer,format,&ctx);CHKERRQ(ierr);
    ierr = PetscObjectDereference((PetscObject)viewer);CHKERRQ(ierr);
    ierr = EPSMonitorSet(eps,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))monitor,ctx,(PetscErrorCode (*)(void**))SlepcConvMonitorDestroy);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   EPSSetFromOptions - Sets EPS options from the options database.
   This routine must be called before EPSSetUp() if the user is to be
   allowed to set the solver type.

   Collective on eps

   Input Parameters:
.  eps - the eigensolver context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner
@*/
PetscErrorCode EPSSetFromOptions(EPS eps)
{
  PetscErrorCode ierr;
  char           type[256];
  PetscBool      set,flg,flg1,flg2,flg3,bval;
  PetscReal      r,array[2]={0,0};
  PetscScalar    s;
  PetscInt       i,j,k;
  PetscDrawLG    lg;
  EPSBalance     bal;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = EPSRegisterAll();CHKERRQ(ierr);
  ierr = PetscObjectOptionsBegin((PetscObject)eps);CHKERRQ(ierr);
    ierr = PetscOptionsFList("-eps_type","Eigensolver method","EPSSetType",EPSList,(char*)(((PetscObject)eps)->type_name?((PetscObject)eps)->type_name:EPSKRYLOVSCHUR),type,256,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = EPSSetType(eps,type);CHKERRQ(ierr);
    } else if (!((PetscObject)eps)->type_name) {
      ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRQ(ierr);
    }

    ierr = PetscOptionsBoolGroupBegin("-eps_hermitian","Hermitian eigenvalue problem","EPSSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_gen_hermitian","Generalized Hermitian eigenvalue problem","EPSSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_non_hermitian","Non-Hermitian eigenvalue problem","EPSSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetProblemType(eps,EPS_NHEP);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_gen_non_hermitian","Generalized non-Hermitian eigenvalue problem","EPSSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetProblemType(eps,EPS_GNHEP);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_pos_gen_non_hermitian","Generalized non-Hermitian eigenvalue problem with positive semi-definite B","EPSSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetProblemType(eps,EPS_PGNHEP);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-eps_gen_indefinite","Generalized Hermitian-indefinite eigenvalue problem","EPSSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetProblemType(eps,EPS_GHIEP);CHKERRQ(ierr); }

    ierr = PetscOptionsBoolGroupBegin("-eps_ritz","Rayleigh-Ritz extraction","EPSSetExtraction",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetExtraction(eps,EPS_RITZ);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_harmonic","Harmonic Ritz extraction","EPSSetExtraction",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetExtraction(eps,EPS_HARMONIC);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_harmonic_relative","Relative harmonic Ritz extraction","EPSSetExtraction",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetExtraction(eps,EPS_HARMONIC_RELATIVE);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_harmonic_right","Right harmonic Ritz extraction","EPSSetExtraction",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetExtraction(eps,EPS_HARMONIC_RIGHT);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_harmonic_largest","Largest harmonic Ritz extraction","EPSSetExtraction",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetExtraction(eps,EPS_HARMONIC_LARGEST);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_refined","Refined Ritz extraction","EPSSetExtraction",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetExtraction(eps,EPS_REFINED);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-eps_refined_harmonic","Refined harmonic Ritz extraction","EPSSetExtraction",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetExtraction(eps,EPS_REFINED_HARMONIC);CHKERRQ(ierr); }

    bal = eps->balance;
    ierr = PetscOptionsEnum("-eps_balance","Balancing method","EPSSetBalance",EPSBalanceTypes,(PetscEnum)bal,(PetscEnum*)&bal,&flg1);CHKERRQ(ierr);
    j = eps->balance_its;
    ierr = PetscOptionsInt("-eps_balance_its","Number of iterations in balancing","EPSSetBalance",eps->balance_its,&j,&flg2);CHKERRQ(ierr);
    r = eps->balance_cutoff;
    ierr = PetscOptionsReal("-eps_balance_cutoff","Cutoff value in balancing","EPSSetBalance",eps->balance_cutoff,&r,&flg3);CHKERRQ(ierr);
    if (flg1 || flg2 || flg3) { ierr = EPSSetBalance(eps,bal,j,r);CHKERRQ(ierr); }

    i = eps->max_it;
    ierr = PetscOptionsInt("-eps_max_it","Maximum number of iterations","EPSSetTolerances",eps->max_it,&i,&flg1);CHKERRQ(ierr);
    r = eps->tol;
    ierr = PetscOptionsReal("-eps_tol","Tolerance","EPSSetTolerances",eps->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL:eps->tol,&r,&flg2);CHKERRQ(ierr);
    if (flg1 || flg2) { ierr = EPSSetTolerances(eps,r,i);CHKERRQ(ierr); }

    ierr = PetscOptionsBoolGroupBegin("-eps_conv_rel","Relative error convergence test","EPSSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetConvergenceTest(eps,EPS_CONV_REL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_conv_norm","Convergence test relative to the eigenvalue and the matrix norms","EPSSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetConvergenceTest(eps,EPS_CONV_NORM);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_conv_abs","Absolute error convergence test","EPSSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetConvergenceTest(eps,EPS_CONV_ABS);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-eps_conv_user","User-defined convergence test","EPSSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetConvergenceTest(eps,EPS_CONV_USER);CHKERRQ(ierr); }

    ierr = PetscOptionsBoolGroupBegin("-eps_stop_basic","Stop iteration if all eigenvalues converged or max_it reached","EPSSetStoppingTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetStoppingTest(eps,EPS_STOP_BASIC);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-eps_stop_user","User-defined stopping test","EPSSetStoppingTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetStoppingTest(eps,EPS_STOP_USER);CHKERRQ(ierr); }

    i = eps->nev;
    ierr = PetscOptionsInt("-eps_nev","Number of eigenvalues to compute","EPSSetDimensions",eps->nev,&i,&flg1);CHKERRQ(ierr);
    j = eps->ncv;
    ierr = PetscOptionsInt("-eps_ncv","Number of basis vectors","EPSSetDimensions",eps->ncv,&j,&flg2);CHKERRQ(ierr);
    k = eps->mpd;
    ierr = PetscOptionsInt("-eps_mpd","Maximum dimension of projected problem","EPSSetDimensions",eps->mpd,&k,&flg3);CHKERRQ(ierr);
    if (flg1 || flg2 || flg3) { ierr = EPSSetDimensions(eps,i,j,k);CHKERRQ(ierr); }

    ierr = PetscOptionsBoolGroupBegin("-eps_largest_magnitude","Compute largest eigenvalues in magnitude","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_smallest_magnitude","Compute smallest eigenvalues in magnitude","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_largest_real","Compute eigenvalues with largest real parts","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_smallest_real","Compute eigenvalues with smallest real parts","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_largest_imaginary","Compute eigenvalues with largest imaginary parts","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_IMAGINARY);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_smallest_imaginary","Compute eigenvalues with smallest imaginary parts","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_IMAGINARY);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_target_magnitude","Compute eigenvalues closest to target","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_target_real","Compute eigenvalues with real parts closest to target","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_REAL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-eps_target_imaginary","Compute eigenvalues with imaginary parts closest to target","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_IMAGINARY);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-eps_all","Compute all eigenvalues in an interval or a region","EPSSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRQ(ierr); }

    ierr = PetscOptionsScalar("-eps_target","Value of the target","EPSSetTarget",eps->target,&s,&flg);CHKERRQ(ierr);
    if (flg) {
      if (eps->which!=EPS_TARGET_REAL && eps->which!=EPS_TARGET_IMAGINARY) {
        ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);CHKERRQ(ierr);
      }
      ierr = EPSSetTarget(eps,s);CHKERRQ(ierr);
    }

    k = 2;
    ierr = PetscOptionsRealArray("-eps_interval","Computational interval (two real values separated with a comma without spaces)","EPSSetInterval",array,&k,&flg);CHKERRQ(ierr);
    if (flg) {
      if (k<2) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_SIZ,"Must pass two values in -eps_interval (comma-separated without spaces)");
      ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRQ(ierr);
      ierr = EPSSetInterval(eps,array[0],array[1]);CHKERRQ(ierr);
    }

    ierr = PetscOptionsBool("-eps_true_residual","Compute true residuals explicitly","EPSSetTrueResidual",eps->trueres,&eps->trueres,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-eps_purify","Postprocess eigenvectors for purification","EPSSetPurify",eps->purify,&bval,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetPurify(eps,bval);CHKERRQ(ierr); }
    ierr = PetscOptionsBool("-eps_two_sided","Use two-sided variant (to compute left eigenvectors)","EPSSetTwoSided",eps->twosided,&bval,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSSetTwoSided(eps,bval);CHKERRQ(ierr); }

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to EPSSetFromOptions()
    */
    ierr = PetscOptionsBool("-eps_monitor_cancel","Remove any hardwired monitor routines","EPSMonitorCancel",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = EPSMonitorCancel(eps);CHKERRQ(ierr);
    }
    /*
      Text monitors
    */
    ierr = EPSMonitorSetFromOptions(eps,"-eps_monitor","Monitor first unconverged approximate eigenvalue and error estimate","EPSMonitorFirst",EPSMonitorFirst,PETSC_FALSE);CHKERRQ(ierr);
    ierr = EPSConvMonitorSetFromOptions(eps,"-eps_monitor_conv","Monitor approximate eigenvalues and error estimates as they converge","EPSMonitorConverged",EPSMonitorConverged);CHKERRQ(ierr);
    ierr = EPSMonitorSetFromOptions(eps,"-eps_monitor_all","Monitor approximate eigenvalues and error estimates","EPSMonitorAll",EPSMonitorAll,PETSC_TRUE);CHKERRQ(ierr);
    /*
      Line graph monitors
    */
    ierr = PetscOptionsBool("-eps_monitor_lg","Monitor first unconverged approximate eigenvalue and error estimate graphically","EPSMonitorSet",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = EPSMonitorLGCreate(PetscObjectComm((PetscObject)eps),NULL,"Error estimates",PETSC_DECIDE,PETSC_DECIDE,300,300,&lg);CHKERRQ(ierr);
      ierr = EPSMonitorSet(eps,EPSMonitorLG,lg,(PetscErrorCode (*)(void**))PetscDrawLGDestroy);CHKERRQ(ierr);
    }
    ierr = PetscOptionsBool("-eps_monitor_lg_all","Monitor error estimates graphically","EPSMonitorSet",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = EPSMonitorLGCreate(PetscObjectComm((PetscObject)eps),NULL,"Error estimates",PETSC_DECIDE,PETSC_DECIDE,300,300,&lg);CHKERRQ(ierr);
      ierr = EPSMonitorSet(eps,EPSMonitorLGAll,lg,(PetscErrorCode (*)(void**))PetscDrawLGDestroy);CHKERRQ(ierr);
      ierr = EPSSetTrackAll(eps,PETSC_TRUE);CHKERRQ(ierr);
    }

    /* -----------------------------------------------------------------------*/
    ierr = PetscOptionsName("-eps_view","Print detailed information on solver used","EPSView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-eps_view_vectors","View computed eigenvectors","EPSVectorsView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-eps_view_values","View computed eigenvalues","EPSValuesView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-eps_converged_reason","Print reason for convergence, and number of iterations","EPSReasonView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-eps_error_absolute","Print absolute errors of each eigenpair","EPSErrorView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-eps_error_relative","Print relative errors of each eigenpair","EPSErrorView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-eps_error_backward","Print backward errors of each eigenpair","EPSErrorView",NULL);CHKERRQ(ierr);

    if (eps->ops->setfromoptions) {
      ierr = (*eps->ops->setfromoptions)(PetscOptionsObject,eps);CHKERRQ(ierr);
    }
    ierr = PetscObjectProcessOptionsHandlers(PetscOptionsObject,(PetscObject)eps);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  if (!eps->V) { ierr = EPSGetBV(eps,&eps->V);CHKERRQ(ierr); }
  ierr = BVSetFromOptions(eps->V);CHKERRQ(ierr);
  if (!eps->rg) { ierr = EPSGetRG(eps,&eps->rg);CHKERRQ(ierr); }
  ierr = RGSetFromOptions(eps->rg);CHKERRQ(ierr);
  if (eps->useds) {
    if (!eps->ds) { ierr = EPSGetDS(eps,&eps->ds);CHKERRQ(ierr); }
    ierr = DSSetFromOptions(eps->ds);CHKERRQ(ierr);
  }
  if (!eps->st) { ierr = EPSGetST(eps,&eps->st);CHKERRQ(ierr); }
  ierr = EPSSetDefaultST(eps);CHKERRQ(ierr);
  ierr = STSetFromOptions(eps->st);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   EPSGetTolerances - Gets the tolerance and maximum iteration count used
   by the EPS convergence tests.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: EPSSetTolerances()
@*/
PetscErrorCode EPSGetTolerances(EPS eps,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (tol)    *tol    = eps->tol;
  if (maxits) *maxits = eps->max_it;
  PetscFunctionReturn(0);
}

/*@
   EPSSetTolerances - Sets the tolerance and maximum iteration count used
   by the EPS convergence tests.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigensolver context
.  tol - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -eps_tol <tol> - Sets the convergence tolerance
-  -eps_max_it <maxits> - Sets the maximum number of iterations allowed

   Notes:
   Use PETSC_DEFAULT for either argument to assign a reasonably good value.

   Level: intermediate

.seealso: EPSGetTolerances()
@*/
PetscErrorCode EPSSetTolerances(EPS eps,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveReal(eps,tol,2);
  PetscValidLogicalCollectiveInt(eps,maxits,3);
  if (tol == PETSC_DEFAULT) {
    eps->tol   = PETSC_DEFAULT;
    eps->state = EPS_STATE_INITIAL;
  } else {
    if (tol <= 0.0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    eps->tol = tol;
  }
  if (maxits == PETSC_DEFAULT || maxits == PETSC_DECIDE) {
    eps->max_it = PETSC_DEFAULT;
    eps->state  = EPS_STATE_INITIAL;
  } else {
    if (maxits <= 0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    eps->max_it = maxits;
  }
  PetscFunctionReturn(0);
}

/*@C
   EPSGetDimensions - Gets the number of eigenvalues to compute
   and the dimension of the subspace.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameters:
+  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Level: intermediate

.seealso: EPSSetDimensions()
@*/
PetscErrorCode EPSGetDimensions(EPS eps,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (nev) *nev = eps->nev;
  if (ncv) *ncv = eps->ncv;
  if (mpd) *mpd = eps->mpd;
  PetscFunctionReturn(0);
}

/*@
   EPSSetDimensions - Sets the number of eigenvalues to compute
   and the dimension of the subspace.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigensolver context
.  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Keys:
+  -eps_nev <nev> - Sets the number of eigenvalues
.  -eps_ncv <ncv> - Sets the dimension of the subspace
-  -eps_mpd <mpd> - Sets the maximum projected dimension

   Notes:
   Use PETSC_DEFAULT for ncv and mpd to assign a reasonably good value, which is
   dependent on the solution method.

   The parameters ncv and mpd are intimately related, so that the user is advised
   to set one of them at most. Normal usage is that
   (a) in cases where nev is small, the user sets ncv (a reasonable default is 2*nev); and
   (b) in cases where nev is large, the user sets mpd.

   The value of ncv should always be between nev and (nev+mpd), typically
   ncv=nev+mpd. If nev is not too large, mpd=nev is a reasonable choice, otherwise
   a smaller value should be used.

   When computing all eigenvalues in an interval, see EPSSetInterval(), these
   parameters lose relevance, and tuning must be done with
   EPSKrylovSchurSetDimensions().

   Level: intermediate

.seealso: EPSGetDimensions(), EPSSetInterval(), EPSKrylovSchurSetDimensions()
@*/
PetscErrorCode EPSSetDimensions(EPS eps,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveInt(eps,nev,2);
  PetscValidLogicalCollectiveInt(eps,ncv,3);
  PetscValidLogicalCollectiveInt(eps,mpd,4);
  if (nev<1) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nev. Must be > 0");
  eps->nev = nev;
  if (ncv == PETSC_DECIDE || ncv == PETSC_DEFAULT) {
    eps->ncv = PETSC_DEFAULT;
  } else {
    if (ncv<1) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    eps->ncv = ncv;
  }
  if (mpd == PETSC_DECIDE || mpd == PETSC_DEFAULT) {
    eps->mpd = PETSC_DEFAULT;
  } else {
    if (mpd<1) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    eps->mpd = mpd;
  }
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSSetWhichEigenpairs - Specifies which portion of the spectrum is
   to be sought.

   Logically Collective on eps

   Input Parameters:
+  eps   - eigensolver context obtained from EPSCreate()
-  which - the portion of the spectrum to be sought

   Possible values:
   The parameter 'which' can have one of these values

+     EPS_LARGEST_MAGNITUDE - largest eigenvalues in magnitude (default)
.     EPS_SMALLEST_MAGNITUDE - smallest eigenvalues in magnitude
.     EPS_LARGEST_REAL - largest real parts
.     EPS_SMALLEST_REAL - smallest real parts
.     EPS_LARGEST_IMAGINARY - largest imaginary parts
.     EPS_SMALLEST_IMAGINARY - smallest imaginary parts
.     EPS_TARGET_MAGNITUDE - eigenvalues closest to the target (in magnitude)
.     EPS_TARGET_REAL - eigenvalues with real part closest to target
.     EPS_TARGET_IMAGINARY - eigenvalues with imaginary part closest to target
.     EPS_ALL - all eigenvalues contained in a given interval or region
-     EPS_WHICH_USER - user defined ordering set with EPSSetEigenvalueComparison()

   Options Database Keys:
+   -eps_largest_magnitude - Sets largest eigenvalues in magnitude
.   -eps_smallest_magnitude - Sets smallest eigenvalues in magnitude
.   -eps_largest_real - Sets largest real parts
.   -eps_smallest_real - Sets smallest real parts
.   -eps_largest_imaginary - Sets largest imaginary parts
.   -eps_smallest_imaginary - Sets smallest imaginary parts
.   -eps_target_magnitude - Sets eigenvalues closest to target
.   -eps_target_real - Sets real parts closest to target
.   -eps_target_imaginary - Sets imaginary parts closest to target
-   -eps_all - Sets all eigenvalues in an interval or region

   Notes:
   Not all eigensolvers implemented in EPS account for all the possible values
   stated above. Also, some values make sense only for certain types of
   problems. If SLEPc is compiled for real numbers EPS_LARGEST_IMAGINARY
   and EPS_SMALLEST_IMAGINARY use the absolute value of the imaginary part
   for eigenvalue selection.

   The target is a scalar value provided with EPSSetTarget().

   The criterion EPS_TARGET_IMAGINARY is available only in case PETSc and
   SLEPc have been built with complex scalars.

   EPS_ALL is intended for use in combination with an interval (see
   EPSSetInterval()), when all eigenvalues within the interval are requested,
   or in the context of the CISS solver for computing all eigenvalues in a region.
   In those cases, the number of eigenvalues is unknown, so the nev parameter
   has a different sense, see EPSSetDimensions().

   Level: intermediate

.seealso: EPSGetWhichEigenpairs(), EPSSetTarget(), EPSSetInterval(),
          EPSSetDimensions(), EPSSetEigenvalueComparison(), EPSWhich
@*/
PetscErrorCode EPSSetWhichEigenpairs(EPS eps,EPSWhich which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,which,2);
  switch (which) {
    case EPS_LARGEST_MAGNITUDE:
    case EPS_SMALLEST_MAGNITUDE:
    case EPS_LARGEST_REAL:
    case EPS_SMALLEST_REAL:
    case EPS_LARGEST_IMAGINARY:
    case EPS_SMALLEST_IMAGINARY:
    case EPS_TARGET_MAGNITUDE:
    case EPS_TARGET_REAL:
#if defined(PETSC_USE_COMPLEX)
    case EPS_TARGET_IMAGINARY:
#endif
    case EPS_ALL:
    case EPS_WHICH_USER:
      if (eps->which != which) {
        eps->state = EPS_STATE_INITIAL;
        eps->which = which;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'which' value");
  }
  PetscFunctionReturn(0);
}

/*@
   EPSGetWhichEigenpairs - Returns which portion of the spectrum is to be
   sought.

   Not Collective

   Input Parameter:
.  eps - eigensolver context obtained from EPSCreate()

   Output Parameter:
.  which - the portion of the spectrum to be sought

   Notes:
   See EPSSetWhichEigenpairs() for possible values of 'which'.

   Level: intermediate

.seealso: EPSSetWhichEigenpairs(), EPSWhich
@*/
PetscErrorCode EPSGetWhichEigenpairs(EPS eps,EPSWhich *which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(which,2);
  *which = eps->which;
  PetscFunctionReturn(0);
}

/*@C
   EPSSetEigenvalueComparison - Specifies the eigenvalue comparison function
   when EPSSetWhichEigenpairs() is set to EPS_WHICH_USER.

   Logically Collective on eps

   Input Parameters:
+  eps  - eigensolver context obtained from EPSCreate()
.  func - a pointer to the comparison function
-  ctx  - a context pointer (the last parameter to the comparison function)

   Calling Sequence of func:
$   func(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx)

+   ar     - real part of the 1st eigenvalue
.   ai     - imaginary part of the 1st eigenvalue
.   br     - real part of the 2nd eigenvalue
.   bi     - imaginary part of the 2nd eigenvalue
.   res    - result of comparison
-   ctx    - optional context, as set by EPSSetEigenvalueComparison()

   Note:
   The returning parameter 'res' can be
+  negative - if the 1st eigenvalue is preferred to the 2st one
.  zero     - if both eigenvalues are equally preferred
-  positive - if the 2st eigenvalue is preferred to the 1st one

   Level: advanced

.seealso: EPSSetWhichEigenpairs(), EPSWhich
@*/
PetscErrorCode EPSSetEigenvalueComparison(EPS eps,PetscErrorCode (*func)(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*),void* ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  eps->sc->comparison    = func;
  eps->sc->comparisonctx = ctx;
  eps->which             = EPS_WHICH_USER;
  PetscFunctionReturn(0);
}

/*@C
   EPSSetArbitrarySelection - Specifies a function intended to look for
   eigenvalues according to an arbitrary selection criterion. This criterion
   can be based on a computation involving the current eigenvector approximation.

   Logically Collective on eps

   Input Parameters:
+  eps  - eigensolver context obtained from EPSCreate()
.  func - a pointer to the evaluation function
-  ctx  - a context pointer (the last parameter to the evaluation function)

   Calling Sequence of func:
$   func(PetscScalar er,PetscScalar ei,Vec xr,Vec xi,PetscScalar *rr,PetscScalar *ri,void *ctx)

+   er     - real part of the current eigenvalue approximation
.   ei     - imaginary part of the current eigenvalue approximation
.   xr     - real part of the current eigenvector approximation
.   xi     - imaginary part of the current eigenvector approximation
.   rr     - result of evaluation (real part)
.   ri     - result of evaluation (imaginary part)
-   ctx    - optional context, as set by EPSSetArbitrarySelection()

   Notes:
   This provides a mechanism to select eigenpairs by evaluating a user-defined
   function. When a function has been provided, the default selection based on
   sorting the eigenvalues is replaced by the sorting of the results of this
   function (with the same sorting criterion given in EPSSetWhichEigenpairs()).

   For instance, suppose you want to compute those eigenvectors that maximize
   a certain computable expression. Then implement the computation using
   the arguments xr and xi, and return the result in rr. Then set the standard
   sorting by magnitude so that the eigenpair with largest value of rr is
   selected.

   This evaluation function is collective, that is, all processes call it and
   it can use collective operations; furthermore, the computed result must
   be the same in all processes.

   The result of func is expressed as a complex number so that it is possible to
   use the standard eigenvalue sorting functions, but normally only rr is used.
   Set ri to zero unless it is meaningful in your application.

   Level: advanced

.seealso: EPSSetWhichEigenpairs()
@*/
PetscErrorCode EPSSetArbitrarySelection(EPS eps,PetscErrorCode (*func)(PetscScalar,PetscScalar,Vec,Vec,PetscScalar*,PetscScalar*,void*),void* ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  eps->arbitrary    = func;
  eps->arbitraryctx = ctx;
  eps->state        = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@C
   EPSSetConvergenceTestFunction - Sets a function to compute the error estimate
   used in the convergence test.

   Logically Collective on eps

   Input Parameters:
+  eps     - eigensolver context obtained from EPSCreate()
.  func    - a pointer to the convergence test function
.  ctx     - context for private data for the convergence routine (may be null)
-  destroy - a routine for destroying the context (may be null)

   Calling Sequence of func:
$   func(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)

+   eps    - eigensolver context obtained from EPSCreate()
.   eigr   - real part of the eigenvalue
.   eigi   - imaginary part of the eigenvalue
.   res    - residual norm associated to the eigenpair
.   errest - (output) computed error estimate
-   ctx    - optional context, as set by EPSSetConvergenceTestFunction()

   Note:
   If the error estimate returned by the convergence test function is less than
   the tolerance, then the eigenvalue is accepted as converged.

   Level: advanced

.seealso: EPSSetConvergenceTest(), EPSSetTolerances()
@*/
PetscErrorCode EPSSetConvergenceTestFunction(EPS eps,PetscErrorCode (*func)(EPS,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*),void* ctx,PetscErrorCode (*destroy)(void*))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (eps->convergeddestroy) {
    ierr = (*eps->convergeddestroy)(eps->convergedctx);CHKERRQ(ierr);
  }
  eps->convergeduser    = func;
  eps->convergeddestroy = destroy;
  eps->convergedctx     = ctx;
  if (func == EPSConvergedRelative) eps->conv = EPS_CONV_REL;
  else if (func == EPSConvergedNorm) eps->conv = EPS_CONV_NORM;
  else if (func == EPSConvergedAbsolute) eps->conv = EPS_CONV_ABS;
  else {
    eps->conv      = EPS_CONV_USER;
    eps->converged = eps->convergeduser;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSSetConvergenceTest - Specifies how to compute the error estimate
   used in the convergence test.

   Logically Collective on eps

   Input Parameters:
+  eps  - eigensolver context obtained from EPSCreate()
-  conv - the type of convergence test

   Options Database Keys:
+  -eps_conv_abs  - Sets the absolute convergence test
.  -eps_conv_rel  - Sets the convergence test relative to the eigenvalue
.  -eps_conv_norm - Sets the convergence test relative to the matrix norms
-  -eps_conv_user - Selects the user-defined convergence test

   Note:
   The parameter 'conv' can have one of these values
+     EPS_CONV_ABS  - absolute error ||r||
.     EPS_CONV_REL  - error relative to the eigenvalue l, ||r||/|l|
.     EPS_CONV_NORM - error relative to the matrix norms, ||r||/(||A||+|l|*||B||)
-     EPS_CONV_USER - function set by EPSSetConvergenceTestFunction()

   Level: intermediate

.seealso: EPSGetConvergenceTest(), EPSSetConvergenceTestFunction(), EPSSetStoppingTest(), EPSConv
@*/
PetscErrorCode EPSSetConvergenceTest(EPS eps,EPSConv conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,conv,2);
  switch (conv) {
    case EPS_CONV_ABS:  eps->converged = EPSConvergedAbsolute; break;
    case EPS_CONV_REL:  eps->converged = EPSConvergedRelative; break;
    case EPS_CONV_NORM: eps->converged = EPSConvergedNorm; break;
    case EPS_CONV_USER:
      if (!eps->convergeduser) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ORDER,"Must call EPSSetConvergenceTestFunction() first");
      eps->converged = eps->convergeduser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'conv' value");
  }
  eps->conv = conv;
  PetscFunctionReturn(0);
}

/*@
   EPSGetConvergenceTest - Gets the method used to compute the error estimate
   used in the convergence test.

   Not Collective

   Input Parameters:
.  eps   - eigensolver context obtained from EPSCreate()

   Output Parameters:
.  conv  - the type of convergence test

   Level: intermediate

.seealso: EPSSetConvergenceTest(), EPSConv
@*/
PetscErrorCode EPSGetConvergenceTest(EPS eps,EPSConv *conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(conv,2);
  *conv = eps->conv;
  PetscFunctionReturn(0);
}

/*@C
   EPSSetStoppingTestFunction - Sets a function to decide when to stop the outer
   iteration of the eigensolver.

   Logically Collective on eps

   Input Parameters:
+  eps     - eigensolver context obtained from EPSCreate()
.  func    - pointer to the stopping test function
.  ctx     - context for private data for the stopping routine (may be null)
-  destroy - a routine for destroying the context (may be null)

   Calling Sequence of func:
$   func(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ctx)

+   eps    - eigensolver context obtained from EPSCreate()
.   its    - current number of iterations
.   max_it - maximum number of iterations
.   nconv  - number of currently converged eigenpairs
.   nev    - number of requested eigenpairs
.   reason - (output) result of the stopping test
-   ctx    - optional context, as set by EPSSetStoppingTestFunction()

   Note:
   Normal usage is to first call the default routine EPSStoppingBasic() and then
   set reason to EPS_CONVERGED_USER if some user-defined conditions have been
   met. To let the eigensolver continue iterating, the result must be left as
   EPS_CONVERGED_ITERATING.

   Level: advanced

.seealso: EPSSetStoppingTest(), EPSStoppingBasic()
@*/
PetscErrorCode EPSSetStoppingTestFunction(EPS eps,PetscErrorCode (*func)(EPS,PetscInt,PetscInt,PetscInt,PetscInt,EPSConvergedReason*,void*),void* ctx,PetscErrorCode (*destroy)(void*))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (eps->stoppingdestroy) {
    ierr = (*eps->stoppingdestroy)(eps->stoppingctx);CHKERRQ(ierr);
  }
  eps->stoppinguser    = func;
  eps->stoppingdestroy = destroy;
  eps->stoppingctx     = ctx;
  if (func == EPSStoppingBasic) eps->stop = EPS_STOP_BASIC;
  else {
    eps->stop     = EPS_STOP_USER;
    eps->stopping = eps->stoppinguser;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSSetStoppingTest - Specifies how to decide the termination of the outer
   loop of the eigensolver.

   Logically Collective on eps

   Input Parameters:
+  eps  - eigensolver context obtained from EPSCreate()
-  stop - the type of stopping test

   Options Database Keys:
+  -eps_stop_basic - Sets the default stopping test
-  -eps_stop_user  - Selects the user-defined stopping test

   Note:
   The parameter 'stop' can have one of these values
+     EPS_STOP_BASIC - default stopping test
-     EPS_STOP_USER  - function set by EPSSetStoppingTestFunction()

   Level: advanced

.seealso: EPSGetStoppingTest(), EPSSetStoppingTestFunction(), EPSSetConvergenceTest(), EPSStop
@*/
PetscErrorCode EPSSetStoppingTest(EPS eps,EPSStop stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,stop,2);
  switch (stop) {
    case EPS_STOP_BASIC: eps->stopping = EPSStoppingBasic; break;
    case EPS_STOP_USER:
      if (!eps->stoppinguser) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ORDER,"Must call EPSSetStoppingTestFunction() first");
      eps->stopping = eps->stoppinguser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'stop' value");
  }
  eps->stop = stop;
  PetscFunctionReturn(0);
}

/*@
   EPSGetStoppingTest - Gets the method used to decide the termination of the outer
   loop of the eigensolver.

   Not Collective

   Input Parameters:
.  eps   - eigensolver context obtained from EPSCreate()

   Output Parameters:
.  stop  - the type of stopping test

   Level: advanced

.seealso: EPSSetStoppingTest(), EPSStop
@*/
PetscErrorCode EPSGetStoppingTest(EPS eps,EPSStop *stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(stop,2);
  *stop = eps->stop;
  PetscFunctionReturn(0);
}

/*@
   EPSSetProblemType - Specifies the type of the eigenvalue problem.

   Logically Collective on eps

   Input Parameters:
+  eps      - the eigensolver context
-  type     - a known type of eigenvalue problem

   Options Database Keys:
+  -eps_hermitian - Hermitian eigenvalue problem
.  -eps_gen_hermitian - generalized Hermitian eigenvalue problem
.  -eps_non_hermitian - non-Hermitian eigenvalue problem
.  -eps_gen_non_hermitian - generalized non-Hermitian eigenvalue problem
-  -eps_pos_gen_non_hermitian - generalized non-Hermitian eigenvalue problem
   with positive semi-definite B

   Notes:
   Allowed values for the problem type are: Hermitian (EPS_HEP), non-Hermitian
   (EPS_NHEP), generalized Hermitian (EPS_GHEP), generalized non-Hermitian
   (EPS_GNHEP), generalized non-Hermitian with positive semi-definite B
   (EPS_PGNHEP), and generalized Hermitian-indefinite (EPS_GHIEP).

   This function must be used to instruct SLEPc to exploit symmetry. If no
   problem type is specified, by default a non-Hermitian problem is assumed
   (either standard or generalized). If the user knows that the problem is
   Hermitian (i.e. A=A^H) or generalized Hermitian (i.e. A=A^H, B=B^H, and
   B positive definite) then it is recommended to set the problem type so
   that eigensolver can exploit these properties.

   Level: intermediate

.seealso: EPSSetOperators(), EPSSetType(), EPSGetProblemType(), EPSProblemType
@*/
PetscErrorCode EPSSetProblemType(EPS eps,EPSProblemType type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,type,2);
  if (type == eps->problem_type) PetscFunctionReturn(0);
  switch (type) {
    case EPS_HEP:
      eps->isgeneralized = PETSC_FALSE;
      eps->ishermitian = PETSC_TRUE;
      eps->ispositive = PETSC_FALSE;
      break;
    case EPS_NHEP:
      eps->isgeneralized = PETSC_FALSE;
      eps->ishermitian = PETSC_FALSE;
      eps->ispositive = PETSC_FALSE;
      break;
    case EPS_GHEP:
      eps->isgeneralized = PETSC_TRUE;
      eps->ishermitian = PETSC_TRUE;
      eps->ispositive = PETSC_TRUE;
      break;
    case EPS_GNHEP:
      eps->isgeneralized = PETSC_TRUE;
      eps->ishermitian = PETSC_FALSE;
      eps->ispositive = PETSC_FALSE;
      break;
    case EPS_PGNHEP:
      eps->isgeneralized = PETSC_TRUE;
      eps->ishermitian = PETSC_FALSE;
      eps->ispositive = PETSC_TRUE;
      break;
    case EPS_GHIEP:
      eps->isgeneralized = PETSC_TRUE;
      eps->ishermitian = PETSC_TRUE;
      eps->ispositive = PETSC_FALSE;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_WRONG,"Unknown eigenvalue problem type");
  }
  eps->problem_type = type;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSGetProblemType - Gets the problem type from the EPS object.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  type - the problem type

   Level: intermediate

.seealso: EPSSetProblemType(), EPSProblemType
@*/
PetscErrorCode EPSGetProblemType(EPS eps,EPSProblemType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(type,2);
  *type = eps->problem_type;
  PetscFunctionReturn(0);
}

/*@
   EPSSetExtraction - Specifies the type of extraction technique to be employed
   by the eigensolver.

   Logically Collective on eps

   Input Parameters:
+  eps  - the eigensolver context
-  extr - a known type of extraction

   Options Database Keys:
+  -eps_ritz - Rayleigh-Ritz extraction
.  -eps_harmonic - harmonic Ritz extraction
.  -eps_harmonic_relative - harmonic Ritz extraction relative to the eigenvalue
.  -eps_harmonic_right - harmonic Ritz extraction for rightmost eigenvalues
.  -eps_harmonic_largest - harmonic Ritz extraction for largest magnitude
   (without target)
.  -eps_refined - refined Ritz extraction
-  -eps_refined_harmonic - refined harmonic Ritz extraction

   Notes:
   Not all eigensolvers support all types of extraction. See the SLEPc
   Users Manual for details.

   By default, a standard Rayleigh-Ritz extraction is used. Other extractions
   may be useful when computing interior eigenvalues.

   Harmonic-type extractions are used in combination with a 'target'.

   Level: advanced

.seealso: EPSSetTarget(), EPSGetExtraction(), EPSExtraction
@*/
PetscErrorCode EPSSetExtraction(EPS eps,EPSExtraction extr)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,extr,2);
  eps->extraction = extr;
  PetscFunctionReturn(0);
}

/*@
   EPSGetExtraction - Gets the extraction type used by the EPS object.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  extr - name of extraction type

   Level: advanced

.seealso: EPSSetExtraction(), EPSExtraction
@*/
PetscErrorCode EPSGetExtraction(EPS eps,EPSExtraction *extr)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(extr,2);
  *extr = eps->extraction;
  PetscFunctionReturn(0);
}

/*@
   EPSSetBalance - Specifies the balancing technique to be employed by the
   eigensolver, and some parameters associated to it.

   Logically Collective on eps

   Input Parameters:
+  eps    - the eigensolver context
.  bal    - the balancing method, one of EPS_BALANCE_NONE, EPS_BALANCE_ONESIDE,
            EPS_BALANCE_TWOSIDE, or EPS_BALANCE_USER
.  its    - number of iterations of the balancing algorithm
-  cutoff - cutoff value

   Options Database Keys:
+  -eps_balance <method> - the balancing method, where <method> is one of
                           'none', 'oneside', 'twoside', or 'user'
.  -eps_balance_its <its> - number of iterations
-  -eps_balance_cutoff <cutoff> - cutoff value

   Notes:
   When balancing is enabled, the solver works implicitly with matrix DAD^-1,
   where D is an appropriate diagonal matrix. This improves the accuracy of
   the computed results in some cases. See the SLEPc Users Manual for details.

   Balancing makes sense only for non-Hermitian problems when the required
   precision is high (i.e. a small tolerance such as 1e-15).

   By default, balancing is disabled. The two-sided method is much more
   effective than the one-sided counterpart, but it requires the system
   matrices to have the MatMultTranspose operation defined.

   The parameter 'its' is the number of iterations performed by the method. The
   cutoff value is used only in the two-side variant. Use PETSC_DEFAULT to assign
   a reasonably good value.

   User-defined balancing is allowed provided that the corresponding matrix
   is set via STSetBalanceMatrix.

   Level: intermediate

.seealso: EPSGetBalance(), EPSBalance, STSetBalanceMatrix()
@*/
PetscErrorCode EPSSetBalance(EPS eps,EPSBalance bal,PetscInt its,PetscReal cutoff)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,bal,2);
  PetscValidLogicalCollectiveInt(eps,its,3);
  PetscValidLogicalCollectiveReal(eps,cutoff,4);
  switch (bal) {
    case EPS_BALANCE_NONE:
    case EPS_BALANCE_ONESIDE:
    case EPS_BALANCE_TWOSIDE:
    case EPS_BALANCE_USER:
      if (eps->balance != bal) {
        eps->state = EPS_STATE_INITIAL;
        eps->balance = bal;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid value of argument 'bal'");
  }
  if (its==PETSC_DECIDE || its==PETSC_DEFAULT) eps->balance_its = 5;
  else {
    if (its <= 0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of its. Must be > 0");
    eps->balance_its = its;
  }
  if (cutoff==PETSC_DECIDE || cutoff==PETSC_DEFAULT) eps->balance_cutoff = 1e-8;
  else {
    if (cutoff <= 0.0) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of cutoff. Must be > 0");
    eps->balance_cutoff = cutoff;
  }
  PetscFunctionReturn(0);
}

/*@C
   EPSGetBalance - Gets the balancing type used by the EPS object, and the
   associated parameters.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameters:
+  bal    - the balancing method
.  its    - number of iterations of the balancing algorithm
-  cutoff - cutoff value

   Level: intermediate

   Note:
   The user can specify NULL for any parameter that is not needed.

.seealso: EPSSetBalance(), EPSBalance
@*/
PetscErrorCode EPSGetBalance(EPS eps,EPSBalance *bal,PetscInt *its,PetscReal *cutoff)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (bal)    *bal = eps->balance;
  if (its)    *its = eps->balance_its;
  if (cutoff) *cutoff = eps->balance_cutoff;
  PetscFunctionReturn(0);
}

/*@
   EPSSetTwoSided - Sets the solver to use a two-sided variant so that left
   eigenvectors are also computed.

   Logically Collective on eps

   Input Parameters:
+  eps      - the eigensolver context
-  twosided - whether the two-sided variant is to be used or not

   Options Database Keys:
.  -eps_two_sided <boolean> - Sets/resets the twosided flag

   Notes:
   If the user sets twosided=PETSC_TRUE then the solver uses a variant of
   the algorithm that computes both right and left eigenvectors. This is
   usually much more costly. This option is not available in all solvers.

   When using two-sided solvers, the problem matrices must have both the
   MatMult and MatMultTranspose operations defined.

   Level: advanced

.seealso: EPSGetTwoSided(), EPSGetLeftEigenvector()
@*/
PetscErrorCode EPSSetTwoSided(EPS eps,PetscBool twosided)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,twosided,2);
  if (twosided!=eps->twosided) {
    eps->twosided = twosided;
    eps->state    = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSGetTwoSided - Returns the flag indicating whether a two-sided variant
   of the algorithm is being used or not.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  twosided - the returned flag

   Level: advanced

.seealso: EPSSetTwoSided()
@*/
PetscErrorCode EPSGetTwoSided(EPS eps,PetscBool *twosided)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(twosided,2);
  *twosided = eps->twosided;
  PetscFunctionReturn(0);
}

/*@
   EPSSetTrueResidual - Specifies if the solver must compute the true residual
   explicitly or not.

   Logically Collective on eps

   Input Parameters:
+  eps     - the eigensolver context
-  trueres - whether true residuals are required or not

   Options Database Keys:
.  -eps_true_residual <boolean> - Sets/resets the boolean flag 'trueres'

   Notes:
   If the user sets trueres=PETSC_TRUE then the solver explicitly computes
   the true residual for each eigenpair approximation, and uses it for
   convergence testing. Computing the residual is usually an expensive
   operation. Some solvers (e.g., Krylov solvers) can avoid this computation
   by using a cheap estimate of the residual norm, but this may sometimes
   give inaccurate results (especially if a spectral transform is being
   used). On the contrary, preconditioned eigensolvers (e.g., Davidson solvers)
   do rely on computing the true residual, so this option is irrelevant for them.

   Level: advanced

.seealso: EPSGetTrueResidual()
@*/
PetscErrorCode EPSSetTrueResidual(EPS eps,PetscBool trueres)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,trueres,2);
  eps->trueres = trueres;
  PetscFunctionReturn(0);
}

/*@
   EPSGetTrueResidual - Returns the flag indicating whether true
   residuals must be computed explicitly or not.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  trueres - the returned flag

   Level: advanced

.seealso: EPSSetTrueResidual()
@*/
PetscErrorCode EPSGetTrueResidual(EPS eps,PetscBool *trueres)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(trueres,2);
  *trueres = eps->trueres;
  PetscFunctionReturn(0);
}

/*@
   EPSSetTrackAll - Specifies if the solver must compute the residual norm of all
   approximate eigenpairs or not.

   Logically Collective on eps

   Input Parameters:
+  eps      - the eigensolver context
-  trackall - whether to compute all residuals or not

   Notes:
   If the user sets trackall=PETSC_TRUE then the solver computes (or estimates)
   the residual norm for each eigenpair approximation. Computing the residual is
   usually an expensive operation and solvers commonly compute only the residual
   associated to the first unconverged eigenpair.

   The options '-eps_monitor_all' and '-eps_monitor_lg_all' automatically
   activate this option.

   Level: developer

.seealso: EPSGetTrackAll()
@*/
PetscErrorCode EPSSetTrackAll(EPS eps,PetscBool trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,trackall,2);
  eps->trackall = trackall;
  PetscFunctionReturn(0);
}

/*@
   EPSGetTrackAll - Returns the flag indicating whether all residual norms must
   be computed or not.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  trackall - the returned flag

   Level: developer

.seealso: EPSSetTrackAll()
@*/
PetscErrorCode EPSGetTrackAll(EPS eps,PetscBool *trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(trackall,2);
  *trackall = eps->trackall;
  PetscFunctionReturn(0);
}

/*@
   EPSSetPurify - Deactivate eigenvector purification (which is activated by default).

   Logically Collective on eps

   Input Parameters:
+  eps    - the eigensolver context
-  purify - whether purification is required or not

   Options Database Keys:
.  -eps_purify <boolean> - Sets/resets the boolean flag 'purify'

   Notes:
   By default, eigenvectors of generalized symmetric eigenproblems are purified
   in order to purge directions in the nullspace of matrix B. If the user knows
   that B is non-singular, then purification can be safely deactivated and some
   computational cost is avoided (this is particularly important in interval computations).

   Level: intermediate

.seealso: EPSGetPurify(), EPSSetInterval()
@*/
PetscErrorCode EPSSetPurify(EPS eps,PetscBool purify)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,purify,2);
  if (purify!=eps->purify) {
    eps->purify = purify;
    eps->state  = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSGetPurify - Returns the flag indicating whether purification is activated
   or not.

   Not Collective

   Input Parameter:
.  eps - the eigensolver context

   Output Parameter:
.  purify - the returned flag

   Level: intermediate

.seealso: EPSSetPurify()
@*/
PetscErrorCode EPSGetPurify(EPS eps,PetscBool *purify)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(purify,2);
  *purify = eps->purify;
  PetscFunctionReturn(0);
}

/*@C
   EPSSetOptionsPrefix - Sets the prefix used for searching for all
   EPS options in the database.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigensolver context
-  prefix - the prefix string to prepend to all EPS option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different EPS contexts, one could call
.vb
      EPSSetOptionsPrefix(eps1,"eig1_")
      EPSSetOptionsPrefix(eps2,"eig2_")
.ve

   Level: advanced

.seealso: EPSAppendOptionsPrefix(), EPSGetOptionsPrefix()
@*/
PetscErrorCode EPSSetOptionsPrefix(EPS eps,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!eps->st) { ierr = EPSGetST(eps,&eps->st);CHKERRQ(ierr); }
  ierr = STSetOptionsPrefix(eps->st,prefix);CHKERRQ(ierr);
  if (!eps->V) { ierr = EPSGetBV(eps,&eps->V);CHKERRQ(ierr); }
  ierr = BVSetOptionsPrefix(eps->V,prefix);CHKERRQ(ierr);
  if (!eps->ds) { ierr = EPSGetDS(eps,&eps->ds);CHKERRQ(ierr); }
  ierr = DSSetOptionsPrefix(eps->ds,prefix);CHKERRQ(ierr);
  if (!eps->rg) { ierr = EPSGetRG(eps,&eps->rg);CHKERRQ(ierr); }
  ierr = RGSetOptionsPrefix(eps->rg,prefix);CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)eps,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   EPSAppendOptionsPrefix - Appends to the prefix used for searching for all
   EPS options in the database.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigensolver context
-  prefix - the prefix string to prepend to all EPS option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: EPSSetOptionsPrefix(), EPSGetOptionsPrefix()
@*/
PetscErrorCode EPSAppendOptionsPrefix(EPS eps,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!eps->st) { ierr = EPSGetST(eps,&eps->st);CHKERRQ(ierr); }
  ierr = STAppendOptionsPrefix(eps->st,prefix);CHKERRQ(ierr);
  if (!eps->V) { ierr = EPSGetBV(eps,&eps->V);CHKERRQ(ierr); }
  ierr = BVAppendOptionsPrefix(eps->V,prefix);CHKERRQ(ierr);
  if (!eps->ds) { ierr = EPSGetDS(eps,&eps->ds);CHKERRQ(ierr); }
  ierr = DSAppendOptionsPrefix(eps->ds,prefix);CHKERRQ(ierr);
  if (!eps->rg) { ierr = EPSGetRG(eps,&eps->rg);CHKERRQ(ierr); }
  ierr = RGAppendOptionsPrefix(eps->rg,prefix);CHKERRQ(ierr);
  ierr = PetscObjectAppendOptionsPrefix((PetscObject)eps,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   EPSGetOptionsPrefix - Gets the prefix used for searching for all
   EPS options in the database.

   Not Collective

   Input Parameters:
.  eps - the eigensolver context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: EPSSetOptionsPrefix(), EPSAppendOptionsPrefix()
@*/
PetscErrorCode EPSGetOptionsPrefix(EPS eps,const char *prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(prefix,2);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)eps,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

