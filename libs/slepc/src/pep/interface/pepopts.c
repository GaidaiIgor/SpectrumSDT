/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   PEP routines related to options that can be set via the command-line
   or procedurally
*/

#include <slepc/private/pepimpl.h>       /*I "slepcpep.h" I*/
#include <petscdraw.h>

/*@C
   PEPMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective on pep

   Input Parameters:
+  pep      - the polynomial eigensolver context
.  name     - the monitor option name
.  help     - message indicating what monitoring is done
.  manual   - manual page for the monitor
.  monitor  - the monitor function, whose context is a PetscViewerAndFormat
-  trackall - whether this monitor tracks all eigenvalues or not

   Level: developer

.seealso: PEPMonitorSet(), PEPSetTrackAll(), PEPConvMonitorSetFromOptions()
@*/
PetscErrorCode PEPMonitorSetFromOptions(PEP pep,const char name[],const char help[],const char manual[],PetscErrorCode (*monitor)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscBool trackall)
{
  PetscErrorCode       ierr;
  PetscBool            flg;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerAndFormat *vf;

  PetscFunctionBegin;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,name,&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerAndFormatCreate(viewer,format,&vf);CHKERRQ(ierr);
    ierr = PetscObjectDereference((PetscObject)viewer);CHKERRQ(ierr);
    ierr = PEPMonitorSet(pep,(PetscErrorCode (*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))monitor,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);
    if (trackall) {
      ierr = PEPSetTrackAll(pep,PETSC_TRUE);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@C
   PEPConvMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user (for monitors that only show iteration numbers of convergence).

   Collective on pep

   Input Parameters:
+  pep      - the polynomial eigensolver context
.  name     - the monitor option name
.  help     - message indicating what monitoring is done
.  manual   - manual page for the monitor
-  monitor  - the monitor function, whose context is a SlepcConvMonitor

   Level: developer

.seealso: PEPMonitorSet(), PEPMonitorSetFromOptions()
@*/
PetscErrorCode PEPConvMonitorSetFromOptions(PEP pep,const char name[],const char help[],const char manual[],PetscErrorCode (*monitor)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,SlepcConvMonitor))
{
  PetscErrorCode    ierr;
  PetscBool         flg;
  PetscViewer       viewer;
  PetscViewerFormat format;
  SlepcConvMonitor  ctx;

  PetscFunctionBegin;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,name,&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = SlepcConvMonitorCreate(viewer,format,&ctx);CHKERRQ(ierr);
    ierr = PetscObjectDereference((PetscObject)viewer);CHKERRQ(ierr);
    ierr = PEPMonitorSet(pep,(PetscErrorCode (*)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))monitor,ctx,(PetscErrorCode (*)(void**))SlepcConvMonitorDestroy);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PEPSetFromOptions - Sets PEP options from the options database.
   This routine must be called before PEPSetUp() if the user is to be
   allowed to set the solver type.

   Collective on pep

   Input Parameters:
.  pep - the polynomial eigensolver context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner
@*/
PetscErrorCode PEPSetFromOptions(PEP pep)
{
  PetscErrorCode  ierr;
  char            type[256];
  PetscBool       set,flg,flg1,flg2,flg3,flg4,flg5;
  PetscReal       r,t,array[2]={0,0};
  PetscScalar     s;
  PetscInt        i,j,k;
  PetscDrawLG     lg;
  PEPScale        scale;
  PEPRefine       refine;
  PEPRefineScheme scheme;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  ierr = PEPRegisterAll();CHKERRQ(ierr);
  ierr = PetscObjectOptionsBegin((PetscObject)pep);CHKERRQ(ierr);
    ierr = PetscOptionsFList("-pep_type","Polynomial eigensolver method","PEPSetType",PEPList,(char*)(((PetscObject)pep)->type_name?((PetscObject)pep)->type_name:PEPTOAR),type,256,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = PEPSetType(pep,type);CHKERRQ(ierr);
    } else if (!((PetscObject)pep)->type_name) {
      ierr = PEPSetType(pep,PEPTOAR);CHKERRQ(ierr);
    }

    ierr = PetscOptionsBoolGroupBegin("-pep_general","General polynomial eigenvalue problem","PEPSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetProblemType(pep,PEP_GENERAL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_hermitian","Hermitian polynomial eigenvalue problem","PEPSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetProblemType(pep,PEP_HERMITIAN);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_hyperbolic","Hyperbolic polynomial eigenvalue problem","PEPSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetProblemType(pep,PEP_HYPERBOLIC);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-pep_gyroscopic","Gyroscopic polynomial eigenvalue problem","PEPSetProblemType",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetProblemType(pep,PEP_GYROSCOPIC);CHKERRQ(ierr); }

    scale = pep->scale;
    ierr = PetscOptionsEnum("-pep_scale","Scaling strategy","PEPSetScale",PEPScaleTypes,(PetscEnum)scale,(PetscEnum*)&scale,&flg1);CHKERRQ(ierr);
    r = pep->sfactor;
    ierr = PetscOptionsReal("-pep_scale_factor","Scale factor","PEPSetScale",pep->sfactor,&r,&flg2);CHKERRQ(ierr);
    if (!flg2 && r==1.0) r = PETSC_DEFAULT;
    j = pep->sits;
    ierr = PetscOptionsInt("-pep_scale_its","Number of iterations in diagonal scaling","PEPSetScale",pep->sits,&j,&flg3);CHKERRQ(ierr);
    t = pep->slambda;
    ierr = PetscOptionsReal("-pep_scale_lambda","Estimate of eigenvalue (modulus) for diagonal scaling","PEPSetScale",pep->slambda,&t,&flg4);CHKERRQ(ierr);
    if (flg1 || flg2 || flg3 || flg4) { ierr = PEPSetScale(pep,scale,r,NULL,NULL,j,t);CHKERRQ(ierr); }

    ierr = PetscOptionsEnum("-pep_extract","Extraction method","PEPSetExtract",PEPExtractTypes,(PetscEnum)pep->extract,(PetscEnum*)&pep->extract,NULL);CHKERRQ(ierr);

    refine = pep->refine;
    ierr = PetscOptionsEnum("-pep_refine","Iterative refinement method","PEPSetRefine",PEPRefineTypes,(PetscEnum)refine,(PetscEnum*)&refine,&flg1);CHKERRQ(ierr);
    i = pep->npart;
    ierr = PetscOptionsInt("-pep_refine_partitions","Number of partitions of the communicator for iterative refinement","PEPSetRefine",pep->npart,&i,&flg2);CHKERRQ(ierr);
    r = pep->rtol;
    ierr = PetscOptionsReal("-pep_refine_tol","Tolerance for iterative refinement","PEPSetRefine",pep->rtol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL/1000:pep->rtol,&r,&flg3);CHKERRQ(ierr);
    j = pep->rits;
    ierr = PetscOptionsInt("-pep_refine_its","Maximum number of iterations for iterative refinement","PEPSetRefine",pep->rits,&j,&flg4);CHKERRQ(ierr);
    scheme = pep->scheme;
    ierr = PetscOptionsEnum("-pep_refine_scheme","Scheme used for linear systems within iterative refinement","PEPSetRefine",PEPRefineSchemes,(PetscEnum)scheme,(PetscEnum*)&scheme,&flg5);CHKERRQ(ierr);
    if (flg1 || flg2 || flg3 || flg4 || flg5) { ierr = PEPSetRefine(pep,refine,i,r,j,scheme);CHKERRQ(ierr); }

    i = pep->max_it;
    ierr = PetscOptionsInt("-pep_max_it","Maximum number of iterations","PEPSetTolerances",pep->max_it,&i,&flg1);CHKERRQ(ierr);
    r = pep->tol;
    ierr = PetscOptionsReal("-pep_tol","Tolerance","PEPSetTolerances",pep->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL:pep->tol,&r,&flg2);CHKERRQ(ierr);
    if (flg1 || flg2) { ierr = PEPSetTolerances(pep,r,i);CHKERRQ(ierr); }

    ierr = PetscOptionsBoolGroupBegin("-pep_conv_rel","Relative error convergence test","PEPSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetConvergenceTest(pep,PEP_CONV_REL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_conv_norm","Convergence test relative to the matrix norms","PEPSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetConvergenceTest(pep,PEP_CONV_NORM);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_conv_abs","Absolute error convergence test","PEPSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetConvergenceTest(pep,PEP_CONV_ABS);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-pep_conv_user","User-defined convergence test","PEPSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetConvergenceTest(pep,PEP_CONV_USER);CHKERRQ(ierr); }

    ierr = PetscOptionsBoolGroupBegin("-pep_stop_basic","Stop iteration if all eigenvalues converged or max_it reached","PEPSetStoppingTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetStoppingTest(pep,PEP_STOP_BASIC);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-pep_stop_user","User-defined stopping test","PEPSetStoppingTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetStoppingTest(pep,PEP_STOP_USER);CHKERRQ(ierr); }

    i = pep->nev;
    ierr = PetscOptionsInt("-pep_nev","Number of eigenvalues to compute","PEPSetDimensions",pep->nev,&i,&flg1);CHKERRQ(ierr);
    j = pep->ncv;
    ierr = PetscOptionsInt("-pep_ncv","Number of basis vectors","PEPSetDimensions",pep->ncv,&j,&flg2);CHKERRQ(ierr);
    k = pep->mpd;
    ierr = PetscOptionsInt("-pep_mpd","Maximum dimension of projected problem","PEPSetDimensions",pep->mpd,&k,&flg3);CHKERRQ(ierr);
    if (flg1 || flg2 || flg3) { ierr = PEPSetDimensions(pep,i,j,k);CHKERRQ(ierr); }

    ierr = PetscOptionsEnum("-pep_basis","Polynomial basis","PEPSetBasis",PEPBasisTypes,(PetscEnum)pep->basis,(PetscEnum*)&pep->basis,NULL);CHKERRQ(ierr);

    ierr = PetscOptionsBoolGroupBegin("-pep_largest_magnitude","Compute largest eigenvalues in magnitude","PEPSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetWhichEigenpairs(pep,PEP_LARGEST_MAGNITUDE);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_smallest_magnitude","Compute smallest eigenvalues in magnitude","PEPSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetWhichEigenpairs(pep,PEP_SMALLEST_MAGNITUDE);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_largest_real","Compute eigenvalues with largest real parts","PEPSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetWhichEigenpairs(pep,PEP_LARGEST_REAL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_smallest_real","Compute eigenvalues with smallest real parts","PEPSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetWhichEigenpairs(pep,PEP_SMALLEST_REAL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_largest_imaginary","Compute eigenvalues with largest imaginary parts","PEPSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetWhichEigenpairs(pep,PEP_LARGEST_IMAGINARY);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_smallest_imaginary","Compute eigenvalues with smallest imaginary parts","PEPSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetWhichEigenpairs(pep,PEP_SMALLEST_IMAGINARY);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_target_magnitude","Compute eigenvalues closest to target","PEPSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetWhichEigenpairs(pep,PEP_TARGET_MAGNITUDE);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-pep_target_real","Compute eigenvalues with real parts closest to target","PEPSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetWhichEigenpairs(pep,PEP_TARGET_REAL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-pep_target_imaginary","Compute eigenvalues with imaginary parts closest to target","PEPSetWhichEigenpairs",&flg);CHKERRQ(ierr);
    if (flg) { ierr = PEPSetWhichEigenpairs(pep,PEP_TARGET_IMAGINARY);CHKERRQ(ierr); }

    ierr = PetscOptionsScalar("-pep_target","Value of the target","PEPSetTarget",pep->target,&s,&flg);CHKERRQ(ierr);
    if (flg) {
      if (pep->which!=PEP_TARGET_REAL && pep->which!=PEP_TARGET_IMAGINARY) {
        ierr = PEPSetWhichEigenpairs(pep,PEP_TARGET_MAGNITUDE);CHKERRQ(ierr);
      }
      ierr = PEPSetTarget(pep,s);CHKERRQ(ierr);
    }

    k = 2;
    ierr = PetscOptionsRealArray("-pep_interval","Computational interval (two real values separated with a comma without spaces)","PEPSetInterval",array,&k,&flg);CHKERRQ(ierr);
    if (flg) {
      if (k<2) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_SIZ,"Must pass two values in -pep_interval (comma-separated without spaces)");
      ierr = PEPSetWhichEigenpairs(pep,PEP_ALL);CHKERRQ(ierr);
      ierr = PEPSetInterval(pep,array[0],array[1]);CHKERRQ(ierr);
    }

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to PEPSetFromOptions()
    */
    ierr = PetscOptionsBool("-pep_monitor_cancel","Remove any hardwired monitor routines","PEPMonitorCancel",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = PEPMonitorCancel(pep);CHKERRQ(ierr);
    }
    /*
      Text monitors
    */
    ierr = PEPMonitorSetFromOptions(pep,"-pep_monitor","Monitor first unconverged approximate eigenvalue and error estimate","PEPMonitorFirst",PEPMonitorFirst,PETSC_FALSE);CHKERRQ(ierr);
    ierr = PEPConvMonitorSetFromOptions(pep,"-pep_monitor_conv","Monitor approximate eigenvalues and error estimates as they converge","PEPMonitorConverged",PEPMonitorConverged);CHKERRQ(ierr);
    ierr = PEPMonitorSetFromOptions(pep,"-pep_monitor_all","Monitor approximate eigenvalues and error estimates","PEPMonitorAll",PEPMonitorAll,PETSC_TRUE);CHKERRQ(ierr);
    /*
      Line graph monitors
    */
    ierr = PetscOptionsBool("-pep_monitor_lg","Monitor first unconverged approximate error estimate graphically","PEPMonitorSet",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = PEPMonitorLGCreate(PetscObjectComm((PetscObject)pep),NULL,"Error estimates",PETSC_DECIDE,PETSC_DECIDE,300,300,&lg);CHKERRQ(ierr);
      ierr = PEPMonitorSet(pep,PEPMonitorLG,lg,(PetscErrorCode (*)(void**))PetscDrawLGDestroy);CHKERRQ(ierr);
    }
    ierr = PetscOptionsBool("-pep_monitor_lg_all","Monitor error estimates graphically","PEPMonitorSet",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = PEPMonitorLGCreate(PetscObjectComm((PetscObject)pep),NULL,"Error estimates",PETSC_DECIDE,PETSC_DECIDE,300,300,&lg);CHKERRQ(ierr);
      ierr = PEPMonitorSet(pep,PEPMonitorLGAll,lg,(PetscErrorCode (*)(void**))PetscDrawLGDestroy);CHKERRQ(ierr);
      ierr = PEPSetTrackAll(pep,PETSC_TRUE);CHKERRQ(ierr);
    }

    /* -----------------------------------------------------------------------*/
    ierr = PetscOptionsName("-pep_view","Print detailed information on solver used","PEPView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-pep_view_vectors","View computed eigenvectors","PEPVectorsView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-pep_view_values","View computed eigenvalues","PEPValuesView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-pep_converged_reason","Print reason for convergence, and number of iterations","PEPReasonView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-pep_error_absolute","Print absolute errors of each eigenpair","PEPErrorView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-pep_error_relative","Print relative errors of each eigenpair","PEPErrorView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-pep_error_backward","Print backward errors of each eigenpair","PEPErrorView",NULL);CHKERRQ(ierr);

    if (pep->ops->setfromoptions) {
      ierr = (*pep->ops->setfromoptions)(PetscOptionsObject,pep);CHKERRQ(ierr);
    }
    ierr = PetscObjectProcessOptionsHandlers(PetscOptionsObject,(PetscObject)pep);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  if (!pep->V) { ierr = PEPGetBV(pep,&pep->V);CHKERRQ(ierr); }
  ierr = BVSetFromOptions(pep->V);CHKERRQ(ierr);
  if (!pep->rg) { ierr = PEPGetRG(pep,&pep->rg);CHKERRQ(ierr); }
  ierr = RGSetFromOptions(pep->rg);CHKERRQ(ierr);
  if (!pep->ds) { ierr = PEPGetDS(pep,&pep->ds);CHKERRQ(ierr); }
  ierr = DSSetFromOptions(pep->ds);CHKERRQ(ierr);
  if (!pep->st) { ierr = PEPGetST(pep,&pep->st);CHKERRQ(ierr); }
  ierr = PEPSetDefaultST(pep);CHKERRQ(ierr);
  ierr = STSetFromOptions(pep->st);CHKERRQ(ierr);
  if (!pep->refineksp) { ierr = PEPRefineGetKSP(pep,&pep->refineksp);CHKERRQ(ierr); }
  ierr = KSPSetFromOptions(pep->refineksp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PEPGetTolerances - Gets the tolerance and maximum iteration count used
   by the PEP convergence tests.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: PEPSetTolerances()
@*/
PetscErrorCode PEPGetTolerances(PEP pep,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (tol)    *tol    = pep->tol;
  if (maxits) *maxits = pep->max_it;
  PetscFunctionReturn(0);
}

/*@
   PEPSetTolerances - Sets the tolerance and maximum iteration count used
   by the PEP convergence tests.

   Logically Collective on pep

   Input Parameters:
+  pep - the polynomial eigensolver context
.  tol - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -pep_tol <tol> - Sets the convergence tolerance
-  -pep_max_it <maxits> - Sets the maximum number of iterations allowed

   Notes:
   Use PETSC_DEFAULT for either argument to assign a reasonably good value.

   Level: intermediate

.seealso: PEPGetTolerances()
@*/
PetscErrorCode PEPSetTolerances(PEP pep,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(pep,tol,2);
  PetscValidLogicalCollectiveInt(pep,maxits,3);
  if (tol == PETSC_DEFAULT) {
    pep->tol   = PETSC_DEFAULT;
    pep->state = PEP_STATE_INITIAL;
  } else {
    if (tol <= 0.0) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    pep->tol = tol;
  }
  if (maxits == PETSC_DEFAULT || maxits == PETSC_DECIDE) {
    pep->max_it = PETSC_DEFAULT;
    pep->state  = PEP_STATE_INITIAL;
  } else {
    if (maxits <= 0) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    pep->max_it = maxits;
  }
  PetscFunctionReturn(0);
}

/*@C
   PEPGetDimensions - Gets the number of eigenvalues to compute
   and the dimension of the subspace.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: PEPSetDimensions()
@*/
PetscErrorCode PEPGetDimensions(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (nev) *nev = pep->nev;
  if (ncv) *ncv = pep->ncv;
  if (mpd) *mpd = pep->mpd;
  PetscFunctionReturn(0);
}

/*@
   PEPSetDimensions - Sets the number of eigenvalues to compute
   and the dimension of the subspace.

   Logically Collective on pep

   Input Parameters:
+  pep - the polynomial eigensolver context
.  nev - number of eigenvalues to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Keys:
+  -pep_nev <nev> - Sets the number of eigenvalues
.  -pep_ncv <ncv> - Sets the dimension of the subspace
-  -pep_mpd <mpd> - Sets the maximum projected dimension

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

   When computing all eigenvalues in an interval, see PEPSetInterval(), these
   parameters lose relevance, and tuning must be done with PEPSTOARSetDimensions().

   Level: intermediate

.seealso: PEPGetDimensions(), PEPSetInterval(), PEPSTOARSetDimensions()
@*/
PetscErrorCode PEPSetDimensions(PEP pep,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveInt(pep,nev,2);
  PetscValidLogicalCollectiveInt(pep,ncv,3);
  PetscValidLogicalCollectiveInt(pep,mpd,4);
  if (nev<1) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nev. Must be > 0");
  pep->nev = nev;
  if (ncv == PETSC_DECIDE || ncv == PETSC_DEFAULT) {
    pep->ncv = PETSC_DEFAULT;
  } else {
    if (ncv<1) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    pep->ncv = ncv;
  }
  if (mpd == PETSC_DECIDE || mpd == PETSC_DEFAULT) {
    pep->mpd = PETSC_DEFAULT;
  } else {
    if (mpd<1) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    pep->mpd = mpd;
  }
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   PEPSetWhichEigenpairs - Specifies which portion of the spectrum is
   to be sought.

   Logically Collective on pep

   Input Parameters:
+  pep   - eigensolver context obtained from PEPCreate()
-  which - the portion of the spectrum to be sought

   Possible values:
   The parameter 'which' can have one of these values

+     PEP_LARGEST_MAGNITUDE - largest eigenvalues in magnitude (default)
.     PEP_SMALLEST_MAGNITUDE - smallest eigenvalues in magnitude
.     PEP_LARGEST_REAL - largest real parts
.     PEP_SMALLEST_REAL - smallest real parts
.     PEP_LARGEST_IMAGINARY - largest imaginary parts
.     PEP_SMALLEST_IMAGINARY - smallest imaginary parts
.     PEP_TARGET_MAGNITUDE - eigenvalues closest to the target (in magnitude)
.     PEP_TARGET_REAL - eigenvalues with real part closest to target
.     PEP_TARGET_IMAGINARY - eigenvalues with imaginary part closest to target
.     PEP_ALL - all eigenvalues contained in a given interval
-     PEP_WHICH_USER - user defined ordering set with PEPSetEigenvalueComparison()

   Options Database Keys:
+   -pep_largest_magnitude - Sets largest eigenvalues in magnitude
.   -pep_smallest_magnitude - Sets smallest eigenvalues in magnitude
.   -pep_largest_real - Sets largest real parts
.   -pep_smallest_real - Sets smallest real parts
.   -pep_largest_imaginary - Sets largest imaginary parts
.   -pep_smallest_imaginary - Sets smallest imaginary parts
.   -pep_target_magnitude - Sets eigenvalues closest to target
.   -pep_target_real - Sets real parts closest to target
.   -pep_target_imaginary - Sets imaginary parts closest to target
-   -pep_all - Sets all eigenvalues in an interval

   Notes:
   Not all eigensolvers implemented in PEP account for all the possible values
   stated above. If SLEPc is compiled for real numbers PEP_LARGEST_IMAGINARY
   and PEP_SMALLEST_IMAGINARY use the absolute value of the imaginary part
   for eigenvalue selection.

   The target is a scalar value provided with PEPSetTarget().

   The criterion PEP_TARGET_IMAGINARY is available only in case PETSc and
   SLEPc have been built with complex scalars.

   PEP_ALL is intended for use in combination with an interval (see
   PEPSetInterval()), when all eigenvalues within the interval are requested.
   In that case, the number of eigenvalues is unknown, so the nev parameter
   has a different sense, see PEPSetDimensions().

   Level: intermediate

.seealso: PEPGetWhichEigenpairs(), PEPSetTarget(), PEPSetInterval(),
          PEPSetDimensions(), PEPSetEigenvalueComparison(), PEPWhich
@*/
PetscErrorCode PEPSetWhichEigenpairs(PEP pep,PEPWhich which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,which,2);
  switch (which) {
    case PEP_LARGEST_MAGNITUDE:
    case PEP_SMALLEST_MAGNITUDE:
    case PEP_LARGEST_REAL:
    case PEP_SMALLEST_REAL:
    case PEP_LARGEST_IMAGINARY:
    case PEP_SMALLEST_IMAGINARY:
    case PEP_TARGET_MAGNITUDE:
    case PEP_TARGET_REAL:
#if defined(PETSC_USE_COMPLEX)
    case PEP_TARGET_IMAGINARY:
#endif
    case PEP_ALL:
    case PEP_WHICH_USER:
      if (pep->which != which) {
        pep->state = PEP_STATE_INITIAL;
        pep->which = which;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'which' value");
  }
  PetscFunctionReturn(0);
}

/*@
    PEPGetWhichEigenpairs - Returns which portion of the spectrum is to be
    sought.

    Not Collective

    Input Parameter:
.   pep - eigensolver context obtained from PEPCreate()

    Output Parameter:
.   which - the portion of the spectrum to be sought

    Notes:
    See PEPSetWhichEigenpairs() for possible values of 'which'.

    Level: intermediate

.seealso: PEPSetWhichEigenpairs(), PEPWhich
@*/
PetscErrorCode PEPGetWhichEigenpairs(PEP pep,PEPWhich *which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(which,2);
  *which = pep->which;
  PetscFunctionReturn(0);
}

/*@C
   PEPSetEigenvalueComparison - Specifies the eigenvalue comparison function
   when PEPSetWhichEigenpairs() is set to PEP_WHICH_USER.

   Logically Collective on pep

   Input Parameters:
+  pep  - eigensolver context obtained from PEPCreate()
.  func - a pointer to the comparison function
-  ctx  - a context pointer (the last parameter to the comparison function)

   Calling Sequence of func:
$   func(PetscScalar ar,PetscScalar ai,PetscScalar br,PetscScalar bi,PetscInt *res,void *ctx)

+   ar     - real part of the 1st eigenvalue
.   ai     - imaginary part of the 1st eigenvalue
.   br     - real part of the 2nd eigenvalue
.   bi     - imaginary part of the 2nd eigenvalue
.   res    - result of comparison
-   ctx    - optional context, as set by PEPSetEigenvalueComparison()

   Note:
   The returning parameter 'res' can be
+  negative - if the 1st eigenvalue is preferred to the 2st one
.  zero     - if both eigenvalues are equally preferred
-  positive - if the 2st eigenvalue is preferred to the 1st one

   Level: advanced

.seealso: PEPSetWhichEigenpairs(), PEPWhich
@*/
PetscErrorCode PEPSetEigenvalueComparison(PEP pep,PetscErrorCode (*func)(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt*,void*),void* ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  pep->sc->comparison    = func;
  pep->sc->comparisonctx = ctx;
  pep->which             = PEP_WHICH_USER;
  PetscFunctionReturn(0);
}

/*@
   PEPSetProblemType - Specifies the type of the polynomial eigenvalue problem.

   Logically Collective on pep

   Input Parameters:
+  pep  - the polynomial eigensolver context
-  type - a known type of polynomial eigenvalue problem

   Options Database Keys:
+  -pep_general - general problem with no particular structure
.  -pep_hermitian - problem whose coefficient matrices are Hermitian
.  -pep_hyperbolic - Hermitian problem that satisfies the definition of hyperbolic
-  -pep_gyroscopic - problem with Hamiltonian structure

   Notes:
   Allowed values for the problem type are: general (PEP_GENERAL), Hermitian
   (PEP_HERMITIAN), hyperbolic (PEP_HYPERBOLIC), and gyroscopic (PEP_GYROSCOPIC).

   This function is used to instruct SLEPc to exploit certain structure in
   the polynomial eigenproblem. By default, no particular structure is assumed.

   If the problem matrices are Hermitian (symmetric in the real case) or
   Hermitian/skew-Hermitian then the solver can exploit this fact to perform
   less operations or provide better stability. Hyperbolic problems are a
   particular case of Hermitian problems, some solvers may treat them simply as
   Hermitian.

   Level: intermediate

.seealso: PEPSetOperators(), PEPSetType(), PEPGetProblemType(), PEPProblemType
@*/
PetscErrorCode PEPSetProblemType(PEP pep,PEPProblemType type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,type,2);
  if (type!=PEP_GENERAL && type!=PEP_HERMITIAN && type!=PEP_HYPERBOLIC && type!=PEP_GYROSCOPIC) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"Unknown eigenvalue problem type");
  if (type != pep->problem_type) {
    pep->problem_type = type;
    pep->state = PEP_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   PEPGetProblemType - Gets the problem type from the PEP object.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameter:
.  type - the problem type

   Level: intermediate

.seealso: PEPSetProblemType(), PEPProblemType
@*/
PetscErrorCode PEPGetProblemType(PEP pep,PEPProblemType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(type,2);
  *type = pep->problem_type;
  PetscFunctionReturn(0);
}

/*@
   PEPSetBasis - Specifies the type of polynomial basis used to describe the
   polynomial eigenvalue problem.

   Logically Collective on pep

   Input Parameters:
+  pep   - the polynomial eigensolver context
-  basis - the type of polynomial basis

   Options Database Key:
.  -pep_basis <basis> - Select the basis type

   Notes:
   By default, the coefficient matrices passed via PEPSetOperators() are
   expressed in the monomial basis, i.e.
   P(lambda) = A_0 + lambda*A_1 + lambda^2*A_2 + ... + lambda^d*A_d.
   Other polynomial bases may have better numerical behaviour, but the user
   must then pass the coefficient matrices accordingly.

   Level: intermediate

.seealso: PEPSetOperators(), PEPGetBasis(), PEPBasis
@*/
PetscErrorCode PEPSetBasis(PEP pep,PEPBasis basis)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,basis,2);
  pep->basis = basis;
  PetscFunctionReturn(0);
}

/*@
   PEPGetBasis - Gets the type of polynomial basis from the PEP object.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameter:
.  basis - the polynomial basis

   Level: intermediate

.seealso: PEPSetBasis(), PEPBasis
@*/
PetscErrorCode PEPGetBasis(PEP pep,PEPBasis *basis)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(basis,2);
  *basis = pep->basis;
  PetscFunctionReturn(0);
}

/*@
   PEPSetTrackAll - Specifies if the solver must compute the residual of all
   approximate eigenpairs or not.

   Logically Collective on pep

   Input Parameters:
+  pep      - the eigensolver context
-  trackall - whether compute all residuals or not

   Notes:
   If the user sets trackall=PETSC_TRUE then the solver explicitly computes
   the residual for each eigenpair approximation. Computing the residual is
   usually an expensive operation and solvers commonly compute the associated
   residual to the first unconverged eigenpair.
   The options '-pep_monitor_all' and '-pep_monitor_lg_all' automatically
   activate this option.

   Level: developer

.seealso: PEPGetTrackAll()
@*/
PetscErrorCode PEPSetTrackAll(PEP pep,PetscBool trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveBool(pep,trackall,2);
  pep->trackall = trackall;
  PetscFunctionReturn(0);
}

/*@
   PEPGetTrackAll - Returns the flag indicating whether all residual norms must
   be computed or not.

   Not Collective

   Input Parameter:
.  pep - the eigensolver context

   Output Parameter:
.  trackall - the returned flag

   Level: developer

.seealso: PEPSetTrackAll()
@*/
PetscErrorCode PEPGetTrackAll(PEP pep,PetscBool *trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidBoolPointer(trackall,2);
  *trackall = pep->trackall;
  PetscFunctionReturn(0);
}

/*@C
   PEPSetConvergenceTestFunction - Sets a function to compute the error estimate
   used in the convergence test.

   Logically Collective on pep

   Input Parameters:
+  pep     - eigensolver context obtained from PEPCreate()
.  func    - a pointer to the convergence test function
.  ctx     - context for private data for the convergence routine (may be null)
-  destroy - a routine for destroying the context (may be null)

   Calling Sequence of func:
$   func(PEP pep,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)

+   pep    - eigensolver context obtained from PEPCreate()
.   eigr   - real part of the eigenvalue
.   eigi   - imaginary part of the eigenvalue
.   res    - residual norm associated to the eigenpair
.   errest - (output) computed error estimate
-   ctx    - optional context, as set by PEPSetConvergenceTestFunction()

   Note:
   If the error estimate returned by the convergence test function is less than
   the tolerance, then the eigenvalue is accepted as converged.

   Level: advanced

.seealso: PEPSetConvergenceTest(), PEPSetTolerances()
@*/
PetscErrorCode PEPSetConvergenceTestFunction(PEP pep,PetscErrorCode (*func)(PEP,PetscScalar,PetscScalar,PetscReal,PetscReal*,void*),void* ctx,PetscErrorCode (*destroy)(void*))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (pep->convergeddestroy) {
    ierr = (*pep->convergeddestroy)(pep->convergedctx);CHKERRQ(ierr);
  }
  pep->convergeduser    = func;
  pep->convergeddestroy = destroy;
  pep->convergedctx     = ctx;
  if (func == PEPConvergedRelative) pep->conv = PEP_CONV_REL;
  else if (func == PEPConvergedNorm) pep->conv = PEP_CONV_NORM;
  else if (func == PEPConvergedAbsolute) pep->conv = PEP_CONV_ABS;
  else {
    pep->conv      = PEP_CONV_USER;
    pep->converged = pep->convergeduser;
  }
  PetscFunctionReturn(0);
}

/*@
   PEPSetConvergenceTest - Specifies how to compute the error estimate
   used in the convergence test.

   Logically Collective on pep

   Input Parameters:
+  pep  - eigensolver context obtained from PEPCreate()
-  conv - the type of convergence test

   Options Database Keys:
+  -pep_conv_abs    - Sets the absolute convergence test
.  -pep_conv_rel    - Sets the convergence test relative to the eigenvalue
.  -pep_conv_norm   - Sets the convergence test relative to the matrix norms
-  -pep_conv_user   - Selects the user-defined convergence test

   Note:
   The parameter 'conv' can have one of these values
+     PEP_CONV_ABS    - absolute error ||r||
.     PEP_CONV_REL    - error relative to the eigenvalue l, ||r||/|l|
.     PEP_CONV_NORM   - error relative matrix norms, ||r||/sum_i(l^i*||A_i||)
-     PEP_CONV_USER   - function set by PEPSetConvergenceTestFunction()

   Level: intermediate

.seealso: PEPGetConvergenceTest(), PEPSetConvergenceTestFunction(), PEPSetStoppingTest(), PEPConv
@*/
PetscErrorCode PEPSetConvergenceTest(PEP pep,PEPConv conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,conv,2);
  switch (conv) {
    case PEP_CONV_ABS:  pep->converged = PEPConvergedAbsolute; break;
    case PEP_CONV_REL:  pep->converged = PEPConvergedRelative; break;
    case PEP_CONV_NORM: pep->converged = PEPConvergedNorm; break;
    case PEP_CONV_USER:
      if (!pep->convergeduser) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ORDER,"Must call PEPSetConvergenceTestFunction() first");
      pep->converged = pep->convergeduser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'conv' value");
  }
  pep->conv = conv;
  PetscFunctionReturn(0);
}

/*@
   PEPGetConvergenceTest - Gets the method used to compute the error estimate
   used in the convergence test.

   Not Collective

   Input Parameters:
.  pep   - eigensolver context obtained from PEPCreate()

   Output Parameters:
.  conv  - the type of convergence test

   Level: intermediate

.seealso: PEPSetConvergenceTest(), PEPConv
@*/
PetscErrorCode PEPGetConvergenceTest(PEP pep,PEPConv *conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(conv,2);
  *conv = pep->conv;
  PetscFunctionReturn(0);
}

/*@C
   PEPSetStoppingTestFunction - Sets a function to decide when to stop the outer
   iteration of the eigensolver.

   Logically Collective on pep

   Input Parameters:
+  pep     - eigensolver context obtained from PEPCreate()
.  func    - pointer to the stopping test function
.  ctx     - context for private data for the stopping routine (may be null)
-  destroy - a routine for destroying the context (may be null)

   Calling Sequence of func:
$   func(PEP pep,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,PEPConvergedReason *reason,void *ctx)

+   pep    - eigensolver context obtained from PEPCreate()
.   its    - current number of iterations
.   max_it - maximum number of iterations
.   nconv  - number of currently converged eigenpairs
.   nev    - number of requested eigenpairs
.   reason - (output) result of the stopping test
-   ctx    - optional context, as set by PEPSetStoppingTestFunction()

   Note:
   Normal usage is to first call the default routine PEPStoppingBasic() and then
   set reason to PEP_CONVERGED_USER if some user-defined conditions have been
   met. To let the eigensolver continue iterating, the result must be left as
   PEP_CONVERGED_ITERATING.

   Level: advanced

.seealso: PEPSetStoppingTest(), PEPStoppingBasic()
@*/
PetscErrorCode PEPSetStoppingTestFunction(PEP pep,PetscErrorCode (*func)(PEP,PetscInt,PetscInt,PetscInt,PetscInt,PEPConvergedReason*,void*),void* ctx,PetscErrorCode (*destroy)(void*))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (pep->stoppingdestroy) {
    ierr = (*pep->stoppingdestroy)(pep->stoppingctx);CHKERRQ(ierr);
  }
  pep->stoppinguser    = func;
  pep->stoppingdestroy = destroy;
  pep->stoppingctx     = ctx;
  if (func == PEPStoppingBasic) pep->stop = PEP_STOP_BASIC;
  else {
    pep->stop     = PEP_STOP_USER;
    pep->stopping = pep->stoppinguser;
  }
  PetscFunctionReturn(0);
}

/*@
   PEPSetStoppingTest - Specifies how to decide the termination of the outer
   loop of the eigensolver.

   Logically Collective on pep

   Input Parameters:
+  pep  - eigensolver context obtained from PEPCreate()
-  stop - the type of stopping test

   Options Database Keys:
+  -pep_stop_basic - Sets the default stopping test
-  -pep_stop_user  - Selects the user-defined stopping test

   Note:
   The parameter 'stop' can have one of these values
+     PEP_STOP_BASIC - default stopping test
-     PEP_STOP_USER  - function set by PEPSetStoppingTestFunction()

   Level: advanced

.seealso: PEPGetStoppingTest(), PEPSetStoppingTestFunction(), PEPSetConvergenceTest(), PEPStop
@*/
PetscErrorCode PEPSetStoppingTest(PEP pep,PEPStop stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,stop,2);
  switch (stop) {
    case PEP_STOP_BASIC: pep->stopping = PEPStoppingBasic; break;
    case PEP_STOP_USER:
      if (!pep->stoppinguser) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ORDER,"Must call PEPSetStoppingTestFunction() first");
      pep->stopping = pep->stoppinguser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'stop' value");
  }
  pep->stop = stop;
  PetscFunctionReturn(0);
}

/*@
   PEPGetStoppingTest - Gets the method used to decide the termination of the outer
   loop of the eigensolver.

   Not Collective

   Input Parameters:
.  pep   - eigensolver context obtained from PEPCreate()

   Output Parameters:
.  stop  - the type of stopping test

   Level: advanced

.seealso: PEPSetStoppingTest(), PEPStop
@*/
PetscErrorCode PEPGetStoppingTest(PEP pep,PEPStop *stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(stop,2);
  *stop = pep->stop;
  PetscFunctionReturn(0);
}

/*@
   PEPSetScale - Specifies the scaling strategy to be used.

   Logically Collective on pep

   Input Parameters:
+  pep    - the eigensolver context
.  scale  - scaling strategy
.  alpha  - the scaling factor used in the scalar strategy
.  Dl     - the left diagonal matrix of the diagonal scaling algorithm
.  Dr     - the right diagonal matrix of the diagonal scaling algorithm
.  its    - number of iterations of the diagonal scaling algorithm
-  lambda - approximation to wanted eigenvalues (modulus)

   Options Database Keys:
+  -pep_scale <type> - scaling type, one of <none,scalar,diagonal,both>
.  -pep_scale_factor <alpha> - the scaling factor
.  -pep_scale_its <its> - number of iterations
-  -pep_scale_lambda <lambda> - approximation to eigenvalues

   Notes:
   There are two non-exclusive scaling strategies: scalar and diagonal.

   In the scalar strategy, scaling is applied to the eigenvalue, that is,
   mu = lambda/alpha is the new eigenvalue and all matrices are scaled
   accordingly. After solving the scaled problem, the original lambda is
   recovered. Parameter 'alpha' must be positive. Use PETSC_DEFAULT to let
   the solver compute a reasonable scaling factor.

   In the diagonal strategy, the solver works implicitly with matrix Dl*A*Dr,
   where Dl and Dr are appropriate diagonal matrices. This improves the accuracy
   of the computed results in some cases. The user may provide the Dr and Dl
   matrices represented as Vec objects storing diagonal elements. If not
   provided, these matrices are computed internally. This option requires
   that the polynomial coefficient matrices are of MATAIJ type.
   The parameter 'its' is the number of iterations performed by the method.
   Parameter 'lambda' must be positive. Use PETSC_DEFAULT or set lambda = 1.0 if
   no information about eigenvalues is available.

   Level: intermediate

.seealso: PEPGetScale()
@*/
PetscErrorCode PEPSetScale(PEP pep,PEPScale scale,PetscReal alpha,Vec Dl,Vec Dr,PetscInt its,PetscReal lambda)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,scale,2);
  pep->scale = scale;
  if (scale==PEP_SCALE_SCALAR || scale==PEP_SCALE_BOTH) {
    PetscValidLogicalCollectiveReal(pep,alpha,3);
    if (alpha == PETSC_DEFAULT || alpha == PETSC_DECIDE) {
      pep->sfactor = 0.0;
      pep->sfactor_set = PETSC_FALSE;
    } else {
      if (alpha<=0.0) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of alpha. Must be > 0");
      pep->sfactor = alpha;
      pep->sfactor_set = PETSC_TRUE;
    }
  }
  if (scale==PEP_SCALE_DIAGONAL || scale==PEP_SCALE_BOTH) {
    if (Dl) {
      PetscValidHeaderSpecific(Dl,VEC_CLASSID,4);
      PetscCheckSameComm(pep,1,Dl,4);
      ierr = PetscObjectReference((PetscObject)Dl);CHKERRQ(ierr);
      ierr = VecDestroy(&pep->Dl);CHKERRQ(ierr);
      pep->Dl = Dl;
    }
    if (Dr) {
      PetscValidHeaderSpecific(Dr,VEC_CLASSID,5);
      PetscCheckSameComm(pep,1,Dr,5);
      ierr = PetscObjectReference((PetscObject)Dr);CHKERRQ(ierr);
      ierr = VecDestroy(&pep->Dr);CHKERRQ(ierr);
      pep->Dr = Dr;
    }
    PetscValidLogicalCollectiveInt(pep,its,6);
    PetscValidLogicalCollectiveReal(pep,lambda,7);
    if (its==PETSC_DECIDE || its==PETSC_DEFAULT) pep->sits = 5;
    else pep->sits = its;
    if (lambda==PETSC_DECIDE || lambda==PETSC_DEFAULT) pep->slambda = 1.0;
    else if (lambda<=0.0) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of lambda. Must be > 0");
    else pep->slambda = lambda;
  }
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@C
   PEPGetScale - Gets the scaling strategy used by the PEP object, and the
   associated parameters.

   Not Collectiv, but vectors are shared by all processors that share the PEP

   Input Parameter:
.  pep - the eigensolver context

   Output Parameters:
+  scale  - scaling strategy
.  alpha  - the scaling factor used in the scalar strategy
.  Dl     - the left diagonal matrix of the diagonal scaling algorithm
.  Dr     - the right diagonal matrix of the diagonal scaling algorithm
.  its    - number of iterations of the diagonal scaling algorithm
-  lambda - approximation to wanted eigenvalues (modulus)

   Level: intermediate

   Note:
   The user can specify NULL for any parameter that is not needed.

   If Dl or Dr were not set by the user, then the ones computed internally are
   returned (or a null pointer if called before PEPSetUp).

.seealso: PEPSetScale(), PEPSetUp()
@*/
PetscErrorCode PEPGetScale(PEP pep,PEPScale *scale,PetscReal *alpha,Vec *Dl,Vec *Dr,PetscInt *its,PetscReal *lambda)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (scale)  *scale  = pep->scale;
  if (alpha)  *alpha  = pep->sfactor;
  if (Dl)     *Dl     = pep->Dl;
  if (Dr)     *Dr     = pep->Dr;
  if (its)    *its    = pep->sits;
  if (lambda) *lambda = pep->slambda;
  PetscFunctionReturn(0);
}

/*@
   PEPSetExtract - Specifies the extraction strategy to be used.

   Logically Collective on pep

   Input Parameters:
+  pep     - the eigensolver context
-  extract - extraction strategy

   Options Database Keys:
.  -pep_extract <type> - extraction type, one of <none,norm,residual,structured>

   Level: intermediate

.seealso: PEPGetExtract()
@*/
PetscErrorCode PEPSetExtract(PEP pep,PEPExtract extract)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,extract,2);
  pep->extract = extract;
  PetscFunctionReturn(0);
}

/*@
   PEPGetExtract - Gets the extraction strategy used by the PEP object.

   Not Collective

   Input Parameter:
.  pep - the eigensolver context

   Output Parameter:
.  extract - extraction strategy

   Level: intermediate

.seealso: PEPSetExtract()
@*/
PetscErrorCode PEPGetExtract(PEP pep,PEPExtract *extract)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (extract) *extract = pep->extract;
  PetscFunctionReturn(0);
}

/*@
   PEPSetRefine - Specifies the refinement type (and options) to be used
   after the solve.

   Logically Collective on pep

   Input Parameters:
+  pep    - the polynomial eigensolver context
.  refine - refinement type
.  npart  - number of partitions of the communicator
.  tol    - the convergence tolerance
.  its    - maximum number of refinement iterations
-  scheme - which scheme to be used for solving the involved linear systems

   Options Database Keys:
+  -pep_refine <type> - refinement type, one of <none,simple,multiple>
.  -pep_refine_partitions <n> - the number of partitions
.  -pep_refine_tol <tol> - the tolerance
.  -pep_refine_its <its> - number of iterations
-  -pep_refine_scheme - to set the scheme for the linear solves

   Notes:
   By default, iterative refinement is disabled, since it may be very
   costly. There are two possible refinement strategies: simple and multiple.
   The simple approach performs iterative refinement on each of the
   converged eigenpairs individually, whereas the multiple strategy works
   with the invariant pair as a whole, refining all eigenpairs simultaneously.
   The latter may be required for the case of multiple eigenvalues.

   In some cases, especially when using direct solvers within the
   iterative refinement method, it may be helpful for improved scalability
   to split the communicator in several partitions. The npart parameter
   indicates how many partitions to use (defaults to 1).

   The tol and its parameters specify the stopping criterion. In the simple
   method, refinement continues until the residual of each eigenpair is
   below the tolerance (tol defaults to the PEP tol, but may be set to a
   different value). In contrast, the multiple method simply performs its
   refinement iterations (just one by default).

   The scheme argument is used to change the way in which linear systems are
   solved. Possible choices are: explicit, mixed block elimination (MBE),
   and Schur complement.

   Level: intermediate

.seealso: PEPGetRefine()
@*/
PetscErrorCode PEPSetRefine(PEP pep,PEPRefine refine,PetscInt npart,PetscReal tol,PetscInt its,PEPRefineScheme scheme)
{
  PetscErrorCode ierr;
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pep,refine,2);
  PetscValidLogicalCollectiveInt(pep,npart,3);
  PetscValidLogicalCollectiveReal(pep,tol,4);
  PetscValidLogicalCollectiveInt(pep,its,5);
  PetscValidLogicalCollectiveEnum(pep,scheme,6);
  pep->refine = refine;
  if (refine) {  /* process parameters only if not REFINE_NONE */
    if (npart!=pep->npart) {
      ierr = PetscSubcommDestroy(&pep->refinesubc);CHKERRQ(ierr);
      ierr = KSPDestroy(&pep->refineksp);CHKERRQ(ierr);
    }
    if (npart == PETSC_DEFAULT || npart == PETSC_DECIDE) {
      pep->npart = 1;
    } else {
      ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pep),&size);CHKERRQ(ierr);
      if (npart<1 || npart>size) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of npart");
      pep->npart = npart;
    }
    if (tol == PETSC_DEFAULT || tol == PETSC_DECIDE) {
      pep->rtol = PETSC_DEFAULT;
    } else {
      if (tol<=0.0) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
      pep->rtol = tol;
    }
    if (its==PETSC_DECIDE || its==PETSC_DEFAULT) {
      pep->rits = PETSC_DEFAULT;
    } else {
      if (its<0) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of its. Must be >= 0");
      pep->rits = its;
    }
    pep->scheme = scheme;
  }
  pep->state = PEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@C
   PEPGetRefine - Gets the refinement strategy used by the PEP object, and the
   associated parameters.

   Not Collective

   Input Parameter:
.  pep - the polynomial eigensolver context

   Output Parameters:
+  refine - refinement type
.  npart  - number of partitions of the communicator
.  tol    - the convergence tolerance
.  its    - maximum number of refinement iterations
-  scheme - the scheme used for solving linear systems

   Level: intermediate

   Note:
   The user can specify NULL for any parameter that is not needed.

.seealso: PEPSetRefine()
@*/
PetscErrorCode PEPGetRefine(PEP pep,PEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,PEPRefineScheme *scheme)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (refine) *refine = pep->refine;
  if (npart)  *npart  = pep->npart;
  if (tol)    *tol    = pep->rtol;
  if (its)    *its    = pep->rits;
  if (scheme) *scheme = pep->scheme;
  PetscFunctionReturn(0);
}

/*@C
   PEPSetOptionsPrefix - Sets the prefix used for searching for all
   PEP options in the database.

   Logically Collective on pep

   Input Parameters:
+  pep - the polynomial eigensolver context
-  prefix - the prefix string to prepend to all PEP option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different PEP contexts, one could call
.vb
      PEPSetOptionsPrefix(pep1,"qeig1_")
      PEPSetOptionsPrefix(pep2,"qeig2_")
.ve

   Level: advanced

.seealso: PEPAppendOptionsPrefix(), PEPGetOptionsPrefix()
@*/
PetscErrorCode PEPSetOptionsPrefix(PEP pep,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!pep->st) { ierr = PEPGetST(pep,&pep->st);CHKERRQ(ierr); }
  ierr = STSetOptionsPrefix(pep->st,prefix);CHKERRQ(ierr);
  if (!pep->V) { ierr = PEPGetBV(pep,&pep->V);CHKERRQ(ierr); }
  ierr = BVSetOptionsPrefix(pep->V,prefix);CHKERRQ(ierr);
  if (!pep->ds) { ierr = PEPGetDS(pep,&pep->ds);CHKERRQ(ierr); }
  ierr = DSSetOptionsPrefix(pep->ds,prefix);CHKERRQ(ierr);
  if (!pep->rg) { ierr = PEPGetRG(pep,&pep->rg);CHKERRQ(ierr); }
  ierr = RGSetOptionsPrefix(pep->rg,prefix);CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)pep,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PEPAppendOptionsPrefix - Appends to the prefix used for searching for all
   PEP options in the database.

   Logically Collective on pep

   Input Parameters:
+  pep - the polynomial eigensolver context
-  prefix - the prefix string to prepend to all PEP option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: PEPSetOptionsPrefix(), PEPGetOptionsPrefix()
@*/
PetscErrorCode PEPAppendOptionsPrefix(PEP pep,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!pep->st) { ierr = PEPGetST(pep,&pep->st);CHKERRQ(ierr); }
  ierr = STAppendOptionsPrefix(pep->st,prefix);CHKERRQ(ierr);
  if (!pep->V) { ierr = PEPGetBV(pep,&pep->V);CHKERRQ(ierr); }
  ierr = BVAppendOptionsPrefix(pep->V,prefix);CHKERRQ(ierr);
  if (!pep->ds) { ierr = PEPGetDS(pep,&pep->ds);CHKERRQ(ierr); }
  ierr = DSAppendOptionsPrefix(pep->ds,prefix);CHKERRQ(ierr);
  if (!pep->rg) { ierr = PEPGetRG(pep,&pep->rg);CHKERRQ(ierr); }
  ierr = RGAppendOptionsPrefix(pep->rg,prefix);CHKERRQ(ierr);
  ierr = PetscObjectAppendOptionsPrefix((PetscObject)pep,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PEPGetOptionsPrefix - Gets the prefix used for searching for all
   PEP options in the database.

   Not Collective

   Input Parameters:
.  pep - the polynomial eigensolver context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: PEPSetOptionsPrefix(), PEPAppendOptionsPrefix()
@*/
PetscErrorCode PEPGetOptionsPrefix(PEP pep,const char *prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(prefix,2);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)pep,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

