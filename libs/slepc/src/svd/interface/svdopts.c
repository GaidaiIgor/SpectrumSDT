/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SVD routines for setting solver options
*/

#include <slepc/private/svdimpl.h>      /*I "slepcsvd.h" I*/
#include <petscdraw.h>

/*@
   SVDSetImplicitTranspose - Indicates how to handle the transpose of the matrix
   associated with the singular value problem.

   Logically Collective on svd

   Input Parameters:
+  svd  - the singular value solver context
-  impl - how to handle the transpose (implicitly or not)

   Options Database Key:
.  -svd_implicittranspose - Activate the implicit transpose mode.

   Notes:
   By default, the transpose of the matrix is explicitly built (if the matrix
   has defined the MatTranspose operation).

   If this flag is set to true, the solver does not build the transpose, but
   handles it implicitly via MatMultTranspose() (or MatMultHermitianTranspose()
   in the complex case) operations. This is likely to be more inefficient
   than the default behaviour, both in sequential and in parallel, but
   requires less storage.

   Level: advanced

.seealso: SVDGetImplicitTranspose(), SVDSolve(), SVDSetOperator()
@*/
PetscErrorCode SVDSetImplicitTranspose(SVD svd,PetscBool impl)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,impl,2);
  if (svd->impltrans!=impl) {
    svd->impltrans = impl;
    svd->state     = SVD_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   SVDGetImplicitTranspose - Gets the mode used to handle the transpose
   of the matrix associated with the singular value problem.

   Not Collective

   Input Parameter:
.  svd  - the singular value solver context

   Output Parameter:
.  impl - how to handle the transpose (implicitly or not)

   Level: advanced

.seealso: SVDSetImplicitTranspose(), SVDSolve(), SVDSetOperator()
@*/
PetscErrorCode SVDGetImplicitTranspose(SVD svd,PetscBool *impl)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidBoolPointer(impl,2);
  *impl = svd->impltrans;
  PetscFunctionReturn(0);
}

/*@
   SVDSetTolerances - Sets the tolerance and maximum
   iteration count used by the default SVD convergence testers.

   Logically Collective on svd

   Input Parameters:
+  svd - the singular value solver context
.  tol - the convergence tolerance
-  maxits - maximum number of iterations to use

   Options Database Keys:
+  -svd_tol <tol> - Sets the convergence tolerance
-  -svd_max_it <maxits> - Sets the maximum number of iterations allowed

   Note:
   Use PETSC_DEFAULT for either argument to assign a reasonably good value.

   Level: intermediate

.seealso: SVDGetTolerances()
@*/
PetscErrorCode SVDSetTolerances(SVD svd,PetscReal tol,PetscInt maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveReal(svd,tol,2);
  PetscValidLogicalCollectiveInt(svd,maxits,3);
  if (tol == PETSC_DEFAULT) {
    svd->tol   = PETSC_DEFAULT;
    svd->state = SVD_STATE_INITIAL;
  } else {
    if (tol <= 0.0) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    svd->tol = tol;
  }
  if (maxits == PETSC_DEFAULT || maxits == PETSC_DECIDE) {
    svd->max_it = PETSC_DEFAULT;
    svd->state  = SVD_STATE_INITIAL;
  } else {
    if (maxits <= 0) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of maxits. Must be > 0");
    svd->max_it = maxits;
  }
  PetscFunctionReturn(0);
}

/*@C
   SVDGetTolerances - Gets the tolerance and maximum
   iteration count used by the default SVD convergence tests.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameters:
+  tol - the convergence tolerance
-  maxits - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: SVDSetTolerances()
@*/
PetscErrorCode SVDGetTolerances(SVD svd,PetscReal *tol,PetscInt *maxits)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (tol)    *tol    = svd->tol;
  if (maxits) *maxits = svd->max_it;
  PetscFunctionReturn(0);
}

/*@
   SVDSetDimensions - Sets the number of singular values to compute
   and the dimension of the subspace.

   Logically Collective on svd

   Input Parameters:
+  svd - the singular value solver context
.  nsv - number of singular values to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Options Database Keys:
+  -svd_nsv <nsv> - Sets the number of singular values
.  -svd_ncv <ncv> - Sets the dimension of the subspace
-  -svd_mpd <mpd> - Sets the maximum projected dimension

   Notes:
   Use PETSC_DEFAULT for ncv and mpd to assign a reasonably good value, which is
   dependent on the solution method and the number of singular values required.

   The parameters ncv and mpd are intimately related, so that the user is advised
   to set one of them at most. Normal usage is that
   (a) in cases where nsv is small, the user sets ncv (a reasonable default is 2*nsv); and
   (b) in cases where nsv is large, the user sets mpd.

   The value of ncv should always be between nsv and (nsv+mpd), typically
   ncv=nsv+mpd. If nsv is not too large, mpd=nsv is a reasonable choice, otherwise
   a smaller value should be used.

   Level: intermediate

.seealso: SVDGetDimensions()
@*/
PetscErrorCode SVDSetDimensions(SVD svd,PetscInt nsv,PetscInt ncv,PetscInt mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveInt(svd,nsv,2);
  PetscValidLogicalCollectiveInt(svd,ncv,3);
  PetscValidLogicalCollectiveInt(svd,mpd,4);
  if (nsv<1) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of nsv. Must be > 0");
  svd->nsv = nsv;
  if (ncv == PETSC_DEFAULT || ncv == PETSC_DECIDE) {
    svd->ncv = PETSC_DEFAULT;
  } else {
    if (ncv<1) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of ncv. Must be > 0");
    svd->ncv = ncv;
  }
  if (mpd == PETSC_DECIDE || mpd == PETSC_DEFAULT) {
    svd->mpd = PETSC_DEFAULT;
  } else {
    if (mpd<1) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of mpd. Must be > 0");
    svd->mpd = mpd;
  }
  svd->state = SVD_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@C
   SVDGetDimensions - Gets the number of singular values to compute
   and the dimension of the subspace.

   Not Collective

   Input Parameter:
.  svd - the singular value context

   Output Parameters:
+  nsv - number of singular values to compute
.  ncv - the maximum dimension of the subspace to be used by the solver
-  mpd - the maximum dimension allowed for the projected problem

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: SVDSetDimensions()
@*/
PetscErrorCode SVDGetDimensions(SVD svd,PetscInt *nsv,PetscInt *ncv,PetscInt *mpd)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (nsv) *nsv = svd->nsv;
  if (ncv) *ncv = svd->ncv;
  if (mpd) *mpd = svd->mpd;
  PetscFunctionReturn(0);
}

/*@
    SVDSetWhichSingularTriplets - Specifies which singular triplets are
    to be sought.

    Logically Collective on svd

    Input Parameter:
.   svd - singular value solver context obtained from SVDCreate()

    Output Parameter:
.   which - which singular triplets are to be sought

    Possible values:
    The parameter 'which' can have one of these values

+     SVD_LARGEST  - largest singular values
-     SVD_SMALLEST - smallest singular values

    Options Database Keys:
+   -svd_largest  - Sets largest singular values
-   -svd_smallest - Sets smallest singular values

    Level: intermediate

.seealso: SVDGetWhichSingularTriplets(), SVDWhich
@*/
PetscErrorCode SVDSetWhichSingularTriplets(SVD svd,SVDWhich which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,which,2);
  switch (which) {
    case SVD_LARGEST:
    case SVD_SMALLEST:
      if (svd->which != which) {
        svd->state = SVD_STATE_INITIAL;
        svd->which = which;
      }
      break;
  default:
    SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'which' parameter");
  }
  PetscFunctionReturn(0);
}

/*@
    SVDGetWhichSingularTriplets - Returns which singular triplets are
    to be sought.

    Not Collective

    Input Parameter:
.   svd - singular value solver context obtained from SVDCreate()

    Output Parameter:
.   which - which singular triplets are to be sought

    Notes:
    See SVDSetWhichSingularTriplets() for possible values of which

    Level: intermediate

.seealso: SVDSetWhichSingularTriplets(), SVDWhich
@*/
PetscErrorCode SVDGetWhichSingularTriplets(SVD svd,SVDWhich *which)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(which,2);
  *which = svd->which;
  PetscFunctionReturn(0);
}

/*@C
   SVDSetConvergenceTestFunction - Sets a function to compute the error estimate
   used in the convergence test.

   Logically Collective on svd

   Input Parameters:
+  svd     - singular value solver context obtained from SVDCreate()
.  func    - a pointer to the convergence test function
.  ctx     - context for private data for the convergence routine (may be null)
-  destroy - a routine for destroying the context (may be null)

   Calling Sequence of func:
$   func(SVD svd,PetscReal sigma,PetscReal res,PetscReal *errest,void *ctx)

+   svd    - singular value solver context obtained from SVDCreate()
.   sigma  - computed singular value
.   res    - residual norm associated to the singular triplet
.   errest - (output) computed error estimate
-   ctx    - optional context, as set by SVDSetConvergenceTestFunction()

   Note:
   If the error estimate returned by the convergence test function is less than
   the tolerance, then the singular value is accepted as converged.

   Level: advanced

.seealso: SVDSetConvergenceTest(), SVDSetTolerances()
@*/
PetscErrorCode SVDSetConvergenceTestFunction(SVD svd,PetscErrorCode (*func)(SVD,PetscReal,PetscReal,PetscReal*,void*),void* ctx,PetscErrorCode (*destroy)(void*))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (svd->convergeddestroy) {
    ierr = (*svd->convergeddestroy)(svd->convergedctx);CHKERRQ(ierr);
  }
  svd->convergeduser    = func;
  svd->convergeddestroy = destroy;
  svd->convergedctx     = ctx;
  if (func == SVDConvergedRelative) svd->conv = SVD_CONV_REL;
  else if (func == SVDConvergedAbsolute) svd->conv = SVD_CONV_ABS;
  else {
    svd->conv      = SVD_CONV_USER;
    svd->converged = svd->convergeduser;
  }
  PetscFunctionReturn(0);
}

/*@
   SVDSetConvergenceTest - Specifies how to compute the error estimate
   used in the convergence test.

   Logically Collective on svd

   Input Parameters:
+  svd  - singular value solver context obtained from SVDCreate()
-  conv - the type of convergence test

   Options Database Keys:
+  -svd_conv_abs  - Sets the absolute convergence test
.  -svd_conv_rel  - Sets the convergence test relative to the singular value
-  -svd_conv_user - Selects the user-defined convergence test

   Note:
   The parameter 'conv' can have one of these values
+     SVD_CONV_ABS  - absolute error ||r||
.     SVD_CONV_REL  - error relative to the singular value l, ||r||/sigma
-     SVD_CONV_USER - function set by SVDSetConvergenceTestFunction()

   Level: intermediate

.seealso: SVDGetConvergenceTest(), SVDSetConvergenceTestFunction(), SVDSetStoppingTest(), SVDConv
@*/
PetscErrorCode SVDSetConvergenceTest(SVD svd,SVDConv conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,conv,2);
  switch (conv) {
    case SVD_CONV_ABS:  svd->converged = SVDConvergedAbsolute; break;
    case SVD_CONV_REL:  svd->converged = SVDConvergedRelative; break;
    case SVD_CONV_USER:
      if (!svd->convergeduser) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ORDER,"Must call SVDSetConvergenceTestFunction() first");
      svd->converged = svd->convergeduser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'conv' value");
  }
  svd->conv = conv;
  PetscFunctionReturn(0);
}

/*@
   SVDGetConvergenceTest - Gets the method used to compute the error estimate
   used in the convergence test.

   Not Collective

   Input Parameters:
.  svd   - singular value solver context obtained from SVDCreate()

   Output Parameters:
.  conv  - the type of convergence test

   Level: intermediate

.seealso: SVDSetConvergenceTest(), SVDConv
@*/
PetscErrorCode SVDGetConvergenceTest(SVD svd,SVDConv *conv)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(conv,2);
  *conv = svd->conv;
  PetscFunctionReturn(0);
}

/*@C
   SVDSetStoppingTestFunction - Sets a function to decide when to stop the outer
   iteration of the singular value solver.

   Logically Collective on svd

   Input Parameters:
+  svd     - singular value solver context obtained from SVDCreate()
.  func    - pointer to the stopping test function
.  ctx     - context for private data for the stopping routine (may be null)
-  destroy - a routine for destroying the context (may be null)

   Calling Sequence of func:
$   func(SVD svd,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nsv,SVDConvergedReason *reason,void *ctx)

+   svd    - singular value solver context obtained from SVDCreate()
.   its    - current number of iterations
.   max_it - maximum number of iterations
.   nconv  - number of currently converged singular triplets
.   nsv    - number of requested singular triplets
.   reason - (output) result of the stopping test
-   ctx    - optional context, as set by SVDSetStoppingTestFunction()

   Note:
   Normal usage is to first call the default routine SVDStoppingBasic() and then
   set reason to SVD_CONVERGED_USER if some user-defined conditions have been
   met. To let the singular value solver continue iterating, the result must be
   left as SVD_CONVERGED_ITERATING.

   Level: advanced

.seealso: SVDSetStoppingTest(), SVDStoppingBasic()
@*/
PetscErrorCode SVDSetStoppingTestFunction(SVD svd,PetscErrorCode (*func)(SVD,PetscInt,PetscInt,PetscInt,PetscInt,SVDConvergedReason*,void*),void* ctx,PetscErrorCode (*destroy)(void*))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (svd->stoppingdestroy) {
    ierr = (*svd->stoppingdestroy)(svd->stoppingctx);CHKERRQ(ierr);
  }
  svd->stoppinguser    = func;
  svd->stoppingdestroy = destroy;
  svd->stoppingctx     = ctx;
  if (func == SVDStoppingBasic) svd->stop = SVD_STOP_BASIC;
  else {
    svd->stop     = SVD_STOP_USER;
    svd->stopping = svd->stoppinguser;
  }
  PetscFunctionReturn(0);
}

/*@
   SVDSetStoppingTest - Specifies how to decide the termination of the outer
   loop of the singular value solver.

   Logically Collective on svd

   Input Parameters:
+  svd  - singular value solver context obtained from SVDCreate()
-  stop - the type of stopping test

   Options Database Keys:
+  -svd_stop_basic - Sets the default stopping test
-  -svd_stop_user  - Selects the user-defined stopping test

   Note:
   The parameter 'stop' can have one of these values
+     SVD_STOP_BASIC - default stopping test
-     SVD_STOP_USER  - function set by SVDSetStoppingTestFunction()

   Level: advanced

.seealso: SVDGetStoppingTest(), SVDSetStoppingTestFunction(), SVDSetConvergenceTest(), SVDStop
@*/
PetscErrorCode SVDSetStoppingTest(SVD svd,SVDStop stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveEnum(svd,stop,2);
  switch (stop) {
    case SVD_STOP_BASIC: svd->stopping = SVDStoppingBasic; break;
    case SVD_STOP_USER:
      if (!svd->stoppinguser) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ORDER,"Must call SVDSetStoppingTestFunction() first");
      svd->stopping = svd->stoppinguser;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Invalid 'stop' value");
  }
  svd->stop = stop;
  PetscFunctionReturn(0);
}

/*@
   SVDGetStoppingTest - Gets the method used to decide the termination of the outer
   loop of the singular value solver.

   Not Collective

   Input Parameters:
.  svd   - singular value solver context obtained from SVDCreate()

   Output Parameters:
.  stop  - the type of stopping test

   Level: advanced

.seealso: SVDSetStoppingTest(), SVDStop
@*/
PetscErrorCode SVDGetStoppingTest(SVD svd,SVDStop *stop)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(stop,2);
  *stop = svd->stop;
  PetscFunctionReturn(0);
}

/*@C
   SVDMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user.

   Collective on svd

   Input Parameters:
+  svd      - the singular value solver context
.  name     - the monitor option name
.  help     - message indicating what monitoring is done
.  manual   - manual page for the monitor
.  monitor  - the monitor function, whose context is a PetscViewerAndFormat
-  trackall - whether this monitor tracks all singular values or not

   Level: developer

.seealso: SVDMonitorSet(), SVDSetTrackAll(), SVDConvMonitorSetFromOptions()
@*/
PetscErrorCode SVDMonitorSetFromOptions(SVD svd,const char name[],const char help[],const char manual[],PetscErrorCode (*monitor)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,PetscViewerAndFormat*),PetscBool trackall)
{
  PetscErrorCode       ierr;
  PetscBool            flg;
  PetscViewer          viewer;
  PetscViewerFormat    format;
  PetscViewerAndFormat *vf;

  PetscFunctionBegin;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)svd),((PetscObject)svd)->options,((PetscObject)svd)->prefix,name,&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerAndFormatCreate(viewer,format,&vf);CHKERRQ(ierr);
    ierr = PetscObjectDereference((PetscObject)viewer);CHKERRQ(ierr);
    ierr = SVDMonitorSet(svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))monitor,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);
    if (trackall) {
      ierr = SVDSetTrackAll(svd,PETSC_TRUE);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@C
   SVDConvMonitorSetFromOptions - Sets a monitor function and viewer appropriate for the type
   indicated by the user (for monitors that only show iteration numbers of convergence).

   Collective on svd

   Input Parameters:
+  svd      - the singular value solver context
.  name     - the monitor option name
.  help     - message indicating what monitoring is done
.  manual   - manual page for the monitor
-  monitor  - the monitor function, whose context is a SlepcConvMonitor

   Level: developer

.seealso: SVDMonitorSet(), SVDMonitorSetFromOptions()
@*/
PetscErrorCode SVDConvMonitorSetFromOptions(SVD svd,const char name[],const char help[],const char manual[],PetscErrorCode (*monitor)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,SlepcConvMonitor))
{
  PetscErrorCode    ierr;
  PetscBool         flg;
  PetscViewer       viewer;
  PetscViewerFormat format;
  SlepcConvMonitor  ctx;

  PetscFunctionBegin;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)svd),((PetscObject)svd)->options,((PetscObject)svd)->prefix,name,&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = SlepcConvMonitorCreate(viewer,format,&ctx);CHKERRQ(ierr);
    ierr = PetscObjectDereference((PetscObject)viewer);CHKERRQ(ierr);
    ierr = SVDMonitorSet(svd,(PetscErrorCode (*)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*))monitor,ctx,(PetscErrorCode (*)(void**))SlepcConvMonitorDestroy);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   SVDSetFromOptions - Sets SVD options from the options database.
   This routine must be called before SVDSetUp() if the user is to be
   allowed to set the solver type.

   Collective on svd

   Input Parameters:
.  svd - the singular value solver context

   Notes:
   To see all options, run your program with the -help option.

   Level: beginner

.seealso:
@*/
PetscErrorCode SVDSetFromOptions(SVD svd)
{
  PetscErrorCode ierr;
  char           type[256];
  PetscBool      set,flg,val,flg1,flg2,flg3;
  PetscInt       i,j,k;
  PetscReal      r;
  PetscDrawLG    lg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  ierr = SVDRegisterAll();CHKERRQ(ierr);
  ierr = PetscObjectOptionsBegin((PetscObject)svd);CHKERRQ(ierr);
    ierr = PetscOptionsFList("-svd_type","SVD solver method","SVDSetType",SVDList,(char*)(((PetscObject)svd)->type_name?((PetscObject)svd)->type_name:SVDCROSS),type,256,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = SVDSetType(svd,type);CHKERRQ(ierr);
    } else if (!((PetscObject)svd)->type_name) {
      ierr = SVDSetType(svd,SVDCROSS);CHKERRQ(ierr);
    }

    ierr = PetscOptionsBool("-svd_implicittranspose","Handle matrix transpose implicitly","SVDSetImplicitTranspose",svd->impltrans,&val,&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDSetImplicitTranspose(svd,val);CHKERRQ(ierr); }

    i = svd->max_it;
    ierr = PetscOptionsInt("-svd_max_it","Maximum number of iterations","SVDSetTolerances",svd->max_it,&i,&flg1);CHKERRQ(ierr);
    r = svd->tol;
    ierr = PetscOptionsReal("-svd_tol","Tolerance","SVDSetTolerances",svd->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL:svd->tol,&r,&flg2);CHKERRQ(ierr);
    if (flg1 || flg2) { ierr = SVDSetTolerances(svd,r,i);CHKERRQ(ierr); }

    ierr = PetscOptionsBoolGroupBegin("-svd_conv_rel","Relative error convergence test","SVDSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDSetConvergenceTest(svd,SVD_CONV_REL);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroup("-svd_conv_abs","Absolute error convergence test","SVDSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDSetConvergenceTest(svd,SVD_CONV_ABS);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-svd_conv_user","User-defined convergence test","SVDSetConvergenceTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDSetConvergenceTest(svd,SVD_CONV_USER);CHKERRQ(ierr); }

    ierr = PetscOptionsBoolGroupBegin("-svd_stop_basic","Stop iteration if all singular values converged or max_it reached","SVDSetStoppingTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDSetStoppingTest(svd,SVD_STOP_BASIC);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-svd_stop_user","User-defined stopping test","SVDSetStoppingTest",&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDSetStoppingTest(svd,SVD_STOP_USER);CHKERRQ(ierr); }

    i = svd->nsv;
    ierr = PetscOptionsInt("-svd_nsv","Number of singular values to compute","SVDSetDimensions",svd->nsv,&i,&flg1);CHKERRQ(ierr);
    j = svd->ncv;
    ierr = PetscOptionsInt("-svd_ncv","Number of basis vectors","SVDSetDimensions",svd->ncv,&j,&flg2);CHKERRQ(ierr);
    k = svd->mpd;
    ierr = PetscOptionsInt("-svd_mpd","Maximum dimension of projected problem","SVDSetDimensions",svd->mpd,&k,&flg3);CHKERRQ(ierr);
    if (flg1 || flg2 || flg3) { ierr = SVDSetDimensions(svd,i,j,k);CHKERRQ(ierr); }

    ierr = PetscOptionsBoolGroupBegin("-svd_largest","Compute largest singular values","SVDSetWhichSingularTriplets",&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDSetWhichSingularTriplets(svd,SVD_LARGEST);CHKERRQ(ierr); }
    ierr = PetscOptionsBoolGroupEnd("-svd_smallest","Compute smallest singular values","SVDSetWhichSingularTriplets",&flg);CHKERRQ(ierr);
    if (flg) { ierr = SVDSetWhichSingularTriplets(svd,SVD_SMALLEST);CHKERRQ(ierr); }

    /* -----------------------------------------------------------------------*/
    /*
      Cancels all monitors hardwired into code before call to SVDSetFromOptions()
    */
    ierr = PetscOptionsBool("-svd_monitor_cancel","Remove any hardwired monitor routines","SVDMonitorCancel",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = SVDMonitorCancel(svd);CHKERRQ(ierr);
    }
    /*
      Text monitors
    */
    ierr = SVDMonitorSetFromOptions(svd,"-svd_monitor","Monitor first unconverged approximate singular value and error estimate","SVDMonitorFirst",SVDMonitorFirst,PETSC_FALSE);CHKERRQ(ierr);
    ierr = SVDConvMonitorSetFromOptions(svd,"-svd_monitor_conv","Monitor approximate singular values and error estimates as they converge","SVDMonitorConverged",SVDMonitorConverged);CHKERRQ(ierr);
    ierr = SVDMonitorSetFromOptions(svd,"-svd_monitor_all","Monitor approximate singular values and error estimates","SVDMonitorAll",SVDMonitorAll,PETSC_TRUE);CHKERRQ(ierr);
    /*
      Line graph monitors
    */
    ierr = PetscOptionsBool("-svd_monitor_lg","Monitor first unconverged approximate singular value and error estimate graphically","SVDMonitorSet",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = SVDMonitorLGCreate(PetscObjectComm((PetscObject)svd),NULL,"Error estimates",PETSC_DECIDE,PETSC_DECIDE,300,300,&lg);CHKERRQ(ierr);
      ierr = SVDMonitorSet(svd,SVDMonitorLG,lg,(PetscErrorCode (*)(void**))PetscDrawLGDestroy);CHKERRQ(ierr);
    }
    ierr = PetscOptionsBool("-svd_monitor_lg_all","Monitor error estimates graphically","SVDMonitorSet",PETSC_FALSE,&flg,&set);CHKERRQ(ierr);
    if (set && flg) {
      ierr = SVDMonitorLGCreate(PetscObjectComm((PetscObject)svd),NULL,"Error estimates",PETSC_DECIDE,PETSC_DECIDE,300,300,&lg);CHKERRQ(ierr);
      ierr = SVDMonitorSet(svd,SVDMonitorLGAll,lg,(PetscErrorCode (*)(void**))PetscDrawLGDestroy);CHKERRQ(ierr);
      ierr = SVDSetTrackAll(svd,PETSC_TRUE);CHKERRQ(ierr);
    }

    /* -----------------------------------------------------------------------*/
    ierr = PetscOptionsName("-svd_view","Print detailed information on solver used","SVDView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-svd_view_vectors","View computed singular vectors","SVDVectorsView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-svd_view_values","View computed singular values","SVDValuesView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-svd_converged_reason","Print reason for convergence, and number of iterations","SVDReasonView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-svd_error_absolute","Print absolute errors of each singular triplet","SVDErrorView",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsName("-svd_error_relative","Print relative errors of each singular triplet","SVDErrorView",NULL);CHKERRQ(ierr);

    if (svd->ops->setfromoptions) {
      ierr = (*svd->ops->setfromoptions)(PetscOptionsObject,svd);CHKERRQ(ierr);
    }
    ierr = PetscObjectProcessOptionsHandlers(PetscOptionsObject,(PetscObject)svd);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  if (!svd->V) { ierr = SVDGetBV(svd,&svd->V,NULL);CHKERRQ(ierr); }
  ierr = BVSetFromOptions(svd->V);CHKERRQ(ierr);
  if (!svd->U) { ierr = SVDGetBV(svd,NULL,&svd->U);CHKERRQ(ierr); }
  ierr = BVSetFromOptions(svd->U);CHKERRQ(ierr);
  if (!svd->ds) { ierr = SVDGetDS(svd,&svd->ds);CHKERRQ(ierr); }
  ierr = DSSetFromOptions(svd->ds);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   SVDSetTrackAll - Specifies if the solver must compute the residual norm of all
   approximate singular value or not.

   Logically Collective on svd

   Input Parameters:
+  svd      - the singular value solver context
-  trackall - whether to compute all residuals or not

   Notes:
   If the user sets trackall=PETSC_TRUE then the solver computes (or estimates)
   the residual norm for each singular value approximation. Computing the residual is
   usually an expensive operation and solvers commonly compute only the residual
   associated to the first unconverged singular value.

   The options '-svd_monitor_all' and '-svd_monitor_lg_all' automatically
   activate this option.

   Level: developer

.seealso: SVDGetTrackAll()
@*/
PetscErrorCode SVDSetTrackAll(SVD svd,PetscBool trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,trackall,2);
  svd->trackall = trackall;
  PetscFunctionReturn(0);
}

/*@
   SVDGetTrackAll - Returns the flag indicating whether all residual norms must
   be computed or not.

   Not Collective

   Input Parameter:
.  svd - the singular value solver context

   Output Parameter:
.  trackall - the returned flag

   Level: developer

.seealso: SVDSetTrackAll()
@*/
PetscErrorCode SVDGetTrackAll(SVD svd,PetscBool *trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidBoolPointer(trackall,2);
  *trackall = svd->trackall;
  PetscFunctionReturn(0);
}


/*@C
   SVDSetOptionsPrefix - Sets the prefix used for searching for all
   SVD options in the database.

   Logically Collective on svd

   Input Parameters:
+  svd - the singular value solver context
-  prefix - the prefix string to prepend to all SVD option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   For example, to distinguish between the runtime options for two
   different SVD contexts, one could call
.vb
      SVDSetOptionsPrefix(svd1,"svd1_")
      SVDSetOptionsPrefix(svd2,"svd2_")
.ve

   Level: advanced

.seealso: SVDAppendOptionsPrefix(), SVDGetOptionsPrefix()
@*/
PetscErrorCode SVDSetOptionsPrefix(SVD svd,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!svd->V) { ierr = SVDGetBV(svd,&svd->V,&svd->U);CHKERRQ(ierr); }
  ierr = BVSetOptionsPrefix(svd->V,prefix);CHKERRQ(ierr);
  ierr = BVSetOptionsPrefix(svd->U,prefix);CHKERRQ(ierr);
  if (!svd->ds) { ierr = SVDGetDS(svd,&svd->ds);CHKERRQ(ierr); }
  ierr = DSSetOptionsPrefix(svd->ds,prefix);CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)svd,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   SVDAppendOptionsPrefix - Appends to the prefix used for searching for all
   SVD options in the database.

   Logically Collective on svd

   Input Parameters:
+  svd - the singular value solver context
-  prefix - the prefix string to prepend to all SVD option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the hyphen.

   Level: advanced

.seealso: SVDSetOptionsPrefix(), SVDGetOptionsPrefix()
@*/
PetscErrorCode SVDAppendOptionsPrefix(SVD svd,const char *prefix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (!svd->V) { ierr = SVDGetBV(svd,&svd->V,&svd->U);CHKERRQ(ierr); }
  ierr = BVAppendOptionsPrefix(svd->V,prefix);CHKERRQ(ierr);
  ierr = BVAppendOptionsPrefix(svd->U,prefix);CHKERRQ(ierr);
  if (!svd->ds) { ierr = SVDGetDS(svd,&svd->ds);CHKERRQ(ierr); }
  ierr = DSAppendOptionsPrefix(svd->ds,prefix);CHKERRQ(ierr);
  ierr = PetscObjectAppendOptionsPrefix((PetscObject)svd,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   SVDGetOptionsPrefix - Gets the prefix used for searching for all
   SVD options in the database.

   Not Collective

   Input Parameters:
.  svd - the singular value solver context

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

   Note:
   On the Fortran side, the user should pass in a string 'prefix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: SVDSetOptionsPrefix(), SVDAppendOptionsPrefix()
@*/
PetscErrorCode SVDGetOptionsPrefix(SVD svd,const char *prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(prefix,2);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)svd,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

