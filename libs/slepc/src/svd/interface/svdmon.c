/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SVD routines related to monitors
*/

#include <slepc/private/svdimpl.h>   /*I "slepcsvd.h" I*/
#include <petscdraw.h>

/*
   Runs the user provided monitor routines, if any.
*/
PetscErrorCode SVDMonitor(SVD svd,PetscInt it,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest)
{
  PetscErrorCode ierr;
  PetscInt       i,n = svd->numbermonitors;

  PetscFunctionBegin;
  for (i=0;i<n;i++) {
    ierr = (*svd->monitor[i])(svd,it,nconv,sigma,errest,nest,svd->monitorcontext[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   SVDMonitorSet - Sets an ADDITIONAL function to be called at every
   iteration to monitor the error estimates for each requested singular triplet.

   Collective on svd

   Input Parameters:
+  svd     - singular value solver context obtained from SVDCreate()
.  monitor - pointer to function (if this is NULL, it turns off monitoring)
-  mctx    - [optional] context for private data for the
             monitor routine (use NULL if no context is desired)

   Calling Sequence of monitor:
$   monitor(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,void *mctx)

+  svd    - singular value solver context obtained from SVDCreate()
.  its    - iteration number
.  nconv  - number of converged singular triplets
.  sigma  - singular values
.  errest - relative error estimates for each singular triplet
.  nest   - number of error estimates
-  mctx   - optional monitoring context, as set by SVDMonitorSet()

   Options Database Keys:
+    -svd_monitor        - print only the first error estimate
.    -svd_monitor_all    - print error estimates at each iteration
.    -svd_monitor_conv   - print the singular value approximations only when
      convergence has been reached
.    -svd_monitor_lg     - sets line graph monitor for the first unconverged
      approximate singular value
.    -svd_monitor_lg_all - sets line graph monitor for all unconverged
      approximate singular values
-    -svd_monitor_cancel - cancels all monitors that have been hardwired into
      a code by calls to SVDMonitorSet(), but does not cancel those set via
      the options database.

   Notes:
   Several different monitoring routines may be set by calling
   SVDMonitorSet() multiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.seealso: SVDMonitorFirst(), SVDMonitorAll(), SVDMonitorCancel()
@*/
PetscErrorCode SVDMonitorSet(SVD svd,PetscErrorCode (*monitor)(SVD,PetscInt,PetscInt,PetscReal*,PetscReal*,PetscInt,void*),void *mctx,PetscErrorCode (*monitordestroy)(void**))
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  if (svd->numbermonitors >= MAXSVDMONITORS) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_ARG_OUTOFRANGE,"Too many SVD monitors set");
  svd->monitor[svd->numbermonitors]           = monitor;
  svd->monitorcontext[svd->numbermonitors]    = (void*)mctx;
  svd->monitordestroy[svd->numbermonitors++]  = monitordestroy;
  PetscFunctionReturn(0);
}

/*@
   SVDMonitorCancel - Clears all monitors for an SVD object.

   Collective on svd

   Input Parameters:
.  svd - singular value solver context obtained from SVDCreate()

   Options Database Key:
.    -svd_monitor_cancel - Cancels all monitors that have been hardwired
      into a code by calls to SVDMonitorSet(),
      but does not cancel those set via the options database.

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDMonitorCancel(SVD svd)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  for (i=0; i<svd->numbermonitors; i++) {
    if (svd->monitordestroy[i]) {
      ierr = (*svd->monitordestroy[i])(&svd->monitorcontext[i]);CHKERRQ(ierr);
    }
  }
  svd->numbermonitors = 0;
  PetscFunctionReturn(0);
}

/*@C
   SVDGetMonitorContext - Gets the monitor context, as set by
   SVDMonitorSet() for the FIRST monitor only.

   Not Collective

   Input Parameter:
.  svd - singular value solver context obtained from SVDCreate()

   Output Parameter:
.  ctx - monitor context

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDGetMonitorContext(SVD svd,void **ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  *ctx = svd->monitorcontext[0];
  PetscFunctionReturn(0);
}

/*@C
   SVDMonitorAll - Print the current approximate values and
   error estimates at each iteration of the singular value solver.

   Collective on svd

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets so far
.  sigma  - singular values
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Level: intermediate

.seealso: SVDMonitorSet(), SVDMonitorFirst(), SVDMonitorConverged()
@*/
PetscErrorCode SVDMonitorAll(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(vf,7);
  viewer = vf->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  ierr = PetscViewerPushFormat(viewer,vf->format);CHKERRQ(ierr);
  ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)svd)->tablevel);CHKERRQ(ierr);
  if (its==1 && ((PetscObject)svd)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Singular value approximations and residual norms for %s solve.\n",((PetscObject)svd)->prefix);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"%3D SVD nconv=%D Values (Errors)",its,nconv);CHKERRQ(ierr);
  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
  for (i=0;i<nest;i++) {
    ierr = PetscViewerASCIIPrintf(viewer," %g (%10.8e)",(double)sigma[i],(double)errest[i]);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)svd)->tablevel);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   SVDMonitorFirst - Print the first unconverged approximate values and
   error estimates at each iteration of the singular value solver.

   Collective on svd

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets so far
.  sigma  - singular values
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Level: intermediate

.seealso: SVDMonitorSet(), SVDMonitorAll(), SVDMonitorConverged()
@*/
PetscErrorCode SVDMonitorFirst(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(vf,7);
  viewer = vf->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  if (its==1 && ((PetscObject)svd)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Singular value approximations and residual norms for %s solve.\n",((PetscObject)svd)->prefix);CHKERRQ(ierr);
  }
  if (nconv<nest) {
    ierr = PetscViewerPushFormat(viewer,vf->format);CHKERRQ(ierr);
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)svd)->tablevel);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%3D SVD nconv=%D first unconverged value (error)",its,nconv);CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer," %g (%10.8e)\n",(double)sigma[nconv],(double)errest[nconv]);CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)svd)->tablevel);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   SVDMonitorConverged - Print the approximate values and error estimates as they converge.

   Collective on svd

   Input Parameters:
+  svd    - singular value solver context
.  its    - iteration number
.  nconv  - number of converged singular triplets so far
.  sigma  - singular values
.  errest - error estimates
.  nest   - number of error estimates to display
-  ctx    - monitor context

   Level: intermediate

.seealso: SVDMonitorSet(), SVDMonitorFirst(), SVDMonitorAll()
@*/
PetscErrorCode SVDMonitorConverged(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,SlepcConvMonitor ctx)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(ctx,7);
  viewer = ctx->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,7);
  if (its==1 && ((PetscObject)svd)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Convergence history for %s solve.\n",((PetscObject)svd)->prefix);CHKERRQ(ierr);
  }
  if (its==1) ctx->oldnconv = 0;
  if (ctx->oldnconv!=nconv) {
    ierr = PetscViewerPushFormat(viewer,ctx->format);CHKERRQ(ierr);
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)svd)->tablevel);CHKERRQ(ierr);
    for (i=ctx->oldnconv;i<nconv;i++) {
      ierr = PetscViewerASCIIPrintf(viewer,"%3D SVD converged value (error) #%D",its,i);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer," %g (%10.8e)\n",(double)sigma[i],(double)errest[i]);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)svd)->tablevel);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ctx->oldnconv = nconv;
  }
  PetscFunctionReturn(0);
}

/*@C
   SVDMonitorLGCreate - Creates a line graph context for use with
   SVD to monitor convergence.

   Collective

   Input Parameters:
+  comm - communicator context
.  host - the X display to open, or null for the local machine
.  label - the title to put in the title bar
.  x, y - the screen coordinates of the upper left coordinate of
          the window
-  m, n - the screen width and height in pixels

   Output Parameter:
.  lgctx - the drawing context

   Options Database Keys:
+  -svd_monitor_lg - Sets line graph monitor for the first residual
-  -svd_monitor_lg_all - Sets line graph monitor for all residuals

   Notes:
   Use PetscDrawLGDestroy() to destroy this line graph.

   Level: intermediate

.seealso: SVDMonitorSet()
@*/
PetscErrorCode SVDMonitorLGCreate(MPI_Comm comm,const char host[],const char label[],int x,int y,int m,int n,PetscDrawLG *lgctx)
{
  PetscDraw      draw;
  PetscDrawLG    lg;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscDrawCreate(comm,host,label,x,y,m,n,&draw);CHKERRQ(ierr);
  ierr = PetscDrawSetFromOptions(draw);CHKERRQ(ierr);
  ierr = PetscDrawLGCreate(draw,1,&lg);CHKERRQ(ierr);
  ierr = PetscDrawLGSetFromOptions(lg);CHKERRQ(ierr);
  ierr = PetscDrawDestroy(&draw);CHKERRQ(ierr);
  *lgctx = lg;
  PetscFunctionReturn(0);
}

PetscErrorCode SVDMonitorLG(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscDrawLG    lg = (PetscDrawLG)ctx;
  PetscReal      x,y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  if (its==1) {
    ierr = PetscDrawLGReset(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSetDimension(lg,1);CHKERRQ(ierr);
    ierr = PetscDrawLGSetLimits(lg,1,1.0,PetscLog10Real(svd->tol)-2,0.0);CHKERRQ(ierr);
  }
  x = (PetscReal)its;
  if (errest[nconv] > 0.0) y = PetscLog10Real(errest[nconv]);
  else y = 0.0;
  ierr = PetscDrawLGAddPoint(lg,&x,&y);CHKERRQ(ierr);
  if (its <= 20 || !(its % 5) || svd->reason) {
    ierr = PetscDrawLGDraw(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSave(lg);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SVDMonitorLGAll(SVD svd,PetscInt its,PetscInt nconv,PetscReal *sigma,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscDrawLG    lg = (PetscDrawLG)ctx;
  PetscInt       i,n = PetscMin(svd->nsv,255);
  PetscReal      *x,*y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  if (its==1) {
    ierr = PetscDrawLGReset(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSetDimension(lg,n);CHKERRQ(ierr);
    ierr = PetscDrawLGSetLimits(lg,1,1.0,PetscLog10Real(svd->tol)-2,0.0);CHKERRQ(ierr);
  }
  ierr = PetscMalloc2(n,&x,n,&y);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
    x[i] = (PetscReal)its;
    if (i < nest && errest[i] > 0.0) y[i] = PetscLog10Real(errest[i]);
    else y[i] = 0.0;
  }
  ierr = PetscDrawLGAddPoint(lg,x,y);CHKERRQ(ierr);
  if (its <= 20 || !(its % 5) || svd->reason) {
    ierr = PetscDrawLGDraw(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSave(lg);CHKERRQ(ierr);
  }
  ierr = PetscFree2(x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

