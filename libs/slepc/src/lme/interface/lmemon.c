/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   LME routines related to monitors
*/

#include <slepc/private/lmeimpl.h>   /*I "slepclme.h" I*/
#include <petscdraw.h>

/*
   Runs the user provided monitor routines, if any.
*/
PetscErrorCode LMEMonitor(LME lme,PetscInt it,PetscReal errest)
{
  PetscErrorCode ierr;
  PetscInt       i,n = lme->numbermonitors;

  PetscFunctionBegin;
  for (i=0;i<n;i++) {
    ierr = (*lme->monitor[i])(lme,it,errest,lme->monitorcontext[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   LMEMonitorSet - Sets an ADDITIONAL function to be called at every
   iteration to monitor convergence.

   Logically Collective on lme

   Input Parameters:
+  lme     - linear matrix equation solver context obtained from LMECreate()
.  monitor - pointer to function (if this is NULL, it turns off monitoring)
.  mctx    - [optional] context for private data for the
             monitor routine (use NULL if no context is desired)
-  monitordestroy - [optional] routine that frees monitor context (may be NULL)

   Calling Sequence of monitor:
$   monitor(LME lme,int its,PetscReal errest,void *mctx)

+  lme    - linear matrix equation solver context obtained from LMECreate()
.  its    - iteration number
.  errest - error estimate
-  mctx   - optional monitoring context, as set by LMEMonitorSet()

   Options Database Keys:
+    -lme_monitor        - print the error estimate
.    -lme_monitor_lg     - sets line graph monitor for the error estimate
-    -lme_monitor_cancel - cancels all monitors that have been hardwired into
      a code by calls to LMEMonitorSet(), but does not cancel those set via
      the options database.

   Notes:
   Several different monitoring routines may be set by calling
   LMEMonitorSet() multiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.seealso: LMEMonitorCancel()
@*/
PetscErrorCode LMEMonitorSet(LME lme,PetscErrorCode (*monitor)(LME,PetscInt,PetscReal,void*),void *mctx,PetscErrorCode (*monitordestroy)(void**))
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  if (lme->numbermonitors >= MAXLMEMONITORS) SETERRQ(PetscObjectComm((PetscObject)lme),PETSC_ERR_ARG_OUTOFRANGE,"Too many LME monitors set");
  lme->monitor[lme->numbermonitors]           = monitor;
  lme->monitorcontext[lme->numbermonitors]    = (void*)mctx;
  lme->monitordestroy[lme->numbermonitors++]  = monitordestroy;
  PetscFunctionReturn(0);
}

/*@
   LMEMonitorCancel - Clears all monitors for an LME object.

   Logically Collective on lme

   Input Parameters:
.  lme - linear matrix equation solver context obtained from LMECreate()

   Options Database Key:
.    -lme_monitor_cancel - Cancels all monitors that have been hardwired
      into a code by calls to LMEMonitorSet(),
      but does not cancel those set via the options database.

   Level: intermediate

.seealso: LMEMonitorSet()
@*/
PetscErrorCode LMEMonitorCancel(LME lme)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  for (i=0; i<lme->numbermonitors; i++) {
    if (lme->monitordestroy[i]) {
      ierr = (*lme->monitordestroy[i])(&lme->monitorcontext[i]);CHKERRQ(ierr);
    }
  }
  lme->numbermonitors = 0;
  PetscFunctionReturn(0);
}

/*@C
   LMEGetMonitorContext - Gets the monitor context, as set by
   LMEMonitorSet() for the FIRST monitor only.

   Not Collective

   Input Parameter:
.  lme - linear matrix equation solver context obtained from LMECreate()

   Output Parameter:
.  ctx - monitor context

   Level: intermediate

.seealso: LMEMonitorSet()
@*/
PetscErrorCode LMEGetMonitorContext(LME lme,void **ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  *ctx = lme->monitorcontext[0];
  PetscFunctionReturn(0);
}

/*@C
   LMEMonitorDefault - Print the error estimate of the current approximation at each
   iteration of the linear matrix equation solver.

   Collective on lme

   Input Parameters:
+  lme    - linear matrix equation solver context
.  its    - iteration number
.  errest - error estimate
-  vf     - viewer and format for monitoring

   Level: intermediate

.seealso: LMEMonitorSet()
@*/
PetscErrorCode LMEMonitorDefault(LME lme,PetscInt its,PetscReal errest,PetscViewerAndFormat *vf)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lme,LME_CLASSID,1);
  PetscValidPointer(vf,4);
  viewer = vf->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,4);
  ierr = PetscViewerPushFormat(viewer,vf->format);CHKERRQ(ierr);
  ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)lme)->tablevel);CHKERRQ(ierr);
  if (its == 1 && ((PetscObject)lme)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Error estimates for %s solve.\n",((PetscObject)lme)->prefix);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"%3D LME Error estimate %14.12e\n",its,(double)errest);CHKERRQ(ierr);
  ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)lme)->tablevel);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   LMEMonitorLGCreate - Creates a line graph context for use with
   LME to monitor convergence.

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
.  -lme_monitor_lg - Sets line graph monitor

   Notes:
   Use PetscDrawLGDestroy() to destroy this line graph.

   Level: intermediate

.seealso: LMEMonitorSet()
@*/
PetscErrorCode LMEMonitorLGCreate(MPI_Comm comm,const char host[],const char label[],int x,int y,int m,int n,PetscDrawLG *lgctx)
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

PetscErrorCode LMEMonitorLG(LME lme,PetscInt its,PetscReal errest,void *ctx)
{
  PetscDrawLG    lg = (PetscDrawLG)ctx;
  PetscReal      x,y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  if (its==1) {
    ierr = PetscDrawLGReset(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSetDimension(lg,1);CHKERRQ(ierr);
    ierr = PetscDrawLGSetLimits(lg,1,1.0,PetscLog10Real(lme->tol)-2,0.0);CHKERRQ(ierr);
  }
  x = (PetscReal)its;
  if (errest > 0.0) y = PetscLog10Real(errest);
  else y = 0.0;
  ierr = PetscDrawLGAddPoint(lg,&x,&y);CHKERRQ(ierr);
  if (its <= 20 || !(its % 5) || lme->reason) {
    ierr = PetscDrawLGDraw(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSave(lg);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

