/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   MFN routines related to monitors
*/

#include <slepc/private/mfnimpl.h>   /*I "slepcmfn.h" I*/
#include <petscdraw.h>

/*
   Runs the user provided monitor routines, if any.
*/
PetscErrorCode MFNMonitor(MFN mfn,PetscInt it,PetscReal errest)
{
  PetscErrorCode ierr;
  PetscInt       i,n = mfn->numbermonitors;

  PetscFunctionBegin;
  for (i=0;i<n;i++) {
    ierr = (*mfn->monitor[i])(mfn,it,errest,mfn->monitorcontext[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   MFNMonitorSet - Sets an ADDITIONAL function to be called at every
   iteration to monitor convergence.

   Logically Collective on mfn

   Input Parameters:
+  mfn     - matrix function context obtained from MFNCreate()
.  monitor - pointer to function (if this is NULL, it turns off monitoring)
.  mctx    - [optional] context for private data for the
             monitor routine (use NULL if no context is desired)
-  monitordestroy - [optional] routine that frees monitor context (may be NULL)

   Calling Sequence of monitor:
$   monitor(MFN mfn,int its,PetscReal errest,void *mctx)

+  mfn    - matrix function context obtained from MFNCreate()
.  its    - iteration number
.  errest - error estimate
-  mctx   - optional monitoring context, as set by MFNMonitorSet()

   Options Database Keys:
+    -mfn_monitor        - print the error estimate
.    -mfn_monitor_lg     - sets line graph monitor for the error estimate
-    -mfn_monitor_cancel - cancels all monitors that have been hardwired into
      a code by calls to MFNMonitorSet(), but does not cancel those set via
      the options database.

   Notes:
   Several different monitoring routines may be set by calling
   MFNMonitorSet() multiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.seealso: MFNMonitorCancel()
@*/
PetscErrorCode MFNMonitorSet(MFN mfn,PetscErrorCode (*monitor)(MFN,PetscInt,PetscReal,void*),void *mctx,PetscErrorCode (*monitordestroy)(void**))
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  if (mfn->numbermonitors >= MAXMFNMONITORS) SETERRQ(PetscObjectComm((PetscObject)mfn),PETSC_ERR_ARG_OUTOFRANGE,"Too many MFN monitors set");
  mfn->monitor[mfn->numbermonitors]           = monitor;
  mfn->monitorcontext[mfn->numbermonitors]    = (void*)mctx;
  mfn->monitordestroy[mfn->numbermonitors++]  = monitordestroy;
  PetscFunctionReturn(0);
}

/*@
   MFNMonitorCancel - Clears all monitors for an MFN object.

   Logically Collective on mfn

   Input Parameters:
.  mfn - matrix function context obtained from MFNCreate()

   Options Database Key:
.    -mfn_monitor_cancel - Cancels all monitors that have been hardwired
      into a code by calls to MFNMonitorSet(),
      but does not cancel those set via the options database.

   Level: intermediate

.seealso: MFNMonitorSet()
@*/
PetscErrorCode MFNMonitorCancel(MFN mfn)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  for (i=0; i<mfn->numbermonitors; i++) {
    if (mfn->monitordestroy[i]) {
      ierr = (*mfn->monitordestroy[i])(&mfn->monitorcontext[i]);CHKERRQ(ierr);
    }
  }
  mfn->numbermonitors = 0;
  PetscFunctionReturn(0);
}

/*@C
   MFNGetMonitorContext - Gets the monitor context, as set by
   MFNMonitorSet() for the FIRST monitor only.

   Not Collective

   Input Parameter:
.  mfn - matrix function context obtained from MFNCreate()

   Output Parameter:
.  ctx - monitor context

   Level: intermediate

.seealso: MFNMonitorSet()
@*/
PetscErrorCode MFNGetMonitorContext(MFN mfn,void **ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  *ctx = mfn->monitorcontext[0];
  PetscFunctionReturn(0);
}

/*@C
   MFNMonitorDefault - Print the error estimate of the current approximation at each
   iteration of the matrix function solver.

   Collective on mfn

   Input Parameters:
+  mfn    - matrix function context
.  its    - iteration number
.  errest - error estimate
-  vf     - viewer and format for monitoring

   Level: intermediate

.seealso: MFNMonitorSet()
@*/
PetscErrorCode MFNMonitorDefault(MFN mfn,PetscInt its,PetscReal errest,PetscViewerAndFormat *vf)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mfn,MFN_CLASSID,1);
  PetscValidPointer(vf,4);
  viewer = vf->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,4);
  ierr = PetscViewerPushFormat(viewer,vf->format);CHKERRQ(ierr);
  ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)mfn)->tablevel);CHKERRQ(ierr);
  if (its == 1 && ((PetscObject)mfn)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Error estimates for %s solve.\n",((PetscObject)mfn)->prefix);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"%3D MFN Error estimate %14.12e\n",its,(double)errest);CHKERRQ(ierr);
  ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)mfn)->tablevel);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   MFNMonitorLGCreate - Creates a line graph context for use with
   MFN to monitor convergence.

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
.  -mfn_monitor_lg - Sets line graph monitor

   Notes:
   Use PetscDrawLGDestroy() to destroy this line graph.

   Level: intermediate

.seealso: MFNMonitorSet()
@*/
PetscErrorCode MFNMonitorLGCreate(MPI_Comm comm,const char host[],const char label[],int x,int y,int m,int n,PetscDrawLG *lgctx)
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

PetscErrorCode MFNMonitorLG(MFN mfn,PetscInt its,PetscReal errest,void *ctx)
{
  PetscDrawLG    lg = (PetscDrawLG)ctx;
  PetscReal      x,y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  if (its==1) {
    ierr = PetscDrawLGReset(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSetDimension(lg,1);CHKERRQ(ierr);
    ierr = PetscDrawLGSetLimits(lg,1,1.0,PetscLog10Real(mfn->tol)-2,0.0);CHKERRQ(ierr);
  }
  x = (PetscReal)its;
  if (errest > 0.0) y = PetscLog10Real(errest);
  else y = 0.0;
  ierr = PetscDrawLGAddPoint(lg,&x,&y);CHKERRQ(ierr);
  if (its <= 20 || !(its % 5) || mfn->reason) {
    ierr = PetscDrawLGDraw(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSave(lg);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

