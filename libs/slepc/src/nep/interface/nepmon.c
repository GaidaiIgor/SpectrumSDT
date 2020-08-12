/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to monitors
*/

#include <slepc/private/nepimpl.h>      /*I "slepcnep.h" I*/
#include <petscdraw.h>

/*
   Runs the user provided monitor routines, if any.
*/
PetscErrorCode NEPMonitor(NEP nep,PetscInt it,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest)
{
  PetscErrorCode ierr;
  PetscInt       i,n = nep->numbermonitors;

  PetscFunctionBegin;
  for (i=0;i<n;i++) {
    ierr = (*nep->monitor[i])(nep,it,nconv,eigr,eigi,errest,nest,nep->monitorcontext[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   NEPMonitorSet - Sets an ADDITIONAL function to be called at every
   iteration to monitor the error estimates for each requested eigenpair.

   Logically Collective on nep

   Input Parameters:
+  nep     - eigensolver context obtained from NEPCreate()
.  monitor - pointer to function (if this is NULL, it turns off monitoring)
.  mctx    - [optional] context for private data for the
             monitor routine (use NULL if no context is desired)
-  monitordestroy - [optional] routine that frees monitor context
             (may be NULL)

   Calling Sequence of monitor:
$   monitor(NEP nep,int its,int nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal* errest,int nest,void *mctx)

+  nep    - nonlinear eigensolver context obtained from NEPCreate()
.  its    - iteration number
.  nconv  - number of converged eigenpairs
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates for each eigenpair
.  nest   - number of error estimates
-  mctx   - optional monitoring context, as set by NEPMonitorSet()

   Options Database Keys:
+    -nep_monitor        - print only the first error estimate
.    -nep_monitor_all    - print error estimates at each iteration
.    -nep_monitor_conv   - print the eigenvalue approximations only when
      convergence has been reached
.    -nep_monitor_lg     - sets line graph monitor for the first unconverged
      approximate eigenvalue
.    -nep_monitor_lg_all - sets line graph monitor for all unconverged
      approximate eigenvalues
-    -nep_monitor_cancel - cancels all monitors that have been hardwired into
      a code by calls to NEPMonitorSet(), but does not cancel those set via
      the options database.

   Notes:
   Several different monitoring routines may be set by calling
   NEPMonitorSet() multiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.seealso: NEPMonitorFirst(), NEPMonitorAll(), NEPMonitorCancel()
@*/
PetscErrorCode NEPMonitorSet(NEP nep,PetscErrorCode (*monitor)(NEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*),void *mctx,PetscErrorCode (*monitordestroy)(void**))
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (nep->numbermonitors >= MAXNEPMONITORS) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Too many NEP monitors set");
  nep->monitor[nep->numbermonitors]           = monitor;
  nep->monitorcontext[nep->numbermonitors]    = (void*)mctx;
  nep->monitordestroy[nep->numbermonitors++]  = monitordestroy;
  PetscFunctionReturn(0);
}

/*@
   NEPMonitorCancel - Clears all monitors for a NEP object.

   Logically Collective on nep

   Input Parameters:
.  nep - eigensolver context obtained from NEPCreate()

   Options Database Key:
.    -nep_monitor_cancel - Cancels all monitors that have been hardwired
      into a code by calls to NEPMonitorSet(),
      but does not cancel those set via the options database.

   Level: intermediate

.seealso: NEPMonitorSet()
@*/
PetscErrorCode NEPMonitorCancel(NEP nep)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  for (i=0; i<nep->numbermonitors; i++) {
    if (nep->monitordestroy[i]) {
      ierr = (*nep->monitordestroy[i])(&nep->monitorcontext[i]);CHKERRQ(ierr);
    }
  }
  nep->numbermonitors = 0;
  PetscFunctionReturn(0);
}

/*@C
   NEPGetMonitorContext - Gets the monitor context, as set by
   NEPMonitorSet() for the FIRST monitor only.

   Not Collective

   Input Parameter:
.  nep - eigensolver context obtained from NEPCreate()

   Output Parameter:
.  ctx - monitor context

   Level: intermediate

.seealso: NEPMonitorSet(), NEPDefaultMonitor()
@*/
PetscErrorCode NEPGetMonitorContext(NEP nep,void **ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  *ctx = nep->monitorcontext[0];
  PetscFunctionReturn(0);
}

/*@C
   NEPMonitorAll - Print the current approximate values and
   error estimates at each iteration of the nonlinear eigensolver.

   Collective on nep

   Input Parameters:
+  nep    - nonlinear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Level: intermediate

.seealso: NEPMonitorSet(), NEPMonitorFirst(), NEPMonitorConverged()
@*/
PetscErrorCode NEPMonitorAll(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(vf,8);
  viewer = vf->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  ierr = PetscViewerPushFormat(viewer,vf->format);CHKERRQ(ierr);
  ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)nep)->tablevel);CHKERRQ(ierr);
  if (its==1 && ((PetscObject)nep)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)nep)->prefix);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"%3D NEP nconv=%D Values (Errors)",its,nconv);CHKERRQ(ierr);
  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
  for (i=0;i<nest;i++) {
#if defined(PETSC_USE_COMPLEX)
    ierr = PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(eigr[i]),(double)PetscImaginaryPart(eigr[i]));CHKERRQ(ierr);
#else
    ierr = PetscViewerASCIIPrintf(viewer," %g",(double)eigr[i]);CHKERRQ(ierr);
    if (eigi[i]!=0.0) { ierr = PetscViewerASCIIPrintf(viewer,"%+gi",(double)eigi[i]);CHKERRQ(ierr); }
#endif
    ierr = PetscViewerASCIIPrintf(viewer," (%10.8e)",(double)errest[i]);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)nep)->tablevel);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   NEPMonitorFirst - Print the first unconverged approximate value and
   error estimate at each iteration of the nonlinear eigensolver.

   Collective on nep

   Input Parameters:
+  nep    - nonlinear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Level: intermediate

.seealso: NEPMonitorSet(), NEPMonitorAll(), NEPMonitorConverged()
@*/
PetscErrorCode NEPMonitorFirst(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(vf,8);
  viewer = vf->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  if (its==1 && ((PetscObject)nep)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)nep)->prefix);CHKERRQ(ierr);
  }
  if (nconv<nest) {
    ierr = PetscViewerPushFormat(viewer,vf->format);CHKERRQ(ierr);
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)nep)->tablevel);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%3D NEP nconv=%D first unconverged value (error)",its,nconv);CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(eigr[nconv]),(double)PetscImaginaryPart(eigr[nconv]));CHKERRQ(ierr);
#else
    ierr = PetscViewerASCIIPrintf(viewer," %g",(double)eigr[nconv]);CHKERRQ(ierr);
    if (eigi[nconv]!=0.0) { ierr = PetscViewerASCIIPrintf(viewer,"%+gi",(double)eigi[nconv]);CHKERRQ(ierr); }
#endif
    ierr = PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[nconv]);CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)nep)->tablevel);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   NEPMonitorConverged - Print the approximate values and
   error estimates as they converge.

   Collective on nep

   Input Parameters:
+  nep    - nonlinear eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  ctx    - monitor context

   Level: intermediate

.seealso: NEPMonitorSet(), NEPMonitorFirst(), NEPMonitorAll()
@*/
PetscErrorCode NEPMonitorConverged(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,SlepcConvMonitor ctx)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(ctx,8);
  viewer = ctx->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  if (its==1 && ((PetscObject)nep)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Convergence history for %s solve.\n",((PetscObject)nep)->prefix);CHKERRQ(ierr);
  }
  if (its==1) ctx->oldnconv = 0;
  if (ctx->oldnconv!=nconv) {
    ierr = PetscViewerPushFormat(viewer,ctx->format);CHKERRQ(ierr);
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)nep)->tablevel);CHKERRQ(ierr);
    for (i=ctx->oldnconv;i<nconv;i++) {
      ierr = PetscViewerASCIIPrintf(viewer,"%3D NEP converged value (error) #%D",its,i);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
      ierr = PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(eigr[i]),(double)PetscImaginaryPart(eigr[i]));CHKERRQ(ierr);
#else
      ierr = PetscViewerASCIIPrintf(viewer," %g",(double)eigr[i]);CHKERRQ(ierr);
      if (eigi[i]!=0.0) { ierr = PetscViewerASCIIPrintf(viewer,"%+gi",(double)eigi[i]);CHKERRQ(ierr); }
#endif
      ierr = PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[i]);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)nep)->tablevel);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ctx->oldnconv = nconv;
  }
  PetscFunctionReturn(0);
}

/*@C
   NEPMonitorLGCreate - Creates a line graph context for use with
   NEP to monitor convergence.

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
+  -nep_monitor_lg - Sets line graph monitor for the first residual
-  -nep_monitor_lg_all - Sets line graph monitor for all residuals

   Notes:
   Use PetscDrawLGDestroy() to destroy this line graph.

   Level: intermediate

.seealso: NEPMonitorSet()
@*/
PetscErrorCode NEPMonitorLGCreate(MPI_Comm comm,const char host[],const char label[],int x,int y,int m,int n,PetscDrawLG *lgctx)
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

PetscErrorCode NEPMonitorLG(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscDrawLG    lg = (PetscDrawLG)ctx;
  PetscReal      x,y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  if (its==1) {
    ierr = PetscDrawLGReset(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSetDimension(lg,1);CHKERRQ(ierr);
    ierr = PetscDrawLGSetLimits(lg,1,1.0,PetscLog10Real(nep->tol)-2,0.0);CHKERRQ(ierr);
  }
  x = (PetscReal)its;
  if (errest[nconv] > 0.0) y = PetscLog10Real(errest[nconv]);
  else y = 0.0;
  ierr = PetscDrawLGAddPoint(lg,&x,&y);CHKERRQ(ierr);
  if (its <= 20 || !(its % 5) || nep->reason) {
    ierr = PetscDrawLGDraw(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSave(lg);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPMonitorLGAll(NEP nep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscDrawLG    lg = (PetscDrawLG)ctx;
  PetscInt       i,n = PetscMin(nep->nev,255);
  PetscReal      *x,*y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  if (its==1) {
    ierr = PetscDrawLGReset(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSetDimension(lg,n);CHKERRQ(ierr);
    ierr = PetscDrawLGSetLimits(lg,1,1.0,PetscLog10Real(nep->tol)-2,0.0);CHKERRQ(ierr);
  }
  ierr = PetscMalloc2(n,&x,n,&y);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
    x[i] = (PetscReal)its;
    if (i < nest && errest[i] > 0.0) y[i] = PetscLog10Real(errest[i]);
    else y[i] = 0.0;
  }
  ierr = PetscDrawLGAddPoint(lg,x,y);CHKERRQ(ierr);
  if (its <= 20 || !(its % 5) || nep->reason) {
    ierr = PetscDrawLGDraw(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSave(lg);CHKERRQ(ierr);
  }
  ierr = PetscFree2(x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

