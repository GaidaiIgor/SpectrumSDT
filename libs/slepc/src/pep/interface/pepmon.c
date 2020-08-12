/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   PEP routines related to monitors
*/

#include <slepc/private/pepimpl.h>      /*I "slepcpep.h" I*/
#include <petscdraw.h>

/*
   Runs the user provided monitor routines, if any.
*/
PetscErrorCode PEPMonitor(PEP pep,PetscInt it,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest)
{
  PetscErrorCode ierr;
  PetscInt       i,n = pep->numbermonitors;

  PetscFunctionBegin;
  for (i=0;i<n;i++) {
    ierr = (*pep->monitor[i])(pep,it,nconv,eigr,eigi,errest,nest,pep->monitorcontext[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   PEPMonitorSet - Sets an ADDITIONAL function to be called at every
   iteration to monitor the error estimates for each requested eigenpair.

   Logically Collective on pep

   Input Parameters:
+  pep     - eigensolver context obtained from PEPCreate()
.  monitor - pointer to function (if this is NULL, it turns off monitoring)
.  mctx    - [optional] context for private data for the
             monitor routine (use NULL if no context is desired)
-  monitordestroy - [optional] routine that frees monitor context (may be NULL)

   Calling Sequence of monitor:
$   monitor(PEP pep,int its,int nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal* errest,int nest,void *mctx)

+  pep    - polynomial eigensolver context obtained from PEPCreate()
.  its    - iteration number
.  nconv  - number of converged eigenpairs
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - relative error estimates for each eigenpair
.  nest   - number of error estimates
-  mctx   - optional monitoring context, as set by PEPMonitorSet()

   Options Database Keys:
+    -pep_monitor        - print only the first error estimate
.    -pep_monitor_all    - print error estimates at each iteration
.    -pep_monitor_conv   - print the eigenvalue approximations only when
      convergence has been reached
.    -pep_monitor_lg     - sets line graph monitor for the first unconverged
      approximate eigenvalue
.    -pep_monitor_lg_all - sets line graph monitor for all unconverged
      approximate eigenvalues
-    -pep_monitor_cancel - cancels all monitors that have been hardwired into
      a code by calls to PEPMonitorSet(), but does not cancel those set via
      the options database.

   Notes:
   Several different monitoring routines may be set by calling
   PEPMonitorSet() multiple times; all will be called in the
   order in which they were set.

   Level: intermediate

.seealso: PEPMonitorFirst(), PEPMonitorAll(), PEPMonitorCancel()
@*/
PetscErrorCode PEPMonitorSet(PEP pep,PetscErrorCode (*monitor)(PEP,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*),void *mctx,PetscErrorCode (*monitordestroy)(void**))
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (pep->numbermonitors >= MAXPEPMONITORS) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_OUTOFRANGE,"Too many PEP monitors set");
  pep->monitor[pep->numbermonitors]           = monitor;
  pep->monitorcontext[pep->numbermonitors]    = (void*)mctx;
  pep->monitordestroy[pep->numbermonitors++]  = monitordestroy;
  PetscFunctionReturn(0);
}

/*@
   PEPMonitorCancel - Clears all monitors for a PEP object.

   Logically Collective on pep

   Input Parameters:
.  pep - eigensolver context obtained from PEPCreate()

   Options Database Key:
.    -pep_monitor_cancel - Cancels all monitors that have been hardwired
      into a code by calls to PEPMonitorSet(),
      but does not cancel those set via the options database.

   Level: intermediate

.seealso: PEPMonitorSet()
@*/
PetscErrorCode PEPMonitorCancel(PEP pep)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  for (i=0; i<pep->numbermonitors; i++) {
    if (pep->monitordestroy[i]) {
      ierr = (*pep->monitordestroy[i])(&pep->monitorcontext[i]);CHKERRQ(ierr);
    }
  }
  pep->numbermonitors = 0;
  PetscFunctionReturn(0);
}

/*@C
   PEPGetMonitorContext - Gets the monitor context, as set by
   PEPMonitorSet() for the FIRST monitor only.

   Not Collective

   Input Parameter:
.  pep - eigensolver context obtained from PEPCreate()

   Output Parameter:
.  ctx - monitor context

   Level: intermediate

.seealso: PEPMonitorSet(), PEPDefaultMonitor()
@*/
PetscErrorCode PEPGetMonitorContext(PEP pep,void **ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  *ctx = pep->monitorcontext[0];
  PetscFunctionReturn(0);
}

/*
   Helper function to compute eigenvalue that must be viewed in monitor
 */
static PetscErrorCode PEPMonitorGetTrueEig(PEP pep,PetscScalar *er,PetscScalar *ei)
{
  PetscErrorCode ierr;
  PetscBool      flg,sinv;

  PetscFunctionBegin;
  ierr = STBackTransform(pep->st,1,er,ei);CHKERRQ(ierr);
  if (pep->sfactor!=1.0) {
    ierr = STGetTransform(pep->st,&flg);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv);CHKERRQ(ierr);
    if (flg && sinv) {
      *er /= pep->sfactor; *ei /= pep->sfactor;
    } else {
      *er *= pep->sfactor; *ei *= pep->sfactor;
    }
  }
  PetscFunctionReturn(0);
}

/*@C
   PEPMonitorAll - Print the current approximate values and
   error estimates at each iteration of the polynomial eigensolver.

   Collective on pep

   Input Parameters:
+  pep    - polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Level: intermediate

.seealso: PEPMonitorSet(), PEPMonitorFirst(), PEPMonitorConverged()
@*/
PetscErrorCode PEPMonitorAll(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    er,ei;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(vf,8);
  viewer = vf->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  ierr = PetscViewerPushFormat(viewer,vf->format);CHKERRQ(ierr);
  ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)pep)->tablevel);CHKERRQ(ierr);
  if (its==1 && ((PetscObject)pep)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)pep)->prefix);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"%3D PEP nconv=%D Values (Errors)",its,nconv);CHKERRQ(ierr);
  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
  for (i=0;i<nest;i++) {
    er = eigr[i]; ei = eigi[i];
    ierr = PEPMonitorGetTrueEig(pep,&er,&ei);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(er),(double)PetscImaginaryPart(er));CHKERRQ(ierr);
#else
    ierr = PetscViewerASCIIPrintf(viewer," %g",(double)er);CHKERRQ(ierr);
    if (eigi[i]!=0.0) { ierr = PetscViewerASCIIPrintf(viewer,"%+gi",(double)ei);CHKERRQ(ierr); }
#endif
    ierr = PetscViewerASCIIPrintf(viewer," (%10.8e)",(double)errest[i]);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)pep)->tablevel);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PEPMonitorFirst - Print the first unconverged approximate value and
   error estimate at each iteration of the polynomial eigensolver.

   Collective on pep

   Input Parameters:
+  pep    - polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  vf     - viewer and format for monitoring

   Level: intermediate

.seealso: PEPMonitorSet(), PEPMonitorAll(), PEPMonitorConverged()
@*/
PetscErrorCode PEPMonitorFirst(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,PetscViewerAndFormat *vf)
{
  PetscErrorCode ierr;
  PetscScalar    er,ei;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(vf,8);
  viewer = vf->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  if (its==1 && ((PetscObject)pep)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Eigenvalue approximations and residual norms for %s solve.\n",((PetscObject)pep)->prefix);CHKERRQ(ierr);
  }
  if (nconv<nest) {
    ierr = PetscViewerPushFormat(viewer,vf->format);CHKERRQ(ierr);
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)pep)->tablevel);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%3D PEP nconv=%D first unconverged value (error)",its,nconv);CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    er = eigr[nconv]; ei = eigi[nconv];
    ierr = PEPMonitorGetTrueEig(pep,&er,&ei);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(er),(double)PetscImaginaryPart(er));CHKERRQ(ierr);
#else
    ierr = PetscViewerASCIIPrintf(viewer," %g",(double)er);CHKERRQ(ierr);
    if (eigi[nconv]!=0.0) { ierr = PetscViewerASCIIPrintf(viewer,"%+gi",(double)ei);CHKERRQ(ierr); }
#endif
    ierr = PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[nconv]);CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)pep)->tablevel);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   PEPMonitorConverged - Print the approximate values and
   error estimates as they converge.

   Collective on pep

   Input Parameters:
+  pep    - polynomial eigensolver context
.  its    - iteration number
.  nconv  - number of converged eigenpairs so far
.  eigr   - real part of the eigenvalues
.  eigi   - imaginary part of the eigenvalues
.  errest - error estimates
.  nest   - number of error estimates to display
-  ctx    - monitor context

   Level: intermediate

.seealso: PEPMonitorSet(), PEPMonitorFirst(), PEPMonitorAll()
@*/
PetscErrorCode PEPMonitorConverged(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,SlepcConvMonitor ctx)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    er,ei;
  PetscViewer    viewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  PetscValidPointer(ctx,8);
  viewer = ctx->viewer;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,8);
  if (its==1 && ((PetscObject)pep)->prefix) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Convergence history for %s solve.\n",((PetscObject)pep)->prefix);CHKERRQ(ierr);
  }
  if (its==1) ctx->oldnconv = 0;
  if (ctx->oldnconv!=nconv) {
    ierr = PetscViewerPushFormat(viewer,ctx->format);CHKERRQ(ierr);
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)pep)->tablevel);CHKERRQ(ierr);
    for (i=ctx->oldnconv;i<nconv;i++) {
      ierr = PetscViewerASCIIPrintf(viewer,"%3D PEP converged value (error) #%D",its,i);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
      er = eigr[i]; ei = eigi[i];
      ierr = PEPMonitorGetTrueEig(pep,&er,&ei);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
      ierr = PetscViewerASCIIPrintf(viewer," %g%+gi",(double)PetscRealPart(er),(double)PetscImaginaryPart(er));CHKERRQ(ierr);
#else
      ierr = PetscViewerASCIIPrintf(viewer," %g",(double)er);CHKERRQ(ierr);
      if (eigi[i]!=0.0) { ierr = PetscViewerASCIIPrintf(viewer,"%+gi",(double)ei);CHKERRQ(ierr); }
#endif
      ierr = PetscViewerASCIIPrintf(viewer," (%10.8e)\n",(double)errest[i]);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)pep)->tablevel);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ctx->oldnconv = nconv;
  }
  PetscFunctionReturn(0);
}

/*@C
   PEPMonitorLGCreate - Creates a line graph context for use with
   PEP to monitor convergence.

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
+  -pep_monitor_lg - Sets line graph monitor for the first residual
-  -pep_monitor_lg_all - Sets line graph monitor for all residuals

   Notes:
   Use PetscDrawLGDestroy() to destroy this line graph.

   Level: intermediate

.seealso: PEPMonitorSet()
@*/
PetscErrorCode PEPMonitorLGCreate(MPI_Comm comm,const char host[],const char label[],int x,int y,int m,int n,PetscDrawLG *lgctx)
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

PetscErrorCode PEPMonitorLG(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscDrawLG    lg = (PetscDrawLG)ctx;
  PetscReal      x,y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  if (its==1) {
    ierr = PetscDrawLGReset(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSetDimension(lg,1);CHKERRQ(ierr);
    ierr = PetscDrawLGSetLimits(lg,1,1.0,PetscLog10Real(pep->tol)-2,0.0);CHKERRQ(ierr);
  }
  x = (PetscReal)its;
  if (errest[nconv] > 0.0) y = PetscLog10Real(errest[nconv]);
  else y = 0.0;
  ierr = PetscDrawLGAddPoint(lg,&x,&y);CHKERRQ(ierr);
  if (its <= 20 || !(its % 5) || pep->reason) {
    ierr = PetscDrawLGDraw(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSave(lg);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPMonitorLGAll(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscDrawLG    lg = (PetscDrawLG)ctx;
  PetscInt       i,n = PetscMin(pep->nev,255);
  PetscReal      *x,*y;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(lg,PETSC_DRAWLG_CLASSID,8);
  if (its==1) {
    ierr = PetscDrawLGReset(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSetDimension(lg,n);CHKERRQ(ierr);
    ierr = PetscDrawLGSetLimits(lg,1,1.0,PetscLog10Real(pep->tol)-2,0.0);CHKERRQ(ierr);
  }
  ierr = PetscMalloc2(n,&x,n,&y);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
    x[i] = (PetscReal)its;
    if (i < nest && errest[i] > 0.0) y[i] = PetscLog10Real(errest[i]);
    else y[i] = 0.0;
  }
  ierr = PetscDrawLGAddPoint(lg,x,y);CHKERRQ(ierr);
  if (its <= 20 || !(its % 5) || pep->reason) {
    ierr = PetscDrawLGDraw(lg);CHKERRQ(ierr);
    ierr = PetscDrawLGSave(lg);CHKERRQ(ierr);
  }
  ierr = PetscFree2(x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

