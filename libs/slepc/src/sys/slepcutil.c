/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/slepcimpl.h>            /*I "slepcsys.h" I*/

/*@C
   SlepcConvMonitorCreate - Creates a SlepcConvMonitor context.

   Collective on viewer

   Input Parameters:
+  viewer - the viewer where the monitor must send data
-  format - the format

   Output Parameter:
.  ctx - the created context

   Notes:
   The created context is used for EPS, SVD, PEP, and NEP monitor functions that just
   print the iteration numbers at which convergence takes place (XXXMonitorConverged).

   This function increases the reference count of the viewer so you can destroy the
   viewer object after this call.

   Level: developer

.seealso: SlepcConvMonitorDestroy()
@*/
PetscErrorCode SlepcConvMonitorCreate(PetscViewer viewer,PetscViewerFormat format,SlepcConvMonitor *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)viewer);CHKERRQ(ierr);
  ierr = PetscNew(ctx);CHKERRQ(ierr);
  (*ctx)->viewer = viewer;
  (*ctx)->format = format;
  PetscFunctionReturn(0);
}

/*@C
   SlepcConvMonitorDestroy - Destroys a SlepcConvMonitor context.

   Collective

   Input Parameters:
.  ctx - the SlepcConvMonitor context to be destroyed.

   Level: developer

.seealso: SlepcConvMonitorCreate()
@*/
PetscErrorCode SlepcConvMonitorDestroy(SlepcConvMonitor *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*ctx) PetscFunctionReturn(0);
  ierr = PetscViewerDestroy(&(*ctx)->viewer);CHKERRQ(ierr);
  ierr = PetscFree(*ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Given n vectors in V, this function gets references of them into W.
   If m<0 then some previous non-processed vectors remain in W and must be freed.
 */
PetscErrorCode SlepcBasisReference_Private(PetscInt n,Vec *V,PetscInt *m,Vec **W)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  for (i=0;i<n;i++) {
    ierr = PetscObjectReference((PetscObject)V[i]);CHKERRQ(ierr);
  }
  ierr = SlepcBasisDestroy_Private(m,W);CHKERRQ(ierr);
  if (n>0) {
    ierr = PetscMalloc1(n,W);CHKERRQ(ierr);
    for (i=0;i<n;i++) (*W)[i] = V[i];
    *m = -n;
  }
  PetscFunctionReturn(0);
}

/*
   Destroys a set of vectors.
   A negative value of m indicates that W contains vectors to be destroyed.
 */
PetscErrorCode SlepcBasisDestroy_Private(PetscInt *m,Vec **W)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  if (*m<0) {
    for (i=0;i<-(*m);i++) {
      ierr = VecDestroy(&(*W)[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree(*W);CHKERRQ(ierr);
  }
  *m = 0;
  PetscFunctionReturn(0);
}

/*@C
   SlepcSNPrintfScalar - Prints a PetscScalar variable to a string of
   given length.

   Not Collective

   Input Parameters:
+  str - the string to print to
.  len - the length of str
.  val - scalar value to be printed
-  exp - to be used within an expression, print leading sign and parentheses
         in case of nonzero imaginary part

   Level: developer
@*/
PetscErrorCode SlepcSNPrintfScalar(char *str,size_t len,PetscScalar val,PetscBool exp)
{
  PetscErrorCode ierr;
#if defined(PETSC_USE_COMPLEX)
  PetscReal      re,im;
#endif

  PetscFunctionBegin;
#if !defined(PETSC_USE_COMPLEX)
  if (exp) {
    ierr = PetscSNPrintf(str,len,"%+g",(double)val);CHKERRQ(ierr);
  } else {
    ierr = PetscSNPrintf(str,len,"%g",(double)val);CHKERRQ(ierr);
  }
#else
  re = PetscRealPart(val);
  im = PetscImaginaryPart(val);
  if (im!=0.0) {
    if (exp) {
      ierr = PetscSNPrintf(str,len,"+(%g%+gi)",(double)re,(double)im);CHKERRQ(ierr);
    } else {
      ierr = PetscSNPrintf(str,len,"%g%+gi",(double)re,(double)im);CHKERRQ(ierr);
    }
  } else {
    if (exp) {
      ierr = PetscSNPrintf(str,len,"%+g",(double)re);CHKERRQ(ierr);
    } else {
      ierr = PetscSNPrintf(str,len,"%g",(double)re);CHKERRQ(ierr);
    }
  }
#endif
  PetscFunctionReturn(0);
}

/*
   SlepcDebugViewMatrix - prints an array as a matrix, to be used from within a debugger.
   Output can be pasted to Matlab.

     nrows, ncols: size of printed matrix
     Xr, Xi: array to be printed (Xi not referenced in complex scalars)
     ldx: leading dimension
     s: name of Matlab variable
     filename: optionally write output to a file
 */
#if defined(PETSC_USE_DEBUG)
PetscErrorCode SlepcDebugViewMatrix(PetscInt nrows,PetscInt ncols,PetscScalar *Xr,PetscScalar *Xi,PetscInt ldx,const char *s,const char *filename)
{
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscViewer    viewer;

  PetscFunctionBegin;
  if (filename) {
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"%s = [\n",s);CHKERRQ(ierr);
  for (i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
#if defined(PETSC_USE_COMPLEX)
      ierr = PetscViewerASCIIPrintf(viewer,"%.18g+%.18gi ",PetscRealPart(Xr[i+j*ldx]),PetscImaginaryPart(Xr[i+j*ldx]));CHKERRQ(ierr);
#else
      if (Xi) { ierr = PetscViewerASCIIPrintf(viewer,"%.18g+%.18gi ",Xr[i+j*ldx],Xi[i+j*ldx]);CHKERRQ(ierr); }
      else { ierr = PetscViewerASCIIPrintf(viewer,"%.18g ",Xr[i+j*ldx]);CHKERRQ(ierr); }
#endif
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"];\n");CHKERRQ(ierr);
  if (filename) { ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}
#endif

/*
   SlepcDebugSetMatlabStdout - sets Matlab format in stdout, to be used from within a debugger.
 */
#if defined(PETSC_USE_DEBUG)
PETSC_UNUSED PetscErrorCode SlepcDebugSetMatlabStdout(void)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;

  PetscFunctionBegin;
  ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

