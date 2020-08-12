/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   NEP routines related to various viewers
*/

#include <slepc/private/nepimpl.h>      /*I "slepcnep.h" I*/
#include <petscdraw.h>

/*@C
   NEPView - Prints the NEP data structure.

   Collective on nep

   Input Parameters:
+  nep - the nonlinear eigenproblem solver context
-  viewer - optional visualization context

   Options Database Key:
.  -nep_view -  Calls NEPView() at end of NEPSolve()

   Note:
   The available visualization contexts include
+     PETSC_VIEWER_STDOUT_SELF - standard output (default)
-     PETSC_VIEWER_STDOUT_WORLD - synchronized standard
         output where only the first processor opens
         the file.  All other processors send their
         data to the first processor to print.

   The user can open an alternative visualization context with
   PetscViewerASCIIOpen() - output to a specified file.

   Level: beginner
@*/
PetscErrorCode NEPView(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  const char     *type=NULL;
  char           str[50];
  PetscInt       i;
  PetscBool      isascii,istrivial;
  PetscViewer    sviewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)nep),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(nep,1,viewer,2);

  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscObjectPrintClassNamePrefixType((PetscObject)nep,viewer);CHKERRQ(ierr);
    if (nep->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*nep->ops->view)(nep,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    if (nep->problem_type) {
      switch (nep->problem_type) {
        case NEP_GENERAL:  type = "general nonlinear eigenvalue problem"; break;
        case NEP_RATIONAL: type = "rational eigenvalue problem"; break;
      }
    } else type = "not yet set";
    ierr = PetscViewerASCIIPrintf(viewer,"  problem type: %s\n",type);CHKERRQ(ierr);
    if (nep->fui) {
      switch (nep->fui) {
      case NEP_USER_INTERFACE_CALLBACK:
        ierr = PetscViewerASCIIPrintf(viewer,"  nonlinear operator from user callbacks\n");CHKERRQ(ierr);
        break;
      case NEP_USER_INTERFACE_SPLIT:
        ierr = PetscViewerASCIIPrintf(viewer,"  nonlinear operator in split form, with %D terms\n",nep->nt);CHKERRQ(ierr);
        break;
      }
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"  nonlinear operator not specified yet\n");CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"  selected portion of the spectrum: ");CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    ierr = SlepcSNPrintfScalar(str,50,nep->target,PETSC_FALSE);CHKERRQ(ierr);
    if (!nep->which) {
      ierr = PetscViewerASCIIPrintf(viewer,"not yet set\n");CHKERRQ(ierr);
    } else switch (nep->which) {
      case NEP_WHICH_USER:
        ierr = PetscViewerASCIIPrintf(viewer,"user defined\n");CHKERRQ(ierr);
        break;
      case NEP_TARGET_MAGNITUDE:
        ierr = PetscViewerASCIIPrintf(viewer,"closest to target: %s (in magnitude)\n",str);CHKERRQ(ierr);
        break;
      case NEP_TARGET_REAL:
        ierr = PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the real axis)\n",str);CHKERRQ(ierr);
        break;
      case NEP_TARGET_IMAGINARY:
        ierr = PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the imaginary axis)\n",str);CHKERRQ(ierr);
        break;
      case NEP_LARGEST_MAGNITUDE:
        ierr = PetscViewerASCIIPrintf(viewer,"largest eigenvalues in magnitude\n");CHKERRQ(ierr);
        break;
      case NEP_SMALLEST_MAGNITUDE:
        ierr = PetscViewerASCIIPrintf(viewer,"smallest eigenvalues in magnitude\n");CHKERRQ(ierr);
        break;
      case NEP_LARGEST_REAL:
        ierr = PetscViewerASCIIPrintf(viewer,"largest real parts\n");CHKERRQ(ierr);
        break;
      case NEP_SMALLEST_REAL:
        ierr = PetscViewerASCIIPrintf(viewer,"smallest real parts\n");CHKERRQ(ierr);
        break;
      case NEP_LARGEST_IMAGINARY:
        ierr = PetscViewerASCIIPrintf(viewer,"largest imaginary parts\n");CHKERRQ(ierr);
        break;
      case NEP_SMALLEST_IMAGINARY:
        ierr = PetscViewerASCIIPrintf(viewer,"smallest imaginary parts\n");CHKERRQ(ierr);
        break;
      case NEP_ALL:
        ierr = PetscViewerASCIIPrintf(viewer,"all eigenvalues in the region\n");CHKERRQ(ierr);
        break;
    }
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    if (nep->twosided) {
      ierr = PetscViewerASCIIPrintf(viewer,"  using two-sided variant (for left eigenvectors)\n");CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"  number of eigenvalues (nev): %D\n",nep->nev);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  number of column vectors (ncv): %D\n",nep->ncv);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum dimension of projected problem (mpd): %D\n",nep->mpd);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum number of iterations: %D\n",nep->max_it);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  tolerance: %g\n",(double)nep->tol);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  convergence test: ");CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    switch (nep->conv) {
    case NEP_CONV_ABS:
      ierr = PetscViewerASCIIPrintf(viewer,"absolute\n");CHKERRQ(ierr);break;
    case NEP_CONV_REL:
      ierr = PetscViewerASCIIPrintf(viewer,"relative to the eigenvalue\n");CHKERRQ(ierr);break;
    case NEP_CONV_NORM:
      ierr = PetscViewerASCIIPrintf(viewer,"relative to the matrix norms\n");CHKERRQ(ierr);
      if (nep->nrma) {
        ierr = PetscViewerASCIIPrintf(viewer,"  computed matrix norms: %g",(double)nep->nrma[0]);CHKERRQ(ierr);
        for (i=1;i<nep->nt;i++) {
          ierr = PetscViewerASCIIPrintf(viewer,", %g",(double)nep->nrma[i]);CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
      }
      break;
    case NEP_CONV_USER:
      ierr = PetscViewerASCIIPrintf(viewer,"user-defined\n");CHKERRQ(ierr);break;
    }
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    if (nep->refine) {
      ierr = PetscViewerASCIIPrintf(viewer,"  iterative refinement: %s, with %s scheme\n",NEPRefineTypes[nep->refine],NEPRefineSchemes[nep->scheme]);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"  refinement stopping criterion: tol=%g, its=%D\n",(double)nep->rtol,nep->rits);CHKERRQ(ierr);
      if (nep->npart>1) {
        ierr = PetscViewerASCIIPrintf(viewer,"  splitting communicator in %D partitions for refinement\n",nep->npart);CHKERRQ(ierr);
      }
    }
    if (nep->nini) {
      ierr = PetscViewerASCIIPrintf(viewer,"  dimension of user-provided initial space: %D\n",PetscAbs(nep->nini));CHKERRQ(ierr);
    }
  } else {
    if (nep->ops->view) {
      ierr = (*nep->ops->view)(nep,viewer);CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
  if (!nep->V) { ierr = NEPGetBV(nep,&nep->V);CHKERRQ(ierr); }
  ierr = BVView(nep->V,viewer);CHKERRQ(ierr);
  if (!nep->rg) { ierr = NEPGetRG(nep,&nep->rg);CHKERRQ(ierr); }
  ierr = RGIsTrivial(nep->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) { ierr = RGView(nep->rg,viewer);CHKERRQ(ierr); }
  if (nep->useds) {
    if (!nep->ds) { ierr = NEPGetDS(nep,&nep->ds);CHKERRQ(ierr); }
    ierr = DSView(nep->ds,viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  if (nep->refine!=NEP_REFINE_NONE) {
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    if (nep->npart>1) {
      if (nep->refinesubc->color==0) {
        ierr = PetscViewerASCIIGetStdout(PetscSubcommChild(nep->refinesubc),&sviewer);CHKERRQ(ierr);
        ierr = KSPView(nep->refineksp,sviewer);CHKERRQ(ierr);
      }
    } else {
      ierr = KSPView(nep->refineksp,viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   NEPViewFromOptions - View from options

   Collective on NEP

   Input Parameters:
+  nep  - the nonlinear eigensolver context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: NEPView(), NEPCreate()
@*/
PetscErrorCode NEPViewFromOptions(NEP nep,PetscObject obj,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  ierr = PetscObjectViewFromOptions((PetscObject)nep,obj,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   NEPReasonView - Displays the reason a NEP solve converged or diverged.

   Collective on nep

   Input Parameters:
+  nep - the nonlinear eigensolver context
-  viewer - the viewer to display the reason

   Options Database Keys:
.  -nep_converged_reason - print reason for convergence, and number of iterations

   Level: intermediate

.seealso: NEPSetConvergenceTest(), NEPSetTolerances(), NEPGetIterationNumber()
@*/
PetscErrorCode NEPReasonView(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isAscii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isAscii);CHKERRQ(ierr);
  if (isAscii) {
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)nep)->tablevel);CHKERRQ(ierr);
    if (nep->reason > 0) {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Nonlinear eigensolve converged (%D eigenpair%s) due to %s; iterations %D\n",((PetscObject)nep)->prefix?((PetscObject)nep)->prefix:"",nep->nconv,(nep->nconv>1)?"s":"",NEPConvergedReasons[nep->reason],nep->its);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Nonlinear eigensolve did not converge due to %s; iterations %D\n",((PetscObject)nep)->prefix?((PetscObject)nep)->prefix:"",NEPConvergedReasons[nep->reason],nep->its);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)nep)->tablevel);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   NEPReasonViewFromOptions - Processes command line options to determine if/how
   the NEP converged reason is to be viewed.

   Collective on nep

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Level: developer
@*/
PetscErrorCode NEPReasonViewFromOptions(NEP nep)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_converged_reason",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = NEPReasonView(nep,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPErrorView_ASCII(NEP nep,NEPErrorType etype,PetscViewer viewer)
{
  PetscReal      error;
  PetscInt       i,j,k,nvals;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  nvals = (nep->which==NEP_ALL)? nep->nconv: nep->nev;
  if (nep->which!=NEP_ALL && nep->nconv<nvals) {
    ierr = PetscViewerASCIIPrintf(viewer," Problem: less than %D eigenvalues converged\n\n",nep->nev);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  if (nep->which==NEP_ALL && !nvals) {
    ierr = PetscViewerASCIIPrintf(viewer," No eigenvalues have been found\n\n");CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  for (i=0;i<nvals;i++) {
    ierr = NEPComputeError(nep,i,etype,&error);CHKERRQ(ierr);
    if (error>=5.0*nep->tol) {
      ierr = PetscViewerASCIIPrintf(viewer," Problem: some of the first %D relative errors are higher than the tolerance\n\n",nvals);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }
  if (nep->which==NEP_ALL) {
    ierr = PetscViewerASCIIPrintf(viewer," Found %D eigenvalues, all of them computed up to the required tolerance:",nvals);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerASCIIPrintf(viewer," All requested eigenvalues computed up to the required tolerance:");CHKERRQ(ierr);
  }
  for (i=0;i<=(nvals-1)/8;i++) {
    ierr = PetscViewerASCIIPrintf(viewer,"\n     ");CHKERRQ(ierr);
    for (j=0;j<PetscMin(8,nvals-8*i);j++) {
      k = nep->perm[8*i+j];
      ierr = SlepcPrintEigenvalueASCII(viewer,nep->eigr[k],nep->eigi[k]);CHKERRQ(ierr);
      if (8*i+j+1<nvals) { ierr = PetscViewerASCIIPrintf(viewer,", ");CHKERRQ(ierr); }
    }
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPErrorView_DETAIL(NEP nep,NEPErrorType etype,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscReal      error,re,im;
  PetscScalar    kr,ki;
  PetscInt       i;
#define EXLEN 30
  char           ex[EXLEN],sep[]=" ---------------------- --------------------\n";

  PetscFunctionBegin;
  if (!nep->nconv) PetscFunctionReturn(0);
  switch (etype) {
    case NEP_ERROR_ABSOLUTE:
      ierr = PetscSNPrintf(ex,EXLEN,"    ||T(k)x||");CHKERRQ(ierr);
      break;
    case NEP_ERROR_RELATIVE:
      ierr = PetscSNPrintf(ex,EXLEN," ||T(k)x||/||kx||");CHKERRQ(ierr);
      break;
    case NEP_ERROR_BACKWARD:
      ierr = PetscSNPrintf(ex,EXLEN,"    eta(x,k)");CHKERRQ(ierr);
      break;
  }
  ierr = PetscViewerASCIIPrintf(viewer,"%s            k             %s\n%s",sep,ex,sep);CHKERRQ(ierr);
  for (i=0;i<nep->nconv;i++) {
    ierr = NEPGetEigenpair(nep,i,&kr,&ki,NULL,NULL);CHKERRQ(ierr);
    ierr = NEPComputeError(nep,i,etype,&error);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(kr);
    im = PetscImaginaryPart(kr);
#else
    re = kr;
    im = ki;
#endif
    if (im!=0.0) {
      ierr = PetscViewerASCIIPrintf(viewer,"  % 9f%+9fi      %12g\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"    % 12f           %12g\n",(double)re,(double)error);CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerASCIIPrintf(viewer,sep);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPErrorView_MATLAB(NEP nep,NEPErrorType etype,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscReal      error;
  PetscInt       i;
  const char     *name;

  PetscFunctionBegin;
  ierr = PetscObjectGetName((PetscObject)nep,&name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Error_%s = [\n",name);CHKERRQ(ierr);
  for (i=0;i<nep->nconv;i++) {
    ierr = NEPComputeError(nep,i,etype,&error);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)error);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"];\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   NEPErrorView - Displays the errors associated with the computed solution
   (as well as the eigenvalues).

   Collective on nep

   Input Parameters:
+  nep    - the nonlinear eigensolver context
.  etype  - error type
-  viewer - optional visualization context

   Options Database Key:
+  -nep_error_absolute - print absolute errors of each eigenpair
.  -nep_error_relative - print relative errors of each eigenpair
-  -nep_error_backward - print backward errors of each eigenpair

   Notes:
   By default, this function checks the error of all eigenpairs and prints
   the eigenvalues if all of them are below the requested tolerance.
   If the viewer has format=PETSC_VIEWER_ASCII_INFO_DETAIL then a table with
   eigenvalues and corresponding errors is printed.

   Level: intermediate

.seealso: NEPSolve(), NEPValuesView(), NEPVectorsView()
@*/
PetscErrorCode NEPErrorView(NEP nep,NEPErrorType etype,PetscViewer viewer)
{
  PetscBool         isascii;
  PetscViewerFormat format;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)nep),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(nep,1,viewer,2);
  NEPCheckSolved(nep,1);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (!isascii) PetscFunctionReturn(0);

  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  switch (format) {
    case PETSC_VIEWER_DEFAULT:
    case PETSC_VIEWER_ASCII_INFO:
      ierr = NEPErrorView_ASCII(nep,etype,viewer);CHKERRQ(ierr);
      break;
    case PETSC_VIEWER_ASCII_INFO_DETAIL:
      ierr = NEPErrorView_DETAIL(nep,etype,viewer);CHKERRQ(ierr);
      break;
    case PETSC_VIEWER_ASCII_MATLAB:
      ierr = NEPErrorView_MATLAB(nep,etype,viewer);CHKERRQ(ierr);
      break;
    default:
      ierr = PetscInfo1(nep,"Unsupported viewer format %s\n",PetscViewerFormats[format]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   NEPErrorViewFromOptions - Processes command line options to determine if/how
   the errors of the computed solution are to be viewed.

   Collective on nep

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Level: developer
@*/
PetscErrorCode NEPErrorViewFromOptions(NEP nep)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_error_absolute",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = NEPErrorView(nep,NEP_ERROR_ABSOLUTE,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_error_relative",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = NEPErrorView(nep,NEP_ERROR_RELATIVE,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_error_backward",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = NEPErrorView(nep,NEP_ERROR_BACKWARD,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPValuesView_DRAW(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscDraw      draw;
  PetscDrawSP    drawsp;
  PetscReal      re,im;
  PetscInt       i,k;

  PetscFunctionBegin;
  if (!nep->nconv) PetscFunctionReturn(0);
  ierr = PetscViewerDrawGetDraw(viewer,0,&draw);CHKERRQ(ierr);
  ierr = PetscDrawSetTitle(draw,"Computed Eigenvalues");CHKERRQ(ierr);
  ierr = PetscDrawSPCreate(draw,1,&drawsp);CHKERRQ(ierr);
  for (i=0;i<nep->nconv;i++) {
    k = nep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(nep->eigr[k]);
    im = PetscImaginaryPart(nep->eigr[k]);
#else
    re = nep->eigr[k];
    im = nep->eigi[k];
#endif
    ierr = PetscDrawSPAddPoint(drawsp,&re,&im);CHKERRQ(ierr);
  }
  ierr = PetscDrawSPDraw(drawsp,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscDrawSPSave(drawsp);CHKERRQ(ierr);
  ierr = PetscDrawSPDestroy(&drawsp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPValuesView_ASCII(NEP nep,PetscViewer viewer)
{
  PetscInt       i,k;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"Eigenvalues = \n");CHKERRQ(ierr);
  for (i=0;i<nep->nconv;i++) {
    k = nep->perm[i];
    ierr = PetscViewerASCIIPrintf(viewer,"   ");CHKERRQ(ierr);
    ierr = SlepcPrintEigenvalueASCII(viewer,nep->eigr[k],nep->eigi[k]);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPValuesView_MATLAB(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  PetscReal      re,im;
  const char     *name;

  PetscFunctionBegin;
  ierr = PetscObjectGetName((PetscObject)nep,&name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Lambda_%s = [\n",name);CHKERRQ(ierr);
  for (i=0;i<nep->nconv;i++) {
    k = nep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(nep->eigr[k]);
    im = PetscImaginaryPart(nep->eigr[k]);
#else
    re = nep->eigr[k];
    im = nep->eigi[k];
#endif
    if (im!=0.0) {
      ierr = PetscViewerASCIIPrintf(viewer,"%18.16e%+18.16ei\n",(double)re,(double)im);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)re);CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerASCIIPrintf(viewer,"];\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   NEPValuesView - Displays the computed eigenvalues in a viewer.

   Collective on nep

   Input Parameters:
+  nep    - the nonlinear eigensolver context
-  viewer - the viewer

   Options Database Key:
.  -nep_view_values - print computed eigenvalues

   Level: intermediate

.seealso: NEPSolve(), NEPVectorsView(), NEPErrorView()
@*/
PetscErrorCode NEPValuesView(NEP nep,PetscViewer viewer)
{
  PetscBool         isascii,isdraw;
  PetscViewerFormat format;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)nep),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(nep,1,viewer,2);
  NEPCheckSolved(nep,1);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isdraw) {
    ierr = NEPValuesView_DRAW(nep,viewer);CHKERRQ(ierr);
  } else if (isascii) {
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    switch (format) {
      case PETSC_VIEWER_DEFAULT:
      case PETSC_VIEWER_ASCII_INFO:
      case PETSC_VIEWER_ASCII_INFO_DETAIL:
        ierr = NEPValuesView_ASCII(nep,viewer);CHKERRQ(ierr);
        break;
      case PETSC_VIEWER_ASCII_MATLAB:
        ierr = NEPValuesView_MATLAB(nep,viewer);CHKERRQ(ierr);
        break;
      default:
        ierr = PetscInfo1(nep,"Unsupported viewer format %s\n",PetscViewerFormats[format]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@
   NEPValuesViewFromOptions - Processes command line options to determine if/how
   the computed eigenvalues are to be viewed.

   Collective on nep

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Level: developer
@*/
PetscErrorCode NEPValuesViewFromOptions(NEP nep)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_view_values",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = NEPValuesView(nep,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
   NEPVectorsView - Outputs computed eigenvectors to a viewer.

   Collective on nep

   Input Parameters:
+  nep    - the nonlinear eigensolver context
-  viewer - the viewer

   Options Database Keys:
.  -nep_view_vectors - output eigenvectors.

   Notes:
   If PETSc was configured with real scalars, complex conjugate eigenvectors
   will be viewed as two separate real vectors, one containing the real part
   and another one containing the imaginary part.

   If left eigenvectors were computed with a two-sided eigensolver, the right
   and left eigenvectors are interleaved, that is, the vectors are output in
   the following order: X0, Y0, X1, Y1, X2, Y2, ...

   Level: intermediate

.seealso: NEPSolve(), NEPValuesView(), NEPErrorView()
@*/
PetscErrorCode NEPVectorsView(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  Vec            x;
#define NMLEN 30
  char           vname[NMLEN];
  const char     *ename;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)nep),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(nep,1,viewer,2);
  NEPCheckSolved(nep,1);
  if (nep->nconv) {
    ierr = PetscObjectGetName((PetscObject)nep,&ename);CHKERRQ(ierr);
    ierr = NEPComputeVectors(nep);CHKERRQ(ierr);
    for (i=0;i<nep->nconv;i++) {
      k = nep->perm[i];
      ierr = PetscSNPrintf(vname,NMLEN,"X%d_%s",(int)i,ename);CHKERRQ(ierr);
      ierr = BVGetColumn(nep->V,k,&x);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)x,vname);CHKERRQ(ierr);
      ierr = VecView(x,viewer);CHKERRQ(ierr);
      ierr = BVRestoreColumn(nep->V,k,&x);CHKERRQ(ierr);
      if (nep->twosided) {
        ierr = PetscSNPrintf(vname,NMLEN,"Y%d_%s",(int)i,ename);CHKERRQ(ierr);
        ierr = BVGetColumn(nep->W,k,&x);CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)x,vname);CHKERRQ(ierr);
        ierr = VecView(x,viewer);CHKERRQ(ierr);
        ierr = BVRestoreColumn(nep->W,k,&x);CHKERRQ(ierr);
      }
    }
  }
  PetscFunctionReturn(0);
}

/*@
   NEPVectorsViewFromOptions - Processes command line options to determine if/how
   the computed eigenvectors are to be viewed.

   Collective on nep

   Input Parameter:
.  nep - the nonlinear eigensolver context

   Level: developer
@*/
PetscErrorCode NEPVectorsViewFromOptions(NEP nep)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg = PETSC_FALSE;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)nep),((PetscObject)nep)->options,((PetscObject)nep)->prefix,"-nep_view_vectors",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = NEPVectorsView(nep,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

