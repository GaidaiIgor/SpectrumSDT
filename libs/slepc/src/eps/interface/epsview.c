/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   EPS routines related to various viewers
*/

#include <slepc/private/epsimpl.h>      /*I "slepceps.h" I*/
#include <petscdraw.h>

/*@C
   EPSView - Prints the EPS data structure.

   Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  viewer - optional visualization context

   Options Database Key:
.  -eps_view -  Calls EPSView() at end of EPSSolve()

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

.seealso: STView()
@*/
PetscErrorCode EPSView(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  const char     *type=NULL,*extr=NULL,*bal=NULL;
  char           str[50];
  PetscBool      isascii,isexternal,istrivial;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)eps),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(eps,1,viewer,2);

#if defined(PETSC_USE_COMPLEX)
#define HERM "hermitian"
#else
#define HERM "symmetric"
#endif
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscObjectPrintClassNamePrefixType((PetscObject)eps,viewer);CHKERRQ(ierr);
    if (eps->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*eps->ops->view)(eps,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    if (eps->problem_type) {
      switch (eps->problem_type) {
        case EPS_HEP:    type = HERM " eigenvalue problem"; break;
        case EPS_GHEP:   type = "generalized " HERM " eigenvalue problem"; break;
        case EPS_NHEP:   type = "non-" HERM " eigenvalue problem"; break;
        case EPS_GNHEP:  type = "generalized non-" HERM " eigenvalue problem"; break;
        case EPS_PGNHEP: type = "generalized non-" HERM " eigenvalue problem with " HERM " positive definite B"; break;
        case EPS_GHIEP:  type = "generalized " HERM "-indefinite eigenvalue problem"; break;
      }
    } else type = "not yet set";
    ierr = PetscViewerASCIIPrintf(viewer,"  problem type: %s\n",type);CHKERRQ(ierr);
    if (eps->extraction) {
      switch (eps->extraction) {
        case EPS_RITZ:              extr = "Rayleigh-Ritz"; break;
        case EPS_HARMONIC:          extr = "harmonic Ritz"; break;
        case EPS_HARMONIC_RELATIVE: extr = "relative harmonic Ritz"; break;
        case EPS_HARMONIC_RIGHT:    extr = "right harmonic Ritz"; break;
        case EPS_HARMONIC_LARGEST:  extr = "largest harmonic Ritz"; break;
        case EPS_REFINED:           extr = "refined Ritz"; break;
        case EPS_REFINED_HARMONIC:  extr = "refined harmonic Ritz"; break;
      }
      ierr = PetscViewerASCIIPrintf(viewer,"  extraction type: %s\n",extr);CHKERRQ(ierr);
    }
    if (!eps->ishermitian && eps->balance!=EPS_BALANCE_NONE) {
      switch (eps->balance) {
        case EPS_BALANCE_NONE:    break;
        case EPS_BALANCE_ONESIDE: bal = "one-sided Krylov"; break;
        case EPS_BALANCE_TWOSIDE: bal = "two-sided Krylov"; break;
        case EPS_BALANCE_USER:    bal = "user-defined matrix"; break;
      }
      ierr = PetscViewerASCIIPrintf(viewer,"  balancing enabled: %s",bal);CHKERRQ(ierr);
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
      if (eps->balance==EPS_BALANCE_ONESIDE || eps->balance==EPS_BALANCE_TWOSIDE) {
        ierr = PetscViewerASCIIPrintf(viewer,", with its=%D",eps->balance_its);CHKERRQ(ierr);
      }
      if (eps->balance==EPS_BALANCE_TWOSIDE && eps->balance_cutoff!=0.0) {
        ierr = PetscViewerASCIIPrintf(viewer," and cutoff=%g",(double)eps->balance_cutoff);CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"  selected portion of the spectrum: ");CHKERRQ(ierr);
    ierr = SlepcSNPrintfScalar(str,50,eps->target,PETSC_FALSE);CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    if (!eps->which) {
      ierr = PetscViewerASCIIPrintf(viewer,"not yet set\n");CHKERRQ(ierr);
    } else switch (eps->which) {
      case EPS_WHICH_USER:
        ierr = PetscViewerASCIIPrintf(viewer,"user defined\n");CHKERRQ(ierr);
        break;
      case EPS_TARGET_MAGNITUDE:
        ierr = PetscViewerASCIIPrintf(viewer,"closest to target: %s (in magnitude)\n",str);CHKERRQ(ierr);
        break;
      case EPS_TARGET_REAL:
        ierr = PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the real axis)\n",str);CHKERRQ(ierr);
        break;
      case EPS_TARGET_IMAGINARY:
        ierr = PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the imaginary axis)\n",str);CHKERRQ(ierr);
        break;
      case EPS_LARGEST_MAGNITUDE:
        ierr = PetscViewerASCIIPrintf(viewer,"largest eigenvalues in magnitude\n");CHKERRQ(ierr);
        break;
      case EPS_SMALLEST_MAGNITUDE:
        ierr = PetscViewerASCIIPrintf(viewer,"smallest eigenvalues in magnitude\n");CHKERRQ(ierr);
        break;
      case EPS_LARGEST_REAL:
        ierr = PetscViewerASCIIPrintf(viewer,"largest real parts\n");CHKERRQ(ierr);
        break;
      case EPS_SMALLEST_REAL:
        ierr = PetscViewerASCIIPrintf(viewer,"smallest real parts\n");CHKERRQ(ierr);
        break;
      case EPS_LARGEST_IMAGINARY:
        ierr = PetscViewerASCIIPrintf(viewer,"largest imaginary parts\n");CHKERRQ(ierr);
        break;
      case EPS_SMALLEST_IMAGINARY:
        ierr = PetscViewerASCIIPrintf(viewer,"smallest imaginary parts\n");CHKERRQ(ierr);
        break;
      case EPS_ALL:
        if (eps->inta || eps->intb) {
          ierr = PetscViewerASCIIPrintf(viewer,"all eigenvalues in interval [%g,%g]\n",(double)eps->inta,(double)eps->intb);CHKERRQ(ierr);
        } else {
          ierr = PetscViewerASCIIPrintf(viewer,"all eigenvalues in the region\n");CHKERRQ(ierr);
        }
        break;
    }
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    if (eps->twosided && eps->problem_type!=EPS_HEP && eps->problem_type!=EPS_GHEP) {
      ierr = PetscViewerASCIIPrintf(viewer,"  using two-sided variant (for left eigenvectors)\n");CHKERRQ(ierr);
    }
    if (eps->purify) {
      ierr = PetscViewerASCIIPrintf(viewer,"  postprocessing eigenvectors with purification\n");CHKERRQ(ierr);
    }
    if (eps->trueres) {
      ierr = PetscViewerASCIIPrintf(viewer,"  computing true residuals explicitly\n");CHKERRQ(ierr);
    }
    if (eps->trackall) {
      ierr = PetscViewerASCIIPrintf(viewer,"  computing all residuals (for tracking convergence)\n");CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"  number of eigenvalues (nev): %D\n",eps->nev);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  number of column vectors (ncv): %D\n",eps->ncv);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum dimension of projected problem (mpd): %D\n",eps->mpd);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum number of iterations: %D\n",eps->max_it);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  tolerance: %g\n",(double)eps->tol);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  convergence test: ");CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    switch (eps->conv) {
    case EPS_CONV_ABS:
      ierr = PetscViewerASCIIPrintf(viewer,"absolute\n");CHKERRQ(ierr);break;
    case EPS_CONV_REL:
      ierr = PetscViewerASCIIPrintf(viewer,"relative to the eigenvalue\n");CHKERRQ(ierr);break;
    case EPS_CONV_NORM:
      ierr = PetscViewerASCIIPrintf(viewer,"relative to the eigenvalue and matrix norms\n");CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"  computed matrix norms: norm(A)=%g",(double)eps->nrma);CHKERRQ(ierr);
      if (eps->isgeneralized) {
        ierr = PetscViewerASCIIPrintf(viewer,", norm(B)=%g",(double)eps->nrmb);CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
      break;
    case EPS_CONV_USER:
      ierr = PetscViewerASCIIPrintf(viewer,"user-defined\n");CHKERRQ(ierr);break;
    }
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    if (eps->nini) {
      ierr = PetscViewerASCIIPrintf(viewer,"  dimension of user-provided initial space: %D\n",PetscAbs(eps->nini));CHKERRQ(ierr);
    }
    if (eps->ninil) {
      ierr = PetscViewerASCIIPrintf(viewer,"  dimension of user-provided left initial space: %D\n",PetscAbs(eps->ninil));CHKERRQ(ierr);
    }
    if (eps->nds) {
      ierr = PetscViewerASCIIPrintf(viewer,"  dimension of user-provided deflation space: %D\n",PetscAbs(eps->nds));CHKERRQ(ierr);
    }
  } else {
    if (eps->ops->view) {
      ierr = (*eps->ops->view)(eps,viewer);CHKERRQ(ierr);
    }
  }
  ierr = PetscObjectTypeCompareAny((PetscObject)eps,&isexternal,EPSARPACK,EPSBLOPEX,EPSBLZPACK,EPSPRIMME,EPSTRLAN,"");CHKERRQ(ierr);
  if (!isexternal) {
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
    if (!eps->V) { ierr = EPSGetBV(eps,&eps->V);CHKERRQ(ierr); }
    ierr = BVView(eps->V,viewer);CHKERRQ(ierr);
    if (eps->rg) {
      ierr = RGIsTrivial(eps->rg,&istrivial);CHKERRQ(ierr);
      if (!istrivial) { ierr = RGView(eps->rg,viewer);CHKERRQ(ierr); }
    }
    if (eps->useds) {
      if (!eps->ds) { ierr = EPSGetDS(eps,&eps->ds);CHKERRQ(ierr); }
      ierr = DSView(eps->ds,viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  }
  if (!eps->st) { ierr = EPSGetST(eps,&eps->st);CHKERRQ(ierr); }
  ierr = STView(eps->st,viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   EPSViewFromOptions - View from options

   Collective on EPS

   Input Parameters:
+  eps  - the eigensolver context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: EPSView(), EPSCreate()
@*/
PetscErrorCode EPSViewFromOptions(EPS eps,PetscObject obj,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  ierr = PetscObjectViewFromOptions((PetscObject)eps,obj,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   EPSReasonView - Displays the reason an EPS solve converged or diverged.

   Collective on eps

   Input Parameters:
+  eps - the eigensolver context
-  viewer - the viewer to display the reason

   Options Database Keys:
.  -eps_converged_reason - print reason for convergence, and number of iterations

   Level: intermediate

.seealso: EPSSetConvergenceTest(), EPSSetTolerances(), EPSGetIterationNumber()
@*/
PetscErrorCode EPSReasonView(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isAscii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isAscii);CHKERRQ(ierr);
  if (isAscii) {
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)eps)->tablevel);CHKERRQ(ierr);
    if (eps->reason > 0) {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Linear eigensolve converged (%D eigenpair%s) due to %s; iterations %D\n",((PetscObject)eps)->prefix?((PetscObject)eps)->prefix:"",eps->nconv,(eps->nconv>1)?"s":"",EPSConvergedReasons[eps->reason],eps->its);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Linear eigensolve did not converge due to %s; iterations %D\n",((PetscObject)eps)->prefix?((PetscObject)eps)->prefix:"",EPSConvergedReasons[eps->reason],eps->its);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)eps)->tablevel);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   EPSReasonViewFromOptions - Processes command line options to determine if/how
   the EPS converged reason is to be viewed.

   Collective on eps

   Input Parameter:
.  eps - the eigensolver context

   Level: developer
@*/
PetscErrorCode EPSReasonViewFromOptions(EPS eps)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)eps),((PetscObject)eps)->options,((PetscObject)eps)->prefix,"-eps_converged_reason",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = EPSReasonView(eps,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSErrorView_ASCII(EPS eps,EPSErrorType etype,PetscViewer viewer)
{
  PetscReal      error;
  PetscInt       i,j,k,nvals;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  nvals = (eps->which==EPS_ALL)? eps->nconv: eps->nev;
  if (eps->which!=EPS_ALL && eps->nconv<nvals) {
    ierr = PetscViewerASCIIPrintf(viewer," Problem: less than %D eigenvalues converged\n\n",eps->nev);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  if (eps->which==EPS_ALL && !nvals) {
    ierr = PetscViewerASCIIPrintf(viewer," No eigenvalues have been found\n\n");CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  for (i=0;i<nvals;i++) {
    ierr = EPSComputeError(eps,i,etype,&error);CHKERRQ(ierr);
    if (error>=5.0*eps->tol) {
      ierr = PetscViewerASCIIPrintf(viewer," Problem: some of the first %D relative errors are higher than the tolerance\n\n",nvals);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }
  if (eps->which==EPS_ALL) {
    ierr = PetscViewerASCIIPrintf(viewer," Found %D eigenvalues, all of them computed up to the required tolerance:",nvals);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerASCIIPrintf(viewer," All requested eigenvalues computed up to the required tolerance:");CHKERRQ(ierr);
  }
  for (i=0;i<=(nvals-1)/8;i++) {
    ierr = PetscViewerASCIIPrintf(viewer,"\n     ");CHKERRQ(ierr);
    for (j=0;j<PetscMin(8,nvals-8*i);j++) {
      k = eps->perm[8*i+j];
      ierr = SlepcPrintEigenvalueASCII(viewer,eps->eigr[k],eps->eigi[k]);CHKERRQ(ierr);
      if (8*i+j+1<nvals) { ierr = PetscViewerASCIIPrintf(viewer,", ");CHKERRQ(ierr); }
    }
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSErrorView_DETAIL(EPS eps,EPSErrorType etype,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscReal      error,re,im;
  PetscScalar    kr,ki;
  PetscInt       i;
#define EXLEN 30
  char           ex[EXLEN],sep[]=" ---------------------- --------------------\n";

  PetscFunctionBegin;
  if (!eps->nconv) PetscFunctionReturn(0);
  switch (etype) {
    case EPS_ERROR_ABSOLUTE:
      ierr = PetscSNPrintf(ex,EXLEN,"   ||Ax-k%sx||",eps->isgeneralized?"B":"");CHKERRQ(ierr);
      break;
    case EPS_ERROR_RELATIVE:
      ierr = PetscSNPrintf(ex,EXLEN,"||Ax-k%sx||/||kx||",eps->isgeneralized?"B":"");CHKERRQ(ierr);
      break;
    case EPS_ERROR_BACKWARD:
      ierr = PetscSNPrintf(ex,EXLEN,"    eta(x,k)");CHKERRQ(ierr);
      break;
  }
  ierr = PetscViewerASCIIPrintf(viewer,"%s            k             %s\n%s",sep,ex,sep);CHKERRQ(ierr);
  for (i=0;i<eps->nconv;i++) {
    ierr = EPSGetEigenpair(eps,i,&kr,&ki,NULL,NULL);CHKERRQ(ierr);
    ierr = EPSComputeError(eps,i,etype,&error);CHKERRQ(ierr);
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

static PetscErrorCode EPSErrorView_MATLAB(EPS eps,EPSErrorType etype,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscReal      error;
  PetscInt       i;
  const char     *name;

  PetscFunctionBegin;
  ierr = PetscObjectGetName((PetscObject)eps,&name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Error_%s = [\n",name);CHKERRQ(ierr);
  for (i=0;i<eps->nconv;i++) {
    ierr = EPSComputeError(eps,i,etype,&error);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)error);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"];\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   EPSErrorView - Displays the errors associated with the computed solution
   (as well as the eigenvalues).

   Collective on eps

   Input Parameters:
+  eps    - the eigensolver context
.  etype  - error type
-  viewer - optional visualization context

   Options Database Key:
+  -eps_error_absolute - print absolute errors of each eigenpair
.  -eps_error_relative - print relative errors of each eigenpair
-  -eps_error_backward - print backward errors of each eigenpair

   Notes:
   By default, this function checks the error of all eigenpairs and prints
   the eigenvalues if all of them are below the requested tolerance.
   If the viewer has format=PETSC_VIEWER_ASCII_INFO_DETAIL then a table with
   eigenvalues and corresponding errors is printed.

   Level: intermediate

.seealso: EPSSolve(), EPSValuesView(), EPSVectorsView()
@*/
PetscErrorCode EPSErrorView(EPS eps,EPSErrorType etype,PetscViewer viewer)
{
  PetscBool         isascii;
  PetscViewerFormat format;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)eps),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(eps,1,viewer,2);
  EPSCheckSolved(eps,1);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (!isascii) PetscFunctionReturn(0);

  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  switch (format) {
    case PETSC_VIEWER_DEFAULT:
    case PETSC_VIEWER_ASCII_INFO:
      ierr = EPSErrorView_ASCII(eps,etype,viewer);CHKERRQ(ierr);
      break;
    case PETSC_VIEWER_ASCII_INFO_DETAIL:
      ierr = EPSErrorView_DETAIL(eps,etype,viewer);CHKERRQ(ierr);
      break;
    case PETSC_VIEWER_ASCII_MATLAB:
      ierr = EPSErrorView_MATLAB(eps,etype,viewer);CHKERRQ(ierr);
      break;
    default:
      ierr = PetscInfo1(eps,"Unsupported viewer format %s\n",PetscViewerFormats[format]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   EPSErrorViewFromOptions - Processes command line options to determine if/how
   the errors of the computed solution are to be viewed.

   Collective on eps

   Input Parameter:
.  eps - the eigensolver context

   Level: developer
@*/
PetscErrorCode EPSErrorViewFromOptions(EPS eps)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)eps),((PetscObject)eps)->options,((PetscObject)eps)->prefix,"-eps_error_absolute",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_ABSOLUTE,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)eps),((PetscObject)eps)->options,((PetscObject)eps)->prefix,"-eps_error_relative",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)eps),((PetscObject)eps)->options,((PetscObject)eps)->prefix,"-eps_error_backward",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_BACKWARD,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSValuesView_DRAW(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscDraw      draw;
  PetscDrawSP    drawsp;
  PetscReal      re,im;
  PetscInt       i,k;

  PetscFunctionBegin;
  if (!eps->nconv) PetscFunctionReturn(0);
  ierr = PetscViewerDrawGetDraw(viewer,0,&draw);CHKERRQ(ierr);
  ierr = PetscDrawSetTitle(draw,"Computed Eigenvalues");CHKERRQ(ierr);
  ierr = PetscDrawSPCreate(draw,1,&drawsp);CHKERRQ(ierr);
  for (i=0;i<eps->nconv;i++) {
    k = eps->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(eps->eigr[k]);
    im = PetscImaginaryPart(eps->eigr[k]);
#else
    re = eps->eigr[k];
    im = eps->eigi[k];
#endif
    ierr = PetscDrawSPAddPoint(drawsp,&re,&im);CHKERRQ(ierr);
  }
  ierr = PetscDrawSPDraw(drawsp,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscDrawSPSave(drawsp);CHKERRQ(ierr);
  ierr = PetscDrawSPDestroy(&drawsp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSValuesView_ASCII(EPS eps,PetscViewer viewer)
{
  PetscInt       i,k;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"Eigenvalues = \n");CHKERRQ(ierr);
  for (i=0;i<eps->nconv;i++) {
    k = eps->perm[i];
    ierr = PetscViewerASCIIPrintf(viewer,"   ");CHKERRQ(ierr);
    ierr = SlepcPrintEigenvalueASCII(viewer,eps->eigr[k],eps->eigi[k]);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSValuesView_MATLAB(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  PetscReal      re,im;
  const char     *name;

  PetscFunctionBegin;
  ierr = PetscObjectGetName((PetscObject)eps,&name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Lambda_%s = [\n",name);CHKERRQ(ierr);
  for (i=0;i<eps->nconv;i++) {
    k = eps->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(eps->eigr[k]);
    im = PetscImaginaryPart(eps->eigr[k]);
#else
    re = eps->eigr[k];
    im = eps->eigi[k];
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
   EPSValuesView - Displays the computed eigenvalues in a viewer.

   Collective on eps

   Input Parameters:
+  eps    - the eigensolver context
-  viewer - the viewer

   Options Database Key:
.  -eps_view_values - print computed eigenvalues

   Level: intermediate

.seealso: EPSSolve(), EPSVectorsView(), EPSErrorView()
@*/
PetscErrorCode EPSValuesView(EPS eps,PetscViewer viewer)
{
  PetscBool         isascii,isdraw;
  PetscViewerFormat format;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)eps),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(eps,1,viewer,2);
  EPSCheckSolved(eps,1);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isdraw) {
    ierr = EPSValuesView_DRAW(eps,viewer);CHKERRQ(ierr);
  } else if (isascii) {
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    switch (format) {
      case PETSC_VIEWER_DEFAULT:
      case PETSC_VIEWER_ASCII_INFO:
      case PETSC_VIEWER_ASCII_INFO_DETAIL:
        ierr = EPSValuesView_ASCII(eps,viewer);CHKERRQ(ierr);
        break;
      case PETSC_VIEWER_ASCII_MATLAB:
        ierr = EPSValuesView_MATLAB(eps,viewer);CHKERRQ(ierr);
        break;
      default:
        ierr = PetscInfo1(eps,"Unsupported viewer format %s\n",PetscViewerFormats[format]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@
   EPSValuesViewFromOptions - Processes command line options to determine if/how
   the computed eigenvalues are to be viewed.

   Collective on eps

   Input Parameters:
.  eps - the eigensolver context

   Level: developer
@*/
PetscErrorCode EPSValuesViewFromOptions(EPS eps)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)eps),((PetscObject)eps)->options,((PetscObject)eps)->prefix,"-eps_view_values",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = EPSValuesView(eps,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
   EPSVectorsView - Outputs computed eigenvectors to a viewer.

   Collective on eps

   Parameter:
+  eps    - the eigensolver context
-  viewer - the viewer

   Options Database Keys:
.  -eps_view_vectors - output eigenvectors.

   Notes:
   If PETSc was configured with real scalars, complex conjugate eigenvectors
   will be viewed as two separate real vectors, one containing the real part
   and another one containing the imaginary part.

   If left eigenvectors were computed with a two-sided eigensolver, the right
   and left eigenvectors are interleaved, that is, the vectors are output in
   the following order: X0, Y0, X1, Y1, X2, Y2, ...

   Level: intermediate

.seealso: EPSSolve(), EPSValuesView(), EPSErrorView()
@*/
PetscErrorCode EPSVectorsView(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  Vec            x;
#define NMLEN 30
  char           vname[NMLEN];
  const char     *ename;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)eps),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(eps,1,viewer,2);
  EPSCheckSolved(eps,1);
  if (eps->nconv) {
    ierr = PetscObjectGetName((PetscObject)eps,&ename);CHKERRQ(ierr);
    ierr = EPSComputeVectors(eps);CHKERRQ(ierr);
    for (i=0;i<eps->nconv;i++) {
      k = eps->perm[i];
      ierr = PetscSNPrintf(vname,NMLEN,"X%d_%s",(int)i,ename);CHKERRQ(ierr);
      ierr = BVGetColumn(eps->V,k,&x);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)x,vname);CHKERRQ(ierr);
      ierr = VecView(x,viewer);CHKERRQ(ierr);
      ierr = BVRestoreColumn(eps->V,k,&x);CHKERRQ(ierr);
      if (eps->twosided) {
        ierr = PetscSNPrintf(vname,NMLEN,"Y%d_%s",(int)i,ename);CHKERRQ(ierr);
        ierr = BVGetColumn(eps->W,k,&x);CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)x,vname);CHKERRQ(ierr);
        ierr = VecView(x,viewer);CHKERRQ(ierr);
        ierr = BVRestoreColumn(eps->W,k,&x);CHKERRQ(ierr);
      }
    }
  }
  PetscFunctionReturn(0);
}

/*@
   EPSVectorsViewFromOptions - Processes command line options to determine if/how
   the computed eigenvectors are to be viewed.

   Collective on eps

   Input Parameter:
.  eps - the eigensolver context

   Level: developer
@*/
PetscErrorCode EPSVectorsViewFromOptions(EPS eps)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg = PETSC_FALSE;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)eps),((PetscObject)eps)->options,((PetscObject)eps)->prefix,"-eps_view_vectors",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = EPSVectorsView(eps,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

