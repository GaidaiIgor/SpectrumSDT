/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   The PEP routines related to various viewers
*/

#include <slepc/private/pepimpl.h>      /*I "slepcpep.h" I*/
#include <petscdraw.h>

/*@C
   PEPView - Prints the PEP data structure.

   Collective on pep

   Input Parameters:
+  pep - the polynomial eigenproblem solver context
-  viewer - optional visualization context

   Options Database Key:
.  -pep_view -  Calls PEPView() at end of PEPSolve()

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
PetscErrorCode PEPView(PEP pep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  const char     *type=NULL;
  char           str[50];
  PetscBool      isascii,islinear,istrivial;
  PetscInt       i;
  PetscViewer    sviewer;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)pep),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(pep,1,viewer,2);

#if defined(PETSC_USE_COMPLEX)
#define HERM "hermitian"
#else
#define HERM "symmetric"
#endif
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscObjectPrintClassNamePrefixType((PetscObject)pep,viewer);CHKERRQ(ierr);
    if (pep->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*pep->ops->view)(pep,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    if (pep->problem_type) {
      switch (pep->problem_type) {
        case PEP_GENERAL:    type = "general polynomial eigenvalue problem"; break;
        case PEP_HERMITIAN:  type = HERM " polynomial eigenvalue problem"; break;
        case PEP_HYPERBOLIC: type = "hyperbolic polynomial eigenvalue problem"; break;
        case PEP_GYROSCOPIC: type = "gyroscopic polynomial eigenvalue problem"; break;
      }
    } else type = "not yet set";
    ierr = PetscViewerASCIIPrintf(viewer,"  problem type: %s\n",type);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  polynomial represented in %s basis\n",PEPBasisTypes[pep->basis]);CHKERRQ(ierr);
    switch (pep->scale) {
      case PEP_SCALE_NONE:
        break;
      case PEP_SCALE_SCALAR:
        ierr = PetscViewerASCIIPrintf(viewer,"  parameter scaling enabled, with scaling factor=%g\n",(double)pep->sfactor);CHKERRQ(ierr);
        break;
      case PEP_SCALE_DIAGONAL:
        ierr = PetscViewerASCIIPrintf(viewer,"  diagonal balancing enabled, with its=%D and lambda=%g\n",pep->sits,(double)pep->slambda);CHKERRQ(ierr);
        break;
      case PEP_SCALE_BOTH:
        ierr = PetscViewerASCIIPrintf(viewer,"  parameter scaling & diagonal balancing enabled, with scaling factor=%g, its=%D and lambda=%g\n",(double)pep->sfactor,pep->sits,(double)pep->slambda);CHKERRQ(ierr);
        break;
    }
    ierr = PetscViewerASCIIPrintf(viewer,"  selected portion of the spectrum: ");CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    ierr = SlepcSNPrintfScalar(str,50,pep->target,PETSC_FALSE);CHKERRQ(ierr);
    if (!pep->which) {
      ierr = PetscViewerASCIIPrintf(viewer,"not yet set\n");CHKERRQ(ierr);
    } else switch (pep->which) {
      case PEP_WHICH_USER:
        ierr = PetscViewerASCIIPrintf(viewer,"user defined\n");CHKERRQ(ierr);
        break;
      case PEP_TARGET_MAGNITUDE:
        ierr = PetscViewerASCIIPrintf(viewer,"closest to target: %s (in magnitude)\n",str);CHKERRQ(ierr);
        break;
      case PEP_TARGET_REAL:
        ierr = PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the real axis)\n",str);CHKERRQ(ierr);
        break;
      case PEP_TARGET_IMAGINARY:
        ierr = PetscViewerASCIIPrintf(viewer,"closest to target: %s (along the imaginary axis)\n",str);CHKERRQ(ierr);
        break;
      case PEP_LARGEST_MAGNITUDE:
        ierr = PetscViewerASCIIPrintf(viewer,"largest eigenvalues in magnitude\n");CHKERRQ(ierr);
        break;
      case PEP_SMALLEST_MAGNITUDE:
        ierr = PetscViewerASCIIPrintf(viewer,"smallest eigenvalues in magnitude\n");CHKERRQ(ierr);
        break;
      case PEP_LARGEST_REAL:
        ierr = PetscViewerASCIIPrintf(viewer,"largest real parts\n");CHKERRQ(ierr);
        break;
      case PEP_SMALLEST_REAL:
        ierr = PetscViewerASCIIPrintf(viewer,"smallest real parts\n");CHKERRQ(ierr);
        break;
      case PEP_LARGEST_IMAGINARY:
        ierr = PetscViewerASCIIPrintf(viewer,"largest imaginary parts\n");CHKERRQ(ierr);
        break;
      case PEP_SMALLEST_IMAGINARY:
        ierr = PetscViewerASCIIPrintf(viewer,"smallest imaginary parts\n");CHKERRQ(ierr);
        break;
      case PEP_ALL:
        if (pep->inta || pep->intb) {
          ierr = PetscViewerASCIIPrintf(viewer,"all eigenvalues in interval [%g,%g]\n",(double)pep->inta,(double)pep->intb);CHKERRQ(ierr);
        } else {
          ierr = PetscViewerASCIIPrintf(viewer,"all eigenvalues in the region\n");CHKERRQ(ierr);
        }
        break;
    }
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  number of eigenvalues (nev): %D\n",pep->nev);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  number of column vectors (ncv): %D\n",pep->ncv);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum dimension of projected problem (mpd): %D\n",pep->mpd);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  maximum number of iterations: %D\n",pep->max_it);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  tolerance: %g\n",(double)pep->tol);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  convergence test: ");CHKERRQ(ierr);
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
    switch (pep->conv) {
    case PEP_CONV_ABS:
      ierr = PetscViewerASCIIPrintf(viewer,"absolute\n");CHKERRQ(ierr);break;
    case PEP_CONV_REL:
      ierr = PetscViewerASCIIPrintf(viewer,"relative to the eigenvalue\n");CHKERRQ(ierr);break;
    case PEP_CONV_NORM:
      ierr = PetscViewerASCIIPrintf(viewer,"relative to the matrix norms\n");CHKERRQ(ierr);
      if (pep->nrma) {
        ierr = PetscViewerASCIIPrintf(viewer,"  computed matrix norms: %g",(double)pep->nrma[0]);CHKERRQ(ierr);
        for (i=1;i<pep->nmat;i++) {
          ierr = PetscViewerASCIIPrintf(viewer,", %g",(double)pep->nrma[i]);CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
      }
      break;
    case PEP_CONV_USER:
      ierr = PetscViewerASCIIPrintf(viewer,"user-defined\n");CHKERRQ(ierr);break;
    }
    ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  extraction type: %s\n",PEPExtractTypes[pep->extract]);CHKERRQ(ierr);
    if (pep->refine) {
      ierr = PetscViewerASCIIPrintf(viewer,"  iterative refinement: %s, with %s scheme\n",PEPRefineTypes[pep->refine],PEPRefineSchemes[pep->scheme]);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"  refinement stopping criterion: tol=%g, its=%D\n",(double)pep->rtol,pep->rits);CHKERRQ(ierr);
      if (pep->npart>1) {
        ierr = PetscViewerASCIIPrintf(viewer,"  splitting communicator in %D partitions for refinement\n",pep->npart);CHKERRQ(ierr);
      }
    }
    if (pep->nini) {
      ierr = PetscViewerASCIIPrintf(viewer,"  dimension of user-provided initial space: %D\n",PetscAbs(pep->nini));CHKERRQ(ierr);
    }
  } else {
    if (pep->ops->view) {
      ierr = (*pep->ops->view)(pep,viewer);CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
  if (!pep->V) { ierr = PEPGetBV(pep,&pep->V);CHKERRQ(ierr); }
  ierr = BVView(pep->V,viewer);CHKERRQ(ierr);
  if (!pep->rg) { ierr = PEPGetRG(pep,&pep->rg);CHKERRQ(ierr); }
  ierr = RGIsTrivial(pep->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) { ierr = RGView(pep->rg,viewer);CHKERRQ(ierr); }
  ierr = PetscObjectTypeCompare((PetscObject)pep,PEPLINEAR,&islinear);CHKERRQ(ierr);
  if (!islinear) {
    if (!pep->ds) { ierr = PEPGetDS(pep,&pep->ds);CHKERRQ(ierr); }
    ierr = DSView(pep->ds,viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  if (!pep->st) { ierr = PEPGetST(pep,&pep->st);CHKERRQ(ierr); }
  ierr = STView(pep->st,viewer);CHKERRQ(ierr);
  if (pep->refine!=PEP_REFINE_NONE) {
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    if (pep->npart>1) {
      if (pep->refinesubc->color==0) {
        ierr = PetscViewerASCIIGetStdout(PetscSubcommChild(pep->refinesubc),&sviewer);CHKERRQ(ierr);
        ierr = KSPView(pep->refineksp,sviewer);CHKERRQ(ierr);
      }
    } else {
      ierr = KSPView(pep->refineksp,viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   PEPViewFromOptions - View from options

   Collective on PEP

   Input Parameters:
+  pep  - the eigensolver context
.  obj  - optional object
-  name - command line option

   Level: intermediate

.seealso: PEPView(), PEPCreate()
@*/
PetscErrorCode PEPViewFromOptions(PEP pep,PetscObject obj,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  ierr = PetscObjectViewFromOptions((PetscObject)pep,obj,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PEPReasonView - Displays the reason a PEP solve converged or diverged.

   Collective on pep

   Input Parameters:
+  pep - the eigensolver context
-  viewer - the viewer to display the reason

   Options Database Keys:
.  -pep_converged_reason - print reason for convergence, and number of iterations

   Level: intermediate

.seealso: PEPSetConvergenceTest(), PEPSetTolerances(), PEPGetIterationNumber()
@*/
PetscErrorCode PEPReasonView(PEP pep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isAscii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isAscii);CHKERRQ(ierr);
  if (isAscii) {
    ierr = PetscViewerASCIIAddTab(viewer,((PetscObject)pep)->tablevel);CHKERRQ(ierr);
    if (pep->reason > 0) {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Polynomial eigensolve converged (%D eigenpair%s) due to %s; iterations %D\n",((PetscObject)pep)->prefix?((PetscObject)pep)->prefix:"",pep->nconv,(pep->nconv>1)?"s":"",PEPConvergedReasons[pep->reason],pep->its);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Polynomial eigensolve did not converge due to %s; iterations %D\n",((PetscObject)pep)->prefix?((PetscObject)pep)->prefix:"",PEPConvergedReasons[pep->reason],pep->its);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIISubtractTab(viewer,((PetscObject)pep)->tablevel);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PEPReasonViewFromOptions - Processes command line options to determine if/how
   the PEP converged reason is to be viewed.

   Collective on pep

   Input Parameter:
.  pep - the eigensolver context

   Level: developer
@*/
PetscErrorCode PEPReasonViewFromOptions(PEP pep)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_converged_reason",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = PEPReasonView(pep,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPErrorView_ASCII(PEP pep,PEPErrorType etype,PetscViewer viewer)
{
  PetscReal      error;
  PetscInt       i,j,k;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (pep->nconv<pep->nev) {
    ierr = PetscViewerASCIIPrintf(viewer," Problem: less than %D eigenvalues converged\n\n",pep->nev);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  for (i=0;i<pep->nev;i++) {
    ierr = PEPComputeError(pep,i,etype,&error);CHKERRQ(ierr);
    if (error>=5.0*pep->tol) {
      ierr = PetscViewerASCIIPrintf(viewer," Problem: some of the first %D relative errors are higher than the tolerance\n\n",pep->nev);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }
  ierr = PetscViewerASCIIPrintf(viewer," All requested eigenvalues computed up to the required tolerance:");CHKERRQ(ierr);
  for (i=0;i<=(pep->nev-1)/8;i++) {
    ierr = PetscViewerASCIIPrintf(viewer,"\n     ");CHKERRQ(ierr);
    for (j=0;j<PetscMin(8,pep->nev-8*i);j++) {
      k = pep->perm[8*i+j];
      ierr = SlepcPrintEigenvalueASCII(viewer,pep->eigr[k],pep->eigi[k]);CHKERRQ(ierr);
      if (8*i+j+1<pep->nev) { ierr = PetscViewerASCIIPrintf(viewer,", ");CHKERRQ(ierr); }
    }
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPErrorView_DETAIL(PEP pep,PEPErrorType etype,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscReal      error,re,im;
  PetscScalar    kr,ki;
  PetscInt       i;
#define EXLEN 30
  char           ex[EXLEN],sep[]=" ---------------------- --------------------\n";

  PetscFunctionBegin;
  if (!pep->nconv) PetscFunctionReturn(0);
  switch (etype) {
    case PEP_ERROR_ABSOLUTE:
      ierr = PetscSNPrintf(ex,EXLEN,"   ||P(k)x||");CHKERRQ(ierr);
      break;
    case PEP_ERROR_RELATIVE:
      ierr = PetscSNPrintf(ex,EXLEN,"||P(k)x||/||kx||");CHKERRQ(ierr);
      break;
    case PEP_ERROR_BACKWARD:
      ierr = PetscSNPrintf(ex,EXLEN,"    eta(x,k)");CHKERRQ(ierr);
      break;
  }
  ierr = PetscViewerASCIIPrintf(viewer,"%s            k             %s\n%s",sep,ex,sep);CHKERRQ(ierr);
  for (i=0;i<pep->nconv;i++) {
    ierr = PEPGetEigenpair(pep,i,&kr,&ki,NULL,NULL);CHKERRQ(ierr);
    ierr = PEPComputeError(pep,i,etype,&error);CHKERRQ(ierr);
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

static PetscErrorCode PEPErrorView_MATLAB(PEP pep,PEPErrorType etype,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscReal      error;
  PetscInt       i;
  const char     *name;

  PetscFunctionBegin;
  ierr = PetscObjectGetName((PetscObject)pep,&name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Error_%s = [\n",name);CHKERRQ(ierr);
  for (i=0;i<pep->nconv;i++) {
    ierr = PEPComputeError(pep,i,etype,&error);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"%18.16e\n",(double)error);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"];\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PEPErrorView - Displays the errors associated with the computed solution
   (as well as the eigenvalues).

   Collective on pep

   Input Parameters:
+  pep    - the eigensolver context
.  etype  - error type
-  viewer - optional visualization context

   Options Database Key:
+  -pep_error_absolute - print absolute errors of each eigenpair
.  -pep_error_relative - print relative errors of each eigenpair
-  -pep_error_backward - print backward errors of each eigenpair

   Notes:
   By default, this function checks the error of all eigenpairs and prints
   the eigenvalues if all of them are below the requested tolerance.
   If the viewer has format=PETSC_VIEWER_ASCII_INFO_DETAIL then a table with
   eigenvalues and corresponding errors is printed.

   Level: intermediate

.seealso: PEPSolve(), PEPValuesView(), PEPVectorsView()
@*/
PetscErrorCode PEPErrorView(PEP pep,PEPErrorType etype,PetscViewer viewer)
{
  PetscBool         isascii;
  PetscViewerFormat format;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)pep),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(pep,1,viewer,2);
  PEPCheckSolved(pep,1);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (!isascii) PetscFunctionReturn(0);

  ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
  switch (format) {
    case PETSC_VIEWER_DEFAULT:
    case PETSC_VIEWER_ASCII_INFO:
      ierr = PEPErrorView_ASCII(pep,etype,viewer);CHKERRQ(ierr);
      break;
    case PETSC_VIEWER_ASCII_INFO_DETAIL:
      ierr = PEPErrorView_DETAIL(pep,etype,viewer);CHKERRQ(ierr);
      break;
    case PETSC_VIEWER_ASCII_MATLAB:
      ierr = PEPErrorView_MATLAB(pep,etype,viewer);CHKERRQ(ierr);
      break;
    default:
      ierr = PetscInfo1(pep,"Unsupported viewer format %s\n",PetscViewerFormats[format]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PEPErrorViewFromOptions - Processes command line options to determine if/how
   the errors of the computed solution are to be viewed.

   Collective on pep

   Input Parameter:
.  pep - the eigensolver context

   Level: developer
@*/
PetscErrorCode PEPErrorViewFromOptions(PEP pep)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_error_absolute",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = PEPErrorView(pep,PEP_ERROR_ABSOLUTE,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_error_relative",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = PEPErrorView(pep,PEP_ERROR_RELATIVE,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_error_backward",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPValuesView_DRAW(PEP pep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscDraw      draw;
  PetscDrawSP    drawsp;
  PetscReal      re,im;
  PetscInt       i,k;

  PetscFunctionBegin;
  if (!pep->nconv) PetscFunctionReturn(0);
  ierr = PetscViewerDrawGetDraw(viewer,0,&draw);CHKERRQ(ierr);
  ierr = PetscDrawSetTitle(draw,"Computed Eigenvalues");CHKERRQ(ierr);
  ierr = PetscDrawSPCreate(draw,1,&drawsp);CHKERRQ(ierr);
  for (i=0;i<pep->nconv;i++) {
    k = pep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(pep->eigr[k]);
    im = PetscImaginaryPart(pep->eigr[k]);
#else
    re = pep->eigr[k];
    im = pep->eigi[k];
#endif
    ierr = PetscDrawSPAddPoint(drawsp,&re,&im);CHKERRQ(ierr);
  }
  ierr = PetscDrawSPDraw(drawsp,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscDrawSPSave(drawsp);CHKERRQ(ierr);
  ierr = PetscDrawSPDestroy(&drawsp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPValuesView_ASCII(PEP pep,PetscViewer viewer)
{
  PetscInt       i,k;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(viewer,"Eigenvalues = \n");CHKERRQ(ierr);
  for (i=0;i<pep->nconv;i++) {
    k = pep->perm[i];
    ierr = PetscViewerASCIIPrintf(viewer,"   ");CHKERRQ(ierr);
    ierr = SlepcPrintEigenvalueASCII(viewer,pep->eigr[k],pep->eigi[k]);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPValuesView_MATLAB(PEP pep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  PetscReal      re,im;
  const char     *name;

  PetscFunctionBegin;
  ierr = PetscObjectGetName((PetscObject)pep,&name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Lambda_%s = [\n",name);CHKERRQ(ierr);
  for (i=0;i<pep->nconv;i++) {
    k = pep->perm[i];
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(pep->eigr[k]);
    im = PetscImaginaryPart(pep->eigr[k]);
#else
    re = pep->eigr[k];
    im = pep->eigi[k];
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
   PEPValuesView - Displays the computed eigenvalues in a viewer.

   Collective on pep

   Input Parameters:
+  pep    - the eigensolver context
-  viewer - the viewer

   Options Database Key:
.  -pep_view_values - print computed eigenvalues

   Level: intermediate

.seealso: PEPSolve(), PEPVectorsView(), PEPErrorView()
@*/
PetscErrorCode PEPValuesView(PEP pep,PetscViewer viewer)
{
  PetscBool         isascii,isdraw;
  PetscViewerFormat format;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)pep),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(pep,1,viewer,2);
  PEPCheckSolved(pep,1);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isdraw) {
    ierr = PEPValuesView_DRAW(pep,viewer);CHKERRQ(ierr);
  } else if (isascii) {
    ierr = PetscViewerGetFormat(viewer,&format);CHKERRQ(ierr);
    switch (format) {
      case PETSC_VIEWER_DEFAULT:
      case PETSC_VIEWER_ASCII_INFO:
      case PETSC_VIEWER_ASCII_INFO_DETAIL:
        ierr = PEPValuesView_ASCII(pep,viewer);CHKERRQ(ierr);
        break;
      case PETSC_VIEWER_ASCII_MATLAB:
        ierr = PEPValuesView_MATLAB(pep,viewer);CHKERRQ(ierr);
        break;
      default:
        ierr = PetscInfo1(pep,"Unsupported viewer format %s\n",PetscViewerFormats[format]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@
   PEPValuesViewFromOptions - Processes command line options to determine if/how
   the computed eigenvalues are to be viewed.

   Collective on pep

   Input Parameter:
.  pep - the eigensolver context

   Level: developer
@*/
PetscErrorCode PEPValuesViewFromOptions(PEP pep)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_view_values",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = PEPValuesView(pep,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@C
   PEPVectorsView - Outputs computed eigenvectors to a viewer.

   Collective on pep

   Input Parameters:
+  pep    - the eigensolver context
-  viewer - the viewer

   Options Database Keys:
.  -pep_view_vectors - output eigenvectors.

   Note:
   If PETSc was configured with real scalars, complex conjugate eigenvectors
   will be viewed as two separate real vectors, one containing the real part
   and another one containing the imaginary part.

   Level: intermediate

.seealso: PEPSolve(), PEPValuesView(), PEPErrorView()
@*/
PetscErrorCode PEPVectorsView(PEP pep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  Vec            x;
#define NMLEN 30
  char           vname[NMLEN];
  const char     *ename;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pep,PEP_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)pep),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(pep,1,viewer,2);
  PEPCheckSolved(pep,1);
  if (pep->nconv) {
    ierr = PetscObjectGetName((PetscObject)pep,&ename);CHKERRQ(ierr);
    ierr = PEPComputeVectors(pep);CHKERRQ(ierr);
    for (i=0;i<pep->nconv;i++) {
      k = pep->perm[i];
      ierr = PetscSNPrintf(vname,NMLEN,"V%d_%s",(int)i,ename);CHKERRQ(ierr);
      ierr = BVGetColumn(pep->V,k,&x);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)x,vname);CHKERRQ(ierr);
      ierr = VecView(x,viewer);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pep->V,k,&x);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@
   PEPVectorsViewFromOptions - Processes command line options to determine if/how
   the computed eigenvectors are to be viewed.

   Collective on pep

   Input Parameter:
.  pep - the eigensolver context

   Level: developer
@*/
PetscErrorCode PEPVectorsViewFromOptions(PEP pep)
{
  PetscErrorCode    ierr;
  PetscViewer       viewer;
  PetscBool         flg = PETSC_FALSE;
  static PetscBool  incall = PETSC_FALSE;
  PetscViewerFormat format;

  PetscFunctionBegin;
  if (incall) PetscFunctionReturn(0);
  incall = PETSC_TRUE;
  ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)pep),((PetscObject)pep)->options,((PetscObject)pep)->prefix,"-pep_view_vectors",&viewer,&format,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerPushFormat(viewer,format);CHKERRQ(ierr);
    ierr = PEPVectorsView(pep,viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  incall = PETSC_FALSE;
  PetscFunctionReturn(0);
}

