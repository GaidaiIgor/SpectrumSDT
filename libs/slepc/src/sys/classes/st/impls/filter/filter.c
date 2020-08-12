/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Filter spectral transformation, to encapsulate polynomial filters
*/

#include <slepc/private/stimpl.h>         /*I "slepcst.h" I*/
#include "filter.h"

/*
   Operator (filter):
               Op               P         M
   if nmat=1:  p(A)             NULL      p(A)
*/
PetscErrorCode STComputeOperator_Filter(ST st)
{
  PetscErrorCode ierr;
  ST_FILTER      *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (st->nmat>1) SETERRQ(PetscObjectComm((PetscObject)st),1,"Only implemented for standard eigenvalue problem");
  if (ctx->intb >= PETSC_MAX_REAL && ctx->inta <= PETSC_MIN_REAL) SETERRQ(PetscObjectComm((PetscObject)st),1,"Must pass an interval with STFilterSetInterval()");
  if (ctx->right == 0.0 && ctx->left == 0.0) SETERRQ(PetscObjectComm((PetscObject)st),1,"Must pass an approximate numerical range with STFilterSetRange()");
  if (ctx->left > ctx->inta || ctx->right < ctx->intb) SETERRQ4(PetscObjectComm((PetscObject)st),1,"The requested interval [%g,%g] must be contained in the numerical range [%g,%g]",(double)ctx->inta,(double)ctx->intb,(double)ctx->left,(double)ctx->right);
  if (!ctx->polyDegree) ctx->polyDegree = 100;
  ctx->frame[0] = ctx->left;
  ctx->frame[1] = ctx->inta;
  ctx->frame[2] = ctx->intb;
  ctx->frame[3] = ctx->right;
  ierr = STFilter_FILTLAN_setFilter(st,&st->T[0]);CHKERRQ(ierr);
  st->M = st->T[0];
  ierr = MatDestroy(&st->P);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STSetUp_Filter(ST st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = STSetWorkVecs(st,4);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STSetFromOptions_Filter(PetscOptionItems *PetscOptionsObject,ST st)
{
  PetscErrorCode ierr;
  PetscReal      array[2]={0,0};
  PetscInt       k;
  PetscBool      flg;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"ST Filter Options");CHKERRQ(ierr);

    k = 2;
    ierr = PetscOptionsRealArray("-st_filter_interval","Interval containing the desired eigenvalues (two real values separated with a comma without spaces)","STFilterSetInterval",array,&k,&flg);CHKERRQ(ierr);
    if (flg) {
      if (k<2) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_SIZ,"Must pass two values in -st_filter_interval (comma-separated without spaces)");
      ierr = STFilterSetInterval(st,array[0],array[1]);CHKERRQ(ierr);
    }
    k = 2;
    ierr = PetscOptionsRealArray("-st_filter_range","Interval containing all eigenvalues (two real values separated with a comma without spaces)","STFilterSetRange",array,&k,&flg);CHKERRQ(ierr);
    if (flg) {
      if (k<2) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_SIZ,"Must pass two values in -st_filter_range (comma-separated without spaces)");
      ierr = STFilterSetRange(st,array[0],array[1]);CHKERRQ(ierr);
    }
    ierr = PetscOptionsInt("-st_filter_degree","Degree of filter polynomial","STFilterSetDegree",100,&k,&flg);CHKERRQ(ierr);
    if (flg) { ierr = STFilterSetDegree(st,k);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STFilterSetInterval_Filter(ST st,PetscReal inta,PetscReal intb)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (inta>intb) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be inta<intb");
  if (ctx->inta != inta || ctx->intb != intb) {
    ctx->inta   = inta;
    ctx->intb   = intb;
    st->state   = ST_STATE_INITIAL;
    st->opready = PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}

/*@
   STFilterSetInterval - Defines the interval containing the desired eigenvalues.

   Logically Collective on st

   Input Parameters:
+  st   - the spectral transformation context
.  inta - left end of the interval
-  intb - right end of the interval

   Options Database Key:
.  -st_filter_interval <a,b> - set [a,b] as the interval of interest

   Level: intermediate

   Notes:
   The filter will be configured to emphasize eigenvalues contained in the given
   interval, and damp out eigenvalues outside it. If the interval is open, then
   the filter is low- or high-pass, otherwise it is mid-pass.

   Common usage is to set the interval in EPS with EPSSetInterval().

   The interval must be contained within the numerical range of the matrix, see
   STFilterSetRange().

.seealso: STFilterGetInterval(), STFilterSetRange(), EPSSetInterval()
@*/
PetscErrorCode STFilterSetInterval(ST st,PetscReal inta,PetscReal intb)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveReal(st,inta,2);
  PetscValidLogicalCollectiveReal(st,intb,3);
  ierr = PetscTryMethod(st,"STFilterSetInterval_C",(ST,PetscReal,PetscReal),(st,inta,intb));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STFilterGetInterval_Filter(ST st,PetscReal *inta,PetscReal *intb)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (inta) *inta = ctx->inta;
  if (intb) *intb = ctx->intb;
  PetscFunctionReturn(0);
}

/*@
   STFilterGetInterval - Gets the interval containing the desired eigenvalues.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
+  inta - left end of the interval
-  intb - right end of the interval

   Level: intermediate

.seealso: STFilterSetInterval()
@*/
PetscErrorCode STFilterGetInterval(ST st,PetscReal *inta,PetscReal *intb)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  ierr = PetscUseMethod(st,"STFilterGetInterval_C",(ST,PetscReal*,PetscReal*),(st,inta,intb));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STFilterSetRange_Filter(ST st,PetscReal left,PetscReal right)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (left>right) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_WRONG,"Badly defined interval, must be left<right");
  if (ctx->left != left || ctx->right != right) {
    ctx->left   = left;
    ctx->right  = right;
    st->state   = ST_STATE_INITIAL;
    st->opready = PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}

/*@
   STFilterSetRange - Defines the numerical range (or field of values) of the matrix, that is,
   the interval containing all eigenvalues.

   Logically Collective on st

   Input Parameters:
+  st    - the spectral transformation context
.  left  - left end of the interval
-  right - right end of the interval

   Options Database Key:
.  -st_filter_range <a,b> - set [a,b] as the numerical range

   Level: intermediate

   Notes:
   The filter will be most effective if the numerical range is tight, that is, left and right
   are good approximations to the leftmost and rightmost eigenvalues, respectively.

.seealso: STFilterGetRange(), STFilterSetInterval()
@*/
PetscErrorCode STFilterSetRange(ST st,PetscReal left,PetscReal right)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveReal(st,left,2);
  PetscValidLogicalCollectiveReal(st,right,3);
  ierr = PetscTryMethod(st,"STFilterSetRange_C",(ST,PetscReal,PetscReal),(st,left,right));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STFilterGetRange_Filter(ST st,PetscReal *left,PetscReal *right)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (left)  *left  = ctx->left;
  if (right) *right = ctx->right;
  PetscFunctionReturn(0);
}

/*@
   STFilterGetRange - Gets the interval containing all eigenvalues.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
+  left  - left end of the interval
-  right - right end of the interval

   Level: intermediate

.seealso: STFilterSetRange()
@*/
PetscErrorCode STFilterGetRange(ST st,PetscReal *left,PetscReal *right)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  ierr = PetscUseMethod(st,"STFilterGetRange_C",(ST,PetscReal*,PetscReal*),(st,left,right));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STFilterSetDegree_Filter(ST st,PetscInt deg)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  if (deg == PETSC_DEFAULT || deg == PETSC_DECIDE) {
    ctx->polyDegree = 0;
    st->state       = ST_STATE_INITIAL;
    st->opready     = PETSC_FALSE;
  } else {
    if (deg<=0) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of degree. Must be > 0");
    if (ctx->polyDegree != deg) {
      ctx->polyDegree = deg;
      st->state       = ST_STATE_INITIAL;
      st->opready     = PETSC_FALSE;
    }
  }
  PetscFunctionReturn(0);
}

/*@
   STFilterSetDegree - Sets the degree of the filter polynomial.

   Logically Collective on st

   Input Parameters:
+  st  - the spectral transformation context
-  deg - polynomial degree

   Options Database Key:
.  -st_filter_degree <deg> - sets the degree of the filter polynomial

   Level: intermediate

.seealso: STFilterGetDegree()
@*/
PetscErrorCode STFilterSetDegree(ST st,PetscInt deg)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidLogicalCollectiveInt(st,deg,2);
  ierr = PetscTryMethod(st,"STFilterSetDegree_C",(ST,PetscInt),(st,deg));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STFilterGetDegree_Filter(ST st,PetscInt *deg)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  *deg = ctx->polyDegree;
  PetscFunctionReturn(0);
}

/*@
   STFilterGetDegree - Gets the degree of the filter polynomial.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
.  deg - polynomial degree

   Level: intermediate

.seealso: STFilterSetDegree()
@*/
PetscErrorCode STFilterGetDegree(ST st,PetscInt *deg)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidIntPointer(deg,2);
  ierr = PetscUseMethod(st,"STFilterGetDegree_C",(ST,PetscInt*),(st,deg));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STFilterGetThreshold_Filter(ST st,PetscReal *gamma)
{
  ST_FILTER *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  *gamma = ctx->filterInfo->yLimit;
  PetscFunctionReturn(0);
}

/*@
   STFilterGetThreshold - Gets the threshold value gamma such that rho(lambda)>=gamma for
   eigenvalues lambda inside the wanted interval and rho(lambda)<gamma for those outside.

   Not Collective

   Input Parameter:
.  st  - the spectral transformation context

   Output Parameter:
.  gamma - the threshold value

   Level: developer
@*/
PetscErrorCode STFilterGetThreshold(ST st,PetscReal *gamma)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidRealPointer(gamma,2);
  ierr = PetscUseMethod(st,"STFilterGetThreshold_C",(ST,PetscReal*),(st,gamma));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STReset_Filter(ST st)
{
  PetscErrorCode ierr;
  ST_FILTER      *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  ctx->left  = 0.0;
  ctx->right = 0.0;
  ierr = MatDestroy(&ctx->T);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode STView_Filter(ST st,PetscViewer viewer)
{
  PetscErrorCode ierr;
  ST_FILTER      *ctx = (ST_FILTER*)st->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  Filter: interval of desired eigenvalues is [%g,%g]\n",(double)ctx->inta,(double)ctx->intb);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  Filter: numerical range is [%g,%g]\n",(double)ctx->left,(double)ctx->right);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  Filter: degree of filter polynomial is %D\n",ctx->polyDegree);CHKERRQ(ierr);
    if (st->state>=ST_STATE_SETUP) {
      ierr = PetscViewerASCIIPrintf(viewer,"  Filter: limit to accept eigenvalues: theta=%g\n",ctx->filterInfo->yLimit);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode STDestroy_Filter(ST st)
{
  PetscErrorCode ierr;
  ST_FILTER      *ctx = (ST_FILTER*)st->data;

  PetscFunctionBegin;
  ierr = PetscFree(ctx->opts);CHKERRQ(ierr);
  ierr = PetscFree(ctx->filterInfo);CHKERRQ(ierr);
  ierr = PetscFree(ctx->baseFilter);CHKERRQ(ierr);
  ierr = PetscFree(st->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterSetInterval_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterGetInterval_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterSetRange_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterGetRange_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterSetDegree_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterGetDegree_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterGetThreshold_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode STCreate_Filter(ST st)
{
  PetscErrorCode ierr;
  ST_FILTER      *ctx;
  FILTLAN_IOP    iop;
  FILTLAN_PFI    pfi;
  PetscFunctionBegin;
  ierr = PetscNewLog(st,&ctx);CHKERRQ(ierr);
  st->data = (void*)ctx;

  st->usesksp = PETSC_FALSE;

  ctx->inta               = PETSC_MIN_REAL;
  ctx->intb               = PETSC_MAX_REAL;
  ctx->left               = 0.0;
  ctx->right              = 0.0;
  ctx->polyDegree         = 0;
  ctx->baseDegree         = 10;

  ierr = PetscNewLog(st,&iop);CHKERRQ(ierr);
  ctx->opts               = iop;
  iop->intervalWeights[0] = 100.0;
  iop->intervalWeights[1] = 1.0;
  iop->intervalWeights[2] = 1.0;
  iop->intervalWeights[3] = 1.0;
  iop->intervalWeights[4] = 100.0;
  iop->transIntervalRatio = 0.6;
  iop->reverseInterval    = PETSC_FALSE;
  iop->initialPlateau     = 0.1;
  iop->plateauShrinkRate  = 1.5;
  iop->initialShiftStep   = 0.01;
  iop->shiftStepExpanRate = 1.5;
  iop->maxInnerIter       = 30;
  iop->yLimitTol          = 0.0001;
  iop->maxOuterIter       = 50;
  iop->numGridPoints      = 200;
  iop->yBottomLine        = 0.001;
  iop->yRippleLimit       = 0.8;

  ierr = PetscNewLog(st,&pfi);CHKERRQ(ierr);
  ctx->filterInfo         = pfi;

  st->ops->apply           = STApply_Generic;
  st->ops->setup           = STSetUp_Filter;
  st->ops->computeoperator = STComputeOperator_Filter;
  st->ops->setfromoptions  = STSetFromOptions_Filter;
  st->ops->destroy         = STDestroy_Filter;
  st->ops->reset           = STReset_Filter;
  st->ops->view            = STView_Filter;

  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterSetInterval_C",STFilterSetInterval_Filter);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterGetInterval_C",STFilterGetInterval_Filter);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterSetRange_C",STFilterSetRange_Filter);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterGetRange_C",STFilterGetRange_Filter);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterSetDegree_C",STFilterSetDegree_Filter);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterGetDegree_C",STFilterGetDegree_Filter);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STFilterGetThreshold_C",STFilterGetThreshold_Filter);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

