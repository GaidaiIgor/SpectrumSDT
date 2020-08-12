/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc nonlinear eigensolver: "interpol"

   Method: Polynomial interpolation

   Algorithm:

       Uses a PEP object to solve the interpolated NEP. Currently supports
       only Chebyshev interpolation on an interval.

   References:

       [1] C. Effenberger and D. Kresser, "Chebyshev interpolation for
           nonlinear eigenvalue problems", BIT 52:933-951, 2012.
*/

#include <slepc/private/nepimpl.h>         /*I "slepcnep.h" I*/
#include <slepc/private/pepimpl.h>         /*I "slepcpep.h" I*/

typedef struct {
  PEP       pep;
  PetscReal tol;       /* tolerance for norm of polynomial coefficients */
  PetscInt  maxdeg;    /* maximum degree of interpolation polynomial */
  PetscInt  deg;       /* actual degree of interpolation polynomial */
} NEP_INTERPOL;

PetscErrorCode NEPSetUp_Interpol(NEP nep)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;
  ST             st;
  RG             rg;
  PetscReal      a,b,c,d,s,tol;
  PetscScalar    zero=0.0;
  PetscBool      flg,istrivial,trackall;
  PetscInt       its,in;

  PetscFunctionBegin;
  ierr = NEPSetDimensions_Default(nep,nep->nev,&nep->ncv,&nep->mpd);CHKERRQ(ierr);
  if (nep->ncv>nep->nev+nep->mpd) SETERRQ(PetscObjectComm((PetscObject)nep),1,"The value of ncv must not be larger than nev+mpd");
  if (nep->max_it==PETSC_DEFAULT) nep->max_it = PetscMax(5000,2*nep->n/nep->ncv);
  if (nep->fui!=NEP_USER_INTERFACE_SPLIT) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"NEPINTERPOL only available for split operator");
  if (nep->stopping!=NEPStoppingBasic) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"This solver does not support user-defined stopping test");

  /* transfer PEP options */
  if (!ctx->pep) { ierr = NEPInterpolGetPEP(nep,&ctx->pep);CHKERRQ(ierr); }
  ierr = PEPSetBasis(ctx->pep,PEP_BASIS_CHEBYSHEV1);CHKERRQ(ierr);
  ierr = PEPSetWhichEigenpairs(ctx->pep,PEP_TARGET_MAGNITUDE);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)ctx->pep,PEPJD,&flg);CHKERRQ(ierr);
  if (!flg) {
    ierr = PEPGetST(ctx->pep,&st);CHKERRQ(ierr);
    ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
  }
  ierr = PEPSetDimensions(ctx->pep,nep->nev,nep->ncv,nep->mpd);CHKERRQ(ierr);
  ierr = PEPGetTolerances(ctx->pep,&tol,&its);CHKERRQ(ierr);
  if (tol==PETSC_DEFAULT) tol = (nep->tol==PETSC_DEFAULT)?SLEPC_DEFAULT_TOL:nep->tol;
  if (ctx->tol==PETSC_DEFAULT) ctx->tol = tol;
  if (its==PETSC_DEFAULT) its = nep->max_it;
  ierr = PEPSetTolerances(ctx->pep,tol,its);CHKERRQ(ierr);
  ierr = NEPGetTrackAll(nep,&trackall);CHKERRQ(ierr);
  ierr = PEPSetTrackAll(ctx->pep,trackall);CHKERRQ(ierr);

  /* transfer region options */
  ierr = RGIsTrivial(nep->rg,&istrivial);CHKERRQ(ierr);
  if (istrivial) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"NEPINTERPOL requires a nontrivial region");
  ierr = PetscObjectTypeCompare((PetscObject)nep->rg,RGINTERVAL,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Only implemented for interval regions");
  ierr = RGIntervalGetEndpoints(nep->rg,&a,&b,&c,&d);CHKERRQ(ierr);
  if (a<=-PETSC_MAX_REAL || b>=PETSC_MAX_REAL) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Only implemented for bounded intervals");
  ierr = PEPGetRG(ctx->pep,&rg);CHKERRQ(ierr);
  ierr = RGSetType(rg,RGINTERVAL);CHKERRQ(ierr);
  if (a==b) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"Only implemented for intervals on the real axis");
  s = 2.0/(b-a);
  c = c*s;
  d = d*s;
  ierr = RGIntervalSetEndpoints(rg,-1.0,1.0,c,d);CHKERRQ(ierr);
  ierr = RGCheckInside(nep->rg,1,&nep->target,&zero,&in);CHKERRQ(ierr);
  if (in<0) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_SUP,"The target is not inside the target set");
  ierr = PEPSetTarget(ctx->pep,(nep->target-(a+b)/2)*s);CHKERRQ(ierr);

  ierr = NEPAllocateSolution(nep,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Input:
    d, number of nodes to compute
    a,b, interval extrems
  Output:
    *x, array containing the d Chebyshev nodes of the interval [a,b]
    *dct2, coefficients to compute a discrete cosine transformation (DCT-II)
*/
static PetscErrorCode ChebyshevNodes(PetscInt d,PetscReal a,PetscReal b,PetscScalar *x,PetscReal *dct2)
{
  PetscInt  j,i;
  PetscReal t;

  PetscFunctionBegin;
  for (j=0;j<d+1;j++) {
    t = ((2*j+1)*PETSC_PI)/(2*(d+1));
    x[j] = (a+b)/2.0+((b-a)/2.0)*PetscCosReal(t);
    for (i=0;i<d+1;i++) dct2[j*(d+1)+i] = PetscCosReal(i*t);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSolve_Interpol(NEP nep)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;
  Mat            *A;
  PetscScalar    *x,*fx,t;
  PetscReal      *cs,a,b,s,aprox,aprox0=1.0,*matnorm;
  PetscInt       i,j,k,deg=ctx->maxdeg;
  PetscBool      hasmnorm=PETSC_FALSE;
  Vec            vr,vi=NULL;

  PetscFunctionBegin;
  ierr = PetscMalloc5(deg+1,&A,(deg+1)*(deg+1),&cs,deg+1,&x,(deg+1)*nep->nt,&fx,nep->nt,&matnorm);CHKERRQ(ierr);
  for  (j=0;j<nep->nt;j++) {
    ierr = MatHasOperation(nep->A[j],MATOP_NORM,&hasmnorm);CHKERRQ(ierr);
    if (!hasmnorm) break;
    ierr = MatNorm(nep->A[j],NORM_INFINITY,matnorm+j);CHKERRQ(ierr);
  }
  if (!hasmnorm) for (j=0;j<nep->nt;j++) matnorm[j] = 1.0;
  ierr = RGIntervalGetEndpoints(nep->rg,&a,&b,NULL,NULL);CHKERRQ(ierr);
  ierr = ChebyshevNodes(deg,a,b,x,cs);CHKERRQ(ierr);
  for (j=0;j<nep->nt;j++) {
    for (i=0;i<=deg;i++) {
      ierr = FNEvaluateFunction(nep->f[j],x[i],&fx[i+j*(deg+1)]);CHKERRQ(ierr);
    }
  }
  /* Polynomial coefficients */
  ctx->deg = deg;
  for (k=0;k<=deg;k++) {
    ierr = MatDuplicate(nep->A[0],MAT_COPY_VALUES,&A[k]);CHKERRQ(ierr);
    t = 0.0;
    for (i=0;i<deg+1;i++) t += fx[i]*cs[i*(deg+1)+k];
    t *= 2.0/(deg+1);
    if (k==0) t /= 2.0;
    aprox = matnorm[0]*PetscAbsScalar(t);
    ierr = MatScale(A[k],t);CHKERRQ(ierr);
    for (j=1;j<nep->nt;j++) {
      t = 0.0;
      for (i=0;i<deg+1;i++) t += fx[i+j*(deg+1)]*cs[i*(deg+1)+k];
      t *= 2.0/(deg+1);
      if (k==0) t /= 2.0;
      aprox += matnorm[j]*PetscAbsScalar(t);
      ierr = MatAXPY(A[k],t,nep->A[j],nep->mstr);CHKERRQ(ierr);
    }
    if (k==0) aprox0 = aprox;
    if (k>1 && aprox/aprox0<ctx->tol) { ctx->deg = k; deg = k; break; }
  }

  ierr = PEPSetOperators(ctx->pep,deg+1,A);CHKERRQ(ierr);
  for (k=0;k<=deg;k++) {
    ierr = MatDestroy(&A[k]);CHKERRQ(ierr);
  }
  ierr = PetscFree5(A,cs,x,fx,matnorm);CHKERRQ(ierr);

  /* Solve polynomial eigenproblem */
  ierr = PEPSolve(ctx->pep);CHKERRQ(ierr);
  ierr = PEPGetConverged(ctx->pep,&nep->nconv);CHKERRQ(ierr);
  ierr = PEPGetIterationNumber(ctx->pep,&nep->its);CHKERRQ(ierr);
  ierr = PEPGetConvergedReason(ctx->pep,(PEPConvergedReason*)&nep->reason);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(nep->V,0,nep->nconv);CHKERRQ(ierr);
  ierr = BVCreateVec(nep->V,&vr);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  ierr = VecDuplicate(vr,&vi);CHKERRQ(ierr);
#endif
  s = 2.0/(b-a);
  for (i=0;i<nep->nconv;i++) {
    ierr = PEPGetEigenpair(ctx->pep,i,&nep->eigr[i],&nep->eigi[i],vr,vi);CHKERRQ(ierr);
    nep->eigr[i] /= s;
    nep->eigr[i] += (a+b)/2.0;
    nep->eigi[i] /= s;
    ierr = BVInsertVec(nep->V,i,vr);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    if (nep->eigi[i]!=0.0) {
      ierr = BVInsertVec(nep->V,++i,vi);CHKERRQ(ierr);
    }
#endif
  }
  ierr = VecDestroy(&vr);CHKERRQ(ierr);
  ierr = VecDestroy(&vi);CHKERRQ(ierr);

  nep->state = NEP_STATE_EIGENVECTORS;
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPMonitor_Interpol(PEP pep,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscInt       i,n;
  NEP            nep = (NEP)ctx;
  PetscReal      a,b,s;
  ST             st;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  n = PetscMin(nest,nep->ncv);
  for (i=0;i<n;i++) {
    nep->eigr[i]   = eigr[i];
    nep->eigi[i]   = eigi[i];
    nep->errest[i] = errest[i];
  }
  ierr = PEPGetST(pep,&st);CHKERRQ(ierr);
  ierr = STBackTransform(st,n,nep->eigr,nep->eigi);CHKERRQ(ierr);
  ierr = RGIntervalGetEndpoints(nep->rg,&a,&b,NULL,NULL);CHKERRQ(ierr);
  s = 2.0/(b-a);
  for (i=0;i<n;i++) {
    nep->eigr[i] /= s;
    nep->eigr[i] += (a+b)/2.0;
    nep->eigi[i] /= s;
  }
  ierr = NEPMonitor(nep,its,nconv,nep->eigr,nep->eigi,nep->errest,nest);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPSetFromOptions_Interpol(PetscOptionItems *PetscOptionsObject,NEP nep)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;
  PetscInt       i;
  PetscBool      flg1,flg2;
  PetscReal      r;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"NEP Interpol Options");CHKERRQ(ierr);

    ierr = NEPInterpolGetInterpolation(nep,&r,&i);CHKERRQ(ierr);
    if (!i) i = PETSC_DEFAULT;
    ierr = PetscOptionsInt("-nep_interpol_interpolation_degree","Maximum degree of polynomial interpolation","NEPInterpolSetInterpolation",i,&i,&flg1);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-nep_interpol_interpolation_tol","Tolerance for interpolation coefficients","NEPInterpolSetInterpolation",r,&r,&flg2);CHKERRQ(ierr);
    if (flg1 || flg2) { ierr = NEPInterpolSetInterpolation(nep,r,i);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);

  if (!ctx->pep) { ierr = NEPInterpolGetPEP(nep,&ctx->pep);CHKERRQ(ierr); }
  ierr = PEPSetFromOptions(ctx->pep);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPInterpolSetInterpolation_Interpol(NEP nep,PetscReal tol,PetscInt degree)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  if (tol == PETSC_DEFAULT) {
    ctx->tol   = PETSC_DEFAULT;
    nep->state = NEP_STATE_INITIAL;
  } else {
    if (tol <= 0.0) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of tol. Must be > 0");
    ctx->tol = tol;
  }
  if (degree == PETSC_DEFAULT || degree == PETSC_DECIDE) {
    ctx->maxdeg = 0;
    if (nep->state) { ierr = NEPReset(nep);CHKERRQ(ierr); }
    nep->state = NEP_STATE_INITIAL;
  } else {
    if (degree <= 0) SETERRQ(PetscObjectComm((PetscObject)nep),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of degree. Must be > 0");
    if (ctx->maxdeg != degree) {
      ctx->maxdeg = degree;
      if (nep->state) { ierr = NEPReset(nep);CHKERRQ(ierr); }
      nep->state = NEP_STATE_INITIAL;
    }
  }
  PetscFunctionReturn(0);
}

/*@
   NEPInterpolSetInterpolation - Sets the tolerance and maximum degree when building
   the interpolation polynomial.

   Collective on nep

   Input Parameters:
+  nep - nonlinear eigenvalue solver
.  tol - tolerance to stop computing polynomial coefficients
-  deg - maximum degree of interpolation

   Options Database Key:
+  -nep_interpol_interpolation_tol <tol> - Sets the tolerance to stop computing polynomial coefficients
-  -nep_interpol_interpolation_degree <degree> - Sets the maximum degree of interpolation

   Notes:
   Use PETSC_DEFAULT for either argument to assign a reasonably good value.

   Level: advanced

.seealso: NEPInterpolGetInterpolation()
@*/
PetscErrorCode NEPInterpolSetInterpolation(NEP nep,PetscReal tol,PetscInt deg)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidLogicalCollectiveReal(nep,tol,2);
  PetscValidLogicalCollectiveInt(nep,deg,3);
  ierr = PetscTryMethod(nep,"NEPInterpolSetInterpolation_C",(NEP,PetscReal,PetscInt),(nep,tol,deg));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPInterpolGetInterpolation_Interpol(NEP nep,PetscReal *tol,PetscInt *deg)
{
  NEP_INTERPOL *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  if (tol) *tol = ctx->tol;
  if (deg) *deg = ctx->maxdeg;
  PetscFunctionReturn(0);
}

/*@
   NEPInterpolGetInterpolation - Gets the tolerance and maximum degree when building
   the interpolation polynomial.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
+  tol - tolerance to stop computing polynomial coefficients
-  deg - maximum degree of interpolation

   Level: advanced

.seealso: NEPInterpolSetInterpolation()
@*/
PetscErrorCode NEPInterpolGetInterpolation(NEP nep,PetscReal *tol,PetscInt *deg)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  ierr = PetscUseMethod(nep,"NEPInterpolGetInterpolation_C",(NEP,PetscReal*,PetscInt*),(nep,tol,deg));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPInterpolSetPEP_Interpol(NEP nep,PEP pep)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)pep);CHKERRQ(ierr);
  ierr = PEPDestroy(&ctx->pep);CHKERRQ(ierr);
  ctx->pep = pep;
  ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->pep);CHKERRQ(ierr);
  nep->state = NEP_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   NEPInterpolSetPEP - Associate a polynomial eigensolver object (PEP) to the
   nonlinear eigenvalue solver.

   Collective on nep

   Input Parameters:
+  nep - nonlinear eigenvalue solver
-  pep - the polynomial eigensolver object

   Level: advanced

.seealso: NEPInterpolGetPEP()
@*/
PetscErrorCode NEPInterpolSetPEP(NEP nep,PEP pep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidHeaderSpecific(pep,PEP_CLASSID,2);
  PetscCheckSameComm(nep,1,pep,2);
  ierr = PetscTryMethod(nep,"NEPInterpolSetPEP_C",(NEP,PEP),(nep,pep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode NEPInterpolGetPEP_Interpol(NEP nep,PEP *pep)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  if (!ctx->pep) {
    ierr = PEPCreate(PetscObjectComm((PetscObject)nep),&ctx->pep);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)ctx->pep,(PetscObject)nep,1);CHKERRQ(ierr);
    ierr = PEPSetOptionsPrefix(ctx->pep,((PetscObject)nep)->prefix);CHKERRQ(ierr);
    ierr = PEPAppendOptionsPrefix(ctx->pep,"nep_interpol_");CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)nep,(PetscObject)ctx->pep);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)ctx->pep,((PetscObject)nep)->options);CHKERRQ(ierr);
    ierr = PEPMonitorSet(ctx->pep,PEPMonitor_Interpol,nep,NULL);CHKERRQ(ierr);
  }
  *pep = ctx->pep;
  PetscFunctionReturn(0);
}

/*@
   NEPInterpolGetPEP - Retrieve the polynomial eigensolver object (PEP)
   associated with the nonlinear eigenvalue solver.

   Not Collective

   Input Parameter:
.  nep - nonlinear eigenvalue solver

   Output Parameter:
.  pep - the polynomial eigensolver object

   Level: advanced

.seealso: NEPInterpolSetPEP()
@*/
PetscErrorCode NEPInterpolGetPEP(NEP nep,PEP *pep)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(nep,NEP_CLASSID,1);
  PetscValidPointer(pep,2);
  ierr = PetscUseMethod(nep,"NEPInterpolGetPEP_C",(NEP,PEP*),(nep,pep));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPView_Interpol(NEP nep,PetscViewer viewer)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    if (!ctx->pep) { ierr = NEPInterpolGetPEP(nep,&ctx->pep);CHKERRQ(ierr); }
    ierr = PetscViewerASCIIPrintf(viewer,"  polynomial degree %D, max=%D\n",ctx->deg,ctx->maxdeg);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  tolerance for norm of polynomial coefficients %g\n",ctx->tol);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = PEPView(ctx->pep,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NEPReset_Interpol(NEP nep)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  ierr = PEPReset(ctx->pep);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode NEPDestroy_Interpol(NEP nep)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx = (NEP_INTERPOL*)nep->data;

  PetscFunctionBegin;
  ierr = PEPDestroy(&ctx->pep);CHKERRQ(ierr);
  ierr = PetscFree(nep->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolSetInterpolation_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolGetInterpolation_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolSetPEP_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolGetPEP_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode NEPCreate_Interpol(NEP nep)
{
  PetscErrorCode ierr;
  NEP_INTERPOL   *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(nep,&ctx);CHKERRQ(ierr);
  nep->data   = (void*)ctx;
  ctx->maxdeg = 5;
  ctx->tol    = PETSC_DEFAULT;

  nep->ops->solve          = NEPSolve_Interpol;
  nep->ops->setup          = NEPSetUp_Interpol;
  nep->ops->setfromoptions = NEPSetFromOptions_Interpol;
  nep->ops->reset          = NEPReset_Interpol;
  nep->ops->destroy        = NEPDestroy_Interpol;
  nep->ops->view           = NEPView_Interpol;

  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolSetInterpolation_C",NEPInterpolSetInterpolation_Interpol);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolGetInterpolation_C",NEPInterpolGetInterpolation_Interpol);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolSetPEP_C",NEPInterpolSetPEP_Interpol);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)nep,"NEPInterpolGetPEP_C",NEPInterpolGetPEP_Interpol);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

