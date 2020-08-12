#include <../src/snes/impls/gs/gsimpl.h>      /*I "petscsnes.h"  I*/

/*@
   SNESNGSSetTolerances - Sets various parameters used in convergence tests.

   Logically Collective on SNES

   Input Parameters:
+  snes - the SNES context
.  abstol - absolute convergence tolerance
.  rtol - relative convergence tolerance
.  stol -  convergence tolerance in terms of the norm of the change in the solution between steps,  || delta x || < stol*|| x ||
-  maxit - maximum number of iterations


   Options Database Keys:
+    -snes_ngs_atol <abstol> - Sets abstol
.    -snes_ngs_rtol <rtol> - Sets rtol
.    -snes_ngs_stol <stol> - Sets stol
-    -snes_max_it <maxit> - Sets maxit

   Level: intermediate

.seealso: SNESSetTrustRegionTolerance()
@*/
PetscErrorCode  SNESNGSSetTolerances(SNES snes,PetscReal abstol,PetscReal rtol,PetscReal stol,PetscInt maxit)
{
  SNES_NGS *gs = (SNES_NGS*)snes->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(snes,SNES_CLASSID,1);

  if (abstol != PETSC_DEFAULT) {
    if (abstol < 0.0) SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_ARG_OUTOFRANGE,"Absolute tolerance %g must be non-negative",(double)abstol);
    gs->abstol = abstol;
  }
  if (rtol != PETSC_DEFAULT) {
    if (rtol < 0.0 || 1.0 <= rtol) SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_ARG_OUTOFRANGE,"Relative tolerance %g must be non-negative and less than 1.0",(double)rtol);
    gs->rtol = rtol;
  }
  if (stol != PETSC_DEFAULT) {
    if (stol < 0.0) SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_ARG_OUTOFRANGE,"Step tolerance %g must be non-negative",(double)stol);
    gs->stol = stol;
  }
  if (maxit != PETSC_DEFAULT) {
    if (maxit < 0) SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_ARG_OUTOFRANGE,"Maximum number of iterations %D must be non-negative",maxit);
    gs->max_its = maxit;
  }
  PetscFunctionReturn(0);
}

/*@
   SNESNGSGetTolerances - Gets various parameters used in convergence tests.

   Not Collective

   Input Parameters:
+  snes - the SNES context
.  atol - absolute convergence tolerance
.  rtol - relative convergence tolerance
.  stol -  convergence tolerance in terms of the norm
           of the change in the solution between steps
-  maxit - maximum number of iterations

   Notes:
   The user can specify NULL for any parameter that is not needed.

   Level: intermediate

.seealso: SNESSetTolerances()
@*/
PetscErrorCode  SNESNGSGetTolerances(SNES snes,PetscReal *atol,PetscReal *rtol,PetscReal *stol,PetscInt *maxit)
{
  SNES_NGS *gs = (SNES_NGS*)snes->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  if (atol) *atol = gs->abstol;
  if (rtol) *rtol = gs->rtol;
  if (stol) *stol = gs->stol;
  if (maxit) *maxit = gs->max_its;
  PetscFunctionReturn(0);
}

/*@
   SNESNGSSetSweeps - Sets the number of sweeps of GS to use.

   Input Parameters:
+  snes   - the SNES context
-  sweeps  - the number of sweeps of GS to perform.

   Level: intermediate

.seealso: SNESSetNGS(), SNESGetNGS(), SNESSetNPC(), SNESNGSGetSweeps()
@*/

PetscErrorCode SNESNGSSetSweeps(SNES snes, PetscInt sweeps)
{
  SNES_NGS *gs = (SNES_NGS*)snes->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  gs->sweeps = sweeps;
  PetscFunctionReturn(0);
}

/*@
   SNESNGSGetSweeps - Gets the number of sweeps GS will use.

   Input Parameters:
.  snes   - the SNES context

   Output Parameters:
.  sweeps  - the number of sweeps of GS to perform.

   Level: intermediate

.seealso: SNESSetNGS(), SNESGetNGS(), SNESSetNPC(), SNESNGSSetSweeps()
@*/
PetscErrorCode SNESNGSGetSweeps(SNES snes, PetscInt * sweeps)
{
  SNES_NGS *gs = (SNES_NGS*)snes->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  *sweeps = gs->sweeps;
  PetscFunctionReturn(0);
}

PetscErrorCode SNESReset_NGS(SNES snes)
{
  SNES_NGS       *gs = (SNES_NGS*)snes->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ISColoringDestroy(&gs->coloring);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SNESDestroy_NGS(SNES snes)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SNESReset_NGS(snes);CHKERRQ(ierr);
  ierr = PetscFree(snes->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SNESSetUp_NGS(SNES snes)
{
  PetscErrorCode ierr;
  PetscErrorCode (*f)(SNES,Vec,Vec,void*);

  PetscFunctionBegin;
  ierr = SNESGetNGS(snes,&f,NULL);CHKERRQ(ierr);
  if (!f) {
    ierr = SNESSetNGS(snes,SNESComputeNGSDefaultSecant,NULL);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SNESSetFromOptions_NGS(PetscOptionItems *PetscOptionsObject,SNES snes)
{
  SNES_NGS       *gs = (SNES_NGS*)snes->data;
  PetscErrorCode ierr;
  PetscInt       sweeps,max_its=PETSC_DEFAULT;
  PetscReal      rtol=PETSC_DEFAULT,atol=PETSC_DEFAULT,stol=PETSC_DEFAULT;
  PetscBool      flg,flg1,flg2,flg3;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"SNES GS options");CHKERRQ(ierr);
  /* GS Options */
  ierr = PetscOptionsInt("-snes_ngs_sweeps","Number of sweeps of GS to apply","SNESComputeGS",gs->sweeps,&sweeps,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = SNESNGSSetSweeps(snes,sweeps);CHKERRQ(ierr);
  }
  ierr = PetscOptionsReal("-snes_ngs_atol","Absolute residual tolerance for GS iteration","SNESComputeGS",gs->abstol,&atol,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-snes_ngs_rtol","Relative residual tolerance for GS iteration","SNESComputeGS",gs->rtol,&rtol,&flg1);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-snes_ngs_stol","Absolute update tolerance for GS iteration","SNESComputeGS",gs->stol,&stol,&flg2);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-snes_ngs_max_it","Maximum number of sweeps of GS to apply","SNESComputeGS",gs->max_its,&max_its,&flg3);CHKERRQ(ierr);
  if (flg || flg1 || flg2 || flg3) {
    ierr = SNESNGSSetTolerances(snes,atol,rtol,stol,max_its);CHKERRQ(ierr);
  }
  flg  = PETSC_FALSE;
  ierr = PetscOptionsBool("-snes_ngs_secant","Use finite difference secant approximation with coloring","",flg,&flg,NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = SNESSetNGS(snes,SNESComputeNGSDefaultSecant,NULL);CHKERRQ(ierr);
    ierr = PetscInfo(snes,"Setting default finite difference secant approximation with coloring\n");CHKERRQ(ierr);
  }
  ierr = PetscOptionsReal("-snes_ngs_secant_h","Differencing parameter for secant search","",gs->h,&gs->h,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-snes_ngs_secant_mat_coloring","Use the graph coloring of the Jacobian for the secant GS","",gs->secant_mat,&gs->secant_mat,&flg);CHKERRQ(ierr);

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SNESView_NGS(SNES snes, PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscErrorCode (*f)(SNES,Vec,Vec,void*);
  SNES_NGS       *gs = (SNES_NGS*)snes->data;
  PetscBool      iascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = DMSNESGetNGS(snes->dm,&f,NULL);CHKERRQ(ierr);
    if (f == SNESComputeNGSDefaultSecant) {
      ierr = PetscViewerASCIIPrintf(viewer,"  Use finite difference secant approximation with coloring with h = %g \n",(double)gs->h);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SNESSolve_NGS(SNES snes)
{
  Vec              F;
  Vec              X;
  Vec              B;
  PetscInt         i;
  PetscReal        fnorm;
  PetscErrorCode   ierr;
  SNESNormSchedule normschedule;

  PetscFunctionBegin;

  if (snes->xl || snes->xu || snes->ops->computevariablebounds) SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_ARG_WRONGSTATE, "SNES solver %s does not support bounds", ((PetscObject)snes)->type_name);

  ierr = PetscCitationsRegister(SNESCitation,&SNEScite);CHKERRQ(ierr);
  X = snes->vec_sol;
  F = snes->vec_func;
  B = snes->vec_rhs;

  ierr         = PetscObjectSAWsTakeAccess((PetscObject)snes);CHKERRQ(ierr);
  snes->iter   = 0;
  snes->norm   = 0.;
  ierr         = PetscObjectSAWsGrantAccess((PetscObject)snes);CHKERRQ(ierr);
  snes->reason = SNES_CONVERGED_ITERATING;

  ierr = SNESGetNormSchedule(snes, &normschedule);CHKERRQ(ierr);
  if (normschedule == SNES_NORM_ALWAYS || normschedule == SNES_NORM_INITIAL_ONLY || normschedule == SNES_NORM_INITIAL_FINAL_ONLY) {
    /* compute the initial function and preconditioned update delX */
    if (!snes->vec_func_init_set) {
      ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
    } else snes->vec_func_init_set = PETSC_FALSE;

    ierr = VecNorm(F, NORM_2, &fnorm);CHKERRQ(ierr); /* fnorm <- ||F||  */
    SNESCheckFunctionNorm(snes,fnorm);
    ierr       = PetscObjectSAWsTakeAccess((PetscObject)snes);CHKERRQ(ierr);
    snes->iter = 0;
    snes->norm = fnorm;
    ierr       = PetscObjectSAWsGrantAccess((PetscObject)snes);CHKERRQ(ierr);
    ierr       = SNESLogConvergenceHistory(snes,snes->norm,0);CHKERRQ(ierr);
    ierr       = SNESMonitor(snes,0,snes->norm);CHKERRQ(ierr);

    /* test convergence */
    ierr = (*snes->ops->converged)(snes,0,0.0,0.0,fnorm,&snes->reason,snes->cnvP);CHKERRQ(ierr);
    if (snes->reason) PetscFunctionReturn(0);
  } else {
    ierr = PetscObjectSAWsGrantAccess((PetscObject)snes);CHKERRQ(ierr);
    ierr = SNESLogConvergenceHistory(snes,snes->norm,0);CHKERRQ(ierr);
  }

  /* Call general purpose update function */
  if (snes->ops->update) {
    ierr = (*snes->ops->update)(snes, snes->iter);CHKERRQ(ierr);
  }

  for (i = 0; i < snes->max_its; i++) {
    ierr = SNESComputeNGS(snes, B, X);CHKERRQ(ierr);
    /* only compute norms if requested or about to exit due to maximum iterations */
    if (normschedule == SNES_NORM_ALWAYS || ((i == snes->max_its - 1) && (normschedule == SNES_NORM_INITIAL_FINAL_ONLY || normschedule == SNES_NORM_FINAL_ONLY))) {
      ierr = SNESComputeFunction(snes,X,F);CHKERRQ(ierr);
      ierr = VecNorm(F, NORM_2, &fnorm);CHKERRQ(ierr); /* fnorm <- ||F||  */
      SNESCheckFunctionNorm(snes,fnorm);
      /* Monitor convergence */
      ierr       = PetscObjectSAWsTakeAccess((PetscObject)snes);CHKERRQ(ierr);
      snes->iter = i+1;
      snes->norm = fnorm;
      ierr       = PetscObjectSAWsGrantAccess((PetscObject)snes);CHKERRQ(ierr);
      ierr       = SNESLogConvergenceHistory(snes,snes->norm,0);CHKERRQ(ierr);
      ierr       = SNESMonitor(snes,snes->iter,snes->norm);CHKERRQ(ierr);
    }
    /* Test for convergence */
    if (normschedule == SNES_NORM_ALWAYS) ierr = (*snes->ops->converged)(snes,snes->iter,0.0,0.0,fnorm,&snes->reason,snes->cnvP);CHKERRQ(ierr);
    if (snes->reason) PetscFunctionReturn(0);
    /* Call general purpose update function */
    if (snes->ops->update) {
      ierr = (*snes->ops->update)(snes, snes->iter);CHKERRQ(ierr);
    }
  }
  if (normschedule == SNES_NORM_ALWAYS) {
    if (i == snes->max_its) {
      ierr = PetscInfo1(snes,"Maximum number of iterations has been reached: %D\n",snes->max_its);CHKERRQ(ierr);
      if (!snes->reason) snes->reason = SNES_DIVERGED_MAX_IT;
    }
  } else if (!snes->reason) snes->reason = SNES_CONVERGED_ITS; /* GS is meant to be used as a preconditioner */
  PetscFunctionReturn(0);
}

/*MC
  SNESNGS - Either calls the user-provided solution routine provided with SNESSetNGS() or does a finite difference secant approximation
            using coloring.

   Level: advanced

  Options Database:
+   -snes_ngs_sweeps <n> - Number of sweeps of GS to apply
.    -snes_ngs_atol <atol> - Absolute residual tolerance for GS iteration
.    -snes_ngs_rtol <rtol> - Relative residual tolerance for GS iteration
.    -snes_ngs_stol <stol> - Absolute update tolerance for GS iteration
.    -snes_ngs_max_it <maxit> - Maximum number of sweeps of GS to apply
.    -snes_ngs_secant - Use pointwise secant local Jacobian approximation with coloring instead of user provided Gauss-Seidel routine, this is
                        used by default if no user provided Gauss-Seidel routine is available. Requires either that a DM that can compute a coloring
                        is available or a Jacobian sparse matrix is provided (from which to get the coloring).
.    -snes_ngs_secant_h <h> - Differencing parameter for secant approximation
.    -snes_ngs_secant_mat_coloring - Use the graph coloring of the Jacobian for the secant GS even if a DM is available.
-    -snes_norm_schedule <none, always, initialonly, finalonly, initalfinalonly> - how often the residual norms are computed

  Notes:
  the Gauss-Seidel smoother is inherited through composition.  If a solver has been created with SNESGetNPC(), it will have
  its parent's Gauss-Seidel routine associated with it.

  By default this routine computes the solution norm at each iteration, this can be time consuming, you can turn this off with SNESSetNormSchedule() or -snes_norm_schedule
   References:
.  1. - Peter R. Brune, Matthew G. Knepley, Barry F. Smith, and Xuemin Tu, "Composing Scalable Nonlinear Algebraic Solvers",
   SIAM Review, 57(4), 2015

.seealso: SNESCreate(), SNES, SNESSetType(), SNESSetNGS(), SNESType (for list of available types), SNESNGSSetSweeps(), SNESNGSSetTolerances(),
          SNESSetNormSchedule()
M*/

PETSC_EXTERN PetscErrorCode SNESCreate_NGS(SNES snes)
{
  SNES_NGS        *gs;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  snes->ops->destroy        = SNESDestroy_NGS;
  snes->ops->setup          = SNESSetUp_NGS;
  snes->ops->setfromoptions = SNESSetFromOptions_NGS;
  snes->ops->view           = SNESView_NGS;
  snes->ops->solve          = SNESSolve_NGS;
  snes->ops->reset          = SNESReset_NGS;

  snes->usesksp = PETSC_FALSE;
  snes->usesnpc = PETSC_FALSE;

  snes->alwayscomputesfinalresidual = PETSC_FALSE;

  if (!snes->tolerancesset) {
    snes->max_its   = 10000;
    snes->max_funcs = 10000;
  }

  ierr = PetscNewLog(snes,&gs);CHKERRQ(ierr);

  gs->sweeps  = 1;
  gs->rtol    = 1e-5;
  gs->abstol  = PETSC_MACHINE_EPSILON;
  gs->stol    = 1000*PETSC_MACHINE_EPSILON;
  gs->max_its = 50;
  gs->h       = PETSC_SQRT_MACHINE_EPSILON;

  snes->data = (void*) gs;
  PetscFunctionReturn(0);
}
