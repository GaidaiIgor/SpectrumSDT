#include <petsc/private/tshistoryimpl.h> /*I "petscts.h" I*/

typedef struct {
  TSHistory hist;
  PetscBool bw;
} TSAdapt_History;

static PetscErrorCode TSAdaptChoose_History(TSAdapt adapt,TS ts,PetscReal h,PetscInt *next_sc,PetscReal *next_h,PetscBool *accept,PetscReal *wlte,PetscReal *wltea,PetscReal *wlter)
{
  PetscErrorCode  ierr;
  PetscInt        step;
  TSAdapt_History *thadapt = (TSAdapt_History*)adapt->data;

  PetscFunctionBegin;
  if (!thadapt->hist) SETERRQ(PetscObjectComm((PetscObject)adapt),PETSC_ERR_ORDER,"Need to call TSAdaptHistorySetHistory() first");
  ierr = TSGetStepNumber(ts,&step);CHKERRQ(ierr);
  ierr = TSHistoryGetTimeStep(thadapt->hist,thadapt->bw,step+1,next_h);CHKERRQ(ierr);
  *accept  = PETSC_TRUE;
  *next_sc = 0;
  *wlte    = -1;
  *wltea   = -1;
  *wlter   = -1;
  PetscFunctionReturn(0);
}

static PetscErrorCode TSAdaptReset_History(TSAdapt adapt)
{
  TSAdapt_History *thadapt = (TSAdapt_History*)adapt->data;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = TSHistoryDestroy(&thadapt->hist);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode TSAdaptDestroy_History(TSAdapt adapt)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = TSAdaptReset_History(adapt);CHKERRQ(ierr);
  ierr = PetscFree(adapt->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* this is not public, as TSHistory is not a public object */
PetscErrorCode TSAdaptHistorySetTSHistory(TSAdapt adapt, TSHistory hist, PetscBool backward)
{
  PetscReal      *hist_t;
  PetscInt       n;
  PetscBool      flg;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(adapt,TSADAPT_CLASSID,1);
  PetscValidLogicalCollectiveBool(adapt,backward,3);
  ierr = PetscObjectTypeCompare((PetscObject)adapt,TSADAPTHISTORY,&flg);CHKERRQ(ierr);
  if (!flg) PetscFunctionReturn(0);
  ierr = TSHistoryGetHistory(hist,&n,(const PetscReal**)&hist_t,NULL,NULL);CHKERRQ(ierr);
  ierr = TSAdaptHistorySetHistory(adapt,n,hist_t,backward);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   TSAdaptHistoryGetStep - Gets time and time step for a given step number in the history

   Logically Collective on TSAdapt

   Input Parameters:
+  adapt    - the TSAdapt context
-  step     - the step number

   Output Parameters:
+  t  - the time corresponding to the requested step (can be NULL)
-  dt - the time step to be taken at the requested step (can be NULL)

   Notes: The time history is internally copied, and the user can free the hist array. The user still needs to specify the termination of the solve via TSSetMaxSteps().

   Level: advanced

.seealso: TSGetAdapt(), TSAdaptSetType(), TSAdaptHistorySetTrajectory(), TSADAPTHISTORY
@*/
PetscErrorCode TSAdaptHistoryGetStep(TSAdapt adapt, PetscInt step, PetscReal *t, PetscReal *dt)
{
  TSAdapt_History *thadapt;
  PetscBool       flg;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(adapt,TSADAPT_CLASSID,1);
  PetscValidLogicalCollectiveInt(adapt,step,2);
  ierr = PetscObjectTypeCompare((PetscObject)adapt,TSADAPTHISTORY,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ1(PetscObjectComm((PetscObject)adapt),PETSC_ERR_SUP,"Not for type %s",((PetscObject)adapt)->type_name);
  thadapt = (TSAdapt_History*)adapt->data;
  ierr = TSHistoryGetTimeStep(thadapt->hist,thadapt->bw,step,dt);CHKERRQ(ierr);
  ierr = TSHistoryGetTime(thadapt->hist,thadapt->bw,step,t);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   TSAdaptHistorySetHistory - Sets a time history in the adaptor

   Logically Collective on TSAdapt

   Input Parameters:
+  adapt    - the TSAdapt context
.  n        - size of the time history
.  hist     - the time history
-  backward - if the time history has to be followed backward

   Notes: The time history is internally copied, and the user can free the hist array. The user still needs to specify the termination of the solve via TSSetMaxSteps().

   Level: advanced

.seealso: TSGetAdapt(), TSAdaptSetType(), TSAdaptHistorySetTrajectory(), TSADAPTHISTORY
@*/
PetscErrorCode TSAdaptHistorySetHistory(TSAdapt adapt, PetscInt n, PetscReal hist[], PetscBool backward)
{
  TSAdapt_History *thadapt;
  PetscBool       flg;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(adapt,TSADAPT_CLASSID,1);
  PetscValidLogicalCollectiveInt(adapt,n,2);
  PetscValidRealPointer(hist,3);
  PetscValidLogicalCollectiveBool(adapt,backward,4);
  ierr = PetscObjectTypeCompare((PetscObject)adapt,TSADAPTHISTORY,&flg);CHKERRQ(ierr);
  if (!flg) PetscFunctionReturn(0);
  thadapt = (TSAdapt_History*)adapt->data;
  ierr = TSHistoryDestroy(&thadapt->hist);CHKERRQ(ierr);
  ierr = TSHistoryCreate(PetscObjectComm((PetscObject)adapt),&thadapt->hist);CHKERRQ(ierr);
  ierr = TSHistorySetHistory(thadapt->hist,n,hist,NULL,PETSC_FALSE);CHKERRQ(ierr);
  thadapt->bw = backward;
  PetscFunctionReturn(0);
}

/*@
   TSAdaptHistorySetTrajectory - Sets a time history in the adaptor from a given TSTrajectory

   Logically Collective on TSAdapt

   Input Parameters:
+  adapt    - the TSAdapt context
.  tj       - the TSTrajectory context
-  backward - if the time history has to be followed backward

   Notes: The time history is internally copied, and the user can destroy the TSTrajectory if not needed. The user still needs to specify the termination of the solve via TSSetMaxSteps().

   Level: advanced

.seealso: TSGetAdapt(), TSAdaptSetType(), TSAdaptHistorySetHistory(), TSADAPTHISTORY
@*/
PetscErrorCode TSAdaptHistorySetTrajectory(TSAdapt adapt, TSTrajectory tj, PetscBool backward)
{
  PetscBool       flg;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(adapt,TSADAPT_CLASSID,1);
  PetscValidHeaderSpecific(tj,TSTRAJECTORY_CLASSID,2);
  PetscValidLogicalCollectiveBool(adapt,backward,3);
  ierr = PetscObjectTypeCompare((PetscObject)adapt,TSADAPTHISTORY,&flg);CHKERRQ(ierr);
  if (!flg) PetscFunctionReturn(0);
  ierr = TSAdaptHistorySetTSHistory(adapt,tj->tsh,backward);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*MC
   TSADAPTHISTORY - Time stepping controller that follows a given time history, used for Tangent Linear Model simulations

   Level: developer

.seealso: TS, TSAdapt, TSGetAdapt(), TSAdaptHistorySetHistory()
M*/
PETSC_EXTERN PetscErrorCode TSAdaptCreate_History(TSAdapt adapt)
{
  PetscErrorCode     ierr;
  TSAdapt_History *thadapt;

  PetscFunctionBegin;
  ierr = PetscNew(&thadapt);CHKERRQ(ierr);
  adapt->matchstepfac[0] = PETSC_SMALL; /* prevent from accumulation errors */
  adapt->matchstepfac[1] = 0.0;         /* we will always match the final step, prevent TSAdaptChoose to mess with it */
  adapt->data            = thadapt;

  adapt->ops->choose  = TSAdaptChoose_History;
  adapt->ops->reset   = TSAdaptReset_History;
  adapt->ops->destroy = TSAdaptDestroy_History;
  PetscFunctionReturn(0);
}
