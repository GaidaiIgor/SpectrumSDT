/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This provides a simple shell interface for programmers to create
   their own spectral transformations without writing much interface code
*/

#include <slepc/private/stimpl.h>        /*I "slepcst.h" I*/

typedef struct {
  void           *ctx;                       /* user provided context */
  PetscErrorCode (*apply)(ST,Vec,Vec);
  PetscErrorCode (*applytrans)(ST,Vec,Vec);
  PetscErrorCode (*backtransform)(ST,PetscInt n,PetscScalar*,PetscScalar*);
} ST_SHELL;

/*@C
   STShellGetContext - Returns the user-provided context associated with a shell ST

   Not Collective

   Input Parameter:
.  st - spectral transformation context

   Output Parameter:
.  ctx - the user provided context

   Level: advanced

   Notes:
   This routine is intended for use within various shell routines

.seealso: STShellSetContext()
@*/
PetscErrorCode STShellGetContext(ST st,void **ctx)
{
  PetscErrorCode ierr;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidPointer(ctx,2);
  ierr = PetscObjectTypeCompare((PetscObject)st,STSHELL,&flg);CHKERRQ(ierr);
  if (!flg) *ctx = 0;
  else      *ctx = ((ST_SHELL*)(st->data))->ctx;
  PetscFunctionReturn(0);
}

/*@
   STShellSetContext - Sets the context for a shell ST

   Logically Collective on st

   Input Parameters:
+  st - the shell ST
-  ctx - the context

   Level: advanced

   Fortran Notes:
   To use this from Fortran you must write a Fortran interface definition
   for this function that tells Fortran the Fortran derived data type that
   you are passing in as the ctx argument.

.seealso: STShellGetContext()
@*/
PetscErrorCode STShellSetContext(ST st,void *ctx)
{
  ST_SHELL       *shell = (ST_SHELL*)st->data;
  PetscErrorCode ierr;
  PetscBool      flg;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  ierr = PetscObjectTypeCompare((PetscObject)st,STSHELL,&flg);CHKERRQ(ierr);
  if (flg) shell->ctx = ctx;
  PetscFunctionReturn(0);
}

PetscErrorCode STApply_Shell(ST st,Vec x,Vec y)
{
  PetscErrorCode   ierr;
  ST_SHELL         *shell = (ST_SHELL*)st->data;
  PetscObjectState instate,outstate;

  PetscFunctionBegin;
  if (!shell->apply) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_USER,"No apply() routine provided to Shell ST");
  ierr = PetscObjectStateGet((PetscObject)y,&instate);CHKERRQ(ierr);
  PetscStackCall("STSHELL user function apply()",ierr = (*shell->apply)(st,x,y);CHKERRQ(ierr));
  ierr = PetscObjectStateGet((PetscObject)y,&outstate);CHKERRQ(ierr);
  if (instate == outstate) {
    /* user forgot to increase the state of the output vector */
    ierr = PetscObjectStateIncrease((PetscObject)y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode STApplyTranspose_Shell(ST st,Vec x,Vec y)
{
  PetscErrorCode ierr;
  ST_SHELL       *shell = (ST_SHELL*)st->data;
  PetscObjectState instate,outstate;

  PetscFunctionBegin;
  if (!shell->applytrans) SETERRQ(PetscObjectComm((PetscObject)st),PETSC_ERR_USER,"No applytranspose() routine provided to Shell ST");
  ierr = PetscObjectStateGet((PetscObject)y,&instate);CHKERRQ(ierr);
  PetscStackCall("STSHELL user function applytrans()",ierr = (*shell->applytrans)(st,x,y);CHKERRQ(ierr));
  ierr = PetscObjectStateGet((PetscObject)y,&outstate);CHKERRQ(ierr);
  if (instate == outstate) {
    /* user forgot to increase the state of the output vector */
    ierr = PetscObjectStateIncrease((PetscObject)y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode STBackTransform_Shell(ST st,PetscInt n,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscErrorCode ierr;
  ST_SHELL       *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  if (shell->backtransform) PetscStackCall("STSHELL user function backtransform()",ierr = (*shell->backtransform)(st,n,eigr,eigi);CHKERRQ(ierr));
  PetscFunctionReturn(0);
}

/*
   STIsInjective_Shell - Check if the user has provided the backtransform operation.
*/
PetscErrorCode STIsInjective_Shell(ST st,PetscBool* is)
{
  ST_SHELL *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  *is = shell->backtransform? PETSC_TRUE: PETSC_FALSE;
  PetscFunctionReturn(0);
}

PetscErrorCode STDestroy_Shell(ST st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(st->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STShellSetApply_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STShellSetApplyTranspose_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STShellSetBackTransform_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STShellSetApply_Shell(ST st,PetscErrorCode (*apply)(ST,Vec,Vec))
{
  ST_SHELL *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  shell->apply = apply;
  PetscFunctionReturn(0);
}

/*@C
   STShellSetApply - Sets routine to use as the application of the
   operator to a vector in the user-defined spectral transformation.

   Logically Collective on st

   Input Parameters:
+  st    - the spectral transformation context
-  apply - the application-provided transformation routine

   Calling sequence of apply:
$  PetscErrorCode apply(ST st,Vec xin,Vec xout)

+  st   - the spectral transformation context
.  xin  - input vector
-  xout - output vector

   Level: advanced

.seealso: STShellSetBackTransform(), STShellSetApplyTranspose()
@*/
PetscErrorCode STShellSetApply(ST st,PetscErrorCode (*apply)(ST,Vec,Vec))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  ierr = PetscTryMethod(st,"STShellSetApply_C",(ST,PetscErrorCode (*)(ST,Vec,Vec)),(st,apply));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STShellSetApplyTranspose_Shell(ST st,PetscErrorCode (*applytrans)(ST,Vec,Vec))
{
  ST_SHELL *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  shell->applytrans = applytrans;
  PetscFunctionReturn(0);
}

/*@C
   STShellSetApplyTranspose - Sets routine to use as the application of the
   transposed operator to a vector in the user-defined spectral transformation.

   Logically Collective on st

   Input Parameters:
+  st    - the spectral transformation context
-  applytrans - the application-provided transformation routine

   Calling sequence of applytrans:
$  PetscErrorCode applytrans(ST st,Vec xin,Vec xout)

+  st   - the spectral transformation context
.  xin  - input vector
-  xout - output vector

   Level: advanced

.seealso: STShellSetApply(), STShellSetBackTransform()
@*/
PetscErrorCode STShellSetApplyTranspose(ST st,PetscErrorCode (*applytrans)(ST,Vec,Vec))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  ierr = PetscTryMethod(st,"STShellSetApplyTranspose_C",(ST,PetscErrorCode (*)(ST,Vec,Vec)),(st,applytrans));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode STShellSetBackTransform_Shell(ST st,PetscErrorCode (*backtr)(ST,PetscInt,PetscScalar*,PetscScalar*))
{
  ST_SHELL *shell = (ST_SHELL*)st->data;

  PetscFunctionBegin;
  shell->backtransform = backtr;
  PetscFunctionReturn(0);
}

/*@C
   STShellSetBackTransform - Sets the routine to be called after the
   eigensolution process has finished in order to transform back the
   computed eigenvalues.

   Logically Collective on st

   Input Parameters:
+  st     - the spectral transformation context
-  backtr - the application-provided backtransform routine

   Calling sequence of backtr:
$  PetscErrorCode backtr(ST st,PetscScalar *eigr,PetscScalar *eigi)

+  st   - the spectral transformation context
.  eigr - pointer ot the real part of the eigenvalue to transform back
-  eigi - pointer ot the imaginary part

   Level: advanced

.seealso: STShellSetApply(), STShellSetApplyTranspose()
@*/
PetscErrorCode STShellSetBackTransform(ST st,PetscErrorCode (*backtr)(ST,PetscInt,PetscScalar*,PetscScalar*))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  ierr = PetscTryMethod(st,"STShellSetBackTransform_C",(ST,PetscErrorCode (*)(ST,PetscInt,PetscScalar*,PetscScalar*)),(st,backtr));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*MC
   STSHELL - Creates a new spectral transformation class.
          This is intended to provide a simple class to use with EPS.
          You should not use this if you plan to make a complete class.

  Level: advanced

  Usage:
$             extern PetscErrorCode (*apply)(void*,Vec,Vec);
$             extern PetscErrorCode (*applytrans)(void*,Vec,Vec);
$             extern PetscErrorCode (*backtr)(void*,PetscScalar*,PetscScalar*);
$
$             STCreate(comm,&st);
$             STSetType(st,STSHELL);
$             STShellSetContext(st,ctx);
$             STShellSetApply(st,apply);
$             STShellSetApplyTranspose(st,applytrans);  (optional)
$             STShellSetBackTransform(st,backtr);       (optional)

M*/

SLEPC_EXTERN PetscErrorCode STCreate_Shell(ST st)
{
  PetscErrorCode ierr;
  ST_SHELL       *ctx;

  PetscFunctionBegin;
  ierr = PetscNewLog(st,&ctx);CHKERRQ(ierr);
  st->data = (void*)ctx;

  st->usesksp = PETSC_FALSE;

  st->ops->apply           = STApply_Shell;
  st->ops->applytrans      = STApplyTranspose_Shell;
  st->ops->backtransform   = STBackTransform_Shell;
  st->ops->destroy         = STDestroy_Shell;

  ierr = PetscObjectComposeFunction((PetscObject)st,"STShellSetApply_C",STShellSetApply_Shell);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STShellSetApplyTranspose_C",STShellSetApplyTranspose_Shell);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)st,"STShellSetBackTransform_C",STShellSetBackTransform_Shell);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

