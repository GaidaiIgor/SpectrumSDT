
/*
    The PC (preconditioner) interface routines, callable by users.
*/
#include <petsc/private/pcimpl.h>            /*I "petscksp.h" I*/
#include <petscdm.h>

/* Logging support */
PetscClassId  PC_CLASSID;
PetscLogEvent PC_SetUp, PC_SetUpOnBlocks, PC_Apply, PC_ApplyCoarse, PC_ApplyMultiple, PC_ApplySymmetricLeft;
PetscLogEvent PC_ApplySymmetricRight, PC_ModifySubMatrices, PC_ApplyOnBlocks, PC_ApplyTransposeOnBlocks;
PetscInt      PetscMGLevelId;

PetscErrorCode PCGetDefaultType_Private(PC pc,const char *type[])
{
  PetscErrorCode ierr;
  PetscMPIInt    size;
  PetscBool      hasop,flg1,flg2,set,flg3;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pc),&size);CHKERRQ(ierr);
  if (pc->pmat) {
    ierr = MatHasOperation(pc->pmat,MATOP_GET_DIAGONAL_BLOCK,&hasop);CHKERRQ(ierr);
    if (size == 1) {
      ierr = MatGetFactorAvailable(pc->pmat,"petsc",MAT_FACTOR_ICC,&flg1);CHKERRQ(ierr);
      ierr = MatGetFactorAvailable(pc->pmat,"petsc",MAT_FACTOR_ILU,&flg2);CHKERRQ(ierr);
      ierr = MatIsSymmetricKnown(pc->pmat,&set,&flg3);CHKERRQ(ierr);
      if (flg1 && (!flg2 || (set && flg3))) {
        *type = PCICC;
      } else if (flg2) {
        *type = PCILU;
      } else if (hasop) { /* likely is a parallel matrix run on one processor */
        *type = PCBJACOBI;
      } else {
        *type = PCNONE;
      }
    } else {
       if (hasop) {
        *type = PCBJACOBI;
      } else {
        *type = PCNONE;
      }
    }
  } else {
    if (size == 1) {
      *type = PCILU;
    } else {
      *type = PCBJACOBI;
    }
  }
  PetscFunctionReturn(0);
}

/*@
   PCReset - Resets a PC context to the pcsetupcalled = 0 state and removes any allocated Vecs and Mats

   Collective on PC

   Input Parameter:
.  pc - the preconditioner context

   Level: developer

   Notes:
    This allows a PC to be reused for a different sized linear system but using the same options that have been previously set in the PC

.seealso: PCCreate(), PCSetUp()
@*/
PetscErrorCode  PCReset(PC pc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  if (pc->ops->reset) {
    ierr = (*pc->ops->reset)(pc);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&pc->diagonalscaleright);CHKERRQ(ierr);
  ierr = VecDestroy(&pc->diagonalscaleleft);CHKERRQ(ierr);
  ierr = MatDestroy(&pc->pmat);CHKERRQ(ierr);
  ierr = MatDestroy(&pc->mat);CHKERRQ(ierr);

  pc->setupcalled = 0;
  PetscFunctionReturn(0);
}

/*@
   PCDestroy - Destroys PC context that was created with PCCreate().

   Collective on PC

   Input Parameter:
.  pc - the preconditioner context

   Level: developer

.seealso: PCCreate(), PCSetUp()
@*/
PetscErrorCode  PCDestroy(PC *pc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*pc) PetscFunctionReturn(0);
  PetscValidHeaderSpecific((*pc),PC_CLASSID,1);
  if (--((PetscObject)(*pc))->refct > 0) {*pc = 0; PetscFunctionReturn(0);}

  ierr = PCReset(*pc);CHKERRQ(ierr);

  /* if memory was published with SAWs then destroy it */
  ierr = PetscObjectSAWsViewOff((PetscObject)*pc);CHKERRQ(ierr);
  if ((*pc)->ops->destroy) {ierr = (*(*pc)->ops->destroy)((*pc));CHKERRQ(ierr);}
  ierr = DMDestroy(&(*pc)->dm);CHKERRQ(ierr);
  ierr = PetscHeaderDestroy(pc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PCGetDiagonalScale - Indicates if the preconditioner applies an additional left and right
      scaling as needed by certain time-stepping codes.

   Logically Collective on PC

   Input Parameter:
.  pc - the preconditioner context

   Output Parameter:
.  flag - PETSC_TRUE if it applies the scaling

   Level: developer

   Notes:
    If this returns PETSC_TRUE then the system solved via the Krylov method is
$           D M A D^{-1} y = D M b  for left preconditioning or
$           D A M D^{-1} z = D b for right preconditioning

.seealso: PCCreate(), PCSetUp(), PCDiagonalScaleLeft(), PCDiagonalScaleRight(), PCSetDiagonalScale()
@*/
PetscErrorCode  PCGetDiagonalScale(PC pc,PetscBool  *flag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidBoolPointer(flag,2);
  *flag = pc->diagonalscale;
  PetscFunctionReturn(0);
}

/*@
   PCSetDiagonalScale - Indicates the left scaling to use to apply an additional left and right
      scaling as needed by certain time-stepping codes.

   Logically Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  s - scaling vector

   Level: intermediate

   Notes:
    The system solved via the Krylov method is
$           D M A D^{-1} y = D M b  for left preconditioning or
$           D A M D^{-1} z = D b for right preconditioning

   PCDiagonalScaleLeft() scales a vector by D. PCDiagonalScaleRight() scales a vector by D^{-1}.

.seealso: PCCreate(), PCSetUp(), PCDiagonalScaleLeft(), PCDiagonalScaleRight(), PCGetDiagonalScale()
@*/
PetscErrorCode  PCSetDiagonalScale(PC pc,Vec s)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(s,VEC_CLASSID,2);
  pc->diagonalscale     = PETSC_TRUE;

  ierr = PetscObjectReference((PetscObject)s);CHKERRQ(ierr);
  ierr = VecDestroy(&pc->diagonalscaleleft);CHKERRQ(ierr);

  pc->diagonalscaleleft = s;

  ierr = VecDuplicate(s,&pc->diagonalscaleright);CHKERRQ(ierr);
  ierr = VecCopy(s,pc->diagonalscaleright);CHKERRQ(ierr);
  ierr = VecReciprocal(pc->diagonalscaleright);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PCDiagonalScaleLeft - Scales a vector by the left scaling as needed by certain time-stepping codes.

   Logically Collective on PC

   Input Parameters:
+  pc - the preconditioner context
.  in - input vector
-  out - scaled vector (maybe the same as in)

   Level: intermediate

   Notes:
    The system solved via the Krylov method is
$           D M A D^{-1} y = D M b  for left preconditioning or
$           D A M D^{-1} z = D b for right preconditioning

   PCDiagonalScaleLeft() scales a vector by D. PCDiagonalScaleRight() scales a vector by D^{-1}.

   If diagonal scaling is turned off and in is not out then in is copied to out

.seealso: PCCreate(), PCSetUp(), PCDiagonalScaleSet(), PCDiagonalScaleRight(), PCDiagonalScale()
@*/
PetscErrorCode  PCDiagonalScaleLeft(PC pc,Vec in,Vec out)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(in,VEC_CLASSID,2);
  PetscValidHeaderSpecific(out,VEC_CLASSID,3);
  if (pc->diagonalscale) {
    ierr = VecPointwiseMult(out,pc->diagonalscaleleft,in);CHKERRQ(ierr);
  } else if (in != out) {
    ierr = VecCopy(in,out);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PCDiagonalScaleRight - Scales a vector by the right scaling as needed by certain time-stepping codes.

   Logically Collective on PC

   Input Parameters:
+  pc - the preconditioner context
.  in - input vector
-  out - scaled vector (maybe the same as in)

   Level: intermediate

   Notes:
    The system solved via the Krylov method is
$           D M A D^{-1} y = D M b  for left preconditioning or
$           D A M D^{-1} z = D b for right preconditioning

   PCDiagonalScaleLeft() scales a vector by D. PCDiagonalScaleRight() scales a vector by D^{-1}.

   If diagonal scaling is turned off and in is not out then in is copied to out

.seealso: PCCreate(), PCSetUp(), PCDiagonalScaleLeft(), PCDiagonalScaleSet(), PCDiagonalScale()
@*/
PetscErrorCode  PCDiagonalScaleRight(PC pc,Vec in,Vec out)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(in,VEC_CLASSID,2);
  PetscValidHeaderSpecific(out,VEC_CLASSID,3);
  if (pc->diagonalscale) {
    ierr = VecPointwiseMult(out,pc->diagonalscaleright,in);CHKERRQ(ierr);
  } else if (in != out) {
    ierr = VecCopy(in,out);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PCSetUseAmat - Sets a flag to indicate that when the preconditioner needs to apply (part of) the
   operator during the preconditioning process it applies the Amat provided to TSSetRHSJacobian(),
   TSSetIJacobian(), SNESSetJacobian(), KSPSetOperator() or PCSetOperator() not the Pmat.

   Logically Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  flg - PETSC_TRUE to use the Amat, PETSC_FALSE to use the Pmat (default is false)

   Options Database Key:
.  -pc_use_amat <true,false>

   Notes:
   For the common case in which the linear system matrix and the matrix used to construct the
   preconditioner are identical, this routine is does nothing.

   Level: intermediate

.seealso: PCGetUseAmat(), PCBJACOBI, PGMG, PCFIELDSPLIT, PCCOMPOSITE
@*/
PetscErrorCode  PCSetUseAmat(PC pc,PetscBool flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  pc->useAmat = flg;
  PetscFunctionReturn(0);
}

/*@
   PCSetErrorIfFailure - Causes PC to generate an error if a FPE, for example a zero pivot, is detected.

   Logically Collective on PC

   Input Parameters:
+  pc - iterative context obtained from PCCreate()
-  flg - PETSC_TRUE indicates you want the error generated

   Level: advanced

   Notes:
    Normally PETSc continues if a linear solver fails due to a failed setup of a preconditioner, you can call KSPGetConvergedReason() after a KSPSolve()
    to determine if it has converged or failed. Or use -ksp_error_if_not_converged to cause the program to terminate as soon as lack of convergence is
    detected.

    This is propagated into KSPs used by this PC, which then propagate it into PCs used by those KSPs

.seealso: PCGetInitialGuessNonzero(), PCSetInitialGuessKnoll(), PCGetInitialGuessKnoll()
@*/
PetscErrorCode  PCSetErrorIfFailure(PC pc,PetscBool flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidLogicalCollectiveBool(pc,flg,2);
  pc->erroriffailure = flg;
  PetscFunctionReturn(0);
}

/*@
   PCGetUseAmat - Gets a flag to indicate that when the preconditioner needs to apply (part of) the
   operator during the preconditioning process it applies the Amat provided to TSSetRHSJacobian(),
   TSSetIJacobian(), SNESSetJacobian(), KSPSetOperator() or PCSetOperator() not the Pmat.

   Logically Collective on PC

   Input Parameter:
.  pc - the preconditioner context

   Output Parameter:
.  flg - PETSC_TRUE to use the Amat, PETSC_FALSE to use the Pmat (default is false)

   Notes:
   For the common case in which the linear system matrix and the matrix used to construct the
   preconditioner are identical, this routine is does nothing.

   Level: intermediate

.seealso: PCSetUseAmat(), PCBJACOBI, PGMG, PCFIELDSPLIT, PCCOMPOSITE
@*/
PetscErrorCode  PCGetUseAmat(PC pc,PetscBool *flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  *flg = pc->useAmat;
  PetscFunctionReturn(0);
}

/*@
   PCCreate - Creates a preconditioner context.

   Collective

   Input Parameter:
.  comm - MPI communicator

   Output Parameter:
.  pc - location to put the preconditioner context

   Notes:
   The default preconditioner for sparse matrices is PCILU or PCICC with 0 fill on one process and block Jacobi with PCILU or PCICC
   in parallel. For dense matrices it is always PCNONE.

   Level: developer

.seealso: PCSetUp(), PCApply(), PCDestroy()
@*/
PetscErrorCode  PCCreate(MPI_Comm comm,PC *newpc)
{
  PC             pc;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidPointer(newpc,1);
  *newpc = 0;
  ierr = PCInitializePackage();CHKERRQ(ierr);

  ierr = PetscHeaderCreate(pc,PC_CLASSID,"PC","Preconditioner","PC",comm,PCDestroy,PCView);CHKERRQ(ierr);

  pc->mat                  = 0;
  pc->pmat                 = 0;
  pc->setupcalled          = 0;
  pc->setfromoptionscalled = 0;
  pc->data                 = 0;
  pc->diagonalscale        = PETSC_FALSE;
  pc->diagonalscaleleft    = 0;
  pc->diagonalscaleright   = 0;

  pc->modifysubmatrices  = 0;
  pc->modifysubmatricesP = 0;

  *newpc = pc;
  PetscFunctionReturn(0);

}

/* -------------------------------------------------------------------------------*/

/*@
   PCApply - Applies the preconditioner to a vector.

   Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  x - input vector

   Output Parameter:
.  y - output vector

   Level: developer

.seealso: PCApplyTranspose(), PCApplyBAorAB()
@*/
PetscErrorCode  PCApply(PC pc,Vec x,Vec y)
{
  PetscErrorCode ierr;
  PetscInt       m,n,mv,nv;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  if (x == y) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  if (pc->erroriffailure) {ierr = VecValidValues(x,2,PETSC_TRUE);CHKERRQ(ierr);}
  /* use pmat to check vector sizes since for KSPLQR the pmat may be of a different size than mat */
  ierr = MatGetLocalSize(pc->pmat,&m,&n);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x,&nv);CHKERRQ(ierr);
  ierr = VecGetLocalSize(y,&mv);CHKERRQ(ierr);
  if (mv != m) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Preconditioner number of local rows %D does not equal resulting vector number of rows %D",m,mv);
  if (nv != n) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Preconditioner number of local columns %D does not equal resulting vector number of rows %D",n,nv);
  ierr = VecSetErrorIfLocked(y,3);CHKERRQ(ierr);

  ierr = PCSetUp(pc);CHKERRQ(ierr);
  if (!pc->ops->apply) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"PC does not have apply");
  ierr = VecLockReadPush(x);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(PC_Apply,pc,x,y,0);CHKERRQ(ierr);
  ierr = (*pc->ops->apply)(pc,x,y);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(PC_Apply,pc,x,y,0);CHKERRQ(ierr);
  if (pc->erroriffailure) {ierr = VecValidValues(y,3,PETSC_FALSE);CHKERRQ(ierr);}
  ierr = VecLockReadPop(x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PCApplySymmetricLeft - Applies the left part of a symmetric preconditioner to a vector.

   Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  x - input vector

   Output Parameter:
.  y - output vector

   Notes:
   Currently, this routine is implemented only for PCICC and PCJACOBI preconditioners.

   Level: developer

.seealso: PCApply(), PCApplySymmetricRight()
@*/
PetscErrorCode  PCApplySymmetricLeft(PC pc,Vec x,Vec y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  if (x == y) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  if (pc->erroriffailure) {ierr = VecValidValues(x,2,PETSC_TRUE);CHKERRQ(ierr);}
  ierr = PCSetUp(pc);CHKERRQ(ierr);
  if (!pc->ops->applysymmetricleft) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"PC does not have left symmetric apply");
  ierr = VecLockReadPush(x);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(PC_ApplySymmetricLeft,pc,x,y,0);CHKERRQ(ierr);
  ierr = (*pc->ops->applysymmetricleft)(pc,x,y);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(PC_ApplySymmetricLeft,pc,x,y,0);CHKERRQ(ierr);
  ierr = VecLockReadPop(x);CHKERRQ(ierr);
  if (pc->erroriffailure) {ierr = VecValidValues(y,3,PETSC_FALSE);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

/*@
   PCApplySymmetricRight - Applies the right part of a symmetric preconditioner to a vector.

   Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  x - input vector

   Output Parameter:
.  y - output vector

   Level: developer

   Notes:
   Currently, this routine is implemented only for PCICC and PCJACOBI preconditioners.

.seealso: PCApply(), PCApplySymmetricLeft()
@*/
PetscErrorCode  PCApplySymmetricRight(PC pc,Vec x,Vec y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  if (x == y) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  if (pc->erroriffailure) {ierr = VecValidValues(x,2,PETSC_TRUE);CHKERRQ(ierr);}
  ierr = PCSetUp(pc);CHKERRQ(ierr);
  if (!pc->ops->applysymmetricright) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"PC does not have left symmetric apply");
  ierr = VecLockReadPush(x);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(PC_ApplySymmetricRight,pc,x,y,0);CHKERRQ(ierr);
  ierr = (*pc->ops->applysymmetricright)(pc,x,y);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(PC_ApplySymmetricRight,pc,x,y,0);CHKERRQ(ierr);
  ierr = VecLockReadPop(x);CHKERRQ(ierr);
  if (pc->erroriffailure) {ierr = VecValidValues(y,3,PETSC_FALSE);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

/*@
   PCApplyTranspose - Applies the transpose of preconditioner to a vector.

   Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  x - input vector

   Output Parameter:
.  y - output vector

   Notes:
    For complex numbers this applies the non-Hermitian transpose.

   Developer Notes:
    We need to implement a PCApplyHermitianTranspose()

   Level: developer

.seealso: PCApply(), PCApplyBAorAB(), PCApplyBAorABTranspose(), PCApplyTransposeExists()
@*/
PetscErrorCode  PCApplyTranspose(PC pc,Vec x,Vec y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  if (x == y) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  if (pc->erroriffailure) {ierr = VecValidValues(x,2,PETSC_TRUE);CHKERRQ(ierr);}
  ierr = PCSetUp(pc);CHKERRQ(ierr);
  if (!pc->ops->applytranspose) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"PC does not have apply transpose");
  ierr = VecLockReadPush(x);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(PC_Apply,pc,x,y,0);CHKERRQ(ierr);
  ierr = (*pc->ops->applytranspose)(pc,x,y);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(PC_Apply,pc,x,y,0);CHKERRQ(ierr);
  ierr = VecLockReadPop(x);CHKERRQ(ierr);
  if (pc->erroriffailure) {ierr = VecValidValues(y,3,PETSC_FALSE);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

/*@
   PCApplyTransposeExists - Test whether the preconditioner has a transpose apply operation

   Collective on PC

   Input Parameters:
.  pc - the preconditioner context

   Output Parameter:
.  flg - PETSC_TRUE if a transpose operation is defined

   Level: developer

.seealso: PCApplyTranspose()
@*/
PetscErrorCode  PCApplyTransposeExists(PC pc,PetscBool  *flg)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidBoolPointer(flg,2);
  if (pc->ops->applytranspose) *flg = PETSC_TRUE;
  else *flg = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@
   PCApplyBAorAB - Applies the preconditioner and operator to a vector. y = B*A*x or y = A*B*x.

   Collective on PC

   Input Parameters:
+  pc - the preconditioner context
.  side - indicates the preconditioner side, one of PC_LEFT, PC_RIGHT, or PC_SYMMETRIC
.  x - input vector
-  work - work vector

   Output Parameter:
.  y - output vector

   Level: developer

   Notes:
    If the PC has had PCSetDiagonalScale() set then D M A D^{-1} for left preconditioning or  D A M D^{-1} is actually applied. Note that the
   specific KSPSolve() method must also be written to handle the post-solve "correction" for the diagonal scaling.

.seealso: PCApply(), PCApplyTranspose(), PCApplyBAorABTranspose()
@*/
PetscErrorCode  PCApplyBAorAB(PC pc,PCSide side,Vec x,Vec y,Vec work)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidLogicalCollectiveEnum(pc,side,2);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  PetscValidHeaderSpecific(y,VEC_CLASSID,4);
  PetscValidHeaderSpecific(work,VEC_CLASSID,5);
  PetscCheckSameComm(pc,1,x,3);
  PetscCheckSameComm(pc,1,y,4);
  PetscCheckSameComm(pc,1,work,5);
  if (x == y) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  if (side != PC_LEFT && side != PC_SYMMETRIC && side != PC_RIGHT) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_OUTOFRANGE,"Side must be right, left, or symmetric");
  if (pc->diagonalscale && side == PC_SYMMETRIC) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"Cannot include diagonal scaling with symmetric preconditioner application");
  if (pc->erroriffailure) {ierr = VecValidValues(x,3,PETSC_TRUE);CHKERRQ(ierr);}

  ierr = PCSetUp(pc);CHKERRQ(ierr);
  if (pc->diagonalscale) {
    if (pc->ops->applyBA) {
      Vec work2; /* this is expensive, but to fix requires a second work vector argument to PCApplyBAorAB() */
      ierr = VecDuplicate(x,&work2);CHKERRQ(ierr);
      ierr = PCDiagonalScaleRight(pc,x,work2);CHKERRQ(ierr);
      ierr = (*pc->ops->applyBA)(pc,side,work2,y,work);CHKERRQ(ierr);
      ierr = PCDiagonalScaleLeft(pc,y,y);CHKERRQ(ierr);
      ierr = VecDestroy(&work2);CHKERRQ(ierr);
    } else if (side == PC_RIGHT) {
      ierr = PCDiagonalScaleRight(pc,x,y);CHKERRQ(ierr);
      ierr = PCApply(pc,y,work);CHKERRQ(ierr);
      ierr = MatMult(pc->mat,work,y);CHKERRQ(ierr);
      ierr = PCDiagonalScaleLeft(pc,y,y);CHKERRQ(ierr);
    } else if (side == PC_LEFT) {
      ierr = PCDiagonalScaleRight(pc,x,y);CHKERRQ(ierr);
      ierr = MatMult(pc->mat,y,work);CHKERRQ(ierr);
      ierr = PCApply(pc,work,y);CHKERRQ(ierr);
      ierr = PCDiagonalScaleLeft(pc,y,y);CHKERRQ(ierr);
    } else if (side == PC_SYMMETRIC) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"Cannot provide diagonal scaling with symmetric application of preconditioner");
  } else {
    if (pc->ops->applyBA) {
      ierr = (*pc->ops->applyBA)(pc,side,x,y,work);CHKERRQ(ierr);
    } else if (side == PC_RIGHT) {
      ierr = PCApply(pc,x,work);CHKERRQ(ierr);
      ierr = MatMult(pc->mat,work,y);CHKERRQ(ierr);
    } else if (side == PC_LEFT) {
      ierr = MatMult(pc->mat,x,work);CHKERRQ(ierr);
      ierr = PCApply(pc,work,y);CHKERRQ(ierr);
    } else if (side == PC_SYMMETRIC) {
      /* There's an extra copy here; maybe should provide 2 work vectors instead? */
      ierr = PCApplySymmetricRight(pc,x,work);CHKERRQ(ierr);
      ierr = MatMult(pc->mat,work,y);CHKERRQ(ierr);
      ierr = VecCopy(y,work);CHKERRQ(ierr);
      ierr = PCApplySymmetricLeft(pc,work,y);CHKERRQ(ierr);
    }
  }
  if (pc->erroriffailure) {ierr = VecValidValues(y,4,PETSC_FALSE);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

/*@
   PCApplyBAorABTranspose - Applies the transpose of the preconditioner
   and operator to a vector. That is, applies tr(B) * tr(A) with left preconditioning,
   NOT tr(B*A) = tr(A)*tr(B).

   Collective on PC

   Input Parameters:
+  pc - the preconditioner context
.  side - indicates the preconditioner side, one of PC_LEFT, PC_RIGHT, or PC_SYMMETRIC
.  x - input vector
-  work - work vector

   Output Parameter:
.  y - output vector


   Notes:
    this routine is used internally so that the same Krylov code can be used to solve A x = b and A' x = b, with a preconditioner
      defined by B'. This is why this has the funny form that it computes tr(B) * tr(A)

    Level: developer

.seealso: PCApply(), PCApplyTranspose(), PCApplyBAorAB()
@*/
PetscErrorCode  PCApplyBAorABTranspose(PC pc,PCSide side,Vec x,Vec y,Vec work)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,3);
  PetscValidHeaderSpecific(y,VEC_CLASSID,4);
  PetscValidHeaderSpecific(work,VEC_CLASSID,5);
  if (x == y) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_IDN,"x and y must be different vectors");
  if (pc->erroriffailure) {ierr = VecValidValues(x,3,PETSC_TRUE);CHKERRQ(ierr);}
  if (pc->ops->applyBAtranspose) {
    ierr = (*pc->ops->applyBAtranspose)(pc,side,x,y,work);CHKERRQ(ierr);
    if (pc->erroriffailure) {ierr = VecValidValues(y,4,PETSC_FALSE);CHKERRQ(ierr);}
    PetscFunctionReturn(0);
  }
  if (side != PC_LEFT && side != PC_RIGHT) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_OUTOFRANGE,"Side must be right or left");

  ierr = PCSetUp(pc);CHKERRQ(ierr);
  if (side == PC_RIGHT) {
    ierr = PCApplyTranspose(pc,x,work);CHKERRQ(ierr);
    ierr = MatMultTranspose(pc->mat,work,y);CHKERRQ(ierr);
  } else if (side == PC_LEFT) {
    ierr = MatMultTranspose(pc->mat,x,work);CHKERRQ(ierr);
    ierr = PCApplyTranspose(pc,work,y);CHKERRQ(ierr);
  }
  /* add support for PC_SYMMETRIC */
  if (pc->erroriffailure) {ierr = VecValidValues(y,4,PETSC_FALSE);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------------*/

/*@
   PCApplyRichardsonExists - Determines whether a particular preconditioner has a
   built-in fast application of Richardson's method.

   Not Collective

   Input Parameter:
.  pc - the preconditioner

   Output Parameter:
.  exists - PETSC_TRUE or PETSC_FALSE

   Level: developer

.seealso: PCApplyRichardson()
@*/
PetscErrorCode  PCApplyRichardsonExists(PC pc,PetscBool  *exists)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidPointer(exists,2);
  if (pc->ops->applyrichardson) *exists = PETSC_TRUE;
  else *exists = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@
   PCApplyRichardson - Applies several steps of Richardson iteration with
   the particular preconditioner. This routine is usually used by the
   Krylov solvers and not the application code directly.

   Collective on PC

   Input Parameters:
+  pc  - the preconditioner context
.  b   - the right hand side
.  w   - one work vector
.  rtol - relative decrease in residual norm convergence criteria
.  abstol - absolute residual norm convergence criteria
.  dtol - divergence residual norm increase criteria
.  its - the number of iterations to apply.
-  guesszero - if the input x contains nonzero initial guess

   Output Parameter:
+  outits - number of iterations actually used (for SOR this always equals its)
.  reason - the reason the apply terminated
-  y - the solution (also contains initial guess if guesszero is PETSC_FALSE

   Notes:
   Most preconditioners do not support this function. Use the command
   PCApplyRichardsonExists() to determine if one does.

   Except for the multigrid PC this routine ignores the convergence tolerances
   and always runs for the number of iterations

   Level: developer

.seealso: PCApplyRichardsonExists()
@*/
PetscErrorCode  PCApplyRichardson(PC pc,Vec b,Vec y,Vec w,PetscReal rtol,PetscReal abstol, PetscReal dtol,PetscInt its,PetscBool guesszero,PetscInt *outits,PCRichardsonConvergedReason *reason)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(b,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  PetscValidHeaderSpecific(w,VEC_CLASSID,4);
  if (b == y) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_IDN,"b and y must be different vectors");
  ierr = PCSetUp(pc);CHKERRQ(ierr);
  if (!pc->ops->applyrichardson) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"PC does not have apply richardson");
  ierr = (*pc->ops->applyrichardson)(pc,b,y,w,rtol,abstol,dtol,its,guesszero,outits,reason);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PCGetFailedReason - Gets the reason a PCSetUp() failed or 0 if it did not fail

   Logically Collective on PC

   Input Parameter:
.  pc - the preconditioner context

   Output Parameter:
.  reason - the reason it failed, currently only -1

   Level: advanced

.seealso: PCCreate(), PCApply(), PCDestroy()
@*/
PetscErrorCode PCGetFailedReason(PC pc,PCFailedReason *reason)
{
  PetscFunctionBegin;
  if (pc->setupcalled < 0) *reason = (PCFailedReason)pc->setupcalled;
  else *reason = pc->failedreason;
  PetscFunctionReturn(0);
}


/*
      a setupcall of 0 indicates never setup,
                     1 indicates has been previously setup
                    -1 indicates a PCSetUp() was attempted and failed
*/
/*@
   PCSetUp - Prepares for the use of a preconditioner.

   Collective on PC

   Input Parameter:
.  pc - the preconditioner context

   Level: developer

.seealso: PCCreate(), PCApply(), PCDestroy()
@*/
PetscErrorCode  PCSetUp(PC pc)
{
  PetscErrorCode   ierr;
  const char       *def;
  PetscObjectState matstate, matnonzerostate;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  if (!pc->mat) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_ARG_WRONGSTATE,"Matrix must be set first");

  if (pc->setupcalled && pc->reusepreconditioner) {
    ierr = PetscInfo(pc,"Leaving PC with identical preconditioner since reuse preconditioner is set\n");CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ierr = PetscObjectStateGet((PetscObject)pc->pmat,&matstate);CHKERRQ(ierr);
  ierr = MatGetNonzeroState(pc->pmat,&matnonzerostate);CHKERRQ(ierr);
  if (!pc->setupcalled) {
    ierr     = PetscInfo(pc,"Setting up PC for first time\n");CHKERRQ(ierr);
    pc->flag = DIFFERENT_NONZERO_PATTERN;
  } else if (matstate == pc->matstate) {
    ierr = PetscInfo(pc,"Leaving PC with identical preconditioner since operator is unchanged\n");CHKERRQ(ierr);
    PetscFunctionReturn(0);
  } else {
    if (matnonzerostate > pc->matnonzerostate) {
       ierr = PetscInfo(pc,"Setting up PC with different nonzero pattern\n");CHKERRQ(ierr);
       pc->flag = DIFFERENT_NONZERO_PATTERN;
    } else {
      ierr = PetscInfo(pc,"Setting up PC with same nonzero pattern\n");CHKERRQ(ierr);
      pc->flag = SAME_NONZERO_PATTERN;
    }
  }
  pc->matstate        = matstate;
  pc->matnonzerostate = matnonzerostate;

  if (!((PetscObject)pc)->type_name) {
    ierr = PCGetDefaultType_Private(pc,&def);CHKERRQ(ierr);
    ierr = PCSetType(pc,def);CHKERRQ(ierr);
  }

  ierr = MatSetErrorIfFailure(pc->pmat,pc->erroriffailure);CHKERRQ(ierr);
  ierr = MatSetErrorIfFailure(pc->mat,pc->erroriffailure);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(PC_SetUp,pc,0,0,0);CHKERRQ(ierr);
  if (pc->ops->setup) {
    ierr = (*pc->ops->setup)(pc);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(PC_SetUp,pc,0,0,0);CHKERRQ(ierr);
  if (!pc->setupcalled) pc->setupcalled = 1;
  PetscFunctionReturn(0);
}

/*@
   PCSetUpOnBlocks - Sets up the preconditioner for each block in
   the block Jacobi, block Gauss-Seidel, and overlapping Schwarz
   methods.

   Collective on PC

   Input Parameters:
.  pc - the preconditioner context

   Level: developer

.seealso: PCCreate(), PCApply(), PCDestroy(), PCSetUp()
@*/
PetscErrorCode  PCSetUpOnBlocks(PC pc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  if (!pc->ops->setuponblocks) PetscFunctionReturn(0);
  ierr = PetscLogEventBegin(PC_SetUpOnBlocks,pc,0,0,0);CHKERRQ(ierr);
  ierr = (*pc->ops->setuponblocks)(pc);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(PC_SetUpOnBlocks,pc,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PCSetModifySubMatrices - Sets a user-defined routine for modifying the
   submatrices that arise within certain subdomain-based preconditioners.
   The basic submatrices are extracted from the preconditioner matrix as
   usual; the user can then alter these (for example, to set different boundary
   conditions for each submatrix) before they are used for the local solves.

   Logically Collective on PC

   Input Parameters:
+  pc - the preconditioner context
.  func - routine for modifying the submatrices
-  ctx - optional user-defined context (may be null)

   Calling sequence of func:
$     func (PC pc,PetscInt nsub,IS *row,IS *col,Mat *submat,void *ctx);

+  row - an array of index sets that contain the global row numbers
         that comprise each local submatrix
.  col - an array of index sets that contain the global column numbers
         that comprise each local submatrix
.  submat - array of local submatrices
-  ctx - optional user-defined context for private data for the
         user-defined func routine (may be null)

   Notes:
   PCSetModifySubMatrices() MUST be called before KSPSetUp() and
   KSPSolve().

   A routine set by PCSetModifySubMatrices() is currently called within
   the block Jacobi (PCBJACOBI) and additive Schwarz (PCASM)
   preconditioners.  All other preconditioners ignore this routine.

   Level: advanced

.seealso: PCModifySubMatrices()
@*/
PetscErrorCode  PCSetModifySubMatrices(PC pc,PetscErrorCode (*func)(PC,PetscInt,const IS[],const IS[],Mat[],void*),void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  pc->modifysubmatrices  = func;
  pc->modifysubmatricesP = ctx;
  PetscFunctionReturn(0);
}

/*@C
   PCModifySubMatrices - Calls an optional user-defined routine within
   certain preconditioners if one has been set with PCSetModifySubMatrices().

   Collective on PC

   Input Parameters:
+  pc - the preconditioner context
.  nsub - the number of local submatrices
.  row - an array of index sets that contain the global row numbers
         that comprise each local submatrix
.  col - an array of index sets that contain the global column numbers
         that comprise each local submatrix
.  submat - array of local submatrices
-  ctx - optional user-defined context for private data for the
         user-defined routine (may be null)

   Output Parameter:
.  submat - array of local submatrices (the entries of which may
            have been modified)

   Notes:
   The user should NOT generally call this routine, as it will
   automatically be called within certain preconditioners (currently
   block Jacobi, additive Schwarz) if set.

   The basic submatrices are extracted from the preconditioner matrix
   as usual; the user can then alter these (for example, to set different
   boundary conditions for each submatrix) before they are used for the
   local solves.

   Level: developer

.seealso: PCSetModifySubMatrices()
@*/
PetscErrorCode  PCModifySubMatrices(PC pc,PetscInt nsub,const IS row[],const IS col[],Mat submat[],void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  if (!pc->modifysubmatrices) PetscFunctionReturn(0);
  ierr = PetscLogEventBegin(PC_ModifySubMatrices,pc,0,0,0);CHKERRQ(ierr);
  ierr = (*pc->modifysubmatrices)(pc,nsub,row,col,submat,ctx);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(PC_ModifySubMatrices,pc,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PCSetOperators - Sets the matrix associated with the linear system and
   a (possibly) different one associated with the preconditioner.

   Logically Collective on PC

   Input Parameters:
+  pc - the preconditioner context
.  Amat - the matrix that defines the linear system
-  Pmat - the matrix to be used in constructing the preconditioner, usually the same as Amat.

   Notes:
    Passing a NULL for Amat or Pmat removes the matrix that is currently used.

    If you wish to replace either Amat or Pmat but leave the other one untouched then
    first call KSPGetOperators() to get the one you wish to keep, call PetscObjectReference()
    on it and then pass it back in in your call to KSPSetOperators().

   More Notes about Repeated Solution of Linear Systems:
   PETSc does NOT reset the matrix entries of either Amat or Pmat
   to zero after a linear solve; the user is completely responsible for
   matrix assembly.  See the routine MatZeroEntries() if desiring to
   zero all elements of a matrix.

   Level: intermediate

.seealso: PCGetOperators(), MatZeroEntries()
 @*/
PetscErrorCode  PCSetOperators(PC pc,Mat Amat,Mat Pmat)
{
  PetscErrorCode   ierr;
  PetscInt         m1,n1,m2,n2;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  if (Amat) PetscValidHeaderSpecific(Amat,MAT_CLASSID,2);
  if (Pmat) PetscValidHeaderSpecific(Pmat,MAT_CLASSID,3);
  if (Amat) PetscCheckSameComm(pc,1,Amat,2);
  if (Pmat) PetscCheckSameComm(pc,1,Pmat,3);
  if (pc->setupcalled && pc->mat && pc->pmat && Amat && Pmat) {
    ierr = MatGetLocalSize(Amat,&m1,&n1);CHKERRQ(ierr);
    ierr = MatGetLocalSize(pc->mat,&m2,&n2);CHKERRQ(ierr);
    if (m1 != m2 || n1 != n2) SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Cannot change local size of Amat after use old sizes %D %D new sizes %D %D",m2,n2,m1,n1);
    ierr = MatGetLocalSize(Pmat,&m1,&n1);CHKERRQ(ierr);
    ierr = MatGetLocalSize(pc->pmat,&m2,&n2);CHKERRQ(ierr);
    if (m1 != m2 || n1 != n2) SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Cannot change local size of Pmat after use old sizes %D %D new sizes %D %D",m2,n2,m1,n1);
  }

  if (Pmat != pc->pmat) {
    /* changing the operator that defines the preconditioner thus reneed to clear current states so new preconditioner is built */
    pc->matnonzerostate = -1;
    pc->matstate        = -1;
  }

  /* reference first in case the matrices are the same */
  if (Amat) {ierr = PetscObjectReference((PetscObject)Amat);CHKERRQ(ierr);}
  ierr = MatDestroy(&pc->mat);CHKERRQ(ierr);
  if (Pmat) {ierr = PetscObjectReference((PetscObject)Pmat);CHKERRQ(ierr);}
  ierr     = MatDestroy(&pc->pmat);CHKERRQ(ierr);
  pc->mat  = Amat;
  pc->pmat = Pmat;
  PetscFunctionReturn(0);
}

/*@
   PCSetReusePreconditioner - reuse the current preconditioner even if the operator in the preconditioner has changed.

   Logically Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  flag - PETSC_TRUE do not compute a new preconditioner, PETSC_FALSE do compute a new preconditioner

    Level: intermediate

.seealso: PCGetOperators(), MatZeroEntries(), PCGetReusePreconditioner(), KSPSetReusePreconditioner()
 @*/
PetscErrorCode  PCSetReusePreconditioner(PC pc,PetscBool flag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  pc->reusepreconditioner = flag;
  PetscFunctionReturn(0);
}

/*@
   PCGetReusePreconditioner - Determines if the PC reuses the current preconditioner even if the operator in the preconditioner has changed.

   Not Collective

   Input Parameter:
.  pc - the preconditioner context

   Output Parameter:
.  flag - PETSC_TRUE do not compute a new preconditioner, PETSC_FALSE do compute a new preconditioner

   Level: intermediate

.seealso: PCGetOperators(), MatZeroEntries(), PCSetReusePreconditioner()
 @*/
PetscErrorCode  PCGetReusePreconditioner(PC pc,PetscBool *flag)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  *flag = pc->reusepreconditioner;
  PetscFunctionReturn(0);
}

/*@
   PCGetOperators - Gets the matrix associated with the linear system and
   possibly a different one associated with the preconditioner.

   Not collective, though parallel Mats are returned if the PC is parallel

   Input Parameter:
.  pc - the preconditioner context

   Output Parameters:
+  Amat - the matrix defining the linear system
-  Pmat - the matrix from which the preconditioner is constructed, usually the same as Amat.

   Level: intermediate

   Notes:
    Does not increase the reference count of the matrices, so you should not destroy them

   Alternative usage: If the operators have NOT been set with KSP/PCSetOperators() then the operators
      are created in PC and returned to the user. In this case, if both operators
      mat and pmat are requested, two DIFFERENT operators will be returned. If
      only one is requested both operators in the PC will be the same (i.e. as
      if one had called KSP/PCSetOperators() with the same argument for both Mats).
      The user must set the sizes of the returned matrices and their type etc just
      as if the user created them with MatCreate(). For example,

$         KSP/PCGetOperators(ksp/pc,&Amat,NULL); is equivalent to
$           set size, type, etc of Amat

$         MatCreate(comm,&mat);
$         KSP/PCSetOperators(ksp/pc,Amat,Amat);
$         PetscObjectDereference((PetscObject)mat);
$           set size, type, etc of Amat

     and

$         KSP/PCGetOperators(ksp/pc,&Amat,&Pmat); is equivalent to
$           set size, type, etc of Amat and Pmat

$         MatCreate(comm,&Amat);
$         MatCreate(comm,&Pmat);
$         KSP/PCSetOperators(ksp/pc,Amat,Pmat);
$         PetscObjectDereference((PetscObject)Amat);
$         PetscObjectDereference((PetscObject)Pmat);
$           set size, type, etc of Amat and Pmat

    The rational for this support is so that when creating a TS, SNES, or KSP the hierarchy
    of underlying objects (i.e. SNES, KSP, PC, Mat) and their livespans can be completely
    managed by the top most level object (i.e. the TS, SNES, or KSP). Another way to look
    at this is when you create a SNES you do not NEED to create a KSP and attach it to
    the SNES object (the SNES object manages it for you). Similarly when you create a KSP
    you do not need to attach a PC to it (the KSP object manages the PC object for you).
    Thus, why should YOU have to create the Mat and attach it to the SNES/KSP/PC, when
    it can be created for you?


.seealso: PCSetOperators(), KSPGetOperators(), KSPSetOperators(), PCGetOperatorsSet()
@*/
PetscErrorCode  PCGetOperators(PC pc,Mat *Amat,Mat *Pmat)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  if (Amat) {
    if (!pc->mat) {
      if (pc->pmat && !Pmat) {  /* Apmat has been set, but user did not request it, so use for Amat */
        pc->mat = pc->pmat;
        ierr    = PetscObjectReference((PetscObject)pc->mat);CHKERRQ(ierr);
      } else {                  /* both Amat and Pmat are empty */
        ierr = MatCreate(PetscObjectComm((PetscObject)pc),&pc->mat);CHKERRQ(ierr);
        if (!Pmat) { /* user did NOT request Pmat, so make same as Amat */
          pc->pmat = pc->mat;
          ierr     = PetscObjectReference((PetscObject)pc->pmat);CHKERRQ(ierr);
        }
      }
    }
    *Amat = pc->mat;
  }
  if (Pmat) {
    if (!pc->pmat) {
      if (pc->mat && !Amat) {    /* Amat has been set but was not requested, so use for pmat */
        pc->pmat = pc->mat;
        ierr     = PetscObjectReference((PetscObject)pc->pmat);CHKERRQ(ierr);
      } else {
        ierr = MatCreate(PetscObjectComm((PetscObject)pc),&pc->pmat);CHKERRQ(ierr);
        if (!Amat) { /* user did NOT request Amat, so make same as Pmat */
          pc->mat = pc->pmat;
          ierr    = PetscObjectReference((PetscObject)pc->mat);CHKERRQ(ierr);
        }
      }
    }
    *Pmat = pc->pmat;
  }
  PetscFunctionReturn(0);
}

/*@C
   PCGetOperatorsSet - Determines if the matrix associated with the linear system and
   possibly a different one associated with the preconditioner have been set in the PC.

   Not collective, though the results on all processes should be the same

   Input Parameter:
.  pc - the preconditioner context

   Output Parameters:
+  mat - the matrix associated with the linear system was set
-  pmat - matrix associated with the preconditioner was set, usually the same

   Level: intermediate

.seealso: PCSetOperators(), KSPGetOperators(), KSPSetOperators(), PCGetOperators()
@*/
PetscErrorCode  PCGetOperatorsSet(PC pc,PetscBool  *mat,PetscBool  *pmat)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  if (mat) *mat = (pc->mat)  ? PETSC_TRUE : PETSC_FALSE;
  if (pmat) *pmat = (pc->pmat) ? PETSC_TRUE : PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*@
   PCFactorGetMatrix - Gets the factored matrix from the
   preconditioner context.  This routine is valid only for the LU,
   incomplete LU, Cholesky, and incomplete Cholesky methods.

   Not Collective on PC though Mat is parallel if PC is parallel

   Input Parameters:
.  pc - the preconditioner context

   Output parameters:
.  mat - the factored matrix

   Level: advanced

   Notes:
    Does not increase the reference count for the matrix so DO NOT destroy it

@*/
PetscErrorCode  PCFactorGetMatrix(PC pc,Mat *mat)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidPointer(mat,2);
  if (pc->ops->getfactoredmatrix) {
    ierr = (*pc->ops->getfactoredmatrix)(pc,mat);CHKERRQ(ierr);
  } else SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"PC type does not support getting factor matrix");
  PetscFunctionReturn(0);
}

/*@C
   PCSetOptionsPrefix - Sets the prefix used for searching for all
   PC options in the database.

   Logically Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  prefix - the prefix string to prepend to all PC option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: PCAppendOptionsPrefix(), PCGetOptionsPrefix()
@*/
PetscErrorCode  PCSetOptionsPrefix(PC pc,const char prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  ierr = PetscObjectSetOptionsPrefix((PetscObject)pc,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PCAppendOptionsPrefix - Appends to the prefix used for searching for all
   PC options in the database.

   Logically Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  prefix - the prefix string to prepend to all PC option requests

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

   Level: advanced

.seealso: PCSetOptionsPrefix(), PCGetOptionsPrefix()
@*/
PetscErrorCode  PCAppendOptionsPrefix(PC pc,const char prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  ierr = PetscObjectAppendOptionsPrefix((PetscObject)pc,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PCGetOptionsPrefix - Gets the prefix used for searching for all
   PC options in the database.

   Not Collective

   Input Parameters:
.  pc - the preconditioner context

   Output Parameters:
.  prefix - pointer to the prefix string used, is returned

   Notes:
    On the fortran side, the user should pass in a string 'prifix' of
   sufficient length to hold the prefix.

   Level: advanced

.seealso: PCSetOptionsPrefix(), PCAppendOptionsPrefix()
@*/
PetscErrorCode  PCGetOptionsPrefix(PC pc,const char *prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidPointer(prefix,2);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)pc,prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Indicates the right hand side will be changed by KSPSolve(), this occurs for a few
  preconditioners including BDDC and Eisentat that transform the equations before applying
  the Krylov methods
*/
PETSC_INTERN PetscErrorCode  PCPreSolveChangeRHS(PC pc,PetscBool *change)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidPointer(change,2);
  *change = PETSC_FALSE;
  ierr = PetscTryMethod(pc,"PCPreSolveChangeRHS_C",(PC,PetscBool*),(pc,change));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PCPreSolve - Optional pre-solve phase, intended for any
   preconditioner-specific actions that must be performed before
   the iterative solve itself.

   Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  ksp - the Krylov subspace context

   Level: developer

   Sample of Usage:
.vb
    PCPreSolve(pc,ksp);
    KSPSolve(ksp,b,x);
    PCPostSolve(pc,ksp);
.ve

   Notes:
   The pre-solve phase is distinct from the PCSetUp() phase.

   KSPSolve() calls this directly, so is rarely called by the user.

.seealso: PCPostSolve()
@*/
PetscErrorCode  PCPreSolve(PC pc,KSP ksp)
{
  PetscErrorCode ierr;
  Vec            x,rhs;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,2);
  pc->presolvedone++;
  if (pc->presolvedone > 2) SETERRQ(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"Cannot embed PCPreSolve() more than twice");
  ierr = KSPGetSolution(ksp,&x);CHKERRQ(ierr);
  ierr = KSPGetRhs(ksp,&rhs);CHKERRQ(ierr);

  if (pc->ops->presolve) {
    ierr = (*pc->ops->presolve)(pc,ksp,rhs,x);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PCPostSolve - Optional post-solve phase, intended for any
   preconditioner-specific actions that must be performed after
   the iterative solve itself.

   Collective on PC

   Input Parameters:
+  pc - the preconditioner context
-  ksp - the Krylov subspace context

   Sample of Usage:
.vb
    PCPreSolve(pc,ksp);
    KSPSolve(ksp,b,x);
    PCPostSolve(pc,ksp);
.ve

   Note:
   KSPSolve() calls this routine directly, so it is rarely called by the user.

   Level: developer

.seealso: PCPreSolve(), KSPSolve()
@*/
PetscErrorCode  PCPostSolve(PC pc,KSP ksp)
{
  PetscErrorCode ierr;
  Vec            x,rhs;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidHeaderSpecific(ksp,KSP_CLASSID,2);
  pc->presolvedone--;
  ierr = KSPGetSolution(ksp,&x);CHKERRQ(ierr);
  ierr = KSPGetRhs(ksp,&rhs);CHKERRQ(ierr);
  if (pc->ops->postsolve) {
    ierr =  (*pc->ops->postsolve)(pc,ksp,rhs,x);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
  PCLoad - Loads a PC that has been stored in binary  with PCView().

  Collective on PetscViewer

  Input Parameters:
+ newdm - the newly loaded PC, this needs to have been created with PCCreate() or
           some related function before a call to PCLoad().
- viewer - binary file viewer, obtained from PetscViewerBinaryOpen()

   Level: intermediate

  Notes:
   The type is determined by the data in the file, any type set into the PC before this call is ignored.

  Notes for advanced users:
  Most users should not need to know the details of the binary storage
  format, since PCLoad() and PCView() completely hide these details.
  But for anyone who's interested, the standard binary matrix storage
  format is
.vb
     has not yet been determined
.ve

.seealso: PetscViewerBinaryOpen(), PCView(), MatLoad(), VecLoad()
@*/
PetscErrorCode  PCLoad(PC newdm, PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isbinary;
  PetscInt       classid;
  char           type[256];

  PetscFunctionBegin;
  PetscValidHeaderSpecific(newdm,PC_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary);CHKERRQ(ierr);
  if (!isbinary) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Invalid viewer; open viewer with PetscViewerBinaryOpen()");

  ierr = PetscViewerBinaryRead(viewer,&classid,1,NULL,PETSC_INT);CHKERRQ(ierr);
  if (classid != PC_FILE_CLASSID) SETERRQ(PetscObjectComm((PetscObject)newdm),PETSC_ERR_ARG_WRONG,"Not PC next in file");
  ierr = PetscViewerBinaryRead(viewer,type,256,NULL,PETSC_CHAR);CHKERRQ(ierr);
  ierr = PCSetType(newdm, type);CHKERRQ(ierr);
  if (newdm->ops->load) {
    ierr = (*newdm->ops->load)(newdm,viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#include <petscdraw.h>
#if defined(PETSC_HAVE_SAWS)
#include <petscviewersaws.h>
#endif

/*@C
   PCViewFromOptions - View from Options

   Collective on PC

   Input Parameters:
+  A - the PC context
.  obj - Optional object
-  name - command line option

   Level: intermediate
.seealso:  PC, PCView, PetscObjectViewFromOptions(), PCCreate()
@*/
PetscErrorCode  PCViewFromOptions(PC A,PetscObject obj,const char name[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(A,PC_CLASSID,1);
  ierr = PetscObjectViewFromOptions((PetscObject)A,obj,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PCView - Prints the PC data structure.

   Collective on PC

   Input Parameters:
+  PC - the PC context
-  viewer - optional visualization context

   Note:
   The available visualization contexts include
+     PETSC_VIEWER_STDOUT_SELF - standard output (default)
-     PETSC_VIEWER_STDOUT_WORLD - synchronized standard
         output where only the first processor opens
         the file.  All other processors send their
         data to the first processor to print.

   The user can open an alternative visualization contexts with
   PetscViewerASCIIOpen() (output to a specified file).

   Level: developer

.seealso: KSPView(), PetscViewerASCIIOpen()
@*/
PetscErrorCode  PCView(PC pc,PetscViewer viewer)
{
  PCType         cstr;
  PetscErrorCode ierr;
  PetscBool      iascii,isstring,isbinary,isdraw;
#if defined(PETSC_HAVE_SAWS)
  PetscBool      issaws;
#endif

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)pc),&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(pc,1,viewer,2);

  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERSTRING,&isstring);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERBINARY,&isbinary);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDRAW,&isdraw);CHKERRQ(ierr);
#if defined(PETSC_HAVE_SAWS)
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERSAWS,&issaws);CHKERRQ(ierr);
#endif

  if (iascii) {
    ierr = PetscObjectPrintClassNamePrefixType((PetscObject)pc,viewer);CHKERRQ(ierr);
    if (!pc->setupcalled) {
      ierr = PetscViewerASCIIPrintf(viewer,"  PC has not been set up so information may be incomplete\n");CHKERRQ(ierr);
    }
    if (pc->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*pc->ops->view)(pc,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    if (pc->mat) {
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
      if (pc->pmat == pc->mat) {
        ierr = PetscViewerASCIIPrintf(viewer,"  linear system matrix = precond matrix:\n");CHKERRQ(ierr);
        ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
        ierr = MatView(pc->mat,viewer);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
      } else {
        if (pc->pmat) {
          ierr = PetscViewerASCIIPrintf(viewer,"  linear system matrix followed by preconditioner matrix:\n");CHKERRQ(ierr);
        } else {
          ierr = PetscViewerASCIIPrintf(viewer,"  linear system matrix:\n");CHKERRQ(ierr);
        }
        ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
        ierr = MatView(pc->mat,viewer);CHKERRQ(ierr);
        if (pc->pmat) {ierr = MatView(pc->pmat,viewer);CHKERRQ(ierr);}
        ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
      }
      ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
    }
  } else if (isstring) {
    ierr = PCGetType(pc,&cstr);CHKERRQ(ierr);
    ierr = PetscViewerStringSPrintf(viewer," PCType: %-7.7s",cstr);CHKERRQ(ierr);
    if (pc->ops->view) {ierr = (*pc->ops->view)(pc,viewer);CHKERRQ(ierr);}
    if (pc->mat) {ierr = MatView(pc->mat,viewer);CHKERRQ(ierr);}
    if (pc->pmat && pc->pmat != pc->mat) {ierr = MatView(pc->pmat,viewer);CHKERRQ(ierr);}
  } else if (isbinary) {
    PetscInt    classid = PC_FILE_CLASSID;
    MPI_Comm    comm;
    PetscMPIInt rank;
    char        type[256];

    ierr = PetscObjectGetComm((PetscObject)pc,&comm);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    if (!rank) {
      ierr = PetscViewerBinaryWrite(viewer,&classid,1,PETSC_INT);CHKERRQ(ierr);
      ierr = PetscStrncpy(type,((PetscObject)pc)->type_name,256);CHKERRQ(ierr);
      ierr = PetscViewerBinaryWrite(viewer,type,256,PETSC_CHAR);CHKERRQ(ierr);
    }
    if (pc->ops->view) {
      ierr = (*pc->ops->view)(pc,viewer);CHKERRQ(ierr);
    }
  } else if (isdraw) {
    PetscDraw draw;
    char      str[25];
    PetscReal x,y,bottom,h;
    PetscInt  n;

    ierr = PetscViewerDrawGetDraw(viewer,0,&draw);CHKERRQ(ierr);
    ierr = PetscDrawGetCurrentPoint(draw,&x,&y);CHKERRQ(ierr);
    if (pc->mat) {
      ierr = MatGetSize(pc->mat,&n,NULL);CHKERRQ(ierr);
      ierr = PetscSNPrintf(str,25,"PC: %s (%D)",((PetscObject)pc)->type_name,n);CHKERRQ(ierr);
    } else {
      ierr = PetscSNPrintf(str,25,"PC: %s",((PetscObject)pc)->type_name);CHKERRQ(ierr);
    }
    ierr   = PetscDrawStringBoxed(draw,x,y,PETSC_DRAW_RED,PETSC_DRAW_BLACK,str,NULL,&h);CHKERRQ(ierr);
    bottom = y - h;
    ierr   = PetscDrawPushCurrentPoint(draw,x,bottom);CHKERRQ(ierr);
    if (pc->ops->view) {
      ierr = (*pc->ops->view)(pc,viewer);CHKERRQ(ierr);
    }
    ierr = PetscDrawPopCurrentPoint(draw);CHKERRQ(ierr);
#if defined(PETSC_HAVE_SAWS)
  } else if (issaws) {
    PetscMPIInt rank;

    ierr = PetscObjectName((PetscObject)pc);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    if (!((PetscObject)pc)->amsmem && !rank) {
      ierr = PetscObjectViewSAWs((PetscObject)pc,viewer);CHKERRQ(ierr);
    }
    if (pc->mat) {ierr = MatView(pc->mat,viewer);CHKERRQ(ierr);}
    if (pc->pmat && pc->pmat != pc->mat) {ierr = MatView(pc->pmat,viewer);CHKERRQ(ierr);}
#endif
  }
  PetscFunctionReturn(0);
}

/*@C
  PCRegister -  Adds a method to the preconditioner package.

   Not collective

   Input Parameters:
+  name_solver - name of a new user-defined solver
-  routine_create - routine to create method context

   Notes:
   PCRegister() may be called multiple times to add several user-defined preconditioners.

   Sample usage:
.vb
   PCRegister("my_solver", MySolverCreate);
.ve

   Then, your solver can be chosen with the procedural interface via
$     PCSetType(pc,"my_solver")
   or at runtime via the option
$     -pc_type my_solver

   Level: advanced

.seealso: PCRegisterAll()
@*/
PetscErrorCode  PCRegister(const char sname[],PetscErrorCode (*function)(PC))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PCInitializePackage();CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&PCList,sname,function);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_PC(Mat A,Vec X,Vec Y)
{
  PC             pc;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,&pc);CHKERRQ(ierr);
  ierr = PCApply(pc,X,Y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
    PCComputeOperator - Computes the explicit preconditioned operator.

    Collective on PC

    Input Parameter:
+   pc - the preconditioner object
-   mattype - the matrix type to be used for the operator

    Output Parameter:
.   mat - the explict preconditioned operator

    Notes:
    This computation is done by applying the operators to columns of the identity matrix.
    This routine is costly in general, and is recommended for use only with relatively small systems.
    Currently, this routine uses a dense matrix format when mattype == NULL

    Level: advanced

.seealso: KSPComputeOperator(), MatType

@*/
PetscErrorCode  PCComputeOperator(PC pc,MatType mattype,Mat *mat)
{
  PetscErrorCode ierr;
  PetscInt       N,M,m,n;
  Mat            A,Apc;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidPointer(mat,3);
  ierr = PCGetOperators(pc,&A,NULL);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&m,&n);CHKERRQ(ierr);
  ierr = MatGetSize(A,&M,&N);CHKERRQ(ierr);
  ierr = MatCreateShell(PetscObjectComm((PetscObject)pc),m,n,M,N,pc,&Apc);CHKERRQ(ierr);
  ierr = MatShellSetOperation(Apc,MATOP_MULT,(void (*)(void))MatMult_PC);CHKERRQ(ierr);
  ierr = MatComputeOperator(Apc,mattype,mat);CHKERRQ(ierr);
  ierr = MatDestroy(&Apc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PCSetCoordinates - sets the coordinates of all the nodes on the local process

   Collective on PC

   Input Parameters:
+  pc - the solver context
.  dim - the dimension of the coordinates 1, 2, or 3
.  nloc - the blocked size of the coordinates array
-  coords - the coordinates array

   Level: intermediate

   Notes:
   coords is an array of the dim coordinates for the nodes on
   the local processor, of size dim*nloc.
   If there are 108 equation on a processor
   for a displacement finite element discretization of elasticity (so
   that there are nloc = 36 = 108/3 nodes) then the array must have 108
   double precision values (ie, 3 * 36).  These x y z coordinates
   should be ordered for nodes 0 to N-1 like so: [ 0.x, 0.y, 0.z, 1.x,
   ... , N-1.z ].

.seealso: MatSetNearNullSpace()
@*/
PetscErrorCode PCSetCoordinates(PC pc, PetscInt dim, PetscInt nloc, PetscReal coords[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidLogicalCollectiveInt(pc,dim,2);
  ierr = PetscTryMethod(pc,"PCSetCoordinates_C",(PC,PetscInt,PetscInt,PetscReal*),(pc,dim,nloc,coords));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PCGetInterpolations - Gets interpolation matrices for all levels (except level 0)

   Logically Collective on PC

   Input Parameters:
+  pc - the precondition context

   Output Parameter:
-  num_levels - the number of levels
.  interpolations - the interpolation matrices (size of num_levels-1)

   Level: advanced

.keywords: MG, GAMG, BoomerAMG, multigrid, interpolation, level

.seealso: PCMGGetRestriction(), PCMGSetInterpolation(), PCMGGetInterpolation(), PCGetCoarseOperators()
@*/
PetscErrorCode PCGetInterpolations(PC pc,PetscInt *num_levels,Mat *interpolations[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidPointer(num_levels,2);
  PetscValidPointer(interpolations,3);
  ierr = PetscUseMethod(pc,"PCGetInterpolations_C",(PC,PetscInt*,Mat*[]),(pc,num_levels,interpolations));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PCGetCoarseOperators - Gets coarse operator matrices for all levels (except the finest level)

   Logically Collective on PC

   Input Parameters:
+  pc - the precondition context

   Output Parameter:
-  num_levels - the number of levels
.  coarseOperators - the coarse operator matrices (size of num_levels-1)

   Level: advanced

.keywords: MG, GAMG, BoomerAMG, get, multigrid, interpolation, level

.seealso: PCMGGetRestriction(), PCMGSetInterpolation(), PCMGGetRScale(), PCMGGetInterpolation(), PCGetInterpolations()
@*/
PetscErrorCode PCGetCoarseOperators(PC pc,PetscInt *num_levels,Mat *coarseOperators[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc,PC_CLASSID,1);
  PetscValidPointer(num_levels,2);
  PetscValidPointer(coarseOperators,3);
  ierr = PetscUseMethod(pc,"PCGetCoarseOperators_C",(PC,PetscInt*,Mat*[]),(pc,num_levels,coarseOperators));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
