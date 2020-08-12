/*
   This provides a thin wrapper around LMVM matrices in order to use their MatLMVMSolve 
   methods as preconditioner applications in KSP solves.
*/

#include <petsc/private/pcimpl.h>        /*I "petscpc.h" I*/
#include <petsc/private/matimpl.h>

typedef struct {
  Vec  xwork, ywork;
  IS   inactive;
  Mat  B;
  PetscBool allocated;
} PC_LMVM;

/*@
   PCLMVMSetMatLMVM - Replaces the LMVM matrix inside the preconditioner with 
   the one provided by the user.
   
   Input Parameters:
+  pc - An LMVM preconditioner
-  B  - An LMVM-type matrix (MATLDFP, MATLBFGS, MATLSR1, MATLBRDN, MATLMBRDN, MATLSBRDN)
  
   Level: intermediate
@*/
PetscErrorCode PCLMVMSetMatLMVM(PC pc, Mat B)
{
  PC_LMVM          *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode   ierr;
  PetscBool        same;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc, PC_CLASSID, 1);
  PetscValidHeaderSpecific(B, MAT_CLASSID, 2);
  ierr = PetscObjectTypeCompare((PetscObject)pc, PCLMVM, &same);CHKERRQ(ierr);
  if (!same) SETERRQ(PetscObjectComm((PetscObject)pc), PETSC_ERR_ARG_WRONG, "PC must be a PCLMVM type.");
  ierr = PetscObjectBaseTypeCompare((PetscObject)B, MATLMVM, &same);CHKERRQ(ierr);
  if (!same) SETERRQ(PetscObjectComm((PetscObject)pc), PETSC_ERR_ARG_WRONG, "Matrix must be an LMVM-type.");
  ierr = MatDestroy(&ctx->B);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)B);CHKERRQ(ierr);
  ctx->B = B;
  PetscFunctionReturn(0);
}

/*@
   PCLMVMGetMatLMVM - Returns a pointer to the underlying LMVM matrix.
   
   Input Parameters:
.  pc - An LMVM preconditioner

   Output Parameters:
.  B - LMVM matrix inside the preconditioner
  
   Level: intermediate
@*/
PetscErrorCode PCLMVMGetMatLMVM(PC pc, Mat *B)
{
  PC_LMVM          *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode   ierr;
  PetscBool        same;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc, PC_CLASSID, 1);
  ierr = PetscObjectTypeCompare((PetscObject)pc, PCLMVM, &same);CHKERRQ(ierr);
  if (!same) SETERRQ(PetscObjectComm((PetscObject)pc), PETSC_ERR_ARG_WRONG, "PC must be a PCLMVM type.");
  *B = ctx->B;
  PetscFunctionReturn(0);
}

/*@
   PCLMVMSetIS - Sets the index sets that reduce the PC application.
   
   Input Parameters:
+  pc - An LMVM preconditioner
-  inactive - Index set defining the variables removed from the problem
  
   Level: intermediate

.seealso:  MatLMVMUpdate()
@*/
PetscErrorCode PCLMVMSetIS(PC pc, IS inactive)
{
  PC_LMVM          *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode   ierr;
  PetscBool        same;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc, PC_CLASSID, 1);
  PetscValidHeaderSpecific(inactive, IS_CLASSID, 2);
  ierr = PetscObjectTypeCompare((PetscObject)pc, PCLMVM, &same);CHKERRQ(ierr);
  if (!same) SETERRQ(PetscObjectComm((PetscObject)pc), PETSC_ERR_ARG_WRONG, "PC must be a PCLMVM type.");
  ierr = PCLMVMClearIS(pc);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)inactive);CHKERRQ(ierr);
  ctx->inactive = inactive;
  PetscFunctionReturn(0);
}

/*@
   PCLMVMClearIS - Removes the inactive variable index set.
   
   Input Parameters:
.  pc - An LMVM preconditioner
  
   Level: intermediate

.seealso:  MatLMVMUpdate()
@*/
PetscErrorCode PCLMVMClearIS(PC pc)
{
  PC_LMVM          *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode   ierr;
  PetscBool        same;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(pc, PC_CLASSID, 1);
  ierr = PetscObjectTypeCompare((PetscObject)pc, PCLMVM, &same);CHKERRQ(ierr);
  if (!same) SETERRQ(PetscObjectComm((PetscObject)pc), PETSC_ERR_ARG_WRONG, "PC must be a PCLMVM type.");
  if (ctx->inactive) {
    ierr = ISDestroy(&ctx->inactive);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PCApply_LMVM(PC pc,Vec x,Vec y)
{
  PC_LMVM          *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode   ierr;
  Vec              xsub, ysub;

  PetscFunctionBegin;
  if (ctx->inactive) {
    ierr = VecZeroEntries(ctx->xwork);CHKERRQ(ierr);
    ierr = VecGetSubVector(ctx->xwork, ctx->inactive, &xsub);CHKERRQ(ierr);
    ierr = VecCopy(x, xsub);CHKERRQ(ierr);
    ierr = VecRestoreSubVector(ctx->xwork, ctx->inactive, &xsub);CHKERRQ(ierr);
  } else {
    ierr = VecCopy(x, ctx->xwork);CHKERRQ(ierr);
  }
  ierr = MatSolve(ctx->B, ctx->xwork, ctx->ywork);CHKERRQ(ierr);
  if (ctx->inactive) {
    ierr = VecGetSubVector(ctx->ywork, ctx->inactive, &ysub);CHKERRQ(ierr);
    ierr = VecCopy(ysub, y);CHKERRQ(ierr);
    ierr = VecRestoreSubVector(ctx->ywork, ctx->inactive, &ysub);CHKERRQ(ierr);
  } else {
    ierr = VecCopy(ctx->ywork, y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PCReset_LMVM(PC pc)
{
  PC_LMVM        *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  if (ctx->xwork) {
    ierr = VecDestroy(&ctx->xwork);CHKERRQ(ierr);
  }
  if (ctx->ywork) {
    ierr = VecDestroy(&ctx->ywork);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PCSetUp_LMVM(PC pc)
{
  PC_LMVM        *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode ierr;
  PetscInt       n, N;
  PetscBool      allocated;

  PetscFunctionBegin;
  ierr = MatLMVMIsAllocated(ctx->B, &allocated);CHKERRQ(ierr);
  if (!allocated) {
    ierr = MatCreateVecs(pc->mat, &ctx->xwork, &ctx->ywork);CHKERRQ(ierr);
    ierr = VecGetLocalSize(ctx->xwork, &n);CHKERRQ(ierr);
    ierr = VecGetSize(ctx->xwork, &N);CHKERRQ(ierr);
    ierr = MatSetSizes(ctx->B, n, n, N, N);CHKERRQ(ierr);
    ierr = MatLMVMAllocate(ctx->B, ctx->xwork, ctx->ywork);CHKERRQ(ierr);
  } else {
    ierr = MatCreateVecs(ctx->B, &ctx->xwork, &ctx->ywork);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PCView_LMVM(PC pc,PetscViewer viewer)
{
  PC_LMVM        *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode ierr;
  PetscBool      iascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  if (iascii && ctx->B->assembled) {
    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
    ierr = MatView(ctx->B, viewer);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode PCSetFromOptions_LMVM(PetscOptionItems* PetscOptionsObject, PC pc)
{
  PC_LMVM        *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = MatSetFromOptions(ctx->B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PCDestroy_LMVM(PC pc)
{
  PC_LMVM        *ctx = (PC_LMVM*)pc->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (ctx->inactive) {
    ierr = ISDestroy(&ctx->inactive);CHKERRQ(ierr);
  }
  if (pc->setupcalled) {
    ierr = VecDestroy(&ctx->xwork);CHKERRQ(ierr);
    ierr = VecDestroy(&ctx->ywork);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&ctx->B);CHKERRQ(ierr);
  ierr = PetscFree(pc->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*MC
   PCLMVM - Creates a preconditioner around an LMVM matrix. Options for the 
            underlying LMVM matrix can be access with the "-pc_lmvm_" prefix.

   Level: intermediate

.seealso:  PCCreate(), PCSetType(), PCType (for list of available types), 
           PC, MATLMVM, PCLMVMUpdate(), PCLMVMSetMatLMVM(), PCLMVMGetMatLMVM()
M*/
PETSC_EXTERN PetscErrorCode PCCreate_LMVM(PC pc)
{
  PetscErrorCode ierr;
  PC_LMVM        *ctx;

  PetscFunctionBegin;
  ierr     = PetscNewLog(pc,&ctx);CHKERRQ(ierr);
  pc->data = (void*)ctx;
  
  pc->ops->reset           = PCReset_LMVM;
  pc->ops->setup           = PCSetUp_LMVM;
  pc->ops->destroy         = PCDestroy_LMVM;
  pc->ops->view            = PCView_LMVM;
  pc->ops->apply           = PCApply_LMVM;
  pc->ops->setfromoptions  = PCSetFromOptions_LMVM;
  pc->ops->applysymmetricleft  = 0;
  pc->ops->applysymmetricright = 0;
  pc->ops->applytranspose  = 0;
  pc->ops->applyrichardson = 0;
  pc->ops->presolve        = 0;
  pc->ops->postsolve       = 0;
  
  ierr = PCSetReusePreconditioner(pc, PETSC_TRUE);CHKERRQ(ierr);
  
  ierr = MatCreate(PetscObjectComm((PetscObject)pc), &ctx->B);CHKERRQ(ierr);
  ierr = MatSetType(ctx->B, MATLMVMBFGS);CHKERRQ(ierr);
  ierr = PetscObjectIncrementTabLevel((PetscObject)ctx->B, (PetscObject)pc, 1);CHKERRQ(ierr);
  ierr = MatSetOptionsPrefix(ctx->B, "pc_lmvm_");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}






