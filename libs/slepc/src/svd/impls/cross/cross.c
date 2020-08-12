/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc singular value solver: "cross"

   Method: Uses a Hermitian eigensolver for A^T*A
*/

#include <slepc/private/svdimpl.h>                /*I "slepcsvd.h" I*/
#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/

typedef struct {
  PetscBool explicitmatrix;
  EPS       eps;
  PetscBool usereps;
  Mat       mat;
  Vec       w,diag;
} SVD_CROSS;

static PetscErrorCode MatMult_Cross(Mat B,Vec x,Vec y)
{
  PetscErrorCode ierr;
  SVD            svd;
  SVD_CROSS      *cross;

  PetscFunctionBegin;
  ierr = MatShellGetContext(B,(void**)&svd);CHKERRQ(ierr);
  cross = (SVD_CROSS*)svd->data;
  ierr = SVDMatMult(svd,PETSC_FALSE,x,cross->w);CHKERRQ(ierr);
  ierr = SVDMatMult(svd,PETSC_TRUE,cross->w,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatCreateVecs_Cross(Mat B,Vec *right,Vec *left)
{
  PetscErrorCode ierr;
  SVD            svd;

  PetscFunctionBegin;
  ierr = MatShellGetContext(B,(void**)&svd);CHKERRQ(ierr);
  if (right) {
    ierr = SVDMatCreateVecs(svd,right,NULL);CHKERRQ(ierr);
    if (left) { ierr = VecDuplicate(*right,left);CHKERRQ(ierr); }
  } else {
    ierr = SVDMatCreateVecs(svd,left,NULL);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatGetDiagonal_Cross(Mat B,Vec d)
{
  PetscErrorCode    ierr;
  SVD               svd;
  SVD_CROSS         *cross;
  PetscMPIInt       len;
  PetscInt          N,n,i,j,start,end,ncols;
  PetscScalar       *work1,*work2,*diag;
  const PetscInt    *cols;
  const PetscScalar *vals;

  PetscFunctionBegin;
  ierr = MatShellGetContext(B,(void**)&svd);CHKERRQ(ierr);
  cross = (SVD_CROSS*)svd->data;
  if (!cross->diag) {
    /* compute diagonal from rows and store in cross->diag */
    ierr = VecDuplicate(d,&cross->diag);CHKERRQ(ierr);
    ierr = SVDMatGetSize(svd,NULL,&N);CHKERRQ(ierr);
    ierr = SVDMatGetLocalSize(svd,NULL,&n);CHKERRQ(ierr);
    ierr = PetscCalloc2(N,&work1,N,&work2);CHKERRQ(ierr);
    if (svd->AT) {
      ierr = MatGetOwnershipRange(svd->AT,&start,&end);CHKERRQ(ierr);
      for (i=start;i<end;i++) {
        ierr = MatGetRow(svd->AT,i,&ncols,NULL,&vals);CHKERRQ(ierr);
        for (j=0;j<ncols;j++)
          work1[i] += vals[j]*vals[j];
        ierr = MatRestoreRow(svd->AT,i,&ncols,NULL,&vals);CHKERRQ(ierr);
      }
    } else {
      ierr = MatGetOwnershipRange(svd->A,&start,&end);CHKERRQ(ierr);
      for (i=start;i<end;i++) {
        ierr = MatGetRow(svd->A,i,&ncols,&cols,&vals);CHKERRQ(ierr);
        for (j=0;j<ncols;j++)
          work1[cols[j]] += vals[j]*vals[j];
        ierr = MatRestoreRow(svd->A,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      }
    }
    ierr = PetscMPIIntCast(N,&len);CHKERRQ(ierr);
    ierr = MPI_Allreduce(work1,work2,len,MPIU_SCALAR,MPIU_SUM,PetscObjectComm((PetscObject)svd));CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(cross->diag,&start,&end);CHKERRQ(ierr);
    ierr = VecGetArray(cross->diag,&diag);CHKERRQ(ierr);
    for (i=start;i<end;i++) diag[i-start] = work2[i];
    ierr = VecRestoreArray(cross->diag,&diag);CHKERRQ(ierr);
    ierr = PetscFree2(work1,work2);CHKERRQ(ierr);
  }
  ierr = VecCopy(cross->diag,d);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSetUp_Cross(SVD svd)
{
  PetscErrorCode ierr;
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  PetscInt       n;
  PetscBool      trackall;

  PetscFunctionBegin;
  if (!cross->mat) {
    if (cross->explicitmatrix) {
      if (svd->A && svd->AT) {  /* explicit transpose */
        ierr = MatProductCreate(svd->AT,svd->A,NULL,&cross->mat);CHKERRQ(ierr);
        ierr = MatProductSetType(cross->mat,MATPRODUCT_AB);CHKERRQ(ierr);
      } else {  /* implicit transpose */
#if defined(PETSC_USE_COMPLEX)
        SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"Must use explicit transpose with complex scalars");
#else
        if (svd->A) {
          ierr = MatProductCreate(svd->A,svd->A,NULL,&cross->mat);CHKERRQ(ierr);
          ierr = MatProductSetType(cross->mat,MATPRODUCT_AtB);CHKERRQ(ierr);
        } else {
          ierr = MatProductCreate(svd->AT,svd->AT,NULL,&cross->mat);CHKERRQ(ierr);
          ierr = MatProductSetType(cross->mat,MATPRODUCT_ABt);CHKERRQ(ierr);
        }
#endif
      }
      ierr = MatProductSetFromOptions(cross->mat);CHKERRQ(ierr);
      ierr = MatProductSymbolic(cross->mat);CHKERRQ(ierr);
      ierr = MatProductNumeric(cross->mat);CHKERRQ(ierr);
    } else {
      ierr = SVDMatGetLocalSize(svd,NULL,&n);CHKERRQ(ierr);
      ierr = MatCreateShell(PetscObjectComm((PetscObject)svd),n,n,PETSC_DETERMINE,PETSC_DETERMINE,svd,&cross->mat);CHKERRQ(ierr);
      ierr = MatShellSetOperation(cross->mat,MATOP_MULT,(void(*)(void))MatMult_Cross);CHKERRQ(ierr);
      ierr = MatShellSetOperation(cross->mat,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_Cross);CHKERRQ(ierr);
      ierr = MatShellSetOperation(cross->mat,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Cross);CHKERRQ(ierr);
      ierr = SVDMatCreateVecs(svd,NULL,&cross->w);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)cross->w);CHKERRQ(ierr);
    }
    ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)cross->mat);CHKERRQ(ierr);
  }

  if (!cross->eps) { ierr = SVDCrossGetEPS(svd,&cross->eps);CHKERRQ(ierr); }
  ierr = EPSSetOperators(cross->eps,cross->mat,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(cross->eps,EPS_HEP);CHKERRQ(ierr);
  if (!cross->usereps) {
    ierr = EPSSetWhichEigenpairs(cross->eps,svd->which==SVD_LARGEST?EPS_LARGEST_REAL:EPS_SMALLEST_REAL);CHKERRQ(ierr);
    ierr = EPSSetDimensions(cross->eps,svd->nsv,svd->ncv,svd->mpd);CHKERRQ(ierr);
    ierr = EPSSetTolerances(cross->eps,svd->tol==PETSC_DEFAULT?SLEPC_DEFAULT_TOL/10.0:svd->tol,svd->max_it);CHKERRQ(ierr);
    switch (svd->conv) {
    case SVD_CONV_ABS:
      ierr = EPSSetConvergenceTest(cross->eps,EPS_CONV_ABS);CHKERRQ(ierr);break;
    case SVD_CONV_REL:
      ierr = EPSSetConvergenceTest(cross->eps,EPS_CONV_REL);CHKERRQ(ierr);break;
    case SVD_CONV_USER:
      SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"User-defined convergence test not supported in this solver");
    }
  }
  if (svd->stop!=SVD_STOP_BASIC) SETERRQ(PetscObjectComm((PetscObject)svd),PETSC_ERR_SUP,"User-defined stopping test not supported in this solver");
  /* Transfer the trackall option from svd to eps */
  ierr = SVDGetTrackAll(svd,&trackall);CHKERRQ(ierr);
  ierr = EPSSetTrackAll(cross->eps,trackall);CHKERRQ(ierr);
  /* Transfer the initial space from svd to eps */
  if (svd->nini<0) {
    ierr = EPSSetInitialSpace(cross->eps,-svd->nini,svd->IS);CHKERRQ(ierr);
    ierr = SlepcBasisDestroy_Private(&svd->nini,&svd->IS);CHKERRQ(ierr);
  }
  ierr = EPSSetUp(cross->eps);CHKERRQ(ierr);
  ierr = EPSGetDimensions(cross->eps,NULL,&svd->ncv,&svd->mpd);CHKERRQ(ierr);
  ierr = EPSGetTolerances(cross->eps,NULL,&svd->max_it);CHKERRQ(ierr);
  if (svd->tol==PETSC_DEFAULT) svd->tol = SLEPC_DEFAULT_TOL;

  svd->leftbasis = PETSC_FALSE;
  ierr = SVDAllocateSolution(svd,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSolve_Cross(SVD svd)
{
  PetscErrorCode ierr;
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  PetscInt       i;
  PetscScalar    lambda;
  PetscReal      sigma;
  Vec            v;

  PetscFunctionBegin;
  ierr = EPSSolve(cross->eps);CHKERRQ(ierr);
  ierr = EPSGetConverged(cross->eps,&svd->nconv);CHKERRQ(ierr);
  ierr = EPSGetIterationNumber(cross->eps,&svd->its);CHKERRQ(ierr);
  ierr = EPSGetConvergedReason(cross->eps,(EPSConvergedReason*)&svd->reason);CHKERRQ(ierr);
  for (i=0;i<svd->nconv;i++) {
    ierr = BVGetColumn(svd->V,i,&v);CHKERRQ(ierr);
    ierr = EPSGetEigenpair(cross->eps,i,&lambda,NULL,v,NULL);CHKERRQ(ierr);
    ierr = BVRestoreColumn(svd->V,i,&v);CHKERRQ(ierr);
    sigma = PetscRealPart(lambda);
    if (sigma<-10*PETSC_MACHINE_EPSILON) SETERRQ1(PetscObjectComm((PetscObject)svd),1,"Negative eigenvalue computed by EPS: %g",sigma);
    if (sigma<0.0) {
      ierr = PetscInfo1(svd,"Negative eigenvalue computed by EPS: %g, resetting to 0\n",sigma);CHKERRQ(ierr);
      sigma = 0.0;
    }
    svd->sigma[i] = PetscSqrtReal(sigma);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSMonitor_Cross(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *ctx)
{
  PetscInt       i;
  SVD            svd = (SVD)ctx;
  PetscScalar    er,ei;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (i=0;i<PetscMin(nest,svd->ncv);i++) {
    er = eigr[i]; ei = eigi[i];
    ierr = STBackTransform(eps->st,1,&er,&ei);CHKERRQ(ierr);
    svd->sigma[i] = PetscSqrtReal(PetscAbsReal(PetscRealPart(er)));
    svd->errest[i] = errest[i];
  }
  ierr = SVDMonitor(svd,its,nconv,svd->sigma,svd->errest,nest);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDSetFromOptions_Cross(PetscOptionItems *PetscOptionsObject,SVD svd)
{
  PetscErrorCode ierr;
  PetscBool      set,val;
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  ST             st;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"SVD Cross Options");CHKERRQ(ierr);

    ierr = PetscOptionsBool("-svd_cross_explicitmatrix","Use cross explicit matrix","SVDCrossSetExplicitMatrix",cross->explicitmatrix,&val,&set);CHKERRQ(ierr);
    if (set) { ierr = SVDCrossSetExplicitMatrix(svd,val);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);

  if (!cross->eps) { ierr = SVDCrossGetEPS(svd,&cross->eps);CHKERRQ(ierr); }
  if (!cross->explicitmatrix && !cross->usereps) {
    /* use as default an ST with shell matrix and Jacobi */
    ierr = EPSGetST(cross->eps,&st);CHKERRQ(ierr);
    ierr = STSetMatMode(st,ST_MATMODE_SHELL);CHKERRQ(ierr);
  }
  ierr = EPSSetFromOptions(cross->eps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDCrossSetExplicitMatrix_Cross(SVD svd,PetscBool explicitmatrix)
{
  SVD_CROSS *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  if (cross->explicitmatrix != explicitmatrix) {
    cross->explicitmatrix = explicitmatrix;
    svd->state = SVD_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   SVDCrossSetExplicitMatrix - Indicate if the eigensolver operator A^T*A must
   be computed explicitly.

   Logically Collective on svd

   Input Parameters:
+  svd      - singular value solver
-  explicit - boolean flag indicating if A^T*A is built explicitly

   Options Database Key:
.  -svd_cross_explicitmatrix <boolean> - Indicates the boolean flag

   Level: advanced

.seealso: SVDCrossGetExplicitMatrix()
@*/
PetscErrorCode SVDCrossSetExplicitMatrix(SVD svd,PetscBool explicitmatrix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidLogicalCollectiveBool(svd,explicitmatrix,2);
  ierr = PetscTryMethod(svd,"SVDCrossSetExplicitMatrix_C",(SVD,PetscBool),(svd,explicitmatrix));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDCrossGetExplicitMatrix_Cross(SVD svd,PetscBool *explicitmatrix)
{
  SVD_CROSS *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  *explicitmatrix = cross->explicitmatrix;
  PetscFunctionReturn(0);
}

/*@
   SVDCrossGetExplicitMatrix - Returns the flag indicating if A^T*A is built explicitly.

   Not Collective

   Input Parameter:
.  svd  - singular value solver

   Output Parameter:
.  explicit - the mode flag

   Level: advanced

.seealso: SVDCrossSetExplicitMatrix()
@*/
PetscErrorCode SVDCrossGetExplicitMatrix(SVD svd,PetscBool *explicitmatrix)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidBoolPointer(explicitmatrix,2);
  ierr = PetscUseMethod(svd,"SVDCrossGetExplicitMatrix_C",(SVD,PetscBool*),(svd,explicitmatrix));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDCrossSetEPS_Cross(SVD svd,EPS eps)
{
  PetscErrorCode ierr;
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)eps);CHKERRQ(ierr);
  ierr = EPSDestroy(&cross->eps);CHKERRQ(ierr);
  cross->eps = eps;
  cross->usereps = PETSC_TRUE;
  ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)cross->eps);CHKERRQ(ierr);
  svd->state = SVD_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   SVDCrossSetEPS - Associate an eigensolver object (EPS) to the
   singular value solver.

   Collective on svd

   Input Parameters:
+  svd - singular value solver
-  eps - the eigensolver object

   Level: advanced

.seealso: SVDCrossGetEPS()
@*/
PetscErrorCode SVDCrossSetEPS(SVD svd,EPS eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidHeaderSpecific(eps,EPS_CLASSID,2);
  PetscCheckSameComm(svd,1,eps,2);
  ierr = PetscTryMethod(svd,"SVDCrossSetEPS_C",(SVD,EPS),(svd,eps));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SVDCrossGetEPS_Cross(SVD svd,EPS *eps)
{
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!cross->eps) {
    ierr = EPSCreate(PetscObjectComm((PetscObject)svd),&cross->eps);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)cross->eps,(PetscObject)svd,1);CHKERRQ(ierr);
    ierr = EPSSetOptionsPrefix(cross->eps,((PetscObject)svd)->prefix);CHKERRQ(ierr);
    ierr = EPSAppendOptionsPrefix(cross->eps,"svd_cross_");CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)svd,(PetscObject)cross->eps);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)cross->eps,((PetscObject)svd)->options);CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(cross->eps,EPS_LARGEST_REAL);CHKERRQ(ierr);
    ierr = EPSMonitorSet(cross->eps,EPSMonitor_Cross,svd,NULL);CHKERRQ(ierr);
  }
  *eps = cross->eps;
  PetscFunctionReturn(0);
}

/*@
   SVDCrossGetEPS - Retrieve the eigensolver object (EPS) associated
   to the singular value solver.

   Not Collective

   Input Parameter:
.  svd - singular value solver

   Output Parameter:
.  eps - the eigensolver object

   Level: advanced

.seealso: SVDCrossSetEPS()
@*/
PetscErrorCode SVDCrossGetEPS(SVD svd,EPS *eps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_CLASSID,1);
  PetscValidPointer(eps,2);
  ierr = PetscUseMethod(svd,"SVDCrossGetEPS_C",(SVD,EPS*),(svd,eps));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDView_Cross(SVD svd,PetscViewer viewer)
{
  PetscErrorCode ierr;
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    if (!cross->eps) { ierr = SVDCrossGetEPS(svd,&cross->eps);CHKERRQ(ierr); }
    ierr = PetscViewerASCIIPrintf(viewer,"  %s matrix\n",cross->explicitmatrix?"explicit":"implicit");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = EPSView(cross->eps,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SVDReset_Cross(SVD svd)
{
  PetscErrorCode ierr;
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  ierr = EPSReset(cross->eps);CHKERRQ(ierr);
  ierr = MatDestroy(&cross->mat);CHKERRQ(ierr);
  ierr = VecDestroy(&cross->w);CHKERRQ(ierr);
  ierr = VecDestroy(&cross->diag);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SVDDestroy_Cross(SVD svd)
{
  PetscErrorCode ierr;
  SVD_CROSS      *cross = (SVD_CROSS*)svd->data;

  PetscFunctionBegin;
  ierr = EPSDestroy(&cross->eps);CHKERRQ(ierr);
  ierr = PetscFree(svd->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetEPS_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetEPS_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetExplicitMatrix_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetExplicitMatrix_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode SVDCreate_Cross(SVD svd)
{
  PetscErrorCode ierr;
  SVD_CROSS      *cross;

  PetscFunctionBegin;
  ierr = PetscNewLog(svd,&cross);CHKERRQ(ierr);
  svd->data = (void*)cross;

  svd->ops->solve          = SVDSolve_Cross;
  svd->ops->setup          = SVDSetUp_Cross;
  svd->ops->setfromoptions = SVDSetFromOptions_Cross;
  svd->ops->destroy        = SVDDestroy_Cross;
  svd->ops->reset          = SVDReset_Cross;
  svd->ops->view           = SVDView_Cross;
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetEPS_C",SVDCrossSetEPS_Cross);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetEPS_C",SVDCrossGetEPS_Cross);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDCrossSetExplicitMatrix_C",SVDCrossSetExplicitMatrix_Cross);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)svd,"SVDCrossGetExplicitMatrix_C",SVDCrossGetExplicitMatrix_Cross);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

