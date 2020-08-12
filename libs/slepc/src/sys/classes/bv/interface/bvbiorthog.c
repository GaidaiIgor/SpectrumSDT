/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV bi-orthogonalization routines
*/

#include <slepc/private/bvimpl.h>          /*I   "slepcbv.h"   I*/

/*
   BVBiorthogonalizeMGS1 - Compute one step of Modified Gram-Schmidt bi-orthogonalization
*/
static PetscErrorCode BVBiorthogonalizeMGS1(BV V,BV W,Vec v,PetscScalar *h,PetscScalar *c)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    dot;
  Vec            vi,wi;

  PetscFunctionBegin;
  for (i=-V->nc;i<V->k;i++) {
    ierr = BVGetColumn(W,i,&wi);CHKERRQ(ierr);
    /* h_i = (v, w_i) */
    ierr = VecDot(v,wi,&dot);CHKERRQ(ierr);
    ierr = BVRestoreColumn(W,i,&wi);CHKERRQ(ierr);
    /* v <- v - h_i v_i */
    ierr = BV_SetValue(V,i,0,c,dot);CHKERRQ(ierr);
    ierr = BVGetColumn(V,i,&vi);CHKERRQ(ierr);
    ierr = VecAXPY(v,-dot,vi);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,i,&vi);CHKERRQ(ierr);
  }
  ierr = BV_AddCoefficients(V,V->k,h,c);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   BVBiorthogonalizeCGS1 - Compute one step of CGS bi-orthogonalization: v = (I-V*W')*v
*/
static PetscErrorCode BVBiorthogonalizeCGS1(BV V,BV W,Vec v,PetscScalar *h,PetscScalar *c)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* h = W'*v */
  ierr = BVDotVec(W,v,c);CHKERRQ(ierr);

  /* v = v - V h */
  ierr = BVMultVec(V,-1.0,1.0,v,c);CHKERRQ(ierr);

  ierr = BV_AddCoefficients(V,V->k,h,c);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define BVBiorthogonalizeGS1(a,b,c,d,e) ((V->orthog_type==BV_ORTHOG_MGS)?BVBiorthogonalizeMGS1:BVBiorthogonalizeCGS1)(a,b,c,d,e)

/*
   BVBiorthogonalizeGS - Orthogonalize with (classical or modified) Gram-Schmidt

   V, W - the two basis vectors objects
   v    - the vector to bi-orthogonalize
*/
static PetscErrorCode BVBiorthogonalizeGS(BV V,BV W,Vec v)
{
  PetscErrorCode ierr;
  PetscScalar    *h,*c;

  PetscFunctionBegin;
  h = V->h;
  c = V->c;
  ierr = BV_CleanCoefficients(V,V->k,h);CHKERRQ(ierr);
  ierr = BVBiorthogonalizeGS1(V,W,v,h,c);CHKERRQ(ierr);
  if (V->orthog_ref!=BV_ORTHOG_REFINE_NEVER) {
    ierr = BVBiorthogonalizeGS1(V,W,v,h,c);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   BVBiorthogonalizeColumn - Bi-orthogonalize a column of two BV objects.

   Collective on V

   Input Parameters:
+  V,W - two basis vectors contexts
-  j   - index of column to be bi-orthonormalized

   Notes:
   This function bi-orthogonalizes vectors V[j],W[j] against W[0..j-1],
   and V[0..j-1], respectively, so that W[0..j]'*V[0..j] = diagonal.

   Level: advanced

.seealso: BVOrthogonalizeColumn(), BVBiorthonormalizeColumn()
@*/
PetscErrorCode BVBiorthogonalizeColumn(BV V,BV W,PetscInt j)
{
  PetscErrorCode ierr;
  PetscInt       ksavev,lsavev,ksavew,lsavew;
  Vec            y,z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidHeaderSpecific(W,BV_CLASSID,2);
  PetscValidLogicalCollectiveInt(V,j,3);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidType(W,2);
  BVCheckSizes(W,2);
  PetscCheckSameTypeAndComm(V,1,W,2);
  if (j<0) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  if (j>=V->m) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%D but V only has %D columns",j,V->m);
  if (j>=W->m) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%D but W only has %D columns",j,W->m);
  if (V->n!=W->n) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension V %D, W %D",V->n,W->n);
  if (V->matrix || W->matrix) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_WRONGSTATE,"V,W must not have an inner product matrix");
  if (V->nc || W->nc) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_WRONGSTATE,"V,W cannot have different number of constraints");
  if (V->ops->gramschmidt || W->ops->gramschmidt) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Object has a special GS function");

  /* bi-orthogonalize */
  ierr = PetscLogEventBegin(BV_OrthogonalizeVec,V,0,0,0);CHKERRQ(ierr);
  ksavev = V->k;
  lsavev = V->l;
  ksavew = W->k;
  lsavew = W->l;
  V->k = j;
  V->l = -V->nc;  /* must also bi-orthogonalize against constraints and leading columns */
  W->k = j;
  W->l = -W->nc;
  ierr = BV_AllocateCoeffs(V);CHKERRQ(ierr);
  ierr = BV_AllocateCoeffs(W);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j,&y);CHKERRQ(ierr);
  ierr = BVBiorthogonalizeGS(V,W,y);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j,&y);CHKERRQ(ierr);
  ierr = BVGetColumn(W,j,&z);CHKERRQ(ierr);
  ierr = BVBiorthogonalizeGS(W,V,z);CHKERRQ(ierr);
  ierr = BVRestoreColumn(W,j,&z);CHKERRQ(ierr);
  V->k = ksavev;
  V->l = lsavev;
  W->k = ksavew;
  W->l = lsavew;
  ierr = PetscLogEventEnd(BV_OrthogonalizeVec,V,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)W);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVBiorthonormalizeColumn - Bi-orthonormalize a column of two BV objects.

   Collective on V

   Input Parameters:
+  V,W - two basis vectors contexts
-  j   - index of column to be bi-orthonormalized

   Output Parameters:
.  delta - (optional) value used for normalization

   Notes:
   This function first bi-orthogonalizes vectors V[j],W[j] against W[0..j-1],
   and V[0..j-1], respectively. Then, it scales the vectors with 1/delta, so
   that the resulting vectors satisfy W[j]'*V[j] = 1.

   Level: advanced

.seealso: BVOrthonormalizeColumn(), BVBiorthogonalizeColumn()
@*/
PetscErrorCode BVBiorthonormalizeColumn(BV V,BV W,PetscInt j,PetscReal *delta)
{
  PetscErrorCode ierr;
  PetscScalar    alpha;
  PetscReal      deltat;
  PetscInt       ksavev,lsavev,ksavew,lsavew;
  Vec            y,z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidHeaderSpecific(W,BV_CLASSID,2);
  PetscValidLogicalCollectiveInt(V,j,3);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidType(W,2);
  BVCheckSizes(W,2);
  PetscCheckSameTypeAndComm(V,1,W,2);
  if (j<0) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  if (j>=V->m) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%D but V only has %D columns",j,V->m);
  if (j>=W->m) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%D but W only has %D columns",j,W->m);
  if (V->n!=W->n) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension V %D, W %D",V->n,W->n);
  if (V->matrix || W->matrix) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_WRONGSTATE,"V,W must not have an inner product matrix");
  if (V->nc || W->nc) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_WRONGSTATE,"V,W cannot have different number of constraints");
  if (V->ops->gramschmidt || W->ops->gramschmidt) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Object has a special GS function");

  /* bi-orthogonalize */
  ierr = PetscLogEventBegin(BV_OrthogonalizeVec,V,0,0,0);CHKERRQ(ierr);
  ksavev = V->k;
  lsavev = V->l;
  ksavew = W->k;
  lsavew = W->l;
  V->k = j;
  V->l = -V->nc;  /* must also bi-orthogonalize against constraints and leading columns */
  W->k = j;
  W->l = -W->nc;
  ierr = BV_AllocateCoeffs(V);CHKERRQ(ierr);
  ierr = BV_AllocateCoeffs(W);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j,&y);CHKERRQ(ierr);
  ierr = BVBiorthogonalizeGS(V,W,y);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j,&y);CHKERRQ(ierr);
  ierr = BVGetColumn(W,j,&z);CHKERRQ(ierr);
  ierr = BVBiorthogonalizeGS(W,V,z);CHKERRQ(ierr);
  ierr = BVRestoreColumn(W,j,&z);CHKERRQ(ierr);
  V->k = ksavev;
  V->l = lsavev;
  W->k = ksavew;
  W->l = lsavew;
  ierr = PetscLogEventEnd(BV_OrthogonalizeVec,V,0,0,0);CHKERRQ(ierr);

  /* scale */
  ierr = PetscLogEventBegin(BV_Scale,V,0,0,0);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j,&y);CHKERRQ(ierr);
  ierr = BVGetColumn(W,j,&z);CHKERRQ(ierr);
  ierr = VecDot(z,y,&alpha);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j,&y);CHKERRQ(ierr);
  ierr = BVRestoreColumn(W,j,&z);CHKERRQ(ierr);
  deltat = PetscSqrtReal(PetscAbsScalar(alpha));
  if (V->n) { ierr = (*V->ops->scale)(V,j,1.0/PetscConj(alpha/deltat));CHKERRQ(ierr); }
  if (W->n) { ierr = (*W->ops->scale)(W,j,1.0/deltat);CHKERRQ(ierr); }
  ierr = PetscLogEventEnd(BV_Scale,V,0,0,0);CHKERRQ(ierr);
  if (delta) *delta = deltat;
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)W);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

