/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV operations involving global communication
*/

#include <slepc/private/bvimpl.h>      /*I "slepcbv.h" I*/

/*
  BVDot for the particular case of non-standard inner product with
  matrix B, which is assumed to be symmetric (or complex Hermitian)
*/
PETSC_STATIC_INLINE PetscErrorCode BVDot_Private(BV X,BV Y,Mat M)
{
  PetscErrorCode ierr;
  PetscObjectId  idx,idy;
  PetscInt       i,j,jend,m;
  PetscScalar    *marray;
  PetscBool      symm=PETSC_FALSE;
  Mat            B;
  Vec            z;

  PetscFunctionBegin;
  BVCheckOp(Y,1,dotvec);
  ierr = MatGetSize(M,&m,NULL);CHKERRQ(ierr);
  ierr = MatDenseGetArray(M,&marray);CHKERRQ(ierr);
  ierr = PetscObjectGetId((PetscObject)X,&idx);CHKERRQ(ierr);
  ierr = PetscObjectGetId((PetscObject)Y,&idy);CHKERRQ(ierr);
  B = Y->matrix;
  Y->matrix = NULL;
  if (idx==idy) symm=PETSC_TRUE;  /* M=X'BX is symmetric */
  jend = X->k;
  for (j=X->l;j<jend;j++) {
    if (symm) Y->k = j+1;
    ierr = BVGetColumn(X->cached,j,&z);CHKERRQ(ierr);
    ierr = (*Y->ops->dotvec)(Y,z,marray+j*m+Y->l);CHKERRQ(ierr);
    ierr = BVRestoreColumn(X->cached,j,&z);CHKERRQ(ierr);
    if (symm) {
      for (i=X->l;i<j;i++)
        marray[j+i*m] = PetscConj(marray[i+j*m]);
    }
  }
  ierr = MatDenseRestoreArray(M,&marray);CHKERRQ(ierr);
  Y->matrix = B;
  PetscFunctionReturn(0);
}

/*@
   BVDot - Computes the 'block-dot' product of two basis vectors objects.

   Collective on X

   Input Parameters:
+  X, Y - basis vectors
-  M    - Mat object where the result must be placed

   Output Parameter:
.  M    - the resulting matrix

   Notes:
   This is the generalization of VecDot() for a collection of vectors, M = Y^H*X.
   The result is a matrix M whose entry m_ij is equal to y_i^H x_j (where y_i^H
   denotes the conjugate transpose of y_i).

   If a non-standard inner product has been specified with BVSetMatrix(),
   then the result is M = Y^H*B*X. In this case, both X and Y must have
   the same associated matrix.

   On entry, M must be a sequential dense Mat with dimensions m,n at least, where
   m is the number of active columns of Y and n is the number of active columns of X.
   Only rows (resp. columns) of M starting from ly (resp. lx) are computed,
   where ly (resp. lx) is the number of leading columns of Y (resp. X).

   X and Y need not be different objects.

   Level: intermediate

.seealso: BVDotVec(), BVDotColumn(), BVSetActiveColumns(), BVSetMatrix()
@*/
PetscErrorCode BVDot(BV X,BV Y,Mat M)
{
  PetscErrorCode ierr;
  PetscBool      match;
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidHeaderSpecific(Y,BV_CLASSID,2);
  PetscValidHeaderSpecific(M,MAT_CLASSID,3);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  PetscValidType(Y,2);
  BVCheckSizes(Y,2);
  PetscValidType(M,3);
  PetscCheckSameTypeAndComm(X,1,Y,2);
  ierr = PetscObjectTypeCompare((PetscObject)M,MATSEQDENSE,&match);CHKERRQ(ierr);
  if (!match) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_SUP,"Mat argument must be of type seqdense");

  ierr = MatGetSize(M,&m,&n);CHKERRQ(ierr);
  if (m<Y->k) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Mat argument has %D rows, should have at least %D",m,Y->k);
  if (n<X->k) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Mat argument has %D columns, should have at least %D",n,X->k);
  if (X->n!=Y->n) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %D, Y %D",X->n,Y->n);
  if (X->matrix!=Y->matrix) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"X and Y must have the same inner product matrix");
  if (X->l==X->k || Y->l==Y->k) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(BV_Dot,X,Y,0,0);CHKERRQ(ierr);
  if (X->matrix) { /* non-standard inner product */
    /* compute BX first */
    ierr = BV_IPMatMultBV(X);CHKERRQ(ierr);
    if (X->vmm==BV_MATMULT_VECS) {
      /* perform computation column by column */
      ierr = BVDot_Private(X,Y,M);CHKERRQ(ierr);
    } else {
      ierr = (*X->ops->dot)(X->cached,Y,M);CHKERRQ(ierr);
    }
  } else {
    ierr = (*X->ops->dot)(X,Y,M);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(BV_Dot,X,Y,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVDotVec - Computes multiple dot products of a vector against all the
   column vectors of a BV.

   Collective on X

   Input Parameters:
+  X - basis vectors
-  y - a vector

   Output Parameter:
.  m - an array where the result must be placed

   Notes:
   This is analogue to VecMDot(), but using BV to represent a collection
   of vectors. The result is m = X^H*y, so m_i is equal to x_j^H y. Note
   that here X is transposed as opposed to BVDot().

   If a non-standard inner product has been specified with BVSetMatrix(),
   then the result is m = X^H*B*y.

   The length of array m must be equal to the number of active columns of X
   minus the number of leading columns, i.e. the first entry of m is the
   product of the first non-leading column with y.

   Level: intermediate

.seealso: BVDot(), BVDotColumn(), BVSetActiveColumns(), BVSetMatrix()
@*/
PetscErrorCode BVDotVec(BV X,Vec y,PetscScalar m[])
{
  PetscErrorCode ierr;
  PetscInt       n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidHeaderSpecific(y,VEC_CLASSID,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  BVCheckOp(X,1,dotvec);
  PetscValidType(y,2);
  PetscCheckSameTypeAndComm(X,1,y,2);

  ierr = VecGetLocalSize(y,&n);CHKERRQ(ierr);
  if (X->n!=n) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %D, y %D",X->n,n);

  ierr = PetscLogEventBegin(BV_DotVec,X,y,0,0);CHKERRQ(ierr);
  ierr = (*X->ops->dotvec)(X,y,m);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_DotVec,X,y,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVDotVecBegin - Starts a split phase dot product computation.

   Input Parameters:
+  X - basis vectors
.  y - a vector
-  m - an array where the result will go (can be NULL)

   Note:
   Each call to BVDotVecBegin() should be paired with a call to BVDotVecEnd().

   Level: advanced

.seealso: BVDotVecEnd(), BVDotVec()
@*/
PetscErrorCode BVDotVecBegin(BV X,Vec y,PetscScalar *m)
{
  PetscErrorCode      ierr;
  PetscInt            i,n,nv;
  PetscSplitReduction *sr;
  MPI_Comm            comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidHeaderSpecific(y,VEC_CLASSID,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  PetscValidType(y,2);
  PetscCheckSameTypeAndComm(X,1,y,2);

  ierr = VecGetLocalSize(y,&n);CHKERRQ(ierr);
  if (X->n!=n) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %D, y %D",X->n,n);

  if (X->ops->dotvec_begin) {
    ierr = (*X->ops->dotvec_begin)(X,y,m);CHKERRQ(ierr);
  } else {
    BVCheckOp(X,1,dotvec_local);
    nv = X->k-X->l;
    ierr = PetscObjectGetComm((PetscObject)X,&comm);CHKERRQ(ierr);
    ierr = PetscSplitReductionGet(comm,&sr);CHKERRQ(ierr);
    if (sr->state != STATE_BEGIN) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ORDER,"Called before all VecxxxEnd() called");
    for (i=0;i<nv;i++) {
      if (sr->numopsbegin+i >= sr->maxops) {
        ierr = PetscSplitReductionExtend(sr);CHKERRQ(ierr);
      }
      sr->reducetype[sr->numopsbegin+i] = PETSC_SR_REDUCE_SUM;
      sr->invecs[sr->numopsbegin+i]     = (void*)X;
    }
    ierr = PetscLogEventBegin(BV_DotVec,X,y,0,0);CHKERRQ(ierr);
    ierr = (*X->ops->dotvec_local)(X,y,sr->lvalues+sr->numopsbegin);CHKERRQ(ierr);
    sr->numopsbegin += nv;
    ierr = PetscLogEventEnd(BV_DotVec,X,y,0,0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   BVDotVecEnd - Ends a split phase dot product computation.

   Input Parameters:
+  X - basis vectors
.  y - a vector
-  m - an array where the result will go

   Note:
   Each call to BVDotVecBegin() should be paired with a call to BVDotVecEnd().

   Level: advanced

.seealso: BVDotVecBegin(), BVDotVec()
@*/
PetscErrorCode BVDotVecEnd(BV X,Vec y,PetscScalar *m)
{
  PetscErrorCode      ierr;
  PetscInt            i,nv;
  PetscSplitReduction *sr;
  MPI_Comm            comm;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidType(X,1);
  BVCheckSizes(X,1);

  if (X->ops->dotvec_end) {
    ierr = (*X->ops->dotvec_end)(X,y,m);CHKERRQ(ierr);
  } else {
    nv = X->k-X->l;
    ierr = PetscObjectGetComm((PetscObject)X,&comm);CHKERRQ(ierr);
    ierr = PetscSplitReductionGet(comm,&sr);CHKERRQ(ierr);
    ierr = PetscSplitReductionEnd(sr);CHKERRQ(ierr);

    if (sr->numopsend >= sr->numopsbegin) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Called VecxxxEnd() more times than VecxxxBegin()");
    if ((void*)X != sr->invecs[sr->numopsend]) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Called BVxxxEnd() in a different order or with a different BV than BVxxxBegin()");
    if (sr->reducetype[sr->numopsend] != PETSC_SR_REDUCE_SUM) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Wrong type of reduction");
    for (i=0;i<nv;i++) m[i] = sr->gvalues[sr->numopsend++];

    /* Finished getting all the results so reset to no outstanding requests */
    if (sr->numopsend == sr->numopsbegin) {
      sr->state       = STATE_BEGIN;
      sr->numopsend   = 0;
      sr->numopsbegin = 0;
    }
  }
  PetscFunctionReturn(0);
}

/*@
   BVDotColumn - Computes multiple dot products of a column against all the
   previous columns of a BV.

   Collective on X

   Input Parameters:
+  X - basis vectors
-  j - the column index

   Output Parameter:
.  q - an array where the result must be placed

   Notes:
   This operation is equivalent to BVDotVec() but it uses column j of X
   rather than taking a Vec as an argument. The number of active columns of
   X is set to j before the computation, and restored afterwards.
   If X has leading columns specified, then these columns do not participate
   in the computation. Therefore, the length of array q must be equal to j
   minus the number of leading columns.

   Developer Notes:
   If q is NULL, then the result is written in position nc+l of the internal
   buffer vector, see BVGetBufferVec().

   Level: advanced

.seealso: BVDot(), BVDotVec(), BVSetActiveColumns(), BVSetMatrix()
@*/
PetscErrorCode BVDotColumn(BV X,PetscInt j,PetscScalar *q)
{
  PetscErrorCode ierr;
  PetscInt       ksave;
  Vec            y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(X,j,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  BVCheckOp(X,1,dotvec);

  if (j<0) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  if (j>=X->m) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%D but BV only has %D columns",j,X->m);

  ierr = PetscLogEventBegin(BV_DotVec,X,0,0,0);CHKERRQ(ierr);
  ksave = X->k;
  X->k = j;
  if (!q && !X->buffer) { ierr = BVGetBufferVec(X,&X->buffer);CHKERRQ(ierr); }
  ierr = BVGetColumn(X,j,&y);CHKERRQ(ierr);
  ierr = (*X->ops->dotvec)(X,y,q);CHKERRQ(ierr);
  ierr = BVRestoreColumn(X,j,&y);CHKERRQ(ierr);
  X->k = ksave;
  ierr = PetscLogEventEnd(BV_DotVec,X,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVDotColumnBegin - Starts a split phase dot product computation.

   Input Parameters:
+  X - basis vectors
-  j - the column index
-  m - an array where the result will go (can be NULL)

   Note:
   Each call to BVDotColumnBegin() should be paired with a call to BVDotColumnEnd().

   Level: advanced

.seealso: BVDotColumnEnd(), BVDotColumn()
@*/
PetscErrorCode BVDotColumnBegin(BV X,PetscInt j,PetscScalar *m)
{
  PetscErrorCode      ierr;
  PetscInt            i,nv,ksave;
  PetscSplitReduction *sr;
  MPI_Comm            comm;
  Vec                 y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(X,j,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);

  if (j<0) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  if (j>=X->m) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%D but BV only has %D columns",j,X->m);
  ksave = X->k;
  X->k = j;
  ierr = BVGetColumn(X,j,&y);CHKERRQ(ierr);

  if (X->ops->dotvec_begin) {
    ierr = (*X->ops->dotvec_begin)(X,y,m);CHKERRQ(ierr);
  } else {
    BVCheckOp(X,1,dotvec_local);
    nv = X->k-X->l;
    ierr = PetscObjectGetComm((PetscObject)X,&comm);CHKERRQ(ierr);
    ierr = PetscSplitReductionGet(comm,&sr);CHKERRQ(ierr);
    if (sr->state != STATE_BEGIN) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ORDER,"Called before all VecxxxEnd() called");
    for (i=0;i<nv;i++) {
      if (sr->numopsbegin+i >= sr->maxops) {
        ierr = PetscSplitReductionExtend(sr);CHKERRQ(ierr);
      }
      sr->reducetype[sr->numopsbegin+i] = PETSC_SR_REDUCE_SUM;
      sr->invecs[sr->numopsbegin+i]     = (void*)X;
    }
    ierr = PetscLogEventBegin(BV_DotVec,X,0,0,0);CHKERRQ(ierr);
    ierr = (*X->ops->dotvec_local)(X,y,sr->lvalues+sr->numopsbegin);CHKERRQ(ierr);
    sr->numopsbegin += nv;
    ierr = PetscLogEventEnd(BV_DotVec,X,0,0,0);CHKERRQ(ierr);
  }
  ierr = BVRestoreColumn(X,j,&y);CHKERRQ(ierr);
  X->k = ksave;
  PetscFunctionReturn(0);
}

/*@
   BVDotColumnEnd - Ends a split phase dot product computation.

   Input Parameters:
+  X - basis vectors
.  j - the column index
-  m - an array where the result will go

   Notes:
   Each call to BVDotColumnBegin() should be paired with a call to BVDotColumnEnd().

   Level: advanced

.seealso: BVDotColumnBegin(), BVDotColumn()
@*/
PetscErrorCode BVDotColumnEnd(BV X,PetscInt j,PetscScalar *m)
{
  PetscErrorCode      ierr;
  PetscInt            i,nv,ksave;
  PetscSplitReduction *sr;
  MPI_Comm            comm;
  Vec                 y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(X,j,2);
  PetscValidType(X,1);
  BVCheckSizes(X,1);

  if (j<0) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  if (j>=X->m) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%D but BV only has %D columns",j,X->m);
  ksave = X->k;
  X->k = j;

  if (X->ops->dotvec_end) {
    ierr = BVGetColumn(X,j,&y);CHKERRQ(ierr);
    ierr = (*X->ops->dotvec_end)(X,y,m);CHKERRQ(ierr);
    ierr = BVRestoreColumn(X,j,&y);CHKERRQ(ierr);
  } else {
    nv = X->k-X->l;
    ierr = PetscObjectGetComm((PetscObject)X,&comm);CHKERRQ(ierr);
    ierr = PetscSplitReductionGet(comm,&sr);CHKERRQ(ierr);
    ierr = PetscSplitReductionEnd(sr);CHKERRQ(ierr);

    if (sr->numopsend >= sr->numopsbegin) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Called VecxxxEnd() more times than VecxxxBegin()");
    if ((void*)X != sr->invecs[sr->numopsend]) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Called BVxxxEnd() in a different order or with a different BV than BVxxxBegin()");
    if (sr->reducetype[sr->numopsend] != PETSC_SR_REDUCE_SUM) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_WRONGSTATE,"Wrong type of reduction");
    for (i=0;i<nv;i++) m[i] = sr->gvalues[sr->numopsend++];

    /* Finished getting all the results so reset to no outstanding requests */
    if (sr->numopsend == sr->numopsbegin) {
      sr->state       = STATE_BEGIN;
      sr->numopsend   = 0;
      sr->numopsbegin = 0;
    }
  }
  X->k = ksave;
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode BVNorm_Private(BV bv,Vec z,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  PetscScalar    p;

  PetscFunctionBegin;
  ierr = BV_IPMatMult(bv,z);CHKERRQ(ierr);
  ierr = VecDot(bv->Bx,z,&p);CHKERRQ(ierr);
  ierr = BV_SafeSqrt(bv,p,val);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode BVNorm_Begin_Private(BV bv,Vec z,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  PetscScalar    p;

  PetscFunctionBegin;
  ierr = BV_IPMatMult(bv,z);CHKERRQ(ierr);
  ierr = VecDotBegin(bv->Bx,z,&p);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode BVNorm_End_Private(BV bv,Vec z,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  PetscScalar    p;

  PetscFunctionBegin;
  ierr = VecDotEnd(bv->Bx,z,&p);CHKERRQ(ierr);
  ierr = BV_SafeSqrt(bv,p,val);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVNorm - Computes the matrix norm of the BV.

   Collective on bv

   Input Parameters:
+  bv   - basis vectors
-  type - the norm type

   Output Parameter:
.  val  - the norm

   Notes:
   All active columns (except the leading ones) are considered as a matrix.
   The allowed norms are NORM_1, NORM_FROBENIUS, and NORM_INFINITY.

   This operation fails if a non-standard inner product has been
   specified with BVSetMatrix().

   Level: intermediate

.seealso: BVNormVec(), BVNormColumn(), BVSetActiveColumns(), BVSetMatrix()
@*/
PetscErrorCode BVNorm(BV bv,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveEnum(bv,type,2);
  PetscValidRealPointer(val,3);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  if (type==NORM_2 || type==NORM_1_AND_2) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");
  if (bv->matrix) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Matrix norm not available for non-standard inner product");

  ierr = PetscLogEventBegin(BV_Norm,bv,0,0,0);CHKERRQ(ierr);
  ierr = (*bv->ops->norm)(bv,-1,type,val);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_Norm,bv,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVNormVec - Computes the norm of a given vector.

   Collective on bv

   Input Parameters:
+  bv   - basis vectors
.  v    - the vector
-  type - the norm type

   Output Parameter:
.  val  - the norm

   Notes:
   This is the analogue of BVNormColumn() but for a vector that is not in the BV.
   If a non-standard inner product has been specified with BVSetMatrix(),
   then the returned value is sqrt(v'*B*v), where B is the inner product
   matrix (argument 'type' is ignored). Otherwise, VecNorm() is called.

   Level: developer

.seealso: BVNorm(), BVNormColumn(), BVSetMatrix()
@*/
PetscErrorCode BVNormVec(BV bv,Vec v,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  PetscInt       n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(v,VEC_CLASSID,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscValidRealPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscValidType(v,2);
  PetscCheckSameComm(bv,1,v,2);

  if (type==NORM_1_AND_2 && !bv->matrix) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  ierr = PetscLogEventBegin(BV_NormVec,bv,0,0,0);CHKERRQ(ierr);
  if (bv->matrix) { /* non-standard inner product */
    ierr = VecGetLocalSize(v,&n);CHKERRQ(ierr);
    if (bv->n!=n) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension bv %D, v %D",bv->n,n);
    ierr = BVNorm_Private(bv,v,type,val);CHKERRQ(ierr);
  } else {
    ierr = VecNorm(v,type,val);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(BV_NormVec,bv,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVNormVecBegin - Starts a split phase norm computation.

   Input Parameters:
+  bv   - basis vectors
.  v    - the vector
.  type - the norm type
-  val  - the norm

   Note:
   Each call to BVNormVecBegin() should be paired with a call to BVNormVecEnd().

   Level: advanced

.seealso: BVNormVecEnd(), BVNormVec()
@*/
PetscErrorCode BVNormVecBegin(BV bv,Vec v,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  PetscInt       n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(v,VEC_CLASSID,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscValidRealPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscValidType(v,2);
  PetscCheckSameTypeAndComm(bv,1,v,2);

  if (type==NORM_1_AND_2 && !bv->matrix) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  ierr = PetscLogEventBegin(BV_NormVec,bv,0,0,0);CHKERRQ(ierr);
  if (bv->matrix) { /* non-standard inner product */
    ierr = VecGetLocalSize(v,&n);CHKERRQ(ierr);
    if (bv->n!=n) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension bv %D, v %D",bv->n,n);
    ierr = BVNorm_Begin_Private(bv,v,type,val);CHKERRQ(ierr);
  } else {
    ierr = VecNormBegin(v,type,val);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(BV_NormVec,bv,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVNormVecEnd - Ends a split phase norm computation.

   Input Parameters:
+  bv   - basis vectors
.  v    - the vector
.  type - the norm type
-  val  - the norm

   Note:
   Each call to BVNormVecBegin() should be paired with a call to BVNormVecEnd().

   Level: advanced

.seealso: BVNormVecBegin(), BVNormVec()
@*/
PetscErrorCode BVNormVecEnd(BV bv,Vec v,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscValidRealPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  if (type==NORM_1_AND_2) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  if (bv->matrix) { /* non-standard inner product */
    ierr = BVNorm_End_Private(bv,v,type,val);CHKERRQ(ierr);
  } else {
    ierr = VecNormEnd(v,type,val);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   BVNormColumn - Computes the vector norm of a selected column.

   Collective on bv

   Input Parameters:
+  bv   - basis vectors
.  j    - column number to be used
-  type - the norm type

   Output Parameter:
.  val  - the norm

   Notes:
   The norm of V[j] is computed (NORM_1, NORM_2, or NORM_INFINITY).
   If a non-standard inner product has been specified with BVSetMatrix(),
   then the returned value is sqrt(V[j]'*B*V[j]),
   where B is the inner product matrix (argument 'type' is ignored).

   Level: intermediate

.seealso: BVNorm(), BVNormVec(), BVSetActiveColumns(), BVSetMatrix()
@*/
PetscErrorCode BVNormColumn(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode ierr;
  Vec            z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscValidRealPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  if (j<0 || j>=bv->m) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Argument j has wrong value %D, the number of columns is %D",j,bv->m);
  if (type==NORM_1_AND_2 && !bv->matrix) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  ierr = PetscLogEventBegin(BV_NormVec,bv,0,0,0);CHKERRQ(ierr);
  if (bv->matrix) { /* non-standard inner product */
    ierr = BVGetColumn(bv,j,&z);CHKERRQ(ierr);
    ierr = BVNorm_Private(bv,z,type,val);CHKERRQ(ierr);
    ierr = BVRestoreColumn(bv,j,&z);CHKERRQ(ierr);
  } else {
    ierr = (*bv->ops->norm)(bv,j,type,val);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(BV_NormVec,bv,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVNormColumnBegin - Starts a split phase norm computation.

   Input Parameters:
+  bv   - basis vectors
.  j    - column number to be used
.  type - the norm type
-  val  - the norm

   Note:
   Each call to BVNormColumnBegin() should be paired with a call to BVNormColumnEnd().

   Level: advanced

.seealso: BVNormColumnEnd(), BVNormColumn()
@*/
PetscErrorCode BVNormColumnBegin(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode      ierr;
  PetscSplitReduction *sr;
  PetscReal           lresult;
  MPI_Comm            comm;
  Vec                 z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscValidRealPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  if (j<0 || j>=bv->m) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Argument j has wrong value %D, the number of columns is %D",j,bv->m);
  if (type==NORM_1_AND_2 && !bv->matrix) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  ierr = PetscLogEventBegin(BV_NormVec,bv,0,0,0);CHKERRQ(ierr);
  ierr = BVGetColumn(bv,j,&z);CHKERRQ(ierr);
  if (bv->matrix) { /* non-standard inner product */
    ierr = BVNorm_Begin_Private(bv,z,type,val);CHKERRQ(ierr);
  } else if (bv->ops->norm_begin) {
    ierr = (*bv->ops->norm_begin)(bv,j,type,val);CHKERRQ(ierr);
  } else {
    BVCheckOp(bv,1,norm_local);
    ierr = PetscObjectGetComm((PetscObject)z,&comm);CHKERRQ(ierr);
    ierr = PetscSplitReductionGet(comm,&sr);CHKERRQ(ierr);
    if (sr->state != STATE_BEGIN) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ORDER,"Called before all VecxxxEnd() called");
    if (sr->numopsbegin >= sr->maxops) {
      ierr = PetscSplitReductionExtend(sr);CHKERRQ(ierr);
    }
    sr->invecs[sr->numopsbegin] = (void*)bv;
    ierr = (*bv->ops->norm_local)(bv,j,type,&lresult);CHKERRQ(ierr);
    if (type == NORM_2) lresult = lresult*lresult;
    if (type == NORM_MAX) sr->reducetype[sr->numopsbegin] = PETSC_SR_REDUCE_MAX;
    else sr->reducetype[sr->numopsbegin] = PETSC_SR_REDUCE_SUM;
    sr->lvalues[sr->numopsbegin++] = lresult;
  }
  ierr = BVRestoreColumn(bv,j,&z);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_NormVec,bv,0,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVNormColumnEnd - Ends a split phase norm computation.

   Input Parameters:
+  bv   - basis vectors
.  j    - column number to be used
.  type - the norm type
-  val  - the norm

   Note:
   Each call to BVNormColumnBegin() should be paired with a call to BVNormColumnEnd().

   Level: advanced

.seealso: BVNormColumnBegin(), BVNormColumn()
@*/
PetscErrorCode BVNormColumnEnd(BV bv,PetscInt j,NormType type,PetscReal *val)
{
  PetscErrorCode      ierr;
  PetscSplitReduction *sr;
  MPI_Comm            comm;
  Vec                 z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidLogicalCollectiveEnum(bv,type,3);
  PetscValidRealPointer(val,4);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  if (type==NORM_1_AND_2) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Requested norm not available");

  ierr = BVGetColumn(bv,j,&z);CHKERRQ(ierr);
  if (bv->matrix) { /* non-standard inner product */
    ierr = BVNorm_End_Private(bv,z,type,val);CHKERRQ(ierr);
  } else if (bv->ops->norm_end) {
    ierr = (*bv->ops->norm_end)(bv,j,type,val);CHKERRQ(ierr);
  } else {
    ierr = PetscObjectGetComm((PetscObject)z,&comm);CHKERRQ(ierr);
    ierr = PetscSplitReductionGet(comm,&sr);CHKERRQ(ierr);
    ierr = PetscSplitReductionEnd(sr);CHKERRQ(ierr);

    if (sr->numopsend >= sr->numopsbegin) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Called VecxxxEnd() more times then VecxxxBegin()");
    if ((void*)bv != sr->invecs[sr->numopsend]) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Called VecxxxEnd() in a different order or with a different vector than VecxxxBegin()");
    if (sr->reducetype[sr->numopsend] != PETSC_SR_REDUCE_MAX && type == NORM_MAX) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Called BVNormEnd(,NORM_MAX,) on a reduction started with VecDotBegin() or NORM_1 or NORM_2");
    *val = PetscRealPart(sr->gvalues[sr->numopsend++]);
    if (type == NORM_2) *val = PetscSqrtReal(*val);
    if (sr->numopsend == sr->numopsbegin) {
      sr->state       = STATE_BEGIN;
      sr->numopsend   = 0;
      sr->numopsbegin = 0;
    }
  }
  ierr = BVRestoreColumn(bv,j,&z);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Compute Y^H*A*X: right part column by column (with MatMult) and bottom
  part row by row (with MatMultHermitianTranspose); result placed in marray[*,ldm]
*/
PETSC_STATIC_INLINE PetscErrorCode BVMatProject_Vec(BV X,Mat A,BV Y,PetscScalar *marray,PetscInt ldm,PetscBool symm)
{
  PetscErrorCode ierr;
  PetscInt       i,j,lx,ly,kx,ky,ulim;
  Vec            z,f;

  PetscFunctionBegin;
  lx = X->l; kx = X->k;
  ly = Y->l; ky = Y->k;
  ierr = BVCreateVec(X,&f);CHKERRQ(ierr);
  BVCheckOp(Y,3,dotvec);
  for (j=lx;j<kx;j++) {
    ierr = BVGetColumn(X,j,&z);CHKERRQ(ierr);
    ierr = MatMult(A,z,f);CHKERRQ(ierr);
    ierr = BVRestoreColumn(X,j,&z);CHKERRQ(ierr);
    ulim = PetscMin(ly+(j-lx)+1,ky);
    Y->l = 0; Y->k = ulim;
    ierr = (*Y->ops->dotvec)(Y,f,marray+j*ldm);CHKERRQ(ierr);
    if (symm) {
      for (i=0;i<j;i++) marray[j+i*ldm] = PetscConj(marray[i+j*ldm]);
    }
  }
  if (!symm) {
    BVCheckOp(X,1,dotvec);
    ierr = BV_AllocateCoeffs(Y);CHKERRQ(ierr);
    for (j=ly;j<ky;j++) {
      ierr = BVGetColumn(Y,j,&z);CHKERRQ(ierr);
      ierr = MatMultHermitianTranspose(A,z,f);CHKERRQ(ierr);
      ierr = BVRestoreColumn(Y,j,&z);CHKERRQ(ierr);
      ulim = PetscMin(lx+(j-ly),kx);
      X->l = 0; X->k = ulim;
      ierr = (*X->ops->dotvec)(X,f,Y->h);CHKERRQ(ierr);
      for (i=0;i<ulim;i++) marray[j+i*ldm] = PetscConj(Y->h[i]);
    }
  }
  ierr = VecDestroy(&f);CHKERRQ(ierr);
  X->l = lx; X->k = kx;
  Y->l = ly; Y->k = ky;
  PetscFunctionReturn(0);
}

/*
  Compute Y^H*A*X= [   --   | Y0'*W1 ]
                   [ Y1'*W0 | Y1'*W1 ]
  Allocates auxiliary BV to store the result of A*X, then one BVDot
  call for top-right part and another one for bottom part;
  result placed in marray[*,ldm]
*/
PETSC_STATIC_INLINE PetscErrorCode BVMatProject_MatMult(BV X,Mat A,BV Y,PetscScalar *marray,PetscInt ldm)
{
  PetscErrorCode ierr;
  PetscInt       j,lx,ly,kx,ky;
  PetscScalar    *harray;
  Mat            H;
  BV             W;

  PetscFunctionBegin;
  lx = X->l; kx = X->k;
  ly = Y->l; ky = Y->k;
  ierr = BVDuplicate(X,&W);CHKERRQ(ierr);
  X->l = 0; X->k = kx;
  W->l = 0; W->k = kx;
  ierr = BVMatMult(X,A,W);CHKERRQ(ierr);

  /* top-right part, Y0'*AX1 */
  if (ly>0 && lx<kx) {
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,ly,kx,NULL,&H);CHKERRQ(ierr);
    W->l = lx; W->k = kx;
    Y->l = 0;  Y->k = ly;
    ierr = BVDot(W,Y,H);CHKERRQ(ierr);
    ierr = MatDenseGetArray(H,&harray);CHKERRQ(ierr);
    for (j=lx;j<kx;j++) {
      ierr = PetscArraycpy(marray+j*ldm,harray+j*ly,ly);CHKERRQ(ierr);
    }
    ierr = MatDenseRestoreArray(H,&harray);CHKERRQ(ierr);
    ierr = MatDestroy(&H);CHKERRQ(ierr);
  }

  /* bottom part, Y1'*AX */
  if (kx>0 && ly<ky) {
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,ky,kx,NULL,&H);CHKERRQ(ierr);
    W->l = 0;  W->k = kx;
    Y->l = ly; Y->k = ky;
    ierr = BVDot(W,Y,H);CHKERRQ(ierr);
    ierr = MatDenseGetArray(H,&harray);CHKERRQ(ierr);
    for (j=0;j<kx;j++) {
      ierr = PetscArraycpy(marray+j*ldm+ly,harray+j*ky+ly,ky-ly);CHKERRQ(ierr);
    }
    ierr = MatDenseRestoreArray(H,&harray);CHKERRQ(ierr);
    ierr = MatDestroy(&H);CHKERRQ(ierr);
  }
  ierr = BVDestroy(&W);CHKERRQ(ierr);
  X->l = lx; X->k = kx;
  Y->l = ly; Y->k = ky;
  PetscFunctionReturn(0);
}

/*
  Compute Y^H*A*X= [   --   | Y0'*W1 ]
                   [ Y1'*W0 | Y1'*W1 ]
  First stage: allocate auxiliary BV to store A*X1, one BVDot for right part;
  Second stage: resize BV to accomodate A'*Y1, then call BVDot for transpose of
  bottom-left part; result placed in marray[*,ldm]
*/
PETSC_STATIC_INLINE PetscErrorCode BVMatProject_MatMult_2(BV X,Mat A,BV Y,PetscScalar *marray,PetscInt ldm,PetscBool symm)
{
  PetscErrorCode ierr;
  PetscInt       i,j,lx,ly,kx,ky;
  PetscScalar    *harray;
  Mat            H;
  BV             W;

  PetscFunctionBegin;
  lx = X->l; kx = X->k;
  ly = Y->l; ky = Y->k;

  /* right part, Y'*AX1 */
  ierr = BVDuplicateResize(X,kx-lx,&W);CHKERRQ(ierr);
  if (ky>0 && lx<kx) {
    ierr = BVMatMult(X,A,W);CHKERRQ(ierr);
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,ky,kx-lx,NULL,&H);CHKERRQ(ierr);
    Y->l = 0; Y->k = ky;
    ierr = BVDot(W,Y,H);CHKERRQ(ierr);
    ierr = MatDenseGetArray(H,&harray);CHKERRQ(ierr);
    for (j=lx;j<kx;j++) {
      ierr = PetscArraycpy(marray+j*ldm,harray+(j-lx)*ky,ky);CHKERRQ(ierr);
    }
    ierr = MatDenseRestoreArray(H,&harray);CHKERRQ(ierr);
    ierr = MatDestroy(&H);CHKERRQ(ierr);
  }

  /* bottom-left part, Y1'*AX0 */
  if (lx>0 && ly<ky) {
    if (symm) {
      /* do not compute, just copy symmetric elements */
      for (i=ly;i<ky;i++) {
        for (j=0;j<lx;j++) marray[i+j*ldm] = PetscConj(marray[j+i*ldm]);
      }
    } else {
      ierr = BVResize(W,ky-ly,PETSC_FALSE);CHKERRQ(ierr);
      Y->l = ly; Y->k = ky;
      ierr = BVMatMultHermitianTranspose(Y,A,W);CHKERRQ(ierr);
      ierr = MatCreateSeqDense(PETSC_COMM_SELF,lx,ky-ly,NULL,&H);CHKERRQ(ierr);
      X->l = 0; X->k = lx;
      ierr = BVDot(W,X,H);CHKERRQ(ierr);
      ierr = MatDenseGetArray(H,&harray);CHKERRQ(ierr);
      for (i=0;i<ky-ly;i++) {
        for (j=0;j<lx;j++) {
          marray[i+j*ldm+ly] = PetscConj(harray[j+i*(ky-ly)]);
        }
      }
      ierr = MatDenseRestoreArray(H,&harray);CHKERRQ(ierr);
      ierr = MatDestroy(&H);CHKERRQ(ierr);
    }
  }
  ierr = BVDestroy(&W);CHKERRQ(ierr);
  X->l = lx; X->k = kx;
  Y->l = ly; Y->k = ky;
  PetscFunctionReturn(0);
}

/*
  Compute Y^H*X = [   --   | Y0'*X1 ]     (X contains A*X):
                  [ Y1'*X0 | Y1'*X1 ]
  one BVDot call for top-right part and another one for bottom part;
  result placed in marray[*,ldm]
*/
PETSC_STATIC_INLINE PetscErrorCode BVMatProject_Dot(BV X,BV Y,PetscScalar *marray,PetscInt ldm)
{
  PetscErrorCode ierr;
  PetscInt       j,lx,ly,kx,ky;
  PetscScalar    *harray;
  Mat            H;

  PetscFunctionBegin;
  lx = X->l; kx = X->k;
  ly = Y->l; ky = Y->k;

  /* top-right part, Y0'*X1 */
  if (ly>0 && lx<kx) {
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,ly,kx,NULL,&H);CHKERRQ(ierr);
    X->l = lx; X->k = kx;
    Y->l = 0;  Y->k = ly;
    ierr = BVDot(X,Y,H);CHKERRQ(ierr);
    ierr = MatDenseGetArray(H,&harray);CHKERRQ(ierr);
    for (j=lx;j<kx;j++) {
      ierr = PetscArraycpy(marray+j*ldm,harray+j*ly,ly);CHKERRQ(ierr);
    }
    ierr = MatDenseRestoreArray(H,&harray);CHKERRQ(ierr);
    ierr = MatDestroy(&H);CHKERRQ(ierr);
  }

  /* bottom part, Y1'*X */
  if (kx>0 && ly<ky) {
    ierr = MatCreateSeqDense(PETSC_COMM_SELF,ky,kx,NULL,&H);CHKERRQ(ierr);
    X->l = 0;  X->k = kx;
    Y->l = ly; Y->k = ky;
    ierr = BVDot(X,Y,H);CHKERRQ(ierr);
    ierr = MatDenseGetArray(H,&harray);CHKERRQ(ierr);
    for (j=0;j<kx;j++) {
      ierr = PetscArraycpy(marray+j*ldm+ly,harray+j*ky+ly,ky-ly);CHKERRQ(ierr);
    }
    ierr = MatDenseRestoreArray(H,&harray);CHKERRQ(ierr);
    ierr = MatDestroy(&H);CHKERRQ(ierr);
  }
  X->l = lx; X->k = kx;
  Y->l = ly; Y->k = ky;
  PetscFunctionReturn(0);
}

/*@
   BVMatProject - Computes the projection of a matrix onto a subspace.

   Collective on X

   Input Parameters:
+  X - basis vectors
.  A - (optional) matrix to be projected
.  Y - left basis vectors, can be equal to X
-  M - Mat object where the result must be placed

   Output Parameter:
.  M - the resulting matrix

   Notes:
   If A=NULL, then it is assumed that X already contains A*X.

   This operation is similar to BVDot(), with important differences.
   The goal is to compute the matrix resulting from the orthogonal projection
   of A onto the subspace spanned by the columns of X, M = X^H*A*X, or the
   oblique projection onto X along Y, M = Y^H*A*X.

   A difference with respect to BVDot() is that the standard inner product
   is always used, regardless of a non-standard inner product being specified
   with BVSetMatrix().

   On entry, M must be a sequential dense Mat with dimensions ky,kx at least,
   where ky (resp. kx) is the number of active columns of Y (resp. X).
   Another difference with respect to BVDot() is that all entries of M are
   computed except the leading ly,lx part, where ly (resp. lx) is the
   number of leading columns of Y (resp. X). Hence, the leading columns of
   X and Y participate in the computation, as opposed to BVDot().
   The leading part of M is assumed to be already available from previous
   computations.

   In the orthogonal projection case, Y=X, some computation can be saved if
   A is real symmetric (or complex Hermitian). In order to exploit this
   property, the symmetry flag of A must be set with MatSetOption().

   Level: intermediate

.seealso: BVDot(), BVSetActiveColumns(), BVSetMatrix()
@*/
PetscErrorCode BVMatProject(BV X,Mat A,BV Y,Mat M)
{
  PetscErrorCode ierr;
  PetscBool      match,set,flg,symm=PETSC_FALSE;
  PetscInt       m,n;
  PetscScalar    *marray;
  Mat            Xmatrix,Ymatrix;
  PetscObjectId  idx,idy;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  if (A) PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidHeaderSpecific(Y,BV_CLASSID,3);
  PetscValidHeaderSpecific(M,MAT_CLASSID,4);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  if (A) {
    PetscValidType(A,2);
    PetscCheckSameComm(X,1,A,2);
  }
  PetscValidType(Y,3);
  BVCheckSizes(Y,3);
  PetscValidType(M,4);
  PetscCheckSameTypeAndComm(X,1,Y,3);
  ierr = PetscObjectTypeCompare((PetscObject)M,MATSEQDENSE,&match);CHKERRQ(ierr);
  if (!match) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_SUP,"Matrix M must be of type seqdense");

  ierr = MatGetSize(M,&m,&n);CHKERRQ(ierr);
  if (m<Y->k) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Matrix M has %D rows, should have at least %D",m,Y->k);
  if (n<X->k) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_SIZ,"Matrix M has %D columns, should have at least %D",n,X->k);
  if (X->n!=Y->n) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %D, Y %D",X->n,Y->n);

  ierr = PetscLogEventBegin(BV_MatProject,X,A,Y,0);CHKERRQ(ierr);
  /* temporarily set standard inner product */
  Xmatrix = X->matrix;
  Ymatrix = Y->matrix;
  X->matrix = Y->matrix = NULL;

  ierr = PetscObjectGetId((PetscObject)X,&idx);CHKERRQ(ierr);
  ierr = PetscObjectGetId((PetscObject)Y,&idy);CHKERRQ(ierr);
  if (A && idx==idy) { /* check symmetry of M=X'AX */
    ierr = MatIsHermitianKnown(A,&set,&flg);CHKERRQ(ierr);
    symm = set? flg: PETSC_FALSE;
  }

  ierr = MatDenseGetArray(M,&marray);CHKERRQ(ierr);

  if (A) {
    if (X->vmm==BV_MATMULT_VECS) {
      /* perform computation column by column */
      ierr = BVMatProject_Vec(X,A,Y,marray,m,symm);CHKERRQ(ierr);
    } else {
      /* use BVMatMult, then BVDot */
      ierr = MatHasOperation(A,MATOP_MULT_TRANSPOSE,&flg);CHKERRQ(ierr);
      if (symm || (flg && X->l>=X->k/2 && Y->l>=Y->k/2)) {
        ierr = BVMatProject_MatMult_2(X,A,Y,marray,m,symm);CHKERRQ(ierr);
      } else {
        ierr = BVMatProject_MatMult(X,A,Y,marray,m);CHKERRQ(ierr);
      }
    }
  } else {
    /* use BVDot on subblocks */
    ierr = BVMatProject_Dot(X,Y,marray,m);CHKERRQ(ierr);
  }

  ierr = MatDenseRestoreArray(M,&marray);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_MatProject,X,A,Y,0);CHKERRQ(ierr);
  /* restore non-standard inner product */
  X->matrix = Xmatrix;
  Y->matrix = Ymatrix;
  PetscFunctionReturn(0);
}

