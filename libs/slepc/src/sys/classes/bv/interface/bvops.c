/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BV operations, except those involving global communication
*/

#include <slepc/private/bvimpl.h>      /*I "slepcbv.h" I*/
#include <slepcds.h>

/*@
   BVMult - Computes Y = beta*Y + alpha*X*Q.

   Logically Collective on Y

   Input Parameters:
+  Y,X        - basis vectors
.  alpha,beta - scalars
-  Q          - (optional) sequential dense matrix

   Output Parameter:
.  Y          - the modified basis vectors

   Notes:
   X and Y must be different objects. The case X=Y can be addressed with
   BVMultInPlace().

   If matrix Q is NULL, then an AXPY operation Y = beta*Y + alpha*X is done
   (i.e. results as if Q = identity). If provided,
   the matrix Q must be a sequential dense Mat, with all entries equal on
   all processes (otherwise each process will compute a different update).
   The dimensions of Q must be at least m,n where m is the number of active
   columns of X and n is the number of active columns of Y.

   The leading columns of Y are not modified. Also, if X has leading
   columns specified, then these columns do not participate in the computation.
   Hence, only rows (resp. columns) of Q starting from lx (resp. ly) are used,
   where lx (resp. ly) is the number of leading columns of X (resp. Y).

   Level: intermediate

.seealso: BVMultVec(), BVMultColumn(), BVMultInPlace(), BVSetActiveColumns()
@*/
PetscErrorCode BVMult(BV Y,PetscScalar alpha,PetscScalar beta,BV X,Mat Q)
{
  PetscErrorCode ierr;
  PetscBool      match;
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(Y,BV_CLASSID,1);
  PetscValidLogicalCollectiveScalar(Y,alpha,2);
  PetscValidLogicalCollectiveScalar(Y,beta,3);
  PetscValidHeaderSpecific(X,BV_CLASSID,4);
  if (Q) PetscValidHeaderSpecific(Q,MAT_CLASSID,5);
  PetscValidType(Y,1);
  BVCheckSizes(Y,1);
  BVCheckOp(Y,1,mult);
  PetscValidType(X,4);
  BVCheckSizes(X,4);
  if (Q) PetscValidType(Q,5);
  PetscCheckSameTypeAndComm(Y,1,X,4);
  if (X==Y) SETERRQ(PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_WRONG,"X and Y arguments must be different");
  if (Q) {
    ierr = PetscObjectTypeCompare((PetscObject)Q,MATSEQDENSE,&match);CHKERRQ(ierr);
    if (!match) SETERRQ(PetscObjectComm((PetscObject)Y),PETSC_ERR_SUP,"Mat argument must be of type seqdense");
    ierr = MatGetSize(Q,&m,&n);CHKERRQ(ierr);
    if (m<X->k) SETERRQ2(PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_SIZ,"Mat argument has %D rows, should have at least %D",m,X->k);
    if (n<Y->k) SETERRQ2(PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_SIZ,"Mat argument has %D columns, should have at least %D",n,Y->k);
  }
  if (X->n!=Y->n) SETERRQ2(PetscObjectComm((PetscObject)Y),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension X %D, Y %D",X->n,Y->n);

  ierr = PetscLogEventBegin(BV_Mult,X,Y,0,0);CHKERRQ(ierr);
  ierr = (*Y->ops->mult)(Y,alpha,beta,X,Q);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_Mult,X,Y,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)Y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMultVec - Computes y = beta*y + alpha*X*q.

   Logically Collective on X

   Input Parameters:
+  X          - a basis vectors object
.  alpha,beta - scalars
.  y          - a vector
-  q          - an array of scalars

   Output Parameter:
.  y          - the modified vector

   Notes:
   This operation is the analogue of BVMult() but with a BV and a Vec,
   instead of two BV. Note that arguments are listed in different order
   with respect to BVMult().

   If X has leading columns specified, then these columns do not participate
   in the computation.

   The length of array q must be equal to the number of active columns of X
   minus the number of leading columns, i.e. the first entry of q multiplies
   the first non-leading column.

   Level: intermediate

.seealso: BVMult(), BVMultColumn(), BVMultInPlace(), BVSetActiveColumns()
@*/
PetscErrorCode BVMultVec(BV X,PetscScalar alpha,PetscScalar beta,Vec y,PetscScalar q[])
{
  PetscErrorCode ierr;
  PetscInt       n,N;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveScalar(X,alpha,2);
  PetscValidLogicalCollectiveScalar(X,beta,3);
  PetscValidHeaderSpecific(y,VEC_CLASSID,4);
  PetscValidScalarPointer(q,5);
  PetscValidType(X,1);
  BVCheckSizes(X,1);
  BVCheckOp(X,1,multvec);
  PetscValidType(y,4);
  PetscCheckSameComm(X,1,y,4);

  ierr = VecGetSize(y,&N);CHKERRQ(ierr);
  ierr = VecGetLocalSize(y,&n);CHKERRQ(ierr);
  if (N!=X->N || n!=X->n) SETERRQ4(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_INCOMP,"Vec sizes (global %D, local %D) do not match BV sizes (global %D, local %D)",N,n,X->N,X->n);

  ierr = PetscLogEventBegin(BV_MultVec,X,y,0,0);CHKERRQ(ierr);
  ierr = (*X->ops->multvec)(X,alpha,beta,y,q);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_MultVec,X,y,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMultColumn - Computes y = beta*y + alpha*X*q, where y is the j-th column
   of X.

   Logically Collective on X

   Input Parameters:
+  X          - a basis vectors object
.  alpha,beta - scalars
.  j          - the column index
-  q          - an array of scalars

   Notes:
   This operation is equivalent to BVMultVec() but it uses column j of X
   rather than taking a Vec as an argument. The number of active columns of
   X is set to j before the computation, and restored afterwards.
   If X has leading columns specified, then these columns do not participate
   in the computation. Therefore, the length of array q must be equal to j
   minus the number of leading columns.

   Developer Notes:
   If q is NULL, then the coefficients are taken from position nc+l of the
   internal buffer vector, see BVGetBufferVec().

   Level: advanced

.seealso: BVMult(), BVMultVec(), BVMultInPlace(), BVSetActiveColumns()
@*/
PetscErrorCode BVMultColumn(BV X,PetscScalar alpha,PetscScalar beta,PetscInt j,PetscScalar *q)
{
  PetscErrorCode ierr;
  PetscInt       ksave;
  Vec            y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  PetscValidLogicalCollectiveScalar(X,alpha,2);
  PetscValidLogicalCollectiveScalar(X,beta,3);
  PetscValidLogicalCollectiveInt(X,j,4);
  PetscValidType(X,1);
  BVCheckSizes(X,1);

  if (j<0) SETERRQ(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  if (j>=X->m) SETERRQ2(PetscObjectComm((PetscObject)X),PETSC_ERR_ARG_OUTOFRANGE,"Index j=%D but BV only has %D columns",j,X->m);

  ierr = PetscLogEventBegin(BV_MultVec,X,0,0,0);CHKERRQ(ierr);
  ksave = X->k;
  X->k = j;
  if (!q && !X->buffer) { ierr = BVGetBufferVec(X,&X->buffer);CHKERRQ(ierr); }
  ierr = BVGetColumn(X,j,&y);CHKERRQ(ierr);
  ierr = (*X->ops->multvec)(X,alpha,beta,y,q);CHKERRQ(ierr);
  ierr = BVRestoreColumn(X,j,&y);CHKERRQ(ierr);
  X->k = ksave;
  ierr = PetscLogEventEnd(BV_MultVec,X,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMultInPlace - Update a set of vectors as V(:,s:e-1) = V*Q(:,s:e-1).

   Logically Collective on V

   Input Parameters:
+  Q - a sequential dense matrix
.  s - first column of V to be overwritten
-  e - first column of V not to be overwritten

   Input/Output Parameter:
.  V - basis vectors

   Notes:
   The matrix Q must be a sequential dense Mat, with all entries equal on
   all processes (otherwise each process will compute a different update).

   This function computes V(:,s:e-1) = V*Q(:,s:e-1), that is, given a set of
   vectors V, columns from s to e-1 are overwritten with columns from s to
   e-1 of the matrix-matrix product V*Q. Only columns s to e-1 of Q are
   referenced.

   Level: intermediate

.seealso: BVMult(), BVMultVec(), BVMultInPlaceTranspose(), BVSetActiveColumns()
@*/
PetscErrorCode BVMultInPlace(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  PetscBool      match;
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidHeaderSpecific(Q,MAT_CLASSID,2);
  PetscValidLogicalCollectiveInt(V,s,3);
  PetscValidLogicalCollectiveInt(V,e,4);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidType(Q,2);
  ierr = PetscObjectTypeCompare((PetscObject)Q,MATSEQDENSE,&match);CHKERRQ(ierr);
  if (!match) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Mat argument must be of type seqdense");

  if (s<V->l || s>V->m) SETERRQ3(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument s has wrong value %D, should be between %D and %D",s,V->l,V->m);
  if (e<V->l || e>V->m) SETERRQ3(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument e has wrong value %D, should be between %D and %D",e,V->l,V->m);
  ierr = MatGetSize(Q,&m,&n);CHKERRQ(ierr);
  if (m<V->k) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat argument has %D rows, should have at least %D",m,V->k);
  if (e>n) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat argument only has %D columns, the requested value of e is larger: %D",n,e);
  if (s>=e) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(BV_MultInPlace,V,Q,0,0);CHKERRQ(ierr);
  ierr = (*V->ops->multinplace)(V,Q,s,e);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_MultInPlace,V,Q,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMultInPlaceTranspose - Update a set of vectors as V(:,s:e-1) = V*Q'(:,s:e-1).

   Logically Collective on V

   Input Parameters:
+  Q - a sequential dense matrix
.  s - first column of V to be overwritten
-  e - first column of V not to be overwritten

   Input/Output Parameter:
.  V - basis vectors

   Notes:
   This is a variant of BVMultInPlace() where the conjugate transpose
   of Q is used.

   Level: intermediate

.seealso: BVMultInPlace()
@*/
PetscErrorCode BVMultInPlaceTranspose(BV V,Mat Q,PetscInt s,PetscInt e)
{
  PetscErrorCode ierr;
  PetscBool      match;
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidHeaderSpecific(Q,MAT_CLASSID,2);
  PetscValidLogicalCollectiveInt(V,s,3);
  PetscValidLogicalCollectiveInt(V,e,4);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidType(Q,2);
  ierr = PetscObjectTypeCompare((PetscObject)Q,MATSEQDENSE,&match);CHKERRQ(ierr);
  if (!match) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Mat argument must be of type seqdense");

  if (s<V->l || s>V->m) SETERRQ3(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument s has wrong value %D, should be between %D and %D",s,V->l,V->m);
  if (e<V->l || e>V->m) SETERRQ3(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Argument e has wrong value %D, should be between %D and %D",e,V->l,V->m);
  ierr = MatGetSize(Q,&m,&n);CHKERRQ(ierr);
  if (n<V->k) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat argument has %D columns, should have at least %D",n,V->k);
  if (e>m) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Mat argument only has %D rows, the requested value of e is larger: %D",m,e);
  if (s>=e || !V->n) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(BV_MultInPlace,V,Q,0,0);CHKERRQ(ierr);
  ierr = (*V->ops->multinplacetrans)(V,Q,s,e);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_MultInPlace,V,Q,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVScale - Multiply the BV entries by a scalar value.

   Logically Collective on bv

   Input Parameters:
+  bv    - basis vectors
-  alpha - scaling factor

   Note:
   All active columns (except the leading ones) are scaled.

   Level: intermediate

.seealso: BVScaleColumn(), BVSetActiveColumns()
@*/
PetscErrorCode BVScale(BV bv,PetscScalar alpha)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveScalar(bv,alpha,2);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  if (alpha == (PetscScalar)1.0) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(BV_Scale,bv,0,0,0);CHKERRQ(ierr);
  if (bv->n) {
    ierr = (*bv->ops->scale)(bv,-1,alpha);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(BV_Scale,bv,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVScaleColumn - Scale one column of a BV.

   Logically Collective on bv

   Input Parameters:
+  bv    - basis vectors
.  j     - column number to be scaled
-  alpha - scaling factor

   Level: intermediate

.seealso: BVScale(), BVSetActiveColumns()
@*/
PetscErrorCode BVScaleColumn(BV bv,PetscInt j,PetscScalar alpha)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidLogicalCollectiveScalar(bv,alpha,3);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  if (j<0 || j>=bv->m) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Argument j has wrong value %D, the number of columns is %D",j,bv->m);
  if (alpha == (PetscScalar)1.0) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(BV_Scale,bv,0,0,0);CHKERRQ(ierr);
  if (bv->n) {
    ierr = (*bv->ops->scale)(bv,j,alpha);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(BV_Scale,bv,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode BVSetRandomColumn_Private(BV bv,PetscInt k)
{
  PetscErrorCode ierr;
  PetscInt       i,low,high;
  PetscScalar    *px,t;
  Vec            x;

  PetscFunctionBegin;
  ierr = BVGetColumn(bv,k,&x);CHKERRQ(ierr);
  if (bv->rrandom) {  /* generate the same vector irrespective of number of processes */
    ierr = VecGetOwnershipRange(x,&low,&high);CHKERRQ(ierr);
    ierr = VecGetArray(x,&px);CHKERRQ(ierr);
    for (i=0;i<bv->N;i++) {
      ierr = PetscRandomGetValue(bv->rand,&t);CHKERRQ(ierr);
      if (i>=low && i<high) px[i-low] = t;
    }
    ierr = VecRestoreArray(x,&px);CHKERRQ(ierr);
  } else {
    ierr = VecSetRandom(x,bv->rand);CHKERRQ(ierr);
  }
  ierr = BVRestoreColumn(bv,k,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVSetRandom - Set the columns of a BV to random numbers.

   Logically Collective on bv

   Input Parameters:
.  bv - basis vectors

   Note:
   All active columns (except the leading ones) are modified.

   Level: advanced

.seealso: BVSetRandomContext(), BVSetRandomColumn(), BVSetRandomCond(), BVSetActiveColumns()
@*/
PetscErrorCode BVSetRandom(BV bv)
{
  PetscErrorCode ierr;
  PetscInt       k;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  ierr = BVGetRandomContext(bv,&bv->rand);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(BV_SetRandom,bv,0,0,0);CHKERRQ(ierr);
  for (k=bv->l;k<bv->k;k++) {
    ierr = BVSetRandomColumn_Private(bv,k);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(BV_SetRandom,bv,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVSetRandomColumn - Set one column of a BV to random numbers.

   Logically Collective on bv

   Input Parameters:
+  bv - basis vectors
-  j  - column number to be set

   Level: advanced

.seealso: BVSetRandomContext(), BVSetRandom(), BVSetActiveColumns()
@*/
PetscErrorCode BVSetRandomColumn(BV bv,PetscInt j)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  if (j<0 || j>=bv->m) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Argument j has wrong value %D, the number of columns is %D",j,bv->m);

  ierr = BVGetRandomContext(bv,&bv->rand);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(BV_SetRandom,bv,0,0,0);CHKERRQ(ierr);
  ierr = BVSetRandomColumn_Private(bv,j);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_SetRandom,bv,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVSetRandomCond - Set the columns of a BV to random numbers, in a way that
   the generated matrix has a given condition number.

   Logically Collective on bv

   Input Parameters:
+  bv    - basis vectors
-  condn - condition number

   Note:
   All active columns (except the leading ones) are modified.

   Level: advanced

.seealso: BVSetRandomContext(), BVSetRandomColumn(), BVSetActiveColumns()
@*/
PetscErrorCode BVSetRandomCond(BV bv,PetscReal condn)
{
  PetscErrorCode ierr;
  PetscInt       k,i;
  PetscScalar    *eig,*d;
  DS             ds;
  Mat            A,X,Xt,M,G;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);

  ierr = BVGetRandomContext(bv,&bv->rand);CHKERRQ(ierr);
  ierr = PetscLogEventBegin(BV_SetRandom,bv,0,0,0);CHKERRQ(ierr);
  /* B = rand(n,k) */
  for (k=bv->l;k<bv->k;k++) {
    ierr = BVSetRandomColumn_Private(bv,k);CHKERRQ(ierr);
  }
  ierr = DSCreate(PetscObjectComm((PetscObject)bv),&ds);CHKERRQ(ierr);
  ierr = DSSetType(ds,DSHEP);CHKERRQ(ierr);
  ierr = DSAllocate(ds,bv->m);CHKERRQ(ierr);
  ierr = DSSetDimensions(ds,bv->k,0,bv->l,bv->k);CHKERRQ(ierr);
  /* [V,S] = eig(B'*B) */
  ierr = DSGetMat(ds,DS_MAT_A,&A);CHKERRQ(ierr);
  ierr = BVDot(bv,bv,A);CHKERRQ(ierr);
  ierr = DSRestoreMat(ds,DS_MAT_A,&A);CHKERRQ(ierr);
  ierr = PetscMalloc1(bv->m,&eig);CHKERRQ(ierr);
  ierr = DSSolve(ds,eig,NULL);CHKERRQ(ierr);
  ierr = DSSynchronize(ds,eig,NULL);CHKERRQ(ierr);
  ierr = DSVectors(ds,DS_MAT_X,NULL,NULL);CHKERRQ(ierr);
  /* M = diag(linspace(1/condn,1,n)./sqrt(diag(S)))' */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,bv->k,bv->k,NULL,&M);CHKERRQ(ierr);
  ierr = MatZeroEntries(M);CHKERRQ(ierr);
  ierr = MatDenseGetArray(M,&d);CHKERRQ(ierr);
  for (i=0;i<bv->k;i++) d[i+i*bv->m] = (1.0/condn+(1.0-1.0/condn)/(bv->k-1)*i)/PetscSqrtScalar(eig[i]);
  ierr = MatDenseRestoreArray(M,&d);CHKERRQ(ierr);
  /* G = X*M*X' */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,bv->k,bv->k,NULL,&Xt);CHKERRQ(ierr);
  ierr = DSGetMat(ds,DS_MAT_X,&X);CHKERRQ(ierr);
  ierr = MatTranspose(X,MAT_REUSE_MATRIX,&Xt);CHKERRQ(ierr);
  ierr = MatProductCreate(Xt,M,NULL,&G);CHKERRQ(ierr);
  ierr = MatProductSetType(G,MATPRODUCT_PtAP);CHKERRQ(ierr);
  ierr = MatProductSetFromOptions(G);CHKERRQ(ierr);
  ierr = MatProductSymbolic(G);CHKERRQ(ierr);
  ierr = MatProductNumeric(G);CHKERRQ(ierr);
  ierr = MatProductClear(G);CHKERRQ(ierr);
  ierr = MatDestroy(&X);CHKERRQ(ierr);
  ierr = MatDestroy(&Xt);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  /* B = B*G */
  ierr = BVMultInPlace(bv,G,bv->l,bv->k);CHKERRQ(ierr);
  ierr = MatDestroy(&G);CHKERRQ(ierr);
  ierr = PetscFree(eig);CHKERRQ(ierr);
  ierr = DSDestroy(&ds);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_SetRandom,bv,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMatMult - Computes the matrix-vector product for each column, Y=A*V.

   Neighbor-wise Collective on A

   Input Parameters:
+  V - basis vectors context
-  A - the matrix

   Output Parameter:
.  Y - the result

   Notes:
   Both V and Y must be distributed in the same manner. Only active columns
   (excluding the leading ones) are processed.
   In the result Y, columns are overwritten starting from the leading ones.
   The number of active columns in V and Y should match, although they need
   not be the same columns.

   It is possible to choose whether the computation is done column by column
   or as a Mat-Mat product, see BVSetMatMultMethod().

   Level: beginner

.seealso: BVCopy(), BVSetActiveColumns(), BVMatMultColumn(), BVMatMultTranspose(), BVMatMultHermitianTranspose(), BVSetMatMultMethod()
@*/
PetscErrorCode BVMatMult(BV V,Mat A,BV Y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  BVCheckOp(V,1,matmult);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidType(A,2);
  PetscValidHeaderSpecific(Y,BV_CLASSID,3);
  PetscValidType(Y,3);
  BVCheckSizes(Y,3);
  PetscCheckSameComm(V,1,A,2);
  PetscCheckSameTypeAndComm(V,1,Y,3);
  if (V->n!=Y->n) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension V %D, Y %D",V->n,Y->n);
  if (V->k-V->l!=Y->k-Y->l) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Y has %D active columns, should match %D active columns in V",Y->k-Y->l,V->k-V->l);

  ierr = PetscLogEventBegin(BV_MatMult,V,A,Y,0);CHKERRQ(ierr);
  ierr = (*V->ops->matmult)(V,A,Y);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_MatMult,V,A,Y,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)Y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMatMultTranspose - Computes the matrix-vector product with the transpose
   of a matrix for each column, Y=A^T*V.

   Neighbor-wise Collective on A

   Input Parameters:
+  V - basis vectors context
-  A - the matrix

   Output Parameter:
.  Y - the result

   Notes:
   Both V and Y must be distributed in the same manner. Only active columns
   (excluding the leading ones) are processed.
   In the result Y, columns are overwritten starting from the leading ones.
   The number of active columns in V and Y should match, although they need
   not be the same columns.

   Currently implemented via MatCreateTranspose().

   Level: beginner

.seealso: BVMatMult(), BVMatMultHermitianTranspose()
@*/
PetscErrorCode BVMatMultTranspose(BV V,Mat A,BV Y)
{
  PetscErrorCode ierr;
  Mat            AT;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidType(A,2);
  PetscValidHeaderSpecific(Y,BV_CLASSID,3);
  PetscValidType(Y,3);
  BVCheckSizes(Y,3);
  PetscCheckSameComm(V,1,A,2);
  PetscCheckSameTypeAndComm(V,1,Y,3);
  if (V->n!=Y->n) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension V %D, Y %D",V->n,Y->n);
  if (V->k-V->l>Y->m-Y->l) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Y has %D non-leading columns, not enough to store %D columns",Y->m-Y->l,V->k-V->l);

  ierr = MatCreateTranspose(A,&AT);CHKERRQ(ierr);
  ierr = BVMatMult(V,AT,Y);CHKERRQ(ierr);
  ierr = MatDestroy(&AT);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMatMultHermitianTranspose - Computes the matrix-vector product with the
   conjugate transpose of a matrix for each column, Y=A^H*V.

   Neighbor-wise Collective on A

   Input Parameters:
+  V - basis vectors context
-  A - the matrix

   Output Parameter:
.  Y - the result

   Note:
   Both V and Y must be distributed in the same manner. Only active columns
   (excluding the leading ones) are processed.
   In the result Y, columns are overwritten starting from the leading ones.
   The number of active columns in V and Y should match, although they need
   not be the same columns.

   Currently implemented via MatCreateHermitianTranspose().

   Level: beginner

.seealso: BVMatMult(), BVMatMultTranspose()
@*/
PetscErrorCode BVMatMultHermitianTranspose(BV V,Mat A,BV Y)
{
  PetscErrorCode ierr;
  Mat            AH;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidType(A,2);
  PetscValidHeaderSpecific(Y,BV_CLASSID,3);
  PetscValidType(Y,3);
  BVCheckSizes(Y,3);
  PetscCheckSameComm(V,1,A,2);
  PetscCheckSameTypeAndComm(V,1,Y,3);
  if (V->n!=Y->n) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension V %D, Y %D",V->n,Y->n);
  if (V->k-V->l>Y->m-Y->l) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"Y has %D non-leading columns, not enough to store %D columns",Y->m-Y->l,V->k-V->l);

  ierr = MatCreateHermitianTranspose(A,&AH);CHKERRQ(ierr);
  ierr = BVMatMult(V,AH,Y);CHKERRQ(ierr);
  ierr = MatDestroy(&AH);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMatMultColumn - Computes the matrix-vector product for a specified
   column, storing the result in the next column: v_{j+1}=A*v_j.

   Neighbor-wise Collective on A

   Input Parameters:
+  V - basis vectors context
.  A - the matrix
-  j - the column

   Level: beginner

.seealso: BVMatMult(), BVMatMultTransposeColumn(), BVMatMultHermitianTransposeColumn()
@*/
PetscErrorCode BVMatMultColumn(BV V,Mat A,PetscInt j)
{
  PetscErrorCode ierr;
  Vec            vj,vj1;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(V,1,A,2);
  PetscValidLogicalCollectiveInt(V,j,3);
  if (j<0) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  if (j+1>=V->m) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Result should go in index j+1=%D but BV only has %D columns",j+1,V->m);

  ierr = PetscLogEventBegin(BV_MatMultVec,V,A,0,0);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j,&vj);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j+1,&vj1);CHKERRQ(ierr);
  ierr = MatMult(A,vj,vj1);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j,&vj);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j+1,&vj1);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_MatMultVec,V,A,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMatMultTransposeColumn - Computes the transpose matrix-vector product for a
   specified column, storing the result in the next column: v_{j+1}=A^T*v_j.

   Neighbor-wise Collective on A

   Input Parameters:
+  V - basis vectors context
.  A - the matrix
-  j - the column

   Level: beginner

.seealso: BVMatMult(), BVMatMultColumn(), BVMatMultHermitianTransposeColumn()
@*/
PetscErrorCode BVMatMultTransposeColumn(BV V,Mat A,PetscInt j)
{
  PetscErrorCode ierr;
  Vec            vj,vj1;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(V,1,A,2);
  PetscValidLogicalCollectiveInt(V,j,3);
  if (j<0) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  if (j+1>=V->m) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Result should go in index j+1=%D but BV only has %D columns",j+1,V->m);

  ierr = PetscLogEventBegin(BV_MatMultVec,V,A,0,0);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j,&vj);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j+1,&vj1);CHKERRQ(ierr);
  ierr = MatMultTranspose(A,vj,vj1);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j,&vj);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j+1,&vj1);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_MatMultVec,V,A,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVMatMultHermitianTransposeColumn - Computes the matrix-vector product for a specified
   column, storing the result in the next column: v_{j+1}=A*v_j.

   Neighbor-wise Collective on A

   Input Parameters:
+  V - basis vectors context
.  A - the matrix
-  j - the column

   Level: beginner

.seealso: BVMatMult(), BVMatMultColumn(), BVMatMultTransposeColumn()
@*/
PetscErrorCode BVMatMultHermitianTransposeColumn(BV V,Mat A,PetscInt j)
{
  PetscErrorCode ierr;
  Vec            vj,vj1;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscCheckSameComm(V,1,A,2);
  PetscValidLogicalCollectiveInt(V,j,3);
  if (j<0) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Index j must be non-negative");
  if (j+1>=V->m) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Result should go in index j+1=%D but BV only has %D columns",j+1,V->m);

  ierr = PetscLogEventBegin(BV_MatMultVec,V,A,0,0);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j,&vj);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j+1,&vj1);CHKERRQ(ierr);
  ierr = MatMultHermitianTranspose(A,vj,vj1);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j,&vj);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j+1,&vj1);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_MatMultVec,V,A,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

