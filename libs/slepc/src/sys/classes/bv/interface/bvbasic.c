/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Basic BV routines
*/

#include <slepc/private/bvimpl.h>      /*I "slepcbv.h" I*/

PetscBool         BVRegisterAllCalled = PETSC_FALSE;
PetscFunctionList BVList = 0;

/*@C
   BVSetType - Selects the type for the BV object.

   Logically Collective on bv

   Input Parameter:
+  bv   - the basis vectors context
-  type - a known type

   Options Database Key:
.  -bv_type <type> - Sets BV type

   Level: intermediate

.seealso: BVGetType()
@*/
PetscErrorCode BVSetType(BV bv,BVType type)
{
  PetscErrorCode ierr,(*r)(BV);
  PetscBool      match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidCharPointer(type,2);

  ierr = PetscObjectTypeCompare((PetscObject)bv,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);
  ierr = PetscStrcmp(type,BVTENSOR,&match);CHKERRQ(ierr);
  if (match) SETERRQ(PetscObjectComm((PetscObject)bv),1,"Use BVCreateTensor() to create a BV of type tensor");

  ierr =  PetscFunctionListFind(BVList,type,&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested BV type %s",type);

  if (bv->ops->destroy) { ierr = (*bv->ops->destroy)(bv);CHKERRQ(ierr); }
  ierr = PetscMemzero(bv->ops,sizeof(struct _BVOps));CHKERRQ(ierr);

  ierr = PetscObjectChangeTypeName((PetscObject)bv,type);CHKERRQ(ierr);
  if (bv->n < 0 && bv->N < 0) {
    bv->ops->create = r;
  } else {
    ierr = PetscLogEventBegin(BV_Create,bv,0,0,0);CHKERRQ(ierr);
    ierr = (*r)(bv);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(BV_Create,bv,0,0,0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@C
   BVGetType - Gets the BV type name (as a string) from the BV context.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  name - name of the type of basis vectors

   Level: intermediate

.seealso: BVSetType()
@*/
PetscErrorCode BVGetType(BV bv,BVType *type)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidPointer(type,2);
  *type = ((PetscObject)bv)->type_name;
  PetscFunctionReturn(0);
}

/*@
   BVSetSizes - Sets the local and global sizes, and the number of columns.

   Collective on bv

   Input Parameters:
+  bv - the basis vectors
.  n  - the local size (or PETSC_DECIDE to have it set)
.  N  - the global size (or PETSC_DECIDE)
-  m  - the number of columns

   Notes:
   n and N cannot be both PETSC_DECIDE.
   If one processor calls this with N of PETSC_DECIDE then all processors must,
   otherwise the program will hang.

   Level: beginner

.seealso: BVSetSizesFromVec(), BVGetSizes(), BVResize()
@*/
PetscErrorCode BVSetSizes(BV bv,PetscInt n,PetscInt N,PetscInt m)
{
  PetscErrorCode ierr;
  PetscInt       ma;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (N >= 0) PetscValidLogicalCollectiveInt(bv,N,3);
  PetscValidLogicalCollectiveInt(bv,m,4);
  if (N >= 0 && n > N) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Local size %D cannot be larger than global size %D",n,N);
  if (m <= 0) SETERRQ1(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Number of columns %D must be positive",m);
  if ((bv->n >= 0 || bv->N >= 0) && (bv->n != n || bv->N != N)) SETERRQ4(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Cannot change/reset vector sizes to %D local %D global after previously setting them to %D local %D global",n,N,bv->n,bv->N);
  if (bv->m > 0 && bv->m != m) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Cannot change the number of columns to %D after previously setting it to %D; use BVResize()",m,bv->m);
  bv->n = n;
  bv->N = N;
  bv->m = m;
  bv->k = m;
  if (!bv->t) {  /* create template vector and get actual dimensions */
    ierr = VecCreate(PetscObjectComm((PetscObject)bv),&bv->t);CHKERRQ(ierr);
    ierr = VecSetSizes(bv->t,bv->n,bv->N);CHKERRQ(ierr);
    ierr = VecSetFromOptions(bv->t);CHKERRQ(ierr);
    ierr = VecGetSize(bv->t,&bv->N);CHKERRQ(ierr);
    ierr = VecGetLocalSize(bv->t,&bv->n);CHKERRQ(ierr);
    if (bv->matrix) {  /* check compatible dimensions of user-provided matrix */
      ierr = MatGetLocalSize(bv->matrix,&ma,NULL);CHKERRQ(ierr);
      if (bv->n!=ma) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Local dimension %D does not match that of matrix given at BVSetMatrix %D",bv->n,ma);
    }
  }
  if (bv->ops->create) {
    ierr = PetscLogEventBegin(BV_Create,bv,0,0,0);CHKERRQ(ierr);
    ierr = (*bv->ops->create)(bv);CHKERRQ(ierr);
    ierr = PetscLogEventEnd(BV_Create,bv,0,0,0);CHKERRQ(ierr);
    bv->ops->create = 0;
    bv->defersfo = PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}

/*@
   BVSetSizesFromVec - Sets the local and global sizes, and the number of columns.
   Local and global sizes are specified indirectly by passing a template vector.

   Collective on bv

   Input Parameters:
+  bv - the basis vectors
.  t  - the template vector
-  m  - the number of columns

   Level: beginner

.seealso: BVSetSizes(), BVGetSizes(), BVResize()
@*/
PetscErrorCode BVSetSizesFromVec(BV bv,Vec t,PetscInt m)
{
  PetscErrorCode ierr;
  PetscInt       ma;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(t,VEC_CLASSID,2);
  PetscCheckSameComm(bv,1,t,2);
  PetscValidLogicalCollectiveInt(bv,m,3);
  if (m <= 0) SETERRQ1(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Number of columns %D must be positive",m);
  if (bv->t) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Template vector was already set by a previous call to BVSetSizes/FromVec");
  ierr = VecGetSize(t,&bv->N);CHKERRQ(ierr);
  ierr = VecGetLocalSize(t,&bv->n);CHKERRQ(ierr);
  if (bv->matrix) {  /* check compatible dimensions of user-provided matrix */
    ierr = MatGetLocalSize(bv->matrix,&ma,NULL);CHKERRQ(ierr);
    if (bv->n!=ma) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Local dimension %D does not match that of matrix given at BVSetMatrix %D",bv->n,ma);
  }
  bv->m = m;
  bv->k = m;
  bv->t = t;
  ierr = PetscObjectReference((PetscObject)t);CHKERRQ(ierr);
  if (bv->ops->create) {
    ierr = (*bv->ops->create)(bv);CHKERRQ(ierr);
    bv->ops->create = 0;
    bv->defersfo = PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}

/*@
   BVGetSizes - Returns the local and global sizes, and the number of columns.

   Not Collective

   Input Parameter:
.  bv - the basis vectors

   Output Parameters:
+  n  - the local size
.  N  - the global size
-  m  - the number of columns

   Note:
   Normal usage requires that bv has already been given its sizes, otherwise
   the call fails. However, this function can also be used to determine if
   a BV object has been initialized completely (sizes and type). For this,
   call with n=NULL and N=NULL, then a return value of m=0 indicates that
   the BV object is not ready for use yet.

   Level: beginner

.seealso: BVSetSizes(), BVSetSizesFromVec()
@*/
PetscErrorCode BVGetSizes(BV bv,PetscInt *n,PetscInt *N,PetscInt *m)
{
  PetscFunctionBegin;
  if (!bv) {
    if (m && !n && !N) *m = 0;
    PetscFunctionReturn(0);
  }
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (n || N) BVCheckSizes(bv,1);
  if (n) *n = bv->n;
  if (N) *N = bv->N;
  if (m) *m = bv->m;
  if (m && !n && !N && !((PetscObject)bv)->type_name) *m = 0;
  PetscFunctionReturn(0);
}

/*@
   BVSetNumConstraints - Set the number of constraints.

   Logically Collective on V

   Input Parameters:
+  V  - basis vectors
-  nc - number of constraints

   Notes:
   This function sets the number of constraints to nc and marks all remaining
   columns as regular. Normal user would call BVInsertConstraints() instead.

   If nc is smaller than the previously set value, then some of the constraints
   are discarded. In particular, using nc=0 removes all constraints preserving
   the content of regular columns.

   Level: developer

.seealso: BVInsertConstraints()
@*/
PetscErrorCode BVSetNumConstraints(BV V,PetscInt nc)
{
  PetscErrorCode ierr;
  PetscInt       total,diff,i;
  Vec            x,y;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(V,nc,2);
  if (nc<0) SETERRQ1(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Number of constraints (given %D) cannot be negative",nc);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  if (V->ci[0]!=-V->nc-1 || V->ci[1]!=-V->nc-1) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_SUP,"Cannot call BVSetNumConstraints after BVGetColumn");

  diff = nc-V->nc;
  if (!diff) PetscFunctionReturn(0);
  total = V->nc+V->m;
  if (total-nc<=0) SETERRQ(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_OUTOFRANGE,"Not enough columns for the given nc value");
  if (diff<0) {  /* lessen constraints, shift contents of BV */
    for (i=0;i<V->m;i++) {
      ierr = BVGetColumn(V,i,&x);CHKERRQ(ierr);
      ierr = BVGetColumn(V,i+diff,&y);CHKERRQ(ierr);
      ierr = VecCopy(x,y);CHKERRQ(ierr);
      ierr = BVRestoreColumn(V,i,&x);CHKERRQ(ierr);
      ierr = BVRestoreColumn(V,i+diff,&y);CHKERRQ(ierr);
    }
  }
  V->nc = nc;
  V->ci[0] = -V->nc-1;
  V->ci[1] = -V->nc-1;
  V->l = 0;
  V->m = total-nc;
  V->k = V->m;
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVGetNumConstraints - Returns the number of constraints.

   Not Collective

   Input Parameter:
.  bv - the basis vectors

   Output Parameters:
.  nc - the number of constraints

   Level: advanced

.seealso: BVGetSizes(), BVInsertConstraints()
@*/
PetscErrorCode BVGetNumConstraints(BV bv,PetscInt *nc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidIntPointer(nc,2);
  *nc = bv->nc;
  PetscFunctionReturn(0);
}

/*@
   BVResize - Change the number of columns.

   Collective on bv

   Input Parameters:
+  bv   - the basis vectors
.  m    - the new number of columns
-  copy - a flag indicating whether current values should be kept

   Note:
   Internal storage is reallocated. If the copy flag is set to true, then
   the contents are copied to the leading part of the new space.

   Level: advanced

.seealso: BVSetSizes(), BVSetSizesFromVec()
@*/
PetscErrorCode BVResize(BV bv,PetscInt m,PetscBool copy)
{
  PetscErrorCode    ierr;
  PetscScalar       *array;
  const PetscScalar *omega;
  Vec               v;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,m,2);
  PetscValidLogicalCollectiveBool(bv,copy,3);
  PetscValidType(bv,1);
  if (m <= 0) SETERRQ1(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Number of columns %D must be positive",m);
  if (bv->nc && !bv->issplit) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Cannot resize a BV with constraints");
  if (bv->m == m) PetscFunctionReturn(0);
  BVCheckOp(bv,1,resize);

  ierr = PetscLogEventBegin(BV_Create,bv,0,0,0);CHKERRQ(ierr);
  ierr = (*bv->ops->resize)(bv,m,copy);CHKERRQ(ierr);
  ierr = VecDestroy(&bv->buffer);CHKERRQ(ierr);
  ierr = BVDestroy(&bv->cached);CHKERRQ(ierr);
  ierr = PetscFree2(bv->h,bv->c);CHKERRQ(ierr);
  if (bv->omega) {
    if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
      ierr = VecCreateSeqCUDA(PETSC_COMM_SELF,m,&v);CHKERRQ(ierr);
#else
      SETERRQ(PetscObjectComm((PetscObject)bv),1,"Something wrong happened");
#endif
    } else {
      ierr = VecCreateSeq(PETSC_COMM_SELF,m,&v);CHKERRQ(ierr);
    }
    ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)v);CHKERRQ(ierr);
    if (copy) {
      ierr = VecGetArray(v,&array);CHKERRQ(ierr);
      ierr = VecGetArrayRead(bv->omega,&omega);CHKERRQ(ierr);
      ierr = PetscArraycpy(array,omega,PetscMin(m,bv->m));CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(bv->omega,&omega);CHKERRQ(ierr);
      ierr = VecRestoreArray(v,&array);CHKERRQ(ierr);
    } else {
      ierr = VecSet(v,1.0);CHKERRQ(ierr);
    }
    ierr = VecDestroy(&bv->omega);CHKERRQ(ierr);
    bv->omega = v;
  }
  bv->m = m;
  bv->k = PetscMin(bv->k,m);
  bv->l = PetscMin(bv->l,m);
  ierr = PetscLogEventEnd(BV_Create,bv,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVSetActiveColumns - Specify the columns that will be involved in operations.

   Logically Collective on bv

   Input Parameters:
+  bv - the basis vectors context
.  l  - number of leading columns
-  k  - number of active columns

   Notes:
   In operations such as BVMult() or BVDot(), only the first k columns are
   considered. This is useful when the BV is filled from left to right, so
   the last m-k columns do not have relevant information.

   Also in operations such as BVMult() or BVDot(), the first l columns are
   normally not included in the computation. See the manpage of each
   operation.

   In orthogonalization operations, the first l columns are treated
   differently: they participate in the orthogonalization but the computed
   coefficients are not stored.

   Level: intermediate

.seealso: BVGetActiveColumns(), BVSetSizes()
@*/
PetscErrorCode BVSetActiveColumns(BV bv,PetscInt l,PetscInt k)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveInt(bv,l,2);
  PetscValidLogicalCollectiveInt(bv,k,3);
  BVCheckSizes(bv,1);
  if (k==PETSC_DECIDE || k==PETSC_DEFAULT) {
    bv->k = bv->m;
  } else {
    if (k<0 || k>bv->m) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of k. Must be between 0 and m");
    bv->k = k;
  }
  if (l==PETSC_DECIDE || l==PETSC_DEFAULT) {
    bv->l = 0;
  } else {
    if (l<0 || l>bv->k) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Illegal value of l. Must be between 0 and k");
    bv->l = l;
  }
  PetscFunctionReturn(0);
}

/*@
   BVGetActiveColumns - Returns the current active dimensions.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
+  l  - number of leading columns
-  k  - number of active columns

   Level: intermediate

.seealso: BVSetActiveColumns()
@*/
PetscErrorCode BVGetActiveColumns(BV bv,PetscInt *l,PetscInt *k)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (l) *l = bv->l;
  if (k) *k = bv->k;
  PetscFunctionReturn(0);
}

/*@
   BVSetMatrix - Specifies the inner product to be used in orthogonalization.

   Collective on bv

   Input Parameters:
+  bv    - the basis vectors context
.  B     - a symmetric matrix (may be NULL)
-  indef - a flag indicating if the matrix is indefinite

   Notes:
   This is used to specify a non-standard inner product, whose matrix
   representation is given by B. Then, all inner products required during
   orthogonalization are computed as (x,y)_B=y^H*B*x rather than the
   standard form (x,y)=y^H*x.

   Matrix B must be real symmetric (or complex Hermitian). A genuine inner
   product requires that B is also positive (semi-)definite. However, we
   also allow for an indefinite B (setting indef=PETSC_TRUE), in which
   case the orthogonalization uses an indefinite inner product.

   This affects operations BVDot(), BVNorm(), BVOrthogonalize(), and variants.

   Setting B=NULL has the same effect as if the identity matrix was passed.

   Level: advanced

.seealso: BVGetMatrix(), BVDot(), BVNorm(), BVOrthogonalize()
@*/
PetscErrorCode BVSetMatrix(BV bv,Mat B,PetscBool indef)
{
  PetscErrorCode ierr;
  PetscInt       m,n;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveBool(bv,indef,3);
  if (B) {
    PetscValidHeaderSpecific(B,MAT_CLASSID,2);
    ierr = MatGetLocalSize(B,&m,&n);CHKERRQ(ierr);
    if (m!=n) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_SIZ,"Matrix must be square");
    if (bv->m && bv->n!=n) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension BV %D, Mat %D",bv->n,n);
  }
  if (B) { ierr = PetscObjectReference((PetscObject)B);CHKERRQ(ierr); }
  ierr = MatDestroy(&bv->matrix);CHKERRQ(ierr);
  bv->matrix = B;
  bv->indef  = indef;
  ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVGetMatrix - Retrieves the matrix representation of the inner product.

   Not collective, though a parallel Mat may be returned

   Input Parameter:
.  bv    - the basis vectors context

   Output Parameter:
+  B     - the matrix of the inner product (may be NULL)
-  indef - the flag indicating if the matrix is indefinite

   Level: advanced

.seealso: BVSetMatrix()
@*/
PetscErrorCode BVGetMatrix(BV bv,Mat *B,PetscBool *indef)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (B)     *B     = bv->matrix;
  if (indef) *indef = bv->indef;
  PetscFunctionReturn(0);
}

/*@
   BVApplyMatrix - Multiplies a vector by the matrix representation of the
   inner product.

   Neighbor-wise Collective on bv

   Input Parameter:
+  bv - the basis vectors context
-  x  - the vector

   Output Parameter:
.  y  - the result

   Note:
   If no matrix was specified this function copies the vector.

   Level: advanced

.seealso: BVSetMatrix(), BVApplyMatrixBV()
@*/
PetscErrorCode BVApplyMatrix(BV bv,Vec x,Vec y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(x,VEC_CLASSID,2);
  PetscValidHeaderSpecific(y,VEC_CLASSID,3);
  if (bv->matrix) {
    ierr = BV_IPMatMult(bv,x);CHKERRQ(ierr);
    ierr = VecCopy(bv->Bx,y);CHKERRQ(ierr);
  } else {
    ierr = VecCopy(x,y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   BVApplyMatrixBV - Multiplies the BV vectors by the matrix representation
   of the inner product.

   Neighbor-wise Collective on X

   Input Parameter:
.  X - the basis vectors context

   Output Parameter:
.  Y - the basis vectors to store the result (optional)

   Note:
   This function computes Y = B*X, where B is the matrix given with
   BVSetMatrix(). This operation is computed as in BVMatMult().
   If no matrix was specified, then it just copies Y = X.

   If no Y is given, the result is stored internally in the cached BV.

   Level: developer

.seealso: BVSetMatrix(), BVApplyMatrix(), BVMatMult(), BVGetCachedBV()
@*/
PetscErrorCode BVApplyMatrixBV(BV X,BV Y)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(X,BV_CLASSID,1);
  if (Y) {
    PetscValidHeaderSpecific(Y,BV_CLASSID,2);
    if (X->matrix) {
      ierr = BVMatMult(X,X->matrix,Y);CHKERRQ(ierr);
    } else {
      ierr = BVCopy(X,Y);CHKERRQ(ierr);
    }
  } else {
    ierr = BV_IPMatMultBV(X);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   BVSetSignature - Sets the signature matrix to be used in orthogonalization.

   Logically Collective on bv

   Input Parameter:
+  bv    - the basis vectors context
-  omega - a vector representing the diagonal of the signature matrix

   Note:
   The signature matrix Omega = V'*B*V is relevant only for an indefinite B.

   Level: developer

.seealso: BVSetMatrix(), BVGetSignature()
@*/
PetscErrorCode BVSetSignature(BV bv,Vec omega)
{
  PetscErrorCode    ierr;
  PetscInt          i,n;
  const PetscScalar *pomega;
  PetscScalar       *intern;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscValidHeaderSpecific(omega,VEC_CLASSID,2);
  PetscValidType(omega,2);

  ierr = VecGetSize(omega,&n);CHKERRQ(ierr);
  if (n!=bv->k) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_SIZ,"Vec argument has %D elements, should be %D",n,bv->k);
  ierr = BV_AllocateSignature(bv);CHKERRQ(ierr);
  if (bv->indef) {
    ierr = VecGetArrayRead(omega,&pomega);CHKERRQ(ierr);
    ierr = VecGetArray(bv->omega,&intern);CHKERRQ(ierr);
    for (i=0;i<n;i++) intern[bv->nc+i] = pomega[i];
    ierr = VecRestoreArray(bv->omega,&intern);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(omega,&pomega);CHKERRQ(ierr);
  } else {
    ierr = PetscInfo(bv,"Ignoring signature because BV is not indefinite\n");CHKERRQ(ierr);
  }
  ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVGetSignature - Retrieves the signature matrix from last orthogonalization.

   Not collective

   Input Parameter:
.  bv    - the basis vectors context

   Output Parameter:
.  omega - a vector representing the diagonal of the signature matrix

   Note:
   The signature matrix Omega = V'*B*V is relevant only for an indefinite B.

   Level: developer

.seealso: BVSetMatrix(), BVSetSignature()
@*/
PetscErrorCode BVGetSignature(BV bv,Vec omega)
{
  PetscErrorCode    ierr;
  PetscInt          i,n;
  PetscScalar       *pomega;
  const PetscScalar *intern;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscValidHeaderSpecific(omega,VEC_CLASSID,2);
  PetscValidType(omega,2);

  ierr = VecGetSize(omega,&n);CHKERRQ(ierr);
  if (n!=bv->k) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_SIZ,"Vec argument has %D elements, should be %D",n,bv->k);
  if (bv->indef && bv->omega) {
    ierr = VecGetArray(omega,&pomega);CHKERRQ(ierr);
    ierr = VecGetArrayRead(bv->omega,&intern);CHKERRQ(ierr);
    for (i=0;i<n;i++) pomega[i] = intern[bv->nc+i];
    ierr = VecRestoreArrayRead(bv->omega,&intern);CHKERRQ(ierr);
    ierr = VecRestoreArray(omega,&pomega);CHKERRQ(ierr);
  } else {
    ierr = VecSet(omega,1.0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   BVSetBufferVec - Attach a vector object to be used as buffer space for
   several operations.

   Collective on bv

   Input Parameters:
+  bv     - the basis vectors context)
-  buffer - the vector

   Notes:
   Use BVGetBufferVec() to retrieve the vector (for example, to free it
   at the end of the computations).

   The vector must be sequential of length (nc+m)*m, where m is the number
   of columns of bv and nc is the number of constraints.

   Level: developer

.seealso: BVGetBufferVec(), BVSetSizes(), BVGetNumConstraints()
@*/
PetscErrorCode BVSetBufferVec(BV bv,Vec buffer)
{
  PetscErrorCode ierr;
  PetscInt       ld,n;
  PetscMPIInt    size;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(buffer,VEC_CLASSID,2);
  BVCheckSizes(bv,1);
  ierr = VecGetSize(buffer,&n);CHKERRQ(ierr);
  ld = bv->m+bv->nc;
  if (n != ld*bv->m) SETERRQ1(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_SIZ,"Buffer size must be %d",ld*bv->m);
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)buffer),&size);CHKERRQ(ierr);
  if (size>1) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Buffer must be a sequential vector");

  ierr = PetscObjectReference((PetscObject)buffer);CHKERRQ(ierr);
  ierr = VecDestroy(&bv->buffer);CHKERRQ(ierr);
  bv->buffer = buffer;
  ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)bv->buffer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVGetBufferVec - Obtain the buffer vector associated with the BV object.

   Not Collective, but Vec returned is parallel if BV is parallel

   Input Parameters:
.  bv - the basis vectors context

   Output Parameter:
.  buffer - vector

   Notes:
   The vector is created if not available previously. It is a sequential vector
   of length (nc+m)*m, where m is the number of columns of bv and nc is the number
   of constraints.

   Developer Notes:
   The buffer vector is viewed as a column-major matrix with leading dimension
   ld=nc+m, and m columns at most. In the most common usage, it has the structure
.vb
      | | C |
      |s|---|
      | | H |
.ve
   where H is an upper Hessenberg matrix of order m x (m-1), C contains coefficients
   related to orthogonalization against constraints (first nc rows), and s is the
   first column that contains scratch values computed during Gram-Schmidt
   orthogonalization. In particular, BVDotColumn() and BVMultColumn() use s to
   store the coefficients.

   Level: developer

.seealso: BVSetBufferVec(), BVSetSizes(), BVGetNumConstraints(), BVDotColumn(), BVMultColumn()
@*/
PetscErrorCode BVGetBufferVec(BV bv,Vec *buffer)
{
  PetscErrorCode ierr;
  PetscInt       ld;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidPointer(buffer,2);
  BVCheckSizes(bv,1);
  if (!bv->buffer) {
    ld = bv->m+bv->nc;
    ierr = VecCreate(PETSC_COMM_SELF,&bv->buffer);CHKERRQ(ierr);
    ierr = VecSetSizes(bv->buffer,PETSC_DECIDE,ld*bv->m);CHKERRQ(ierr);
    ierr = VecSetType(bv->buffer,((PetscObject)bv->t)->type_name);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)bv->buffer);CHKERRQ(ierr);
  }
  *buffer = bv->buffer;
  PetscFunctionReturn(0);
}

/*@
   BVSetRandomContext - Sets the PetscRandom object associated with the BV,
   to be used in operations that need random numbers.

   Collective on bv

   Input Parameters:
+  bv   - the basis vectors context
-  rand - the random number generator context

   Level: advanced

.seealso: BVGetRandomContext(), BVSetRandom(), BVSetRandomColumn()
@*/
PetscErrorCode BVSetRandomContext(BV bv,PetscRandom rand)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidHeaderSpecific(rand,PETSC_RANDOM_CLASSID,2);
  PetscCheckSameComm(bv,1,rand,2);
  ierr = PetscObjectReference((PetscObject)rand);CHKERRQ(ierr);
  ierr = PetscRandomDestroy(&bv->rand);CHKERRQ(ierr);
  bv->rand = rand;
  ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)bv->rand);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVGetRandomContext - Gets the PetscRandom object associated with the BV.

   Not Collective

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  rand - the random number generator context

   Level: advanced

.seealso: BVSetRandomContext(), BVSetRandom(), BVSetRandomColumn()
@*/
PetscErrorCode BVGetRandomContext(BV bv,PetscRandom* rand)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidPointer(rand,2);
  if (!bv->rand) {
    ierr = PetscRandomCreate(PetscObjectComm((PetscObject)bv),&bv->rand);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)bv->rand);CHKERRQ(ierr);
    if (bv->rrandom) {
      ierr = PetscRandomSetSeed(bv->rand,0x12345678);CHKERRQ(ierr);
      ierr = PetscRandomSeed(bv->rand);CHKERRQ(ierr);
    }
  }
  *rand = bv->rand;
  PetscFunctionReturn(0);
}

/*@
   BVSetFromOptions - Sets BV options from the options database.

   Collective on bv

   Input Parameter:
.  bv - the basis vectors context

   Level: beginner
@*/
PetscErrorCode BVSetFromOptions(BV bv)
{
  PetscErrorCode     ierr;
  char               type[256];
  PetscBool          flg1,flg2,flg3,flg4;
  PetscReal          r;
  BVOrthogType       otype;
  BVOrthogRefineType orefine;
  BVOrthogBlockType  oblock;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  ierr = BVRegisterAll();CHKERRQ(ierr);
  ierr = PetscObjectOptionsBegin((PetscObject)bv);CHKERRQ(ierr);
    ierr = PetscOptionsFList("-bv_type","Basis Vectors type","BVSetType",BVList,(char*)(((PetscObject)bv)->type_name?((PetscObject)bv)->type_name:BVSVEC),type,256,&flg1);CHKERRQ(ierr);
    if (flg1) {
      ierr = BVSetType(bv,type);CHKERRQ(ierr);
    } else if (!((PetscObject)bv)->type_name) {
      ierr = BVSetType(bv,BVSVEC);CHKERRQ(ierr);
    }

    otype = bv->orthog_type;
    ierr = PetscOptionsEnum("-bv_orthog_type","Orthogonalization method","BVSetOrthogonalization",BVOrthogTypes,(PetscEnum)otype,(PetscEnum*)&otype,&flg1);CHKERRQ(ierr);
    orefine = bv->orthog_ref;
    ierr = PetscOptionsEnum("-bv_orthog_refine","Iterative refinement mode during orthogonalization","BVSetOrthogonalization",BVOrthogRefineTypes,(PetscEnum)orefine,(PetscEnum*)&orefine,&flg2);CHKERRQ(ierr);
    r = bv->orthog_eta;
    ierr = PetscOptionsReal("-bv_orthog_eta","Parameter of iterative refinement during orthogonalization","BVSetOrthogonalization",r,&r,&flg3);CHKERRQ(ierr);
    oblock = bv->orthog_block;
    ierr = PetscOptionsEnum("-bv_orthog_block","Block orthogonalization method","BVSetOrthogonalization",BVOrthogBlockTypes,(PetscEnum)oblock,(PetscEnum*)&oblock,&flg4);CHKERRQ(ierr);
    if (flg1 || flg2 || flg3 || flg4) { ierr = BVSetOrthogonalization(bv,otype,orefine,r,oblock);CHKERRQ(ierr); }

    ierr = PetscOptionsEnum("-bv_matmult","Method for BVMatMult","BVSetMatMultMethod",BVMatMultTypes,(PetscEnum)bv->vmm,(PetscEnum*)&bv->vmm,NULL);CHKERRQ(ierr);

    /* undocumented option to generate random vectors that are independent of the number of processes */
    ierr = PetscOptionsGetBool(NULL,NULL,"-bv_reproducible_random",&bv->rrandom,NULL);CHKERRQ(ierr);

    if (bv->ops->create) bv->defersfo = PETSC_TRUE;   /* defer call to setfromoptions */
    else if (bv->ops->setfromoptions) {
      ierr = (*bv->ops->setfromoptions)(PetscOptionsObject,bv);CHKERRQ(ierr);
    }
    ierr = PetscObjectProcessOptionsHandlers(PetscOptionsObject,(PetscObject)bv);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  if (!bv->rand) { ierr = BVGetRandomContext(bv,&bv->rand);CHKERRQ(ierr); }
  ierr = PetscRandomSetFromOptions(bv->rand);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVSetOrthogonalization - Specifies the method used for the orthogonalization of
   vectors (classical or modified Gram-Schmidt with or without refinement), and
   for the block-orthogonalization (simultaneous orthogonalization of a set of
   vectors).

   Logically Collective on bv

   Input Parameters:
+  bv     - the basis vectors context
.  type   - the method of vector orthogonalization
.  refine - type of refinement
.  eta    - parameter for selective refinement
-  block  - the method of block orthogonalization

   Options Database Keys:
+  -bv_orthog_type <type> - Where <type> is cgs for Classical Gram-Schmidt orthogonalization
                         (default) or mgs for Modified Gram-Schmidt orthogonalization
.  -bv_orthog_refine <ref> - Where <ref> is one of never, ifneeded (default) or always
.  -bv_orthog_eta <eta> -  For setting the value of eta
-  -bv_orthog_block <block> - Where <block> is the block-orthogonalization method

   Notes:
   The default settings work well for most problems.

   The parameter eta should be a real value between 0 and 1 (or PETSC_DEFAULT).
   The value of eta is used only when the refinement type is "ifneeded".

   When using several processors, MGS is likely to result in bad scalability.

   If the method set for block orthogonalization is GS, then the computation
   is done column by column with the vector orthogonalization.

   Level: advanced

.seealso: BVOrthogonalizeColumn(), BVGetOrthogonalization(), BVOrthogType, BVOrthogRefineType, BVOrthogBlockType
@*/
PetscErrorCode BVSetOrthogonalization(BV bv,BVOrthogType type,BVOrthogRefineType refine,PetscReal eta,BVOrthogBlockType block)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveEnum(bv,type,2);
  PetscValidLogicalCollectiveEnum(bv,refine,3);
  PetscValidLogicalCollectiveReal(bv,eta,4);
  PetscValidLogicalCollectiveEnum(bv,block,5);
  switch (type) {
    case BV_ORTHOG_CGS:
    case BV_ORTHOG_MGS:
      bv->orthog_type = type;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Unknown orthogonalization type");
  }
  switch (refine) {
    case BV_ORTHOG_REFINE_NEVER:
    case BV_ORTHOG_REFINE_IFNEEDED:
    case BV_ORTHOG_REFINE_ALWAYS:
      bv->orthog_ref = refine;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Unknown refinement type");
  }
  if (eta == PETSC_DEFAULT) {
    bv->orthog_eta = 0.7071;
  } else {
    if (eta <= 0.0 || eta > 1.0) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"Invalid eta value");
    bv->orthog_eta = eta;
  }
  switch (block) {
    case BV_ORTHOG_BLOCK_GS:
    case BV_ORTHOG_BLOCK_CHOL:
    case BV_ORTHOG_BLOCK_TSQR:
    case BV_ORTHOG_BLOCK_TSQRCHOL:
    case BV_ORTHOG_BLOCK_SVQB:
      bv->orthog_block = block;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Unknown block orthogonalization type");
  }
  PetscFunctionReturn(0);
}

/*@
   BVGetOrthogonalization - Gets the orthogonalization settings from the BV object.

   Not Collective

   Input Parameter:
.  bv - basis vectors context

   Output Parameter:
+  type   - the method of vector orthogonalization
.  refine - type of refinement
.  eta    - parameter for selective refinement
-  block  - the method of block orthogonalization

   Level: advanced

.seealso: BVOrthogonalizeColumn(), BVSetOrthogonalization(), BVOrthogType, BVOrthogRefineType, BVOrthogBlockType
@*/
PetscErrorCode BVGetOrthogonalization(BV bv,BVOrthogType *type,BVOrthogRefineType *refine,PetscReal *eta,BVOrthogBlockType *block)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  if (type)   *type   = bv->orthog_type;
  if (refine) *refine = bv->orthog_ref;
  if (eta)    *eta    = bv->orthog_eta;
  if (block)  *block  = bv->orthog_block;
  PetscFunctionReturn(0);
}

/*@
   BVSetMatMultMethod - Specifies the method used for the BVMatMult() operation.

   Logically Collective on bv

   Input Parameters:
+  bv     - the basis vectors context
-  method - the method for the BVMatMult() operation

   Options Database Keys:
.  -bv_matmult <meth> - choose one of the methods: vecs, mat

   Notes:
   Allowed values are
+  BV_MATMULT_VECS - perform a matrix-vector multiply per each column
.  BV_MATMULT_MAT - carry out a Mat-Mat product with a dense matrix
-  BV_MATMULT_MAT_SAVE - this case is deprecated

   The default is BV_MATMULT_MAT except in the case of BVVECS.

   Level: advanced

.seealso: BVMatMult(), BVGetMatMultMethod(), BVMatMultType
@*/
PetscErrorCode BVSetMatMultMethod(BV bv,BVMatMultType method)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidLogicalCollectiveEnum(bv,method,2);
  switch (method) {
    case BV_MATMULT_VECS:
    case BV_MATMULT_MAT:
      bv->vmm = method;
      break;
    case BV_MATMULT_MAT_SAVE:
      ierr = PetscInfo(bv,"BV_MATMULT_MAT_SAVE is deprecated, using BV_MATMULT_MAT\n");CHKERRQ(ierr);
      bv->vmm = BV_MATMULT_MAT;
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Unknown matmult method");
  }
  PetscFunctionReturn(0);
}

/*@
   BVGetMatMultMethod - Gets the method used for the BVMatMult() operation.

   Not Collective

   Input Parameter:
.  bv - basis vectors context

   Output Parameter:
.  method - the method for the BVMatMult() operation

   Level: advanced

.seealso: BVMatMult(), BVSetMatMultMethod(), BVMatMultType
@*/
PetscErrorCode BVGetMatMultMethod(BV bv,BVMatMultType *method)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidPointer(bv,method);
  *method = bv->vmm;
  PetscFunctionReturn(0);
}

/*@
   BVGetColumn - Returns a Vec object that contains the entries of the
   requested column of the basis vectors object.

   Logically Collective on bv

   Input Parameters:
+  bv - the basis vectors context
-  j  - the index of the requested column

   Output Parameter:
.  v  - vector containing the jth column

   Notes:
   The returned Vec must be seen as a reference (not a copy) of the BV
   column, that is, modifying the Vec will change the BV entries as well.

   The returned Vec must not be destroyed. BVRestoreColumn() must be
   called when it is no longer needed. At most, two columns can be fetched,
   that is, this function can only be called twice before the corresponding
   BVRestoreColumn() is invoked.

   A negative index j selects the i-th constraint, where i=-j. Constraints
   should not be modified.

   Level: beginner

.seealso: BVRestoreColumn(), BVInsertConstraints()
@*/
PetscErrorCode BVGetColumn(BV bv,PetscInt j,Vec *v)
{
  PetscErrorCode ierr;
  PetscInt       l;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  BVCheckOp(bv,1,getcolumn);
  PetscValidLogicalCollectiveInt(bv,j,2);
  if (j<0 && -j>bv->nc) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"You requested constraint %D but only %D are available",-j,bv->nc);
  if (j>=bv->m) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"You requested column %D but only %D are available",j,bv->m);
  if (j==bv->ci[0] || j==bv->ci[1]) SETERRQ1(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Column %D already fetched in a previous call to BVGetColumn",j);
  l = BVAvailableVec;
  if (l==-1) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_SUP,"Too many requested columns; you must call BVRestoreColumn for one of the previously fetched columns");
  ierr = (*bv->ops->getcolumn)(bv,j,v);CHKERRQ(ierr);
  bv->ci[l] = j;
  ierr = PetscObjectStateGet((PetscObject)bv->cv[l],&bv->st[l]);CHKERRQ(ierr);
  ierr = PetscObjectGetId((PetscObject)bv->cv[l],&bv->id[l]);CHKERRQ(ierr);
  *v = bv->cv[l];
  PetscFunctionReturn(0);
}

/*@
   BVRestoreColumn - Restore a column obtained with BVGetColumn().

   Logically Collective on bv

   Input Parameters:
+  bv - the basis vectors context
.  j  - the index of the column
-  v  - vector obtained with BVGetColumn()

   Note:
   The arguments must match the corresponding call to BVGetColumn().

   Level: beginner

.seealso: BVGetColumn()
@*/
PetscErrorCode BVRestoreColumn(BV bv,PetscInt j,Vec *v)
{
  PetscErrorCode   ierr;
  PetscObjectId    id;
  PetscObjectState st;
  PetscInt         l;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  PetscValidLogicalCollectiveInt(bv,j,2);
  PetscValidPointer(v,3);
  PetscValidHeaderSpecific(*v,VEC_CLASSID,3);
  if (j<0 && -j>bv->nc) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"You requested constraint %D but only %D are available",-j,bv->nc);
  if (j>=bv->m) SETERRQ2(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_OUTOFRANGE,"You requested column %D but only %D are available",j,bv->m);
  if (j!=bv->ci[0] && j!=bv->ci[1]) SETERRQ1(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Column %D has not been fetched with a call to BVGetColumn",j);
  l = (j==bv->ci[0])? 0: 1;
  ierr = PetscObjectGetId((PetscObject)*v,&id);CHKERRQ(ierr);
  if (id!=bv->id[l]) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Argument 3 is not the same Vec that was obtained with BVGetColumn");
  ierr = PetscObjectStateGet((PetscObject)*v,&st);CHKERRQ(ierr);
  if (st!=bv->st[l]) {
    ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  }
  if (bv->ops->restorecolumn) {
    ierr = (*bv->ops->restorecolumn)(bv,j,v);CHKERRQ(ierr);
  } else bv->cv[l] = NULL;
  bv->ci[l] = -bv->nc-1;
  bv->st[l] = -1;
  bv->id[l] = 0;
  *v = NULL;
  PetscFunctionReturn(0);
}

/*@C
   BVGetArray - Returns a pointer to a contiguous array that contains this
   processor's portion of the BV data.

   Logically Collective on bv

   Input Parameters:
.  bv - the basis vectors context

   Output Parameter:
.  a  - location to put pointer to the array

   Notes:
   BVRestoreArray() must be called when access to the array is no longer needed.
   This operation may imply a data copy, for BV types that do not store
   data contiguously in memory.

   The pointer will normally point to the first entry of the first column,
   but if the BV has constraints then these go before the regular columns.

   Level: advanced

.seealso: BVRestoreArray(), BVInsertConstraints()
@*/
PetscErrorCode BVGetArray(BV bv,PetscScalar **a)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  BVCheckOp(bv,1,getarray);
  ierr = (*bv->ops->getarray)(bv,a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   BVRestoreArray - Restore the BV object after BVGetArray() has been called.

   Logically Collective on bv

   Input Parameters:
+  bv - the basis vectors context
-  a  - location of pointer to array obtained from BVGetArray()

   Note:
   This operation may imply a data copy, for BV types that do not store
   data contiguously in memory.

   Level: advanced

.seealso: BVGetColumn()
@*/
PetscErrorCode BVRestoreArray(BV bv,PetscScalar **a)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  if (bv->ops->restorearray) {
    ierr = (*bv->ops->restorearray)(bv,a);CHKERRQ(ierr);
  }
  if (a) *a = NULL;
  ierr = PetscObjectStateIncrease((PetscObject)bv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   BVGetArrayRead - Returns a read-only pointer to a contiguous array that
   contains this processor's portion of the BV data.

   Not Collective

   Input Parameters:
.  bv - the basis vectors context

   Output Parameter:
.  a  - location to put pointer to the array

   Notes:
   BVRestoreArrayRead() must be called when access to the array is no
   longer needed. This operation may imply a data copy, for BV types that
   do not store data contiguously in memory.

   The pointer will normally point to the first entry of the first column,
   but if the BV has constraints then these go before the regular columns.

   Level: advanced

.seealso: BVRestoreArray(), BVInsertConstraints()
@*/
PetscErrorCode BVGetArrayRead(BV bv,const PetscScalar **a)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  BVCheckOp(bv,1,getarrayread);
  ierr = (*bv->ops->getarrayread)(bv,a);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   BVRestoreArrayRead - Restore the BV object after BVGetArrayRead() has
   been called.

   Not Collective

   Input Parameters:
+  bv - the basis vectors context
-  a  - location of pointer to array obtained from BVGetArrayRead()

   Level: advanced

.seealso: BVGetColumn()
@*/
PetscErrorCode BVRestoreArrayRead(BV bv,const PetscScalar **a)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  if (bv->ops->restorearrayread) {
    ierr = (*bv->ops->restorearrayread)(bv,a);CHKERRQ(ierr);
  }
  if (a) *a = NULL;
  PetscFunctionReturn(0);
}

/*@
   BVCreateVec - Creates a new Vec object with the same type and dimensions
   as the columns of the basis vectors object.

   Collective on bv

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  v  - the new vector

   Note:
   The user is responsible of destroying the returned vector.

   Level: beginner

.seealso: BVCreateMat()
@*/
PetscErrorCode BVCreateVec(BV bv,Vec *v)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscValidPointer(v,2);
  ierr = VecDuplicate(bv->t,v);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVCreateMat - Creates a new Mat object of dense type and copies the contents
   of the BV object.

   Collective on bv

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  A  - the new matrix

   Notes:
   The user is responsible of destroying the returned matrix.

   The matrix contains all columns of the BV, not just the active columns.

   Level: intermediate

.seealso: BVCreateFromMat(), BVCreateVec(), BVGetMat()
@*/
PetscErrorCode BVCreateMat(BV bv,Mat *A)
{
  PetscErrorCode    ierr;
  PetscScalar       *aa;
  const PetscScalar *vv;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscValidPointer(A,2);

  ierr = MatCreateDense(PetscObjectComm((PetscObject)bv->t),bv->n,PETSC_DECIDE,bv->N,bv->m,NULL,A);CHKERRQ(ierr);
  ierr = MatDenseGetArray(*A,&aa);CHKERRQ(ierr);
  ierr = BVGetArrayRead(bv,&vv);CHKERRQ(ierr);
  ierr = PetscArraycpy(aa,vv,bv->m*bv->n);CHKERRQ(ierr);
  ierr = BVRestoreArrayRead(bv,&vv);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(*A,&aa);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVGetMat - Returns a Mat object of dense type that shares the memory of
   the BV object.

   Collective on bv

   Input Parameter:
.  bv - the basis vectors context

   Output Parameter:
.  A  - the matrix

   Notes:
   The returned matrix contains only the active columns. If the content of
   the Mat is modified, these changes are also done in the BV object. The
   user must call BVRestoreMat() when no longer needed.

   This operation implies a call to BVGetArray(), which may result in data
   copies.

   Level: advanced

.seealso: BVRestoreMat(), BVCreateMat(), BVGetArray()
@*/
PetscErrorCode BVGetMat(BV bv,Mat *A)
{
  PetscErrorCode ierr;
  PetscScalar    *vv,*aa;
  PetscBool      create=PETSC_FALSE;
  PetscInt       m,cols;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscValidPointer(A,2);
  m = bv->k-bv->l;
  if (!bv->Aget) create=PETSC_TRUE;
  else {
    ierr = MatDenseGetArray(bv->Aget,&aa);CHKERRQ(ierr);
    if (aa) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"BVGetMat already called on this BV");
    ierr = MatGetSize(bv->Aget,NULL,&cols);CHKERRQ(ierr);
    if (cols!=m) {
      ierr = MatDestroy(&bv->Aget);CHKERRQ(ierr);
      create=PETSC_TRUE;
    }
  }
  ierr = BVGetArray(bv,&vv);CHKERRQ(ierr);
  if (create) {
    ierr = MatCreateDense(PetscObjectComm((PetscObject)bv),bv->n,PETSC_DECIDE,bv->N,m,vv,&bv->Aget);CHKERRQ(ierr); /* pass a pointer to avoid allocation of storage */
    ierr = MatDensePlaceArray(bv->Aget,NULL);CHKERRQ(ierr);  /* replace with a null pointer, the value after BVRestoreMat */
    ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)bv->Aget);CHKERRQ(ierr);
  }
  ierr = MatDensePlaceArray(bv->Aget,vv+(bv->nc+bv->l)*bv->n);CHKERRQ(ierr);  /* set the actual pointer */
  *A = bv->Aget;
  PetscFunctionReturn(0);
}

/*@
   BVRestoreMat - Restores the Mat obtained with BVGetMat().

   Logically Collective on bv

   Input Parameters:
+  bv - the basis vectors context
-  A  - the fetched matrix

   Note:
   A call to this function must match a previous call of BVGetMat().
   The effect is that the contents of the Mat are copied back to the
   BV internal data structures.

   Level: advanced

.seealso: BVGetMat(), BVRestoreArray()
@*/
PetscErrorCode BVRestoreMat(BV bv,Mat *A)
{
  PetscErrorCode ierr;
  PetscScalar    *vv,*aa;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  BVCheckSizes(bv,1);
  PetscValidPointer(A,2);
  if (!bv->Aget) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"BVRestoreMat must match a previous call to BVGetMat");
  if (bv->Aget!=*A) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Mat argument is not the same as the one obtained with BVGetMat");
  ierr = MatDenseGetArray(bv->Aget,&aa);CHKERRQ(ierr);
  vv = aa-(bv->nc+bv->l)*bv->n;
  ierr = MatDenseResetArray(bv->Aget);CHKERRQ(ierr);
  ierr = BVRestoreArray(bv,&vv);CHKERRQ(ierr);
  *A = NULL;
  PetscFunctionReturn(0);
}

/*
   Copy all user-provided attributes of V to another BV object W
 */
PETSC_STATIC_INLINE PetscErrorCode BVDuplicate_Private(BV V,BV W)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = BVSetType(W,((PetscObject)V)->type_name);CHKERRQ(ierr);
  W->orthog_type  = V->orthog_type;
  W->orthog_ref   = V->orthog_ref;
  W->orthog_eta   = V->orthog_eta;
  W->orthog_block = V->orthog_block;
  if (V->matrix) { ierr = PetscObjectReference((PetscObject)V->matrix);CHKERRQ(ierr); }
  W->matrix       = V->matrix;
  W->indef        = V->indef;
  W->vmm          = V->vmm;
  W->rrandom      = V->rrandom;
  if (V->rand) { ierr = PetscObjectReference((PetscObject)V->rand);CHKERRQ(ierr); }
  W->rand         = V->rand;
  if (V->ops->duplicate) { ierr = (*V->ops->duplicate)(V,W);CHKERRQ(ierr); }
  ierr = PetscObjectStateIncrease((PetscObject)W);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVDuplicate - Creates a new basis vector object of the same type and
   dimensions as an existing one.

   Collective on V

   Input Parameter:
.  V - basis vectors context

   Output Parameter:
.  W - location to put the new BV

   Notes:
   The new BV has the same type and dimensions as V, and it shares the same
   template vector. Also, the inner product matrix and orthogonalization
   options are copied.

   BVDuplicate() DOES NOT COPY the entries, but rather allocates storage
   for the new basis vectors. Use BVCopy() to copy the contents.

   Level: intermediate

.seealso: BVDuplicateResize(), BVCreate(), BVSetSizesFromVec(), BVCopy()
@*/
PetscErrorCode BVDuplicate(BV V,BV *W)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidPointer(W,2);
  ierr = BVCreate(PetscObjectComm((PetscObject)V),W);CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(*W,V->t,V->m);CHKERRQ(ierr);
  ierr = BVDuplicate_Private(V,*W);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVDuplicateResize - Creates a new basis vector object of the same type and
   dimensions as an existing one, but with possibly different number of columns.

   Collective on V

   Input Parameter:
+  V - basis vectors context
-  m - the new number of columns

   Output Parameter:
.  W - location to put the new BV

   Note:
   This is equivalent of a call to BVDuplicate() followed by BVResize(). The
   contents of V are not copied to W.

   Level: intermediate

.seealso: BVDuplicate(), BVResize()
@*/
PetscErrorCode BVDuplicateResize(BV V,PetscInt m,BV *W)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidLogicalCollectiveInt(V,m,2);
  PetscValidPointer(W,3);
  ierr = BVCreate(PetscObjectComm((PetscObject)V),W);CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(*W,V->t,m);CHKERRQ(ierr);
  ierr = BVDuplicate_Private(V,*W);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVGetCachedBV - Returns a BV object stored internally that holds the
   result of B*X after a call to BVApplyMatrixBV().

   Not collective

   Input Parameter:
.  bv    - the basis vectors context

   Output Parameter:
.  cached - the cached BV

   Note:
   The cached BV is created if not available previously.

   Level: developer

.seealso: BVApplyMatrixBV()
@*/
PetscErrorCode BVGetCachedBV(BV bv,BV *cached)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidPointer(cached,2);
  BVCheckSizes(bv,1);
  if (!bv->cached) {
    ierr = BVCreate(PetscObjectComm((PetscObject)bv),&bv->cached);CHKERRQ(ierr);
    ierr = BVSetSizesFromVec(bv->cached,bv->t,bv->m);CHKERRQ(ierr);
    ierr = BVDuplicate_Private(bv,bv->cached);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)bv->cached);CHKERRQ(ierr);
  }
  *cached = bv->cached;
  PetscFunctionReturn(0);
}

/*@
   BVCopy - Copies a basis vector object into another one, W <- V.

   Logically Collective on V

   Input Parameter:
.  V - basis vectors context

   Output Parameter:
.  W - the copy

   Note:
   Both V and W must be distributed in the same manner; local copies are
   done. Only active columns (excluding the leading ones) are copied.
   In the destination W, columns are overwritten starting from the leading ones.
   Constraints are not copied.

   Level: beginner

.seealso: BVCopyVec(), BVCopyColumn(), BVDuplicate(), BVSetActiveColumns()
@*/
PetscErrorCode BVCopy(BV V,BV W)
{
  PetscErrorCode    ierr;
  PetscScalar       *womega;
  const PetscScalar *vomega;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  BVCheckOp(V,1,copy);
  PetscValidHeaderSpecific(W,BV_CLASSID,2);
  PetscValidType(W,2);
  BVCheckSizes(W,2);
  PetscCheckSameTypeAndComm(V,1,W,2);
  if (V->n!=W->n) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Mismatching local dimension V %D, W %D",V->n,W->n);
  if (V->k-V->l>W->m-W->l) SETERRQ2(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_SIZ,"W has %D non-leading columns, not enough to store %D columns",W->m-W->l,V->k-V->l);
  if (!V->n) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(BV_Copy,V,W,0,0);CHKERRQ(ierr);
  if (V->indef && V->matrix && V->indef==W->indef && V->matrix==W->matrix) {
    /* copy signature */
    ierr = BV_AllocateSignature(W);CHKERRQ(ierr);
    ierr = VecGetArrayRead(V->omega,&vomega);CHKERRQ(ierr);
    ierr = VecGetArray(W->omega,&womega);CHKERRQ(ierr);
    ierr = PetscArraycpy(womega+W->nc+W->l,vomega+V->nc+V->l,V->k-V->l);CHKERRQ(ierr);
    ierr = VecRestoreArray(W->omega,&womega);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(V->omega,&vomega);CHKERRQ(ierr);
  }
  ierr = (*V->ops->copy)(V,W);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_Copy,V,W,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)W);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVCopyVec - Copies one of the columns of a basis vectors object into a Vec.

   Logically Collective on V

   Input Parameter:
+  V - basis vectors context
-  j - the column number to be copied

   Output Parameter:
.  w - the copied column

   Note:
   Both V and w must be distributed in the same manner; local copies are done.

   Level: beginner

.seealso: BVCopy(), BVCopyColumn(), BVInsertVec()
@*/
PetscErrorCode BVCopyVec(BV V,PetscInt j,Vec w)
{
  PetscErrorCode ierr;
  PetscInt       n,N;
  Vec            z;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidLogicalCollectiveInt(V,j,2);
  PetscValidHeaderSpecific(w,VEC_CLASSID,3);
  PetscCheckSameComm(V,1,w,3);

  ierr = VecGetSize(w,&N);CHKERRQ(ierr);
  ierr = VecGetLocalSize(w,&n);CHKERRQ(ierr);
  if (N!=V->N || n!=V->n) SETERRQ4(PetscObjectComm((PetscObject)V),PETSC_ERR_ARG_INCOMP,"Vec sizes (global %D, local %D) do not match BV sizes (global %D, local %D)",N,n,V->N,V->n);

  ierr = PetscLogEventBegin(BV_Copy,V,w,0,0);CHKERRQ(ierr);
  ierr = BVGetColumn(V,j,&z);CHKERRQ(ierr);
  ierr = VecCopy(z,w);CHKERRQ(ierr);
  ierr = BVRestoreColumn(V,j,&z);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(BV_Copy,V,w,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   BVCopyColumn - Copies the values from one of the columns to another one.

   Logically Collective on V

   Input Parameter:
+  V - basis vectors context
.  j - the number of the source column
-  i - the number of the destination column

   Level: beginner

.seealso: BVCopy(), BVCopyVec()
@*/
PetscErrorCode BVCopyColumn(BV V,PetscInt j,PetscInt i)
{
  PetscErrorCode ierr;
  Vec            z,w;
  PetscScalar    *omega;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(V,BV_CLASSID,1);
  PetscValidType(V,1);
  BVCheckSizes(V,1);
  PetscValidLogicalCollectiveInt(V,j,2);
  PetscValidLogicalCollectiveInt(V,i,3);
  if (j==i) PetscFunctionReturn(0);

  ierr = PetscLogEventBegin(BV_Copy,V,0,0,0);CHKERRQ(ierr);
  if (V->omega) {
    ierr = VecGetArray(V->omega,&omega);CHKERRQ(ierr);
    omega[i] = omega[j];
    ierr = VecRestoreArray(V->omega,&omega);CHKERRQ(ierr);
  }
  if (V->ops->copycolumn) {
    ierr = (*V->ops->copycolumn)(V,j,i);CHKERRQ(ierr);
  } else {
    ierr = BVGetColumn(V,j,&z);CHKERRQ(ierr);
    ierr = BVGetColumn(V,i,&w);CHKERRQ(ierr);
    ierr = VecCopy(z,w);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,j,&z);CHKERRQ(ierr);
    ierr = BVRestoreColumn(V,i,&w);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(BV_Copy,V,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectStateIncrease((PetscObject)V);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode BVGetSplit_Private(BV bv,PetscBool left,BV *split)
{
  PetscErrorCode ierr;
  PetscInt       ncols;

  PetscFunctionBegin;
  ncols = left? bv->nc+bv->l: bv->m-bv->l;
  if (!*split) {
    ierr = BVCreate(PetscObjectComm((PetscObject)bv),split);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)*split);CHKERRQ(ierr);
    (*split)->issplit = left? 1: 2;
    (*split)->splitparent = bv;
    ierr = BVSetSizesFromVec(*split,bv->t,ncols);CHKERRQ(ierr);
    ierr = BVDuplicate_Private(bv,*split);CHKERRQ(ierr);
  } else {
    ierr = BVResize(*split,ncols,PETSC_FALSE);CHKERRQ(ierr);
  }
  (*split)->l  = 0;
  (*split)->k  = left? bv->l: bv->k-bv->l;
  (*split)->nc = left? bv->nc: 0;
  (*split)->m  = ncols-(*split)->nc;
  if ((*split)->nc) {
    (*split)->ci[0] = -(*split)->nc-1;
    (*split)->ci[1] = -(*split)->nc-1;
  }
  if (left) {
    ierr = PetscObjectStateGet((PetscObject)*split,&bv->lstate);CHKERRQ(ierr);
  } else {
    ierr = PetscObjectStateGet((PetscObject)*split,&bv->rstate);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   BVGetSplit - Splits the BV object into two BV objects that share the
   internal data, one of them containing the leading columns and the other
   one containing the remaining columns.

   Logically Collective on bv

   Input Parameters:
.  bv - the basis vectors context

   Output Parameters:
+  L - left BV containing leading columns (can be NULL)
-  R - right BV containing remaining columns (can be NULL)

   Notes:
   The columns are split in two sets. The leading columns (including the
   constraints) are assigned to the left BV and the remaining columns
   are assigned to the right BV. The number of leading columns, as
   specified with BVSetActiveColumns(), must be between 1 and m-1 (to
   guarantee that both L and R have at least one column).

   The returned BV's must be seen as references (not copies) of the input
   BV, that is, modifying them will change the entries of bv as well.
   The returned BV's must not be destroyed. BVRestoreSplit() must be called
   when they are no longer needed.

   Pass NULL for any of the output BV's that is not needed.

   Level: advanced

.seealso: BVRestoreSplit(), BVSetActiveColumns(), BVSetNumConstraints()
@*/
PetscErrorCode BVGetSplit(BV bv,BV *L,BV *R)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  if (!bv->l) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Must indicate the number of leading columns with BVSetActiveColumns()");
  if (bv->lsplit) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Cannot get the split BV's twice before restoring them with BVRestoreSplit()");
  bv->lsplit = bv->nc+bv->l;
  ierr = BVGetSplit_Private(bv,PETSC_TRUE,&bv->L);CHKERRQ(ierr);
  ierr = BVGetSplit_Private(bv,PETSC_FALSE,&bv->R);CHKERRQ(ierr);
  if (L) *L = bv->L;
  if (R) *R = bv->R;
  PetscFunctionReturn(0);
}

/*@
   BVRestoreSplit - Restore the BV objects obtained with BVGetSplit().

   Logically Collective on bv

   Input Parameters:
+  bv - the basis vectors context
.  L  - left BV obtained with BVGetSplit()
-  R  - right BV obtained with BVGetSplit()

   Note:
   The arguments must match the corresponding call to BVGetSplit().

   Level: advanced

.seealso: BVGetSplit()
@*/
PetscErrorCode BVRestoreSplit(BV bv,BV *L,BV *R)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(bv,BV_CLASSID,1);
  PetscValidType(bv,1);
  BVCheckSizes(bv,1);
  if (!bv->lsplit) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Must call BVGetSplit first");
  if (L && *L!=bv->L) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Argument 2 is not the same BV that was obtained with BVGetSplit");
  if (R && *R!=bv->R) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONG,"Argument 3 is not the same BV that was obtained with BVGetSplit");
  if (L && ((*L)->ci[0]>(*L)->nc-1 || (*L)->ci[1]>(*L)->nc-1)) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Argument 2 has unrestored columns, use BVRestoreColumn()");
  if (R && ((*R)->ci[0]>(*R)->nc-1 || (*R)->ci[1]>(*R)->nc-1)) SETERRQ(PetscObjectComm((PetscObject)bv),PETSC_ERR_ARG_WRONGSTATE,"Argument 3 has unrestored columns, use BVRestoreColumn()");

  if (bv->ops->restoresplit) {
    ierr = (*bv->ops->restoresplit)(bv,L,R);CHKERRQ(ierr);
  }
  bv->lsplit = 0;
  if (L) *L = NULL;
  if (R) *R = NULL;
  PetscFunctionReturn(0);
}

