/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/slepcimpl.h>            /*I "slepcsys.h" I*/

static PetscErrorCode MatCreateTile_SeqAIJ(PetscScalar a,Mat A,PetscScalar b,Mat B,PetscScalar c,Mat C,PetscScalar d,Mat D,Mat G)
{
  PetscErrorCode    ierr;
  PetscInt          i,j,M1,M2,N1,N2,*nnz,ncols,*scols,bs;
  PetscScalar       *svals,*buf;
  const PetscInt    *cols;
  const PetscScalar *vals;

  PetscFunctionBegin;
  ierr = MatGetSize(A,&M1,&N1);CHKERRQ(ierr);
  ierr = MatGetSize(D,&M2,&N2);CHKERRQ(ierr);
  ierr = MatGetBlockSize(A,&bs);CHKERRQ(ierr);

  ierr = PetscCalloc1((M1+M2)/bs,&nnz);CHKERRQ(ierr);
  /* Preallocate for A */
  if (a!=0.0) {
    for (i=0;i<(M1+bs-1)/bs;i++) {
      ierr = MatGetRow(A,i*bs,&ncols,NULL,NULL);CHKERRQ(ierr);
      nnz[i] += ncols/bs;
      ierr = MatRestoreRow(A,i*bs,&ncols,NULL,NULL);CHKERRQ(ierr);
    }
  }
  /* Preallocate for B */
  if (b!=0.0) {
    for (i=0;i<(M1+bs-1)/bs;i++) {
      ierr = MatGetRow(B,i*bs,&ncols,NULL,NULL);CHKERRQ(ierr);
      nnz[i] += ncols/bs;
      ierr = MatRestoreRow(B,i*bs,&ncols,NULL,NULL);CHKERRQ(ierr);
    }
  }
  /* Preallocate for C */
  if (c!=0.0) {
    for (i=0;i<(M2+bs-1)/bs;i++) {
      ierr = MatGetRow(C,i*bs,&ncols,NULL,NULL);CHKERRQ(ierr);
      nnz[i+M1/bs] += ncols/bs;
      ierr = MatRestoreRow(C,i*bs,&ncols,NULL,NULL);CHKERRQ(ierr);
    }
  }
  /* Preallocate for D */
  if (d!=0.0) {
    for (i=0;i<(M2+bs-1)/bs;i++) {
      ierr = MatGetRow(D,i*bs,&ncols,NULL,NULL);CHKERRQ(ierr);
      nnz[i+M1/bs] += ncols/bs;
      ierr = MatRestoreRow(D,i*bs,&ncols,NULL,NULL);CHKERRQ(ierr);
    }
  }
  ierr = MatXAIJSetPreallocation(G,bs,nnz,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscFree(nnz);CHKERRQ(ierr);

  ierr = PetscMalloc2(PetscMax(N1,N2),&buf,PetscMax(N1,N2),&scols);CHKERRQ(ierr);
  /* Transfer A */
  if (a!=0.0) {
    for (i=0;i<M1;i++) {
      ierr = MatGetRow(A,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      if (a!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*a;
      } else svals=(PetscScalar*)vals;
      ierr = MatSetValues(G,1,&i,ncols,cols,svals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatRestoreRow(A,i,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }
  /* Transfer B */
  if (b!=0.0) {
    for (i=0;i<M1;i++) {
      ierr = MatGetRow(B,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      if (b!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*b;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = cols[j]+N1;
      ierr = MatSetValues(G,1,&i,ncols,scols,svals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatRestoreRow(B,i,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }
  /* Transfer C */
  if (c!=0.0) {
    for (i=0;i<M2;i++) {
      ierr = MatGetRow(C,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      if (c!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*c;
      } else svals=(PetscScalar*)vals;
      j = i+M1;
      ierr = MatSetValues(G,1,&j,ncols,cols,svals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatRestoreRow(C,i,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }
  /* Transfer D */
  if (d!=0.0) {
    for (i=0;i<M2;i++) {
      ierr = MatGetRow(D,i,&ncols,&cols,&vals);CHKERRQ(ierr);
      if (d!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*d;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = cols[j]+N1;
      j = i+M1;
      ierr = MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatRestoreRow(D,i,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }
  ierr = PetscFree2(buf,scols);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatCreateTile_MPIAIJ(PetscScalar a,Mat A,PetscScalar b,Mat B,PetscScalar c,Mat C,PetscScalar d,Mat D,Mat G)
{
  PetscErrorCode ierr;
  PetscMPIInt       np;
  PetscInt          p,i,j,N1,N2,m1,m2,n1,n2,*map1,*map2;
  PetscInt          *dnz,*onz,ncols,*scols,start,gstart;
  PetscScalar       *svals,*buf;
  const PetscInt    *cols,*mapptr1,*mapptr2;
  const PetscScalar *vals;

  PetscFunctionBegin;
  ierr = MatGetSize(A,NULL,&N1);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&m1,&n1);CHKERRQ(ierr);
  ierr = MatGetSize(D,NULL,&N2);CHKERRQ(ierr);
  ierr = MatGetLocalSize(D,&m2,&n2);CHKERRQ(ierr);

  /* Create mappings */
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject)G),&np);CHKERRQ(ierr);
  ierr = MatGetOwnershipRangesColumn(A,&mapptr1);CHKERRQ(ierr);
  ierr = MatGetOwnershipRangesColumn(B,&mapptr2);CHKERRQ(ierr);
  ierr = PetscMalloc4(PetscMax(N1,N2),&buf,PetscMax(N1,N2),&scols,N1,&map1,N2,&map2);CHKERRQ(ierr);
  for (p=0;p<np;p++) {
    for (i=mapptr1[p];i<mapptr1[p+1];i++) map1[i] = i+mapptr2[p];
  }
  for (p=0;p<np;p++) {
    for (i=mapptr2[p];i<mapptr2[p+1];i++) map2[i] = i+mapptr1[p+1];
  }

  ierr = MatPreallocateInitialize(PetscObjectComm((PetscObject)G),m1+m2,n1+n2,dnz,onz);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(G,&gstart,NULL);CHKERRQ(ierr);
  /* Preallocate for A */
  if (a!=0.0) {
    ierr = MatGetOwnershipRange(A,&start,NULL);CHKERRQ(ierr);
    for (i=0;i<m1;i++) {
      ierr = MatGetRow(A,i+start,&ncols,&cols,NULL);CHKERRQ(ierr);
      for (j=0;j<ncols;j++) scols[j] = map1[cols[j]];
      ierr = MatPreallocateSet(gstart+i,ncols,scols,dnz,onz);CHKERRQ(ierr);
      ierr = MatRestoreRow(A,i+start,&ncols,&cols,NULL);CHKERRQ(ierr);
    }
  }
  /* Preallocate for B */
  if (b!=0.0) {
    ierr = MatGetOwnershipRange(B,&start,NULL);CHKERRQ(ierr);
    for (i=0;i<m1;i++) {
      ierr = MatGetRow(B,i+start,&ncols,&cols,NULL);CHKERRQ(ierr);
      for (j=0;j<ncols;j++) scols[j] = map2[cols[j]];
      ierr = MatPreallocateSet(gstart+i,ncols,scols,dnz,onz);CHKERRQ(ierr);
      ierr = MatRestoreRow(B,i+start,&ncols,&cols,NULL);CHKERRQ(ierr);
    }
  }
  /* Preallocate for C */
  if (c!=0.0) {
    ierr = MatGetOwnershipRange(C,&start,NULL);CHKERRQ(ierr);
    for (i=0;i<m2;i++) {
      ierr = MatGetRow(C,i+start,&ncols,&cols,NULL);CHKERRQ(ierr);
      for (j=0;j<ncols;j++) scols[j] = map1[cols[j]];
      ierr = MatPreallocateSet(gstart+m1+i,ncols,scols,dnz,onz);CHKERRQ(ierr);
      ierr = MatRestoreRow(C,i+start,&ncols,&cols,NULL);CHKERRQ(ierr);
    }
  }
  /* Preallocate for D */
  if (d!=0.0) {
    ierr = MatGetOwnershipRange(D,&start,NULL);CHKERRQ(ierr);
    for (i=0;i<m2;i++) {
      ierr = MatGetRow(D,i+start,&ncols,&cols,NULL);CHKERRQ(ierr);
      for (j=0;j<ncols;j++) scols[j] = map2[cols[j]];
      ierr = MatPreallocateSet(gstart+m1+i,ncols,scols,dnz,onz);CHKERRQ(ierr);
      ierr = MatRestoreRow(D,i+start,&ncols,&cols,NULL);CHKERRQ(ierr);
    }
  }
  ierr = MatMPIAIJSetPreallocation(G,0,dnz,0,onz);CHKERRQ(ierr);
  ierr = MatPreallocateFinalize(dnz,onz);CHKERRQ(ierr);

  /* Transfer A */
  if (a!=0.0) {
    ierr = MatGetOwnershipRange(A,&start,NULL);CHKERRQ(ierr);
    for (i=0;i<m1;i++) {
      ierr = MatGetRow(A,i+start,&ncols,&cols,&vals);CHKERRQ(ierr);
      if (a!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*a;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = map1[cols[j]];
      j = gstart+i;
      ierr = MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatRestoreRow(A,i+start,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }
  /* Transfer B */
  if (b!=0.0) {
    ierr = MatGetOwnershipRange(B,&start,NULL);CHKERRQ(ierr);
    for (i=0;i<m1;i++) {
      ierr = MatGetRow(B,i+start,&ncols,&cols,&vals);CHKERRQ(ierr);
      if (b!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*b;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = map2[cols[j]];
      j = gstart+i;
      ierr = MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatRestoreRow(B,i+start,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }
  /* Transfer C */
  if (c!=0.0) {
    ierr = MatGetOwnershipRange(C,&start,NULL);CHKERRQ(ierr);
    for (i=0;i<m2;i++) {
      ierr = MatGetRow(C,i+start,&ncols,&cols,&vals);CHKERRQ(ierr);
      if (c!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*c;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = map1[cols[j]];
      j = gstart+m1+i;
      ierr = MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatRestoreRow(C,i+start,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }
  /* Transfer D */
  if (d!=0.0) {
    ierr = MatGetOwnershipRange(D,&start,NULL);CHKERRQ(ierr);
    for (i=0;i<m2;i++) {
      ierr = MatGetRow(D,i+start,&ncols,&cols,&vals);CHKERRQ(ierr);
      if (d!=1.0) {
        svals=buf;
        for (j=0;j<ncols;j++) svals[j] = vals[j]*d;
      } else svals=(PetscScalar*)vals;
      for (j=0;j<ncols;j++) scols[j] = map2[cols[j]];
      j = gstart+m1+i;
      ierr = MatSetValues(G,1,&j,ncols,scols,svals,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatRestoreRow(D,i+start,&ncols,&cols,&vals);CHKERRQ(ierr);
    }
  }
  ierr = PetscFree4(buf,scols,map1,map2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   MatCreateTile - Explicitly build a matrix from four blocks, G = [ a*A b*B; c*C d*D ].

   Collective on A

   Input parameters:
+  A - matrix for top-left block
.  a - scaling factor for block A
.  B - matrix for top-right block
.  b - scaling factor for block B
.  C - matrix for bottom-left block
.  c - scaling factor for block C
.  D - matrix for bottom-right block
-  d - scaling factor for block D

   Output parameter:
.  G  - the resulting matrix

   Notes:
   In the case of a parallel matrix, a permuted version of G is returned. The permutation
   is a perfect shuffle such that the local parts of A, B, C, D remain in the local part of
   G for the same process.

   Matrix G must be destroyed by the user.

   Level: developer
@*/
PetscErrorCode MatCreateTile(PetscScalar a,Mat A,PetscScalar b,Mat B,PetscScalar c,Mat C,PetscScalar d,Mat D,Mat *G)
{
  PetscErrorCode ierr;
  PetscInt       M1,M2,N1,N2,M,N,m1,m2,n1,n2,m,n,bs;
  PetscBool      flg1;
  MatType        Atype;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(A,MAT_CLASSID,2);
  PetscValidHeaderSpecific(B,MAT_CLASSID,4);
  PetscValidHeaderSpecific(C,MAT_CLASSID,6);
  PetscValidHeaderSpecific(D,MAT_CLASSID,8);
  PetscCheckSameTypeAndComm(A,2,B,4);
  PetscCheckSameTypeAndComm(A,2,C,6);
  PetscCheckSameTypeAndComm(A,2,D,8);
  PetscValidLogicalCollectiveScalar(A,a,1);
  PetscValidLogicalCollectiveScalar(A,b,3);
  PetscValidLogicalCollectiveScalar(A,c,5);
  PetscValidLogicalCollectiveScalar(A,d,7);
  PetscValidPointer(G,9);

  /* check row 1 */
  ierr = MatGetSize(A,&M1,NULL);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,&m1,NULL);CHKERRQ(ierr);
  ierr = MatGetSize(B,&M,NULL);CHKERRQ(ierr);
  ierr = MatGetLocalSize(B,&m,NULL);CHKERRQ(ierr);
  if (M!=M1 || m!=m1) SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible dimensions");
  /* check row 2 */
  ierr = MatGetSize(C,&M2,NULL);CHKERRQ(ierr);
  ierr = MatGetLocalSize(C,&m2,NULL);CHKERRQ(ierr);
  ierr = MatGetSize(D,&M,NULL);CHKERRQ(ierr);
  ierr = MatGetLocalSize(D,&m,NULL);CHKERRQ(ierr);
  if (M!=M2 || m!=m2) SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible dimensions");
  /* check column 1 */
  ierr = MatGetSize(A,NULL,&N1);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A,NULL,&n1);CHKERRQ(ierr);
  ierr = MatGetSize(C,NULL,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(C,NULL,&n);CHKERRQ(ierr);
  if (N!=N1 || n!=n1) SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible dimensions");
  /* check column 2 */
  ierr = MatGetSize(B,NULL,&N2);CHKERRQ(ierr);
  ierr = MatGetLocalSize(B,NULL,&n2);CHKERRQ(ierr);
  ierr = MatGetSize(D,NULL,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(D,NULL,&n);CHKERRQ(ierr);
  if (N!=N2 || n!=n2) SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_ARG_INCOMP,"Incompatible dimensions");

  ierr = MatGetType(A,&Atype);CHKERRQ(ierr);
  ierr = MatGetBlockSize(A,&bs);CHKERRQ(ierr);
  ierr = MatCreate(PetscObjectComm((PetscObject)A),G);CHKERRQ(ierr);
  ierr = MatSetSizes(*G,m1+m2,n1+n2,M1+M2,N1+N2);CHKERRQ(ierr);
  ierr = MatSetType(*G,Atype);CHKERRQ(ierr);
  ierr = MatSetBlockSize(*G,bs);CHKERRQ(ierr);
  ierr = MatSetUp(*G);CHKERRQ(ierr);

  ierr = PetscObjectTypeCompareAny((PetscObject)A,&flg1,MATMPIAIJ,MATMPIAIJCUSPARSE,"");CHKERRQ(ierr);
  if (flg1) {
    ierr = MatCreateTile_MPIAIJ(a,A,b,B,c,C,d,D,*G);CHKERRQ(ierr);
  } else {
    ierr = PetscObjectTypeCompareAny((PetscObject)A,&flg1,MATSEQAIJ,MATSEQAIJCUSPARSE,MATSEQBAIJ,"");CHKERRQ(ierr);
    if (flg1) {
      ierr = MatCreateTile_SeqAIJ(a,A,b,B,c,C,d,D,*G);CHKERRQ(ierr);
    } else SETERRQ(PetscObjectComm((PetscObject)A),PETSC_ERR_SUP,"Not implemented for this matrix type");
  }
  ierr = MatAssemblyBegin(*G,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*G,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   MatCreateVecsEmpty - Get vector(s) compatible with the matrix, i.e. with the same
   parallel layout, but without internal array.

   Collective on mat

   Input Parameter:
.  mat - the matrix

   Output Parameters:
+  right - (optional) vector that the matrix can be multiplied against
-  left - (optional) vector that the matrix vector product can be stored in

   Note:
   This is similar to MatCreateVecs(), but the new vectors do not have an internal
   array, so the intended usage is with VecPlaceArray().

   Level: developer
@*/
PetscErrorCode MatCreateVecsEmpty(Mat mat,Vec *right,Vec *left)
{
  PetscErrorCode ierr;
  PetscBool      standard,cuda=PETSC_FALSE,skip=PETSC_FALSE;
  PetscInt       M,N,mloc,nloc,rbs,cbs;
  PetscMPIInt    size;
  Vec            v;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(mat,MAT_CLASSID,1);
  PetscValidType(mat,1);

  ierr = PetscObjectTypeCompareAny((PetscObject)mat,&standard,MATSEQAIJ,MATMPIAIJ,MATSEQBAIJ,MATMPIBAIJ,MATSEQSBAIJ,MATMPISBAIJ,MATSEQDENSE,MATMPIDENSE,"");CHKERRQ(ierr);
  ierr = PetscObjectTypeCompareAny((PetscObject)mat,&cuda,MATSEQAIJCUSPARSE,MATMPIAIJCUSPARSE,"");CHKERRQ(ierr);
  if (!standard && !cuda) {
    ierr = MatCreateVecs(mat,right,left);CHKERRQ(ierr);
    v = right? *right: *left;
    if (v) {
      ierr = PetscObjectTypeCompareAny((PetscObject)v,&standard,VECSEQ,VECMPI,"");CHKERRQ(ierr);
      ierr = PetscObjectTypeCompareAny((PetscObject)v,&cuda,VECSEQCUDA,VECMPICUDA,"");CHKERRQ(ierr);
    }
    if (!standard && !cuda) skip = PETSC_TRUE;
    else {
      if (right) { ierr = VecDestroy(right);CHKERRQ(ierr); }
      if (left) { ierr = VecDestroy(left);CHKERRQ(ierr); }
    }
  }
  if (!skip) {
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)mat),&size);CHKERRQ(ierr);
    ierr = MatGetLocalSize(mat,&mloc,&nloc);CHKERRQ(ierr);
    ierr = MatGetSize(mat,&M,&N);CHKERRQ(ierr);
    ierr = MatGetBlockSizes(mat,&rbs,&cbs);CHKERRQ(ierr);
    if (right) {
      if (cuda) {
#if defined(PETSC_HAVE_CUDA)
        if (size>1) {
          ierr = VecCreateMPICUDAWithArray(PetscObjectComm((PetscObject)mat),cbs,nloc,N,NULL,right);CHKERRQ(ierr);
        } else {
          ierr = VecCreateSeqCUDAWithArray(PetscObjectComm((PetscObject)mat),cbs,N,NULL,right);CHKERRQ(ierr);
        }
#endif
      } else {
        if (size>1) {
          ierr = VecCreateMPIWithArray(PetscObjectComm((PetscObject)mat),cbs,nloc,N,NULL,right);CHKERRQ(ierr);
        } else {
          ierr = VecCreateSeqWithArray(PetscObjectComm((PetscObject)mat),cbs,N,NULL,right);CHKERRQ(ierr);
        }
      }
    }
    if (left) {
      if (cuda) {
#if defined(PETSC_HAVE_CUDA)
        if (size>1) {
          ierr = VecCreateMPICUDAWithArray(PetscObjectComm((PetscObject)mat),rbs,mloc,M,NULL,left);CHKERRQ(ierr);
        } else {
          ierr = VecCreateSeqCUDAWithArray(PetscObjectComm((PetscObject)mat),rbs,M,NULL,left);CHKERRQ(ierr);
        }
#endif
      } else {
        if (size>1) {
          ierr = VecCreateMPIWithArray(PetscObjectComm((PetscObject)mat),rbs,mloc,M,NULL,left);CHKERRQ(ierr);
        } else {
          ierr = VecCreateSeqWithArray(PetscObjectComm((PetscObject)mat),rbs,M,NULL,left);CHKERRQ(ierr);
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

