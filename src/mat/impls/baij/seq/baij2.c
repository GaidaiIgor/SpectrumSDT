
#include <../src/mat/impls/baij/seq/baij.h>
#include <../src/mat/impls/dense/seq/dense.h>
#include <petsc/private/kernels/blockinvert.h>
#include <petscbt.h>
#include <petscblaslapack.h>

#if defined(PETSC_HAVE_IMMINTRIN_H) && defined(__AVX2__) && defined(__FMA__) && defined(PETSC_USE_REAL_DOUBLE) && !defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_64BIT_INDICES)
#include <immintrin.h>
#endif

PetscErrorCode MatIncreaseOverlap_SeqBAIJ(Mat A,PetscInt is_max,IS is[],PetscInt ov)
{
  Mat_SeqBAIJ    *a = (Mat_SeqBAIJ*)A->data;
  PetscErrorCode ierr;
  PetscInt       row,i,j,k,l,m,n,*nidx,isz,val,ival;
  const PetscInt *idx;
  PetscInt       start,end,*ai,*aj,bs,*nidx2;
  PetscBT        table;

  PetscFunctionBegin;
  m  = a->mbs;
  ai = a->i;
  aj = a->j;
  bs = A->rmap->bs;

  if (ov < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Negative overlap specified");

  ierr = PetscBTCreate(m,&table);CHKERRQ(ierr);
  ierr = PetscMalloc1(m+1,&nidx);CHKERRQ(ierr);
  ierr = PetscMalloc1(A->rmap->N+1,&nidx2);CHKERRQ(ierr);

  for (i=0; i<is_max; i++) {
    /* Initialise the two local arrays */
    isz  = 0;
    ierr = PetscBTMemzero(m,table);CHKERRQ(ierr);

    /* Extract the indices, assume there can be duplicate entries */
    ierr = ISGetIndices(is[i],&idx);CHKERRQ(ierr);
    ierr = ISGetLocalSize(is[i],&n);CHKERRQ(ierr);

    /* Enter these into the temp arrays i.e mark table[row], enter row into new index */
    for (j=0; j<n; ++j) {
      ival = idx[j]/bs; /* convert the indices into block indices */
      if (ival>=m) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"index greater than mat-dim");
      if (!PetscBTLookupSet(table,ival)) nidx[isz++] = ival;
    }
    ierr = ISRestoreIndices(is[i],&idx);CHKERRQ(ierr);
    ierr = ISDestroy(&is[i]);CHKERRQ(ierr);

    k = 0;
    for (j=0; j<ov; j++) { /* for each overlap*/
      n = isz;
      for (; k<n; k++) {  /* do only those rows in nidx[k], which are not done yet */
        row   = nidx[k];
        start = ai[row];
        end   = ai[row+1];
        for (l = start; l<end; l++) {
          val = aj[l];
          if (!PetscBTLookupSet(table,val)) nidx[isz++] = val;
        }
      }
    }
    /* expand the Index Set */
    for (j=0; j<isz; j++) {
      for (k=0; k<bs; k++) nidx2[j*bs+k] = nidx[j]*bs+k;
    }
    ierr = ISCreateGeneral(PETSC_COMM_SELF,isz*bs,nidx2,PETSC_COPY_VALUES,is+i);CHKERRQ(ierr);
  }
  ierr = PetscBTDestroy(&table);CHKERRQ(ierr);
  ierr = PetscFree(nidx);CHKERRQ(ierr);
  ierr = PetscFree(nidx2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatCreateSubMatrix_SeqBAIJ_Private(Mat A,IS isrow,IS iscol,MatReuse scall,Mat *B)
{
  Mat_SeqBAIJ    *a = (Mat_SeqBAIJ*)A->data,*c;
  PetscErrorCode ierr;
  PetscInt       *smap,i,k,kstart,kend,oldcols = a->nbs,*lens;
  PetscInt       row,mat_i,*mat_j,tcol,*mat_ilen;
  const PetscInt *irow,*icol;
  PetscInt       nrows,ncols,*ssmap,bs=A->rmap->bs,bs2=a->bs2;
  PetscInt       *aj = a->j,*ai = a->i;
  MatScalar      *mat_a;
  Mat            C;
  PetscBool      flag;

  PetscFunctionBegin;
  ierr = ISGetIndices(isrow,&irow);CHKERRQ(ierr);
  ierr = ISGetIndices(iscol,&icol);CHKERRQ(ierr);
  ierr = ISGetLocalSize(isrow,&nrows);CHKERRQ(ierr);
  ierr = ISGetLocalSize(iscol,&ncols);CHKERRQ(ierr);

  ierr  = PetscCalloc1(1+oldcols,&smap);CHKERRQ(ierr);
  ssmap = smap;
  ierr  = PetscMalloc1(1+nrows,&lens);CHKERRQ(ierr);
  for (i=0; i<ncols; i++) smap[icol[i]] = i+1;
  /* determine lens of each row */
  for (i=0; i<nrows; i++) {
    kstart  = ai[irow[i]];
    kend    = kstart + a->ilen[irow[i]];
    lens[i] = 0;
    for (k=kstart; k<kend; k++) {
      if (ssmap[aj[k]]) lens[i]++;
    }
  }
  /* Create and fill new matrix */
  if (scall == MAT_REUSE_MATRIX) {
    c = (Mat_SeqBAIJ*)((*B)->data);

    if (c->mbs!=nrows || c->nbs!=ncols || (*B)->rmap->bs!=bs) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Submatrix wrong size");
    ierr = PetscArraycmp(c->ilen,lens,c->mbs,&flag);CHKERRQ(ierr);
    if (!flag) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Cannot reuse matrix. wrong no of nonzeros");
    ierr = PetscArrayzero(c->ilen,c->mbs);CHKERRQ(ierr);
    C    = *B;
  } else {
    ierr = MatCreate(PetscObjectComm((PetscObject)A),&C);CHKERRQ(ierr);
    ierr = MatSetSizes(C,nrows*bs,ncols*bs,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
    ierr = MatSetType(C,((PetscObject)A)->type_name);CHKERRQ(ierr);
    ierr = MatSeqBAIJSetPreallocation(C,bs,0,lens);CHKERRQ(ierr);
  }
  c = (Mat_SeqBAIJ*)(C->data);
  for (i=0; i<nrows; i++) {
    row      = irow[i];
    kstart   = ai[row];
    kend     = kstart + a->ilen[row];
    mat_i    = c->i[i];
    mat_j    = c->j + mat_i;
    mat_a    = c->a + mat_i*bs2;
    mat_ilen = c->ilen + i;
    for (k=kstart; k<kend; k++) {
      if ((tcol=ssmap[a->j[k]])) {
        *mat_j++ = tcol - 1;
        ierr     = PetscArraycpy(mat_a,a->a+k*bs2,bs2);CHKERRQ(ierr);
        mat_a   += bs2;
        (*mat_ilen)++;
      }
    }
  }
  /* sort */
  {
    MatScalar *work;
    ierr = PetscMalloc1(bs2,&work);CHKERRQ(ierr);
    for (i=0; i<nrows; i++) {
      PetscInt ilen;
      mat_i = c->i[i];
      mat_j = c->j + mat_i;
      mat_a = c->a + mat_i*bs2;
      ilen  = c->ilen[i];
      ierr  = PetscSortIntWithDataArray(ilen,mat_j,mat_a,bs2*sizeof(MatScalar),work);CHKERRQ(ierr);
    }
    ierr = PetscFree(work);CHKERRQ(ierr);
  }

  /* Free work space */
  ierr = ISRestoreIndices(iscol,&icol);CHKERRQ(ierr);
  ierr = PetscFree(smap);CHKERRQ(ierr);
  ierr = PetscFree(lens);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = ISRestoreIndices(isrow,&irow);CHKERRQ(ierr);
  *B   = C;
  PetscFunctionReturn(0);
}

PetscErrorCode MatCreateSubMatrix_SeqBAIJ(Mat A,IS isrow,IS iscol,MatReuse scall,Mat *B)
{
  Mat_SeqBAIJ    *a = (Mat_SeqBAIJ*)A->data;
  IS             is1,is2;
  PetscErrorCode ierr;
  PetscInt       *vary,*iary,nrows,ncols,i,bs=A->rmap->bs,count,maxmnbs,j;
  const PetscInt *irow,*icol;

  PetscFunctionBegin;
  ierr = ISGetIndices(isrow,&irow);CHKERRQ(ierr);
  ierr = ISGetIndices(iscol,&icol);CHKERRQ(ierr);
  ierr = ISGetLocalSize(isrow,&nrows);CHKERRQ(ierr);
  ierr = ISGetLocalSize(iscol,&ncols);CHKERRQ(ierr);

  /* Verify if the indices corespond to each element in a block
   and form the IS with compressed IS */
  maxmnbs = PetscMax(a->mbs,a->nbs);
  ierr = PetscMalloc2(maxmnbs,&vary,maxmnbs,&iary);CHKERRQ(ierr);
  ierr = PetscArrayzero(vary,a->mbs);CHKERRQ(ierr);
  for (i=0; i<nrows; i++) vary[irow[i]/bs]++;
  for (i=0; i<a->mbs; i++) {
    if (vary[i]!=0 && vary[i]!=bs) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Index set does not match blocks");
  }
  count = 0;
  for (i=0; i<nrows; i++) {
    j = irow[i] / bs;
    if ((vary[j]--)==bs) iary[count++] = j;
  }
  ierr = ISCreateGeneral(PETSC_COMM_SELF,count,iary,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);

  ierr = PetscArrayzero(vary,a->nbs);CHKERRQ(ierr);
  for (i=0; i<ncols; i++) vary[icol[i]/bs]++;
  for (i=0; i<a->nbs; i++) {
    if (vary[i]!=0 && vary[i]!=bs) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Internal error in PETSc");
  }
  count = 0;
  for (i=0; i<ncols; i++) {
    j = icol[i] / bs;
    if ((vary[j]--)==bs) iary[count++] = j;
  }
  ierr = ISCreateGeneral(PETSC_COMM_SELF,count,iary,PETSC_COPY_VALUES,&is2);CHKERRQ(ierr);
  ierr = ISRestoreIndices(isrow,&irow);CHKERRQ(ierr);
  ierr = ISRestoreIndices(iscol,&icol);CHKERRQ(ierr);
  ierr = PetscFree2(vary,iary);CHKERRQ(ierr);

  ierr = MatCreateSubMatrix_SeqBAIJ_Private(A,is1,is2,scall,B);CHKERRQ(ierr);
  ierr = ISDestroy(&is1);CHKERRQ(ierr);
  ierr = ISDestroy(&is2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatDestroySubMatrix_SeqBAIJ(Mat C)
{
  PetscErrorCode ierr;
  Mat_SeqBAIJ    *c = (Mat_SeqBAIJ*)C->data;
  Mat_SubSppt    *submatj = c->submatis1;

  PetscFunctionBegin;
  ierr = (*submatj->destroy)(C);CHKERRQ(ierr);
  ierr = MatDestroySubMatrix_Private(submatj);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatDestroySubMatrices_SeqBAIJ(PetscInt n,Mat *mat[])
{
  PetscErrorCode ierr;
  PetscInt       i;
  Mat            C;
  Mat_SeqBAIJ    *c;
  Mat_SubSppt    *submatj;

  PetscFunctionBegin;
  for (i=0; i<n; i++) {
    C       = (*mat)[i];
    c       = (Mat_SeqBAIJ*)C->data;
    submatj = c->submatis1;
    if (submatj) {
      ierr = (*submatj->destroy)(C);CHKERRQ(ierr);
      ierr = MatDestroySubMatrix_Private(submatj);CHKERRQ(ierr);
      ierr = PetscFree(C->defaultvectype);CHKERRQ(ierr);
      ierr = PetscLayoutDestroy(&C->rmap);CHKERRQ(ierr);
      ierr = PetscLayoutDestroy(&C->cmap);CHKERRQ(ierr);
      ierr = PetscHeaderDestroy(&C);CHKERRQ(ierr);
    } else {
      ierr = MatDestroy(&C);CHKERRQ(ierr);
    }
  }
  ierr = PetscFree(*mat);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatCreateSubMatrices_SeqBAIJ(Mat A,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *B[])
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  if (scall == MAT_INITIAL_MATRIX) {
    ierr = PetscCalloc1(n+1,B);CHKERRQ(ierr);
  }

  for (i=0; i<n; i++) {
    ierr = MatCreateSubMatrix_SeqBAIJ(A,irow[i],icol[i],scall,&(*B)[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* -------------------------------------------------------*/
/* Should check that shapes of vectors and matrices match */
/* -------------------------------------------------------*/

PetscErrorCode MatMult_SeqBAIJ_1(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z,sum;
  const PetscScalar *x;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,n;
  const PetscInt    *idx,*ii,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&z);CHKERRQ(ierr);

  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(z,a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
  }

  for (i=0; i<mbs; i++) {
    n   = ii[1] - ii[0];
    v   = a->a + ii[0];
    idx = a->j + ii[0];
    ii++;
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+1*n,1*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    sum = 0.0;
    PetscSparseDensePlusDot(sum,x,v,idx,n);
    if (usecprow) {
      z[ridx[i]] = sum;
    } else {
      z[i]       = sum;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&z);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*a->nz - a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMult_SeqBAIJ_2(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,*zarray;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,*idx,*ii,j,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,2*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    sum1        = 0.0; sum2 = 0.0;
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+4*n,4*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 2*(*idx++); x1 = xb[0]; x2 = xb[1];
      sum1 += v[0]*x1 + v[2]*x2;
      sum2 += v[1]*x1 + v[3]*x2;
      v    += 4;
    }
    if (usecprow) z = zarray + 2*ridx[i];
    z[0] = sum1; z[1] = sum2;
    if (!usecprow) z += 2;
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(8.0*a->nz - 2.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMult_SeqBAIJ_3(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,x1,x2,x3,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,*idx,*ii,j,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

#if defined(PETSC_HAVE_PRAGMA_DISJOINT)
#pragma disjoint(*v,*z,*xb)
#endif

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,3*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    sum1        = 0.0; sum2 = 0.0; sum3 = 0.0;
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+9*n,9*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb = x + 3*(*idx++);
      x1 = xb[0];
      x2 = xb[1];
      x3 = xb[2];

      sum1 += v[0]*x1 + v[3]*x2 + v[6]*x3;
      sum2 += v[1]*x1 + v[4]*x2 + v[7]*x3;
      sum3 += v[2]*x1 + v[5]*x2 + v[8]*x3;
      v    += 9;
    }
    if (usecprow) z = zarray + 3*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3;
    if (!usecprow) z += 3;
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(18.0*a->nz - 3.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMult_SeqBAIJ_4(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4,x1,x2,x3,x4,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,*idx,*ii,j,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,4*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0];
    ii++;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;

    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);     /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+16*n,16*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 4*(*idx++);
      x1    = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];
      sum1 += v[0]*x1 + v[4]*x2 + v[8]*x3  + v[12]*x4;
      sum2 += v[1]*x1 + v[5]*x2 + v[9]*x3  + v[13]*x4;
      sum3 += v[2]*x1 + v[6]*x2 + v[10]*x3 + v[14]*x4;
      sum4 += v[3]*x1 + v[7]*x2 + v[11]*x3 + v[15]*x4;
      v    += 16;
    }
    if (usecprow) z = zarray + 4*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4;
    if (!usecprow) z += 4;
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(32.0*a->nz - 4.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMult_SeqBAIJ_5(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       sum1,sum2,sum3,sum4,sum5,x1,x2,x3,x4,x5,*z = 0,*zarray;
  const PetscScalar *xb,*x;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  const PetscInt    *idx,*ii,*ridx=NULL;
  PetscInt          mbs,i,j,n;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,5*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    sum1        = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0;
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);     /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+25*n,25*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 5*(*idx++);
      x1    = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4];
      sum1 += v[0]*x1 + v[5]*x2 + v[10]*x3 + v[15]*x4 + v[20]*x5;
      sum2 += v[1]*x1 + v[6]*x2 + v[11]*x3 + v[16]*x4 + v[21]*x5;
      sum3 += v[2]*x1 + v[7]*x2 + v[12]*x3 + v[17]*x4 + v[22]*x5;
      sum4 += v[3]*x1 + v[8]*x2 + v[13]*x3 + v[18]*x4 + v[23]*x5;
      sum5 += v[4]*x1 + v[9]*x2 + v[14]*x3 + v[19]*x4 + v[24]*x5;
      v    += 25;
    }
    if (usecprow) z = zarray + 5*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5;
    if (!usecprow) z += 5;
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(50.0*a->nz - 5.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMult_SeqBAIJ_6(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4,sum5,sum6;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2,x3,x4,x5,x6,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,*idx,*ii,j,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,6*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n  = ii[1] - ii[0];
    ii++;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    sum6 = 0.0;

    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);     /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+36*n,36*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 6*(*idx++);
      x1    = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5];
      sum1 += v[0]*x1 + v[6]*x2  + v[12]*x3 + v[18]*x4 + v[24]*x5 + v[30]*x6;
      sum2 += v[1]*x1 + v[7]*x2  + v[13]*x3 + v[19]*x4 + v[25]*x5 + v[31]*x6;
      sum3 += v[2]*x1 + v[8]*x2  + v[14]*x3 + v[20]*x4 + v[26]*x5 + v[32]*x6;
      sum4 += v[3]*x1 + v[9]*x2  + v[15]*x3 + v[21]*x4 + v[27]*x5 + v[33]*x6;
      sum5 += v[4]*x1 + v[10]*x2 + v[16]*x3 + v[22]*x4 + v[28]*x5 + v[34]*x6;
      sum6 += v[5]*x1 + v[11]*x2 + v[17]*x3 + v[23]*x4 + v[29]*x5 + v[35]*x6;
      v    += 36;
    }
    if (usecprow) z = zarray + 6*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6;
    if (!usecprow) z += 6;
  }

  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(72.0*a->nz - 6.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMult_SeqBAIJ_7(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2,x3,x4,x5,x6,x7,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,*idx,*ii,j,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,7*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n  = ii[1] - ii[0];
    ii++;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    sum6 = 0.0;
    sum7 = 0.0;

    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);     /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+49*n,49*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 7*(*idx++);
      x1    = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5]; x7 = xb[6];
      sum1 += v[0]*x1 + v[7]*x2  + v[14]*x3 + v[21]*x4 + v[28]*x5 + v[35]*x6 + v[42]*x7;
      sum2 += v[1]*x1 + v[8]*x2  + v[15]*x3 + v[22]*x4 + v[29]*x5 + v[36]*x6 + v[43]*x7;
      sum3 += v[2]*x1 + v[9]*x2  + v[16]*x3 + v[23]*x4 + v[30]*x5 + v[37]*x6 + v[44]*x7;
      sum4 += v[3]*x1 + v[10]*x2 + v[17]*x3 + v[24]*x4 + v[31]*x5 + v[38]*x6 + v[45]*x7;
      sum5 += v[4]*x1 + v[11]*x2 + v[18]*x3 + v[25]*x4 + v[32]*x5 + v[39]*x6 + v[46]*x7;
      sum6 += v[5]*x1 + v[12]*x2 + v[19]*x3 + v[26]*x4 + v[33]*x5 + v[40]*x6 + v[47]*x7;
      sum7 += v[6]*x1 + v[13]*x2 + v[20]*x3 + v[27]*x4 + v[34]*x5 + v[41]*x6 + v[48]*x7;
      v    += 49;
    }
    if (usecprow) z = zarray + 7*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    if (!usecprow) z += 7;
  }

  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(98.0*a->nz - 7.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_IMMINTRIN_H) && defined(__AVX2__) && defined(__FMA__) && defined(PETSC_USE_REAL_DOUBLE) && !defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_64BIT_INDICES)
PetscErrorCode MatMult_SeqBAIJ_9_AVX2(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ    *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar    *z = 0,*work,*workt,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar   *v;
  PetscErrorCode ierr;
  PetscInt       mbs,i,bs=A->rmap->bs,j,n,bs2=a->bs2;
  const PetscInt *idx,*ii,*ridx=NULL;
  PetscInt       k;
  PetscBool      usecprow=a->compressedrow.use;

  __m256d a0,a1,a2,a3,a4,a5;
  __m256d w0,w1,w2,w3;
  __m256d z0,z1,z2;
  __m256i mask1 = _mm256_set_epi64x(0LL, 0LL, 0LL, 1LL<<63);

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,bs*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  if (!a->mult_work) {
    k    = PetscMax(A->rmap->n,A->cmap->n);
    ierr = PetscMalloc1(k+1,&a->mult_work);CHKERRQ(ierr);
  }

  work = a->mult_work;
  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    workt       = work;
    for (j=0; j<n; j++) {
      xb = x + bs*(*idx++);
      for (k=0; k<bs; k++) workt[k] = xb[k];
      workt += bs;
    }
    if (usecprow) z = zarray + bs*ridx[i];

    z0 = _mm256_setzero_pd(); z1 = _mm256_setzero_pd(); z2 = _mm256_setzero_pd();

    for (j=0; j<n; j++) {
      /* first column of a */
      w0 = _mm256_set1_pd(work[j*9  ]);
      a0 = _mm256_loadu_pd(&v[j*81  ]); z0 = _mm256_fmadd_pd(a0,w0,z0);
      a1 = _mm256_loadu_pd(&v[j*81+4]); z1 = _mm256_fmadd_pd(a1,w0,z1);
      a2 = _mm256_loadu_pd(&v[j*81+8]); z2 = _mm256_fmadd_pd(a2,w0,z2);

      /* second column of a */
      w1 = _mm256_set1_pd(work[j*9+ 1]);
      a0 = _mm256_loadu_pd(&v[j*81+ 9]); z0 = _mm256_fmadd_pd(a0,w1,z0);
      a1 = _mm256_loadu_pd(&v[j*81+13]); z1 = _mm256_fmadd_pd(a1,w1,z1);
      a2 = _mm256_loadu_pd(&v[j*81+17]); z2 = _mm256_fmadd_pd(a2,w1,z2);

      /* third column of a */
      w2 = _mm256_set1_pd(work[j*9 +2]);
      a3 = _mm256_loadu_pd(&v[j*81+18]); z0 = _mm256_fmadd_pd(a3,w2,z0);
      a4 = _mm256_loadu_pd(&v[j*81+22]); z1 = _mm256_fmadd_pd(a4,w2,z1);
      a5 = _mm256_loadu_pd(&v[j*81+26]); z2 = _mm256_fmadd_pd(a5,w2,z2);

      /* fourth column of a */
      w3 = _mm256_set1_pd(work[j*9+ 3]);
      a0 = _mm256_loadu_pd(&v[j*81+27]); z0 = _mm256_fmadd_pd(a0,w3,z0);
      a1 = _mm256_loadu_pd(&v[j*81+31]); z1 = _mm256_fmadd_pd(a1,w3,z1);
      a2 = _mm256_loadu_pd(&v[j*81+35]); z2 = _mm256_fmadd_pd(a2,w3,z2);

      /* fifth column of a */
      w0 = _mm256_set1_pd(work[j*9+ 4]);
      a3 = _mm256_loadu_pd(&v[j*81+36]); z0 = _mm256_fmadd_pd(a3,w0,z0);
      a4 = _mm256_loadu_pd(&v[j*81+40]); z1 = _mm256_fmadd_pd(a4,w0,z1);
      a5 = _mm256_loadu_pd(&v[j*81+44]); z2 = _mm256_fmadd_pd(a5,w0,z2);

      /* sixth column of a */
      w1 = _mm256_set1_pd(work[j*9+ 5]);
      a0 = _mm256_loadu_pd(&v[j*81+45]); z0 = _mm256_fmadd_pd(a0,w1,z0);
      a1 = _mm256_loadu_pd(&v[j*81+49]); z1 = _mm256_fmadd_pd(a1,w1,z1);
      a2 = _mm256_loadu_pd(&v[j*81+53]); z2 = _mm256_fmadd_pd(a2,w1,z2);

      /* seventh column of a */
      w2 = _mm256_set1_pd(work[j*9+ 6]);
      a0 = _mm256_loadu_pd(&v[j*81+54]); z0 = _mm256_fmadd_pd(a0,w2,z0);
      a1 = _mm256_loadu_pd(&v[j*81+58]); z1 = _mm256_fmadd_pd(a1,w2,z1);
      a2 = _mm256_loadu_pd(&v[j*81+62]); z2 = _mm256_fmadd_pd(a2,w2,z2);

      /* eigth column of a */
      w3 = _mm256_set1_pd(work[j*9+ 7]);
      a3 = _mm256_loadu_pd(&v[j*81+63]); z0 = _mm256_fmadd_pd(a3,w3,z0);
      a4 = _mm256_loadu_pd(&v[j*81+67]); z1 = _mm256_fmadd_pd(a4,w3,z1);
      a5 = _mm256_loadu_pd(&v[j*81+71]); z2 = _mm256_fmadd_pd(a5,w3,z2);

      /* ninth column of a */
      w0 = _mm256_set1_pd(work[j*9+ 8]);
      a0 = _mm256_loadu_pd(&v[j*81+72]); z0 = _mm256_fmadd_pd(a0,w0,z0);
      a1 = _mm256_loadu_pd(&v[j*81+76]); z1 = _mm256_fmadd_pd(a1,w0,z1);
      a2 = _mm256_maskload_pd(&v[j*81+80],mask1); z2 = _mm256_fmadd_pd(a2,w0,z2);
    }

    _mm256_storeu_pd(&z[ 0], z0); _mm256_storeu_pd(&z[ 4], z1); _mm256_maskstore_pd(&z[8], mask1, z2);

    v += n*bs2;
    if (!usecprow) z += bs;
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*a->nz*bs2 - bs*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

PetscErrorCode MatMult_SeqBAIJ_11(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11;
  const PetscScalar *x,*xb;
  PetscScalar       *zarray,xv;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  const PetscInt    *ii,*ij=a->j,*idx;
  PetscInt          mbs,i,j,k,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  v = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,11*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n    = ii[i+1] - ii[i];
    idx  = ij + ii[i];
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    sum8 = 0.0; sum9 = 0.0; sum10 = 0.0; sum11 = 0.0;

    for (j=0; j<n; j++) {
      xb = x + 11*(idx[j]);

      for (k=0; k<11; k++) {
        xv     =  xb[k];
        sum1  += v[0]*xv;
        sum2  += v[1]*xv;
        sum3  += v[2]*xv;
        sum4  += v[3]*xv;
        sum5  += v[4]*xv;
        sum6  += v[5]*xv;
        sum7  += v[6]*xv;
        sum8  += v[7]*xv;
        sum9  += v[8]*xv;
        sum10 += v[9]*xv;
        sum11 += v[10]*xv;
        v     += 11;
      }
    }
    if (usecprow) z = zarray + 11*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum8; z[8] = sum9; z[9] = sum10; z[10] = sum11;

    if (!usecprow) z += 11;
  }

  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(242.0*a->nz - 11.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* MatMult_SeqBAIJ_15 version 1: Columns in the block are accessed one at a time */
/* Default MatMult for block size 15 */
PetscErrorCode MatMult_SeqBAIJ_15_ver1(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15;
  const PetscScalar *x,*xb;
  PetscScalar       *zarray,xv;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  const PetscInt    *ii,*ij=a->j,*idx;
  PetscInt          mbs,i,j,k,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  v = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,15*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n    = ii[i+1] - ii[i];
    idx  = ij + ii[i];
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    sum8 = 0.0; sum9 = 0.0; sum10 = 0.0; sum11 = 0.0; sum12 = 0.0; sum13 = 0.0; sum14 = 0.0;sum15 = 0.0;

    for (j=0; j<n; j++) {
      xb = x + 15*(idx[j]);

      for (k=0; k<15; k++) {
        xv     =  xb[k];
        sum1  += v[0]*xv;
        sum2  += v[1]*xv;
        sum3  += v[2]*xv;
        sum4  += v[3]*xv;
        sum5  += v[4]*xv;
        sum6  += v[5]*xv;
        sum7  += v[6]*xv;
        sum8  += v[7]*xv;
        sum9  += v[8]*xv;
        sum10 += v[9]*xv;
        sum11 += v[10]*xv;
        sum12 += v[11]*xv;
        sum13 += v[12]*xv;
        sum14 += v[13]*xv;
        sum15 += v[14]*xv;
        v     += 15;
      }
    }
    if (usecprow) z = zarray + 15*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum8; z[8] = sum9; z[9] = sum10; z[10] = sum11; z[11] = sum12; z[12] = sum13; z[13] = sum14;z[14] = sum15;

    if (!usecprow) z += 15;
  }

  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(450.0*a->nz - 15.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* MatMult_SeqBAIJ_15_ver2 : Columns in the block are accessed in sets of 4,4,4,3 */
PetscErrorCode MatMult_SeqBAIJ_15_ver2(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2,x3,x4,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  const PetscInt    *ii,*ij=a->j,*idx;
  PetscInt          mbs,i,j,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  v = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,15*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n    = ii[i+1] - ii[i];
    idx  = ij + ii[i];
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    sum8 = 0.0; sum9 = 0.0; sum10 = 0.0; sum11 = 0.0; sum12 = 0.0; sum13 = 0.0; sum14 = 0.0;sum15 = 0.0;

    for (j=0; j<n; j++) {
      xb = x + 15*(idx[j]);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];

      sum1  += v[0]*x1 + v[15]*x2 + v[30]*x3   + v[45]*x4;
      sum2  += v[1]*x1 + v[16]*x2 + v[31]*x3   + v[46]*x4;
      sum3  += v[2]*x1 + v[17]*x2 + v[32]*x3  + v[47]*x4;
      sum4  += v[3]*x1 + v[18]*x2 + v[33]*x3  + v[48]*x4;
      sum5  += v[4]*x1 + v[19]*x2 + v[34]*x3   + v[49]*x4;
      sum6  += v[5]*x1 + v[20]*x2 + v[35]*x3   + v[50]*x4;
      sum7  += v[6]*x1 + v[21]*x2 + v[36]*x3  + v[51]*x4;
      sum8  += v[7]*x1 + v[22]*x2 + v[37]*x3  + v[52]*x4;
      sum9  += v[8]*x1 + v[23]*x2 + v[38]*x3   + v[53]*x4;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3   + v[54]*x4;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3  + v[55]*x4;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3  + v[56]*x4;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3   + v[57]*x4;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3   + v[58]*x4;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3  + v[59]*x4;

      v += 60;

      x1 = xb[4]; x2 = xb[5]; x3 = xb[6]; x4 = xb[7];

      sum1  += v[0]*x1 + v[15]*x2 + v[30]*x3   + v[45]*x4;
      sum2  += v[1]*x1 + v[16]*x2 + v[31]*x3   + v[46]*x4;
      sum3  += v[2]*x1 + v[17]*x2 + v[32]*x3  + v[47]*x4;
      sum4  += v[3]*x1 + v[18]*x2 + v[33]*x3  + v[48]*x4;
      sum5  += v[4]*x1 + v[19]*x2 + v[34]*x3   + v[49]*x4;
      sum6  += v[5]*x1 + v[20]*x2 + v[35]*x3   + v[50]*x4;
      sum7  += v[6]*x1 + v[21]*x2 + v[36]*x3  + v[51]*x4;
      sum8  += v[7]*x1 + v[22]*x2 + v[37]*x3  + v[52]*x4;
      sum9  += v[8]*x1 + v[23]*x2 + v[38]*x3   + v[53]*x4;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3   + v[54]*x4;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3  + v[55]*x4;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3  + v[56]*x4;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3   + v[57]*x4;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3   + v[58]*x4;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3  + v[59]*x4;
      v     += 60;

      x1     = xb[8]; x2 = xb[9]; x3 = xb[10]; x4 = xb[11];
      sum1  += v[0]*x1 + v[15]*x2 + v[30]*x3   + v[45]*x4;
      sum2  += v[1]*x1 + v[16]*x2 + v[31]*x3   + v[46]*x4;
      sum3  += v[2]*x1 + v[17]*x2 + v[32]*x3  + v[47]*x4;
      sum4  += v[3]*x1 + v[18]*x2 + v[33]*x3  + v[48]*x4;
      sum5  += v[4]*x1 + v[19]*x2 + v[34]*x3   + v[49]*x4;
      sum6  += v[5]*x1 + v[20]*x2 + v[35]*x3   + v[50]*x4;
      sum7  += v[6]*x1 + v[21]*x2 + v[36]*x3  + v[51]*x4;
      sum8  += v[7]*x1 + v[22]*x2 + v[37]*x3  + v[52]*x4;
      sum9  += v[8]*x1 + v[23]*x2 + v[38]*x3   + v[53]*x4;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3   + v[54]*x4;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3  + v[55]*x4;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3  + v[56]*x4;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3   + v[57]*x4;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3   + v[58]*x4;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3  + v[59]*x4;
      v     += 60;

      x1     = xb[12]; x2 = xb[13]; x3 = xb[14];
      sum1  += v[0]*x1 + v[15]*x2 + v[30]*x3;
      sum2  += v[1]*x1 + v[16]*x2 + v[31]*x3;
      sum3  += v[2]*x1 + v[17]*x2 + v[32]*x3;
      sum4  += v[3]*x1 + v[18]*x2 + v[33]*x3;
      sum5  += v[4]*x1 + v[19]*x2 + v[34]*x3;
      sum6  += v[5]*x1 + v[20]*x2 + v[35]*x3;
      sum7  += v[6]*x1 + v[21]*x2 + v[36]*x3;
      sum8  += v[7]*x1 + v[22]*x2 + v[37]*x3;
      sum9  += v[8]*x1 + v[23]*x2 + v[38]*x3;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3;
      v     += 45;
    }
    if (usecprow) z = zarray + 15*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum8; z[8] = sum9; z[9] = sum10; z[10] = sum11; z[11] = sum12; z[12] = sum13; z[13] = sum14;z[14] = sum15;

    if (!usecprow) z += 15;
  }

  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(450.0*a->nz - 15.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* MatMult_SeqBAIJ_15_ver3 : Columns in the block are accessed in sets of 8,7 */
PetscErrorCode MatMult_SeqBAIJ_15_ver3(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2,x3,x4,x5,x6,x7,x8,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  const PetscInt    *ii,*ij=a->j,*idx;
  PetscInt          mbs,i,j,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  v = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,15*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n    = ii[i+1] - ii[i];
    idx  = ij + ii[i];
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    sum8 = 0.0; sum9 = 0.0; sum10 = 0.0; sum11 = 0.0; sum12 = 0.0; sum13 = 0.0; sum14 = 0.0;sum15 = 0.0;

    for (j=0; j<n; j++) {
      xb = x + 15*(idx[j]);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5]; x7 = xb[6];
      x8 = xb[7];

      sum1  += v[0]*x1 + v[15]*x2  + v[30]*x3  + v[45]*x4 + v[60]*x5 + v[75]*x6 + v[90]*x7 + v[105]*x8;
      sum2  += v[1]*x1 + v[16]*x2  + v[31]*x3  + v[46]*x4 + v[61]*x5 + v[76]*x6 + v[91]*x7 + v[106]*x8;
      sum3  += v[2]*x1 + v[17]*x2  + v[32]*x3  + v[47]*x4 + v[62]*x5 + v[77]*x6 + v[92]*x7 + v[107]*x8;
      sum4  += v[3]*x1 + v[18]*x2 + v[33]*x3  + v[48]*x4 + v[63]*x5 + v[78]*x6 + v[93]*x7 + v[108]*x8;
      sum5  += v[4]*x1 + v[19]*x2 + v[34]*x3  + v[49]*x4 + v[64]*x5 + v[79]*x6 + v[94]*x7 + v[109]*x8;
      sum6  += v[5]*x1 + v[20]*x2 + v[35]*x3  + v[50]*x4 + v[65]*x5 + v[80]*x6 + v[95]*x7 + v[110]*x8;
      sum7  += v[6]*x1 + v[21]*x2 + v[36]*x3  + v[51]*x4 + v[66]*x5 + v[81]*x6 + v[96]*x7 + v[111]*x8;
      sum8  += v[7]*x1 + v[22]*x2  + v[37]*x3  + v[52]*x4 + v[67]*x5 + v[82]*x6 + v[97]*x7 + v[112]*x8;
      sum9  += v[8]*x1 + v[23]*x2  + v[38]*x3  + v[53]*x4 + v[68]*x5 + v[83]*x6 + v[98]*x7 + v[113]*x8;
      sum10 += v[9]*x1 + v[24]*x2  + v[39]*x3  + v[54]*x4 + v[69]*x5 + v[84]*x6 + v[99]*x7 + v[114]*x8;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3  + v[55]*x4 + v[70]*x5 + v[85]*x6 + v[100]*x7 + v[115]*x8;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3  + v[56]*x4 + v[71]*x5 + v[86]*x6 + v[101]*x7 + v[116]*x8;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3  + v[57]*x4 + v[72]*x5 + v[87]*x6 + v[102]*x7 + v[117]*x8;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3  + v[58]*x4 + v[73]*x5 + v[88]*x6 + v[103]*x7 + v[118]*x8;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3  + v[59]*x4 + v[74]*x5 + v[89]*x6 + v[104]*x7 + v[119]*x8;
      v     += 120;

      x1 = xb[8]; x2 = xb[9]; x3 = xb[10]; x4 = xb[11]; x5 = xb[12]; x6 = xb[13]; x7 = xb[14];

      sum1  += v[0]*x1 + v[15]*x2  + v[30]*x3  + v[45]*x4 + v[60]*x5 + v[75]*x6 + v[90]*x7;
      sum2  += v[1]*x1 + v[16]*x2  + v[31]*x3  + v[46]*x4 + v[61]*x5 + v[76]*x6 + v[91]*x7;
      sum3  += v[2]*x1 + v[17]*x2  + v[32]*x3  + v[47]*x4 + v[62]*x5 + v[77]*x6 + v[92]*x7;
      sum4  += v[3]*x1 + v[18]*x2 + v[33]*x3  + v[48]*x4 + v[63]*x5 + v[78]*x6 + v[93]*x7;
      sum5  += v[4]*x1 + v[19]*x2 + v[34]*x3  + v[49]*x4 + v[64]*x5 + v[79]*x6 + v[94]*x7;
      sum6  += v[5]*x1 + v[20]*x2 + v[35]*x3  + v[50]*x4 + v[65]*x5 + v[80]*x6 + v[95]*x7;
      sum7  += v[6]*x1 + v[21]*x2 + v[36]*x3  + v[51]*x4 + v[66]*x5 + v[81]*x6 + v[96]*x7;
      sum8  += v[7]*x1 + v[22]*x2  + v[37]*x3  + v[52]*x4 + v[67]*x5 + v[82]*x6 + v[97]*x7;
      sum9  += v[8]*x1 + v[23]*x2  + v[38]*x3  + v[53]*x4 + v[68]*x5 + v[83]*x6 + v[98]*x7;
      sum10 += v[9]*x1 + v[24]*x2  + v[39]*x3  + v[54]*x4 + v[69]*x5 + v[84]*x6 + v[99]*x7;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3  + v[55]*x4 + v[70]*x5 + v[85]*x6 + v[100]*x7;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3  + v[56]*x4 + v[71]*x5 + v[86]*x6 + v[101]*x7;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3  + v[57]*x4 + v[72]*x5 + v[87]*x6 + v[102]*x7;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3  + v[58]*x4 + v[73]*x5 + v[88]*x6 + v[103]*x7;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3  + v[59]*x4 + v[74]*x5 + v[89]*x6 + v[104]*x7;
      v     += 105;
    }
    if (usecprow) z = zarray + 15*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum8; z[8] = sum9; z[9] = sum10; z[10] = sum11; z[11] = sum12; z[12] = sum13; z[13] = sum14;z[14] = sum15;

    if (!usecprow) z += 15;
  }

  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(450.0*a->nz - 15.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* MatMult_SeqBAIJ_15_ver4 : All columns in the block are accessed at once */
PetscErrorCode MatMult_SeqBAIJ_15_ver4(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  const PetscInt    *ii,*ij=a->j,*idx;
  PetscInt          mbs,i,j,n,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  v = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,15*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  for (i=0; i<mbs; i++) {
    n    = ii[i+1] - ii[i];
    idx  = ij + ii[i];
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    sum8 = 0.0; sum9 = 0.0; sum10 = 0.0; sum11 = 0.0; sum12 = 0.0; sum13 = 0.0; sum14 = 0.0;sum15 = 0.0;

    for (j=0; j<n; j++) {
      xb = x + 15*(idx[j]);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5]; x7 = xb[6];
      x8 = xb[7]; x9 = xb[8]; x10 = xb[9]; x11 = xb[10]; x12 = xb[11]; x13 = xb[12]; x14 = xb[13];x15 = xb[14];

      sum1  +=  v[0]*x1  + v[15]*x2 + v[30]*x3 + v[45]*x4 + v[60]*x5 + v[75]*x6 + v[90]*x7  + v[105]*x8 + v[120]*x9 + v[135]*x10 + v[150]*x11 + v[165]*x12 + v[180]*x13 + v[195]*x14 + v[210]*x15;
      sum2  +=  v[1]*x1  + v[16]*x2 + v[31]*x3 + v[46]*x4 + v[61]*x5 + v[76]*x6 + v[91]*x7  + v[106]*x8 + v[121]*x9 + v[136]*x10 + v[151]*x11 + v[166]*x12 + v[181]*x13 + v[196]*x14 + v[211]*x15;
      sum3  +=  v[2]*x1  + v[17]*x2 + v[32]*x3 + v[47]*x4 + v[62]*x5 + v[77]*x6 + v[92]*x7  + v[107]*x8 + v[122]*x9 + v[137]*x10 + v[152]*x11 + v[167]*x12 + v[182]*x13 + v[197]*x14 + v[212]*x15;
      sum4  +=  v[3]*x1  + v[18]*x2 + v[33]*x3 + v[48]*x4 + v[63]*x5 + v[78]*x6 + v[93]*x7  + v[108]*x8 + v[123]*x9 + v[138]*x10 + v[153]*x11 + v[168]*x12 + v[183]*x13 + v[198]*x14 + v[213]*x15;
      sum5  += v[4]*x1  + v[19]*x2 + v[34]*x3 + v[49]*x4 + v[64]*x5 + v[79]*x6 + v[94]*x7  + v[109]*x8 + v[124]*x9 + v[139]*x10 + v[154]*x11 + v[169]*x12 + v[184]*x13 + v[199]*x14 + v[214]*x15;
      sum6  += v[5]*x1  + v[20]*x2 + v[35]*x3 + v[50]*x4 + v[65]*x5 + v[80]*x6 + v[95]*x7  + v[110]*x8 + v[125]*x9 + v[140]*x10 + v[155]*x11 + v[170]*x12 + v[185]*x13 + v[200]*x14 + v[215]*x15;
      sum7  += v[6]*x1  + v[21]*x2 + v[36]*x3 + v[51]*x4 + v[66]*x5 + v[81]*x6 + v[96]*x7  + v[111]*x8 + v[126]*x9 + v[141]*x10 + v[156]*x11 + v[171]*x12 + v[186]*x13 + v[201]*x14 + v[216]*x15;
      sum8  += v[7]*x1  + v[22]*x2 + v[37]*x3 + v[52]*x4 + v[67]*x5 + v[82]*x6 + v[97]*x7  + v[112]*x8 + v[127]*x9 + v[142]*x10 + v[157]*x11 + v[172]*x12 + v[187]*x13 + v[202]*x14 + v[217]*x15;
      sum9  += v[8]*x1  + v[23]*x2 + v[38]*x3 + v[53]*x4 + v[68]*x5 + v[83]*x6 + v[98]*x7  + v[113]*x8 + v[128]*x9 + v[143]*x10 + v[158]*x11 + v[173]*x12 + v[188]*x13 + v[203]*x14 + v[218]*x15;
      sum10 += v[9]*x1  + v[24]*x2 + v[39]*x3 + v[54]*x4 + v[69]*x5 + v[84]*x6 + v[99]*x7  + v[114]*x8 + v[129]*x9 + v[144]*x10 + v[159]*x11 + v[174]*x12 + v[189]*x13 + v[204]*x14 + v[219]*x15;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3 + v[55]*x4 + v[70]*x5 + v[85]*x6 + v[100]*x7 + v[115]*x8 + v[130]*x9 + v[145]*x10 + v[160]*x11 + v[175]*x12 + v[190]*x13 + v[205]*x14 + v[220]*x15;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3 + v[56]*x4 + v[71]*x5 + v[86]*x6 + v[101]*x7 + v[116]*x8 + v[131]*x9 + v[146]*x10 + v[161]*x11 + v[176]*x12 + v[191]*x13 + v[206]*x14 + v[221]*x15;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3 + v[57]*x4 + v[72]*x5 + v[87]*x6 + v[102]*x7 + v[117]*x8 + v[132]*x9 + v[147]*x10 + v[162]*x11 + v[177]*x12 + v[192]*x13 + v[207]*x14 + v[222]*x15;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3 + v[58]*x4 + v[73]*x5 + v[88]*x6 + v[103]*x7 + v[118]*x8 + v[133]*x9 + v[148]*x10 + v[163]*x11 + v[178]*x12 + v[193]*x13 + v[208]*x14 + v[223]*x15;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3 + v[59]*x4 + v[74]*x5 + v[89]*x6 + v[104]*x7 + v[119]*x8 + v[134]*x9 + v[149]*x10 + v[164]*x11 + v[179]*x12 + v[194]*x13 + v[209]*x14 + v[224]*x15;
      v     += 225;
    }
    if (usecprow) z = zarray + 15*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum8; z[8] = sum9; z[9] = sum10; z[10] = sum11; z[11] = sum12; z[12] = sum13; z[13] = sum14;z[14] = sum15;

    if (!usecprow) z += 15;
  }

  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(450.0*a->nz - 15.0*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    This will not work with MatScalar == float because it calls the BLAS
*/
PetscErrorCode MatMult_SeqBAIJ_N(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,*work,*workt,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,bs=A->rmap->bs,j,n,bs2=a->bs2;
  const PetscInt    *idx,*ii,*ridx=NULL;
  PetscInt          ncols,k;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscArrayzero(zarray,bs*a->mbs);CHKERRQ(ierr);
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  if (!a->mult_work) {
    k    = PetscMax(A->rmap->n,A->cmap->n);
    ierr = PetscMalloc1(k+1,&a->mult_work);CHKERRQ(ierr);
  }
  work = a->mult_work;
  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    ncols       = n*bs;
    workt       = work;
    for (j=0; j<n; j++) {
      xb = x + bs*(*idx++);
      for (k=0; k<bs; k++) workt[k] = xb[k];
      workt += bs;
    }
    if (usecprow) z = zarray + bs*ridx[i];
    PetscKernel_w_gets_Ar_times_v(bs,ncols,work,v,z);
    /* BLASgemv_("N",&bs,&ncols,&_DOne,v,&bs,work,&_One,&_DZero,z,&_One); */
    v += n*bs2;
    if (!usecprow) z += bs;
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*a->nz*bs2 - bs*a->nonzerorowcnt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultAdd_SeqBAIJ_1(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  const PetscScalar *x;
  PetscScalar       *y,*z,sum;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs=a->mbs,i,n,*ridx=NULL;
  const PetscInt    *idx,*ii;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayPair(yy,zz,&y,&z);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    if (zz != yy) {
      ierr = PetscArraycpy(z,y,mbs);CHKERRQ(ierr);
    }
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0];
    ii++;
    if (!usecprow) {
      sum         = y[i];
    } else {
      sum = y[ridx[i]];
    }
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA); /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Entries for the next row */
    PetscSparseDensePlusDot(sum,x,v,idx,n);
    v   += n;
    idx += n;
    if (usecprow) {
      z[ridx[i]] = sum;
    } else {
      z[i] = sum;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayPair(yy,zz,&y,&z);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*a->nz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultAdd_SeqBAIJ_2(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *y = 0,*z = 0,sum1,sum2;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2,*yarray,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs = a->mbs,i,n,j;
  const PetscInt    *idx,*ii,*ridx = NULL;
  PetscBool         usecprow = a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    if (zz != yy) {
      ierr = PetscArraycpy(zarray,yarray,2*mbs);CHKERRQ(ierr);
    }
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y  = yarray;
    z  = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow) {
      z = zarray + 2*ridx[i];
      y = yarray + 2*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1];
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+4*n,4*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb = x + 2*(*idx++);
      x1 = xb[0];
      x2 = xb[1];

      sum1 += v[0]*x1 + v[2]*x2;
      sum2 += v[1]*x1 + v[3]*x2;
      v    += 4;
    }
    z[0] = sum1; z[1] = sum2;
    if (!usecprow) {
      z += 2; y += 2;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(4.0*a->nz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultAdd_SeqBAIJ_3(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *y = 0,*z = 0,sum1,sum2,sum3,x1,x2,x3,*yarray,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs = a->mbs,i,j,n;
  const PetscInt    *idx,*ii,*ridx = NULL;
  PetscBool         usecprow = a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    if (zz != yy) {
      ierr = PetscArraycpy(zarray,yarray,3*mbs);CHKERRQ(ierr);
    }
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y  = yarray;
    z  = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow) {
      z = zarray + 3*ridx[i];
      y = yarray + 3*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2];
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+9*n,9*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 3*(*idx++); x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      sum1 += v[0]*x1 + v[3]*x2 + v[6]*x3;
      sum2 += v[1]*x1 + v[4]*x2 + v[7]*x3;
      sum3 += v[2]*x1 + v[5]*x2 + v[8]*x3;
      v    += 9;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3;
    if (!usecprow) {
      z += 3; y += 3;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(18.0*a->nz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultAdd_SeqBAIJ_4(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *y = 0,*z = 0,sum1,sum2,sum3,sum4,x1,x2,x3,x4,*yarray,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs = a->mbs,i,j,n;
  const PetscInt    *idx,*ii,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    if (zz != yy) {
      ierr = PetscArraycpy(zarray,yarray,4*mbs);CHKERRQ(ierr);
    }
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y  = yarray;
    z  = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow) {
      z = zarray + 4*ridx[i];
      y = yarray + 4*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2]; sum4 = y[3];
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);     /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+16*n,16*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 4*(*idx++);
      x1    = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];
      sum1 += v[0]*x1 + v[4]*x2 + v[8]*x3   + v[12]*x4;
      sum2 += v[1]*x1 + v[5]*x2 + v[9]*x3   + v[13]*x4;
      sum3 += v[2]*x1 + v[6]*x2 + v[10]*x3  + v[14]*x4;
      sum4 += v[3]*x1 + v[7]*x2 + v[11]*x3  + v[15]*x4;
      v    += 16;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4;
    if (!usecprow) {
      z += 4; y += 4;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(32.0*a->nz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultAdd_SeqBAIJ_5(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *y = 0,*z = 0,sum1,sum2,sum3,sum4,sum5,x1,x2,x3,x4,x5;
  const PetscScalar *x,*xb;
  PetscScalar       *yarray,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs = a->mbs,i,j,n;
  const PetscInt    *idx,*ii,*ridx = NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    if (zz != yy) {
      ierr = PetscArraycpy(zarray,yarray,5*mbs);CHKERRQ(ierr);
    }
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y  = yarray;
    z  = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow) {
      z = zarray + 5*ridx[i];
      y = yarray + 5*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2]; sum4 = y[3]; sum5 = y[4];
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);     /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+25*n,25*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 5*(*idx++);
      x1    = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4];
      sum1 += v[0]*x1 + v[5]*x2 + v[10]*x3  + v[15]*x4 + v[20]*x5;
      sum2 += v[1]*x1 + v[6]*x2 + v[11]*x3  + v[16]*x4 + v[21]*x5;
      sum3 += v[2]*x1 + v[7]*x2 + v[12]*x3  + v[17]*x4 + v[22]*x5;
      sum4 += v[3]*x1 + v[8]*x2 + v[13]*x3  + v[18]*x4 + v[23]*x5;
      sum5 += v[4]*x1 + v[9]*x2 + v[14]*x3  + v[19]*x4 + v[24]*x5;
      v    += 25;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5;
    if (!usecprow) {
      z += 5; y += 5;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(50.0*a->nz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultAdd_SeqBAIJ_6(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *y = 0,*z = 0,sum1,sum2,sum3,sum4,sum5,sum6;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2,x3,x4,x5,x6,*yarray,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs = a->mbs,i,j,n;
  const PetscInt    *idx,*ii,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    if (zz != yy) {
      ierr = PetscArraycpy(zarray,yarray,6*mbs);CHKERRQ(ierr);
    }
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y  = yarray;
    z  = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow) {
      z = zarray + 6*ridx[i];
      y = yarray + 6*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2]; sum4 = y[3]; sum5 = y[4]; sum6 = y[5];
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);     /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+36*n,36*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 6*(*idx++);
      x1    = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5];
      sum1 += v[0]*x1 + v[6]*x2  + v[12]*x3  + v[18]*x4 + v[24]*x5 + v[30]*x6;
      sum2 += v[1]*x1 + v[7]*x2  + v[13]*x3  + v[19]*x4 + v[25]*x5 + v[31]*x6;
      sum3 += v[2]*x1 + v[8]*x2  + v[14]*x3  + v[20]*x4 + v[26]*x5 + v[32]*x6;
      sum4 += v[3]*x1 + v[9]*x2  + v[15]*x3  + v[21]*x4 + v[27]*x5 + v[33]*x6;
      sum5 += v[4]*x1 + v[10]*x2 + v[16]*x3  + v[22]*x4 + v[28]*x5 + v[34]*x6;
      sum6 += v[5]*x1 + v[11]*x2 + v[17]*x3  + v[23]*x4 + v[29]*x5 + v[35]*x6;
      v    += 36;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6;
    if (!usecprow) {
      z += 6; y += 6;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(72.0*a->nz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultAdd_SeqBAIJ_7(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *y = 0,*z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2,x3,x4,x5,x6,x7,*yarray,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs = a->mbs,i,j,n;
  const PetscInt    *idx,*ii,*ridx = NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    if (zz != yy) {
      ierr = PetscArraycpy(zarray,yarray,7*mbs);CHKERRQ(ierr);
    }
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y  = yarray;
    z  = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow) {
      z = zarray + 7*ridx[i];
      y = yarray + 7*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2]; sum4 = y[3]; sum5 = y[4]; sum6 = y[5]; sum7 = y[6];
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);     /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+49*n,49*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 7*(*idx++);
      x1    = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5]; x7 = xb[6];
      sum1 += v[0]*x1 + v[7]*x2  + v[14]*x3  + v[21]*x4 + v[28]*x5 + v[35]*x6 + v[42]*x7;
      sum2 += v[1]*x1 + v[8]*x2  + v[15]*x3  + v[22]*x4 + v[29]*x5 + v[36]*x6 + v[43]*x7;
      sum3 += v[2]*x1 + v[9]*x2  + v[16]*x3  + v[23]*x4 + v[30]*x5 + v[37]*x6 + v[44]*x7;
      sum4 += v[3]*x1 + v[10]*x2 + v[17]*x3  + v[24]*x4 + v[31]*x5 + v[38]*x6 + v[45]*x7;
      sum5 += v[4]*x1 + v[11]*x2 + v[18]*x3  + v[25]*x4 + v[32]*x5 + v[39]*x6 + v[46]*x7;
      sum6 += v[5]*x1 + v[12]*x2 + v[19]*x3  + v[26]*x4 + v[33]*x5 + v[40]*x6 + v[47]*x7;
      sum7 += v[6]*x1 + v[13]*x2 + v[20]*x3  + v[27]*x4 + v[34]*x5 + v[41]*x6 + v[48]*x7;
      v    += 49;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    if (!usecprow) {
      z += 7; y += 7;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(98.0*a->nz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_IMMINTRIN_H) && defined(__AVX2__) && defined(__FMA__) && defined(PETSC_USE_REAL_DOUBLE) && !defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_64BIT_INDICES)
PetscErrorCode MatMultAdd_SeqBAIJ_9_AVX2(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ    *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar    *z = 0,*work,*workt,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar   *v;
  PetscErrorCode ierr;
  PetscInt       mbs,i,bs=A->rmap->bs,j,n,bs2=a->bs2;
  PetscInt       k;
  PetscBool      usecprow=a->compressedrow.use;
  const PetscInt *idx,*ii,*ridx=NULL;

  __m256d a0,a1,a2,a3,a4,a5;
  __m256d w0,w1,w2,w3;
  __m256d z0,z1,z2;
  __m256i mask1 = _mm256_set_epi64x(0LL, 0LL, 0LL, 1LL<<63);

  PetscFunctionBegin;
  ierr = VecCopy(yy,zz);CHKERRQ(ierr);
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  if (!a->mult_work) {
    k    = PetscMax(A->rmap->n,A->cmap->n);
    ierr = PetscMalloc1(k+1,&a->mult_work);CHKERRQ(ierr);
  }

  work = a->mult_work;
  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    workt       = work;
    for (j=0; j<n; j++) {
      xb = x + bs*(*idx++);
      for (k=0; k<bs; k++) workt[k] = xb[k];
      workt += bs;
    }
    if (usecprow) z = zarray + bs*ridx[i];

    z0 = _mm256_loadu_pd(&z[ 0]); z1 = _mm256_loadu_pd(&z[ 4]); z2 = _mm256_set1_pd(z[ 8]);

    for (j=0; j<n; j++) {
      /* first column of a */
      w0 = _mm256_set1_pd(work[j*9  ]);
      a0 = _mm256_loadu_pd(&v[j*81  ]); z0 = _mm256_fmadd_pd(a0,w0,z0);
      a1 = _mm256_loadu_pd(&v[j*81+4]); z1 = _mm256_fmadd_pd(a1,w0,z1);
      a2 = _mm256_loadu_pd(&v[j*81+8]); z2 = _mm256_fmadd_pd(a2,w0,z2);

      /* second column of a */
      w1 = _mm256_set1_pd(work[j*9+ 1]);
      a0 = _mm256_loadu_pd(&v[j*81+ 9]); z0 = _mm256_fmadd_pd(a0,w1,z0);
      a1 = _mm256_loadu_pd(&v[j*81+13]); z1 = _mm256_fmadd_pd(a1,w1,z1);
      a2 = _mm256_loadu_pd(&v[j*81+17]); z2 = _mm256_fmadd_pd(a2,w1,z2);

      /* third column of a */
      w2 = _mm256_set1_pd(work[j*9+ 2]);
      a3 = _mm256_loadu_pd(&v[j*81+18]); z0 = _mm256_fmadd_pd(a3,w2,z0);
      a4 = _mm256_loadu_pd(&v[j*81+22]); z1 = _mm256_fmadd_pd(a4,w2,z1);
      a5 = _mm256_loadu_pd(&v[j*81+26]); z2 = _mm256_fmadd_pd(a5,w2,z2);

      /* fourth column of a */
      w3 = _mm256_set1_pd(work[j*9+ 3]);
      a0 = _mm256_loadu_pd(&v[j*81+27]); z0 = _mm256_fmadd_pd(a0,w3,z0);
      a1 = _mm256_loadu_pd(&v[j*81+31]); z1 = _mm256_fmadd_pd(a1,w3,z1);
      a2 = _mm256_loadu_pd(&v[j*81+35]); z2 = _mm256_fmadd_pd(a2,w3,z2);

      /* fifth column of a */
      w0 = _mm256_set1_pd(work[j*9+ 4]);
      a3 = _mm256_loadu_pd(&v[j*81+36]); z0 = _mm256_fmadd_pd(a3,w0,z0);
      a4 = _mm256_loadu_pd(&v[j*81+40]); z1 = _mm256_fmadd_pd(a4,w0,z1);
      a5 = _mm256_loadu_pd(&v[j*81+44]); z2 = _mm256_fmadd_pd(a5,w0,z2);

      /* sixth column of a */
      w1 = _mm256_set1_pd(work[j*9+ 5]);
      a0 = _mm256_loadu_pd(&v[j*81+45]); z0 = _mm256_fmadd_pd(a0,w1,z0);
      a1 = _mm256_loadu_pd(&v[j*81+49]); z1 = _mm256_fmadd_pd(a1,w1,z1);
      a2 = _mm256_loadu_pd(&v[j*81+53]); z2 = _mm256_fmadd_pd(a2,w1,z2);

      /* seventh column of a */
      w2 = _mm256_set1_pd(work[j*9+ 6]);
      a0 = _mm256_loadu_pd(&v[j*81+54]); z0 = _mm256_fmadd_pd(a0,w2,z0);
      a1 = _mm256_loadu_pd(&v[j*81+58]); z1 = _mm256_fmadd_pd(a1,w2,z1);
      a2 = _mm256_loadu_pd(&v[j*81+62]); z2 = _mm256_fmadd_pd(a2,w2,z2);

      /* eigth column of a */
      w3 = _mm256_set1_pd(work[j*9+ 7]);
      a3 = _mm256_loadu_pd(&v[j*81+63]); z0 = _mm256_fmadd_pd(a3,w3,z0);
      a4 = _mm256_loadu_pd(&v[j*81+67]); z1 = _mm256_fmadd_pd(a4,w3,z1);
      a5 = _mm256_loadu_pd(&v[j*81+71]); z2 = _mm256_fmadd_pd(a5,w3,z2);

      /* ninth column of a */
      w0 = _mm256_set1_pd(work[j*9+ 8]);
      a0 = _mm256_loadu_pd(&v[j*81+72]); z0 = _mm256_fmadd_pd(a0,w0,z0);
      a1 = _mm256_loadu_pd(&v[j*81+76]); z1 = _mm256_fmadd_pd(a1,w0,z1);
      a2 = _mm256_maskload_pd(&v[j*81+80],mask1); z2 = _mm256_fmadd_pd(a2,w0,z2);
    }

    _mm256_storeu_pd(&z[ 0], z0); _mm256_storeu_pd(&z[ 4], z1); _mm256_maskstore_pd(&z[8], mask1, z2);

    v += n*bs2;
    if (!usecprow) z += bs;
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(162.0*a->nz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

PetscErrorCode MatMultAdd_SeqBAIJ_11(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *y = 0,*z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11;
  const PetscScalar *x,*xb;
  PetscScalar       x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,*yarray,*zarray;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs = a->mbs,i,j,n;
  const PetscInt    *idx,*ii,*ridx = NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    if (zz != yy) {
      ierr = PetscArraycpy(zarray,yarray,7*mbs);CHKERRQ(ierr);
    }
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y  = yarray;
    z  = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow) {
      z = zarray + 11*ridx[i];
      y = yarray + 11*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2]; sum4 = y[3]; sum5 = y[4]; sum6 = y[5]; sum7 = y[6];
    sum8 = y[7]; sum9 = y[8]; sum10 = y[9]; sum11 = y[10];
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);     /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+121*n,121*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    for (j=0; j<n; j++) {
      xb    = x + 11*(*idx++);
      x1    = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5]; x7 = xb[6];x8 = xb[7]; x9 = xb[8]; x10 = xb[9]; x11 = xb[10];
      sum1 += v[  0]*x1 + v[  11]*x2  + v[2*11]*x3  + v[3*11]*x4 + v[4*11]*x5 + v[5*11]*x6 + v[6*11]*x7+ v[7*11]*x8 + v[8*11]*x9  + v[9*11]*x10  + v[10*11]*x11;
      sum2 += v[1+0]*x1 + v[1+11]*x2  + v[1+2*11]*x3  + v[1+3*11]*x4 + v[1+4*11]*x5 + v[1+5*11]*x6 + v[1+6*11]*x7+ v[1+7*11]*x8 + v[1+8*11]*x9  + v[1+9*11]*x10  + v[1+10*11]*x11;
      sum3 += v[2+0]*x1 + v[2+11]*x2  + v[2+2*11]*x3  + v[2+3*11]*x4 + v[2+4*11]*x5 + v[2+5*11]*x6 + v[2+6*11]*x7+ v[2+7*11]*x8 + v[2+8*11]*x9  + v[2+9*11]*x10  + v[2+10*11]*x11;
      sum4 += v[3+0]*x1 + v[3+11]*x2  + v[3+2*11]*x3  + v[3+3*11]*x4 + v[3+4*11]*x5 + v[3+5*11]*x6 + v[3+6*11]*x7+ v[3+7*11]*x8 + v[3+8*11]*x9  + v[3+9*11]*x10  + v[3+10*11]*x11;
      sum5 += v[4+0]*x1 + v[4+11]*x2  + v[4+2*11]*x3  + v[4+3*11]*x4 + v[4+4*11]*x5 + v[4+5*11]*x6 + v[4+6*11]*x7+ v[4+7*11]*x8 + v[4+8*11]*x9  + v[4+9*11]*x10  + v[4+10*11]*x11;
      sum6 += v[5+0]*x1 + v[5+11]*x2  + v[5+2*11]*x3  + v[5+3*11]*x4 + v[5+4*11]*x5 + v[5+5*11]*x6 + v[5+6*11]*x7+ v[5+7*11]*x8 + v[5+8*11]*x9  + v[5+9*11]*x10  + v[5+10*11]*x11;
      sum7 += v[6+0]*x1 + v[6+11]*x2  + v[6+2*11]*x3  + v[6+3*11]*x4 + v[6+4*11]*x5 + v[6+5*11]*x6 + v[6+6*11]*x7+ v[6+7*11]*x8 + v[6+8*11]*x9  + v[6+9*11]*x10  + v[6+10*11]*x11;
      sum8 += v[7+0]*x1 + v[7+11]*x2  + v[7+2*11]*x3  + v[7+3*11]*x4 + v[7+4*11]*x5 + v[7+5*11]*x6 + v[7+6*11]*x7+ v[7+7*11]*x8 + v[7+8*11]*x9  + v[7+9*11]*x10  + v[7+10*11]*x11;
      sum9 += v[8+0]*x1 + v[8+11]*x2  + v[8+2*11]*x3  + v[8+3*11]*x4 + v[8+4*11]*x5 + v[8+5*11]*x6 + v[8+6*11]*x7+ v[8+7*11]*x8 + v[8+8*11]*x9  + v[8+9*11]*x10  + v[8+10*11]*x11;
      sum10 += v[9+0]*x1 + v[9+11]*x2  + v[9+2*11]*x3  + v[9+3*11]*x4 + v[9+4*11]*x5 + v[9+5*11]*x6 + v[9+6*11]*x7+ v[9+7*11]*x8 + v[9+8*11]*x9  + v[9+9*11]*x10  + v[9+10*11]*x11;
      sum11 += v[10+0]*x1 + v[10+11]*x2  + v[10+2*11]*x3  + v[10+3*11]*x4 + v[10+4*11]*x5 + v[10+5*11]*x6 + v[10+6*11]*x7+ v[10+7*11]*x8 + v[10+8*11]*x9  + v[10+9*11]*x10  + v[10+10*11]*x11;
      v    += 121;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum6; z[8] = sum7; z[9] = sum8; z[10] = sum9; z[11] = sum10;
    if (!usecprow) {
      z += 11; y += 11;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayPair(yy,zz,&yarray,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(242.0*a->nz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultAdd_SeqBAIJ_N(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,*work,*workt,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,bs=A->rmap->bs,j,n,bs2=a->bs2;
  PetscInt          ncols,k;
  const PetscInt    *ridx = NULL,*idx,*ii;
  PetscBool         usecprow = a->compressedrow.use;

  PetscFunctionBegin;
  ierr = VecCopy(yy,zz);CHKERRQ(ierr);
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&zarray);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = zarray;
  }

  if (!a->mult_work) {
    k    = PetscMax(A->rmap->n,A->cmap->n);
    ierr = PetscMalloc1(k+1,&a->mult_work);CHKERRQ(ierr);
  }
  work = a->mult_work;
  for (i=0; i<mbs; i++) {
    n     = ii[1] - ii[0]; ii++;
    ncols = n*bs;
    workt = work;
    for (j=0; j<n; j++) {
      xb = x + bs*(*idx++);
      for (k=0; k<bs; k++) workt[k] = xb[k];
      workt += bs;
    }
    if (usecprow) z = zarray + bs*ridx[i];
    PetscKernel_w_gets_w_plus_Ar_times_v(bs,ncols,work,v,z);
    /* BLASgemv_("N",&bs,&ncols,&_DOne,v,&bs,work,&_One,&_DOne,z,&_One); */
    v += n*bs2;
    if (!usecprow) z += bs;
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&zarray);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*a->nz*bs2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultHermitianTranspose_SeqBAIJ(Mat A,Vec xx,Vec zz)
{
  PetscScalar    zero = 0.0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecSet(zz,zero);CHKERRQ(ierr);
  ierr = MatMultHermitianTransposeAdd_SeqBAIJ(A,xx,zz,zz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultTranspose_SeqBAIJ(Mat A,Vec xx,Vec zz)
{
  PetscScalar    zero = 0.0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecSet(zz,zero);CHKERRQ(ierr);
  ierr = MatMultTransposeAdd_SeqBAIJ(A,xx,zz,zz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultHermitianTransposeAdd_SeqBAIJ(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z,x1,x2,x3,x4,x5;
  const PetscScalar *x,*xb = NULL;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,rval,bs=A->rmap->bs,j,n;
  const PetscInt    *idx,*ii,*ib,*ridx = NULL;
  Mat_CompressedRow cprow = a->compressedrow;
  PetscBool         usecprow = cprow.use;

  PetscFunctionBegin;
  if (yy != zz) { ierr = VecCopy(yy,zz);CHKERRQ(ierr); }
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&z);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = cprow.nrows;
    ii   = cprow.i;
    ridx = cprow.rindex;
  } else {
    mbs=a->mbs;
    ii = a->i;
    xb = x;
  }

  switch (bs) {
  case 1:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + ridx[i];
      x1 = xb[0];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval     = ib[j];
        z[rval] += PetscConj(*v) * x1;
        v++;
      }
      if (!usecprow) xb++;
    }
    break;
  case 2:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 2*ridx[i];
      x1 = xb[0]; x2 = xb[1];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval       = ib[j]*2;
        z[rval++] += PetscConj(v[0])*x1 + PetscConj(v[1])*x2;
        z[rval++] += PetscConj(v[2])*x1 + PetscConj(v[3])*x2;
        v         += 4;
      }
      if (!usecprow) xb += 2;
    }
    break;
  case 3:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 3*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval       = ib[j]*3;
        z[rval++] += PetscConj(v[0])*x1 + PetscConj(v[1])*x2 + PetscConj(v[2])*x3;
        z[rval++] += PetscConj(v[3])*x1 + PetscConj(v[4])*x2 + PetscConj(v[5])*x3;
        z[rval++] += PetscConj(v[6])*x1 + PetscConj(v[7])*x2 + PetscConj(v[8])*x3;
        v         += 9;
      }
      if (!usecprow) xb += 3;
    }
    break;
  case 4:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 4*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval       = ib[j]*4;
        z[rval++] +=  PetscConj(v[0])*x1 + PetscConj(v[1])*x2  + PetscConj(v[2])*x3  + PetscConj(v[3])*x4;
        z[rval++] +=  PetscConj(v[4])*x1 + PetscConj(v[5])*x2  + PetscConj(v[6])*x3  + PetscConj(v[7])*x4;
        z[rval++] +=  PetscConj(v[8])*x1 + PetscConj(v[9])*x2  + PetscConj(v[10])*x3 + PetscConj(v[11])*x4;
        z[rval++] += PetscConj(v[12])*x1 + PetscConj(v[13])*x2 + PetscConj(v[14])*x3 + PetscConj(v[15])*x4;
        v         += 16;
      }
      if (!usecprow) xb += 4;
    }
    break;
  case 5:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 5*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      x4 = xb[3]; x5 = xb[4];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval       = ib[j]*5;
        z[rval++] +=  PetscConj(v[0])*x1 +  PetscConj(v[1])*x2 +  PetscConj(v[2])*x3 +  PetscConj(v[3])*x4 +  PetscConj(v[4])*x5;
        z[rval++] +=  PetscConj(v[5])*x1 +  PetscConj(v[6])*x2 +  PetscConj(v[7])*x3 +  PetscConj(v[8])*x4 +  PetscConj(v[9])*x5;
        z[rval++] += PetscConj(v[10])*x1 + PetscConj(v[11])*x2 + PetscConj(v[12])*x3 + PetscConj(v[13])*x4 + PetscConj(v[14])*x5;
        z[rval++] += PetscConj(v[15])*x1 + PetscConj(v[16])*x2 + PetscConj(v[17])*x3 + PetscConj(v[18])*x4 + PetscConj(v[19])*x5;
        z[rval++] += PetscConj(v[20])*x1 + PetscConj(v[21])*x2 + PetscConj(v[22])*x3 + PetscConj(v[23])*x4 + PetscConj(v[24])*x5;
        v         += 25;
      }
      if (!usecprow) xb += 5;
    }
    break;
  default: /* block sizes larger than 5 by 5 are handled by BLAS */
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"block size larger than 5 is not supported yet");
#if 0
    {
      PetscInt          ncols,k,bs2=a->bs2;
      PetscScalar       *work,*workt,zb;
      const PetscScalar *xtmp;
      if (!a->mult_work) {
        k    = PetscMax(A->rmap->n,A->cmap->n);
        ierr = PetscMalloc1(k+1,&a->mult_work);CHKERRQ(ierr);
      }
      work = a->mult_work;
      xtmp = x;
      for (i=0; i<mbs; i++) {
        n     = ii[1] - ii[0]; ii++;
        ncols = n*bs;
        ierr  = PetscArrayzero(work,ncols);CHKERRQ(ierr);
        if (usecprow) xtmp = x + bs*ridx[i];
        PetscKernel_w_gets_w_plus_trans_Ar_times_v(bs,ncols,xtmp,v,work);
        /* BLASgemv_("T",&bs,&ncols,&_DOne,v,&bs,xtmp,&_One,&_DOne,work,&_One); */
        v += n*bs2;
        if (!usecprow) xtmp += bs;
        workt = work;
        for (j=0; j<n; j++) {
          zb = z + bs*(*idx++);
          for (k=0; k<bs; k++) zb[k] += workt[k] ;
          workt += bs;
        }
      }
    }
#endif
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&z);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*a->nz*a->bs2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMultTransposeAdd_SeqBAIJ(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *zb,*z,x1,x2,x3,x4,x5;
  const PetscScalar *x,*xb = 0;
  const MatScalar   *v;
  PetscErrorCode    ierr;
  PetscInt          mbs,i,rval,bs=A->rmap->bs,j,n,bs2=a->bs2;
  const PetscInt    *idx,*ii,*ib,*ridx = NULL;
  Mat_CompressedRow cprow   = a->compressedrow;
  PetscBool         usecprow=cprow.use;

  PetscFunctionBegin;
  if (yy != zz) { ierr = VecCopy(yy,zz);CHKERRQ(ierr); }
  ierr = VecGetArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecGetArray(zz,&z);CHKERRQ(ierr);

  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = cprow.nrows;
    ii   = cprow.i;
    ridx = cprow.rindex;
  } else {
    mbs=a->mbs;
    ii = a->i;
    xb = x;
  }

  switch (bs) {
  case 1:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + ridx[i];
      x1 = xb[0];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval     = ib[j];
        z[rval] += *v * x1;
        v++;
      }
      if (!usecprow) xb++;
    }
    break;
  case 2:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 2*ridx[i];
      x1 = xb[0]; x2 = xb[1];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval       = ib[j]*2;
        z[rval++] += v[0]*x1 + v[1]*x2;
        z[rval++] += v[2]*x1 + v[3]*x2;
        v         += 4;
      }
      if (!usecprow) xb += 2;
    }
    break;
  case 3:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 3*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval       = ib[j]*3;
        z[rval++] += v[0]*x1 + v[1]*x2 + v[2]*x3;
        z[rval++] += v[3]*x1 + v[4]*x2 + v[5]*x3;
        z[rval++] += v[6]*x1 + v[7]*x2 + v[8]*x3;
        v         += 9;
      }
      if (!usecprow) xb += 3;
    }
    break;
  case 4:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 4*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval       = ib[j]*4;
        z[rval++] +=  v[0]*x1 +  v[1]*x2 +  v[2]*x3 +  v[3]*x4;
        z[rval++] +=  v[4]*x1 +  v[5]*x2 +  v[6]*x3 +  v[7]*x4;
        z[rval++] +=  v[8]*x1 +  v[9]*x2 + v[10]*x3 + v[11]*x4;
        z[rval++] += v[12]*x1 + v[13]*x2 + v[14]*x3 + v[15]*x4;
        v         += 16;
      }
      if (!usecprow) xb += 4;
    }
    break;
  case 5:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 5*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      x4 = xb[3]; x5 = xb[4];
      ib = idx + ii[0];
      n  = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval       = ib[j]*5;
        z[rval++] +=  v[0]*x1 +  v[1]*x2 +  v[2]*x3 +  v[3]*x4 +  v[4]*x5;
        z[rval++] +=  v[5]*x1 +  v[6]*x2 +  v[7]*x3 +  v[8]*x4 +  v[9]*x5;
        z[rval++] += v[10]*x1 + v[11]*x2 + v[12]*x3 + v[13]*x4 + v[14]*x5;
        z[rval++] += v[15]*x1 + v[16]*x2 + v[17]*x3 + v[18]*x4 + v[19]*x5;
        z[rval++] += v[20]*x1 + v[21]*x2 + v[22]*x3 + v[23]*x4 + v[24]*x5;
        v         += 25;
      }
      if (!usecprow) xb += 5;
    }
    break;
  default: {      /* block sizes larger then 5 by 5 are handled by BLAS */
    PetscInt          ncols,k;
    PetscScalar       *work,*workt;
    const PetscScalar *xtmp;
    if (!a->mult_work) {
      k    = PetscMax(A->rmap->n,A->cmap->n);
      ierr = PetscMalloc1(k+1,&a->mult_work);CHKERRQ(ierr);
    }
    work = a->mult_work;
    xtmp = x;
    for (i=0; i<mbs; i++) {
      n     = ii[1] - ii[0]; ii++;
      ncols = n*bs;
      ierr  = PetscArrayzero(work,ncols);CHKERRQ(ierr);
      if (usecprow) xtmp = x + bs*ridx[i];
      PetscKernel_w_gets_w_plus_trans_Ar_times_v(bs,ncols,xtmp,v,work);
      /* BLASgemv_("T",&bs,&ncols,&_DOne,v,&bs,xtmp,&_One,&_DOne,work,&_One); */
      v += n*bs2;
      if (!usecprow) xtmp += bs;
      workt = work;
      for (j=0; j<n; j++) {
        zb = z + bs*(*idx++);
        for (k=0; k<bs; k++) zb[k] += workt[k];
        workt += bs;
      }
    }
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(zz,&z);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*a->nz*a->bs2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatScale_SeqBAIJ(Mat inA,PetscScalar alpha)
{
  Mat_SeqBAIJ    *a      = (Mat_SeqBAIJ*)inA->data;
  PetscInt       totalnz = a->bs2*a->nz;
  PetscScalar    oalpha  = alpha;
  PetscErrorCode ierr;
  PetscBLASInt   one = 1,tnz;

  PetscFunctionBegin;
  ierr = PetscBLASIntCast(totalnz,&tnz);CHKERRQ(ierr);
  PetscStackCallBLAS("BLASscal",BLASscal_(&tnz,&oalpha,a->a,&one));
  ierr = PetscLogFlops(totalnz);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatNorm_SeqBAIJ(Mat A,NormType type,PetscReal *norm)
{
  PetscErrorCode ierr;
  Mat_SeqBAIJ    *a  = (Mat_SeqBAIJ*)A->data;
  MatScalar      *v  = a->a;
  PetscReal      sum = 0.0;
  PetscInt       i,j,k,bs=A->rmap->bs,nz=a->nz,bs2=a->bs2,k1;

  PetscFunctionBegin;
  if (type == NORM_FROBENIUS) {
#if defined(PETSC_USE_REAL___FP16)
    PetscBLASInt one = 1,cnt = bs2*nz;
    *norm = BLASnrm2_(&cnt,v,&one);
#else
    for (i=0; i<bs2*nz; i++) {
      sum += PetscRealPart(PetscConj(*v)*(*v)); v++;
    }
#endif
    *norm = PetscSqrtReal(sum);
    ierr = PetscLogFlops(2.0*bs2*nz);CHKERRQ(ierr);
  } else if (type == NORM_1) { /* maximum column sum */
    PetscReal *tmp;
    PetscInt  *bcol = a->j;
    ierr = PetscCalloc1(A->cmap->n+1,&tmp);CHKERRQ(ierr);
    for (i=0; i<nz; i++) {
      for (j=0; j<bs; j++) {
        k1 = bs*(*bcol) + j; /* column index */
        for (k=0; k<bs; k++) {
          tmp[k1] += PetscAbsScalar(*v); v++;
        }
      }
      bcol++;
    }
    *norm = 0.0;
    for (j=0; j<A->cmap->n; j++) {
      if (tmp[j] > *norm) *norm = tmp[j];
    }
    ierr = PetscFree(tmp);CHKERRQ(ierr);
    ierr = PetscLogFlops(PetscMax(bs2*nz-1,0));CHKERRQ(ierr);
  } else if (type == NORM_INFINITY) { /* maximum row sum */
    *norm = 0.0;
    for (k=0; k<bs; k++) {
      for (j=0; j<a->mbs; j++) {
        v   = a->a + bs2*a->i[j] + k;
        sum = 0.0;
        for (i=0; i<a->i[j+1]-a->i[j]; i++) {
          for (k1=0; k1<bs; k1++) {
            sum += PetscAbsScalar(*v);
            v   += bs;
          }
        }
        if (sum > *norm) *norm = sum;
      }
    }
    ierr = PetscLogFlops(PetscMax(bs2*nz-1,0));CHKERRQ(ierr);
  } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"No support for this norm yet");
  PetscFunctionReturn(0);
}


PetscErrorCode MatEqual_SeqBAIJ(Mat A,Mat B,PetscBool * flg)
{
  Mat_SeqBAIJ    *a = (Mat_SeqBAIJ*)A->data,*b = (Mat_SeqBAIJ*)B->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* If the  matrix/block dimensions are not equal, or no of nonzeros or shift */
  if ((A->rmap->N != B->rmap->N) || (A->cmap->n != B->cmap->n) || (A->rmap->bs != B->rmap->bs)|| (a->nz != b->nz)) {
    *flg = PETSC_FALSE;
    PetscFunctionReturn(0);
  }

  /* if the a->i are the same */
  ierr = PetscArraycmp(a->i,b->i,a->mbs+1,flg);CHKERRQ(ierr);
  if (!*flg) PetscFunctionReturn(0);

  /* if a->j are the same */
  ierr = PetscArraycmp(a->j,b->j,a->nz,flg);CHKERRQ(ierr);
  if (!*flg) PetscFunctionReturn(0);

  /* if a->a are the same */
  ierr = PetscArraycmp(a->a,b->a,(a->nz)*(A->rmap->bs)*(B->rmap->bs),flg);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}

PetscErrorCode MatGetDiagonal_SeqBAIJ(Mat A,Vec v)
{
  Mat_SeqBAIJ    *a = (Mat_SeqBAIJ*)A->data;
  PetscErrorCode ierr;
  PetscInt       i,j,k,n,row,bs,*ai,*aj,ambs,bs2;
  PetscScalar    *x,zero = 0.0;
  MatScalar      *aa,*aa_j;

  PetscFunctionBegin;
  if (A->factortype) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Not for factored matrix");
  bs   = A->rmap->bs;
  aa   = a->a;
  ai   = a->i;
  aj   = a->j;
  ambs = a->mbs;
  bs2  = a->bs2;

  ierr = VecSet(v,zero);CHKERRQ(ierr);
  ierr = VecGetArray(v,&x);CHKERRQ(ierr);
  ierr = VecGetLocalSize(v,&n);CHKERRQ(ierr);
  if (n != A->rmap->N) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Nonconforming matrix and vector");
  for (i=0; i<ambs; i++) {
    for (j=ai[i]; j<ai[i+1]; j++) {
      if (aj[j] == i) {
        row  = i*bs;
        aa_j = aa+j*bs2;
        for (k=0; k<bs2; k+=(bs+1),row++) x[row] = aa_j[k];
        break;
      }
    }
  }
  ierr = VecRestoreArray(v,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatDiagonalScale_SeqBAIJ(Mat A,Vec ll,Vec rr)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  const PetscScalar *l,*r,*li,*ri;
  PetscScalar       x;
  MatScalar         *aa, *v;
  PetscErrorCode    ierr;
  PetscInt          i,j,k,lm,rn,M,m,n,mbs,tmp,bs,bs2,iai;
  const PetscInt    *ai,*aj;

  PetscFunctionBegin;
  ai  = a->i;
  aj  = a->j;
  aa  = a->a;
  m   = A->rmap->n;
  n   = A->cmap->n;
  bs  = A->rmap->bs;
  mbs = a->mbs;
  bs2 = a->bs2;
  if (ll) {
    ierr = VecGetArrayRead(ll,&l);CHKERRQ(ierr);
    ierr = VecGetLocalSize(ll,&lm);CHKERRQ(ierr);
    if (lm != m) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Left scaling vector wrong length");
    for (i=0; i<mbs; i++) { /* for each block row */
      M  = ai[i+1] - ai[i];
      li = l + i*bs;
      v  = aa + bs2*ai[i];
      for (j=0; j<M; j++) { /* for each block */
        for (k=0; k<bs2; k++) {
          (*v++) *= li[k%bs];
        }
      }
    }
    ierr = VecRestoreArrayRead(ll,&l);CHKERRQ(ierr);
    ierr = PetscLogFlops(a->nz);CHKERRQ(ierr);
  }

  if (rr) {
    ierr = VecGetArrayRead(rr,&r);CHKERRQ(ierr);
    ierr = VecGetLocalSize(rr,&rn);CHKERRQ(ierr);
    if (rn != n) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Right scaling vector wrong length");
    for (i=0; i<mbs; i++) { /* for each block row */
      iai = ai[i];
      M   = ai[i+1] - iai;
      v   = aa + bs2*iai;
      for (j=0; j<M; j++) { /* for each block */
        ri = r + bs*aj[iai+j];
        for (k=0; k<bs; k++) {
          x = ri[k];
          for (tmp=0; tmp<bs; tmp++) v[tmp] *= x;
          v += bs;
        }
      }
    }
    ierr = VecRestoreArrayRead(rr,&r);CHKERRQ(ierr);
    ierr = PetscLogFlops(a->nz);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


PetscErrorCode MatGetInfo_SeqBAIJ(Mat A,MatInfoType flag,MatInfo *info)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;

  PetscFunctionBegin;
  info->block_size   = a->bs2;
  info->nz_allocated = a->bs2*a->maxnz;
  info->nz_used      = a->bs2*a->nz;
  info->nz_unneeded  = info->nz_allocated - info->nz_used;
  info->assemblies   = A->num_ass;
  info->mallocs      = A->info.mallocs;
  info->memory       = ((PetscObject)A)->mem;
  if (A->factortype) {
    info->fill_ratio_given  = A->info.fill_ratio_given;
    info->fill_ratio_needed = A->info.fill_ratio_needed;
    info->factor_mallocs    = A->info.factor_mallocs;
  } else {
    info->fill_ratio_given  = 0;
    info->fill_ratio_needed = 0;
    info->factor_mallocs    = 0;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MatZeroEntries_SeqBAIJ(Mat A)
{
  Mat_SeqBAIJ    *a = (Mat_SeqBAIJ*)A->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscArrayzero(a->a,a->bs2*a->i[a->mbs]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMatMultSymbolic_SeqBAIJ_SeqDense(Mat A,Mat B,PetscReal fill,Mat C)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatMatMultSymbolic_SeqDense_SeqDense(A,B,0.0,C);CHKERRQ(ierr);

  C->ops->matmultnumeric = MatMatMultNumeric_SeqBAIJ_SeqDense;
  PetscFunctionReturn(0);
}

PetscErrorCode MatMatMult_SeqBAIJ_1_Private(Mat A,PetscScalar* b,PetscInt bm,PetscScalar* c,PetscInt cm,PetscInt cn)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1;
  const PetscScalar *xb;
  PetscScalar       x1;
  const MatScalar   *v,*vv;
  PetscInt          mbs,i,*idx,*ii,j,*jj,n,k,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = c;
  }

  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+n,n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    if (usecprow) z = c + ridx[i];
    jj = idx;
    vv = v;
    for (k=0; k<cn; k++) {
      idx = jj;
      v = vv;
      sum1    = 0.0;
      for (j=0; j<n; j++) {
        xb    = b + (*idx++); x1 = xb[0+k*bm];
        sum1 += v[0]*x1;
        v    += 1;
      }
      z[0+k*cm] = sum1;
    }
    if (!usecprow) z += 1;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MatMatMult_SeqBAIJ_2_Private(Mat A,PetscScalar* b,PetscInt bm,PetscScalar* c,PetscInt cm,PetscInt cn)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2;
  const PetscScalar *xb;
  PetscScalar       x1,x2;
  const MatScalar   *v,*vv;
  PetscInt          mbs,i,*idx,*ii,j,*jj,n,k,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = c;
  }

  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+4*n,4*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    if (usecprow) z = c + 2*ridx[i];
    jj = idx;
    vv = v;
    for (k=0; k<cn; k++) {
      idx = jj;
      v = vv;
      sum1    = 0.0; sum2 = 0.0;
      for (j=0; j<n; j++) {
        xb    = b + 2*(*idx++); x1 = xb[0+k*bm]; x2 = xb[1+k*bm];
        sum1 += v[0]*x1 + v[2]*x2;
        sum2 += v[1]*x1 + v[3]*x2;
        v    += 4;
      }
      z[0+k*cm] = sum1; z[1+k*cm] = sum2;
    }
    if (!usecprow) z += 2;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MatMatMult_SeqBAIJ_3_Private(Mat A,PetscScalar* b,PetscInt bm,PetscScalar* c,PetscInt cm,PetscInt cn)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3;
  const PetscScalar *xb;
  PetscScalar       x1,x2,x3;
  const MatScalar   *v,*vv;
  PetscInt          mbs,i,*idx,*ii,j,*jj,n,k,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = c;
  }

  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+9*n,9*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    if (usecprow) z = c + 3*ridx[i];
    jj = idx;
    vv = v;
    for (k=0; k<cn; k++) {
      idx = jj;
      v = vv;
      sum1    = 0.0; sum2 = 0.0; sum3 = 0.0;
      for (j=0; j<n; j++) {
        xb    = b + 3*(*idx++); x1 = xb[0+k*bm]; x2 = xb[1+k*bm]; x3 = xb[2+k*bm];
        sum1 += v[0]*x1 + v[3]*x2 + v[6]*x3;
        sum2 += v[1]*x1 + v[4]*x2 + v[7]*x3;
        sum3 += v[2]*x1 + v[5]*x2 + v[8]*x3;
        v    += 9;
      }
      z[0+k*cm] = sum1; z[1+k*cm] = sum2; z[2+k*cm] = sum3;
    }
    if (!usecprow) z += 3;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MatMatMult_SeqBAIJ_4_Private(Mat A,PetscScalar* b,PetscInt bm,PetscScalar* c,PetscInt cm,PetscInt cn)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4;
  const PetscScalar *xb;
  PetscScalar       x1,x2,x3,x4;
  const MatScalar   *v,*vv;
  PetscInt          mbs,i,*idx,*ii,j,*jj,n,k,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = c;
  }

  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+16*n,16*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    if (usecprow) z = c + 4*ridx[i];
    jj = idx;
    vv = v;
    for (k=0; k<cn; k++) {
      idx = jj;
      v = vv;
      sum1    = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0;
      for (j=0; j<n; j++) {
        xb    = b + 4*(*idx++); x1 = xb[0+k*bm]; x2 = xb[1+k*bm]; x3 = xb[2+k*bm]; x4 = xb[3+k*bm];
        sum1 += v[0]*x1 + v[4]*x2 + v[8]*x3   + v[12]*x4;
        sum2 += v[1]*x1 + v[5]*x2 + v[9]*x3   + v[13]*x4;
        sum3 += v[2]*x1 + v[6]*x2 + v[10]*x3  + v[14]*x4;
        sum4 += v[3]*x1 + v[7]*x2 + v[11]*x3  + v[15]*x4;
        v    += 16;
      }
      z[0+k*cm] = sum1; z[1+k*cm] = sum2; z[2+k*cm] = sum3; z[3+k*cm] = sum4;
    }
    if (!usecprow) z += 4;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MatMatMult_SeqBAIJ_5_Private(Mat A,PetscScalar* b,PetscInt bm,PetscScalar* c,PetscInt cm,PetscInt cn)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar       *z = 0,sum1,sum2,sum3,sum4,sum5;
  const PetscScalar *xb;
  PetscScalar       x1,x2,x3,x4,x5;
  const MatScalar   *v,*vv;
  PetscInt          mbs,i,*idx,*ii,j,*jj,n,k,*ridx=NULL;
  PetscBool         usecprow=a->compressedrow.use;

  PetscFunctionBegin;
  idx = a->j;
  v   = a->a;
  if (usecprow) {
    mbs  = a->compressedrow.nrows;
    ii   = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii  = a->i;
    z   = c;
  }

  for (i=0; i<mbs; i++) {
    n           = ii[1] - ii[0]; ii++;
    PetscPrefetchBlock(idx+n,n,0,PETSC_PREFETCH_HINT_NTA);   /* Indices for the next row (assumes same size as this one) */
    PetscPrefetchBlock(v+25*n,25*n,0,PETSC_PREFETCH_HINT_NTA); /* Entries for the next row */
    if (usecprow) z = c + 5*ridx[i];
    jj = idx;
    vv = v;
    for (k=0; k<cn; k++) {
      idx = jj;
      v = vv;
      sum1    = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0;
      for (j=0; j<n; j++) {
        xb    = b + 5*(*idx++); x1 = xb[0+k*bm]; x2 = xb[1+k*bm]; x3 = xb[2+k*bm]; x4 = xb[3+k*bm]; x5 = xb[4+k*bm];
        sum1 += v[0]*x1 + v[5]*x2 + v[10]*x3  + v[15]*x4 + v[20]*x5;
        sum2 += v[1]*x1 + v[6]*x2 + v[11]*x3  + v[16]*x4 + v[21]*x5;
        sum3 += v[2]*x1 + v[7]*x2 + v[12]*x3  + v[17]*x4 + v[22]*x5;
        sum4 += v[3]*x1 + v[8]*x2 + v[13]*x3  + v[18]*x4 + v[23]*x5;
        sum5 += v[4]*x1 + v[9]*x2 + v[14]*x3  + v[19]*x4 + v[24]*x5;
        v    += 25;
      }
      z[0+k*cm] = sum1; z[1+k*cm] = sum2; z[2+k*cm] = sum3; z[3+k*cm] = sum4; z[4+k*cm] = sum5;
    }
    if (!usecprow) z += 5;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MatMatMultNumeric_SeqBAIJ_SeqDense(Mat A,Mat B,Mat C)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  Mat_SeqDense      *bd = (Mat_SeqDense*)B->data;
  Mat_SeqDense      *cd = (Mat_SeqDense*)C->data;
  PetscInt          cm=cd->lda,cn=B->cmap->n,bm=bd->lda;
  PetscInt          mbs,i,bs=A->rmap->bs,j,n,bs2=a->bs2;
  PetscBLASInt      bbs,bcn,bbm,bcm;
  PetscScalar       *z = 0;
  PetscScalar       *c,*b;
  const MatScalar   *v;
  const PetscInt    *idx,*ii,*ridx=NULL;
  PetscScalar       _DZero=0.0,_DOne=1.0;
  PetscBool         usecprow=a->compressedrow.use;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  if (!cm || !cn) PetscFunctionReturn(0);
  if (B->rmap->n != A->cmap->n) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Number columns in A %D not equal rows in B %D\n",A->cmap->n,B->rmap->n);
  if (A->rmap->n != C->rmap->n) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Number rows in C %D not equal rows in A %D\n",C->rmap->n,A->rmap->n);
  if (B->cmap->n != C->cmap->n) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Number columns in B %D not equal columns in C %D\n",B->cmap->n,C->cmap->n);
  b = bd->v;
  if (a->nonzerorowcnt != A->rmap->n) {
    ierr = MatZeroEntries(C);CHKERRQ(ierr);
  }
  ierr = MatDenseGetArray(C,&c);CHKERRQ(ierr);
  switch (bs) {
  case 1:
    ierr = MatMatMult_SeqBAIJ_1_Private(A, b, bm, c, cm, cn);
    break;
  case 2:
    ierr = MatMatMult_SeqBAIJ_2_Private(A, b, bm, c, cm, cn);
    break;
  case 3:
    ierr = MatMatMult_SeqBAIJ_3_Private(A, b, bm, c, cm, cn);
    break;
  case 4:
    ierr = MatMatMult_SeqBAIJ_4_Private(A, b, bm, c, cm, cn);
    break;
  case 5:
    ierr = MatMatMult_SeqBAIJ_5_Private(A, b, bm, c, cm, cn);
    break;
  default: /* block sizes larger than 5 by 5 are handled by BLAS */
    ierr = PetscBLASIntCast(bs,&bbs);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(cn,&bcn);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(bm,&bbm);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(cm,&bcm);CHKERRQ(ierr);
    idx = a->j;
    v   = a->a;
    if (usecprow) {
      mbs  = a->compressedrow.nrows;
      ii   = a->compressedrow.i;
      ridx = a->compressedrow.rindex;
    } else {
      mbs = a->mbs;
      ii  = a->i;
      z   = c;
    }
    for (i=0; i<mbs; i++) {
      n = ii[1] - ii[0]; ii++;
      if (usecprow) z = c + bs*ridx[i];
      if (n) {
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&bbs,&bcn,&bbs,&_DOne,v,&bbs,b+bs*(*idx++),&bbm,&_DZero,z,&bcm));
        v += bs2;
      }
      for (j=1; j<n; j++) {
        PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&bbs,&bcn,&bbs,&_DOne,v,&bbs,b+bs*(*idx++),&bbm,&_DOne,z,&bcm));
        v += bs2;
      }
      if (!usecprow) z += bs;
    }
  }
  ierr = MatDenseRestoreArray(C,&c);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = PetscLogFlops((2.0*a->nz*bs2 - bs*a->nonzerorowcnt)*cn);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
