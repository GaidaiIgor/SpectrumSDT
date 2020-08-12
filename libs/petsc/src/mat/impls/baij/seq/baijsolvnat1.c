#include <../src/mat/impls/baij/seq/baij.h>
#include <petsc/private/kernels/blockinvert.h>


/*
      Special case where the matrix was ILU(0) factored in the natural
   ordering. This eliminates the need for the column and row permutation.
*/
PetscErrorCode MatSolve_SeqBAIJ_1_NaturalOrdering_inplace(Mat A,Vec bb,Vec xx)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  const PetscInt    n  = a->mbs,*vi,*ai=a->i,*aj=a->j,*diag=a->diag;
  PetscErrorCode    ierr;
  const MatScalar   *aa=a->a,*v;
  PetscScalar       *x;
  const PetscScalar *b;
  PetscScalar       s1,x1;
  PetscInt          jdx,idt,idx,nz,i;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);

  /* forward solve the lower triangular */
  idx  = 0;
  x[0] = b[0];
  for (i=1; i<n; i++) {
    v    =  aa      + ai[i];
    vi   =  aj      + ai[i];
    nz   =  diag[i] - ai[i];
    idx +=  1;
    s1   =  b[idx];
    while (nz--) {
      jdx = *vi++;
      x1  = x[jdx];
      s1 -= v[0]*x1;
      v  += 1;
    }
    x[idx] = s1;
  }
  /* backward solve the upper triangular */
  for (i=n-1; i>=0; i--) {
    v   = aa + diag[i] + 1;
    vi  = aj + diag[i] + 1;
    nz  = ai[i+1] - diag[i] - 1;
    idt = i;
    s1  = x[idt];
    while (nz--) {
      idx = *vi++;
      x1  = x[idx];
      s1 -= v[0]*x1;
      v  += 1;
    }
    v      = aa +  diag[i];
    x[idt] = v[0]*s1;
  }
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*(a->nz) - A->cmap->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode MatForwardSolve_SeqBAIJ_1_NaturalOrdering(Mat A,Vec bb,Vec xx)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscErrorCode    ierr;
  const PetscInt    n = a->mbs,*ai = a->i,*aj = a->j,*vi;
  PetscScalar       *x,sum;
  const PetscScalar *b;
  const MatScalar   *aa = a->a,*v;
  PetscInt          i,nz;

  PetscFunctionBegin;
  if (!n) PetscFunctionReturn(0);

  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);

  /* forward solve the lower triangular */
  x[0] = b[0];
  v    = aa;
  vi   = aj;
  for (i=1; i<n; i++) {
    nz  = ai[i+1] - ai[i];
    sum = b[i];
    PetscSparseDenseMinusDot(sum,x,v,vi,nz);
    v   += nz;
    vi  += nz;
    x[i] = sum;
  }
  ierr = PetscLogFlops(a->nz - A->cmap->n);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatBackwardSolve_SeqBAIJ_1_NaturalOrdering(Mat A,Vec bb,Vec xx)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscErrorCode    ierr;
  const PetscInt    n = a->mbs,*aj = a->j,*adiag = a->diag,*vi;
  PetscScalar       *x,sum;
  const PetscScalar *b;
  const MatScalar   *aa = a->a,*v;
  PetscInt          i,nz;

  PetscFunctionBegin;
  if (!n) PetscFunctionReturn(0);

  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);

  /* backward solve the upper triangular */
  for (i=n-1; i>=0; i--) {
    v   = aa + adiag[i+1] + 1;
    vi  = aj + adiag[i+1] + 1;
    nz  = adiag[i] - adiag[i+1]-1;
    sum = b[i];
    PetscSparseDenseMinusDot(sum,x,v,vi,nz);
    x[i] = sum*v[nz]; /* x[i]=aa[adiag[i]]*sum; v++; */
  }

  ierr = PetscLogFlops(2.0*a->nz - A->cmap->n);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatSolve_SeqBAIJ_1_NaturalOrdering(Mat A,Vec bb,Vec xx)
{
  Mat_SeqBAIJ       *a = (Mat_SeqBAIJ*)A->data;
  PetscErrorCode    ierr;
  const PetscInt    n = a->mbs,*ai = a->i,*aj = a->j,*adiag = a->diag,*vi;
  PetscScalar       *x,sum;
  const PetscScalar *b;
  const MatScalar   *aa = a->a,*v;
  PetscInt          i,nz;

  PetscFunctionBegin;
  if (!n) PetscFunctionReturn(0);

  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);

  /* forward solve the lower triangular */
  x[0] = b[0];
  v    = aa;
  vi   = aj;
  for (i=1; i<n; i++) {
    nz  = ai[i+1] - ai[i];
    sum = b[i];
    PetscSparseDenseMinusDot(sum,x,v,vi,nz);
    v   += nz;
    vi  += nz;
    x[i] = sum;
  }

  /* backward solve the upper triangular */
  for (i=n-1; i>=0; i--) {
    v   = aa + adiag[i+1] + 1;
    vi  = aj + adiag[i+1] + 1;
    nz  = adiag[i] - adiag[i+1]-1;
    sum = x[i];
    PetscSparseDenseMinusDot(sum,x,v,vi,nz);
    x[i] = sum*v[nz]; /* x[i]=aa[adiag[i]]*sum; v++; */
  }

  ierr = PetscLogFlops(2.0*a->nz - A->cmap->n);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

