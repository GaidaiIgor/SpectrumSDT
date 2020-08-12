#include <../src/mat/impls/baij/seq/baij.h>
#include <petsc/private/kernels/blockinvert.h>

/* Block operations are done by accessing one column at at time */
/* Default MatSolve for block size 14 */

PetscErrorCode MatSolve_SeqBAIJ_14_NaturalOrdering(Mat A,Vec bb,Vec xx)
{
  Mat_SeqBAIJ       *a=(Mat_SeqBAIJ*)A->data;
  PetscErrorCode    ierr;
  const PetscInt    n=a->mbs,*ai=a->i,*aj=a->j,*adiag=a->diag,*vi,bs=A->rmap->bs,bs2=a->bs2;
  PetscInt          i,k,nz,idx,idt,m;
  const MatScalar   *aa=a->a,*v;
  PetscScalar       s[14];
  PetscScalar       *x,xv;
  const PetscScalar *b;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);

  /* forward solve the lower triangular */
  for (i=0; i<n; i++) {
    v         = aa + bs2*ai[i];
    vi        = aj + ai[i];
    nz        = ai[i+1] - ai[i];
    idt       = bs*i;
    x[idt]    = b[idt];    x[1+idt]  = b[1+idt];  x[2+idt]  = b[2+idt];  x[3+idt]  = b[3+idt];  x[4+idt]  = b[4+idt];
    x[5+idt]  = b[5+idt];  x[6+idt]  = b[6+idt];  x[7+idt]  = b[7+idt];  x[8+idt]  = b[8+idt];  x[9+idt] = b[9+idt];
    x[10+idt] = b[10+idt]; x[11+idt] = b[11+idt]; x[12+idt] = b[12+idt]; x[13+idt] = b[13+idt];
    for (m=0; m<nz; m++) {
      idx = bs*vi[m];
      for (k=0; k<bs; k++) {
        xv         = x[k + idx];
        x[idt]    -= v[0]*xv;
        x[1+idt]  -= v[1]*xv;
        x[2+idt]  -= v[2]*xv;
        x[3+idt]  -= v[3]*xv;
        x[4+idt]  -= v[4]*xv;
        x[5+idt]  -= v[5]*xv;
        x[6+idt]  -= v[6]*xv;
        x[7+idt]  -= v[7]*xv;
        x[8+idt]  -= v[8]*xv;
        x[9+idt]  -= v[9]*xv;
        x[10+idt] -= v[10]*xv;
        x[11+idt] -= v[11]*xv;
        x[12+idt] -= v[12]*xv;
        x[13+idt] -= v[13]*xv;
        v         += 14;
      }
    }
  }
  /* backward solve the upper triangular */
  for (i=n-1; i>=0; i--) {
    v     = aa + bs2*(adiag[i+1]+1);
    vi    = aj + adiag[i+1]+1;
    nz    = adiag[i] - adiag[i+1] - 1;
    idt   = bs*i;
    s[0]  = x[idt];    s[1]  = x[1+idt];  s[2]  = x[2+idt];  s[3]  = x[3+idt];  s[4]  = x[4+idt];
    s[5]  = x[5+idt];  s[6]  = x[6+idt];  s[7]  = x[7+idt];  s[8]  = x[8+idt];  s[9]  = x[9+idt];
    s[10] = x[10+idt]; s[11] = x[11+idt]; s[12] = x[12+idt]; s[13] = x[13+idt];

    for (m=0; m<nz; m++) {
      idx = bs*vi[m];
      for (k=0; k<bs; k++) {
        xv     = x[k + idx];
        s[0]  -= v[0]*xv;
        s[1]  -= v[1]*xv;
        s[2]  -= v[2]*xv;
        s[3]  -= v[3]*xv;
        s[4]  -= v[4]*xv;
        s[5]  -= v[5]*xv;
        s[6]  -= v[6]*xv;
        s[7]  -= v[7]*xv;
        s[8]  -= v[8]*xv;
        s[9]  -= v[9]*xv;
        s[10] -= v[10]*xv;
        s[11] -= v[11]*xv;
        s[12] -= v[12]*xv;
        s[13] -= v[13]*xv;
        v     += 14;
      }
    }
    ierr = PetscArrayzero(x+idt,bs);CHKERRQ(ierr);
    for (k=0; k<bs; k++) {
      x[idt]    += v[0]*s[k];
      x[1+idt]  += v[1]*s[k];
      x[2+idt]  += v[2]*s[k];
      x[3+idt]  += v[3]*s[k];
      x[4+idt]  += v[4]*s[k];
      x[5+idt]  += v[5]*s[k];
      x[6+idt]  += v[6]*s[k];
      x[7+idt]  += v[7]*s[k];
      x[8+idt]  += v[8]*s[k];
      x[9+idt]  += v[9]*s[k];
      x[10+idt] += v[10]*s[k];
      x[11+idt] += v[11]*s[k];
      x[12+idt] += v[12]*s[k];
      x[13+idt] += v[13]*s[k];
      v         += 14;
    }
  }
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*bs2*(a->nz) - bs*A->cmap->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Block operations are done by accessing one column at at time */
/* Default MatSolve for block size 13 */

PetscErrorCode MatSolve_SeqBAIJ_13_NaturalOrdering(Mat A,Vec bb,Vec xx)
{
  Mat_SeqBAIJ       *a=(Mat_SeqBAIJ*)A->data;
  PetscErrorCode    ierr;
  const PetscInt    n=a->mbs,*ai=a->i,*aj=a->j,*adiag=a->diag,*vi,bs=A->rmap->bs,bs2=a->bs2;
  PetscInt          i,k,nz,idx,idt,m;
  const MatScalar   *aa=a->a,*v;
  PetscScalar       s[13];
  PetscScalar       *x,xv;
  const PetscScalar *b;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);

  /* forward solve the lower triangular */
  for (i=0; i<n; i++) {
    v         = aa + bs2*ai[i];
    vi        = aj + ai[i];
    nz        = ai[i+1] - ai[i];
    idt       = bs*i;
    x[idt]    = b[idt];    x[1+idt]  = b[1+idt];  x[2+idt]  = b[2+idt];  x[3+idt]  = b[3+idt];  x[4+idt]  = b[4+idt];
    x[5+idt]  = b[5+idt];  x[6+idt]  = b[6+idt];  x[7+idt]  = b[7+idt];  x[8+idt]  = b[8+idt];  x[9+idt] = b[9+idt];
    x[10+idt] = b[10+idt]; x[11+idt] = b[11+idt]; x[12+idt] = b[12+idt];
    for (m=0; m<nz; m++) {
      idx = bs*vi[m];
      for (k=0; k<bs; k++) {
        xv         = x[k + idx];
        x[idt]    -= v[0]*xv;
        x[1+idt]  -= v[1]*xv;
        x[2+idt]  -= v[2]*xv;
        x[3+idt]  -= v[3]*xv;
        x[4+idt]  -= v[4]*xv;
        x[5+idt]  -= v[5]*xv;
        x[6+idt]  -= v[6]*xv;
        x[7+idt]  -= v[7]*xv;
        x[8+idt]  -= v[8]*xv;
        x[9+idt]  -= v[9]*xv;
        x[10+idt] -= v[10]*xv;
        x[11+idt] -= v[11]*xv;
        x[12+idt] -= v[12]*xv;
        v         += 13;
      }
    }
  }
  /* backward solve the upper triangular */
  for (i=n-1; i>=0; i--) {
    v     = aa + bs2*(adiag[i+1]+1);
    vi    = aj + adiag[i+1]+1;
    nz    = adiag[i] - adiag[i+1] - 1;
    idt   = bs*i;
    s[0]  = x[idt];    s[1]  = x[1+idt];  s[2]  = x[2+idt];  s[3]  = x[3+idt];  s[4]  = x[4+idt];
    s[5]  = x[5+idt];  s[6]  = x[6+idt];  s[7]  = x[7+idt];  s[8]  = x[8+idt];  s[9]  = x[9+idt];
    s[10] = x[10+idt]; s[11] = x[11+idt]; s[12] = x[12+idt];

    for (m=0; m<nz; m++) {
      idx = bs*vi[m];
      for (k=0; k<bs; k++) {
        xv     = x[k + idx];
        s[0]  -= v[0]*xv;
        s[1]  -= v[1]*xv;
        s[2]  -= v[2]*xv;
        s[3]  -= v[3]*xv;
        s[4]  -= v[4]*xv;
        s[5]  -= v[5]*xv;
        s[6]  -= v[6]*xv;
        s[7]  -= v[7]*xv;
        s[8]  -= v[8]*xv;
        s[9]  -= v[9]*xv;
        s[10] -= v[10]*xv;
        s[11] -= v[11]*xv;
        s[12] -= v[12]*xv;
        v     += 13;
      }
    }
    ierr = PetscArrayzero(x+idt,bs);CHKERRQ(ierr);
    for (k=0; k<bs; k++) {
      x[idt]    += v[0]*s[k];
      x[1+idt]  += v[1]*s[k];
      x[2+idt]  += v[2]*s[k];
      x[3+idt]  += v[3]*s[k];
      x[4+idt]  += v[4]*s[k];
      x[5+idt]  += v[5]*s[k];
      x[6+idt]  += v[6]*s[k];
      x[7+idt]  += v[7]*s[k];
      x[8+idt]  += v[8]*s[k];
      x[9+idt]  += v[9]*s[k];
      x[10+idt] += v[10]*s[k];
      x[11+idt] += v[11]*s[k];
      x[12+idt] += v[12]*s[k];
      v         += 13;
    }
  }
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*bs2*(a->nz) - bs*A->cmap->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Block operations are done by accessing one column at at time */
/* Default MatSolve for block size 12 */

PetscErrorCode MatSolve_SeqBAIJ_12_NaturalOrdering(Mat A,Vec bb,Vec xx)
{
  Mat_SeqBAIJ       *a=(Mat_SeqBAIJ*)A->data;
  PetscErrorCode    ierr;
  const PetscInt    n=a->mbs,*ai=a->i,*aj=a->j,*adiag=a->diag,*vi,bs=A->rmap->bs,bs2=a->bs2;
  PetscInt          i,k,nz,idx,idt,m;
  const MatScalar   *aa=a->a,*v;
  PetscScalar       s[12];
  PetscScalar       *x,xv;
  const PetscScalar *b;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecGetArray(xx,&x);CHKERRQ(ierr);

  /* forward solve the lower triangular */
  for (i=0; i<n; i++) {
    v         = aa + bs2*ai[i];
    vi        = aj + ai[i];
    nz        = ai[i+1] - ai[i];
    idt       = bs*i;
    x[idt]    = b[idt];    x[1+idt]  = b[1+idt];  x[2+idt]  = b[2+idt];  x[3+idt]  = b[3+idt];  x[4+idt]  = b[4+idt];
    x[5+idt]  = b[5+idt];  x[6+idt]  = b[6+idt];  x[7+idt]  = b[7+idt];  x[8+idt]  = b[8+idt];  x[9+idt] = b[9+idt];
    x[10+idt] = b[10+idt]; x[11+idt] = b[11+idt];
    for (m=0; m<nz; m++) {
      idx = bs*vi[m];
      for (k=0; k<bs; k++) {
        xv         = x[k + idx];
        x[idt]    -= v[0]*xv;
        x[1+idt]  -= v[1]*xv;
        x[2+idt]  -= v[2]*xv;
        x[3+idt]  -= v[3]*xv;
        x[4+idt]  -= v[4]*xv;
        x[5+idt]  -= v[5]*xv;
        x[6+idt]  -= v[6]*xv;
        x[7+idt]  -= v[7]*xv;
        x[8+idt]  -= v[8]*xv;
        x[9+idt]  -= v[9]*xv;
        x[10+idt] -= v[10]*xv;
        x[11+idt] -= v[11]*xv;
        v         += 12;
      }
    }
  }
  /* backward solve the upper triangular */
  for (i=n-1; i>=0; i--) {
    v     = aa + bs2*(adiag[i+1]+1);
    vi    = aj + adiag[i+1]+1;
    nz    = adiag[i] - adiag[i+1] - 1;
    idt   = bs*i;
    s[0]  = x[idt];    s[1]  = x[1+idt];  s[2]  = x[2+idt];  s[3]  = x[3+idt];  s[4]  = x[4+idt];
    s[5]  = x[5+idt];  s[6]  = x[6+idt];  s[7]  = x[7+idt];  s[8]  = x[8+idt];  s[9]  = x[9+idt];
    s[10] = x[10+idt]; s[11] = x[11+idt];

    for (m=0; m<nz; m++) {
      idx = bs*vi[m];
      for (k=0; k<bs; k++) {
        xv     = x[k + idx];
        s[0]  -= v[0]*xv;
        s[1]  -= v[1]*xv;
        s[2]  -= v[2]*xv;
        s[3]  -= v[3]*xv;
        s[4]  -= v[4]*xv;
        s[5]  -= v[5]*xv;
        s[6]  -= v[6]*xv;
        s[7]  -= v[7]*xv;
        s[8]  -= v[8]*xv;
        s[9]  -= v[9]*xv;
        s[10] -= v[10]*xv;
        s[11] -= v[11]*xv;
        v     += 12;
      }
    }
    ierr = PetscArrayzero(x+idt,bs);CHKERRQ(ierr);
    for (k=0; k<bs; k++) {
      x[idt]    += v[0]*s[k];
      x[1+idt]  += v[1]*s[k];
      x[2+idt]  += v[2]*s[k];
      x[3+idt]  += v[3]*s[k];
      x[4+idt]  += v[4]*s[k];
      x[5+idt]  += v[5]*s[k];
      x[6+idt]  += v[6]*s[k];
      x[7+idt]  += v[7]*s[k];
      x[8+idt]  += v[8]*s[k];
      x[9+idt]  += v[9]*s[k];
      x[10+idt] += v[10]*s[k];
      x[11+idt] += v[11]*s[k];
      v         += 12;
    }
  }
  ierr = VecRestoreArrayRead(bb,&b);CHKERRQ(ierr);
  ierr = VecRestoreArray(xx,&x);CHKERRQ(ierr);
  ierr = PetscLogFlops(2.0*bs2*(a->nz) - bs*A->cmap->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
