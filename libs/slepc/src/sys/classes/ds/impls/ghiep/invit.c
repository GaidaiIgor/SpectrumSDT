/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/dsimpl.h>
#include <slepcblaslapack.h>

struct HRtr
{
  PetscScalar *data;
  PetscInt    m;
  PetscInt    idx[2];
  PetscInt    n[2];
  PetscScalar tau[2];
  PetscReal   alpha;
  PetscReal   cs;
  PetscReal   sn;
  PetscInt    type;
};

/*
  Generates a hyperbolic rotation
    if x1*x1 - x2*x2 != 0
      r = sqrt(|x1*x1 - x2*x2|)
      c = x1/r  s = x2/r

      | c -s||x1|   |d*r|
      |-s  c||x2| = | 0 |
      where d = 1 for type==1 and -1 for type==2
  Returns the condition number of the reduction
*/
static PetscErrorCode HRGen(PetscReal x1,PetscReal x2,PetscInt *type,PetscReal *c,PetscReal *s,PetscReal *r,PetscReal *cond)
{
  PetscReal t,n2,xa,xb;
  PetscInt  type_;

  PetscFunctionBegin;
  if (x2==0.0) {
    *r = PetscAbsReal(x1); *c = (x1>=0.0)?1.0:-1.0; *s = 0.0;
    if (type) *type = 1;
    PetscFunctionReturn(0);
  }
  if (PetscAbsReal(x1) == PetscAbsReal(x2)) {
    /* hyperbolic rotation doesn't exist */
    *c = *s = *r = 0.0;
    if (type) *type = 0;
    *cond = PETSC_MAX_REAL;
    PetscFunctionReturn(0);
  }

  if (PetscAbsReal(x1)>PetscAbsReal(x2)) {
    xa = x1; xb = x2; type_ = 1;
  } else {
    xa = x2; xb = x1; type_ = 2;
  }
  t = xb/xa;
  n2 = PetscAbsReal(1 - t*t);
  *r = PetscSqrtReal(n2)*PetscAbsReal(xa);
  *c = x1/(*r);
  *s = x2/(*r);
  if (type_ == 2) *r *= -1;
  if (type) *type = type_;
  if (cond) *cond = (PetscAbsReal(*c) + PetscAbsReal(*s))/PetscAbsReal(PetscAbsReal(*c) - PetscAbsReal(*s));
  PetscFunctionReturn(0);
}

/*
                                |c  s|
  Applies an hyperbolic rotator |s  c|
           |c  s|
    [x1 x2]|s  c|
*/
static PetscErrorCode HRApply(PetscInt n,PetscScalar *x1,PetscInt inc1,PetscScalar *x2,PetscInt inc2,PetscReal c,PetscReal s)
{
  PetscInt    i;
  PetscReal   t;
  PetscScalar tmp;

  PetscFunctionBegin;
  if (PetscAbsReal(c)>PetscAbsReal(s)) { /* Type I */
    t = s/c;
    for (i=0;i<n;i++) {
      x1[i*inc1] = c*x1[i*inc1] + s*x2[i*inc2];
      x2[i*inc2] = t*x1[i*inc1] + x2[i*inc2]/c;
    }
  } else { /* Type II */
    t = c/s;
    for (i=0;i<n;i++) {
      tmp = x1[i*inc1];
      x1[i*inc1] = c*x1[i*inc1] + s*x2[i*inc2];
      x2[i*inc2] = t*x1[i*inc1] + tmp/s;
    }
  }
  PetscFunctionReturn(0);
}

/*
  Reduction to tridiagonal-diagonal form (see F. Tisseur, SIMAX 26(1), 2004).

  Input:
    A symmetric (only lower triangular part is referred)
    s vector +1 and -1 (signature matrix)
  Output:
    d,e
    s
    Q s-orthogonal matrix with Q^T*A*Q = T (symmetric tridiagonal matrix)
*/
static PetscErrorCode TridiagDiag_HHR(PetscInt n,PetscScalar *A,PetscInt lda,PetscReal *s,PetscScalar* Q,PetscInt ldq,PetscBool flip,PetscReal *d,PetscReal *e,PetscInt *perm_,PetscScalar *work,PetscReal *rwork,PetscBLASInt *iwork)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,*ii,*jj,i0=0,ik=0,tmp,type;
  PetscInt       nwu=0;
  PetscReal      *ss,cond=1.0,cs,sn,r;
  PetscScalar    tau,t,*AA;
  PetscBLASInt   n0,n1,ni,inc=1,m,n_,lda_,ldq_,*perm;
  PetscBool      breakdown = PETSC_TRUE;

  PetscFunctionBegin;
  if (n<3) {
    if (n==1) Q[0]=1;
    if (n==2) {
      Q[0] = Q[1+ldq] = 1;
      Q[1] = Q[ldq] = 0;
    }
    PetscFunctionReturn(0);
  }
  ierr = PetscBLASIntCast(lda,&lda_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldq,&ldq_);CHKERRQ(ierr);
  ss = rwork;
  perm = iwork;
  AA = work;
  for (i=0;i<n;i++) {
    ierr = PetscArraycpy(AA+i*n,A+i*lda,n);CHKERRQ(ierr);
  }
  nwu += n*n;
  k=0;
  while (breakdown && k<n) {
    breakdown = PETSC_FALSE;
    /* Classify (and flip) A and s according to sign */
    if (flip) {
      for (i=0;i<n;i++) {
        perm[i] = n-1-perm_[i];
        if (perm[i]==0) i0 = i;
        if (perm[i]==k) ik = i;
      }
    } else {
      for (i=0;i<n;i++) {
        perm[i] = perm_[i];
        if (perm[i]==0) i0 = i;
        if (perm[i]==k) ik = i;
      }
    }
    perm[ik] = 0;
    perm[i0] = k;
    i=1;
    while (i<n-1 && s[perm[i-1]]==s[perm[0]]) {
      if (s[perm[i]]!=s[perm[0]]) {
        j=i+1;
        while (j<n-1 && s[perm[j]]!=s[perm[0]])j++;
        tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
      }
      i++;
    }
    for (i=0;i<n;i++) {
      ss[i] = s[perm[i]];
    }
    if (flip) {
      ii = &j;
      jj = &i;
    } else {
      ii = &i;
      jj = &j;
    }
    for (i=0;i<n;i++)
      for (j=0;j<n;j++)
        A[i+j*lda] = AA[perm[*ii]+perm[*jj]*n];
    /* Initialize Q */
    for (i=0;i<n;i++) {
      ierr = PetscArrayzero(Q+i*ldq,n);CHKERRQ(ierr);
      Q[perm[i]+i*ldq] = 1.0;
    }
    for (ni=1;ni<n && ss[ni]==ss[0]; ni++);
    n0 = ni-1;
    n1 = n_-ni;
    for (j=0;j<n-2;j++) {
      ierr = PetscBLASIntCast(n-j-1,&m);CHKERRQ(ierr);
      /* Forming and applying reflectors */
      if (n0 > 1) {
        PetscStackCallBLAS("LAPACKlarfg",LAPACKlarfg_(&n0,A+ni-n0+j*lda,A+ni-n0+j*lda+1,&inc,&tau));
        /* Apply reflector */
        if (PetscAbsScalar(tau) != 0.0) {
          t=*(A+ni-n0+j*lda);  *(A+ni-n0+j*lda)=1.0;
          PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("R",&m,&n0,A+ni-n0+j*lda,&inc,&tau,A+j+1+(j+1)*lda,&lda_,work+nwu));
          PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n0,&m,A+ni-n0+j*lda,&inc,&tau,A+j+1+(j+1)*lda,&lda_,work+nwu));
          /* Update Q */
          PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("R",&n_,&n0,A+ni-n0+j*lda,&inc,&tau,Q+(j+1)*ldq,&ldq_,work+nwu));
          *(A+ni-n0+j*lda) = t;
          for (i=1;i<n0;i++) {
            *(A+ni-n0+j*lda+i) = 0.0;  *(A+j+(ni-n0+i)*lda) = 0.0;
          }
          *(A+j+(ni-n0)*lda) = *(A+ni-n0+j*lda);
        }
      }
      if (n1 > 1) {
        PetscStackCallBLAS("LAPACKlarfg",LAPACKlarfg_(&n1,A+n-n1+j*lda,A+n-n1+j*lda+1,&inc,&tau));
        /* Apply reflector */
        if (PetscAbsScalar(tau) != 0.0) {
          t=*(A+n-n1+j*lda);  *(A+n-n1+j*lda)=1.0;
          PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("R",&m,&n1,A+n-n1+j*lda,&inc,&tau,A+j+1+(n-n1)*lda,&lda_,work+nwu));
          PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n1,&m,A+n-n1+j*lda,&inc,&tau,A+n-n1+(j+1)*lda,&lda_,work+nwu));
          /* Update Q */
          PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("R",&n_,&n1,A+n-n1+j*lda,&inc,&tau,Q+(n-n1)*ldq,&ldq_,work+nwu));
          *(A+n-n1+j*lda) = t;
          for (i=1;i<n1;i++) {
            *(A+n-n1+i+j*lda) = 0.0;  *(A+j+(n-n1+i)*lda) = 0.0;
          }
          *(A+j+(n-n1)*lda) = *(A+n-n1+j*lda);
        }
      }
      /* Hyperbolic rotation */
      if (n0 > 0 && n1 > 0) {
        ierr = HRGen(PetscRealPart(A[ni-n0+j*lda]),PetscRealPart(A[n-n1+j*lda]),&type,&cs,&sn,&r,&cond);CHKERRQ(ierr);
        /* Check condition number */
        if (cond > 1.0/(10*PETSC_SQRT_MACHINE_EPSILON)) {
          breakdown = PETSC_TRUE;
          k++;
          if (k==n || flip) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Breakdown in construction of hyperbolic transformation");
          break;
        }
        A[ni-n0+j*lda] = r; A[n-n1+j*lda] = 0.0;
        A[j+(ni-n0)*lda] = r; A[j+(n-n1)*lda] = 0.0;
        /* Apply to A */
        ierr = HRApply(m,A+j+1+(ni-n0)*lda,1,A+j+1+(n-n1)*lda,1,cs,-sn);CHKERRQ(ierr);
        ierr = HRApply(m,A+ni-n0+(j+1)*lda,lda,A+n-n1+(j+1)*lda,lda,cs,-sn);CHKERRQ(ierr);

        /* Update Q */
        ierr = HRApply(n,Q+(ni-n0)*ldq,1,Q+(n-n1)*ldq,1,cs,-sn);CHKERRQ(ierr);
        if (type==2) {
          ss[ni-n0] = -ss[ni-n0]; ss[n-n1] = -ss[n-n1];
          n0++;ni++;n1--;
        }
      }
      if (n0>0) n0--;
      else n1--;
    }
  }

  /* flip matrices */
  if (flip) {
    for (i=0;i<n-1;i++) {
      d[i] = PetscRealPart(A[n-i-1+(n-i-1)*lda]);
      e[i] = PetscRealPart(A[n-i-1+(n-i-2)*lda]);
      s[i] = ss[n-i-1];
    }
    s[n-1] = ss[0];
    d[n-1] = PetscRealPart(A[0]);
    for (i=0;i<n;i++) {
      ierr=PetscArraycpy(work+i*n,Q+i*ldq,n);CHKERRQ(ierr);
    }
    for (i=0;i<n;i++)
      for (j=0;j<n;j++)
        Q[i+j*ldq] = work[i+(n-j-1)*n];
  } else {
    for (i=0;i<n-1;i++) {
      d[i] = PetscRealPart(A[i+i*lda]);
      e[i] = PetscRealPart(A[i+1+i*lda]);
      s[i] = ss[i];
    }
    s[n-1] = ss[n-1];
    d[n-1] = PetscRealPart(A[n-1 + (n-1)*lda]);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MadeHRtr(PetscInt sz,PetscInt n,PetscInt idx0,PetscInt n0,PetscInt idx1,PetscInt n1,struct HRtr *tr1,struct HRtr *tr2,PetscReal *ncond,PetscScalar *work)
{
  PetscErrorCode ierr;
  PetscScalar    *x,*y;
  PetscReal      ncond2=1.0;
  PetscBLASInt   n0_,n1_,inc=1;

  PetscFunctionBegin;
  /* Hyperbolic transformation to make zeros in x */
  x = tr1->data;
  tr1->n[0] = n0;
  tr1->n[1] = n1;
  tr1->idx[0] = idx0;
  tr1->idx[1] = idx1;
  ierr = PetscBLASIntCast(tr1->n[0],&n0_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(tr1->n[1],&n1_);CHKERRQ(ierr);
  if (tr1->n[0] > 1) {
    PetscStackCallBLAS("LAPACKlarfg",LAPACKlarfg_(&n0_,x+tr1->idx[0],x+tr1->idx[0]+1,&inc,tr1->tau));
  }
  if (tr1->n[1]> 1) {
    PetscStackCallBLAS("LAPACKlarfg",LAPACKlarfg_(&n1_,x+tr1->idx[1],x+tr1->idx[1]+1,&inc,tr1->tau+1));
  }
  if (tr1->idx[0]<tr1->idx[1]) {
    ierr = HRGen(PetscRealPart(x[tr1->idx[0]]),PetscRealPart(x[tr1->idx[1]]),&(tr1->type),&(tr1->cs),&(tr1->sn),&(tr1->alpha),ncond);CHKERRQ(ierr);
  } else {
    tr1->alpha = PetscRealPart(x[tr1->idx[0]]);
    *ncond = 1.0;
  }
  if (sz==2) {
    y = tr2->data;
    /* Apply first transformation to second column */
    if (tr1->n[0] > 1 && PetscAbsScalar(tr1->tau[0])!=0.0) {
      x[tr1->idx[0]] = 1.0;
      PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n0_,&inc,x+tr1->idx[0],&inc,tr1->tau,y+tr1->idx[0],&n0_,work));
    }
    if (tr1->n[1] > 1 && PetscAbsScalar(tr1->tau[1])!=0.0) {
      x[tr1->idx[1]] = 1.0;
      PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n1_,&inc,x+tr1->idx[1],&inc,tr1->tau+1,y+tr1->idx[1],&n1_,work));
    }
    if (tr1->idx[0]<tr1->idx[1]) {
      ierr = HRApply(1,y+tr1->idx[0],1,y+tr1->idx[1],1,tr1->cs,-tr1->sn);CHKERRQ(ierr);
    }
    tr2->n[0] = tr1->n[0];
    tr2->n[1] = tr1->n[1];
    tr2->idx[0] = tr1->idx[0];
    tr2->idx[1] = tr1->idx[1];
    if (tr1->idx[0]<tr1->idx[1] && tr1->type==2) {
      tr2->idx[1]++; tr2->n[1]--; tr2->n[0]++;
    }
    if (tr2->n[0]>0) {
      tr2->n[0]--; tr2->idx[0]++;
      if (tr2->n[1]==0) tr2->idx[1] = tr2->idx[0];
    } else {
      tr2->n[1]--; tr2->idx[1]++; tr2->idx[0] = tr2->idx[1];
    }
    /* Hyperbolic transformation to make zeros in y */
    ierr = PetscBLASIntCast(tr2->n[0],&n0_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(tr2->n[1],&n1_);CHKERRQ(ierr);
    if (tr2->n[0] > 1) {
      PetscStackCallBLAS("LAPACKlarfg",LAPACKlarfg_(&n0_,y+tr2->idx[0],y+tr2->idx[0]+1,&inc,tr2->tau));
    }
    if (tr2->n[1]> 1) {
      PetscStackCallBLAS("LAPACKlarfg",LAPACKlarfg_(&n1_,y+tr2->idx[1],y+tr2->idx[1]+1,&inc,tr2->tau+1));
    }
    if (tr2->idx[0]<tr2->idx[1]) {
      ierr = HRGen(PetscRealPart(y[tr2->idx[0]]),PetscRealPart(y[tr2->idx[1]]),&(tr2->type),&(tr2->cs),&(tr2->sn),&(tr2->alpha),&ncond2);CHKERRQ(ierr);
    } else {
      tr2->alpha = PetscRealPart(y[tr2->idx[0]]);
      ncond2 = 1.0;
    }
    if (ncond2>*ncond) *ncond = ncond2;
  }
  PetscFunctionReturn(0);
}

/*
  Auxiliary function to try perform one iteration of hr routine,
  checking condition number. If it is < tolD, apply the
  transformation to H and R, if not, ok=false and it do nothing
  tolE, tolerance to exchange complex pairs to improve conditioning
*/
static PetscErrorCode TryHRIt(PetscInt n,PetscInt j,PetscInt sz,PetscScalar *H,PetscInt ldh,PetscScalar *R,PetscInt ldr,PetscReal *s,PetscBool *exg,PetscBool *ok,PetscInt *n0,PetscInt *n1,PetscInt *idx0,PetscInt *idx1,PetscReal *cond,PetscScalar *work)
{
  PetscErrorCode ierr;
  struct HRtr    *tr1,*tr2,tr1_t,tr2_t,tr1_te,tr2_te;
  PetscScalar    *x,*y;
  PetscReal      ncond,ncond_e;
  PetscInt       nwu=0,i,d=1;
  PetscBLASInt   n0_,n1_,inc=1,mh,mr,n_,ldr_,ldh_;
  PetscReal      tolD = 1e+5;

  PetscFunctionBegin;
  if (cond) *cond = 1.0;
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldr,&ldr_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(ldh,&ldh_);CHKERRQ(ierr);
  x = work+nwu;
  nwu += n;
  ierr = PetscArraycpy(x,R+j*ldr,n);CHKERRQ(ierr);
  *exg = PETSC_FALSE;
  *ok = PETSC_TRUE;
  tr1_t.data = x;
  if (sz==1) {
    /* Hyperbolic transformation to make zeros in x */
    ierr = MadeHRtr(sz,n,*idx0,*n0,*idx1,*n1,&tr1_t,NULL,&ncond,work+nwu);CHKERRQ(ierr);
    /* Check condition number to single column*/
    if (ncond>tolD) *ok = PETSC_FALSE;
    tr1 = &tr1_t;
    tr2 = &tr2_t;
  } else {
    y = work+nwu;
    nwu += n;
    ierr = PetscArraycpy(y,R+(j+1)*ldr,n);CHKERRQ(ierr);
    tr2_t.data = y;
    ierr = MadeHRtr(sz,n,*idx0,*n0,*idx1,*n1,&tr1_t,&tr2_t,&ncond,work+nwu);CHKERRQ(ierr);
    /* Computing hyperbolic transformations also for exchanged vectors */
    tr1_te.data = work+nwu;
    nwu += n;
    ierr = PetscArraycpy(tr1_te.data,R+(j+1)*ldr,n);CHKERRQ(ierr);
    tr2_te.data = work+nwu;
    nwu += n;
    ierr = PetscArraycpy(tr2_te.data,R+j*ldr,n);CHKERRQ(ierr);
    ierr = MadeHRtr(sz,n,*idx0,*n0,*idx1,*n1,&tr1_te,&tr2_te,&ncond_e,work+nwu);CHKERRQ(ierr);
    if (ncond > d*ncond_e) {
      *exg = PETSC_TRUE;
      tr1 = &tr1_te;
      tr2 = &tr2_te;
      ncond = ncond_e;
    } else {
      tr1 = &tr1_t;
      tr2 = &tr2_t;
    }
    if (ncond>tolD) *ok = PETSC_FALSE;
  }
  if (*ok) {
    /* Everything is OK, apply transformations to R and H */
    /* First column */
    if (cond && *cond<ncond) *cond = ncond;
    x = tr1->data;
    ierr = PetscBLASIntCast(tr1->n[0],&n0_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(tr1->n[1],&n1_);CHKERRQ(ierr);
    ierr = PetscBLASIntCast(n-j-sz,&mr);CHKERRQ(ierr);
    if (tr1->n[0] > 1 && PetscAbsScalar(tr1->tau[0])!=0.0) {
      x[tr1->idx[0]] = 1.0;
      PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n0_,&mr,x+tr1->idx[0],&inc,tr1->tau,R+(j+sz)*ldr+tr1->idx[0],&ldr_,work+nwu));
      PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("R",&n_,&n0_,x+tr1->idx[0],&inc,tr1->tau,H+(tr1->idx[0])*ldh,&ldh_,work+nwu));
    }
    if (tr1->n[1] > 1 && PetscAbsScalar(tr1->tau[1])!=0.0) {
      x[tr1->idx[1]] = 1.0;
      PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n1_,&mr,x+tr1->idx[1],&inc,tr1->tau+1,R+(j+sz)*ldr+tr1->idx[1],&ldr_,work+nwu));
      PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("R",&n_,&n1_,x+tr1->idx[1],&inc,tr1->tau+1,H+(tr1->idx[1])*ldh,&ldh_,work+nwu));
    }
    if (tr1->idx[0]<tr1->idx[1]) {
      ierr = HRApply(mr,R+(j+sz)*ldr+tr1->idx[0],ldr,R+(j+sz)*ldr+tr1->idx[1],ldr,tr1->cs,-tr1->sn);CHKERRQ(ierr);
      if (tr1->type==1) {
        ierr = HRApply(n,H+(tr1->idx[0])*ldh,1,H+(tr1->idx[1])*ldh,1,tr1->cs,tr1->sn);CHKERRQ(ierr);
      } else {
        ierr = HRApply(n,H+(tr1->idx[0])*ldh,1,H+(tr1->idx[1])*ldh,1,-tr1->cs,-tr1->sn);CHKERRQ(ierr);
        s[tr1->idx[0]] = -s[tr1->idx[0]];
        s[tr1->idx[1]] = -s[tr1->idx[1]];
      }
    }
    for (i=0;i<tr1->idx[0];i++) *(R+j*ldr+i) = x[i];
    for (i=tr1->idx[0]+1;i<n;i++) *(R+j*ldr+i) = 0.0;
    *(R+j*ldr+tr1->idx[0]) = tr1->alpha;
    if (sz==2) {
      y = tr2->data;
      /* Second column */
      ierr = PetscBLASIntCast(tr2->n[0],&n0_);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(tr2->n[1],&n1_);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(n-j-sz,&mr);CHKERRQ(ierr);
      ierr = PetscBLASIntCast(n-tr2->idx[0],&mh);CHKERRQ(ierr);
      if (tr2->n[0] > 1 && PetscAbsScalar(tr2->tau[0])!=0.0) {
        y[tr2->idx[0]] = 1.0;
        PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n0_,&mr,y+tr2->idx[0],&inc,tr2->tau,R+(j+2)*ldr+tr2->idx[0],&ldr_,work+nwu));
        PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("R",&n_,&n0_,y+tr2->idx[0],&inc,tr2->tau,H+(tr2->idx[0])*ldh,&ldh_,work+nwu));
      }
      if (tr2->n[1] > 1 && PetscAbsScalar(tr2->tau[1])!=0.0) {
        y[tr2->idx[1]] = 1.0;
        PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("L",&n1_,&mr,y+tr2->idx[1],&inc,tr2->tau+1,R+(j+2)*ldr+tr2->idx[1],&ldr_,work+nwu));
        PetscStackCallBLAS("LAPACKlarf",LAPACKlarf_("R",&n_,&n1_,y+tr2->idx[1],&inc,tr2->tau+1,H+(tr2->idx[1])*ldh,&ldh_,work+nwu));
      }
      if (tr2->idx[0]<tr2->idx[1]) {
        ierr = HRApply(mr,R+(j+2)*ldr+tr2->idx[0],ldr,R+(j+2)*ldr+tr2->idx[1],ldr,tr2->cs,-tr2->sn);CHKERRQ(ierr);
        if (tr2->type==1) {
          ierr = HRApply(n,H+(tr2->idx[0])*ldh,1,H+(tr2->idx[1])*ldh,1,tr2->cs,tr2->sn);CHKERRQ(ierr);
        } else {
          ierr = HRApply(n,H+(tr2->idx[0])*ldh,1,H+(tr2->idx[1])*ldh,1,-tr2->cs,-tr2->sn);CHKERRQ(ierr);
          s[tr2->idx[0]] = -s[tr2->idx[0]];
          s[tr2->idx[1]] = -s[tr2->idx[1]];
        }
      }
      for (i=0;i<tr2->idx[0]-1;i++) *(R+(j+1)*ldr+i) = y[i];
      *(R+(j+1)*ldr+tr2->idx[0]-1) = y[tr2->idx[0]-1];
      for (i=tr2->idx[0]+1;i<n;i++) *(R+(j+1)*ldr+i) = 0.0;
      *(R+(j+1)*ldr+tr2->idx[0]) = tr2->alpha;
      *n0 = tr2->n[0];
      *n1 = tr2->n[1];
      *idx0 = tr2->idx[0];
      *idx1 = tr2->idx[1];
      if (tr2->idx[0]<tr2->idx[1] && tr2->type==2) {
        (*idx1)++; (*n1)--; (*n0)++;
      }
    } else {
      *n0 = tr1->n[0];
      *n1 = tr1->n[1];
      *idx0 = tr1->idx[0];
      *idx1 = tr1->idx[1];
      if (tr1->idx[0]<tr1->idx[1] && tr1->type==2) {
        (*idx1)++; (*n1)--; (*n0)++;
      }
    }
    if (*n0>0) {
      (*n0)--; (*idx0)++;
      if (*n1==0) *idx1 = *idx0;
    } else {
      (*n1)--; (*idx1)++; *idx0 = *idx1;
    }
  }
  PetscFunctionReturn(0);
}

/*
  compute V = HR whit H s-orthogonal and R upper triangular
*/
static PetscErrorCode PseudoOrthog_HR(PetscInt *nv,PetscScalar *V,PetscInt ldv,PetscReal *s,PetscScalar *R,PetscInt ldr,PetscBLASInt *perm,PetscBLASInt *cmplxEig,PetscBool *breakdown,PetscScalar *work)
{
  PetscErrorCode ierr;
  PetscInt       i,j,n,n0,n1,np,idx0,idx1,sz=1,k=0,t1,t2,nwu=0;
  PetscScalar    *col1,*col2;
  PetscBool      exg=PETSC_FALSE,ok=PETSC_FALSE;

  PetscFunctionBegin;
  n = *nv;
  col1 = work+nwu;
  nwu += n;
  col2 = work+nwu;
  nwu += n;
  /* Sort R and s according to sing(s) */
  np = 0;
  for (i=0;i<n;i++) if (s[i]>0) np++;
  if (s[0]>0) n1 = np;
  else n1 = n-np;
  n0 = 0;
  for (i=0;i<n;i++) {
    if (s[i]==s[0]) {
      s[n0] = s[0];
      perm[n0++] = i;
    } else perm[n1++] = i;
  }
  for (i=n0;i<n;i++) s[i] = -s[0];
  n1 -= n0;
  idx0 = 0;
  idx1 = n0;
  if (idx1==n) idx1=idx0;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) R[j*ldr+i] = V[j*ldv+perm[i]];
  }
  /* Initialize H */
  for (i=0;i<n;i++) {
    ierr = PetscArrayzero(V+i*ldv,n);CHKERRQ(ierr);
    V[perm[i]+i*ldv] = 1.0;
  }
  for (i=0;i<n;i++) perm[i] = i;
  j = 0;
  while (j<n-k) {
    if (cmplxEig[j]==0) sz=1;
    else sz=2;
    ierr = TryHRIt(n,j,sz,V,ldv,R,ldr,s,&exg,&ok,&n0,&n1,&idx0,&idx1,NULL,work+nwu);CHKERRQ(ierr);
    if (ok) {
      if (exg) cmplxEig[j] = -cmplxEig[j];
      j = j+sz;
    } else { /* to be discarded */
      k = k+1;
      if (cmplxEig[j]==0) {
        if (j<n) {
          t1 = perm[j];
          for (i=j;i<n-1;i++) perm[i] = perm[i+1];
          perm[n-1] = t1;
          t1 = cmplxEig[j];
          for (i=j;i<n-1;i++) cmplxEig[i] = cmplxEig[i+1];
          cmplxEig[n-1] = t1;
          ierr = PetscArraycpy(col1,R+j*ldr,n*sizeof(PetscScalar));CHKERRQ(ierr);
          for (i=j;i<n-1;i++) {
            ierr = PetscArraycpy(R+i*ldr,R+(i+1)*ldr,n*sizeof(PetscScalar));CHKERRQ(ierr);
          }
          ierr = PetscArraycpy(R+(n-1)*ldr,col1,n*sizeof(PetscScalar));CHKERRQ(ierr);
        }
      } else {
        k = k+1;
        if (j<n-1) {
          t1 = perm[j]; t2 = perm[j+1];
          for (i=j;i<n-2;i++) perm[i] = perm[i+2];
          perm[n-2] = t1; perm[n-1] = t2;
          t1 = cmplxEig[j]; t2 = cmplxEig[j+1];
          for (i=j;i<n-2;i++) cmplxEig[i] = cmplxEig[i+2];
          cmplxEig[n-2] = t1; cmplxEig[n-1] = t2;
          ierr = PetscArraycpy(col1,R+j*ldr,n);CHKERRQ(ierr);
          ierr = PetscArraycpy(col2,R+(j+1)*ldr,n);CHKERRQ(ierr);
          for (i=j;i<n-2;i++) {
            ierr = PetscArraycpy(R+i*ldr,R+(i+2)*ldr,n);CHKERRQ(ierr);
          }
          ierr = PetscArraycpy(R+(n-2)*ldr,col1,n);CHKERRQ(ierr);
          ierr = PetscArraycpy(R+(n-1)*ldr,col2,n);CHKERRQ(ierr);
        }
      }
    }
  }
  if (k!=0) {
    if (breakdown) *breakdown = PETSC_TRUE;
    *nv = n-k;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DSGHIEPOrthogEigenv(DS ds,DSMatType mat,PetscScalar *wr,PetscScalar *wi,PetscBool accum)
{
  PetscErrorCode ierr;
  PetscInt       lws,nwus=0,nwui=0,lwi;
  PetscInt       off,n,nv,ld,i,ldr,l;
  PetscScalar    *W,*X,*R,*ts,zeroS=0.0,oneS=1.0;
  PetscReal      *s,vi,vr,tr,*d,*e;
  PetscBLASInt   ld_,n_,nv_,*perm,*cmplxEig;

  PetscFunctionBegin;
  l = ds->l;
  n = ds->n-l;
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ld = ds->ld;
  ierr = PetscBLASIntCast(ld,&ld_);CHKERRQ(ierr);
  off = l*ld+l;
  s = ds->rmat[DS_MAT_D];
  if (!ds->compact) {
    for (i=l;i<ds->n;i++) s[i] = PetscRealPart(*(ds->mat[DS_MAT_B]+i*ld+i));
  }
  lws = n*n+7*n;
  lwi = 2*n;
  ierr = DSAllocateWork_Private(ds,lws,0,lwi);CHKERRQ(ierr);
  R = ds->work+nwus;
  nwus += n*n;
  ldr = n;
  perm = ds->iwork + nwui;
  nwui += n;
  cmplxEig = ds->iwork+nwui;
  X = ds->mat[mat];
  for (i=0;i<n;i++) {
#if defined(PETSC_USE_COMPLEX)
    vi = PetscImaginaryPart(wr[l+i]);
#else
    vi = PetscRealPart(wi[l+i]);
#endif
    if (vi!=0) {
      cmplxEig[i] = 1;
      cmplxEig[i+1] = 2;
      i++;
    } else cmplxEig[i] = 0;
  }
  nv = n;

  /* Perform HR decomposition */
  /* Hyperbolic rotators */
  ierr = PseudoOrthog_HR(&nv,X+off,ld,s+l,R,ldr,perm,cmplxEig,NULL,ds->work+nwus);CHKERRQ(ierr);
  /* Sort wr,wi perm */
  ts = ds->work+nwus;
  ierr = PetscArraycpy(ts,wr+l,n);CHKERRQ(ierr);
  for (i=0;i<n;i++) wr[i+l] = ts[perm[i]];
#if !defined(PETSC_USE_COMPLEX)
  ierr = PetscArraycpy(ts,wi+l,n);CHKERRQ(ierr);
  for (i=0;i<n;i++) wi[i+l] = ts[perm[i]];
#endif
  /* Projected Matrix */
  ierr = PetscArrayzero(ds->rmat[DS_MAT_T]+2*ld,ld);CHKERRQ(ierr);
  d = ds->rmat[DS_MAT_T];
  e = d+ld;
  for (i=0;i<nv;i++) {
    if (cmplxEig[i]==0) { /* Real */
      d[l+i] = PetscRealPart(wr[l+i]*s[l+i]);
      e[l+i] = 0.0;
    } else {
      vr = PetscRealPart(wr[l+i]);
#if defined(PETSC_USE_COMPLEX)
      vi = PetscImaginaryPart(wr[l+i]);
#else
      vi = PetscRealPart(wi[l+i]);
#endif
      if (cmplxEig[i]==-1) vi = -vi;
      tr = PetscRealPart((R[i+(i+1)*ldr]/R[i+i*ldr]))*vi;
      d[l+i] = (vr-tr)*s[l+i];
      d[l+i+1] = (vr+tr)*s[l+i+1];
      e[l+i] = PetscRealPart(s[l+i]*(R[(i+1)+(i+1)*ldr]/R[i+i*ldr])*vi);
      e[l+i+1] = 0.0;
      i++;
    }
  }
  /* accumulate previous Q */
  if (accum) {
    ierr = PetscBLASIntCast(nv,&nv_);CHKERRQ(ierr);
    ierr = DSAllocateMat_Private(ds,DS_MAT_W);CHKERRQ(ierr);
    W = ds->mat[DS_MAT_W];
    ierr = DSCopyMatrix_Private(ds,DS_MAT_W,DS_MAT_Q);CHKERRQ(ierr);
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&nv_,&n_,&oneS,W+off,&ld_,X+off,&ld_,&zeroS,ds->mat[DS_MAT_Q]+off,&ld_));
  } else {
    ierr = PetscArrayzero(ds->mat[DS_MAT_Q],ld*ld);CHKERRQ(ierr);
    for (i=0;i<ds->l;i++) *(ds->mat[DS_MAT_Q]+i+i*ld) = 1.0;
    for (i=0;i<n;i++) { ierr = PetscArraycpy(ds->mat[DS_MAT_Q]+off+i*ld,X+off+i*ld,n);CHKERRQ(ierr); }
  }
  ds->t = nv+l;
  if (!ds->compact) { ierr = DSSwitchFormat_GHIEP(ds,PETSC_FALSE);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   Reduce to tridiagonal-diagonal pair by means of TridiagDiag_HHR.
*/
PetscErrorCode DSIntermediate_GHIEP(DS ds)
{
  PetscErrorCode ierr;
  PetscInt       i,ld,off;
  PetscInt       nwall,nwallr,nwalli;
  PetscScalar    *A,*B,*Q;
  PetscReal      *d,*e,*s;

  PetscFunctionBegin;
  ld = ds->ld;
  A = ds->mat[DS_MAT_A];
  B = ds->mat[DS_MAT_B];
  Q = ds->mat[DS_MAT_Q];
  d = ds->rmat[DS_MAT_T];
  e = ds->rmat[DS_MAT_T]+ld;
  s = ds->rmat[DS_MAT_D];
  off = ds->l+ds->l*ld;
  ierr = PetscArrayzero(Q,ld*ld);CHKERRQ(ierr);
  nwall = ld*ld+ld;
  nwallr = ld;
  nwalli = ld;
  ierr = DSAllocateWork_Private(ds,nwall,nwallr,nwalli);CHKERRQ(ierr);
  for (i=0;i<ds->n;i++) Q[i+i*ld]=1.0;
  for (i=0;i<ds->n-ds->l;i++) *(ds->perm+i)=i;
  if (ds->compact) {
    if (ds->state < DS_STATE_INTERMEDIATE) {
      ierr = DSSwitchFormat_GHIEP(ds,PETSC_FALSE);CHKERRQ(ierr);
      ierr = TridiagDiag_HHR(ds->k-ds->l+1,A+off,ld,s+ds->l,Q+off,ld,PETSC_TRUE,d+ds->l,e+ds->l,ds->perm,ds->work,ds->rwork,ds->iwork);CHKERRQ(ierr);
      ds->k = ds->l;
      ierr = PetscArrayzero(d+2*ld+ds->l,ds->n-ds->l);CHKERRQ(ierr);
    }
  } else {
    if (ds->state < DS_STATE_INTERMEDIATE) {
      for (i=0;i<ds->n;i++) s[i] = PetscRealPart(B[i+i*ld]);
      ierr = TridiagDiag_HHR(ds->n-ds->l,A+off,ld,s+ds->l,Q+off,ld,PETSC_FALSE,d+ds->l,e+ds->l,ds->perm,ds->work,ds->rwork,ds->iwork);CHKERRQ(ierr);
      ierr = PetscArrayzero(d+2*ld,ds->n);CHKERRQ(ierr);
      ds->k = ds->l;
      ierr = DSSwitchFormat_GHIEP(ds,PETSC_FALSE);CHKERRQ(ierr);
    } else {
      ierr = DSSwitchFormat_GHIEP(ds,PETSC_TRUE);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

