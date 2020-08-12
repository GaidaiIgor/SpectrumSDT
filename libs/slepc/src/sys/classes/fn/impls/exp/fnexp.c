/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Exponential function  exp(x)
*/

#include <slepc/private/fnimpl.h>      /*I "slepcfn.h" I*/
#include <slepcblaslapack.h>

PetscErrorCode FNEvaluateFunction_Exp(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscFunctionBegin;
  *y = PetscExpScalar(x);
  PetscFunctionReturn(0);
}

PetscErrorCode FNEvaluateDerivative_Exp(FN fn,PetscScalar x,PetscScalar *y)
{
  PetscFunctionBegin;
  *y = PetscExpScalar(x);
  PetscFunctionReturn(0);
}

#define MAX_PADE 6
#define SWAP(a,b,t) {t=a;a=b;b=t;}

PetscErrorCode FNEvaluateFunctionMat_Exp_Pade(FN fn,Mat A,Mat B)
{
  PetscErrorCode ierr;
  PetscBLASInt   n,ld,ld2,*ipiv,info,inc=1;
  PetscInt       m,j,k,sexp;
  PetscBool      odd;
  const PetscInt p=MAX_PADE;
  PetscReal      c[MAX_PADE+1],s,*rwork;
  PetscScalar    scale,mone=-1.0,one=1.0,two=2.0,zero=0.0;
  PetscScalar    *Aa,*Ba,*As,*A2,*Q,*P,*W,*aux;

  PetscFunctionBegin;
  ierr = MatDenseGetArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseGetArray(B,&Ba);CHKERRQ(ierr);
  ierr = MatGetSize(A,&m,NULL);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(m,&n);CHKERRQ(ierr);
  ld  = n;
  ld2 = ld*ld;
  P   = Ba;
  ierr = PetscMalloc6(m*m,&Q,m*m,&W,m*m,&As,m*m,&A2,ld,&rwork,ld,&ipiv);CHKERRQ(ierr);
  ierr = PetscArraycpy(As,Aa,ld2);CHKERRQ(ierr);

  /* Pade' coefficients */
  c[0] = 1.0;
  for (k=1;k<=p;k++) c[k] = c[k-1]*(p+1-k)/(k*(2*p+1-k));

  /* Scaling */
  s = LAPACKlange_("I",&n,&n,As,&ld,rwork);
  ierr = PetscLogFlops(1.0*n*n);CHKERRQ(ierr);
  if (s>0.5) {
    sexp = PetscMax(0,(int)(PetscLogReal(s)/PetscLogReal(2.0))+2);
    scale = PetscPowRealInt(2.0,-sexp);
    PetscStackCallBLAS("BLASscal",BLASscal_(&ld2,&scale,As,&inc));
    ierr = PetscLogFlops(1.0*n*n);CHKERRQ(ierr);
  } else sexp = 0;

  /* Horner evaluation */
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&one,As,&ld,As,&ld,&zero,A2,&ld));
  ierr = PetscLogFlops(2.0*n*n*n);CHKERRQ(ierr);
  ierr = PetscArrayzero(Q,ld2);CHKERRQ(ierr);
  ierr = PetscArrayzero(P,ld2);CHKERRQ(ierr);
  for (j=0;j<n;j++) {
    Q[j+j*ld] = c[p];
    P[j+j*ld] = c[p-1];
  }

  odd = PETSC_TRUE;
  for (k=p-1;k>0;k--) {
    if (odd) {
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&one,Q,&ld,A2,&ld,&zero,W,&ld));
      SWAP(Q,W,aux);
      for (j=0;j<n;j++) Q[j+j*ld] += c[k-1];
      odd = PETSC_FALSE;
    } else {
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&one,P,&ld,A2,&ld,&zero,W,&ld));
      SWAP(P,W,aux);
      for (j=0;j<n;j++) P[j+j*ld] += c[k-1];
      odd = PETSC_TRUE;
    }
    ierr = PetscLogFlops(2.0*n*n*n);CHKERRQ(ierr);
  }
  /*if (odd) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&one,Q,&ld,As,&ld,&zero,W,&ld));
    SWAP(Q,W,aux);
    PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&ld2,&mone,P,&inc,Q,&inc));
    PetscStackCallBLAS("LAPACKgesv",LAPACKgesv_(&n,&n,Q,&ld,ipiv,P,&ld,&info));
    SlepcCheckLapackInfo("gesv",info);
    PetscStackCallBLAS("BLASscal",BLASscal_(&ld2,&two,P,&inc));
    for (j=0;j<n;j++) P[j+j*ld] += 1.0;
    PetscStackCallBLAS("BLASscal",BLASscal_(&ld2,&mone,P,&inc));
  } else {*/
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&one,P,&ld,As,&ld,&zero,W,&ld));
    SWAP(P,W,aux);
    PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&ld2,&mone,P,&inc,Q,&inc));
    PetscStackCallBLAS("LAPACKgesv",LAPACKgesv_(&n,&n,Q,&ld,ipiv,P,&ld,&info));
    SlepcCheckLapackInfo("gesv",info);
    PetscStackCallBLAS("BLASscal",BLASscal_(&ld2,&two,P,&inc));
    for (j=0;j<n;j++) P[j+j*ld] += 1.0;
  /*}*/
  ierr = PetscLogFlops(2.0*n*n*n+2.0*n*n*n/3.0+4.0*n*n);CHKERRQ(ierr);

  for (k=1;k<=sexp;k++) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&one,P,&ld,P,&ld,&zero,W,&ld));
    ierr = PetscArraycpy(P,W,ld2);CHKERRQ(ierr);
  }
  if (P!=Ba) { ierr = PetscArraycpy(Ba,P,ld2);CHKERRQ(ierr); }
  ierr = PetscLogFlops(2.0*n*n*n*sexp);CHKERRQ(ierr);

  ierr = PetscFree6(Q,W,As,A2,rwork,ipiv);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(B,&Ba);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
 * Set scaling factor (s) and Pade degree (k,m)
 */
static PetscErrorCode sexpm_params(PetscReal nrm,PetscInt *s,PetscInt *k,PetscInt *m)
{
  PetscFunctionBegin;
  if (nrm>1) {
    if      (nrm<200)  {*s = 4; *k = 5; *m = *k-1;}
    else if (nrm<1e4)  {*s = 4; *k = 4; *m = *k+1;}
    else if (nrm<1e6)  {*s = 4; *k = 3; *m = *k+1;}
    else if (nrm<1e9)  {*s = 3; *k = 3; *m = *k+1;}
    else if (nrm<1e11) {*s = 2; *k = 3; *m = *k+1;}
    else if (nrm<1e12) {*s = 2; *k = 2; *m = *k+1;}
    else if (nrm<1e14) {*s = 2; *k = 1; *m = *k+1;}
    else               {*s = 1; *k = 1; *m = *k+1;}
  } else { /* nrm<1 */
    if       (nrm>0.5)  {*s = 4; *k = 4; *m = *k-1;}
    else  if (nrm>0.3)  {*s = 3; *k = 4; *m = *k-1;}
    else  if (nrm>0.15) {*s = 2; *k = 4; *m = *k-1;}
    else  if (nrm>0.07) {*s = 1; *k = 4; *m = *k-1;}
    else  if (nrm>0.01) {*s = 0; *k = 4; *m = *k-1;}
    else  if (nrm>3e-4) {*s = 0; *k = 3; *m = *k-1;}
    else  if (nrm>1e-5) {*s = 0; *k = 3; *m = 0;}
    else  if (nrm>1e-8) {*s = 0; *k = 2; *m = 0;}
    else                {*s = 0; *k = 1; *m = 0;}
  }
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_COMPLEX)
/*
 * Partial fraction form coefficients.
 * If query, the function returns the size necessary to store the coefficients.
 */
static PetscErrorCode getcoeffs(PetscInt k,PetscInt m,PetscComplex *r,PetscComplex *q,PetscComplex *remain,PetscBool query)
{
  PetscInt i;
  const PetscComplex /* m == k+1 */
    p1r4[5] = {-1.582680186458572e+01 - 2.412564578224361e+01*PETSC_i,
               -1.582680186458572e+01 + 2.412564578224361e+01*PETSC_i,
                1.499984465975511e+02 + 6.804227952202417e+01*PETSC_i,
                1.499984465975511e+02 - 6.804227952202417e+01*PETSC_i,
               -2.733432894659307e+02                                },
    p1q4[5] = { 3.655694325463550e+00 + 6.543736899360086e+00*PETSC_i,
                3.655694325463550e+00 - 6.543736899360086e+00*PETSC_i,
                5.700953298671832e+00 + 3.210265600308496e+00*PETSC_i,
                5.700953298671832e+00 - 3.210265600308496e+00*PETSC_i,
                6.286704751729261e+00                               },
    p1r3[4] = {-1.130153999597152e+01 + 1.247167585025031e+01*PETSC_i,
               -1.130153999597152e+01 - 1.247167585025031e+01*PETSC_i,
                1.330153999597152e+01 - 6.007173273704750e+01*PETSC_i,
                1.330153999597152e+01 + 6.007173273704750e+01*PETSC_i},
    p1q3[4] = { 3.212806896871536e+00 + 4.773087433276636e+00*PETSC_i,
                3.212806896871536e+00 - 4.773087433276636e+00*PETSC_i,
                4.787193103128464e+00 + 1.567476416895212e+00*PETSC_i,
                4.787193103128464e+00 - 1.567476416895212e+00*PETSC_i},
    p1r2[3] = { 7.648749087422928e+00 + 4.171640244747463e+00*PETSC_i,
                7.648749087422928e+00 - 4.171640244747463e+00*PETSC_i,
               -1.829749817484586e+01                                },
    p1q2[3] = { 2.681082873627756e+00 + 3.050430199247411e+00*PETSC_i,
                2.681082873627756e+00 - 3.050430199247411e+00*PETSC_i,
                3.637834252744491e+00                                },
    p1r1[2] = { 1.000000000000000e+00 - 3.535533905932738e+00*PETSC_i,
                1.000000000000000e+00 + 3.535533905932738e+00*PETSC_i},
    p1q1[2] = { 2.000000000000000e+00 + 1.414213562373095e+00*PETSC_i,
                2.000000000000000e+00 - 1.414213562373095e+00*PETSC_i};
  const PetscComplex /* m == k-1 */
    m1r5[4] = {-1.423367961376821e+02 - 1.385465094833037e+01*PETSC_i,
               -1.423367961376821e+02 + 1.385465094833037e+01*PETSC_i,
                2.647367961376822e+02 - 4.814394493714596e+02*PETSC_i,
                2.647367961376822e+02 + 4.814394493714596e+02*PETSC_i},
    m1q5[4] = { 5.203941240131764e+00 + 5.805856841805367e+00*PETSC_i,
                5.203941240131764e+00 - 5.805856841805367e+00*PETSC_i,
                6.796058759868242e+00 + 1.886649260140217e+00*PETSC_i,
                6.796058759868242e+00 - 1.886649260140217e+00*PETSC_i},
    m1r4[3] = { 2.484269593165883e+01 + 7.460342395992306e+01*PETSC_i,
                2.484269593165883e+01 - 7.460342395992306e+01*PETSC_i,
               -1.734353918633177e+02                                },
    m1q4[3] = { 4.675757014491557e+00 + 3.913489560603711e+00*PETSC_i,
                4.675757014491557e+00 - 3.913489560603711e+00*PETSC_i,
                5.648485971016893e+00                                },
    m1r3[2] = { 2.533333333333333e+01 - 2.733333333333333e+01*PETSC_i,
                2.533333333333333e+01 + 2.733333333333333e+01*PETSC_i},
    m1q3[2] = { 4.000000000000000e+00 + 2.000000000000000e+00*PETSC_i,
                4.000000000000000e+00 - 2.000000000000000e+00*PETSC_i};
  const PetscScalar /* m == k-1 */
    m1remain5[2] = { 2.000000000000000e-01,  9.800000000000000e+00},
    m1remain4[2] = {-2.500000000000000e-01, -7.750000000000000e+00},
    m1remain3[2] = { 3.333333333333333e-01,  5.666666666666667e+00},
    m1remain2[2] = {-0.5,                   -3.5},
    remain3[4] = {1.0/6.0, 1.0/2.0, 1, 1},
    remain2[3] = {1.0/2.0, 1, 1};

  PetscFunctionBegin;
  if (query) { /* query about buffer's size */
    if (m==k+1) {
      *remain = 0;
      *r = *q = k+1;
      PetscFunctionReturn(0); /* quick return */
    }
    if (m==k-1) {
      *remain = 2;
      if (k==5) *r = *q = 4;
      else if (k==4) *r = *q = 3;
      else if (k==3) *r = *q = 2;
      else if (k==2) *r = *q = 1;
    }
    if (m==0) {
      *r = *q = 0;
      *remain = k+1;
    }
  } else {
    if (m==k+1) {
      if (k==4) {
        for (i=0;i<5;i++) { r[i] = p1r4[i]; q[i] = p1q4[i]; }
      } else if (k==3) {
        for (i=0;i<4;i++) { r[i] = p1r3[i]; q[i] = p1q3[i]; }
      } else if (k==2) {
        for (i=0;i<3;i++) { r[i] = p1r2[i]; q[i] = p1q2[i]; }
      } else if (k==1) {
        for (i=0;i<2;i++) { r[i] = p1r1[i]; q[i] = p1q1[i]; }
      }
      PetscFunctionReturn(0); /* quick return */
    }
    if (m==k-1) {
      if (k==5) {
        for (i=0;i<4;i++) { r[i] = m1r5[i]; q[i] = m1q5[i]; }
        for (i=0;i<2;i++) remain[i] = m1remain5[i];
      } else if (k==4) {
        for (i=0;i<3;i++) { r[i] = m1r4[i]; q[i] = m1q4[i]; }
        for (i=0;i<2;i++) remain[i] = m1remain4[i];
      } else if (k==3) {
        for (i=0;i<2;i++) { r[i] = m1r3[i]; q[i] = m1q3[i]; remain[i] = m1remain3[i]; }
      } else if (k==2) {
        r[0] = -13.5; q[0] = 3;
        for (i=0;i<2;i++) remain[i] = m1remain2[i];
      }
    }
    if (m==0) {
      r = q = 0;
      if (k==3) {
        for (i=0;i<4;i++) remain[i] = remain3[i];
      } else if (k==2) {
        for (i=0;i<3;i++) remain[i] = remain2[i];
      }
    }
  }
  PetscFunctionReturn(0);
}

/*
 * Product form coefficients.
 * If query, the function returns the size necessary to store the coefficients.
 */
static PetscErrorCode getcoeffsproduct(PetscInt k,PetscInt m,PetscComplex *p,PetscComplex *q,PetscComplex *mult,PetscBool query)
{
  PetscInt i;
  const PetscComplex /* m == k+1 */
  p1p4[4] = {-5.203941240131764e+00 + 5.805856841805367e+00*PETSC_i,
             -5.203941240131764e+00 - 5.805856841805367e+00*PETSC_i,
             -6.796058759868242e+00 + 1.886649260140217e+00*PETSC_i,
             -6.796058759868242e+00 - 1.886649260140217e+00*PETSC_i},
  p1q4[5] = { 3.655694325463550e+00 + 6.543736899360086e+00*PETSC_i,
              3.655694325463550e+00 - 6.543736899360086e+00*PETSC_i,
              6.286704751729261e+00                                ,
              5.700953298671832e+00 + 3.210265600308496e+00*PETSC_i,
              5.700953298671832e+00 - 3.210265600308496e+00*PETSC_i},
  p1p3[3] = {-4.675757014491557e+00 + 3.913489560603711e+00*PETSC_i,
             -4.675757014491557e+00 - 3.913489560603711e+00*PETSC_i,
             -5.648485971016893e+00                                },
  p1q3[4] = { 3.212806896871536e+00 + 4.773087433276636e+00*PETSC_i,
              3.212806896871536e+00 - 4.773087433276636e+00*PETSC_i,
              4.787193103128464e+00 + 1.567476416895212e+00*PETSC_i,
              4.787193103128464e+00 - 1.567476416895212e+00*PETSC_i},
  p1p2[2] = {-4.00000000000000e+00  + 2.000000000000000e+00*PETSC_i,
             -4.00000000000000e+00  - 2.000000000000000e+00*PETSC_i},
  p1q2[3] = { 2.681082873627756e+00 + 3.050430199247411e+00*PETSC_i,
              2.681082873627756e+00 - 3.050430199247411e+00*PETSC_i,
              3.637834252744491e+00                               },
  p1q1[2] = { 2.000000000000000e+00 + 1.414213562373095e+00*PETSC_i,
              2.000000000000000e+00 - 1.414213562373095e+00*PETSC_i};
  const PetscComplex /* m == k-1 */
  m1p5[5] = {-3.655694325463550e+00 + 6.543736899360086e+00*PETSC_i,
             -3.655694325463550e+00 - 6.543736899360086e+00*PETSC_i,
             -6.286704751729261e+00                                ,
             -5.700953298671832e+00 + 3.210265600308496e+00*PETSC_i,
             -5.700953298671832e+00 - 3.210265600308496e+00*PETSC_i},
  m1q5[4] = { 5.203941240131764e+00 + 5.805856841805367e+00*PETSC_i,
              5.203941240131764e+00 - 5.805856841805367e+00*PETSC_i,
              6.796058759868242e+00 + 1.886649260140217e+00*PETSC_i,
              6.796058759868242e+00 - 1.886649260140217e+00*PETSC_i},
  m1p4[4] = {-3.212806896871536e+00 + 4.773087433276636e+00*PETSC_i,
             -3.212806896871536e+00 - 4.773087433276636e+00*PETSC_i,
             -4.787193103128464e+00 + 1.567476416895212e+00*PETSC_i,
             -4.787193103128464e+00 - 1.567476416895212e+00*PETSC_i},
  m1q4[3] = { 4.675757014491557e+00 + 3.913489560603711e+00*PETSC_i,
              4.675757014491557e+00 - 3.913489560603711e+00*PETSC_i,
              5.648485971016893e+00                                },
  m1p3[3] = {-2.681082873627756e+00 + 3.050430199247411e+00*PETSC_i,
             -2.681082873627756e+00 - 3.050430199247411e+00*PETSC_i,
             -3.637834252744491e+00                                },
  m1q3[2] = { 4.000000000000000e+00 + 2.000000000000000e+00*PETSC_i,
              4.000000000000000e+00 - 2.000000000000001e+00*PETSC_i},
  m1p2[2] = {-2.000000000000000e+00 + 1.414213562373095e+00*PETSC_i,
             -2.000000000000000e+00 - 1.414213562373095e+00*PETSC_i};

  PetscFunctionBegin;
  if (query) {
    if (m == k+1) {
      *mult = 1;
      *p = k;
      *q = k+1;
      PetscFunctionReturn(0);
    }
    if (m==k-1) {
      *mult = 1;
      *p = k;
      *q = k-1;
    }
  } else {
    if (m == k+1) {
      *mult = PetscPowInt(-1,m);
      *mult *= m;
      if (k==4) {
        for (i=0;i<4;i++) { p[i] = p1p4[i]; q[i] = p1q4[i]; }
        q[4] = p1q4[4];
      } else if (k==3) {
        for (i=0;i<3;i++) { p[i] = p1p3[i]; q[i] = p1q3[i]; }
        q[3] = p1q3[3];
      } else if (k==2) {
        for (i=0;i<2;i++) { p[i] = p1p2[i]; q[i] = p1q2[i]; }
        q[2] = p1q2[2];
      } else if (k==1) {
        p[0] = -3;
        for (i=0;i<2;i++) q[i] = p1q1[i];
      }
      PetscFunctionReturn(0);
    }
    if (m==k-1) {
      *mult = PetscPowInt(-1,m);
      *mult /= k;
      if (k==5) {
        for (i=0;i<4;i++) { p[i] = m1p5[i]; q[i] = m1q5[i]; }
        p[4] = m1p5[4];
      } else if (k==4) {
        for (i=0;i<3;i++) { p[i] = m1p4[i]; q[i] = m1q4[i]; }
        p[3] = m1p4[3];
      } else if (k==3) {
        for (i=0;i<2;i++) { p[i] = m1p3[i]; q[i] = m1q3[i]; }
        p[2] = m1p3[2];
      } else if (k==2) {
        for (i=0;i<2;i++) p[i] = m1p2[i];
        q[0] = 3;
      }
    }
  }
  PetscFunctionReturn(0);
}
#endif /* PETSC_HAVE_COMPLEX */

#if defined(PETSC_USE_COMPLEX)
static PetscErrorCode getisreal(PetscInt n,PetscComplex *a,PetscBool *result)
{
  PetscInt i;

  PetscFunctionBegin;
  *result=PETSC_TRUE;
  for (i=0;i<n&&*result;i++) {
    if (PetscImaginaryPartComplex(a[i])) *result=PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}
#endif

/*
 * Matrix exponential implementation based on algorithm and matlab code by Stefan Guettel
 * and Yuji Nakatsukasa
 *
 *     Stefan Guettel and Yuji Nakatsukasa, "Scaled and Squared Subdiagonal Pade
 *     Approximation for the Matrix Exponential",
 *     SIAM J. Matrix Anal. Appl. 37(1):145-170, 2016.
 *     https://doi.org/10.1137/15M1027553
 */
PetscErrorCode FNEvaluateFunctionMat_Exp_GuettelNakatsukasa(FN fn,Mat A,Mat B)
{
#if !defined(PETSC_HAVE_COMPLEX)
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,1,"This function requires C99 or C++ complex support");
#else
  PetscInt       i,j,n_,s,k,m,mod;
  PetscBLASInt   n,n2,irsize,rsizediv2,ipsize,iremainsize,info,*piv,minlen,lwork,one=1;
  PetscReal      nrm,shift;
#if defined(PETSC_USE_COMPLEX) || defined(PETSC_HAVE_ESSL)
  PetscReal      *rwork=NULL;
#endif
  PetscComplex   *As,*RR,*RR2,*expmA,*expmA2,*Maux,*Maux2,rsize,*r,psize,*p,remainsize,*remainterm,*rootp,*rootq,mult=0.0,scale,cone=1.0,czero=0.0,*aux;
  PetscScalar    *Aa,*Ba,*Ba2,*sMaux,*wr,*wi,expshift,sone=1.0,szero=0.0,*saux;
  PetscErrorCode ierr;
  PetscBool      isreal;
#if defined(PETSC_HAVE_ESSL)
  PetscScalar    sdummy,*wri;
  PetscBLASInt   idummy,io=0;
#else
  PetscBLASInt   query=-1;
  PetscScalar    work1,*work;
#endif

  PetscFunctionBegin;
  ierr = MatGetSize(A,&n_,NULL);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n_,&n);CHKERRQ(ierr);
  ierr = MatDenseGetArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseGetArray(B,&Ba);CHKERRQ(ierr);
  Ba2 = Ba;
  ierr = PetscBLASIntCast(n*n,&n2);CHKERRQ(ierr);

  ierr = PetscMalloc2(n2,&sMaux,n2,&Maux);CHKERRQ(ierr);
  Maux2 = Maux;
  ierr = PetscMalloc2(n,&wr,n,&wi);CHKERRQ(ierr);
  ierr = PetscArraycpy(sMaux,Aa,n2);CHKERRQ(ierr);
  /* estimate rightmost eigenvalue and shift A with it */
#if !defined(PETSC_HAVE_ESSL)
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_("N","N",&n,sMaux,&n,wr,wi,NULL,&n,NULL,&n,&work1,&query,&info));
  SlepcCheckLapackInfo("geev",info);
  ierr = PetscBLASIntCast((PetscInt)PetscRealPart(work1),&lwork);CHKERRQ(ierr);
  ierr = PetscMalloc1(lwork,&work);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_("N","N",&n,sMaux,&n,wr,wi,NULL,&n,NULL,&n,work,&lwork,&info));
  ierr = PetscFree(work);CHKERRQ(ierr);
#else
  ierr = PetscArraycpy(Maux,Aa,n2);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_("N","N",&n,Maux,&n,wr,NULL,&n,NULL,&n,&work1,&query,rwork,&info));
  SlepcCheckLapackInfo("geev",info);
  ierr = PetscBLASIntCast((PetscInt)PetscRealPart(work1),&lwork);CHKERRQ(ierr);
  ierr = PetscMalloc2(2*n,&rwork,lwork,&work);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_("N","N",&n,Maux,&n,wr,NULL,&n,NULL,&n,work,&lwork,rwork,&info));
  ierr = PetscFree2(rwork,work);CHKERRQ(ierr);
#endif
  SlepcCheckLapackInfo("geev",info);
#else /* defined(PETSC_HAVE_ESSL) */
  ierr = PetscBLASIntCast(4*n,&lwork);CHKERRQ(ierr);
  ierr = PetscMalloc2(lwork,&rwork,2*n,&wri);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_(&io,sMaux,&n,wri,&sdummy,&idummy,&idummy,&n,rwork,&lwork));
  for (i=0;i<n;i++) {
    wr[i] = wri[2*i];
    wi[i] = wri[2*i+1];
  }
#else
  PetscStackCallBLAS("LAPACKgeev",LAPACKgeev_(&io,Maux,&n,wri,&sdummy,&idummy,&idummy,&n,rwork,&lwork));
  for (i=0;i<n;i++) wr[i] = wri[i];
#endif
  ierr = PetscFree2(rwork,wri);CHKERRQ(ierr);
#endif
  ierr = PetscLogFlops(25.0*n*n*n+(n*n*n)/3.0+1.0*n*n*n);CHKERRQ(ierr);

  shift = PetscRealPart(wr[0]);
  for (i=1;i<n;i++) {
    if (PetscRealPart(wr[i]) > shift) shift = PetscRealPart(wr[i]);
  }
  ierr = PetscFree2(wr,wi);CHKERRQ(ierr);
  /* shift so that largest real part is (about) 0 */
  ierr = PetscArraycpy(sMaux,Aa,n2);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
    sMaux[i+i*n] -= shift;
  }
  ierr = PetscLogFlops(1.0*n);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr = PetscArraycpy(Maux,Aa,n2);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
    Maux[i+i*n] -= shift;
  }
  ierr = PetscLogFlops(1.0*n);CHKERRQ(ierr);
#endif

  /* estimate norm(A) and select the scaling factor */
  nrm = LAPACKlange_("O",&n,&n,sMaux,&n,NULL);
  ierr = PetscLogFlops(1.0*n*n);CHKERRQ(ierr);
  ierr = sexpm_params(nrm,&s,&k,&m);CHKERRQ(ierr);
  if (s==0 && k==1 && m==0) { /* exp(A) = I+A to eps! */
    expshift = PetscExpReal(shift);
    for (i=0;i<n;i++) sMaux[i+i*n] += 1.0;
    PetscStackCallBLAS("BLASscal",BLASscal_(&n2,&expshift,sMaux,&one));
    ierr = PetscLogFlops(1.0*(n+n2));CHKERRQ(ierr);
    ierr = PetscArraycpy(Ba,sMaux,n2);CHKERRQ(ierr);
    ierr = PetscFree2(sMaux,Maux);CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(A,&Aa);CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(B,&Ba);CHKERRQ(ierr);
    PetscFunctionReturn(0); /* quick return */
  }

  ierr = PetscMalloc4(n2,&expmA,n2,&As,n2,&RR,n,&piv);CHKERRQ(ierr);
  expmA2 = expmA; RR2 = RR;
  /* scale matrix */
#if !defined(PETSC_USE_COMPLEX)
  for (i=0;i<n2;i++) {
    As[i] = sMaux[i];
  }
#else
  ierr = PetscArraycpy(As,sMaux,n2);CHKERRQ(ierr);
#endif
  scale = 1.0/PetscPowRealInt(2.0,s);
  PetscStackCallBLAS("BLASCOMPLEXscal",BLASCOMPLEXscal_(&n2,&scale,As,&one));
  ierr = SlepcLogFlopsComplex(1.0*n2);CHKERRQ(ierr);

  /* evaluate Pade approximant (partial fraction or product form) */
  if (fn->method==3 || !m) { /* partial fraction */
    ierr = getcoeffs(k,m,&rsize,&psize,&remainsize,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscBLASIntCast((PetscInt)PetscRealPartComplex(rsize),&irsize);CHKERRQ(ierr);
    ierr = PetscBLASIntCast((PetscInt)PetscRealPartComplex(psize),&ipsize);CHKERRQ(ierr);
    ierr = PetscBLASIntCast((PetscInt)PetscRealPartComplex(remainsize),&iremainsize);CHKERRQ(ierr);
    ierr = PetscMalloc3(irsize,&r,ipsize,&p,iremainsize,&remainterm);CHKERRQ(ierr);
    ierr = getcoeffs(k,m,r,p,remainterm,PETSC_FALSE);CHKERRQ(ierr);

    ierr = PetscArrayzero(expmA,n2);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
    isreal = PETSC_TRUE;
#else
    ierr = getisreal(n2,Maux,&isreal);CHKERRQ(ierr);
#endif
    if (isreal) {
      rsizediv2 = irsize/2;
      for (i=0;i<rsizediv2;i++) { /* use partial fraction to get R(As) */
        ierr = PetscArraycpy(Maux,As,n2);CHKERRQ(ierr);
        ierr = PetscArrayzero(RR,n2);CHKERRQ(ierr);
        for (j=0;j<n;j++) {
          Maux[j+j*n] -= p[2*i];
          RR[j+j*n] = r[2*i];
        }
        PetscStackCallBLAS("LAPACKCOMPLEXgesv",LAPACKCOMPLEXgesv_(&n,&n,Maux,&n,piv,RR,&n,&info));
        SlepcCheckLapackInfo("gesv",info);
        for (j=0;j<n2;j++) {
          expmA[j] += RR[j] + PetscConj(RR[j]);
        }
        /* loop(n) + gesv + loop(n2) */
        ierr = SlepcLogFlopsComplex(1.0*n+(2.0*n*n*n/3.0+2.0*n*n*n)+2.0*n2);CHKERRQ(ierr);
      }

      mod = ipsize % 2;
      if (mod) {
        ierr = PetscArraycpy(Maux,As,n2);CHKERRQ(ierr);
        ierr = PetscArrayzero(RR,n2);CHKERRQ(ierr);
        for (j=0;j<n;j++) {
          Maux[j+j*n] -= p[ipsize-1];
          RR[j+j*n] = r[irsize-1];
        }
        PetscStackCallBLAS("LAPACKCOMPLEXgesv",LAPACKCOMPLEXgesv_(&n,&n,Maux,&n,piv,RR,&n,&info));
        SlepcCheckLapackInfo("gesv",info);
        for (j=0;j<n2;j++) {
          expmA[j] += RR[j];
        }
        ierr = SlepcLogFlopsComplex(1.0*n+(2.0*n*n*n/3.0+2.0*n*n*n)+1.0*n2);CHKERRQ(ierr);
      }
    } else { /* complex */
      for (i=0;i<irsize;i++) { /* use partial fraction to get R(As) */
        ierr = PetscArraycpy(Maux,As,n2);CHKERRQ(ierr);
        ierr = PetscArrayzero(RR,n2);CHKERRQ(ierr);
        for (j=0;j<n;j++) {
          Maux[j+j*n] -= p[i];
          RR[j+j*n] = r[i];
        }
        PetscStackCallBLAS("LAPACKCOMPLEXgesv",LAPACKCOMPLEXgesv_(&n,&n,Maux,&n,piv,RR,&n,&info));
        SlepcCheckLapackInfo("gesv",info);
        for (j=0;j<n2;j++) {
          expmA[j] += RR[j];
        }
        ierr = SlepcLogFlopsComplex(1.0*n+(2.0*n*n*n/3.0+2.0*n*n*n)+1.0*n2);CHKERRQ(ierr);
      }
    }
    for (i=0;i<iremainsize;i++) {
      if (!i) {
        ierr = PetscArrayzero(RR,n2);CHKERRQ(ierr);
        for (j=0;j<n;j++) {
          RR[j+j*n] = remainterm[iremainsize-1];
        }
      } else {
        ierr = PetscArraycpy(RR,As,n2);CHKERRQ(ierr);
        for (j=1;j<i;j++) {
          PetscStackCallBLAS("BLASCOMPLEXgemm",BLASCOMPLEXgemm_("N","N",&n,&n,&n,&cone,RR,&n,RR,&n,&czero,Maux,&n));
          SWAP(RR,Maux,aux);
          ierr = SlepcLogFlopsComplex(2.0*n*n*n);CHKERRQ(ierr);
        }
        PetscStackCallBLAS("BLASCOMPLEXscal",BLASCOMPLEXscal_(&n2,&remainterm[iremainsize-1-i],RR,&one));
        ierr = SlepcLogFlopsComplex(1.0*n2);CHKERRQ(ierr);
      }
      for (j=0;j<n2;j++) {
        expmA[j] += RR[j];
      }
      ierr = SlepcLogFlopsComplex(1.0*n2);CHKERRQ(ierr);
    }
    ierr = PetscFree3(r,p,remainterm);CHKERRQ(ierr);
  } else { /* product form, default */
    ierr = getcoeffsproduct(k,m,&rsize,&psize,&mult,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscBLASIntCast((PetscInt)PetscRealPartComplex(rsize),&irsize);CHKERRQ(ierr);
    ierr = PetscBLASIntCast((PetscInt)PetscRealPartComplex(psize),&ipsize);CHKERRQ(ierr);
    ierr = PetscMalloc2(irsize,&rootp,ipsize,&rootq);CHKERRQ(ierr);
    ierr = getcoeffsproduct(k,m,rootp,rootq,&mult,PETSC_FALSE);CHKERRQ(ierr);

    ierr = PetscArrayzero(expmA,n2);CHKERRQ(ierr);
    for (i=0;i<n;i++) { /* initialize */
      expmA[i+i*n] = 1.0;
    }
    minlen = PetscMin(irsize,ipsize);
    for (i=0;i<minlen;i++) {
      ierr = PetscArraycpy(RR,As,n2);CHKERRQ(ierr);
      for (j=0;j<n;j++) {
        RR[j+j*n] -= rootp[i];
      }
      PetscStackCallBLAS("BLASCOMPLEXgemm",BLASCOMPLEXgemm_("N","N",&n,&n,&n,&cone,RR,&n,expmA,&n,&czero,Maux,&n));
      SWAP(expmA,Maux,aux);
      ierr = PetscArraycpy(RR,As,n2);CHKERRQ(ierr);
      for (j=0;j<n;j++) {
        RR[j+j*n] -= rootq[i];
      }
      PetscStackCallBLAS("LAPACKCOMPLEXgesv",LAPACKCOMPLEXgesv_(&n,&n,RR,&n,piv,expmA,&n,&info));
      SlepcCheckLapackInfo("gesv",info);
      /* loop(n) + gemm + loop(n) + gesv */
      ierr = SlepcLogFlopsComplex(1.0*n+(2.0*n*n*n)+1.0*n+(2.0*n*n*n/3.0+2.0*n*n*n));CHKERRQ(ierr);
    }
    /* extra numerator */
    for (i=minlen;i<irsize;i++) {
      ierr = PetscArraycpy(RR,As,n2);CHKERRQ(ierr);
      for (j=0;j<n;j++) {
        RR[j+j*n] -= rootp[i];
      }
      PetscStackCallBLAS("BLASCOMPLEXgemm",BLASCOMPLEXgemm_("N","N",&n,&n,&n,&cone,RR,&n,expmA,&n,&czero,Maux,&n));
      SWAP(expmA,Maux,aux);
      ierr = SlepcLogFlopsComplex(1.0*n+2.0*n*n*n);CHKERRQ(ierr);
    }
    /* extra denominator */
    for (i=minlen;i<ipsize;i++) {
      ierr = PetscArraycpy(RR,As,n2);CHKERRQ(ierr);
      for (j=0;j<n;j++) RR[j+j*n] -= rootq[i];
      PetscStackCallBLAS("LAPACKCOMPLEXgesv",LAPACKCOMPLEXgesv_(&n,&n,RR,&n,piv,expmA,&n,&info));
      SlepcCheckLapackInfo("gesv",info);
      ierr = SlepcLogFlopsComplex(1.0*n+(2.0*n*n*n/3.0+2.0*n*n*n));CHKERRQ(ierr);
    }
    PetscStackCallBLAS("BLASCOMPLEXscal",BLASCOMPLEXscal_(&n2,&mult,expmA,&one));
    ierr = SlepcLogFlopsComplex(1.0*n2);CHKERRQ(ierr);
    ierr = PetscFree2(rootp,rootq);CHKERRQ(ierr);
  }

#if !defined(PETSC_USE_COMPLEX)
  for (i=0;i<n2;i++) {
    Ba2[i] = PetscRealPartComplex(expmA[i]);
  }
#else
  ierr = PetscArraycpy(Ba2,expmA,n2);CHKERRQ(ierr);
#endif

  /* perform repeated squaring */
  for (i=0;i<s;i++) { /* final squaring */
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n,&n,&n,&sone,Ba2,&n,Ba2,&n,&szero,sMaux,&n));
    SWAP(Ba2,sMaux,saux);
    ierr = PetscLogFlops(2.0*n*n*n);CHKERRQ(ierr);
  }
  if (Ba2!=Ba) {
    ierr = PetscArraycpy(Ba,Ba2,n2);CHKERRQ(ierr);
    sMaux = Ba2;
  }
  expshift = PetscExpReal(shift);
  PetscStackCallBLAS("BLASscal",BLASscal_(&n2,&expshift,Ba,&one));
  ierr = PetscLogFlops(1.0*n2);CHKERRQ(ierr);

  /* restore pointers */
  Maux = Maux2; expmA = expmA2; RR = RR2;
  ierr = PetscFree2(sMaux,Maux);CHKERRQ(ierr);
  ierr = PetscFree4(expmA,As,RR,piv);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(B,&Ba);CHKERRQ(ierr);
  PetscFunctionReturn(0);
#endif
}

#define SMALLN 100

/*
 * Function needed to compute optimal parameters (required workspace is 3*n*n)
 */
static PetscInt ell(PetscBLASInt n,PetscScalar *A,PetscReal coeff,PetscInt m,PetscScalar *work,PetscRandom rand)
{
  PetscScalar    *Ascaled=work;
  PetscReal      nrm,alpha,beta,rwork[1];
  PetscInt       t;
  PetscBLASInt   i,j;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  beta = PetscPowReal(coeff,1.0/(2*m+1));
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      Ascaled[i+j*n] = beta*PetscAbsScalar(A[i+j*n]);
  nrm = LAPACKlange_("O",&n,&n,A,&n,rwork);
  ierr = PetscLogFlops(2.0*n*n);CHKERRQ(ierr);
  ierr = SlepcNormAm(n,Ascaled,2*m+1,work+n*n,rand,&alpha);CHKERRQ(ierr);
  alpha /= nrm;
  t = PetscMax((PetscInt)PetscCeilReal(PetscLogReal(2.0*alpha/PETSC_MACHINE_EPSILON)/PetscLogReal(2.0)/(2*m)),0);
  PetscFunctionReturn(t);
}

/*
 * Compute scaling parameter (s) and order of Pade approximant (m)  (required workspace is 4*n*n)
 */
static PetscErrorCode expm_params(PetscInt n,PetscScalar **Apowers,PetscInt *s,PetscInt *m,PetscScalar *work)
{
  PetscErrorCode  ierr;
  PetscScalar     sfactor,sone=1.0,szero=0.0,*A=Apowers[0],*Ascaled;
  PetscReal       d4,d6,d8,d10,eta1,eta3,eta4,eta5,rwork[1];
  PetscBLASInt    n_,n2,one=1;
  PetscRandom     rand;
  const PetscReal coeff[5] = { 9.92063492063492e-06, 9.94131285136576e-11,  /* backward error function */
                               2.22819456055356e-16, 1.69079293431187e-22, 8.82996160201868e-36 };
  const PetscReal theta[5] = { 1.495585217958292e-002,    /* m = 3  */
                               2.539398330063230e-001,    /* m = 5  */
                               9.504178996162932e-001,    /* m = 7  */
                               2.097847961257068e+000,    /* m = 9  */
                               5.371920351148152e+000 };  /* m = 13 */

  PetscFunctionBegin;
  *s = 0;
  *m = 13;
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_SELF,&rand);CHKERRQ(ierr);
  d4 = PetscPowReal(LAPACKlange_("O",&n_,&n_,Apowers[2],&n_,rwork),1.0/4.0);
  if (d4==0.0) { /* safeguard for the case A = 0 */
    *m = 3;
    goto done;
  }
  d6 = PetscPowReal(LAPACKlange_("O",&n_,&n_,Apowers[3],&n_,rwork),1.0/6.0);
  ierr = PetscLogFlops(2.0*n*n);CHKERRQ(ierr);
  eta1 = PetscMax(d4,d6);
  if (eta1<=theta[0] && !ell(n_,A,coeff[0],3,work,rand)) {
    *m = 3;
    goto done;
  }
  if (eta1<=theta[1] && !ell(n_,A,coeff[1],5,work,rand)) {
    *m = 5;
    goto done;
  }
  if (n<SMALLN) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[2],&n_,Apowers[2],&n_,&szero,work,&n_));
    d8 = PetscPowReal(LAPACKlange_("O",&n_,&n_,work,&n_,rwork),1.0/8.0);
    ierr = PetscLogFlops(2.0*n*n*n+1.0*n*n);CHKERRQ(ierr);
  } else {
    ierr = SlepcNormAm(n_,Apowers[2],2,work,rand,&d8);CHKERRQ(ierr);
    d8 = PetscPowReal(d8,1.0/8.0);
  }
  eta3 = PetscMax(d6,d8);
  if (eta3<=theta[2] && !ell(n_,A,coeff[2],7,work,rand)) {
    *m = 7;
    goto done;
  }
  if (eta3<=theta[3] && !ell(n_,A,coeff[3],9,work,rand)) {
    *m = 9;
    goto done;
  }
  if (n<SMALLN) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[2],&n_,Apowers[3],&n_,&szero,work,&n_));
    d10 = PetscPowReal(LAPACKlange_("O",&n_,&n_,work,&n_,rwork),1.0/10.0);
    ierr = PetscLogFlops(2.0*n*n*n+1.0*n*n);CHKERRQ(ierr);
  } else {
    ierr = SlepcNormAm(n_,Apowers[1],5,work,rand,&d10);CHKERRQ(ierr);
    d10 = PetscPowReal(d10,1.0/10.0);
  }
  eta4 = PetscMax(d8,d10);
  eta5 = PetscMin(eta3,eta4);
  *s = PetscMax((PetscInt)PetscCeilReal(PetscLogReal(eta5/theta[4])/PetscLogReal(2.0)),0);
  if (*s) {
    Ascaled = work+3*n*n;
    n2 = n_*n_;
    PetscStackCallBLAS("BLAScopy",BLAScopy_(&n2,A,&one,Ascaled,&one));
    sfactor = PetscPowRealInt(2.0,-(*s));
    PetscStackCallBLAS("BLASscal",BLASscal_(&n2,&sfactor,Ascaled,&one));
    ierr = PetscLogFlops(1.0*n*n);CHKERRQ(ierr);
  } else Ascaled = A;
  *s += ell(n_,Ascaled,coeff[4],13,work,rand);
done:
  ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
 * Matrix exponential implementation based on algorithm and matlab code by N. Higham and co-authors
 *
 *     N. J. Higham, "The scaling and squaring method for the matrix exponential
 *     revisited", SIAM J. Matrix Anal. Appl. 26(4):1179-1193, 2005.
 */
PetscErrorCode FNEvaluateFunctionMat_Exp_Higham(FN fn,Mat A,Mat B)
{
  PetscErrorCode    ierr;
  PetscBLASInt      n_,n2,*ipiv,info,one=1;
  PetscInt          n,m,j,s;
  PetscScalar       scale,smone=-1.0,sone=1.0,stwo=2.0,szero=0.0;
  PetscScalar       *Aa,*Ba,*Apowers[5],*Q,*P,*W,*work,*aux;
  const PetscScalar *c;
  const PetscScalar c3[4]   = { 120, 60, 12, 1 };
  const PetscScalar c5[6]   = { 30240, 15120, 3360, 420, 30, 1 };
  const PetscScalar c7[8]   = { 17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1 };
  const PetscScalar c9[10]  = { 17643225600, 8821612800, 2075673600, 302702400, 30270240,
                                2162160, 110880, 3960, 90, 1 };
  const PetscScalar c13[14] = { 64764752532480000, 32382376266240000, 7771770303897600,
                                1187353796428800,  129060195264000,   10559470521600,
                                670442572800,      33522128640,       1323241920,
                                40840800,          960960,            16380,  182,  1 };

  PetscFunctionBegin;
  ierr = MatDenseGetArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseGetArray(B,&Ba);CHKERRQ(ierr);
  ierr = MatGetSize(A,&n,NULL);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(n,&n_);CHKERRQ(ierr);
  n2 = n_*n_;
  ierr = PetscMalloc2(8*n*n,&work,n,&ipiv);CHKERRQ(ierr);

  /* Matrix powers */
  Apowers[0] = work;                  /* Apowers[0] = A   */
  Apowers[1] = Apowers[0] + n*n;      /* Apowers[1] = A^2 */
  Apowers[2] = Apowers[1] + n*n;      /* Apowers[2] = A^4 */
  Apowers[3] = Apowers[2] + n*n;      /* Apowers[3] = A^6 */
  Apowers[4] = Apowers[3] + n*n;      /* Apowers[4] = A^8 */

  ierr = PetscArraycpy(Apowers[0],Aa,n2);CHKERRQ(ierr);
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[0],&n_,Apowers[0],&n_,&szero,Apowers[1],&n_));
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[1],&n_,Apowers[1],&n_,&szero,Apowers[2],&n_));
  PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[1],&n_,Apowers[2],&n_,&szero,Apowers[3],&n_));
  ierr = PetscLogFlops(6.0*n*n*n);CHKERRQ(ierr);

  /* Compute scaling parameter and order of Pade approximant */
  ierr = expm_params(n,Apowers,&s,&m,Apowers[4]);CHKERRQ(ierr);

  if (s) { /* rescale */
    for (j=0;j<4;j++) {
      scale = PetscPowRealInt(2.0,-PetscMax(2*j,1)*s);
      PetscStackCallBLAS("BLASscal",BLASscal_(&n2,&scale,Apowers[j],&one));
    }
    ierr = PetscLogFlops(4.0*n*n);CHKERRQ(ierr);
  }

  /* Evaluate the Pade approximant */
  switch (m) {
    case 3:  c = c3;  break;
    case 5:  c = c5;  break;
    case 7:  c = c7;  break;
    case 9:  c = c9;  break;
    case 13: c = c13; break;
    default: SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong value of m %d",m);
  }
  P = Ba;
  Q = Apowers[4] + n*n;
  W = Q + n*n;
  switch (m) {
    case 3:
    case 5:
    case 7:
    case 9:
      if (m==9) PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[1],&n_,Apowers[3],&n_,&szero,Apowers[4],&n_));
      ierr = PetscArrayzero(P,n2);CHKERRQ(ierr);
      ierr = PetscArrayzero(Q,n2);CHKERRQ(ierr);
      for (j=0;j<n;j++) {
        P[j+j*n] = c[1];
        Q[j+j*n] = c[0];
      }
      for (j=m;j>=3;j-=2) {
        PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[j],Apowers[(j+1)/2-1],&one,P,&one));
        PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[j-1],Apowers[(j+1)/2-1],&one,Q,&one));
        ierr = PetscLogFlops(4.0*n*n);CHKERRQ(ierr);
      }
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[0],&n_,P,&n_,&szero,W,&n_));
      ierr = PetscLogFlops(2.0*n*n*n);CHKERRQ(ierr);
      SWAP(P,W,aux);
      break;
    case 13:
      /*  P = A*(Apowers[3]*(c[13]*Apowers[3] + c[11]*Apowers[2] + c[9]*Apowers[1])
              + c[7]*Apowers[3] + c[5]*Apowers[2] + c[3]*Apowers[1] + c[1]*I)       */
      PetscStackCallBLAS("BLAScopy",BLAScopy_(&n2,Apowers[3],&one,P,&one));
      PetscStackCallBLAS("BLASscal",BLASscal_(&n2,&c[13],P,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[11],Apowers[2],&one,P,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[9],Apowers[1],&one,P,&one));
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[3],&n_,P,&n_,&szero,W,&n_));
      ierr = PetscLogFlops(5.0*n*n+2.0*n*n*n);CHKERRQ(ierr);
      ierr = PetscArrayzero(P,n2);CHKERRQ(ierr);
      for (j=0;j<n;j++) P[j+j*n] = c[1];
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[7],Apowers[3],&one,P,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[5],Apowers[2],&one,P,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[3],Apowers[1],&one,P,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&sone,P,&one,W,&one));
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[0],&n_,W,&n_,&szero,P,&n_));
      ierr = PetscLogFlops(7.0*n*n+2.0*n*n*n);CHKERRQ(ierr);
      /*  Q = Apowers[3]*(c[12]*Apowers[3] + c[10]*Apowers[2] + c[8]*Apowers[1])
              + c[6]*Apowers[3] + c[4]*Apowers[2] + c[2]*Apowers[1] + c[0]*I        */
      PetscStackCallBLAS("BLAScopy",BLAScopy_(&n2,Apowers[3],&one,Q,&one));
      PetscStackCallBLAS("BLASscal",BLASscal_(&n2,&c[12],Q,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[10],Apowers[2],&one,Q,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[8],Apowers[1],&one,Q,&one));
      PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,Apowers[3],&n_,Q,&n_,&szero,W,&n_));
      ierr = PetscLogFlops(5.0*n*n+2.0*n*n*n);CHKERRQ(ierr);
      ierr = PetscArrayzero(Q,n2);CHKERRQ(ierr);
      for (j=0;j<n;j++) Q[j+j*n] = c[0];
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[6],Apowers[3],&one,Q,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[4],Apowers[2],&one,Q,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&c[2],Apowers[1],&one,Q,&one));
      PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&sone,W,&one,Q,&one));
      ierr = PetscLogFlops(7.0*n*n);CHKERRQ(ierr);
      break;
    default: SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong value of m %d",m);
  }
  PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&n2,&smone,P,&one,Q,&one));
  PetscStackCallBLAS("LAPACKgesv",LAPACKgesv_(&n_,&n_,Q,&n_,ipiv,P,&n_,&info));
  SlepcCheckLapackInfo("gesv",info);
  PetscStackCallBLAS("BLASscal",BLASscal_(&n2,&stwo,P,&one));
  for (j=0;j<n;j++) P[j+j*n] += 1.0;
  ierr = PetscLogFlops(2.0*n*n*n/3.0+4.0*n*n);CHKERRQ(ierr);

  /* Squaring */
  for (j=1;j<=s;j++) {
    PetscStackCallBLAS("BLASgemm",BLASgemm_("N","N",&n_,&n_,&n_,&sone,P,&n_,P,&n_,&szero,W,&n_));
    SWAP(P,W,aux);
  }
  if (P!=Ba) { ierr = PetscArraycpy(Ba,P,n2);CHKERRQ(ierr); }
  ierr = PetscLogFlops(2.0*n*n*n*s);CHKERRQ(ierr);

  ierr = PetscFree2(work,ipiv);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(A,&Aa);CHKERRQ(ierr);
  ierr = MatDenseRestoreArray(B,&Ba);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FNView_Exp(FN fn,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscBool      isascii;
  char           str[50];
  const char     *methodname[] = {
                  "scaling & squaring, [m/m] Pade approximant (Higham)",
                  "scaling & squaring, [6/6] Pade approximant",
                  "scaling & squaring, subdiagonal Pade approximant (product form)",
                  "scaling & squaring, subdiagonal Pade approximant (partial fraction)"
  };
  const int      nmeth=sizeof(methodname)/sizeof(methodname[0]);

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    if (fn->beta==(PetscScalar)1.0) {
      if (fn->alpha==(PetscScalar)1.0) {
        ierr = PetscViewerASCIIPrintf(viewer,"  Exponential: exp(x)\n");CHKERRQ(ierr);
      } else {
        ierr = SlepcSNPrintfScalar(str,50,fn->alpha,PETSC_TRUE);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"  Exponential: exp(%s*x)\n",str);CHKERRQ(ierr);
      }
    } else {
      ierr = SlepcSNPrintfScalar(str,50,fn->beta,PETSC_TRUE);CHKERRQ(ierr);
      if (fn->alpha==(PetscScalar)1.0) {
        ierr = PetscViewerASCIIPrintf(viewer,"  Exponential: %s*exp(x)\n",str);CHKERRQ(ierr);
      } else {
        ierr = PetscViewerASCIIPrintf(viewer,"  Exponential: %s",str);CHKERRQ(ierr);
        ierr = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);CHKERRQ(ierr);
        ierr = SlepcSNPrintfScalar(str,50,fn->alpha,PETSC_TRUE);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer,"*exp(%s*x)\n",str);CHKERRQ(ierr);
        ierr = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);CHKERRQ(ierr);
      }
    }
    if (fn->method<nmeth) {
      ierr = PetscViewerASCIIPrintf(viewer,"  computing matrix functions with: %s\n",methodname[fn->method]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode FNCreate_Exp(FN fn)
{
  PetscFunctionBegin;
  fn->ops->evaluatefunction       = FNEvaluateFunction_Exp;
  fn->ops->evaluatederivative     = FNEvaluateDerivative_Exp;
  fn->ops->evaluatefunctionmat[0] = FNEvaluateFunctionMat_Exp_Higham;
  fn->ops->evaluatefunctionmat[1] = FNEvaluateFunctionMat_Exp_Pade;
  fn->ops->evaluatefunctionmat[2] = FNEvaluateFunctionMat_Exp_GuettelNakatsukasa; /* product form */
  fn->ops->evaluatefunctionmat[3] = FNEvaluateFunctionMat_Exp_GuettelNakatsukasa; /* partial fraction */
  fn->ops->view                   = FNView_Exp;
  PetscFunctionReturn(0);
}

