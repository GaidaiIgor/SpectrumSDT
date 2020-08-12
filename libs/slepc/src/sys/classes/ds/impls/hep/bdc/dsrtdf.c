/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   BDC - Block-divide and conquer (see description in README file)
*/

#include <slepc/private/dsimpl.h>
#include <slepcblaslapack.h>

PetscErrorCode BDC_dsrtdf_(PetscBLASInt *k,PetscBLASInt n,PetscBLASInt n1,
        PetscReal *d,PetscReal *q,PetscBLASInt ldq,PetscBLASInt *indxq,
        PetscReal *rho,PetscReal *z,PetscReal *dlamda,PetscReal *w,
        PetscReal *q2,PetscBLASInt *indx,PetscBLASInt *indxc,PetscBLASInt *indxp,
        PetscBLASInt *coltyp,PetscReal reltol,PetscBLASInt *dz,PetscBLASInt *de,
        PetscBLASInt *info)
{
/*  -- Routine written in LAPACK Version 3.0 style -- */
/* *************************************************** */
/*     Written by */
/*     Michael Moldaschl and Wilfried Gansterer */
/*     University of Vienna */
/*     last modification: March 16, 2014 */

/*     Small adaptations of original code written by */
/*     Wilfried Gansterer and Bob Ward, */
/*     Department of Computer Science, University of Tennessee */
/*     see https://doi.org/10.1137/S1064827501399432 */
/* *************************************************** */

/*  Purpose */
/*  ======= */

/*  DSRTDF merges the two sets of eigenvalues of a rank-one modified */
/*  diagonal matrix D+ZZ^T together into a single sorted set. Then it tries */
/*  to deflate the size of the problem. */
/*  There are two ways in which deflation can occur:  when two or more */
/*  eigenvalues of D are close together or if there is a tiny entry in the */
/*  Z vector.  For each such occurrence the order of the related secular */
/*  equation problem is reduced by one. */

/*  Arguments */
/*  ========= */

/*  K      (output) INTEGER */
/*         The number of non-deflated eigenvalues, and the order of the */
/*         related secular equation. 0 <= K <=N. */

/*  N      (input) INTEGER */
/*         The dimension of the diagonal matrix.  N >= 0. */

/*  N1     (input) INTEGER */
/*         The location of the last eigenvalue in the leading sub-matrix. */
/*         min(1,N) <= N1 <= max(1,N). */

/*  D      (input/output) DOUBLE PRECISION array, dimension (N) */
/*         On entry, D contains the eigenvalues of the two submatrices to */
/*         be combined. */
/*         On exit, D contains the trailing (N-K) updated eigenvalues */
/*         (those which were deflated) sorted into increasing order. */

/*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N) */
/*         On entry, Q contains the eigenvectors of two submatrices in */
/*         the two square blocks with corners at (1,1), (N1,N1) */
/*         and (N1+1, N1+1), (N,N). */
/*         On exit, Q contains the trailing (N-K) updated eigenvectors */
/*         (those which were deflated) in its last N-K columns. */

/*  LDQ    (input) INTEGER */
/*         The leading dimension of the array Q.  LDQ >= max(1,N). */

/*  INDXQ  (input/output) INTEGER array, dimension (N) */
/*         The permutation which separately sorts the two sub-problems */
/*         in D into ascending order.  Note that elements in the second */
/*         half of this permutation must first have N1 added to their */
/*         values. Destroyed on exit. */

/*  RHO    (input/output) DOUBLE PRECISION */
/*         On entry, the off-diagonal element associated with the rank-1 */
/*         cut which originally split the two submatrices which are now */
/*         being recombined. */
/*         On exit, RHO has been modified to the value required by */
/*         DLAED3M (made positive and multiplied by two, such that */
/*         the modification vector has norm one). */

/*  Z      (input/output) DOUBLE PRECISION array, dimension (N) */
/*         On entry, Z contains the updating vector. */
/*         On exit, the contents of Z has been destroyed by the updating */
/*         process. */

/*  DLAMDA (output) DOUBLE PRECISION array, dimension (N) */
/*         A copy of the first K non-deflated eigenvalues which */
/*         will be used by DLAED3M to form the secular equation. */

/*  W      (output) DOUBLE PRECISION array, dimension (N) */
/*         The first K values of the final deflation-altered z-vector */
/*         which will be passed to DLAED3M. */

/*  Q2     (output) DOUBLE PRECISION array, dimension */
/*         (N*N) (if everything is deflated) or */
/*         (N1*(COLTYP(1)+COLTYP(2)) + (N-N1)*(COLTYP(2)+COLTYP(3))) */
/*         (if not everything is deflated) */
/*         If everything is deflated, then N*N intermediate workspace */
/*         is needed in Q2. */
/*         If not everything is deflated, Q2 contains on exit a copy of the */
/*         first K non-deflated eigenvectors which will be used by */
/*         DLAED3M in a matrix multiply (DGEMM) to accumulate the new */
/*         eigenvectors, followed by the N-K deflated eigenvectors. */

/*  INDX   (workspace) INTEGER array, dimension (N) */
/*         The permutation used to sort the contents of DLAMDA into */
/*         ascending order. */

/*  INDXC  (output) INTEGER array, dimension (N) */
/*         The permutation used to arrange the columns of the deflated */
/*         Q matrix into three groups:  columns in the first group contain */
/*         non-zero elements only at and above N1 (type 1), columns in the */
/*         second group are dense (type 2), and columns in the third group */
/*         contain non-zero elements only below N1 (type 3). */

/*  INDXP  (workspace) INTEGER array, dimension (N) */
/*         The permutation used to place deflated values of D at the end */
/*         of the array.  INDXP(1:K) points to the nondeflated D-values */
/*         and INDXP(K+1:N) points to the deflated eigenvalues. */

/*  COLTYP (workspace/output) INTEGER array, dimension (N) */
/*         During execution, a label which will indicate which of the */
/*         following types a column in the Q2 matrix is: */
/*         1 : non-zero in the upper half only; */
/*         2 : dense; */
/*         3 : non-zero in the lower half only; */
/*         4 : deflated. */
/*         On exit, COLTYP(i) is the number of columns of type i, */
/*         for i=1 to 4 only. */

/*  RELTOL (input) DOUBLE PRECISION */
/*         User specified deflation tolerance. If RELTOL.LT.20*EPS, then */
/*         the standard value used in the original LAPACK routines is used. */

/*  DZ     (output) INTEGER, DZ.GE.0 */
/*         counts the deflation due to a small component in the modification */
/*         vector. */

/*  DE     (output) INTEGER, DE.GE.0 */
/*         counts the deflation due to close eigenvalues. */

/*  INFO   (output) INTEGER */
/*          = 0:  successful exit. */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

/*  Further Details */
/*  =============== */

/*  Based on code written by */
/*  Wilfried Gansterer and Bob Ward, */
/*  Department of Computer Science, University of Tennessee */

/*  Based on the design of the LAPACK code DLAED2 with modifications */
/*  to allow a block divide and conquer scheme */

/*  DLAED2 was written by Jeff Rutter, Computer Science Division, University */
/*  of California at Berkeley, USA, and modified by Francoise Tisseur, */
/*  University of Tennessee. */

/*  ===================================================================== */

  PetscReal    c, s, t, eps, tau, tol, dmax, dmone = -1.;
  PetscBLASInt i, j, i1, k2, n2, ct, nj, pj=0, js, iq1, iq2;
  PetscBLASInt psm[4], imax, jmax, ctot[4], factmp=1, one=1;

  PetscFunctionBegin;
  *info = 0;

  if (n < 0) *info = -2;
  else if (n1 < PetscMin(1,n) || n1 > PetscMax(1,n)) *info = -3;
  else if (ldq < PetscMax(1,n)) *info = -6;
  if (*info) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong argument %d in DSRTDF",-(*info));

  /* Initialize deflation counters */

  *dz = 0;
  *de = 0;

/* **************************************************************************** */

  /* Quick return if possible */

  if (n == 0) PetscFunctionReturn(0);

/* **************************************************************************** */

  n2 = n - n1;

  /* Modify Z so that RHO is positive. */

  if (*rho < 0.) PetscStackCallBLAS("BLASscal",BLASscal_(&n2, &dmone, &z[n1], &one));

  /* Normalize z so that norm2(z) = 1.  Since z is the concatenation of */
  /* two normalized vectors, norm2(z) = sqrt(2). (NOTE that this holds also */
  /* here in the approximate block-tridiagonal D&C because the two vectors are */
  /* singular vectors, but it would not necessarily hold in a D&C for a */
  /* general banded matrix !) */

  t = 1. / PetscSqrtReal(2.);
  PetscStackCallBLAS("BLASscal",BLASscal_(&n, &t, z, &one));

  /* NOTE: at this point the value of RHO is modified in order to */
  /*       normalize Z:    RHO = ABS( norm2(z)**2 * RHO */
  /*       in our case:    norm2(z) = sqrt(2), and therefore: */

  *rho = PetscAbs(*rho * 2.);

  /* Sort the eigenvalues into increasing order */

  for (i = n1; i < n; ++i) indxq[i] += n1;

  /* re-integrate the deflated parts from the last pass */

  for (i = 0; i < n; ++i) dlamda[i] = d[indxq[i]-1];
  PetscStackCallBLAS("LAPACKlamrg",LAPACKlamrg_(&n1, &n2, dlamda, &one, &one, indxc));
  for (i = 0; i < n; ++i) indx[i] = indxq[indxc[i]-1];
  for (i = 0; i < n; ++i) indxq[i] = 0;

  /* Calculate the allowable deflation tolerance */

  /* imax = BLASamax_(&n, z, &one); */
  imax = 1;
  dmax = PetscAbsReal(z[0]);
  for (i=1;i<n;i++) {
    if (PetscAbsReal(z[i])>dmax) {
      imax = i+1;
      dmax = PetscAbsReal(z[i]);
    }
  }
  /* jmax = BLASamax_(&n, d, &one); */
  jmax = 1;
  dmax = PetscAbsReal(d[0]);
  for (i=1;i<n;i++) {
    if (PetscAbsReal(d[i])>dmax) {
      jmax = i+1;
      dmax = PetscAbsReal(d[i]);
    }
  }

  eps = LAPACKlamch_("Epsilon");
  if (reltol < eps * 20) {
    /* use the standard deflation tolerance from the original LAPACK code */
    tol = eps * 8. * PetscMax(PetscAbs(d[jmax-1]),PetscAbs(z[imax-1]));
  } else {
    /* otherwise set TOL to the input parameter RELTOL ! */
    tol = reltol * PetscMax(PetscAbs(d[jmax-1]),PetscAbs(z[imax-1]));
  }

  /* If the rank-1 modifier is small enough, no more needs to be done */
  /* except to reorganize Q so that its columns correspond with the */
  /* elements in D. */

  if (*rho * PetscAbs(z[imax-1]) <= tol) {

    /* complete deflation because of small rank-one modifier */
    /* NOTE: in this case we need N*N space in the array Q2 ! */

    *dz = n; *k = 0;
    iq2 = 1;
    for (j = 0; j < n; ++j) {
      i = indx[j]; indxq[j] = i;
      PetscStackCallBLAS("BLAScopy",BLAScopy_(&n, &q[(i-1)*ldq], &one, &q2[iq2-1], &one));
      dlamda[j] = d[i-1];
      iq2 += n;
    }
    for (j=0;j<n;j++) for (i=0;i<n;i++) q[i+j*ldq] = q2[i+j*n];
    PetscStackCallBLAS("BLAScopy",BLAScopy_(&n, dlamda, &one, d, &one));
    PetscFunctionReturn(0);
  }

  /* If there are multiple eigenvalues then the problem deflates.  Here */
  /* the number of equal eigenvalues is found.  As each equal */
  /* eigenvalue is found, an elementary reflector is computed to rotate */
  /* the corresponding eigensubspace so that the corresponding */
  /* components of Z are zero in this new basis. */

  /* initialize the column types */

  /* first N1 columns are initially (before deflation) of type 1 */
  for (i = 0; i < n1; ++i) coltyp[i] = 1;
  /* columns N1+1 to N are initially (before deflation) of type 3 */
  for (i = n1; i < n; ++i) coltyp[i] = 3;

  *k = 0;
  k2 = n + 1;
  for (j = 0; j < n; ++j) {
    nj = indx[j]-1;
    if (*rho * PetscAbs(z[nj]) <= tol) {

      /* Deflate due to small z component. */
      ++(*dz);
      --k2;

      /* deflated eigenpair, therefore column type 4 */
      coltyp[nj] = 4;
      indxp[k2-1] = nj+1;
      indxq[k2-1] = nj+1;
      if (j+1 == n) goto L100;
    } else {

      /* not deflated */
      pj = nj;
      goto L80;
    }
  }
  factmp = 1;
L80:
  ++j;
  nj = indx[j]-1;
  if (j+1 > n) goto L100;
  if (*rho * PetscAbs(z[nj]) <= tol) {

    /* Deflate due to small z component. */
    ++(*dz);
    --k2;
    coltyp[nj] = 4;
    indxp[k2-1] = nj+1;
    indxq[k2-1] = nj+1;
  } else {

    /* Check if eigenvalues are close enough to allow deflation. */
    s = z[pj];
    c = z[nj];

    /* Find sqrt(a**2+b**2) without overflow or */
    /* destructive underflow. */

    tau = SlepcAbs(c, s);
    t = d[nj] - d[pj];
    c /= tau;
    s = -s / tau;
    if (PetscAbs(t * c * s) <= tol) {

      /* Deflate due to close eigenvalues. */
      ++(*de);
      z[nj] = tau;
      z[pj] = 0.;
      if (coltyp[nj] != coltyp[pj]) coltyp[nj] = 2;

        /* if deflation happens across the two eigenvector blocks */
        /* (eigenvalues corresponding to different blocks), then */
        /* on column in the eigenvector matrix fills up completely */
        /* (changes from type 1 or 3 to type 2) */

        /* eigenpair PJ is deflated, therefore the column type changes */
        /* to 4 */

        coltyp[pj] = 4;
        PetscStackCallBLAS("BLASrot",BLASrot_(&n, &q[pj*ldq], &one, &q[nj*ldq], &one, &c, &s));
        t = d[pj] * (c * c) + d[nj] * (s * s);
        d[nj] = d[pj] * (s * s) + d[nj] * (c * c);
        d[pj] = t;
        --k2;
        i = 1;
L90:
        if (k2 + i <= n) {
          if (d[pj] < d[indxp[k2 + i-1]-1]) {
            indxp[k2 + i - 2] = indxp[k2 + i - 1];
            indxq[k2 + i - 2] = indxq[k2 + i - 1];
            indxp[k2 + i - 1] = pj + 1;
            indxq[k2 + i - 2] = pj + 1;
            ++i;
            goto L90;
          } else {
            indxp[k2 + i - 2] = pj + 1;
            indxq[k2 + i - 2] = pj + 1;
          }
        } else {
          indxp[k2 + i - 2] = pj + 1;
          indxq[k2 + i - 2] = pj + 1;
        }
        pj = nj;
        factmp = -1;
      } else {

      /* do not deflate */
      ++(*k);
      dlamda[*k-1] = d[pj];
      w[*k-1] = z[pj];
      indxp[*k-1] = pj + 1;
      indxq[*k-1] = pj + 1;
      factmp = 1;
      pj = nj;
    }
  }
  goto L80;
L100:

  /* Record the last eigenvalue. */
  ++(*k);
  dlamda[*k-1] = d[pj];
  w[*k-1] = z[pj];
  indxp[*k-1] = pj+1;
  indxq[*k-1] = (pj+1) * factmp;

  /* Count up the total number of the various types of columns, then */
  /* form a permutation which positions the four column types into */
  /* four uniform groups (although one or more of these groups may be */
  /* empty). */

  for (j = 0; j < 4; ++j) ctot[j] = 0;
  for (j = 0; j < n; ++j) {
    ct = coltyp[j];
    ++ctot[ct-1];
  }

  /* PSM(*) = Position in SubMatrix (of types 1 through 4) */
  psm[0] = 1;
  psm[1] = ctot[0] + 1;
  psm[2] = psm[1] + ctot[1];
  psm[3] = psm[2] + ctot[2];
  *k = n - ctot[3];

  /* Fill out the INDXC array so that the permutation which it induces */
  /* will place all type-1 columns first, all type-2 columns next, */
  /* then all type-3's, and finally all type-4's. */

  for (j = 0; j < n; ++j) {
    js = indxp[j];
    ct = coltyp[js-1];
    indx[psm[ct - 1]-1] = js;
    indxc[psm[ct - 1]-1] = j+1;
    ++psm[ct - 1];
  }

  /* Sort the eigenvalues and corresponding eigenvectors into DLAMDA */
  /* and Q2 respectively.  The eigenvalues/vectors which were not */
  /* deflated go into the first K slots of DLAMDA and Q2 respectively, */
  /* while those which were deflated go into the last N - K slots. */

  i = 0;
  iq1 = 0;
  iq2 = (ctot[0] + ctot[1]) * n1;
  for (j = 0; j < ctot[0]; ++j) {
    js = indx[i];
    PetscStackCallBLAS("BLAScopy",BLAScopy_(&n1, &q[(js-1)*ldq], &one, &q2[iq1], &one));
    z[i] = d[js-1];
    ++i;
    iq1 += n1;
  }

  for (j = 0; j < ctot[1]; ++j) {
    js = indx[i];
    PetscStackCallBLAS("BLAScopy",BLAScopy_(&n1, &q[(js-1)*ldq], &one, &q2[iq1], &one));
    PetscStackCallBLAS("BLAScopy",BLAScopy_(&n2, &q[n1+(js-1)*ldq], &one, &q2[iq2], &one));
    z[i] = d[js-1];
    ++i;
    iq1 += n1;
    iq2 += n2;
  }

  for (j = 0; j < ctot[2]; ++j) {
    js = indx[i];
    PetscStackCallBLAS("BLAScopy",BLAScopy_(&n2, &q[n1+(js-1)*ldq], &one, &q2[iq2], &one));
    z[i] = d[js-1];
    ++i;
    iq2 += n2;
  }

  iq1 = iq2;
  for (j = 0; j < ctot[3]; ++j) {
    js = indx[i];
    PetscStackCallBLAS("BLAScopy",BLAScopy_(&n, &q[(js-1)*ldq], &one, &q2[iq2], &one));
    iq2 += n;
    z[i] = d[js-1];
    ++i;
  }

  /* The deflated eigenvalues and their corresponding vectors go back */
  /* into the last N - K slots of D and Q respectively. */

  for (j=0;j<ctot[3];j++) for (i=0;i<n;i++) q[i+(j+*k)*ldq] = q2[iq1+i+j*n];
  i1 = n - *k;
  PetscStackCallBLAS("BLAScopy",BLAScopy_(&i1, &z[*k], &one, &d[*k], &one));

  /* Copy CTOT into COLTYP for referencing in DLAED3M. */

  for (j = 0; j < 4; ++j) coltyp[j] = ctot[j];
  PetscFunctionReturn(0);
}

