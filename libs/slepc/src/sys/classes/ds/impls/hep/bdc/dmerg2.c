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

PetscErrorCode BDC_dmerg2_(const char *jobz,PetscBLASInt rkct,PetscBLASInt n,
        PetscReal *ev,PetscReal *q,PetscBLASInt ldq,PetscBLASInt *indxq,
        PetscReal *rho,PetscReal *u,PetscBLASInt sbrkp1,PetscReal *v,
        PetscBLASInt sbrk,PetscBLASInt cutpnt,PetscReal *work,PetscBLASInt lwork,
        PetscBLASInt *iwork,PetscReal tol,PetscBLASInt *info,PetscBLASInt jobz_len)
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

/*  DMERG2 computes the updated eigensystem of a diagonal matrix after */
/*  modification by a rank-one symmetric matrix.  The diagonal matrix */
/*  consists of two diagonal submatrices, and the vectors defining the */
/*  rank-1 matrix similarly have two underlying subvectors each. */
/*  The dimension of the first subproblem is CUTPNT, the dimension of */
/*  the second subproblem is N-CUTPNT. */

/*  T = Q(in) (EV(in) + RHO * Z*Z') Q'(in) = Q(out) * EV(out) * Q'(out) */

/*     where Z = Q'[V U']', where V is a row vector and U is a column */
/*     vector with dimensions corresponding to the two underlying */
/*     subproblems. */

/*     The eigenvectors of the original matrix are stored in Q, and the */
/*     eigenvalues in EV.  The algorithm consists of three stages: */

/*        The first stage consists of deflating the size of the problem */
/*        when there are multiple eigenvalues or if there is a zero in */
/*        the Z vector.  For each such occurrence the dimension of the */
/*        secular equation problem is reduced by one.  This stage is */
/*        performed by the routine DSRTDF. */

/*        The second stage consists of calculating the updated */
/*        eigenvalues. This is done by finding the roots of the secular */
/*        equation via the routine DLAED4 (as called by DLAED3M). */
/*        This routine also calculates the eigenvectors of the current */
/*        problem. */

/*        If(JOBZ.EQ.'D') then the final stage consists */
/*        of computing the updated eigenvectors directly using the updated */
/*        eigenvalues. The eigenvectors for the current problem are multiplied */
/*        with the eigenvectors from the overall problem. */

/*  Arguments */
/*  ========= */

/*  JOBZ   (input) CHARACTER*1 */
/*          = 'N': Compute eigenvalues only (not implemented); */
/*          = 'D': Compute eigenvalues and eigenvectors. */
/*                 Eigenvectors are accumulated in the divide-and-conquer */
/*                 process. */

/*  RKCT   (input) INTEGER */
/*         The number of the rank modification which is accounted for */
/*         (RKCT >= 1). Required parameter, because the update operation of the */
/*         modification vector can be performed much more efficiently */
/*         if RKCT.EQ.1. In that case, the eigenvector matrix is still */
/*         block-diagonal. For RKCT.GE.2 the eigenvector matrix for the update */
/*         operation has filled up and is a full matrix. */

/*  N      (input) INTEGER */
/*         The dimension of the symmetric block tridiagonal matrix. */
/*         N >= 0. */

/*  EV     (input/output) DOUBLE PRECISION array, dimension (N) */
/*         On entry, the eigenvalues (=diagonal values) of the */
/*         rank-1-perturbed matrix. */
/*         On exit, the eigenvalues of the repaired matrix. */

/*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*         On entry, the eigenvectors of the rank-1-perturbed matrix. */
/*         On exit, the eigenvectors of the repaired tridiagonal matrix. */

/*  LDQ    (input) INTEGER */
/*         The leading dimension of the array Q.  LDQ >= max(1,N). */

/*  INDXQ  (input/output) INTEGER array, dimension (N) */
/*         On entry, the permutation which separately sorts the two */
/*         subproblems in EV into ascending order. */
/*         On exit, the permutation which will reintegrate the */
/*         subproblems back into sorted order, */
/*         i.e. EV(INDXQ(I = 1, N)) will be in ascending order. */

/*  RHO    (input/output) DOUBLE PRECISION */
/*         The scalar in the rank-1 perturbation. Modified (multiplied */
/*         by 2) in DSRTDF. */

/*  U      (input) DOUBLE PRECISION array; dimension (SBRKP1), where SBRKP1 */
/*         is the size of the first (original) block after CUTPNT. */
/*         The column vector of the rank-1 subdiagonal connecting the */
/*         two diagonal subproblems. */
/*         Theoretically, zero entries might have to be appended after U */
/*         in order to make it have dimension (N-CUTPNT). However, this */
/*         is not required because it can be accounted for in the */
/*         matrix-vector product using the argument SBRKP1. */

/*  SBRKP1 (input) INTEGER */
/*         Dimension of the relevant (non-zero) part of the vector U. */
/*         Equal to the size of the first original block after the */
/*         breakpoint. */

/*  V      (input) DOUBLE PRECISION array; dimension (SBRK), where SBRK */
/*         is the size of the last (original) block before CUTPNT. */
/*         The row vector of the rank-1 subdiagonal connecting the two */
/*         diagonal subproblems. */
/*         Theoretically, zero entries might have to be inserted in front */
/*         of V in order to make it have dimension (CUTPNT). However, this */
/*         is not required because it can be accounted for in the */
/*         matrix-vector product using the argument SBRK. */

/*  SBRK   (input) INTEGER */
/*         Dimension of the relevant (non-zero) part of the vector V. */
/*         Equal to the size of the last original block before the */
/*         breakpoint. */

/*  CUTPNT (input) INTEGER */
/*         The location of the last eigenvalue of the leading diagonal */
/*         block.  min(1,N) <= CUTPNT <= max(1,N). */

/*  WORK   (workspace) DOUBLE PRECISION array, dimension (LWORK) */

/*  LWORK  (input) INTEGER */
/*         The dimension of the array WORK. */
/*         In order to guarantee correct results in all cases, */
/*         LWORK must be at least (2*N**2 + 3*N). In many cases, */
/*         less workspace is required. The absolute minimum required is */
/*         (N**2 + 3*N). */
/*         If the workspace provided is not sufficient, the routine will */
/*         return a corresponding error code and report how much workspace */
/*         was missing (see INFO). */
/*         NOTE: This parameter is needed for determining whether enough */
/*               workspace is provided, and, if not, for computing how */
/*               much workspace is needed. */

/*  IWORK  (workspace) INTEGER array, dimension (4*N) */

/*  TOL    (input) DOUBLE PRECISION */
/*         User specified deflation tolerance for the routine DSRTDF. */

/*  INFO   (output) INTEGER */
/*          = 0:  successful exit. */
/*          < -200: not enough workspace */
/*                ABS(INFO + 200) numbers have to be stored in addition */
/*                to the workspace provided, otherwise some eigenvectors */
/*                will be incorrect. */
/*          < 0, > -99:  if INFO.EQ.-i, the i-th argument had an */
/*                       illegal value. */
/*          > 0:  if INFO.EQ.1, an eigenvalue did not converge */
/*                if INFO.EQ.2, the deflation counters DZ and DE do not sum */
/*                              up to the total number N-K of components */
/*                              deflated */

/*  Further Details */
/*  =============== */

/*  Based on code written by */
/*     Wilfried Gansterer and Bob Ward, */
/*     Department of Computer Science, University of Tennessee */

/*  Based on the design of the LAPACK code Dlaed1.f written by Jeff */
/*  Rutter, Computer Science Division, University of California at */
/*  Berkeley, and modified by Francoise Tisseur, University of Tennessee. */

/*  ===================================================================== */

  PetscBLASInt   i, k, n1, n2, de, is, dz, iw, iz, iq2, nmc, cpp1;
  PetscBLASInt   indx, indxc, indxp, lwmin, idlmda;
  PetscBLASInt   spneed, coltyp, tmpcut, i__1, i__2, one=1, mone=-1;
  char           defl[1];
  PetscReal      done = 1.0, dzero = 0.0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *info = 0;
  lwmin = n*n + n * 3;

  if (n < 0) *info = -3;
  else if (ldq < PetscMax(1,n)) *info = -6;
  else if (cutpnt < PetscMin(1,n) || cutpnt > PetscMax(1,n)) *info = -13;
  else if (lwork < lwmin) *info = -15;
  if (*info) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong argument %d in DMERG2",-(*info));

/* **************************************************************************** */

  /* Quick return if possible */

  if (n == 0) PetscFunctionReturn(0);

/* **************************************************************************** */

  /* The following values are integer pointers which indicate */
  /* the portion of the workspace used by a particular array in DSRTDF */
  /* and DLAED3M. */

  iz = 0;
  idlmda = iz + n;
  iw = idlmda + n;
  iq2 = iw + n;
  is = iq2 + n * n;

  /* After the pointer IS the matrix S is stored and read in WORK */
  /* in the routine DLAED3M. */

  indx = 0;
  indxc = indx + n;
  coltyp = indxc + n;
  indxp = coltyp + n;

  /* If eigenvectors are to be accumulated in the divide-and-conquer */
  /* process (JOBZ.EQ.'D') form the z-vector which consists of */
  /* Q_1^T * V and Q_2^T * U. */

  cpp1 = cutpnt + 1;
  if (rkct == 1) {

    /* for the first rank modification the eigenvector matrix has */
    /* special block-diagonal structure and therefore Q_1^T * V and */
    /* Q_2^T * U can be formed separately. */

    PetscStackCallBLAS("BLASgemv",BLASgemv_("T", &sbrk, &cutpnt, &done,
              &q[cutpnt - sbrk], &ldq, v, &one, &dzero, &work[iz], &one));
    nmc = n - cutpnt;
    PetscStackCallBLAS("BLASgemv",BLASgemv_("T", &sbrkp1, &nmc, &done,
              &q[cpp1-1 + (cpp1-1)*ldq], &ldq, u,
              &one, &dzero, &work[iz + cutpnt], &one));

  } else {

    /* for the higher rank modifications, the vectors V and U */
    /* have to be multiplied with the full eigenvector matrix */

    PetscStackCallBLAS("BLASgemv",BLASgemv_("T", &sbrk, &n, &done,
              &q[cutpnt - sbrk], &ldq, v, &one, &dzero, &work[iz], &one));
    PetscStackCallBLAS("BLASgemv",BLASgemv_("T", &sbrkp1, &n, &done, &q[cpp1-1],
              &ldq, u, &one, &done, &work[iz], &one));

  }

/* **************************************************************************** */

  /* Deflate eigenvalues. */

  if (rkct == 1) {

    /* for the first rank modification we need the actual cutpoint */

    n1 = cutpnt;
    tmpcut = cutpnt;
  } else {

    /* for the later rank modifications there is no actual cutpoint any more */

    n1 = n;

    /* The original value of CUTPNT has to be preserved for the next time */
    /* this subroutine is called (therefore, CUTPNT is an INPUT parameter */
    /* and not to be changed). Thus, assign N to TMPCUT and use the local */
    /* variable TMPCUT from now on for the cut point. */

    tmpcut = n;
  }

  /* initialize the flag DEFL (indicates whether deflation occurred - */
  /* this information is needed later in DLAED3M) */

  *(unsigned char *)defl = '0';

  /* call DSRTDF for deflation */

  ierr = BDC_dsrtdf_(&k, n, n1, ev, q, ldq, indxq, rho, &work[iz],
          &work[idlmda], &work[iw], &work[iq2], &iwork[indx],
          &iwork[indxc], &iwork[indxp], &iwork[coltyp], tol, &dz, &de, info);
          CHKERRQ(ierr);
  if (*info) SETERRQ1(PETSC_COMM_SELF,1,"dmerg2: error in dsrtdf, info = %d",*info);

  if (k < n) {

   /* ....some deflation occurred in dsrtdf, set the flag DEFL */
   /*     (needed in DLAED3M.f, since Givens rotations need to be */
   /*     applied to the eigenvector matrix only if some deflation */
   /*     happened) */

    *(unsigned char *)defl = '1';
  }

/* **************************************************************************** */

  /* Solve the Secular Equation. */

  if (k != 0 || k == 0) {

    /* ....not everything was deflated. */

    /* ....check whether enough workspace is available: */

    /* Note that the following (upper) bound SPNEED for the workspace */
    /* requirements should also hold in the extreme case TMPCUT=N, */
    /* which happens for every rank modification after the first one. */

    i__1 = (iwork[coltyp] + iwork[coltyp+1]) * k;
    i__2 = (iwork[coltyp+1] + iwork[coltyp + 2]) * k;
    spneed = is + PetscMax(i__1,i__2) - 1;

    if (spneed > lwork) SETERRQ1(PETSC_COMM_SELF,1,"dmerg2: Workspace needed exceeds the workspace provided by %d numbers",spneed-lwork);

    /* calling DLAED3M for solving the secular equation. */

    ierr = BDC_dlaed3m_(jobz, defl, k, n, tmpcut, ev, q, ldq,
                *rho, &work[idlmda], &work[iq2], &iwork[indxc], &iwork[coltyp],
                &work[iw], &work[is], info, 1, 1);CHKERRQ(ierr);
    if (*info) SETERRQ1(PETSC_COMM_SELF,1,"dmerg2: error in dlaed3m, info = %d",*info);

    /* Prepare the INDXQ sorting permutation. */

    n1 = k;
    n2 = n - k;
    PetscStackCallBLAS("LAPACKlamrg",LAPACKlamrg_(&n1, &n2, ev, &one, &mone, indxq));
    if (k == 0) for (i = 0; i < n; ++i) indxq[i] = i+1;

  } else {

    /* ....everything was deflated (K.EQ.0) */

    for (i = 0; i < n; ++i) indxq[i] = i+1;
  }
  PetscFunctionReturn(0);
}

