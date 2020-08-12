#include <petsc/private/petscimpl.h>
#include <petsc/private/dtimpl.h> /*I "petscdt.h" I*/

/*MC
   PetscDTAltV - An interface for common operations on k-forms, also known as alternating algebraic forms or alternating k-linear maps.
   The name of the interface comes from the notation "Alt V" for the algebra of all k-forms acting vectors in the space V, also known as the exterior algebra of V*.

   A recommended reference for this material is Section 2 "Exterior algebra and exterior calculus" in "Finite element
   exterior calculus, homological techniques, and applications", by Arnold, Falk, & Winther (2006, doi:10.1017/S0962492906210018).

   A k-form w (k is called the "form degree" of w) is an alternating k-linear map acting on tuples (v_1, ..., v_k) of
   vectors from a vector space V and producing a real number:
   - alternating: swapping any two vectors in a tuple reverses the sign of the result, e.g. w(v_1, v_2, ..., v_k) = -w(v_2, v_1, ..., v_k)
   - k-linear: w acts linear in each vector separately, e.g. w(a*v + b*y, v_2, ..., v_k) = a*w(v,v_2,...,v_k) + b*w(y,v_2,...,v_k)
   This action is implemented as PetscDTAltVApply.

   The k-forms on a vector space form a vector space themselves, Alt^k V.  The dimension of Alt^k V, if V is N dimensional, is N choose k.  (This
   shows that for an N dimensional space, only 0 <= k <= N are valid form degrees.)
   The standard basis for Alt^k V, used in PetscDTAltV, has one basis k-form for each ordered subset of k coordinates of the N dimensional space:
   For example, if the coordinate directions of a four dimensional space are (t, x, y, z), then there are 4 choose 2 = 6 ordered subsets of two coordinates.
   They are, in lexicographic order, (t, x), (t, y), (t, z), (x, y), (x, z) and (y, z).  PetscDTAltV also orders the basis of Alt^k V lexicographically
   by the associated subsets.

   The unit basis k-form associated with coordinates (c_1, ..., c_k) acts on a set of k vectors (v_1, ..., v_k) by creating a square matrix V where
   V[i,j] = v_i[c_j] and taking the determinant of V.

   If j + k <= N, then a j-form f and a k-form g can be multiplied to create a (j+k)-form using the wedge or exterior product, (f wedge g).
   This is an anticommutative product, (f wedge g) = -(g wedge f).  It is sufficient to describe the wedge product of two basis forms.
   Let f be the basis j-form associated with coordinates (f_1,...,f_j) and g be the basis k-form associated with coordinates (g_1,...,g_k):
   - If there is any coordinate in both sets, then (f wedge g) = 0.
   - Otherwise, (f wedge g) is a multiple of the basis (j+k)-form h associated with (f_1,...,f_j,g_1,...,g_k).
   - In fact it is equal to either h or -h depending on how (f_1,...,f_j,g_1,...,g_k) compares to the same list of coordinates given in ascending order: if it is an even permutation of that list, then (f wedge g) = h, otherwise (f wedge g) = -h.
   The wedge product is implemented for either two inputs (f and g) in PetscDTAltVWedge, or for one (just f, giving a
   matrix to multiply against multiple choices of g) in PetscDTAltVWedgeMatrix.

   If k > 0, a k-form w and a vector v can combine to make a (k-1)-formm through the interior product, (w int v),
   defined by (w int v)(v_1,...,v_{k-1}) = w(v,v_1,...,v_{k-1}).

   The interior product is implemented for either two inputs (w and v) in PetscDTAltVInterior, for one (just v, giving a
   matrix to multiply against multiple choices of w) in PetscDTAltVInteriorMatrix,
   or for no inputs (giving the sparsity pattern of PetscDTAltVInteriorMatrix) in PetscDTAltVInteriorPattern.

   When there is a linear map L: V -> W from an N dimensional vector space to an M dimensional vector space,
   it induces the linear pullback map L^* : Alt^k W -> Alt^k V, defined by L^* w(v_1,...,v_k) = w(L v_1, ..., L v_k).
   The pullback is implemented as PetscDTAltVPullback (acting on a known w) and PetscDTAltVPullbackMatrix (creating a matrix that computes the actin of L^*).

   Alt^k V and Alt^(N-k) V have the same dimension, and the Hodge star operator maps between them.  We note that Alt^N V is a one dimensional space, and its
   basis vector is sometime called vol.  The Hodge star operator has the property that (f wedge (star g)) = (f,g) vol, where (f,g) is the simple inner product
   of the basis coefficients of f and g.
   Powers of the Hodge star operator can be applied with PetscDTAltVStar

   Level: intermediate

.seealso: PetscDTAltVApply(), PetscDTAltVWedge(), PetscDTAltVInterior(), PetscDTAltVPullback(), PetscDTAltVStar()
M*/

/*@
   PetscDTAltVApply - Apply an a k-form (an alternating k-linear map) to a set of k N-dimensional vectors

   Input Arguments:
+  N - the dimension of the vector space, N >= 0
.  k - the degree k of the k-form w, 0 <= k <= N
.  w - a k-form, size [N choose k] (each degree of freedom of a k-form is associated with a subset of k coordinates of the N-dimensional vectors: the degrees of freedom are ordered lexicographically by their associated subsets)
-  v - a set of k vectors of size N, size [k x N], each vector stored contiguously

   Output Arguments:
.  wv - w(v_1,...,v_k) = \sum_i w_i * det(V_i): the degree of freedom w_i is associated with coordinates [s_{i,1},...,s_{i,k}], and the square matrix V_i has entry (j,k) given by the s_{i,k}'th coordinate of v_j

   Level: intermediate

.seealso: PetscDTAltV, PetscDTAltVPullback(), PetscDTAltVPullbackMatrix()
@*/
PetscErrorCode PetscDTAltVApply(PetscInt N, PetscInt k, const PetscReal *w, const PetscReal *v, PetscReal *wv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (N < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid dimension");
  if (k < 0 || k > N) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid form degree");
  if (N <= 3) {
    if (!k) {
      *wv = w[0];
    } else {
      if (N == 1)        {*wv = w[0] * v[0];}
      else if (N == 2) {
        if (k == 1)      {*wv = w[0] * v[0] + w[1] * v[1];}
        else             {*wv = w[0] * (v[0] * v[3] - v[1] * v[2]);}
      } else {
        if (k == 1)      {*wv = w[0] * v[0] + w[1] * v[1] + w[2] * v[2];}
        else if (k == 2) {
          *wv = w[0] * (v[0] * v[4] - v[1] * v[3]) +
                w[1] * (v[0] * v[5] - v[2] * v[3]) +
                w[2] * (v[1] * v[5] - v[2] * v[4]);
        } else {
          *wv = w[0] * (v[0] * (v[4] * v[8] - v[5] * v[7]) +
                        v[1] * (v[5] * v[6] - v[3] * v[8]) +
                        v[2] * (v[3] * v[7] - v[4] * v[6]));
        }
      }
    }
  } else {
    PetscInt Nk, Nf;
    PetscInt *subset, *perm;
    PetscInt i, j, l;
    PetscReal sum = 0.;

    ierr = PetscDTFactorialInt(k, &Nf);CHKERRQ(ierr);
    ierr = PetscDTBinomialInt(N, k, &Nk);CHKERRQ(ierr);
    ierr = PetscMalloc2(k, &subset, k, &perm);CHKERRQ(ierr);
    for (i = 0; i < Nk; i++) {
      PetscReal subsum = 0.;

      ierr = PetscDTEnumSubset(N, k, i, subset);CHKERRQ(ierr);
      for (j = 0; j < Nf; j++) {
        PetscBool permOdd;
        PetscReal prod;

        ierr = PetscDTEnumPerm(k, j, perm, &permOdd);CHKERRQ(ierr);
        prod = permOdd ? -1. : 1.;
        for (l = 0; l < k; l++) {
          prod *= v[perm[l] * N + subset[l]];
        }
        subsum += prod;
      }
      sum += w[i] * subsum;
    }
    ierr = PetscFree2(subset, perm);CHKERRQ(ierr);
    *wv = sum;
  }
  PetscFunctionReturn(0);
}

/*@
   PetscDTAltVWedge - Compute the wedge product of a j-form and a k-form, giving a (j+k) form

   Input Arguments:
+  N - the dimension of the vector space, N >= 0
.  j - the degree j of the j-form a, 0 <= j <= N
.  k - the degree k of the k-form b, 0 <= k <= N and 0 <= j+k <= N
.  a - a j-form, size [N choose j]
-  b - a k-form, size [N choose k]

   Output Arguments:
.  awedgeb - the (j+k)-form a wedge b, size [N choose (j+k)]: (a wedge b)(v_1,...,v_{j+k}) = \sum_{s} sign(s) a(v_{s_1},...,v_{s_j}) b(v_{s_{j+1}},...,v_{s_{j+k}}),
             where the sum is over permutations s such that s_1 < s_2 < ... < s_j and s_{j+1} < s_{j+2} < ... < s_{j+k}.

   Level: intermediate

.seealso: PetscDTAltV, PetscDTAltVWedgeMatrix(), PetscDTAltVPullback(), PetscDTAltVPullbackMatrix()
@*/
PetscErrorCode PetscDTAltVWedge(PetscInt N, PetscInt j, PetscInt k, const PetscReal *a, const PetscReal *b, PetscReal *awedgeb)
{
  PetscInt       i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (N < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid dimension");
  if (j < 0 || k < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "negative form degree");
  if (j + k > N) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Wedge greater than dimension");
  if (N <= 3) {
    PetscInt Njk;

    ierr = PetscDTBinomialInt(N, j+k, &Njk);CHKERRQ(ierr);
    if (!j)      {for (i = 0; i < Njk; i++) {awedgeb[i] = a[0] * b[i];}}
    else if (!k) {for (i = 0; i < Njk; i++) {awedgeb[i] = a[i] * b[0];}}
    else {
      if (N == 2) {awedgeb[0] = a[0] * b[1] - a[1] * b[0];}
      else {
        if (j+k == 2) {
          awedgeb[0] = a[0] * b[1] - a[1] * b[0];
          awedgeb[1] = a[0] * b[2] - a[2] * b[0];
          awedgeb[2] = a[1] * b[2] - a[2] * b[1];
        } else {
          awedgeb[0] = a[0] * b[2] - a[1] * b[1] + a[2] * b[0];
        }
      }
    }
  } else {
    PetscInt  Njk;
    PetscInt  JKj;
    PetscInt *subset, *subsetjk, *subsetj, *subsetk;
    PetscInt  i;

    ierr = PetscDTBinomialInt(N, j+k, &Njk);CHKERRQ(ierr);
    ierr = PetscDTBinomialInt(j+k, j, &JKj);CHKERRQ(ierr);
    ierr = PetscMalloc4(j+k, &subset, j+k, &subsetjk, j, &subsetj, k, &subsetk);CHKERRQ(ierr);
    for (i = 0; i < Njk; i++) {
      PetscReal sum = 0.;
      PetscInt  l;

      ierr = PetscDTEnumSubset(N, j+k, i, subset);CHKERRQ(ierr);
      for (l = 0; l < JKj; l++) {
        PetscBool jkOdd;
        PetscInt  m, jInd, kInd;

        ierr = PetscDTEnumSplit(j+k, j, l, subsetjk, &jkOdd);CHKERRQ(ierr);
        for (m = 0; m < j; m++) {
          subsetj[m] = subset[subsetjk[m]];
        }
        for (m = 0; m < k; m++) {
          subsetk[m] = subset[subsetjk[j+m]];
        }
        ierr = PetscDTSubsetIndex(N, j, subsetj, &jInd);CHKERRQ(ierr);
        ierr = PetscDTSubsetIndex(N, k, subsetk, &kInd);CHKERRQ(ierr);
        sum += jkOdd ? -(a[jInd] * b[kInd]) : (a[jInd] * b[kInd]);
      }
      awedgeb[i] = sum;
    }
    ierr = PetscFree4(subset, subsetjk, subsetj, subsetk);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PetscDTAltVWedgeMatrix - Compute the matrix defined by the wedge product with a given j-form that maps k-forms to (j+k)-forms

   Input Arguments:
+  N - the dimension of the vector space, N >= 0
.  j - the degree j of the j-form a, 0 <= j <= N
.  k - the degree k of the k-forms that (a wedge) will be applied to, 0 <= k <= N and 0 <= j+k <= N
-  a - a j-form, size [N choose j]

   Output Arguments:
.  awedge - (a wedge), an [(N choose j+k) x (N choose k)] matrix in row-major order, such that (a wedge) * b = a wedge b

   Level: intermediate

.seealso: PetscDTAltV, PetscDTAltVPullback(), PetscDTAltVPullbackMatrix()
@*/
PetscErrorCode PetscDTAltVWedgeMatrix(PetscInt N, PetscInt j, PetscInt k, const PetscReal *a, PetscReal *awedgeMat)
{
  PetscInt       i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (N < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid dimension");
  if (j < 0 || k < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "negative form degree");
  if (j + k > N) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Wedge greater than dimension");
  if (N <= 3) {
    PetscInt Njk;

    ierr = PetscDTBinomialInt(N, j+k, &Njk);CHKERRQ(ierr);
    if (!j) {
      for (i = 0; i < Njk * Njk; i++) {awedgeMat[i] = 0.;}
      for (i = 0; i < Njk; i++) {awedgeMat[i * (Njk + 1)] = a[0];}
    } else if (!k) {
      for (i = 0; i < Njk; i++) {awedgeMat[i] = a[i];}
    } else {
      if (N == 2) {
        awedgeMat[0] = -a[1]; awedgeMat[1] =  a[0];
      } else {
        if (j+k == 2) {
          awedgeMat[0] = -a[1]; awedgeMat[1] =  a[0]; awedgeMat[2] =    0.;
          awedgeMat[3] = -a[2]; awedgeMat[4] =    0.; awedgeMat[5] =  a[0];
          awedgeMat[6] =    0.; awedgeMat[7] = -a[2]; awedgeMat[8] =  a[1];
        } else {
          awedgeMat[0] =  a[2]; awedgeMat[1] = -a[1]; awedgeMat[2] =  a[0];
        }
      }
    }
  } else {
    PetscInt  Njk;
    PetscInt  Nk;
    PetscInt  JKj, i;
    PetscInt *subset, *subsetjk, *subsetj, *subsetk;

    ierr = PetscDTBinomialInt(N,   k,   &Nk);CHKERRQ(ierr);
    ierr = PetscDTBinomialInt(N,   j+k, &Njk);CHKERRQ(ierr);
    ierr = PetscDTBinomialInt(j+k, j,   &JKj);CHKERRQ(ierr);
    ierr = PetscMalloc4(j+k, &subset, j+k, &subsetjk, j, &subsetj, k, &subsetk);CHKERRQ(ierr);
    for (i = 0; i < Njk * Nk; i++) awedgeMat[i] = 0.;
    for (i = 0; i < Njk; i++) {
      PetscInt  l;

      ierr = PetscDTEnumSubset(N, j+k, i, subset);CHKERRQ(ierr);
      for (l = 0; l < JKj; l++) {
        PetscBool jkOdd;
        PetscInt  m, jInd, kInd;

        ierr = PetscDTEnumSplit(j+k, j, l, subsetjk, &jkOdd);CHKERRQ(ierr);
        for (m = 0; m < j; m++) {
          subsetj[m] = subset[subsetjk[m]];
        }
        for (m = 0; m < k; m++) {
          subsetk[m] = subset[subsetjk[j+m]];
        }
        ierr = PetscDTSubsetIndex(N, j, subsetj, &jInd);CHKERRQ(ierr);
        ierr = PetscDTSubsetIndex(N, k, subsetk, &kInd);CHKERRQ(ierr);
        awedgeMat[i * Nk + kInd] += jkOdd ? - a[jInd] : a[jInd];
      }
    }
    ierr = PetscFree4(subset, subsetjk, subsetj, subsetk);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PetscDTAltVPullback - Compute the pullback of a k-form under a linear transformation of the coordinate space

   Input Arguments:
+  N - the dimension of the origin vector space of the linear transformation, M >= 0
.  M - the dimension of the image vector space of the linear transformation, N >= 0
.  L - a linear transformation, an [M x N] matrix in row-major format
.  k - the *signed* degree k of the |k|-form w, -(min(M,N)) <= k <= min(M,N).  A negative form degree indicates that the pullback should be conjugated by the Hodge star operator (see note).
-  w - a |k|-form in the image space, size [M choose |k|]

   Output Arguments:
.  Lstarw - the pullback of w to a |k|-form in the origin space, size [N choose |k|]: (Lstarw)(v_1,...v_k) = w(L*v_1,...,L*v_k).

   Level: intermediate

   Note: negative form degrees accomodate, e.g., H-div conforming vector fields.  An H-div conforming vector field stores its degrees of freedom as (dx, dy, dz), like a 1-form,
   but its normal trace is integrated on faces, like a 2-form.  The correct pullback then is to apply the Hodge star transformation from (M-2)-form to 2-form, pullback as a 2-form,
   then the inverse Hodge star transformation.

.seealso: PetscDTAltV, PetscDTAltVPullbackMatrix(), PetscDTAltVStar()
@*/
PetscErrorCode PetscDTAltVPullback(PetscInt N, PetscInt M, const PetscReal *L, PetscInt k, const PetscReal *w, PetscReal *Lstarw)
{
  PetscInt         i, j, Nk, Mk;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  if (N < 0 || M < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid dimensions");
  if (PetscAbsInt(k) > N || PetscAbsInt(k) > M) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid form degree");
  if (N <= 3 && M <= 3) {

    ierr = PetscDTBinomialInt(M, PetscAbsInt(k), &Mk);CHKERRQ(ierr);
    ierr = PetscDTBinomialInt(N, PetscAbsInt(k), &Nk);CHKERRQ(ierr);
    if (!k) {
      Lstarw[0] = w[0];
    } else if (k == 1) {
      for (i = 0; i < Nk; i++) {
        PetscReal sum = 0.;

        for (j = 0; j < Mk; j++) {sum += L[j * Nk + i] * w[j];}
        Lstarw[i] = sum;
      }
    } else if (k == -1) {
      PetscReal mult[3] = {1., -1., 1.};

      for (i = 0; i < Nk; i++) {
        PetscReal sum = 0.;

        for (j = 0; j < Mk; j++) {
          sum += L[(Mk - 1 - j) * Nk + (Nk - 1 - i)] * w[j] * mult[j];
        }
        Lstarw[i] = mult[i] * sum;
      }
    } else if (k == 2) {
      PetscInt pairs[3][2] = {{0,1},{0,2},{1,2}};

      for (i = 0; i < Nk; i++) {
        PetscReal sum = 0.;
        for (j = 0; j < Mk; j++) {
          sum += (L[pairs[j][0] * N + pairs[i][0]] * L[pairs[j][1] * N + pairs[i][1]] -
                  L[pairs[j][1] * N + pairs[i][0]] * L[pairs[j][0] * N + pairs[i][1]]) * w[j];
        }
        Lstarw[i] = sum;
      }
    } else if (k == -2) {
      PetscInt  pairs[3][2] = {{1,2},{2,0},{0,1}};
      PetscInt  offi = (N == 2) ? 2 : 0;
      PetscInt  offj = (M == 2) ? 2 : 0;

      for (i = 0; i < Nk; i++) {
        PetscReal sum   = 0.;

        for (j = 0; j < Mk; j++) {
          sum += (L[pairs[offj + j][0] * N + pairs[offi + i][0]] *
                  L[pairs[offj + j][1] * N + pairs[offi + i][1]] -
                  L[pairs[offj + j][1] * N + pairs[offi + i][0]] *
                  L[pairs[offj + j][0] * N + pairs[offi + i][1]]) * w[j];

        }
        Lstarw[i] = sum;
      }
    } else {
      PetscReal detL = L[0] * (L[4] * L[8] - L[5] * L[7]) +
                       L[1] * (L[5] * L[6] - L[3] * L[8]) +
                       L[2] * (L[3] * L[7] - L[4] * L[6]);

      for (i = 0; i < Nk; i++) {Lstarw[i] = detL * w[i];}
    }
  } else {
    PetscInt         Nf, l, p;
    PetscReal       *Lw, *Lwv;
    PetscInt        *subsetw, *subsetv;
    PetscInt        *perm;
    PetscReal       *walloc = NULL;
    const PetscReal *ww = NULL;
    PetscBool        negative = PETSC_FALSE;

    ierr = PetscDTBinomialInt(M, PetscAbsInt(k), &Mk);CHKERRQ(ierr);
    ierr = PetscDTBinomialInt(N, PetscAbsInt(k), &Nk);CHKERRQ(ierr);
    ierr = PetscDTFactorialInt(PetscAbsInt(k), &Nf);CHKERRQ(ierr);
    if (k < 0) {
      negative = PETSC_TRUE;
      k = -k;
      ierr = PetscMalloc1(Mk, &walloc);CHKERRQ(ierr);
      ierr = PetscDTAltVStar(M, M - k, 1, w, walloc);CHKERRQ(ierr);
      ww = walloc;
    } else {
      ww = w;
    }
    ierr = PetscMalloc5(k, &subsetw, k, &subsetv, k, &perm, N * k, &Lw, k * k, &Lwv);CHKERRQ(ierr);
    for (i = 0; i < Nk; i++) Lstarw[i] = 0.;
    for (i = 0; i < Mk; i++) {
      ierr = PetscDTEnumSubset(M, k, i, subsetw);CHKERRQ(ierr);
      for (j = 0; j < Nk; j++) {
        ierr = PetscDTEnumSubset(N, k, j, subsetv);CHKERRQ(ierr);
        for (p = 0; p < Nf; p++) {
          PetscReal prod;
          PetscBool isOdd;

          ierr = PetscDTEnumPerm(k, p, perm, &isOdd);CHKERRQ(ierr);
          prod = isOdd ? -ww[i] : ww[i];
          for (l = 0; l < k; l++) {
            prod *= L[subsetw[perm[l]] * N + subsetv[l]];
          }
          Lstarw[j] += prod;
        }
      }
    }
    if (negative) {
      PetscReal *sLsw;

      ierr = PetscMalloc1(Nk, &sLsw);CHKERRQ(ierr);
      ierr = PetscDTAltVStar(N, N - k, -1,  Lstarw, sLsw);CHKERRQ(ierr);
      for (i = 0; i < Nk; i++) Lstarw[i] = sLsw[i];
      ierr = PetscFree(sLsw);CHKERRQ(ierr);
    }
    ierr = PetscFree5(subsetw, subsetv, perm, Lw, Lwv);CHKERRQ(ierr);
    ierr = PetscFree(walloc);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PetscDTAltVPullbackMatrix - Compute the pullback matrix for k-forms under a linear transformation

   Input Arguments:
+  N - the dimension of the origin vector space of the linear transformation, N >= 0
.  M - the dimension of the image vector space of the linear transformation, M >= 0
.  L - a linear transformation, an [M x N] matrix in row-major format
-  k - the *signed* degree k of the |k|-forms on which Lstar acts, -(min(M,N)) <= k <= min(M,N).  A negative form degree indicates that the pullback should be conjugated by the Hodge star operator (see note in PetscDTAltvPullback())

   Output Arguments:
.  Lstar - the pullback matrix, an [(N choose |k|) x (M choose |k|)] matrix in row-major format such that Lstar * w = L^* w

   Level: intermediate

.seealso: PetscDTAltV, PetscDTAltVPullback(), PetscDTAltVStar()
@*/
PetscErrorCode PetscDTAltVPullbackMatrix(PetscInt N, PetscInt M, const PetscReal *L, PetscInt k, PetscReal *Lstar)
{
  PetscInt        Nk, Mk, Nf, i, j, l, p;
  PetscReal      *Lw, *Lwv;
  PetscInt       *subsetw, *subsetv;
  PetscInt       *perm;
  PetscBool       negative = PETSC_FALSE;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (N < 0 || M < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid dimensions");
  if (PetscAbsInt(k) > N || PetscAbsInt(k) > M) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid form degree");
  if (N <= 3 && M <= 3) {
    PetscReal mult[3] = {1., -1., 1.};

    ierr = PetscDTBinomialInt(M, PetscAbsInt(k), &Mk);CHKERRQ(ierr);
    ierr = PetscDTBinomialInt(N, PetscAbsInt(k), &Nk);CHKERRQ(ierr);
    if (!k) {
      Lstar[0] = 1.;
    } else if (k == 1) {
      for (i = 0; i < Nk; i++) {for (j = 0; j < Mk; j++) {Lstar[i * Mk + j] = L[j * Nk + i];}}
    } else if (k == -1) {
      for (i = 0; i < Nk; i++) {
        for (j = 0; j < Mk; j++) {
          Lstar[i * Mk + j] = L[(Mk - 1 - j) * Nk + (Nk - 1 - i)] * mult[i] * mult[j];
        }
      }
    } else if (k == 2) {
      PetscInt pairs[3][2] = {{0,1},{0,2},{1,2}};

      for (i = 0; i < Nk; i++) {
        for (j = 0; j < Mk; j++) {
          Lstar[i * Mk + j] = L[pairs[j][0] * N + pairs[i][0]] *
                              L[pairs[j][1] * N + pairs[i][1]] -
                              L[pairs[j][1] * N + pairs[i][0]] *
                              L[pairs[j][0] * N + pairs[i][1]];
        }
      }
    } else if (k == -2) {
      PetscInt  pairs[3][2] = {{1,2},{2,0},{0,1}};
      PetscInt  offi = (N == 2) ? 2 : 0;
      PetscInt  offj = (M == 2) ? 2 : 0;

      for (i = 0; i < Nk; i++) {
        for (j = 0; j < Mk; j++) {
          Lstar[i * Mk + j] = L[pairs[offj + j][0] * N + pairs[offi + i][0]] *
                              L[pairs[offj + j][1] * N + pairs[offi + i][1]] -
                              L[pairs[offj + j][1] * N + pairs[offi + i][0]] *
                              L[pairs[offj + j][0] * N + pairs[offi + i][1]];
        }
      }
    } else {
      PetscReal detL = L[0] * (L[4] * L[8] - L[5] * L[7]) +
                       L[1] * (L[5] * L[6] - L[3] * L[8]) +
                       L[2] * (L[3] * L[7] - L[4] * L[6]);

      for (i = 0; i < Nk; i++) {Lstar[i] = detL;}
    }
  } else {
    if (k < 0) {
      negative = PETSC_TRUE;
      k = -k;
    }
    ierr = PetscDTBinomialInt(M, PetscAbsInt(k), &Mk);CHKERRQ(ierr);
    ierr = PetscDTBinomialInt(N, PetscAbsInt(k), &Nk);CHKERRQ(ierr);
    ierr = PetscDTFactorialInt(PetscAbsInt(k), &Nf);CHKERRQ(ierr);
    ierr = PetscMalloc5(M, &subsetw, N, &subsetv, k, &perm, N * k, &Lw, k * k, &Lwv);CHKERRQ(ierr);
    for (i = 0; i < Nk * Mk; i++) Lstar[i] = 0.;
    for (i = 0; i < Mk; i++) {
      PetscBool iOdd;
      PetscInt  iidx, jidx;

      ierr = PetscDTEnumSplit(M, k, i, subsetw, &iOdd);CHKERRQ(ierr);
      iidx = negative ? Mk - 1 - i : i;
      iOdd = negative ? (PetscBool) (iOdd ^ ((k * (M-k)) & 1)) : PETSC_FALSE;
      for (j = 0; j < Nk; j++) {
        PetscBool jOdd;

        ierr = PetscDTEnumSplit(N, k, j, subsetv, &jOdd);CHKERRQ(ierr);
        jidx = negative ? Nk - 1 - j : j;
        jOdd = negative ? (PetscBool) (iOdd ^ jOdd ^ ((k * (N-k)) & 1)) : PETSC_FALSE;
        for (p = 0; p < Nf; p++) {
          PetscReal prod;
          PetscBool isOdd;

          ierr = PetscDTEnumPerm(k, p, perm, &isOdd);CHKERRQ(ierr);
          isOdd = (PetscBool) (isOdd ^ jOdd);
          prod = isOdd ? -1. : 1.;
          for (l = 0; l < k; l++) {
            prod *= L[subsetw[perm[l]] * N + subsetv[l]];
          }
          Lstar[jidx * Mk + iidx] += prod;
        }
      }
    }
    ierr = PetscFree5(subsetw, subsetv, perm, Lw, Lwv);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PetscDTAltVInterior - Compute the interior product of a k-form with a vector

   Input Arguments:
+  N - the dimension of the vector space, N >= 0
.  k - the degree k of the k-form w, 0 <= k <= N
.  w - a k-form, size [N choose k]
-  v - an N dimensional vector

   Output Arguments:
.  wIntv - the (k-1)-form (w int v), size [N choose (k-1)]: (w int v) is defined by its action on (k-1) vectors {v_1, ..., v_{k-1}} as (w inv v)(v_1, ..., v_{k-1}) = w(v, v_1, ..., v_{k-1}).

   Level: intermediate

.seealso: PetscDTAltV, PetscDTAltVInteriorMatrix(), PetscDTAltVInteriorPattern(), PetscDTAltVPullback(), PetscDTAltVPullbackMatrix()
@*/
PetscErrorCode PetscDTAltVInterior(PetscInt N, PetscInt k, const PetscReal *w, const PetscReal *v, PetscReal *wIntv)
{
  PetscInt        i, Nk, Nkm;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (k <= 0 || k > N) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid form degree");
  ierr = PetscDTBinomialInt(N, k,   &Nk);CHKERRQ(ierr);
  ierr = PetscDTBinomialInt(N, k-1, &Nkm);CHKERRQ(ierr);
  if (N <= 3) {
    if (k == 1) {
      PetscReal sum = 0.;

      for (i = 0; i < N; i++) {
        sum += w[i] * v[i];
      }
      wIntv[0] = sum;
    } else if (k == N) {
      PetscReal mult[3] = {1., -1., 1.};

      for (i = 0; i < N; i++) {
        wIntv[N - 1 - i] = w[0] * v[i] * mult[i];
      }
    } else {
      wIntv[0] = - w[0]*v[1] - w[1]*v[2];
      wIntv[1] =   w[0]*v[0] - w[2]*v[2];
      wIntv[2] =   w[1]*v[0] + w[2]*v[1];
    }
  } else {
    PetscInt       *subset, *work;

    ierr = PetscMalloc2(k, &subset, k, &work);CHKERRQ(ierr);
    for (i = 0; i < Nkm; i++) wIntv[i] = 0.;
    for (i = 0; i < Nk; i++) {
      PetscInt  j, l, m;

      ierr = PetscDTEnumSubset(N, k, i, subset);CHKERRQ(ierr);
      for (j = 0; j < k; j++) {
        PetscInt  idx;
        PetscBool flip = (PetscBool) (j & 1);

        for (l = 0, m = 0; l < k; l++) {
          if (l != j) work[m++] = subset[l];
        }
        ierr = PetscDTSubsetIndex(N, k - 1, work, &idx);CHKERRQ(ierr);
        wIntv[idx] += flip ? -(w[i] * v[subset[j]]) :  (w[i] * v[subset[j]]);
      }
    }
    ierr = PetscFree2(subset, work);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PetscDTAltVInteriorMatrix - Compute the matrix of the linear transformation induced on a k-form by the interior product with a vector

   Input Arguments:
+  N - the dimension of the vector space, N >= 0
.  k - the degree k of the k-forms on which intvMat acts, 0 <= k <= N
-  v - an N dimensional vector

   Output Arguments:
.  intvMat - an [(N choose (k-1)) x (N choose k)] matrix, row-major: (intvMat) * w = (w int v)

   Level: intermediate

.seealso: PetscDTAltV, PetscDTAltVInterior(), PetscDTAltVInteriorPattern(), PetscDTAltVPullback(), PetscDTAltVPullbackMatrix()
@*/
PetscErrorCode PetscDTAltVInteriorMatrix(PetscInt N, PetscInt k, const PetscReal *v, PetscReal *intvMat)
{
  PetscInt        i, Nk, Nkm;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (k <= 0 || k > N) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid form degree");
  ierr = PetscDTBinomialInt(N, k,   &Nk);CHKERRQ(ierr);
  ierr = PetscDTBinomialInt(N, k-1, &Nkm);CHKERRQ(ierr);
  if (N <= 3) {
    if (k == 1) {
      for (i = 0; i < N; i++) intvMat[i] = v[i];
    } else if (k == N) {
      PetscReal mult[3] = {1., -1., 1.};

      for (i = 0; i < N; i++) intvMat[N - 1 - i] = v[i] * mult[i];
    } else {
      intvMat[0] = -v[1]; intvMat[1] = -v[2]; intvMat[2] =    0.;
      intvMat[3] =  v[0]; intvMat[4] =    0.; intvMat[5] = -v[2];
      intvMat[6] =    0.; intvMat[7] =  v[0]; intvMat[8] =  v[1];
    }
  } else {
    PetscInt       *subset, *work;

    ierr = PetscMalloc2(k, &subset, k, &work);CHKERRQ(ierr);
    for (i = 0; i < Nk * Nkm; i++) intvMat[i] = 0.;
    for (i = 0; i < Nk; i++) {
      PetscInt  j, l, m;

      ierr = PetscDTEnumSubset(N, k, i, subset);CHKERRQ(ierr);
      for (j = 0; j < k; j++) {
        PetscInt  idx;
        PetscBool flip = (PetscBool) (j & 1);

        for (l = 0, m = 0; l < k; l++) {
          if (l != j) work[m++] = subset[l];
        }
        ierr = PetscDTSubsetIndex(N, k - 1, work, &idx);CHKERRQ(ierr);
        intvMat[idx * Nk + i] += flip ? -v[subset[j]] :  v[subset[j]];
      }
    }
    ierr = PetscFree2(subset, work);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PetscDTAltVInteriorPattern - compute the sparsity and sign pattern of the interior product matrix computed in PetscDTAltVInteriorMatrix()

   Input Arguments:
+  N - the dimension of the vector space, N >= 0
-  k - the degree of the k-forms on which intvMat from PetscDTAltVInteriorMatrix() acts, 0 <= k <= N.

   Output Arguments:
.  indices - The interior product matrix intvMat has size [(N choose (k-1)) x (N choose k)] and has (N choose k) * k
             non-zeros.  indices[i][0] and indices[i][1] are the row and column of a non-zero, and its value is equal to the vector
             coordinate v[j] if indices[i][2] = j, or -v[j] if indices[i][2] = -(j+1)

   Level: intermediate

   Note: this function is useful when the interior product needs to be computed at multiple locations, as when computing the Koszul differential

.seealso: PetscDTAltV, PetscDTAltVInterior(), PetscDTAltVInteriorMatrix(), PetscDTAltVPullback(), PetscDTAltVPullbackMatrix()
@*/
PetscErrorCode PetscDTAltVInteriorPattern(PetscInt N, PetscInt k, PetscInt (*indices)[3])
{
  PetscInt        i, Nk, Nkm;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (k <= 0 || k > N) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid form degree");
  ierr = PetscDTBinomialInt(N, k,   &Nk);CHKERRQ(ierr);
  ierr = PetscDTBinomialInt(N, k-1, &Nkm);CHKERRQ(ierr);
  if (N <= 3) {
    if (k == 1) {
      for (i = 0; i < N; i++) {
        indices[i][0] = 0;
        indices[i][1] = i;
        indices[i][2] = i;
      }
    } else if (k == N) {
      PetscInt val[3] = {0, -2, 2};

      for (i = 0; i < N; i++) {
        indices[i][0] = N - 1 - i;
        indices[i][1] = 0;
        indices[i][2] = val[i];
      }
    } else {
      indices[0][0] = 0; indices[0][1] = 0; indices[0][2] = -(1 + 1);
      indices[1][0] = 0; indices[1][1] = 1; indices[1][2] = -(2 + 1);
      indices[2][0] = 1; indices[2][1] = 0; indices[2][2] = 0;
      indices[3][0] = 1; indices[3][1] = 2; indices[3][2] = -(2 + 1);
      indices[4][0] = 2; indices[4][1] = 1; indices[4][2] = 0;
      indices[5][0] = 2; indices[5][1] = 2; indices[5][2] = 1;
    }
  } else {
    PetscInt       *subset, *work;

    ierr = PetscMalloc2(k, &subset, k, &work);CHKERRQ(ierr);
    for (i = 0; i < Nk; i++) {
      PetscInt  j, l, m;

      ierr = PetscDTEnumSubset(N, k, i, subset);CHKERRQ(ierr);
      for (j = 0; j < k; j++) {
        PetscInt  idx;
        PetscBool flip = (PetscBool) (j & 1);

        for (l = 0, m = 0; l < k; l++) {
          if (l != j) work[m++] = subset[l];
        }
        ierr = PetscDTSubsetIndex(N, k - 1, work, &idx);CHKERRQ(ierr);
        indices[i * k + j][0] = idx;
        indices[i * k + j][1] = i;
        indices[i * k + j][2] = flip ? -(subset[j] + 1) : subset[j];
      }
    }
    ierr = PetscFree2(subset, work);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*@
   PetscDTAltVStar - Apply a power of the Hodge star operator, which maps k-forms to (N-k) forms, to a k-form

   Input Arguments:
+  N - the dimension of the vector space, N >= 0
.  k - the degree k of the k-form w, 0 <= k <= N
.  pow - the number of times to apply the Hodge star operator: pow < 0 indicates that the inverse of the Hodge star operator should be applied |pow| times.
-  w - a k-form, size [N choose k]

   Output Arguments:
.  starw = (star)^pow w.  Each degree of freedom of a k-form is associated with a subset S of k coordinates of the N dimensional vector space: the Hodge start operator (star) maps that degree of freedom to the degree of freedom associated with S', the complement of S, with a sign change if the permutation of coordinates {S[0], ... S[k-1], S'[0], ... S'[N-k- 1]} is an odd permutation.  This implies (star)^2 w = (-1)^{k(N-k)} w, and (star)^4 w = w.

   Level: intermediate

.seealso: PetscDTAltV, PetscDTAltVPullback(), PetscDTAltVPullbackMatrix()
@*/
PetscErrorCode PetscDTAltVStar(PetscInt N, PetscInt k, PetscInt pow, const PetscReal *w, PetscReal *starw)
{
  PetscInt        Nk, i;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (k < 0 || k > N) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "invalid form degree");
  ierr = PetscDTBinomialInt(N, k, &Nk);CHKERRQ(ierr);
  pow = pow % 4;
  pow = (pow + 4) % 4; /* make non-negative */
  /* pow is now 0, 1, 2, 3 */
  if (N <= 3) {
    if (pow & 1) {
      PetscReal mult[3] = {1., -1., 1.};

      for (i = 0; i < Nk; i++) starw[Nk - 1 - i] = w[i] * mult[i];
    } else {
      for (i = 0; i < Nk; i++) starw[i] = w[i];
    }
    if (pow > 1 && ((k * (N - k)) & 1)) {
      for (i = 0; i < Nk; i++) starw[i] = -starw[i];
    }
  } else {
    PetscInt       *subset;

    ierr = PetscMalloc1(N, &subset);CHKERRQ(ierr);
    if (pow % 2) {
      PetscInt l = (pow == 1) ? k : N - k;
      for (i = 0; i < Nk; i++) {
        PetscBool sOdd;
        PetscInt  j, idx;

        ierr = PetscDTEnumSplit(N, l, i, subset, &sOdd);CHKERRQ(ierr);
        ierr = PetscDTSubsetIndex(N, l, subset, &idx);CHKERRQ(ierr);
        ierr = PetscDTSubsetIndex(N, N-l, &subset[l], &j);CHKERRQ(ierr);
        starw[j] = sOdd ? -w[idx] : w[idx];
      }
    } else {
      for (i = 0; i < Nk; i++) starw[i] = w[i];
    }
    /* star^2 = -1^(k * (N - k)) */
    if (pow > 1 && (k * (N - k)) % 2) {
      for (i = 0; i < Nk; i++) starw[i] = -starw[i];
    }
    ierr = PetscFree(subset);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
