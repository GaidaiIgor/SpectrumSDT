/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc polynomial eigensolver: "stoar"

   Method: S-TOAR with spectrum slicing for symmetric quadratic eigenproblems

   Algorithm:

       Symmetric Two-Level Orthogonal Arnoldi.

   References:

       [1] C. Campos and J.E. Roman, "Inertia-based spectrum slicing
           for symmetric quadratic eigenvalue problems", Numer. Linear
           Algebra Appl. (in press), 2020.
*/

#include <slepc/private/pepimpl.h>         /*I "slepcpep.h" I*/
#include "../src/pep/impls/krylov/pepkrylov.h"
#include <slepcblaslapack.h>

static PetscBool  cited = PETSC_FALSE;
static const char citation[] =
  "@Article{slepc-slice-qep,\n"
  "   author = \"C. Campos and J. E. Roman\",\n"
  "   title = \"Inertia-based spectrum slicing for symmetric quadratic eigenvalue problems\",\n"
  "   journal = \"Numer. Linear Algebra Appl.\",\n"
  "   volume = \"IP\",\n"
  "   number = \"x\",\n"
  "   pages = \"xx--xx\",\n"
  "   year = \"2020,\"\n"
  "   doi = \"https://doi.org/10.1002/nla.2293\"\n"
  "}\n";

#define SLICE_PTOL PETSC_SQRT_MACHINE_EPSILON

static PetscErrorCode PEPQSliceResetSR(PEP pep)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr=ctx->sr;
  PEP_shift      s;
  PetscInt       i;

  PetscFunctionBegin;
  if (sr) {
    /* Reviewing list of shifts to free memory */
    s = sr->s0;
    if (s) {
      while (s->neighb[1]) {
        s = s->neighb[1];
        ierr = PetscFree(s->neighb[0]);CHKERRQ(ierr);
      }
      ierr = PetscFree(s);CHKERRQ(ierr);
    }
    ierr = PetscFree(sr->S);CHKERRQ(ierr);
    for (i=0;i<pep->nconv;i++) {ierr = PetscFree(sr->qinfo[i].q);CHKERRQ(ierr);}
    ierr = PetscFree(sr->qinfo);CHKERRQ(ierr);
    for (i=0;i<3;i++) {ierr = VecDestroy(&sr->v[i]);CHKERRQ(ierr);}
    ierr = EPSDestroy(&sr->eps);CHKERRQ(ierr);
    ierr = PetscFree(sr);CHKERRQ(ierr);
  }
  ctx->sr = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode PEPReset_STOAR_QSlice(PEP pep)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  ierr = PEPQSliceResetSR(pep);CHKERRQ(ierr);
  ierr = PetscFree(ctx->inertias);CHKERRQ(ierr);
  ierr = PetscFree(ctx->shifts);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PEPQSliceAllocateSolution - Allocate memory storage for common variables such
  as eigenvalues and eigenvectors.
*/
static PetscErrorCode PEPQSliceAllocateSolution(PEP pep)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PetscInt       k;
  PetscLogDouble cnt;
  BVType         type;
  Vec            t;
  PEP_SR         sr = ctx->sr;

  PetscFunctionBegin;
  /* allocate space for eigenvalues and friends */
  k = PetscMax(1,sr->numEigs);
  ierr = PetscFree4(sr->eigr,sr->eigi,sr->errest,sr->perm);CHKERRQ(ierr);
  ierr = PetscCalloc4(k,&sr->eigr,k,&sr->eigi,k,&sr->errest,k,&sr->perm);CHKERRQ(ierr);
  ierr = PetscFree(sr->qinfo);CHKERRQ(ierr);
  ierr = PetscCalloc1(k,&sr->qinfo);CHKERRQ(ierr);
  cnt = 2*k*sizeof(PetscScalar) + 2*k*sizeof(PetscReal) + k*sizeof(PetscInt);
  ierr = PetscLogObjectMemory((PetscObject)pep,cnt);CHKERRQ(ierr);

  /* allocate sr->V and transfer options from pep->V */
  ierr = BVDestroy(&sr->V);CHKERRQ(ierr);
  ierr = BVCreate(PetscObjectComm((PetscObject)pep),&sr->V);CHKERRQ(ierr);
  ierr = PetscLogObjectParent((PetscObject)pep,(PetscObject)sr->V);CHKERRQ(ierr);
  if (!pep->V) { ierr = PEPGetBV(pep,&pep->V);CHKERRQ(ierr); }
  if (!((PetscObject)(pep->V))->type_name) {
    ierr = BVSetType(sr->V,BVSVEC);CHKERRQ(ierr);
  } else {
    ierr = BVGetType(pep->V,&type);CHKERRQ(ierr);
    ierr = BVSetType(sr->V,type);CHKERRQ(ierr);
  }
  ierr = STMatCreateVecsEmpty(pep->st,&t,NULL);CHKERRQ(ierr);
  ierr = BVSetSizesFromVec(sr->V,t,k+1);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  sr->ld = k;
  ierr = PetscFree(sr->S);CHKERRQ(ierr);
  ierr = PetscMalloc1((k+1)*sr->ld*(pep->nmat-1),&sr->S);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Convergence test to compute positive Ritz values */
static PetscErrorCode ConvergedPositive(EPS eps,PetscScalar eigr,PetscScalar eigi,PetscReal res,PetscReal *errest,void *ctx)
{
  PetscFunctionBegin;
  *errest = (PetscRealPart(eigr)>0.0)?0.0:res;
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPQSliceMatGetInertia(PEP pep,PetscReal shift,PetscInt *inertia,PetscInt *zeros)
{
  KSP            ksp,kspr;
  PC             pc;
  Mat            F;
  PetscBool      flg;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!pep->solvematcoeffs) {
    ierr = PetscMalloc1(pep->nmat,&pep->solvematcoeffs);CHKERRQ(ierr);
  }
  if (shift==PETSC_MAX_REAL) { /* Inertia of matrix A[2] */
    pep->solvematcoeffs[0] = 0.0; pep->solvematcoeffs[1] = 0.0; pep->solvematcoeffs[2] = 1.0;
  } else {
    ierr = PEPEvaluateBasis(pep,shift,0,pep->solvematcoeffs,NULL);CHKERRQ(ierr);
  }
  ierr = STMatSetUp(pep->st,pep->sfactor,pep->solvematcoeffs);CHKERRQ(ierr);
  ierr = STGetKSP(pep->st,&ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)pc,PCREDUNDANT,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PCRedundantGetKSP(pc,&kspr);CHKERRQ(ierr);
    ierr = KSPGetPC(kspr,&pc);CHKERRQ(ierr);
  }
  ierr = PCFactorGetMatrix(pc,&F);CHKERRQ(ierr);
  ierr = MatGetInertia(F,inertia,zeros,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPQSliceGetInertia(PEP pep,PetscReal shift,PetscInt *inertia,PetscInt *zeros,PetscInt correction)
{
  PetscErrorCode ierr;
  KSP            ksp;
  Mat            P;
  PetscReal      nzshift=0.0;
  PetscScalar    dot;
  PetscRandom    rand;
  PetscInt       nconv;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr=ctx->sr;

  PetscFunctionBegin;
  if (shift >= PETSC_MAX_REAL) { /* Right-open interval */
    *inertia = 0;
  } else if (shift <= PETSC_MIN_REAL) {
    *inertia = 0;
    if (zeros) *zeros = 0;
  } else {
    /* If the shift is zero, perturb it to a very small positive value.
       The goal is that the nonzero pattern is the same in all cases and reuse
       the symbolic factorizations */
    nzshift = (shift==0.0)? 10.0/PETSC_MAX_REAL: shift;
    ierr = PEPQSliceMatGetInertia(pep,nzshift,inertia,zeros);CHKERRQ(ierr);
    ierr = STSetShift(pep->st,nzshift);CHKERRQ(ierr);
  }
  if (!correction) {
    if (shift >= PETSC_MAX_REAL) *inertia = 2*pep->n;
    else if (shift>PETSC_MIN_REAL) {
      ierr = STGetKSP(pep->st,&ksp);CHKERRQ(ierr);
      ierr = KSPGetOperators(ksp,&P,NULL);CHKERRQ(ierr);
      if (*inertia!=pep->n && !sr->v[0]) {
        ierr = MatCreateVecs(P,&sr->v[0],NULL);CHKERRQ(ierr);
        ierr = VecDuplicate(sr->v[0],&sr->v[1]);CHKERRQ(ierr);
        ierr = VecDuplicate(sr->v[0],&sr->v[2]);CHKERRQ(ierr);
        ierr = BVGetRandomContext(pep->V,&rand);CHKERRQ(ierr);
        ierr = VecSetRandom(sr->v[0],rand);CHKERRQ(ierr);
      }
      if (*inertia<pep->n && *inertia>0) {
        if (!sr->eps) {
          ierr = EPSCreate(PetscObjectComm((PetscObject)pep),&sr->eps);CHKERRQ(ierr);
          ierr = EPSSetProblemType(sr->eps,EPS_HEP);CHKERRQ(ierr);
          ierr = EPSSetWhichEigenpairs(sr->eps,EPS_LARGEST_REAL);CHKERRQ(ierr);
        }
        ierr = EPSSetConvergenceTestFunction(sr->eps,ConvergedPositive,NULL,NULL);CHKERRQ(ierr);
        ierr = EPSSetOperators(sr->eps,P,NULL);CHKERRQ(ierr);
        ierr = EPSSolve(sr->eps);CHKERRQ(ierr);
        ierr = EPSGetConverged(sr->eps,&nconv);CHKERRQ(ierr);
        if (!nconv) SETERRQ1(((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Inertia computation fails in %g",nzshift);
        ierr = EPSGetEigenpair(sr->eps,0,NULL,NULL,sr->v[0],sr->v[1]);CHKERRQ(ierr);
      }
      if (*inertia!=pep->n) {
        ierr = MatMult(pep->A[1],sr->v[0],sr->v[1]);CHKERRQ(ierr);
        ierr = MatMult(pep->A[2],sr->v[0],sr->v[2]);CHKERRQ(ierr);
        ierr = VecAXPY(sr->v[1],2*nzshift,sr->v[2]);CHKERRQ(ierr);
        ierr = VecDot(sr->v[1],sr->v[0],&dot);CHKERRQ(ierr);
        if (PetscRealPart(dot)>0.0) *inertia = 2*pep->n-*inertia;
      }
    }
  } else if (correction<0) *inertia = 2*pep->n-*inertia;
  PetscFunctionReturn(0);
}

/*
   Check eigenvalue type - used only in non-hyperbolic problems.
   All computed eigenvalues must have the same definite type (positive or negative).
   If ini=TRUE the type is available in omega, otherwise we compute an eigenvalue
   closest to shift and determine its type.
 */
static PetscErrorCode PEPQSliceCheckEigenvalueType(PEP pep,PetscReal shift,PetscReal omega,PetscBool ini)
{
  PetscErrorCode ierr;
  PEP            pep2;
  ST             st;
  PetscInt       nconv;
  PetscScalar    lambda,dot;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr=ctx->sr;

  PetscFunctionBegin;
  if (!ini) {
    if (-(omega/(shift*ctx->alpha+ctx->beta))*sr->type<0) SETERRQ1(((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Different positive/negative type detected in eigenvalue %g",(double)shift);
  } else {
    ierr = PEPCreate(PetscObjectComm((PetscObject)pep),&pep2);CHKERRQ(ierr);
    ierr = PEPSetOptionsPrefix(pep2,((PetscObject)pep)->prefix);CHKERRQ(ierr);
    ierr = PEPAppendOptionsPrefix(pep2,"pep_eigenvalue_type_");CHKERRQ(ierr);
    ierr = PEPSetTolerances(pep2,PETSC_DEFAULT,pep->max_it/4);CHKERRQ(ierr);
    ierr = PEPSetType(pep2,PEPTOAR);CHKERRQ(ierr);
    ierr = PEPSetOperators(pep2,pep->nmat,pep->A);CHKERRQ(ierr);
    ierr = PEPSetWhichEigenpairs(pep2,PEP_TARGET_MAGNITUDE);CHKERRQ(ierr);
    ierr = PEPGetRG(pep2,&pep2->rg);CHKERRQ(ierr);
    ierr = RGSetType(pep2->rg,RGINTERVAL);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    ierr = RGIntervalSetEndpoints(pep2->rg,pep->inta,pep->intb,-PETSC_SQRT_MACHINE_EPSILON,PETSC_SQRT_MACHINE_EPSILON);CHKERRQ(ierr);
#else
    ierr = RGIntervalSetEndpoints(pep2->rg,pep->inta,pep->intb,0.0,0.0);CHKERRQ(ierr);
#endif
    pep2->target = shift;
    st = pep2->st;
    pep2->st = pep->st;
    ierr = PEPSolve(pep2);CHKERRQ(ierr);
    ierr = PEPGetConverged(pep2,&nconv);CHKERRQ(ierr);
    if (nconv) {
      ierr = PEPGetEigenpair(pep2,0,&lambda,NULL,pep2->work[0],NULL);CHKERRQ(ierr);
      ierr = MatMult(pep->A[1],pep2->work[0],pep2->work[1]);CHKERRQ(ierr);
      ierr = MatMult(pep->A[2],pep2->work[0],pep2->work[2]);CHKERRQ(ierr);
      ierr = VecAXPY(pep2->work[1],2.0*lambda*pep->sfactor,pep2->work[2]);CHKERRQ(ierr);
      ierr = VecDot(pep2->work[1],pep2->work[0],&dot);CHKERRQ(ierr);
      ierr = PetscInfo2(pep,"lambda=%g, %s type\n",(double)PetscRealPart(lambda),(PetscRealPart(dot)>0.0)?"positive":"negative");CHKERRQ(ierr);
      if (!sr->type) sr->type = (PetscRealPart(dot)>0.0)?1:-1;
      else {
        if (sr->type*PetscRealPart(dot)<0.0) SETERRQ1(((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Different positive/negative type detected in eigenvalue %g",(double)PetscRealPart(lambda));
      }
    }
    pep2->st = st;
    ierr = PEPDestroy(&pep2);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode PEPQSliceDiscriminant(PEP pep,Vec u,Vec w,PetscReal *d,PetscReal *smas,PetscReal *smenos)
{
  PetscReal      ap,bp,cp,dis;
  PetscScalar    ts;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatMult(pep->A[0],u,w);CHKERRQ(ierr);
  ierr = VecDot(w,u,&ts);CHKERRQ(ierr);
  cp = PetscRealPart(ts);
  ierr = MatMult(pep->A[1],u,w);CHKERRQ(ierr);
  ierr = VecDot(w,u,&ts);CHKERRQ(ierr);
  bp = PetscRealPart(ts);
  ierr = MatMult(pep->A[2],u,w);CHKERRQ(ierr);
  ierr = VecDot(w,u,&ts);CHKERRQ(ierr);
  ap = PetscRealPart(ts);
  dis = bp*bp-4*ap*cp;
  if (dis>=0.0 && smas) {
    if (ap>0) *smas = (-bp+PetscSqrtReal(dis))/(2*ap);
    else if (ap<0) *smas = (-bp-PetscSqrtReal(dis))/(2*ap);
    else {
      if (bp >0) *smas = -cp/bp;
      else *smas = PETSC_MAX_REAL;
    }
  }
  if (dis>=0.0 && smenos) {
    if (ap>0) *smenos = (-bp-PetscSqrtReal(dis))/(2*ap);
    else if (ap<0) *smenos = (-bp+PetscSqrtReal(dis))/(2*ap);
    else {
      if (bp<0) *smenos = -cp/bp;
      else *smenos = PETSC_MAX_REAL;
    }
  }
  if (d) *d = dis;
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode PEPQSliceEvaluateQEP(PEP pep,PetscScalar x,Mat M,MatStructure str)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatCopy(pep->A[0],M,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatAXPY(M,x,pep->A[1],str);CHKERRQ(ierr);
  ierr = MatAXPY(M,x*x,pep->A[2],str);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
   PEPCheckDefiniteQEP - Determines if a symmetric/Hermitian quadratic eigenvalue problem
   is definite or not.

   Logically Collective on pep

   Input Parameter:
.  pep  - eigensolver context

   Output Parameters:
+  xi - first computed parameter
.  mu - second computed parameter
.  definite - flag indicating that the problem is definite
-  hyperbolic - flag indicating that the problem is hyperbolic

   Notes:
   This function is intended for quadratic eigenvalue problems, Q(lambda)=A*lambda^2+B*lambda+C,
   with symmetric (or Hermitian) coefficient matrices A,B,C.

   On output, the flag 'definite' may have the values -1 (meaning that the QEP is not
   definite), 1 (if the problem is definite), or 0 if the algorithm was not able to
   determine whether the problem is definite or not.

   If definite=1, the output flag 'hyperbolic' informs in a similar way about whether the
   problem is hyperbolic or not.

   If definite=1, the computed values xi and mu satisfy Q(xi)<0 and Q(mu)>0, as
   obtained via the method proposed in [Niendorf and Voss, LAA 2010]. Furthermore, if
   hyperbolic=1 then only xi is computed.

   Level: advanced
@*/
PetscErrorCode PEPCheckDefiniteQEP(PEP pep,PetscReal *xi,PetscReal *mu,PetscInt *definite,PetscInt *hyperbolic)
{
  PetscErrorCode ierr;
  PetscRandom    rand;
  Vec            u,w;
  PetscReal      d=0.0,s=0.0,sp,mut=0.0,omg=0.0,omgp;
  PetscInt       k,its=10,hyp=0,check=0,nconv,inertia,n;
  Mat            M=NULL;
  MatStructure   str;
  EPS            eps;
  PetscBool      transform,ptypehyp;

  PetscFunctionBegin;
  if (pep->problem_type!=PEP_HERMITIAN && pep->problem_type!=PEP_HYPERBOLIC) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Only available for Hermitian (or hyperbolic) problems");
  ptypehyp = (pep->problem_type==PEP_HYPERBOLIC)? PETSC_TRUE: PETSC_FALSE;
  if (!pep->st) { ierr = PEPGetST(pep,&pep->st);CHKERRQ(ierr); }
  ierr = PEPSetDefaultST(pep);CHKERRQ(ierr);
  ierr = STSetMatrices(pep->st,pep->nmat,pep->A);CHKERRQ(ierr);
  ierr = MatGetSize(pep->A[0],&n,NULL);CHKERRQ(ierr);
  ierr = STGetTransform(pep->st,&transform);CHKERRQ(ierr);
  ierr = STSetTransform(pep->st,PETSC_FALSE);CHKERRQ(ierr);
  ierr = STSetUp(pep->st);CHKERRQ(ierr);
  ierr = MatCreateVecs(pep->A[0],&u,&w);CHKERRQ(ierr);
  ierr = PEPGetBV(pep,&pep->V);CHKERRQ(ierr);
  ierr = BVGetRandomContext(pep->V,&rand);CHKERRQ(ierr);
  ierr = VecSetRandom(u,rand);CHKERRQ(ierr);
  ierr = VecNormalize(u,NULL);CHKERRQ(ierr);
  ierr = PEPQSliceDiscriminant(pep,u,w,&d,&s,NULL);CHKERRQ(ierr);
  if (d<0.0) check = -1;
  if (!check) {
    ierr = EPSCreate(PetscObjectComm((PetscObject)pep),&eps);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL);CHKERRQ(ierr);
    ierr = EPSSetTolerances(eps,PetscSqrtReal(PETSC_SQRT_MACHINE_EPSILON),PETSC_DECIDE);CHKERRQ(ierr);
    ierr = MatDuplicate(pep->A[0],MAT_DO_NOT_COPY_VALUES,&M);CHKERRQ(ierr);
    ierr = STGetMatStructure(pep->st,&str);CHKERRQ(ierr);
  }
  for (k=0;k<its&&!check;k++) {
    ierr = PEPQSliceEvaluateQEP(pep,s,M,str);CHKERRQ(ierr);
    ierr = EPSSetOperators(eps,M,NULL);CHKERRQ(ierr);
    ierr = EPSSolve(eps);CHKERRQ(ierr);
    ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    if (!nconv) break;
    ierr = EPSGetEigenpair(eps,0,NULL,NULL,u,w);CHKERRQ(ierr);
    sp = s;
    ierr = PEPQSliceDiscriminant(pep,u,w,&d,&s,&omg);CHKERRQ(ierr);
    if (d<0.0) {check = -1; break;}
    if (PetscAbsReal((s-sp)/s)<100*PETSC_MACHINE_EPSILON) break;
    if (s>sp) {hyp = -1;}
    mut = 2*s-sp;
    ierr =  PEPQSliceMatGetInertia(pep,mut,&inertia,NULL);CHKERRQ(ierr);
    if (inertia == n) {check = 1; break;}
  }
  for (;k<its&&!check;k++) {
    mut = (s-omg)/2;
    ierr =  PEPQSliceMatGetInertia(pep,mut,&inertia,NULL);CHKERRQ(ierr);
    if (inertia == n) {check = 1; break;}
    if (PetscAbsReal((s-omg)/omg)<100*PETSC_MACHINE_EPSILON) break;
    ierr = PEPQSliceEvaluateQEP(pep,omg,M,str);CHKERRQ(ierr);
    ierr = EPSSetOperators(eps,M,NULL);CHKERRQ(ierr);
    ierr = EPSSolve(eps);CHKERRQ(ierr);
    ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    if (!nconv) break;
    ierr = EPSGetEigenpair(eps,0,NULL,NULL,u,w);CHKERRQ(ierr);
    omgp = omg;
    ierr = PEPQSliceDiscriminant(pep,u,w,&d,NULL,&omg);CHKERRQ(ierr);
    if (d<0.0) {check = -1; break;}
    if (omg<omgp) hyp = -1;
  }
  if (check==1) *xi = mut;
  if (hyp==-1 && ptypehyp) SETERRQ(PetscObjectComm((PetscObject)pep),1,"Problem does not satisfy hyperbolic test; consider removing the hyperbolicity flag");
  if (check==1 && hyp==0) {
    ierr =  PEPQSliceMatGetInertia(pep,PETSC_MAX_REAL,&inertia,NULL);CHKERRQ(ierr);
    if (inertia == 0) hyp = 1;
    else hyp = -1;
  }
  if (check==1 && hyp!=1) {
    check = 0;
    ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
    for (;k<its&&!check;k++) {
      ierr = PEPQSliceEvaluateQEP(pep,s,M,str);CHKERRQ(ierr);
      ierr = EPSSetOperators(eps,M,NULL);CHKERRQ(ierr);
      ierr = EPSSolve(eps);CHKERRQ(ierr);
      ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
      if (!nconv) break;
      ierr = EPSGetEigenpair(eps,0,NULL,NULL,u,w);CHKERRQ(ierr);
      sp = s;
      ierr = PEPQSliceDiscriminant(pep,u,w,&d,&s,&omg);CHKERRQ(ierr);
      if (d<0.0) {check = -1; break;}
      if (PetscAbsReal((s-sp)/s)<100*PETSC_MACHINE_EPSILON) break;
      mut = 2*s-sp;
      ierr =  PEPQSliceMatGetInertia(pep,mut,&inertia,NULL);CHKERRQ(ierr);
      if (inertia == 0) {check = 1; break;}
    }
    for (;k<its&&!check;k++) {
      mut = (s-omg)/2;
      ierr =  PEPQSliceMatGetInertia(pep,mut,&inertia,NULL);CHKERRQ(ierr);
      if (inertia == 0) {check = 1; break;}
      if (PetscAbsReal((s-omg)/omg)<100*PETSC_MACHINE_EPSILON) break;
      ierr = PEPQSliceEvaluateQEP(pep,omg,M,str);CHKERRQ(ierr);
      ierr = EPSSetOperators(eps,M,NULL);CHKERRQ(ierr);
      ierr = EPSSolve(eps);CHKERRQ(ierr);
      ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
      if (!nconv) break;
      ierr = EPSGetEigenpair(eps,0,NULL,NULL,u,w);CHKERRQ(ierr);
      ierr = PEPQSliceDiscriminant(pep,u,w,&d,NULL,&omg);CHKERRQ(ierr);
      if (d<0.0) {check = -1; break;}
    }
  }
  if (check==1) *mu = mut;
  *definite = check;
  *hyperbolic = hyp;
  if (M) { ierr = MatDestroy(&M);CHKERRQ(ierr); }
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&w);CHKERRQ(ierr);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = STSetTransform(pep->st,transform);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Dummy backtransform operation
 */
static PetscErrorCode PEPBackTransform_Skip(PEP pep)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSetUp_STOAR_QSlice(PEP pep)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr;
  PetscInt       ld,i,zeros=0;
  SlepcSC        sc;
  PetscBool      issinv;
  PetscReal      r;

  PetscFunctionBegin;
  if (pep->intb >= PETSC_MAX_REAL && pep->inta <= PETSC_MIN_REAL) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"The defined computational interval should have at least one of their sides bounded");
  ierr = PetscObjectTypeCompareAny((PetscObject)pep->st,&issinv,STSINVERT,STCAYLEY,"");CHKERRQ(ierr);
  if (!issinv) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"Shift-and-invert or Cayley ST is needed for spectrum slicing");
  if (pep->tol==PETSC_DEFAULT) pep->tol = SLEPC_DEFAULT_TOL*1e-2;  /* use tighter tolerance */
  if (ctx->nev==1) ctx->nev = PetscMin(20,pep->n);  /* nev not set, use default value */
  if (pep->n>10 && ctx->nev<10) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONG,"nev cannot be less than 10 in spectrum slicing runs");
  pep->ops->backtransform = PEPBackTransform_Skip;
  if (pep->max_it==PETSC_DEFAULT) pep->max_it = 100;

  /* create spectrum slicing context and initialize it */
  ierr = PEPQSliceResetSR(pep);CHKERRQ(ierr);
  ierr = PetscNewLog(pep,&sr);CHKERRQ(ierr);
  ctx->sr   = sr;
  sr->itsKs = 0;
  sr->nleap = 0;
  sr->sPres = NULL;

  if (pep->solvematcoeffs) { ierr = PetscFree(pep->solvematcoeffs);CHKERRQ(ierr); }
  ierr = PetscMalloc1(pep->nmat,&pep->solvematcoeffs);CHKERRQ(ierr);
  if (!pep->st) { ierr = PEPGetST(pep,&pep->st);CHKERRQ(ierr); }
  ierr = STSetTransform(pep->st,PETSC_FALSE);CHKERRQ(ierr);
  ierr = STSetUp(pep->st);CHKERRQ(ierr);

  ctx->hyperbolic = (pep->problem_type==PEP_HYPERBOLIC)? PETSC_TRUE: PETSC_FALSE;

  /* check presence of ends and finding direction */
  if (pep->inta > PETSC_MIN_REAL || pep->intb >= PETSC_MAX_REAL) {
    sr->int0 = pep->inta;
    sr->int1 = pep->intb;
    sr->dir = 1;
    if (pep->intb >= PETSC_MAX_REAL) { /* Right-open interval */
      sr->hasEnd = PETSC_FALSE;
    } else sr->hasEnd = PETSC_TRUE;
  } else {
    sr->int0 = pep->intb;
    sr->int1 = pep->inta;
    sr->dir = -1;
    sr->hasEnd = PetscNot(pep->inta <= PETSC_MIN_REAL);
  }

  /* compute inertia0 */
  ierr = PEPQSliceGetInertia(pep,sr->int0,&sr->inertia0,ctx->detect?&zeros:NULL,ctx->hyperbolic?0:1);CHKERRQ(ierr);
  if (zeros && (sr->int0==pep->inta || sr->int0==pep->intb)) SETERRQ(((PetscObject)pep)->comm,PETSC_ERR_USER,"Found singular matrix for the transformed problem in the interval endpoint");
  if (!ctx->hyperbolic && ctx->checket) {
    ierr = PEPQSliceCheckEigenvalueType(pep,sr->int0,0.0,PETSC_TRUE);CHKERRQ(ierr);
  }

  /* compute inertia1 */
  ierr = PEPQSliceGetInertia(pep,sr->int1,&sr->inertia1,ctx->detect?&zeros:NULL,ctx->hyperbolic?0:1);CHKERRQ(ierr);
  if (zeros) SETERRQ(((PetscObject)pep)->comm,PETSC_ERR_USER,"Found singular matrix for the transformed problem in an interval endpoint defined by user");
  if (!ctx->hyperbolic && ctx->checket && sr->hasEnd) {
    ierr = PEPQSliceCheckEigenvalueType(pep,sr->int1,0.0,PETSC_TRUE);CHKERRQ(ierr);
    if (!sr->type && (sr->inertia1-sr->inertia0)) SETERRQ(((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"No information of eigenvalue type in Interval");
    if (sr->type && !(sr->inertia1-sr->inertia0)) SETERRQ(((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Different positive/negative type detected");
    if (sr->dir*(sr->inertia1-sr->inertia0)<0) {
      sr->intcorr = -1;
      sr->inertia0 = 2*pep->n-sr->inertia0;
      sr->inertia1 = 2*pep->n-sr->inertia1;
    } else sr->intcorr = 1;
  } else {
    if (sr->inertia0<=pep->n && sr->inertia1<=pep->n) sr->intcorr = 1;
    else if (sr->inertia0>=pep->n && sr->inertia1>=pep->n) sr->intcorr = -1;
  }

  if (sr->hasEnd) {
    sr->dir = -sr->dir; r = sr->int0; sr->int0 = sr->int1; sr->int1 = r;
    i = sr->inertia0; sr->inertia0 = sr->inertia1; sr->inertia1 = i;
  }

  /* number of eigenvalues in interval */
  sr->numEigs = (sr->dir)*(sr->inertia1 - sr->inertia0);
  ierr = PetscInfo3(pep,"QSlice setup: allocating for %D eigenvalues in [%g,%g]\n",sr->numEigs,(double)pep->inta,(double)pep->intb);CHKERRQ(ierr);
  if (sr->numEigs) {
    ierr = PEPQSliceAllocateSolution(pep);CHKERRQ(ierr);
    ierr = PEPSetDimensions_Default(pep,ctx->nev,&ctx->ncv,&ctx->mpd);CHKERRQ(ierr);
    pep->nev = ctx->nev; pep->ncv = ctx->ncv; pep->mpd = ctx->mpd;
    ld   = ctx->ncv+2;
    ierr = DSSetType(pep->ds,DSGHIEP);CHKERRQ(ierr);
    ierr = DSSetCompact(pep->ds,PETSC_TRUE);CHKERRQ(ierr);
    ierr = DSAllocate(pep->ds,ld);CHKERRQ(ierr);
    ierr = DSGetSlepcSC(pep->ds,&sc);CHKERRQ(ierr);
    sc->rg            = NULL;
    sc->comparison    = SlepcCompareLargestMagnitude;
    sc->comparisonctx = NULL;
    sc->map           = NULL;
    sc->mapobj        = NULL;
  }
  PetscFunctionReturn(0);
}

/*
   Fills the fields of a shift structure
*/
static PetscErrorCode PEPCreateShift(PEP pep,PetscReal val,PEP_shift neighb0,PEP_shift neighb1)
{
  PetscErrorCode ierr;
  PEP_shift      s,*pending2;
  PetscInt       i;
  PEP_SR         sr;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  sr = ctx->sr;
  ierr = PetscNewLog(pep,&s);CHKERRQ(ierr);
  s->value = val;
  s->neighb[0] = neighb0;
  if (neighb0) neighb0->neighb[1] = s;
  s->neighb[1] = neighb1;
  if (neighb1) neighb1->neighb[0] = s;
  s->comp[0] = PETSC_FALSE;
  s->comp[1] = PETSC_FALSE;
  s->index = -1;
  s->neigs = 0;
  s->nconv[0] = s->nconv[1] = 0;
  s->nsch[0] = s->nsch[1]=0;
  /* Inserts in the stack of pending shifts */
  /* If needed, the array is resized */
  if (sr->nPend >= sr->maxPend) {
    sr->maxPend *= 2;
    ierr = PetscMalloc1(sr->maxPend,&pending2);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)pep,sizeof(PEP_shift));CHKERRQ(ierr);
    for (i=0;i<sr->nPend;i++) pending2[i] = sr->pending[i];
    ierr = PetscFree(sr->pending);CHKERRQ(ierr);
    sr->pending = pending2;
  }
  sr->pending[sr->nPend++]=s;
  PetscFunctionReturn(0);
}

/* Provides next shift to be computed */
static PetscErrorCode PEPExtractShift(PEP pep)
{
  PetscErrorCode ierr;
  PetscInt       iner,zeros=0;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr;
  PetscReal      newShift,aux;
  PEP_shift      sPres;

  PetscFunctionBegin;
  sr = ctx->sr;
  if (sr->nPend > 0) {
    if (sr->dirch) {
      aux = sr->int1; sr->int1 = sr->int0; sr->int0 = aux;
      iner = sr->inertia1; sr->inertia1 = sr->inertia0; sr->inertia0 = iner;
      sr->dir *= -1;
      ierr = PetscFree(sr->s0->neighb[1]);CHKERRQ(ierr);
      ierr = PetscFree(sr->s0);CHKERRQ(ierr);
      sr->nPend--;
      ierr = PEPCreateShift(pep,sr->int0,NULL,NULL);CHKERRQ(ierr);
      sr->sPrev = NULL;
      sr->sPres = sr->pending[--sr->nPend];
      pep->target = sr->sPres->value;
      sr->s0 = sr->sPres;
      pep->reason = PEP_CONVERGED_ITERATING;
    } else {
      sr->sPrev = sr->sPres;
      sr->sPres = sr->pending[--sr->nPend];
    }
    sPres = sr->sPres;
    ierr = PEPQSliceGetInertia(pep,sPres->value,&iner,ctx->detect?&zeros:NULL,sr->intcorr);CHKERRQ(ierr);
    if (zeros) {
      newShift = sPres->value*(1.0+SLICE_PTOL);
      if (sr->dir*(sPres->neighb[0] && newShift-sPres->neighb[0]->value) < 0) newShift = (sPres->value+sPres->neighb[0]->value)/2;
      else if (sPres->neighb[1] && sr->dir*(sPres->neighb[1]->value-newShift) < 0) newShift = (sPres->value+sPres->neighb[1]->value)/2;
      ierr = PEPQSliceGetInertia(pep,newShift,&iner,&zeros,sr->intcorr);CHKERRQ(ierr);
      if (zeros) SETERRQ1(((PetscObject)pep)->comm,PETSC_ERR_CONV_FAILED,"Inertia computation fails in %g",newShift);
      sPres->value = newShift;
    }
    sr->sPres->inertia = iner;
    pep->target = sr->sPres->value;
    pep->reason = PEP_CONVERGED_ITERATING;
    pep->its = 0;
  } else sr->sPres = NULL;
  PetscFunctionReturn(0);
}

/*
  Obtains value of subsequent shift
*/
static PetscErrorCode PEPGetNewShiftValue(PEP pep,PetscInt side,PetscReal *newS)
{
  PetscReal lambda,d_prev;
  PetscInt  i,idxP;
  PEP_SR    sr;
  PEP_shift sPres,s;
  PEP_STOAR *ctx=(PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  sr = ctx->sr;
  sPres = sr->sPres;
  if (sPres->neighb[side]) {
  /* Completing a previous interval */
    if (!sPres->neighb[side]->neighb[side] && sPres->neighb[side]->nconv[side]==0) { /* One of the ends might be too far from eigenvalues */
      if (side) *newS = (sPres->value + PetscRealPart(sr->eigr[sr->perm[sr->indexEig-1]]))/2;
      else *newS = (sPres->value + PetscRealPart(sr->eigr[sr->perm[0]]))/2;
    } else *newS=(sPres->value + sPres->neighb[side]->value)/2;
  } else { /* (Only for side=1). Creating a new interval. */
    if (sPres->neigs==0) {/* No value has been accepted*/
      if (sPres->neighb[0]) {
        /* Multiplying by 10 the previous distance */
        *newS = sPres->value + 10*(sr->dir)*PetscAbsReal(sPres->value - sPres->neighb[0]->value);
        sr->nleap++;
        /* Stops when the interval is open and no values are found in the last 5 shifts (there might be infinite eigenvalues) */
        if (!sr->hasEnd && sr->nleap > 5) SETERRQ(PetscObjectComm((PetscObject)pep),1,"Unable to compute the wanted eigenvalues with open interval");
      } else { /* First shift */
        if (pep->nconv != 0) {
          /* Unaccepted values give information for next shift */
          idxP=0;/* Number of values left from shift */
          for (i=0;i<pep->nconv;i++) {
            lambda = PetscRealPart(pep->eigr[i]);
            if ((sr->dir)*(lambda - sPres->value) <0) idxP++;
            else break;
          }
          /* Avoiding subtraction of eigenvalues (might be the same).*/
          if (idxP>0) {
            d_prev = PetscAbsReal(sPres->value - PetscRealPart(pep->eigr[0]))/(idxP+0.3);
          } else {
            d_prev = PetscAbsReal(sPres->value - PetscRealPart(pep->eigr[pep->nconv-1]))/(pep->nconv+0.3);
          }
          *newS = sPres->value + ((sr->dir)*d_prev*pep->nev)/2;
          sr->dirch = PETSC_FALSE;
        } else { /* No values found, no information for next shift */
          if (!sr->dirch) {
            sr->dirch = PETSC_TRUE;
            *newS = sr->int1;
          } else SETERRQ(PetscObjectComm((PetscObject)pep),1,"First shift renders no information");
        }
      }
    } else { /* Accepted values found */
      sr->dirch = PETSC_FALSE;
      sr->nleap = 0;
      /* Average distance of values in previous subinterval */
      s = sPres->neighb[0];
      while (s && PetscAbs(s->inertia - sPres->inertia)==0) {
        s = s->neighb[0];/* Looking for previous shifts with eigenvalues within */
      }
      if (s) {
        d_prev = PetscAbsReal((sPres->value - s->value)/(sPres->inertia - s->inertia));
      } else { /* First shift. Average distance obtained with values in this shift */
        /* first shift might be too far from first wanted eigenvalue (no values found outside the interval)*/
        if ((sr->dir)*(PetscRealPart(sr->eigr[0])-sPres->value)>0 && PetscAbsReal((PetscRealPart(sr->eigr[sr->indexEig-1]) - PetscRealPart(sr->eigr[0]))/PetscRealPart(sr->eigr[0])) > PetscSqrtReal(pep->tol)) {
          d_prev =  PetscAbsReal((PetscRealPart(sr->eigr[sr->indexEig-1]) - PetscRealPart(sr->eigr[0])))/(sPres->neigs+0.3);
        } else {
          d_prev = PetscAbsReal(PetscRealPart(sr->eigr[sr->indexEig-1]) - sPres->value)/(sPres->neigs+0.3);
        }
      }
      /* Average distance is used for next shift by adding it to value on the right or to shift */
      if ((sr->dir)*(PetscRealPart(sr->eigr[sPres->index + sPres->neigs -1]) - sPres->value)>0) {
        *newS = PetscRealPart(sr->eigr[sPres->index + sPres->neigs -1])+ ((sr->dir)*d_prev*(pep->nev))/2;
      } else { /* Last accepted value is on the left of shift. Adding to shift */
        *newS = sPres->value + ((sr->dir)*d_prev*(pep->nev))/2;
      }
    }
    /* End of interval can not be surpassed */
    if ((sr->dir)*(sr->int1 - *newS) < 0) *newS = sr->int1;
  }/* of neighb[side]==null */
  PetscFunctionReturn(0);
}

/*
  Function for sorting an array of real values
*/
static PetscErrorCode sortRealEigenvalues(PetscScalar *r,PetscInt *perm,PetscInt nr,PetscBool prev,PetscInt dir)
{
  PetscReal re;
  PetscInt  i,j,tmp;

  PetscFunctionBegin;
  if (!prev) for (i=0;i<nr;i++) perm[i] = i;
  /* Insertion sort */
  for (i=1;i<nr;i++) {
    re = PetscRealPart(r[perm[i]]);
    j = i-1;
    while (j>=0 && dir*(re - PetscRealPart(r[perm[j]])) <= 0) {
      tmp = perm[j]; perm[j] = perm[j+1]; perm[j+1] = tmp; j--;
    }
  }
  PetscFunctionReturn(0);
}

/* Stores the pairs obtained since the last shift in the global arrays */
static PetscErrorCode PEPStoreEigenpairs(PEP pep)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PetscReal      lambda,err,*errest;
  PetscInt       i,*aux,count=0,ndef,ld,nconv=pep->nconv,d=pep->nmat-1,idx;
  PetscBool      iscayley,divide=PETSC_FALSE;
  PEP_SR         sr = ctx->sr;
  PEP_shift      sPres;
  Vec            w,vomega;
  Mat            MS;
  BV             tV;
  PetscScalar    *S,*eigr,*tS,*omega;

  PetscFunctionBegin;
  sPres = sr->sPres;
  sPres->index = sr->indexEig;

  if (nconv>sr->ndef0+sr->ndef1) {
    /* Back-transform */
    ierr = STBackTransform(pep->st,nconv,pep->eigr,pep->eigi);CHKERRQ(ierr);
    for (i=0;i<nconv;i++) {
#if defined(PETSC_USE_COMPLEX)
      if (PetscImaginaryPart(pep->eigr[i])) pep->eigr[i] = sr->int0-sr->dir;
#else
      if (pep->eigi[i]) pep->eigr[i] = sr->int0-sr->dir;
#endif
    }
    ierr = PetscObjectTypeCompare((PetscObject)pep->st,STCAYLEY,&iscayley);CHKERRQ(ierr);
    /* Sort eigenvalues */
    ierr = sortRealEigenvalues(pep->eigr,pep->perm,nconv,PETSC_FALSE,sr->dir);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,nconv,&vomega);CHKERRQ(ierr);
    ierr = BVGetSignature(ctx->V,vomega);CHKERRQ(ierr);
    ierr = VecGetArray(vomega,&omega);CHKERRQ(ierr);
    ierr = BVGetSizes(pep->V,NULL,NULL,&ld);CHKERRQ(ierr);
    ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
    /* Values stored in global array */
    ierr = PetscCalloc4(nconv,&eigr,nconv,&errest,nconv*nconv*d,&tS,nconv,&aux);CHKERRQ(ierr);
    ndef = sr->ndef0+sr->ndef1;
    for (i=0;i<nconv;i++) {
      lambda = PetscRealPart(pep->eigr[pep->perm[i]]);
      err = pep->errest[pep->perm[i]];
      if ((sr->dir)*(lambda - sPres->ext[0]) > 0 && (sr->dir)*(sPres->ext[1] - lambda) > 0) {/* Valid value */
        if (sr->indexEig+count-ndef>=sr->numEigs) SETERRQ(PetscObjectComm((PetscObject)pep),1,"Unexpected error in Spectrum Slicing");
        ierr = PEPQSliceCheckEigenvalueType(pep,lambda,PetscRealPart(omega[pep->perm[i]]),PETSC_FALSE);CHKERRQ(ierr);
        eigr[count] = lambda;
        errest[count] = err;
        if (((sr->dir)*(sPres->value - lambda) > 0) && ((sr->dir)*(lambda - sPres->ext[0]) > 0)) sPres->nconv[0]++;
        if (((sr->dir)*(lambda - sPres->value) > 0) && ((sr->dir)*(sPres->ext[1] - lambda) > 0)) sPres->nconv[1]++;
        ierr = PetscArraycpy(tS+count*(d*nconv),S+pep->perm[i]*(d*ld),nconv);CHKERRQ(ierr);
        ierr = PetscArraycpy(tS+count*(d*nconv)+nconv,S+pep->perm[i]*(d*ld)+ld,nconv);CHKERRQ(ierr);
        count++;
      }
    }
    ierr = VecRestoreArray(vomega,&omega);CHKERRQ(ierr);
    ierr = VecDestroy(&vomega);CHKERRQ(ierr);
    for (i=0;i<count;i++) {
      ierr = PetscArraycpy(S+i*(d*ld),tS+i*nconv*d,nconv);CHKERRQ(ierr);
      ierr = PetscArraycpy(S+i*(d*ld)+ld,tS+i*nconv*d+nconv,nconv);CHKERRQ(ierr);
    }
    ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
    ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(ctx->V,0,count);CHKERRQ(ierr);
    ierr = BVTensorCompress(ctx->V,count);CHKERRQ(ierr);
    if (sr->sPres->nconv[0] && sr->sPres->nconv[1]) {
      divide = PETSC_TRUE;
      ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
      ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
      ierr = PetscArrayzero(tS,nconv*nconv*d);CHKERRQ(ierr);
      for (i=0;i<count;i++) {
        ierr = PetscArraycpy(tS+i*nconv*d,S+i*(d*ld),count);CHKERRQ(ierr);
        ierr = PetscArraycpy(tS+i*nconv*d+nconv,S+i*(d*ld)+ld,count);CHKERRQ(ierr);
      }
      ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
      ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(pep->V,0,count);CHKERRQ(ierr);
      ierr = BVDuplicateResize(pep->V,count,&tV);CHKERRQ(ierr);
      ierr = BVCopy(pep->V,tV);CHKERRQ(ierr);
    }
    if (sr->sPres->nconv[0]) {
      if (divide) {
        ierr = BVSetActiveColumns(ctx->V,0,sr->sPres->nconv[0]);CHKERRQ(ierr);
        ierr = BVTensorCompress(ctx->V,sr->sPres->nconv[0]);CHKERRQ(ierr);
      }
      for (i=0;i<sr->ndef0;i++) aux[i] = sr->idxDef0[i];
      for (i=sr->ndef0;i<sr->sPres->nconv[0];i++) aux[i] = sr->indexEig+i-sr->ndef0;
      ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
      ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
      for (i=0;i<sr->sPres->nconv[0];i++) {
        sr->eigr[aux[i]] = eigr[i];
        sr->errest[aux[i]] = errest[i];
        ierr = BVGetColumn(pep->V,i,&w);CHKERRQ(ierr);
        ierr = BVInsertVec(sr->V,aux[i],w);CHKERRQ(ierr);
        ierr = BVRestoreColumn(pep->V,i,&w);CHKERRQ(ierr);
        idx = sr->ld*d*aux[i];
        ierr = PetscArrayzero(sr->S+idx,sr->ld*d);CHKERRQ(ierr);
        ierr = PetscArraycpy(sr->S+idx,S+i*(ld*d),sr->sPres->nconv[0]);CHKERRQ(ierr);
        ierr = PetscArraycpy(sr->S+idx+sr->ld,S+i*(ld*d)+ld,sr->sPres->nconv[0]);CHKERRQ(ierr);
        ierr = PetscFree(sr->qinfo[aux[i]].q);CHKERRQ(ierr);
        ierr = PetscMalloc1(sr->sPres->nconv[0],&sr->qinfo[aux[i]].q);CHKERRQ(ierr);
        ierr = PetscArraycpy(sr->qinfo[aux[i]].q,aux,sr->sPres->nconv[0]);CHKERRQ(ierr);
        sr->qinfo[aux[i]].nq = sr->sPres->nconv[0];
      }
      ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
      ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    }

    if (sr->sPres->nconv[1]) {
      if (divide) {
        ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
        ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
        for (i=0;i<sr->sPres->nconv[1];i++) {
          ierr = PetscArraycpy(S+i*(d*ld),tS+(sr->sPres->nconv[0]+i)*nconv*d,count);CHKERRQ(ierr);
          ierr = PetscArraycpy(S+i*(d*ld)+ld,tS+(sr->sPres->nconv[0]+i)*nconv*d+nconv,count);CHKERRQ(ierr);
        }
        ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
        ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
        ierr = BVSetActiveColumns(pep->V,0,count);CHKERRQ(ierr);
        ierr = BVCopy(tV,pep->V);CHKERRQ(ierr);
        ierr = BVSetActiveColumns(ctx->V,0,sr->sPres->nconv[1]);CHKERRQ(ierr);
        ierr = BVTensorCompress(ctx->V,sr->sPres->nconv[1]);CHKERRQ(ierr);
      }
      for (i=0;i<sr->ndef1;i++) aux[i] = sr->idxDef1[i];
      for (i=sr->ndef1;i<sr->sPres->nconv[1];i++) aux[i] = sr->indexEig+sr->sPres->nconv[0]-sr->ndef0+i-sr->ndef1;
      ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
      ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
      for (i=0;i<sr->sPres->nconv[1];i++) {
        sr->eigr[aux[i]] = eigr[sr->sPres->nconv[0]+i];
        sr->errest[aux[i]] = errest[sr->sPres->nconv[0]+i];
        ierr = BVGetColumn(pep->V,i,&w);CHKERRQ(ierr);
        ierr = BVInsertVec(sr->V,aux[i],w);CHKERRQ(ierr);
        ierr = BVRestoreColumn(pep->V,i,&w);CHKERRQ(ierr);
        idx = sr->ld*d*aux[i];
        ierr = PetscArrayzero(sr->S+idx,sr->ld*d);CHKERRQ(ierr);
        ierr = PetscArraycpy(sr->S+idx,S+i*(ld*d),sr->sPres->nconv[1]);CHKERRQ(ierr);
        ierr = PetscArraycpy(sr->S+idx+sr->ld,S+i*(ld*d)+ld,sr->sPres->nconv[1]);CHKERRQ(ierr);
        ierr = PetscFree(sr->qinfo[aux[i]].q);CHKERRQ(ierr);
        ierr = PetscMalloc1(sr->sPres->nconv[1],&sr->qinfo[aux[i]].q);CHKERRQ(ierr);
        ierr = PetscArraycpy(sr->qinfo[aux[i]].q,aux,sr->sPres->nconv[1]);CHKERRQ(ierr);
        sr->qinfo[aux[i]].nq = sr->sPres->nconv[1];
      }
      ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
      ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    }
    sPres->neigs = count-sr->ndef0-sr->ndef1;
    sr->indexEig += sPres->neigs;
    sPres->nconv[0]-= sr->ndef0;
    sPres->nconv[1]-= sr->ndef1;
    ierr = PetscFree4(eigr,errest,tS,aux);CHKERRQ(ierr);
  } else {
    sPres->neigs = 0;
    sPres->nconv[0]= 0;
    sPres->nconv[1]= 0;
  }
  /* Global ordering array updating */
  ierr = sortRealEigenvalues(sr->eigr,sr->perm,sr->indexEig,PETSC_FALSE,sr->dir);CHKERRQ(ierr);
  /* Check for completion */
  sPres->comp[0] = PetscNot(sPres->nconv[0] < sPres->nsch[0]);
  sPres->comp[1] = PetscNot(sPres->nconv[1] < sPres->nsch[1]);
  if (sPres->nconv[0] > sPres->nsch[0] || sPres->nconv[1] > sPres->nsch[1]) SETERRQ(PetscObjectComm((PetscObject)pep),1,"Mismatch between number of values found and information from inertia");
  if (divide) { ierr = BVDestroy(&tV);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPLookForDeflation(PEP pep)
{
  PetscReal val;
  PetscInt  i,count0=0,count1=0;
  PEP_shift sPres;
  PetscInt  ini,fin;
  PEP_SR    sr;
  PEP_STOAR *ctx=(PEP_STOAR*)pep->data;

  PetscFunctionBegin;
  sr = ctx->sr;
  sPres = sr->sPres;

  if (sPres->neighb[0]) ini = (sr->dir)*(sPres->neighb[0]->inertia - sr->inertia0);
  else ini = 0;
  fin = sr->indexEig;
  /* Selection of ends for searching new values */
  if (!sPres->neighb[0]) sPres->ext[0] = sr->int0;/* First shift */
  else sPres->ext[0] = sPres->neighb[0]->value;
  if (!sPres->neighb[1]) {
    if (sr->hasEnd) sPres->ext[1] = sr->int1;
    else sPres->ext[1] = (sr->dir > 0)?PETSC_MAX_REAL:PETSC_MIN_REAL;
  } else sPres->ext[1] = sPres->neighb[1]->value;
  /* Selection of values between right and left ends */
  for (i=ini;i<fin;i++) {
    val=PetscRealPart(sr->eigr[sr->perm[i]]);
    /* Values to the right of left shift */
    if ((sr->dir)*(val - sPres->ext[1]) < 0) {
      if ((sr->dir)*(val - sPres->value) < 0) count0++;
      else count1++;
    } else break;
  }
  /* The number of values on each side are found */
  if (sPres->neighb[0]) {
    sPres->nsch[0] = (sr->dir)*(sPres->inertia - sPres->neighb[0]->inertia)-count0;
    if (sPres->nsch[0]<0) SETERRQ(PetscObjectComm((PetscObject)pep),1,"Mismatch between number of values found and information from inertia");
  } else sPres->nsch[0] = 0;

  if (sPres->neighb[1]) {
    sPres->nsch[1] = (sr->dir)*(sPres->neighb[1]->inertia - sPres->inertia) - count1;
    if (sPres->nsch[1]<0) SETERRQ(PetscObjectComm((PetscObject)pep),1,"Mismatch between number of values found and information from inertia");
  } else sPres->nsch[1] = (sr->dir)*(sr->inertia1 - sPres->inertia);

  /* Completing vector of indexes for deflation */
  for (i=0;i<count0;i++) sr->idxDef0[i] = sr->perm[ini+i];
  sr->ndef0 = count0;
  for (i=0;i<count1;i++) sr->idxDef1[i] = sr->perm[ini+count0+i];
  sr->ndef1 = count1;
  PetscFunctionReturn(0);
}

/*
  Compute a run of Lanczos iterations
*/
static PetscErrorCode PEPSTOARrun_QSlice(PEP pep,PetscReal *a,PetscReal *b,PetscReal *omega,PetscInt k,PetscInt *M,PetscBool *breakdown,PetscBool *symmlost,Vec *t_)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscInt       i,j,m=*M,l,lock;
  PetscInt       lds,d,ld,offq,nqt,ldds;
  Vec            v=t_[0],t=t_[1],q=t_[2];
  PetscReal      norm,sym=0.0,fro=0.0,*f;
  PetscScalar    *y,*S,sigma;
  PetscBLASInt   j_,one=1;
  PetscBool      lindep;
  Mat            MS;

  PetscFunctionBegin;
  ierr = PetscMalloc1(*M,&y);CHKERRQ(ierr);
  ierr = BVGetSizes(pep->V,NULL,NULL,&ld);CHKERRQ(ierr);
  ierr = BVTensorGetDegree(ctx->V,&d);CHKERRQ(ierr);
  ierr = BVGetActiveColumns(pep->V,&lock,&nqt);CHKERRQ(ierr);
  lds = d*ld;
  offq = ld;
  ierr = DSGetLeadingDimension(pep->ds,&ldds);CHKERRQ(ierr);

  *breakdown = PETSC_FALSE; /* ----- */
  ierr = STGetShift(pep->st,&sigma);CHKERRQ(ierr);
  ierr = DSGetDimensions(pep->ds,NULL,NULL,&l,NULL,NULL);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(ctx->V,0,m);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(pep->V,0,nqt);CHKERRQ(ierr);
  for (j=k;j<m;j++) {
    /* apply operator */
    ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
    ierr = BVGetColumn(pep->V,nqt,&t);CHKERRQ(ierr);
    ierr = BVMultVec(pep->V,1.0,0.0,v,S+j*lds);CHKERRQ(ierr);
    ierr = MatMult(pep->A[1],v,q);CHKERRQ(ierr);
    ierr = MatMult(pep->A[2],v,t);CHKERRQ(ierr);
    ierr = VecAXPY(q,sigma*pep->sfactor,t);CHKERRQ(ierr);
    ierr = VecScale(q,pep->sfactor);CHKERRQ(ierr);
    ierr = BVMultVec(pep->V,1.0,0.0,v,S+offq+j*lds);CHKERRQ(ierr);
    ierr = MatMult(pep->A[2],v,t);CHKERRQ(ierr);
    ierr = VecAXPY(q,pep->sfactor*pep->sfactor,t);CHKERRQ(ierr);
    ierr = STMatSolve(pep->st,q,t);CHKERRQ(ierr);
    ierr = VecScale(t,-1.0);CHKERRQ(ierr);
    ierr = BVRestoreColumn(pep->V,nqt,&t);CHKERRQ(ierr);

    /* orthogonalize */
    ierr = BVOrthogonalizeColumn(pep->V,nqt,S+(j+1)*lds,&norm,&lindep);CHKERRQ(ierr);
    if (!lindep) {
      *(S+(j+1)*lds+nqt) = norm;
      ierr = BVScaleColumn(pep->V,nqt,1.0/norm);CHKERRQ(ierr);
      nqt++;
    }
    for (i=0;i<nqt;i++) *(S+(j+1)*lds+offq+i) = *(S+j*lds+i)+sigma*(*(S+(j+1)*lds+i));
    ierr = BVSetActiveColumns(pep->V,0,nqt);CHKERRQ(ierr);
    ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
    ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);

    /* level-2 orthogonalization */
    ierr = BVOrthogonalizeColumn(ctx->V,j+1,y,&norm,&lindep);CHKERRQ(ierr);
    a[j] = PetscRealPart(y[j]);
    omega[j+1] = (norm > 0)?1.0:-1.0;
    ierr = BVScaleColumn(ctx->V,j+1,1.0/norm);CHKERRQ(ierr);
    b[j] = PetscAbsReal(norm);

    /* check symmetry */
    ierr = DSGetArrayReal(pep->ds,DS_MAT_T,&f);CHKERRQ(ierr);
    if (j==k) {
      for (i=l;i<j-1;i++) y[i] = PetscAbsScalar(y[i])-PetscAbsReal(f[2*ldds+i]);
      for (i=0;i<l;i++) y[i] = 0.0;
    }
    ierr = DSRestoreArrayReal(pep->ds,DS_MAT_T,&f);CHKERRQ(ierr);
    if (j>0) y[j-1] = PetscAbsScalar(y[j-1])-PetscAbsReal(b[j-1]);
    ierr = PetscBLASIntCast(j,&j_);CHKERRQ(ierr);
    sym = SlepcAbs(BLASnrm2_(&j_,y,&one),sym);
    fro = SlepcAbs(fro,SlepcAbs(a[j],b[j]));
    if (j>0) fro = SlepcAbs(fro,b[j-1]);
    if (sym/fro>PetscMax(PETSC_SQRT_MACHINE_EPSILON,10*pep->tol)) {
      *symmlost = PETSC_TRUE;
      *M=j;
      break;
    }
  }
  ierr = BVSetActiveColumns(pep->V,lock,nqt);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(ctx->V,0,*M);CHKERRQ(ierr);
  ierr = PetscFree(y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PEPSTOAR_QSlice(PEP pep,Mat B)
{
  PetscErrorCode ierr;
  PEP_STOAR      *ctx = (PEP_STOAR*)pep->data;
  PetscInt       j,k,l,nv=0,ld,ldds,t,nq=0,idx;
  PetscInt       nconv=0,deg=pep->nmat-1,count0=0,count1=0;
  PetscScalar    *Q,*om,sigma,*back,*S,*pQ;
  PetscReal      beta,norm=1.0,*omega,*a,*b,*r,eta,lambda;
  PetscBool      breakdown,symmlost=PETSC_FALSE,sinv,falselock=PETSC_TRUE;
  Mat            MS,MQ;
  Vec            v,vomega;
  PEP_SR         sr;
  BVOrthogType   otype;
  BVOrthogBlockType obtype;

  PetscFunctionBegin;
  /* Resize if needed for deflating vectors  */
  sr = ctx->sr;
  sigma = sr->sPres->value;
  k = sr->ndef0+sr->ndef1;
  pep->ncv = ctx->ncv+k;
  pep->nev = ctx->nev+k;
  ierr = PEPAllocateSolution(pep,3);CHKERRQ(ierr);
  ierr = BVDestroy(&ctx->V);CHKERRQ(ierr);
  ierr = BVCreateTensor(pep->V,pep->nmat-1,&ctx->V);CHKERRQ(ierr);
  ierr = BVGetOrthogonalization(pep->V,&otype,NULL,&eta,&obtype);CHKERRQ(ierr);
  ierr = BVSetOrthogonalization(ctx->V,otype,BV_ORTHOG_REFINE_ALWAYS,eta,obtype);CHKERRQ(ierr);
  ierr = DSAllocate(pep->ds,pep->ncv+2);CHKERRQ(ierr);
  ierr = PetscMalloc1(pep->ncv,&back);CHKERRQ(ierr);
  ierr = DSGetLeadingDimension(pep->ds,&ldds);CHKERRQ(ierr);
  ierr = BVSetMatrix(ctx->V,B,PETSC_TRUE);CHKERRQ(ierr);
  if (ctx->lock) {
    /* undocumented option to use a cheaper locking instead of the true locking */
    ierr = PetscOptionsGetBool(NULL,NULL,"-pep_stoar_falselocking",&falselock,NULL);CHKERRQ(ierr);
  } else SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_SUP,"A locking variant is needed for spectrum slicing");
  ierr = PetscObjectTypeCompare((PetscObject)pep->st,STSINVERT,&sinv);CHKERRQ(ierr);
  ierr = RGPushScale(pep->rg,sinv?pep->sfactor:1.0/pep->sfactor);CHKERRQ(ierr);
  ierr = STScaleShift(pep->st,sinv?pep->sfactor:1.0/pep->sfactor);CHKERRQ(ierr);

  /* Get the starting Arnoldi vector */
  ierr = BVSetActiveColumns(pep->V,0,1);CHKERRQ(ierr);
  ierr = BVTensorBuildFirstColumn(ctx->V,pep->nini);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(ctx->V,0,1);CHKERRQ(ierr);
  if (k) {
    /* Insert deflated vectors */
    ierr = BVSetActiveColumns(pep->V,0,0);CHKERRQ(ierr);
    idx = sr->ndef0?sr->idxDef0[0]:sr->idxDef1[0];
    for (j=0;j<k;j++) {
      ierr = BVGetColumn(pep->V,j,&v);CHKERRQ(ierr);
      ierr = BVCopyVec(sr->V,sr->qinfo[idx].q[j],v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(pep->V,j,&v);CHKERRQ(ierr);
    }
    /* Update innerproduct matrix */
    ierr = BVSetActiveColumns(ctx->V,0,0);CHKERRQ(ierr);
    ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(pep->V,0,k);CHKERRQ(ierr);
    ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);

    ierr = BVGetSizes(pep->V,NULL,NULL,&ld);CHKERRQ(ierr);
    ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
    for (j=0;j<sr->ndef0;j++) {
      ierr = PetscArrayzero(S+j*ld*deg,ld*deg);CHKERRQ(ierr);
      ierr = PetscArraycpy(S+j*ld*deg,sr->S+sr->idxDef0[j]*sr->ld*deg,k);CHKERRQ(ierr);
      ierr = PetscArraycpy(S+j*ld*deg+ld,sr->S+sr->idxDef0[j]*sr->ld*deg+sr->ld,k);CHKERRQ(ierr);
      pep->eigr[j] = sr->eigr[sr->idxDef0[j]];
      pep->errest[j] = sr->errest[sr->idxDef0[j]];
    }
    for (j=0;j<sr->ndef1;j++) {
      ierr = PetscArrayzero(S+(j+sr->ndef0)*ld*deg,ld*deg);CHKERRQ(ierr);
      ierr = PetscArraycpy(S+(j+sr->ndef0)*ld*deg,sr->S+sr->idxDef1[j]*sr->ld*deg,k);CHKERRQ(ierr);
      ierr = PetscArraycpy(S+(j+sr->ndef0)*ld*deg+ld,sr->S+sr->idxDef1[j]*sr->ld*deg+sr->ld,k);CHKERRQ(ierr);
      pep->eigr[j+sr->ndef0] = sr->eigr[sr->idxDef1[j]];
      pep->errest[j+sr->ndef0] = sr->errest[sr->idxDef1[j]];
    }
    ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
    ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(ctx->V,0,k+1);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,k+1,&vomega);CHKERRQ(ierr);
    ierr = VecGetArray(vomega,&om);CHKERRQ(ierr);
    for (j=0;j<k;j++) {
      ierr = BVOrthogonalizeColumn(ctx->V,j,NULL,&norm,NULL);CHKERRQ(ierr);
      ierr = BVScaleColumn(ctx->V,j,1/norm);CHKERRQ(ierr);
      om[j] = (norm>=0.0)?1.0:-1.0;
    }
    ierr = BVTensorGetFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = MatDenseGetArray(MS,&S);CHKERRQ(ierr);
    for (j=0;j<deg;j++) {
      ierr = BVSetRandomColumn(pep->V,k+j);CHKERRQ(ierr);
      ierr = BVOrthogonalizeColumn(pep->V,k+j,S+k*ld*deg+j*ld,&norm,NULL);CHKERRQ(ierr);
      ierr = BVScaleColumn(pep->V,k+j,1.0/norm);CHKERRQ(ierr);
      S[k*ld*deg+j*ld+k+j] = norm;
    }
    ierr = MatDenseRestoreArray(MS,&S);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(pep->V,0,k+deg);CHKERRQ(ierr);
    ierr = BVTensorRestoreFactors(ctx->V,NULL,&MS);CHKERRQ(ierr);
    ierr = BVOrthogonalizeColumn(ctx->V,k,NULL,&norm,NULL);CHKERRQ(ierr);
    ierr = BVScaleColumn(ctx->V,k,1.0/norm);CHKERRQ(ierr);
    om[k] = (norm>=0.0)?1.0:-1.0;
    ierr = VecRestoreArray(vomega,&om);CHKERRQ(ierr);
    ierr = BVSetSignature(ctx->V,vomega);CHKERRQ(ierr);
    ierr = DSGetArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
    ierr = VecGetArray(vomega,&om);CHKERRQ(ierr);
    for (j=0;j<k;j++) a[j] = PetscRealPart(om[j]/(pep->eigr[j]-sigma));
    ierr = VecRestoreArray(vomega,&om);CHKERRQ(ierr);
    ierr = VecDestroy(&vomega);CHKERRQ(ierr);
    ierr = DSRestoreArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
    ierr = DSGetArray(pep->ds,DS_MAT_Q,&pQ);CHKERRQ(ierr);
    ierr = PetscArrayzero(pQ,ldds*k);CHKERRQ(ierr);
    for (j=0;j<k;j++) pQ[j+j*ldds] = 1.0;
    ierr = DSRestoreArray(pep->ds,DS_MAT_Q,&pQ);CHKERRQ(ierr);
  }
  ierr = BVSetActiveColumns(ctx->V,0,k+1);CHKERRQ(ierr);
  ierr = DSGetArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF,k+1,&vomega);CHKERRQ(ierr);
  ierr = BVGetSignature(ctx->V,vomega);CHKERRQ(ierr);
  ierr = VecGetArray(vomega,&om);CHKERRQ(ierr);
  for (j=0;j<k+1;j++) omega[j] = PetscRealPart(om[j]);
  ierr = VecRestoreArray(vomega,&om);CHKERRQ(ierr);
  ierr = DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
  ierr = VecDestroy(&vomega);CHKERRQ(ierr);

  ierr = PetscInfo7(pep,"Start STOAR: sigma=%g in [%g,%g], for deflation: left=%D right=%D, searching: left=%D right=%D\n",(double)sr->sPres->value,(double)(sr->sPres->neighb[0]?sr->sPres->neighb[0]->value:sr->int0),(double)(sr->sPres->neighb[1]?sr->sPres->neighb[1]->value:sr->int1),sr->ndef0,sr->ndef1,sr->sPres->nsch[0],sr->sPres->nsch[1]);CHKERRQ(ierr);

  /* Restart loop */
  l = 0;
  pep->nconv = k;
  while (pep->reason == PEP_CONVERGED_ITERATING) {
    pep->its++;
    ierr = DSGetArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
    b = a+ldds;
    ierr = DSGetArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);

    /* Compute an nv-step Lanczos factorization */
    nv = PetscMin(pep->nconv+pep->mpd,pep->ncv);
    ierr = PEPSTOARrun_QSlice(pep,a,b,omega,pep->nconv+l,&nv,&breakdown,&symmlost,pep->work);CHKERRQ(ierr);
    beta = b[nv-1];
    if (symmlost && nv==pep->nconv+l) {
      pep->reason = PEP_DIVERGED_SYMMETRY_LOST;
      pep->nconv = nconv;
      ierr = PetscInfo2(pep,"Symmetry lost in STOAR sigma=%g nconv=%D\n",(double)sr->sPres->value,nconv);CHKERRQ(ierr);
      if (falselock || !ctx->lock) {
        ierr = BVSetActiveColumns(ctx->V,0,pep->nconv);CHKERRQ(ierr);
        ierr = BVTensorCompress(ctx->V,0);CHKERRQ(ierr);
      }
      break;
    }
    ierr = DSRestoreArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
    ierr = DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
    ierr = DSSetDimensions(pep->ds,nv,0,pep->nconv,pep->nconv+l);CHKERRQ(ierr);
    if (l==0) {
      ierr = DSSetState(pep->ds,DS_STATE_INTERMEDIATE);CHKERRQ(ierr);
    } else {
      ierr = DSSetState(pep->ds,DS_STATE_RAW);CHKERRQ(ierr);
    }

    /* Solve projected problem */
    ierr = DSSolve(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);
    ierr = DSSort(pep->ds,pep->eigr,pep->eigi,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DSSynchronize(pep->ds,pep->eigr,pep->eigi);CHKERRQ(ierr);

    /* Check convergence */
    /* ierr = PEPSTOARpreKConvergence(pep,nv,&norm,pep->work);CHKERRQ(ierr);*/
    norm = 1.0;
    ierr = DSGetDimensions(pep->ds,NULL,NULL,NULL,NULL,&t);CHKERRQ(ierr);
    ierr = PEPKrylovConvergence(pep,PETSC_FALSE,pep->nconv,t-pep->nconv,PetscAbsReal(beta)*norm,&k);CHKERRQ(ierr);
    ierr = (*pep->stopping)(pep,pep->its,pep->max_it,k,pep->nev,&pep->reason,pep->stoppingctx);CHKERRQ(ierr);
    for (j=0;j<k;j++) back[j] = pep->eigr[j];
    ierr = STBackTransform(pep->st,k,back,pep->eigi);CHKERRQ(ierr);
    count0=count1=0;
    for (j=0;j<k;j++) {
      lambda = PetscRealPart(back[j]);
      if (((sr->dir)*(sr->sPres->value - lambda) > 0) && ((sr->dir)*(lambda - sr->sPres->ext[0]) > 0)) count0++;
      if (((sr->dir)*(lambda - sr->sPres->value) > 0) && ((sr->dir)*(sr->sPres->ext[1] - lambda) > 0)) count1++;
    }
    if ((count0-sr->ndef0 >= sr->sPres->nsch[0]) && (count1-sr->ndef1 >= sr->sPres->nsch[1])) pep->reason = PEP_CONVERGED_TOL;
    /* Update l */
    if (pep->reason != PEP_CONVERGED_ITERATING || breakdown) l = 0;
    else {
      l = PetscMax(1,(PetscInt)((nv-k)/2));
      l = PetscMin(l,t);
      if (!breakdown) {
        ierr = DSGetArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
        if (*(a+ldds+k+l-1)!=0) {
          if (k+l<nv-1) l = l+1;
          else l = l-1;
        }
        /* Prepare the Rayleigh quotient for restart */
        ierr = DSGetArray(pep->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
        ierr = DSGetArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
        r = a + 2*ldds;
        for (j=k;j<k+l;j++) {
          r[j] = PetscRealPart(Q[nv-1+j*ldds]*beta);
        }
        b = a+ldds;
        b[k+l-1] = r[k+l-1];
        omega[k+l] = omega[nv];
        ierr = DSRestoreArray(pep->ds,DS_MAT_Q,&Q);CHKERRQ(ierr);
        ierr = DSRestoreArrayReal(pep->ds,DS_MAT_T,&a);CHKERRQ(ierr);
        ierr = DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
      }
    }
    nconv = k;
    if (!ctx->lock && pep->reason == PEP_CONVERGED_ITERATING && !breakdown) { l += k; k = 0; } /* non-locking variant: reset no. of converged pairs */

    /* Update S */
    ierr = DSGetMat(pep->ds,DS_MAT_Q,&MQ);CHKERRQ(ierr);
    ierr = BVMultInPlace(ctx->V,MQ,pep->nconv,k+l);CHKERRQ(ierr);
    ierr = MatDestroy(&MQ);CHKERRQ(ierr);

    /* Copy last column of S */
    ierr = BVCopyColumn(ctx->V,nv,k+l);CHKERRQ(ierr);
    ierr = DSGetArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,k+l,&vomega);CHKERRQ(ierr);
    ierr = VecGetArray(vomega,&om);CHKERRQ(ierr);
    for (j=0;j<k+l;j++) om[j] = omega[j];
    ierr = VecRestoreArray(vomega,&om);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(ctx->V,0,k+l);CHKERRQ(ierr);
    ierr = BVSetSignature(ctx->V,vomega);CHKERRQ(ierr);
    ierr = VecDestroy(&vomega);CHKERRQ(ierr);
    ierr = DSRestoreArrayReal(pep->ds,DS_MAT_D,&omega);CHKERRQ(ierr);

    if (breakdown && pep->reason == PEP_CONVERGED_ITERATING) {
      /* stop if breakdown */
      ierr = PetscInfo2(pep,"Breakdown TOAR method (it=%D norm=%g)\n",pep->its,(double)beta);CHKERRQ(ierr);
      pep->reason = PEP_DIVERGED_BREAKDOWN;
    }
    if (pep->reason != PEP_CONVERGED_ITERATING) l--;
    ierr = BVGetActiveColumns(pep->V,NULL,&nq);CHKERRQ(ierr);
    if (k+l+deg<=nq) {
      ierr = BVSetActiveColumns(ctx->V,pep->nconv,k+l+1);CHKERRQ(ierr);
      if (!falselock && ctx->lock) {
        ierr = BVTensorCompress(ctx->V,k-pep->nconv);CHKERRQ(ierr);
      } else {
        ierr = BVTensorCompress(ctx->V,0);CHKERRQ(ierr);
      }
    }
    pep->nconv = k;
    ierr = PEPMonitor(pep,pep->its,nconv,pep->eigr,pep->eigi,pep->errest,nv);CHKERRQ(ierr);
  }
  sr->itsKs += pep->its;
  if (pep->nconv>0) {
    ierr = BVSetActiveColumns(ctx->V,0,pep->nconv);CHKERRQ(ierr);
    ierr = BVGetActiveColumns(pep->V,NULL,&nq);CHKERRQ(ierr);
    ierr = BVSetActiveColumns(pep->V,0,nq);CHKERRQ(ierr);
    if (nq>pep->nconv) {
      ierr = BVTensorCompress(ctx->V,pep->nconv);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(pep->V,0,pep->nconv);CHKERRQ(ierr);
    }
    for (j=0;j<pep->nconv;j++) {
      pep->eigr[j] *= pep->sfactor;
      pep->eigi[j] *= pep->sfactor;
    }
  }
  ierr = PetscInfo4(pep,"Finished STOAR: nconv=%D (deflated=%D, left=%D, right=%D)\n",pep->nconv,sr->ndef0+sr->ndef1,count0-sr->ndef0,count1-sr->ndef1);CHKERRQ(ierr);
  ierr = STScaleShift(pep->st,sinv?1.0/pep->sfactor:pep->sfactor);CHKERRQ(ierr);
  ierr = RGPopScale(pep->rg);CHKERRQ(ierr);

  if (pep->reason == PEP_DIVERGED_SYMMETRY_LOST && nconv<sr->ndef0+sr->ndef1) SETERRQ1(PetscObjectComm((PetscObject)pep),1,"Symmetry lost at sigma=%g",(double)sr->sPres->value);
  if (pep->reason == PEP_DIVERGED_SYMMETRY_LOST && nconv==sr->ndef0+sr->ndef1) {
    if (++sr->symmlost>10) SETERRQ1(PetscObjectComm((PetscObject)pep),1,"Symmetry lost at sigma=%g",(double)sr->sPres->value);
  } else sr->symmlost = 0;

  /* truncate Schur decomposition and change the state to raw so that
     DSVectors() computes eigenvectors from scratch */
  ierr = DSSetDimensions(pep->ds,pep->nconv,0,0,0);CHKERRQ(ierr);
  ierr = DSSetState(pep->ds,DS_STATE_RAW);CHKERRQ(ierr);
  ierr = PetscFree(back);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define SWAP(a,b,t) {t=a;a=b;b=t;}

static PetscErrorCode PEPQSliceGetInertias(PEP pep,PetscInt *n,PetscReal **shifts,PetscInt **inertias)
{
  PetscErrorCode  ierr;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR          sr=ctx->sr;
  PetscInt        i=0,j,tmpi;
  PetscReal       v,tmpr;
  PEP_shift       s;

  PetscFunctionBegin;
  if (!pep->state) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONGSTATE,"Must call PEPSetUp() first");
  if (!sr) SETERRQ(PetscObjectComm((PetscObject)pep),PETSC_ERR_ARG_WRONGSTATE,"Only available in interval computations, see PEPSetInterval()");
  if (!sr->s0) {  /* PEPSolve not called yet */
    *n = 2;
  } else {
    *n = 1;
    s = sr->s0;
    while (s) {
      (*n)++;
      s = s->neighb[1];
    }
  }
  ierr = PetscMalloc1(*n,shifts);CHKERRQ(ierr);
  ierr = PetscMalloc1(*n,inertias);CHKERRQ(ierr);
  if (!sr->s0) {  /* PEPSolve not called yet */
    (*shifts)[0]   = sr->int0;
    (*shifts)[1]   = sr->int1;
    (*inertias)[0] = sr->inertia0;
    (*inertias)[1] = sr->inertia1;
  } else {
    s = sr->s0;
    while (s) {
      (*shifts)[i]     = s->value;
      (*inertias)[i++] = s->inertia;
      s = s->neighb[1];
    }
    (*shifts)[i]   = sr->int1;
    (*inertias)[i] = sr->inertia1;
  }
  /* remove possible duplicate in last position */
  if ((*shifts)[(*n)-1]==(*shifts)[(*n)-2]) (*n)--;
  /* sort result */
  for (i=0;i<*n;i++) {
    v = (*shifts)[i];
    for (j=i+1;j<*n;j++) {
      if (v > (*shifts)[j]) {
        SWAP((*shifts)[i],(*shifts)[j],tmpr);
        SWAP((*inertias)[i],(*inertias)[j],tmpi);
        v = (*shifts)[i];
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PEPSolve_STOAR_QSlice(PEP pep)
{
  PetscErrorCode ierr;
  PetscInt       i,j,ti,deg=pep->nmat-1;
  PetscReal      newS;
  PEP_STOAR      *ctx=(PEP_STOAR*)pep->data;
  PEP_SR         sr=ctx->sr;
  Mat            S,B;
  PetscScalar    *pS;

  PetscFunctionBegin;
  ierr = PetscCitationsRegister(citation,&cited);CHKERRQ(ierr);

  /* Only with eigenvalues present in the interval ...*/
  if (sr->numEigs==0) {
    pep->reason = PEP_CONVERGED_TOL;
    PetscFunctionReturn(0);
  }

  /* Inner product matrix */
  ierr = PEPSTOARSetUpInnerMatrix(pep,&B);CHKERRQ(ierr);

  /* Array of pending shifts */
  sr->maxPend = 100; /* Initial size */
  sr->nPend = 0;
  ierr = PetscMalloc1(sr->maxPend,&sr->pending);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)pep,(sr->maxPend)*sizeof(PEP_shift));CHKERRQ(ierr);
  ierr = PEPCreateShift(pep,sr->int0,NULL,NULL);CHKERRQ(ierr);
  /* extract first shift */
  sr->sPrev = NULL;
  sr->sPres = sr->pending[--sr->nPend];
  sr->sPres->inertia = sr->inertia0;
  pep->target = sr->sPres->value;
  sr->s0 = sr->sPres;
  sr->indexEig = 0;

  /* Memory reservation for auxiliary variables */
  ierr = PetscLogObjectMemory((PetscObject)pep,(sr->numEigs+2*pep->ncv)*sizeof(PetscScalar));CHKERRQ(ierr);
  for (i=0;i<sr->numEigs;i++) {
    sr->eigr[i]   = 0.0;
    sr->eigi[i]   = 0.0;
    sr->errest[i] = 0.0;
    sr->perm[i]   = i;
  }
  /* Vectors for deflation */
  ierr = PetscMalloc2(sr->numEigs,&sr->idxDef0,sr->numEigs,&sr->idxDef1);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)pep,sr->numEigs*sizeof(PetscInt));CHKERRQ(ierr);
  sr->indexEig = 0;
  while (sr->sPres) {
    /* Search for deflation */
    ierr = PEPLookForDeflation(pep);CHKERRQ(ierr);
    /* KrylovSchur */
    ierr = PEPSTOAR_QSlice(pep,B);CHKERRQ(ierr);

    ierr = PEPStoreEigenpairs(pep);CHKERRQ(ierr);
    /* Select new shift */
    if (!sr->sPres->comp[1]) {
      ierr = PEPGetNewShiftValue(pep,1,&newS);CHKERRQ(ierr);
      ierr = PEPCreateShift(pep,newS,sr->sPres,sr->sPres->neighb[1]);CHKERRQ(ierr);
    }
    if (!sr->sPres->comp[0]) {
      /* Completing earlier interval */
      ierr = PEPGetNewShiftValue(pep,0,&newS);CHKERRQ(ierr);
      ierr = PEPCreateShift(pep,newS,sr->sPres->neighb[0],sr->sPres);CHKERRQ(ierr);
    }
    /* Preparing for a new search of values */
    ierr = PEPExtractShift(pep);CHKERRQ(ierr);
  }

  /* Updating pep values prior to exit */
  ierr = PetscFree2(sr->idxDef0,sr->idxDef1);CHKERRQ(ierr);
  ierr = PetscFree(sr->pending);CHKERRQ(ierr);
  pep->nconv  = sr->indexEig;
  pep->reason = PEP_CONVERGED_TOL;
  pep->its    = sr->itsKs;
  pep->nev    = sr->indexEig;
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,pep->nconv,pep->nconv,NULL,&S);CHKERRQ(ierr);
  ierr = MatDenseGetArray(S,&pS);CHKERRQ(ierr);
  for (i=0;i<pep->nconv;i++) {
    for (j=0;j<sr->qinfo[i].nq;j++) pS[i*pep->nconv+sr->qinfo[i].q[j]] = *(sr->S+i*sr->ld*deg+j);
  }
  ierr = MatDenseRestoreArray(S,&pS);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(sr->V,0,pep->nconv);CHKERRQ(ierr);
  ierr = BVMultInPlace(sr->V,S,0,pep->nconv);CHKERRQ(ierr);
  ierr = MatDestroy(&S);CHKERRQ(ierr);
  ierr = BVDestroy(&pep->V);CHKERRQ(ierr);
  pep->V = sr->V;
  ierr = PetscFree4(pep->eigr,pep->eigi,pep->errest,pep->perm);CHKERRQ(ierr);
  pep->eigr   = sr->eigr;
  pep->eigi   = sr->eigi;
  pep->perm   = sr->perm;
  pep->errest = sr->errest;
  if (sr->dir<0) {
    for (i=0;i<pep->nconv/2;i++) {
      ti = sr->perm[i]; sr->perm[i] = sr->perm[pep->nconv-1-i]; sr->perm[pep->nconv-1-i] = ti;
    }
  }
  ierr = PetscFree(ctx->inertias);CHKERRQ(ierr);
  ierr = PetscFree(ctx->shifts);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = PEPQSliceGetInertias(pep,&ctx->nshifts,&ctx->shifts,&ctx->inertias);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

