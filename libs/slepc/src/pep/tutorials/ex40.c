/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Checking the definite property in quadratic symmetric eigenproblem.\n\n"
  "The command line options are:\n"
  "  -n <n> ... dimension of the matrices.\n"
  "  -transform... whether to transform to a hyperbolic problem or not.\n"
  "  -nonhyperbolic... to test with a modified (definite) problem that is not hyperbolic.\n\n";

#include <slepcpep.h>

/*
  This example is based on spring.c, for fixed values mu=1,tau=10,kappa=5

  The transformations are based on the method proposed in [Niendorf and Voss, LAA 2010].
*/

PetscErrorCode QEPDefiniteTransformGetMatrices(PEP,PetscBool,PetscReal,PetscReal,Mat[3]);
PetscErrorCode QEPDefiniteTransformMap(PetscBool,PetscReal,PetscReal,PetscInt,PetscScalar*,PetscBool);
PetscErrorCode QEPDefiniteCheckError(Mat*,PEP,PetscBool,PetscReal,PetscReal);
PetscErrorCode TransformMatricesMoebius(Mat[3],MatStructure,PetscReal,PetscReal,PetscReal,PetscReal,Mat[3]);

int main(int argc,char **argv)
{
  Mat            M,C,K,*Op,A[3],At[3],B[3]; /* problem matrices */
  PEP            pep;        /* polynomial eigenproblem solver context */
  ST             st;         /* spectral transformation context */
  KSP            ksp;
  PC             pc;
  PEPProblemType type;
  PetscBool      terse,transform=PETSC_FALSE,nohyp=PETSC_FALSE;
  PetscInt       n=100,Istart,Iend,i,def=0,hyp;
  PetscReal      muu=1,tau=10,kappa=5,inta,intb;
  PetscReal      alpha,beta,xi,mu,at[2]={0.0,0.0},c=.857,s;
  PetscScalar    target,targett,ats[2];
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPEP example that checks definite property, n=%D\n\n",n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the matrices that define the eigensystem, (k^2*M+k*C+K)x=0
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* K is a tridiagonal */
  ierr = MatCreate(PETSC_COMM_WORLD,&K);CHKERRQ(ierr);
  ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(K);CHKERRQ(ierr);
  ierr = MatSetUp(K);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(K,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) {
      ierr = MatSetValue(K,i,i-1,-kappa,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatSetValue(K,i,i,kappa*3.0,INSERT_VALUES);CHKERRQ(ierr);
    if (i<n-1) {
      ierr = MatSetValue(K,i,i+1,-kappa,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* C is a tridiagonal */
  ierr = MatCreate(PETSC_COMM_WORLD,&C);CHKERRQ(ierr);
  ierr = MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(C);CHKERRQ(ierr);
  ierr = MatSetUp(C);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(C,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    if (i>0) {
      ierr = MatSetValue(C,i,i-1,-tau,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatSetValue(C,i,i,tau*3.0,INSERT_VALUES);CHKERRQ(ierr);
    if (i<n-1) {
      ierr = MatSetValue(C,i,i+1,-tau,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* M is a diagonal matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&M);CHKERRQ(ierr);
  ierr = MatSetSizes(M,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(M);CHKERRQ(ierr);
  ierr = MatSetUp(M);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(M,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(M,i,i,muu,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(NULL,NULL,"-nonhyperbolic",&nohyp,NULL);CHKERRQ(ierr);
  A[0] = K; A[1] = C; A[2] = M;
  if (nohyp) {
    s = c*.6;
    ierr = TransformMatricesMoebius(A,DIFFERENT_NONZERO_PATTERN,c,s,-s,c,At);CHKERRQ(ierr);
    for (i=0;i<3;i++) { ierr = MatDestroy(&A[i]);CHKERRQ(ierr); }
    Op = At;
  } else Op = A;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
  ierr = PEPSetProblemType(pep,PEP_HERMITIAN);CHKERRQ(ierr);
  ierr = PEPSetType(pep,PEPSTOAR);CHKERRQ(ierr);
  /*
     Set operators and set problem type
  */
  ierr = PEPSetOperators(pep,3,Op);CHKERRQ(ierr);

  /*
     Set shift-and-invert with Cholesky; select MUMPS if available
  */
  ierr = PEPGetST(pep,&st);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCCHOLESKY);CHKERRQ(ierr);
#if defined(PETSC_HAVE_MUMPS)
#if defined(PETSC_USE_COMPLEX)
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Spectrum slicing with MUMPS is not available for complex scalars");
#endif
  ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);CHKERRQ(ierr);
  /*
     Add several MUMPS options (see ex43.c for a better way of setting them in program):
     '-mat_mumps_icntl_13 1': turn off ScaLAPACK for matrix inertia
  */
  ierr = PetscOptionsInsertString(NULL,"-mat_mumps_icntl_13 1 -mat_mumps_icntl_24 1 -mat_mumps_cntl_3 1e-12");CHKERRQ(ierr);
#endif

  /*
     Set solver parameters at runtime
  */
  ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(NULL,NULL,"-transform",&transform,NULL);CHKERRQ(ierr);
  if (transform) {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Check if the problem is definite
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = PEPCheckDefiniteQEP(pep,&xi,&mu,&def,&hyp);CHKERRQ(ierr);
    switch (def) {
      case 1:
        if (hyp==1) {ierr = PetscPrintf(PETSC_COMM_WORLD,"Hyperbolic Problem xi=%g\n",xi);CHKERRQ(ierr);}
        else {ierr = PetscPrintf(PETSC_COMM_WORLD,"Definite Problem xi=%g mu=%g\n",xi,mu);CHKERRQ(ierr);}
        break;
      case -1:
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Not Definite Problem\n");CHKERRQ(ierr);
        break;
      default:
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Cannot determine definiteness\n");CHKERRQ(ierr);
        break;
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Transform the QEP to have a definite inner product in the linearization
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    if (def==1) {
      ierr = QEPDefiniteTransformGetMatrices(pep,hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu,B);CHKERRQ(ierr);
      ierr = PEPSetOperators(pep,3,B);CHKERRQ(ierr);
      ierr = PEPGetTarget(pep,&target);CHKERRQ(ierr);
      targett = target;
      ierr = QEPDefiniteTransformMap(hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu,1,&targett,PETSC_FALSE);CHKERRQ(ierr);
      ierr = PEPSetTarget(pep,targett);CHKERRQ(ierr);
      ierr = PEPGetProblemType(pep,&type);CHKERRQ(ierr);
      ierr = PEPSetProblemType(pep,PEP_HYPERBOLIC);CHKERRQ(ierr);
      ierr = PEPSTOARGetLinearization(pep,&alpha,&beta);CHKERRQ(ierr);
      ierr = PEPSTOARSetLinearization(pep,1.0,0.0);CHKERRQ(ierr);
      ierr = PEPGetInterval(pep,&inta,&intb);CHKERRQ(ierr);
      if (inta!=intb) {
        ats[0] = inta; ats[1] = intb;
        ierr = QEPDefiniteTransformMap(hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu,2,ats,PETSC_FALSE);CHKERRQ(ierr);
        at[0] = PetscRealPart(ats[0]); at[1] = PetscRealPart(ats[1]);
        if (at[0]<at[1]) { ierr = PEPSetInterval(pep,at[0],at[1]);CHKERRQ(ierr); }
        else { ierr = PEPSetInterval(pep,PETSC_MIN_REAL,at[1]);CHKERRQ(ierr); }
      }
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PEPSolve(pep);CHKERRQ(ierr);

  /* show detailed info unless -terse option is given by user */
  if (def!=1) {
    ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
    if (terse) {
      ierr = PEPErrorView(pep,PEP_ERROR_BACKWARD,NULL);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
      ierr = PEPReasonView(pep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      ierr = PEPErrorView(pep,PEP_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
  } else {
    /* Map the solution */
    ierr = PEPReasonView(pep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = QEPDefiniteCheckError(Op,pep,hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu);CHKERRQ(ierr);
    for (i=0;i<3;i++) {ierr = MatDestroy(B+i);CHKERRQ(ierr);}
  }
  if (at[0]>at[1]) {
    ierr = PEPSetInterval(pep,at[0],PETSC_MAX_REAL);CHKERRQ(ierr);
    ierr = PEPSolve(pep);CHKERRQ(ierr);
    ierr = PEPReasonView(pep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    /* Map the solution */
    ierr = QEPDefiniteCheckError(Op,pep,hyp==1?PETSC_TRUE:PETSC_FALSE,xi,mu);CHKERRQ(ierr);
  }
  if (def==1) {
    ierr = PEPSetTarget(pep,target);CHKERRQ(ierr);
    ierr = PEPSetProblemType(pep,type);CHKERRQ(ierr);
    ierr = PEPSTOARSetLinearization(pep,alpha,beta);CHKERRQ(ierr);
    if (inta!=intb) {
      ierr = PEPSetInterval(pep,inta,intb);CHKERRQ(ierr);
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PEPDestroy(&pep);CHKERRQ(ierr);
  for (i=0;i<3;i++) { ierr = MatDestroy(Op+i);CHKERRQ(ierr); }
  ierr = SlepcFinalize();
  return ierr;
}

/* ------------------------------------------------------------------- */
/*
  QEPDefiniteTransformMap_Initial - map a scalar value with a certain Moebius transform

                   a theta + b
         lambda = --------------
                   c theta + d

  Input:
    xi,mu: real values such that Q(xi)<0 and Q(mu)>0
    hyperbolic: if true the problem is assumed hyperbolic (mu is not used)
  Input/Output:
    val (array of length n)
    if backtransform=true returns lambda from theta, else returns theta from lambda
*/
static PetscErrorCode QEPDefiniteTransformMap_Initial(PetscBool hyperbolic,PetscReal xi,PetscReal mu,PetscInt n,PetscScalar *val,PetscBool backtransform)
{
  PetscInt  i;
  PetscReal a,b,c,d,s;

  PetscFunctionBegin;
  if (hyperbolic) { a = 1.0; b = xi; c =0.0; d = 1.0; }
  else { a = mu; b = mu*xi-1; c = 1.0; d = xi+mu; }
  if (!backtransform) { s = a; a = -d; d = -s; }
  for (i=0;i<n;i++) {
    if (PetscRealPart(val[i]) >= PETSC_MAX_REAL || PetscRealPart(val[i]) <= PETSC_MIN_REAL) val[i] = a/c;
    else if (val[i] == -d/c) val[i] = PETSC_MAX_REAL;
    else val[i] = (a*val[i]+b)/(c*val[i]+d);
  }
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
  QEPDefiniteTransformMap - perform the mapping if the problem is hyperbolic, otherwise
  modify the value of xi in advance
*/
PetscErrorCode QEPDefiniteTransformMap(PetscBool hyperbolic,PetscReal xi,PetscReal mu,PetscInt n,PetscScalar *val,PetscBool backtransform)
{
  PetscReal      xit;
  PetscScalar    alpha;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  xit = xi;
  if (!hyperbolic) {
    alpha = xi;
    ierr = QEPDefiniteTransformMap_Initial(PETSC_FALSE,0.0,mu,1,&alpha,PETSC_FALSE);CHKERRQ(ierr);
    xit = PetscRealPart(alpha);
  }
  ierr = QEPDefiniteTransformMap_Initial(hyperbolic,xit,mu,n,val,backtransform);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
  TransformMatricesMoebius - transform the coefficient matrices of a QEP

  Input:
    A: coefficient matrices of the original QEP
    a,b,c,d: parameters of the Moebius transform
    str: structure flag for MatAXPY operations
  Output:
    B: transformed matrices
*/
PetscErrorCode TransformMatricesMoebius(Mat A[3],MatStructure str,PetscReal a,PetscReal b,PetscReal c,PetscReal d,Mat B[3])
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  PetscReal      cf[9];

  PetscFunctionBegin;
  for (i=0;i<3;i++) {
    ierr = MatDuplicate(A[2],MAT_COPY_VALUES,&B[i]);CHKERRQ(ierr);
  }
  /* Ct = b*b*A+b*d*B+d*d*C */
  cf[0] = d*d; cf[1] = b*d; cf[2] = b*b;
  /* Bt = 2*a*b*A+(b*c+a*d)*B+2*c*d*C*/
  cf[3] = 2*c*d; cf[4] = b*c+a*d; cf[5] = 2*a*b;
  /* At = a*a*A+a*c*B+c*c*C */
  cf[6] = c*c; cf[7] = a*c; cf[8] = a*a;
  for (k=0;k<3;k++) {
    ierr = MatScale(B[k],cf[k*3+2]);CHKERRQ(ierr);
    for (i=0;i<2;i++) {
      ierr = MatAXPY(B[k],cf[3*k+i],A[i],str);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
  QEPDefiniteTransformGetMatrices - given a PEP of degree 2, transform the three
  matrices with TransformMatricesMoebius

  Input:
    pep: polynomial eigenproblem to be transformed, with Q(.) being the quadratic polynomial
    xi,mu: real values such that Q(xi)<0 and Q(mu)>0
    hyperbolic: if true the problem is assumed hyperbolic (mu is not used)
  Output:
    T: coefficient matrices of the transformed polynomial
*/
PetscErrorCode QEPDefiniteTransformGetMatrices(PEP pep,PetscBool hyperbolic,PetscReal xi,PetscReal mu,Mat T[3])
{
  PetscErrorCode ierr;
  MatStructure   str;
  ST             st;
  PetscInt       i;
  PetscReal      a,b,c,d;
  PetscScalar    xit;
  Mat            A[3];

  PetscFunctionBegin;
  for (i=2;i>=0;i--) {
    ierr = PEPGetOperators(pep,i,&A[i]);CHKERRQ(ierr);
  }
  if (hyperbolic) { a = 1.0; b = xi; c =0.0; d = 1.0; }
  else {
    xit = xi;
    ierr = QEPDefiniteTransformMap_Initial(PETSC_FALSE,0.0,mu,1,&xit,PETSC_FALSE);CHKERRQ(ierr);
    a = mu; b = mu*PetscRealPart(xit)-1.0; c = 1.0; d = PetscRealPart(xit)+mu;
  }
  ierr = PEPGetST(pep,&st);CHKERRQ(ierr);
  ierr = STGetMatStructure(st,&str);CHKERRQ(ierr);
  ierr = TransformMatricesMoebius(A,str,a,b,c,d,T);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
  Auxiliary funciton to compute the residual norm of an eigenpair of a QEP defined
  by coefficient matrices A
*/
static PetscErrorCode PEPResidualNorm(Mat *A,PetscScalar kr,PetscScalar ki,Vec xr,Vec xi,Vec *z,PetscReal *norm)
{
  PetscErrorCode ierr;
  PetscInt       i,nmat=3;
  PetscScalar    vals[3];
  Vec            u,w;
#if !defined(PETSC_USE_COMPLEX)
  Vec            ui,wi;
  PetscReal      ni;
  PetscBool      imag;
  PetscScalar    ivals[3];
#endif

  PetscFunctionBegin;
  u = z[0]; w = z[1];
  ierr = VecSet(u,0.0);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  ui = z[2]; wi = z[3];
#endif
  vals[0] = 1.0;
  vals[1] = kr;
  vals[2] = kr*kr-ki*ki;
#if !defined(PETSC_USE_COMPLEX)
  ivals[0] = 0.0;
  ivals[1] = ki;
  ivals[2] = 2.0*kr*ki;
  if (ki == 0 || PetscAbsScalar(ki) < PetscAbsScalar(kr*PETSC_MACHINE_EPSILON))
    imag = PETSC_FALSE;
  else {
    imag = PETSC_TRUE;
    ierr = VecSet(ui,0.0);CHKERRQ(ierr);
  }
#endif
  for (i=0;i<nmat;i++) {
    if (vals[i]!=0.0) {
      ierr = MatMult(A[i],xr,w);CHKERRQ(ierr);
      ierr = VecAXPY(u,vals[i],w);CHKERRQ(ierr);
    }
#if !defined(PETSC_USE_COMPLEX)
    if (imag) {
      if (ivals[i]!=0 || vals[i]!=0) {
        ierr = MatMult(A[i],xi,wi);CHKERRQ(ierr);
        if (vals[i]==0) {
          ierr = MatMult(A[i],xr,w);CHKERRQ(ierr);
        }
      }
      if (ivals[i]!=0){
        ierr = VecAXPY(u,-ivals[i],wi);CHKERRQ(ierr);
        ierr = VecAXPY(ui,ivals[i],w);CHKERRQ(ierr);
      }
      if (vals[i]!=0) {
        ierr = VecAXPY(ui,vals[i],wi);CHKERRQ(ierr);
      }
    }
#endif
  }
  ierr = VecNorm(u,NORM_2,norm);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  if (imag) {
    ierr = VecNorm(ui,NORM_2,&ni);CHKERRQ(ierr);
    *norm = SlepcAbsEigenvalue(*norm,ni);
  }
#endif
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
  QEPDefiniteCheckError - check and print the residual norm of a transformed PEP

  Input:
    A: coefficient matrices of the original problem
    pep: solver containing the computed solution of the transformed problem
    xi,mu,hyperbolic: parameters used in transformation
*/
PetscErrorCode QEPDefiniteCheckError(Mat *A,PEP pep,PetscBool hyperbolic,PetscReal xi,PetscReal mu)
{
  PetscErrorCode ierr;
  PetscScalar    er,ei;
  PetscReal      re,im,error;
  Vec            vr,vi,w[4];
  PetscInt       i,nconv;
  BV             bv;
#define EXLEN 30
  char           ex[EXLEN],sep[]=" ---------------------- --------------------\n";

  PetscFunctionBegin;
  ierr = PetscSNPrintf(ex,EXLEN,"||P(k)x||/||kx||");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s            k             %s\n%s",sep,ex,sep);CHKERRQ(ierr);
  ierr = PEPGetConverged(pep,&nconv);CHKERRQ(ierr);
  ierr = PEPGetBV(pep,&bv);CHKERRQ(ierr);
  ierr = BVCreateVec(bv,w);CHKERRQ(ierr);
  ierr = VecDuplicate(w[0],&vr);CHKERRQ(ierr);
  ierr = VecDuplicate(w[0],&vi);CHKERRQ(ierr);
  for (i=1;i<4;i++) {ierr = VecDuplicate(w[0],w+i);CHKERRQ(ierr);}
  for (i=0;i<nconv;i++) {
    ierr = PEPGetEigenpair(pep,i,&er,&ei,vr,vi);CHKERRQ(ierr);
    ierr = QEPDefiniteTransformMap(hyperbolic,xi,mu,1,&er,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PEPResidualNorm(A,er,0.0,vr,vi,w,&error);CHKERRQ(ierr);
    error /= SlepcAbsEigenvalue(er,0.0);
#if defined(PETSC_USE_COMPLEX)
    re = PetscRealPart(er);
    im = PetscImaginaryPart(ei);
#else
    re = er;
    im = ei;
#endif
    if (im!=0.0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  % 9f%+9fi      %12g\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"    % 12f           %12g\n",(double)re,(double)error);CHKERRQ(ierr);
    }
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s",sep);CHKERRQ(ierr);
  for (i=0;i<4;i++) {ierr = VecDestroy(w+i);CHKERRQ(ierr);}
  ierr = VecDestroy(&vi);CHKERRQ(ierr);
  ierr = VecDestroy(&vr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   testset:
      requires: !single
      args: -pep_nev 3 -nonhyperbolic -pep_target 2
      output_file: output/ex40_1.out
      filter: grep -v "Definite" | sed -e "s/iterations [0-9]\([0-9]*\)/iterations xx/g" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      test:
         suffix: 1
         requires: !complex
      test:
         suffix: 1_complex
         requires: complex !mumps
      test:
         suffix: 1_transform
         requires: !complex
         args: -transform
      test:
         suffix: 1_transform_complex
         requires: complex !mumps
         args: -transform

TEST*/
