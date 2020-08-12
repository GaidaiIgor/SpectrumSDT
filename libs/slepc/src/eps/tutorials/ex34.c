/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   This is a nonlinear eigenvalue problem. When p=2, it is reduced to a linear Laplace eigenvalue
   problem.

   -\nabla\cdot(|\nabla u|^{p-2} \nabla u) = k |u|^{p-2} u in (0,1)x(0,1),

   u = 0 on the entire boundary.

   The code is implemented based on DMPlex using Q1 FEM on a quadrilateral mesh. In this code, we consider p=3.

   Contributed  by Fande Kong fdkong.jd@gmail.com
*/

static char help[] = "Nonlinear inverse iteration for A(x)*x=lambda*B(x)*x.\n\n";


#include <slepceps.h>
#include <petscdmplex.h>
#include <petscds.h>

PetscErrorCode CreateSquareMesh(MPI_Comm,DM*);
PetscErrorCode SetupDiscretization(DM);
PetscErrorCode FormJacobianA(SNES,Vec,Mat,Mat,void*);
PetscErrorCode FormFunctionA(SNES,Vec,Vec,void*);
PetscErrorCode MatMult_A(Mat A,Vec x,Vec y);
PetscErrorCode FormJacobianB(SNES,Vec,Mat,Mat,void*);
PetscErrorCode FormFunctionB(SNES,Vec,Vec,void*);
PetscErrorCode MatMult_B(Mat A,Vec x,Vec y);
PetscErrorCode FormFunctionAB(SNES,Vec,Vec,Vec,void*);
PetscErrorCode BoundaryGlobalIndex(DM,const char*,IS*);

typedef struct {
  IS    bdis; /* global indices for boundary DoFs */
  SNES  snes;
} AppCtx;

int main(int argc,char **argv)
{
  DM             dm;
  MPI_Comm       comm;
  AppCtx         user;
  EPS            eps;  /* eigenproblem solver context */
  ST             st;
  EPSType        type;
  Mat            A,B,P;
  Vec            v0;
  PetscContainer container;
  PetscInt       nev,nconv,m,n,M,N;
  PetscBool      nonlin,flg=PETSC_FALSE,update;
  SNES           snes;
  PetscReal      tol,relerr;
  PetscBool      use_shell_matrix=PETSC_FALSE,test_init_sol=PETSC_FALSE;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  comm = PETSC_COMM_WORLD;
  /* Create a quadrilateral mesh on domain (0,1)x(0,1) */
  ierr = CreateSquareMesh(comm,&dm);CHKERRQ(ierr);
  /* Setup basis function */
  ierr = SetupDiscretization(dm);CHKERRQ(ierr);
  ierr = BoundaryGlobalIndex(dm,"marker",&user.bdis);CHKERRQ(ierr);
  /* Check if we are going to use shell matrices */
  ierr = PetscOptionsGetBool(NULL,NULL,"-use_shell_matrix",&use_shell_matrix,NULL);CHKERRQ(ierr);
  if (use_shell_matrix) {
    ierr = DMCreateMatrix(dm,&P);CHKERRQ(ierr);
    ierr = MatGetLocalSize(P,&m,&n);CHKERRQ(ierr);
    ierr = MatGetSize(P,&M,&N);CHKERRQ(ierr);
    ierr = MatCreateShell(comm,m,n,M,N,&user,&A);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatMult_A);CHKERRQ(ierr);
    ierr = MatCreateShell(comm,m,n,M,N,&user,&B);CHKERRQ(ierr);
    ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))MatMult_B);CHKERRQ(ierr);
  } else {
    ierr = DMCreateMatrix(dm,&A);CHKERRQ(ierr);
    ierr = MatDuplicate(A,MAT_COPY_VALUES,&B);CHKERRQ(ierr);
  }

  /*
     Compose callback functions and context that will be needed by the solver
  */
  ierr = PetscObjectComposeFunction((PetscObject)A,"formFunction",FormFunctionA);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-form_function_ab",&flg,NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscObjectComposeFunction((PetscObject)A,"formFunctionAB",FormFunctionAB);CHKERRQ(ierr);
  }
  ierr = PetscObjectComposeFunction((PetscObject)A,"formJacobian",FormJacobianA);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)B,"formFunction",FormFunctionB);CHKERRQ(ierr);
  ierr = PetscContainerCreate(comm,&container);CHKERRQ(ierr);
  ierr = PetscContainerSetPointer(container,&user);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)A,"formFunctionCtx",(PetscObject)container);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)A,"formJacobianCtx",(PetscObject)container);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)B,"formFunctionCtx",(PetscObject)container);CHKERRQ(ierr);
  ierr = PetscContainerDestroy(&container);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(comm,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_GNHEP);CHKERRQ(ierr);
  /*
     Use nonlinear inverse iteration
  */
  ierr = EPSSetType(eps,EPSPOWER);CHKERRQ(ierr);
  ierr = EPSPowerSetNonlinear(eps,PETSC_TRUE);CHKERRQ(ierr);
  /*
    Attach DM to SNES
  */
  ierr = EPSPowerGetSNES(eps,&snes);CHKERRQ(ierr);
  user.snes = snes;
  ierr = SNESSetDM(snes,dm);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* Set a preconditioning matrix to ST */
  if (use_shell_matrix) {
    ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
    ierr = STPrecondSetMatForPC(st,P);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-test_init_sol",&test_init_sol,NULL);CHKERRQ(ierr);
  if (nconv && test_init_sol) {
    PetscScalar   k;
    PetscReal     norm0;
    PetscInt      nits;

    ierr = MatCreateVecs(A,&v0,NULL);CHKERRQ(ierr);
    ierr = EPSGetEigenpair(eps,0,&k,NULL,v0,NULL);CHKERRQ(ierr);
    ierr = EPSSetInitialSpace(eps,1,&v0);CHKERRQ(ierr);
    ierr = VecDestroy(&v0);CHKERRQ(ierr);
    /* Norm of the previous residual */
    ierr = SNESGetFunctionNorm(snes,&norm0);CHKERRQ(ierr);
    /* Make the tolerance smaller than the last residual
       SNES will converge right away if the initial is setup correctly */
    ierr = SNESSetTolerances(snes,norm0*1.2,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    ierr = EPSSolve(eps);CHKERRQ(ierr);
    /* Number of Newton iterations supposes to be zero */
    ierr = SNESGetIterationNumber(snes,&nits);CHKERRQ(ierr);
    if (nits) {
      ierr = PetscPrintf(comm," Number of Newtoniterations %D should be zero \n",nits);CHKERRQ(ierr);
    }
  }

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,NULL);CHKERRQ(ierr);
  ierr = EPSPowerGetNonlinear(eps,&nonlin);CHKERRQ(ierr);
  ierr = EPSPowerGetUpdate(eps,&update);CHKERRQ(ierr);
  ierr = PetscPrintf(comm," Solution method: %s%s\n\n",type,nonlin?(update?" (nonlinear with monolithic update)":" (nonlinear)"):"");CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(comm," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

  /* print eigenvalue and error */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  if (nconv>0) {
    PetscScalar   k;
    PetscReal     na,nb;
    Vec           a,b,eigen;
    ierr = DMCreateGlobalVector(dm,&a);CHKERRQ(ierr);
    ierr = VecDuplicate(a,&b);CHKERRQ(ierr);
    ierr = VecDuplicate(a,&eigen);CHKERRQ(ierr);
    ierr = EPSGetEigenpair(eps,0,&k,NULL,eigen,NULL);CHKERRQ(ierr);
    ierr = FormFunctionA(snes,eigen,a,&user);CHKERRQ(ierr);
    ierr = FormFunctionB(snes,eigen,b,&user);CHKERRQ(ierr);
    ierr = VecAXPY(a,-k,b);CHKERRQ(ierr);
    ierr = VecNorm(a,NORM_2,&na);CHKERRQ(ierr);
    ierr = VecNorm(b,NORM_2,&nb);CHKERRQ(ierr);
    relerr = na/(nb*PetscAbsScalar(k));
    if (relerr<10*tol) {
      ierr = PetscPrintf(comm,"k: %g, relative error below tol\n",(double)PetscRealPart(k));CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(comm,"k: %g, relative error: %g\n",(double)PetscRealPart(k),(double)relerr);CHKERRQ(ierr);
    }
    ierr = VecDestroy(&a);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    ierr = VecDestroy(&eigen);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(comm,"Solver did not converge\n");CHKERRQ(ierr);
  }

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  if (use_shell_matrix) {
    ierr = MatDestroy(&P);CHKERRQ(ierr);
  }
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = ISDestroy(&user.bdis);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/* <|u|u, v> */
static void f0_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                 const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                 PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscScalar cof = PetscAbsScalar(u[0]);

  f0[0] = cof*u[0];
}

/* <|\nabla u| \nabla u, \nabla v> */
static void f1_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                 const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                 PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt    d;
  PetscScalar cof = 0;
  for (d = 0; d < dim; ++d)  cof += u_x[d]*u_x[d];

  cof = PetscSqrtScalar(cof);

  for (d = 0; d < dim; ++d) f1[d] = u_x[d]*cof;
}

/* approximate  Jacobian for   <|u|u, v> */
static void g0_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = 1.0*PetscAbsScalar(u[0]);
}

/* approximate  Jacobian for   <|\nabla u| \nabla u, \nabla v> */
static void g3_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d;

  for (d = 0; d < dim; ++d) g3[d*dim+d] = 1.0;
}

PetscErrorCode SetupDiscretization(DM dm)
{
  PetscFE        fe;
  MPI_Comm       comm;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* Create finite element */
  ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
  ierr = PetscFECreateDefault(comm,2,1,PETSC_FALSE,NULL,-1,&fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)fe,"u");CHKERRQ(ierr);
  ierr = DMSetField(dm,0,NULL,(PetscObject)fe);CHKERRQ(ierr);
  ierr = DMCreateDS(dm);CHKERRQ(ierr);
  ierr = PetscFEDestroy(&fe);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode CreateSquareMesh(MPI_Comm comm,DM *dm)
{
  PetscInt       cells[] = {8,8};
  PetscInt       dim = 2;
  DM             pdm;
  PetscMPIInt    size;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexCreateBoxMesh(comm,dim,PETSC_FALSE,cells,NULL,NULL,NULL,PETSC_TRUE,dm);CHKERRQ(ierr);
  ierr = DMSetFromOptions(*dm);CHKERRQ(ierr);
  ierr = DMSetUp(*dm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  if (size > 1) {
    ierr = DMPlexDistribute(*dm,0,NULL,&pdm);CHKERRQ(ierr);
    ierr = DMDestroy(dm);CHKERRQ(ierr);
    *dm = pdm;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode BoundaryGlobalIndex(DM dm,const char labelname[],IS *bdis)
{
  IS             bdpoints;
  PetscInt       nindices,*indices,numDof,offset,npoints,i,j;
  const PetscInt *bdpoints_indices;
  DMLabel        bdmarker;
  PetscSection   gsection;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetGlobalSection(dm,&gsection);CHKERRQ(ierr);
  ierr = DMGetLabel(dm,labelname,&bdmarker);CHKERRQ(ierr);
  ierr = DMLabelGetStratumIS(bdmarker,1,&bdpoints);CHKERRQ(ierr);
  ierr = ISGetLocalSize(bdpoints,&npoints);CHKERRQ(ierr);
  ierr = ISGetIndices(bdpoints,&bdpoints_indices);CHKERRQ(ierr);
  nindices = 0;
  for (i=0;i<npoints;i++) {
    ierr = PetscSectionGetDof(gsection,bdpoints_indices[i],&numDof);CHKERRQ(ierr);
    if (numDof<=0) continue;
    nindices += numDof;
  }
  ierr = PetscCalloc1(nindices,&indices);CHKERRQ(ierr);
  nindices = 0;
  for (i=0;i<npoints;i++) {
    ierr = PetscSectionGetDof(gsection,bdpoints_indices[i],&numDof);CHKERRQ(ierr);
    if (numDof<=0) continue;
    ierr = PetscSectionGetOffset(gsection,bdpoints_indices[i],&offset);CHKERRQ(ierr);
    for (j=0;j<numDof;j++) indices[nindices++] = offset+j;
  }
  ierr = ISRestoreIndices(bdpoints,&bdpoints_indices);CHKERRQ(ierr);
  ierr = ISDestroy(&bdpoints);CHKERRQ(ierr);
  ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm),nindices,indices,PETSC_OWN_POINTER,bdis);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode FormJacobian(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
  DM             dm;
  Vec            Xloc;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&Xloc);CHKERRQ(ierr);
  ierr = VecZeroEntries(Xloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,X,INSERT_VALUES,Xloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,X,INSERT_VALUES,Xloc);CHKERRQ(ierr);
  CHKMEMQ;
  ierr = DMPlexSNESComputeJacobianFEM(dm,Xloc,A,B,ctx);CHKERRQ(ierr);
  if (A!=B) {
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  CHKMEMQ;
  ierr = DMRestoreLocalVector(dm,&Xloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FormJacobianA(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
  PetscErrorCode ierr;
  DM             dm;
  PetscDS        prob;
  AppCtx         *userctx = (AppCtx *)ctx;

  PetscFunctionBegin;
  ierr = MatSetOption(B,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr = DMGetDS(dm,&prob);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob,0,0,NULL,NULL,NULL,g3_uu);CHKERRQ(ierr);
  ierr = FormJacobian(snes,X,A,B,ctx);CHKERRQ(ierr);
  ierr = MatZeroRowsIS(B,userctx->bdis,1.0,NULL,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FormJacobianB(SNES snes,Vec X,Mat A,Mat B,void *ctx)
{
  PetscErrorCode ierr;
  DM             dm;
  PetscDS        prob;
  AppCtx         *userctx = (AppCtx *)ctx;

  PetscFunctionBegin;
  ierr = MatSetOption(B,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr = DMGetDS(dm,&prob);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob,0,0,g0_uu,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = FormJacobian(snes,X,A,B,ctx);CHKERRQ(ierr);
  ierr = MatZeroRowsIS(B,userctx->bdis,0.0,NULL,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FormFunctionAB(SNES snes,Vec x,Vec Ax,Vec Bx,void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /*
   * In real applications, users should have a generic formFunctionAB which
   * forms Ax and Bx simultaneously for an more efficient calculation.
   * In this example, we just call FormFunctionA+FormFunctionB to mimic how
   * to use FormFunctionAB
   */
  ierr = FormFunctionA(snes,x,Ax,ctx);CHKERRQ(ierr);
  ierr = FormFunctionB(snes,x,Bx,ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


static PetscErrorCode FormFunction(SNES snes,Vec X,Vec F,void *ctx)
{
  DM             dm;
  Vec            Xloc,Floc;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&Xloc);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&Floc);CHKERRQ(ierr);
  ierr = VecZeroEntries(Xloc);CHKERRQ(ierr);
  ierr = VecZeroEntries(Floc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,X,INSERT_VALUES,Xloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,X,INSERT_VALUES,Xloc);CHKERRQ(ierr);
  CHKMEMQ;
  ierr = DMPlexSNESComputeResidualFEM(dm,Xloc,Floc,ctx);CHKERRQ(ierr);
  CHKMEMQ;
  ierr = VecZeroEntries(F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,Floc,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,Floc,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Xloc);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Floc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FormFunctionA(SNES snes,Vec X,Vec F,void *ctx)
{
  PetscErrorCode ierr;
  DM             dm;
  PetscDS        prob;
  PetscInt       nindices,iStart,iEnd,i;
  AppCtx         *userctx = (AppCtx *)ctx;
  PetscScalar    *array,value;
  const PetscInt *indices;
  PetscInt       vecstate;

  PetscFunctionBegin;
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr = DMGetDS(dm,&prob);CHKERRQ(ierr);
  /* hook functions */
  ierr = PetscDSSetResidual(prob,0,NULL,f1_u);CHKERRQ(ierr);
  ierr = FormFunction(snes,X,F,ctx);CHKERRQ(ierr);
  /* Boundary condition */
  ierr = VecLockGet(X,&vecstate);CHKERRQ(ierr);
  if (vecstate>0) {
    ierr = VecLockReadPop(X);CHKERRQ(ierr);
  }
  ierr = VecGetOwnershipRange(X,&iStart,&iEnd);CHKERRQ(ierr);
  ierr = VecGetArray(X,&array);CHKERRQ(ierr);
  ierr = ISGetLocalSize(userctx->bdis,&nindices);CHKERRQ(ierr);
  ierr = ISGetIndices(userctx->bdis,&indices);CHKERRQ(ierr);
  for (i=0;i<nindices;i++) {
    value = array[indices[i]-iStart] - 0.0;
    ierr = VecSetValue(F,indices[i],value,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = ISRestoreIndices(userctx->bdis,&indices);CHKERRQ(ierr);
  ierr = VecRestoreArray(X,&array);CHKERRQ(ierr);
  if (vecstate>0) {
    ierr = VecLockReadPush(X);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(F);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMult_A(Mat A,Vec x,Vec y)
{
  PetscErrorCode ierr;
  AppCtx         *userctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,&userctx);CHKERRQ(ierr);
  ierr = FormFunctionA(userctx->snes,x,y,userctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FormFunctionB(SNES snes,Vec X,Vec F,void *ctx)
{
  PetscErrorCode ierr;
  DM             dm;
  PetscDS        prob;
  PetscInt       nindices,iStart,iEnd,i;
  AppCtx         *userctx = (AppCtx *)ctx;
  PetscScalar    value;
  const PetscInt *indices;

  PetscFunctionBegin;
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  ierr = DMGetDS(dm,&prob);CHKERRQ(ierr);
  /* hook functions */
  ierr = PetscDSSetResidual(prob,0,f0_u,NULL);CHKERRQ(ierr);
  ierr = FormFunction(snes,X,F,ctx);CHKERRQ(ierr);
  /* Boundary condition */
  ierr = VecGetOwnershipRange(F,&iStart,&iEnd);CHKERRQ(ierr);
  ierr = ISGetLocalSize(userctx->bdis,&nindices);CHKERRQ(ierr);
  ierr = ISGetIndices(userctx->bdis,&indices);CHKERRQ(ierr);
  for (i=0;i<nindices;i++) {
    value = 0.0;
    ierr = VecSetValue(F,indices[i],value,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = ISRestoreIndices(userctx->bdis,&indices);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatMult_B(Mat B,Vec x,Vec y)
{
  PetscErrorCode ierr;
  AppCtx         *userctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(B,&userctx);CHKERRQ(ierr);
  ierr = FormFunctionB(userctx->snes,x,y,userctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   testset:
      requires: double
      args: -petscspace_degree 1 -petscspace_poly_tensor
      output_file: output/ex34_1.out
      test:
         suffix: 1
      test:
         suffix: 2
         args: -eps_power_update -form_function_ab {{0 1}}
         filter: sed -e "s/ with monolithic update//"
      test:
         suffix: 3
         args: -use_shell_matrix -eps_power_snes_mf_operator 1
      test:
         suffix: 4
         args: -use_shell_matrix -eps_power_update -init_eps_power_snes_mf_operator 1 -eps_power_snes_mf_operator 1 -form_function_ab {{0 1}}
         filter: sed -e "s/ with monolithic update//"
      test:
         suffix: 5
         args: -use_shell_matrix -eps_power_update -init_eps_power_snes_mf_operator 1 -eps_power_snes_mf_operator 1 -form_function_ab {{0 1}} -test_init_sol 1
         filter: sed -e "s/ with monolithic update//"

TEST*/
