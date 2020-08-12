
static char help[] ="Solves a simple data assimilation problem with one dimensional advection diffusion equation using TSAdjoint\n\n";

/*

    Not yet tested in parallel

*/
/*
   Concepts: TS^time-dependent linear problems
   Concepts: TS^heat equation
   Concepts: TS^diffusion equation
   Concepts: adjoints
   Processors: n
*/

/* ------------------------------------------------------------------------

   This program uses the one-dimensional advection-diffusion equation),
       u_t = mu*u_xx - a u_x,
   on the domain 0 <= x <= 1, with periodic boundary conditions

   to demonstrate solving a data assimilation problem of finding the initial conditions
   to produce a given solution at a fixed time.

   The operators are discretized with the spectral element method

  ------------------------------------------------------------------------- */

/*
   Include "petscts.h" so that we can use TS solvers.  Note that this file
   automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h  - vectors
     petscmat.h  - matrices
     petscis.h     - index sets            petscksp.h  - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h   - preconditioners
     petscksp.h   - linear solvers        petscsnes.h - nonlinear solvers
*/

#include <petsctao.h>
#include <petscts.h>
#include <petscdt.h>
#include <petscdraw.h>
#include <petscdmda.h>

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines.
*/

typedef struct {
  PetscInt  n;                /* number of nodes */
  PetscReal *nodes;           /* GLL nodes */
  PetscReal *weights;         /* GLL weights */
} PetscGLL;

typedef struct {
  PetscInt    N;             /* grid points per elements*/
  PetscInt    E;              /* number of elements */
  PetscReal   tol_L2,tol_max; /* error norms */
  PetscInt    steps;          /* number of timesteps */
  PetscReal   Tend;           /* endtime */
  PetscReal   mu;             /* viscosity */
  PetscReal   a;              /* advection speed */
  PetscReal   L;              /* total length of domain */
  PetscReal   Le;
  PetscReal   Tadj;
} PetscParam;

typedef struct {
  Vec         reference;               /* desired end state */
  Vec         grid;              /* total grid */
  Vec         grad;
  Vec         ic;
  Vec         curr_sol;
  Vec         joe;
  Vec         true_solution;     /* actual initial conditions for the final solution */
} PetscData;

typedef struct {
  Vec         grid;              /* total grid */
  Vec         mass;              /* mass matrix for total integration */
  Mat         stiff;             /* stifness matrix */
  Mat         advec;
  Mat         keptstiff;
  PetscGLL    gll;
} PetscSEMOperators;

typedef struct {
  DM                da;                /* distributed array data structure */
  PetscSEMOperators SEMop;
  PetscParam        param;
  PetscData         dat;
  TS                ts;
  PetscReal         initial_dt;
  PetscReal         *solutioncoefficients;
  PetscInt          ncoeff;
} AppCtx;

/*
   User-defined routines
*/
extern PetscErrorCode FormFunctionGradient(Tao,Vec,PetscReal*,Vec,void*);
extern PetscErrorCode RHSLaplacian(TS,PetscReal,Vec,Mat,Mat,void*);
extern PetscErrorCode RHSAdvection(TS,PetscReal,Vec,Mat,Mat,void*);
extern PetscErrorCode InitialConditions(Vec,AppCtx*);
extern PetscErrorCode ComputeReference(TS,PetscReal,Vec,AppCtx*);
extern PetscErrorCode MonitorError(Tao,void*);
extern PetscErrorCode MonitorDestroy(void**);
extern PetscErrorCode ComputeSolutionCoefficients(AppCtx*);
extern PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);
extern PetscErrorCode RHSJacobian(TS,PetscReal,Vec,Mat,Mat,void*);

int main(int argc,char **argv)
{
  AppCtx         appctx;                 /* user-defined application context */
  Tao            tao;
  Vec            u;                      /* approximate solution vector */
  PetscErrorCode ierr;
  PetscInt       i, xs, xm, ind, j, lenglob;
  PetscReal      x, *wrk_ptr1, *wrk_ptr2;
  MatNullSpace   nsp;

   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program and set problem parameters
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscFunctionBegin;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /*initialize parameters */
  appctx.param.N    = 10;  /* order of the spectral element */
  appctx.param.E    = 8;  /* number of elements */
  appctx.param.L    = 1.0;  /* length of the domain */
  appctx.param.mu   = 0.00001; /* diffusion coefficient */
  appctx.param.a    = 0.0;     /* advection speed */
  appctx.initial_dt = 1e-4;
  appctx.param.steps = PETSC_MAX_INT;
  appctx.param.Tend  = 0.01;
  appctx.ncoeff      = 2;

  ierr = PetscOptionsGetInt(NULL,NULL,"-N",&appctx.param.N,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-E",&appctx.param.E,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ncoeff",&appctx.ncoeff,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-Tend",&appctx.param.Tend,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-mu",&appctx.param.mu,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-a",&appctx.param.a,NULL);CHKERRQ(ierr);
  appctx.param.Le = appctx.param.L/appctx.param.E;


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create GLL data structures
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscMalloc2(appctx.param.N,&appctx.SEMop.gll.nodes,appctx.param.N,&appctx.SEMop.gll.weights);CHKERRQ(ierr);
  ierr = PetscDTGaussLobattoLegendreQuadrature(appctx.param.N,PETSCGAUSSLOBATTOLEGENDRE_VIA_LINEAR_ALGEBRA,appctx.SEMop.gll.nodes,appctx.SEMop.gll.weights);CHKERRQ(ierr);
  appctx.SEMop.gll.n = appctx.param.N;
  lenglob  = appctx.param.E*(appctx.param.N-1);

  /*
     Create distributed array (DMDA) to manage parallel grid and vectors
     and to set up the ghost point communication pattern.  There are E*(Nl-1)+1
     total grid values spread equally among all the processors, except first and last
  */

  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_PERIODIC,lenglob,1,1,NULL,&appctx.da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(appctx.da);CHKERRQ(ierr);
  ierr = DMSetUp(appctx.da);CHKERRQ(ierr);

  /*
     Extract global and local vectors from DMDA; we use these to store the
     approximate solution.  Then duplicate these for remaining vectors that
     have the same types.
  */

  ierr = DMCreateGlobalVector(appctx.da,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&appctx.dat.ic);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&appctx.dat.true_solution);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&appctx.dat.reference);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&appctx.SEMop.grid);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&appctx.SEMop.mass);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&appctx.dat.curr_sol);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&appctx.dat.joe);CHKERRQ(ierr);

  ierr = DMDAGetCorners(appctx.da,&xs,NULL,NULL,&xm,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(appctx.da,appctx.SEMop.grid,&wrk_ptr1);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(appctx.da,appctx.SEMop.mass,&wrk_ptr2);CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */

    xs=xs/(appctx.param.N-1);
    xm=xm/(appctx.param.N-1);

  /*
     Build total grid and mass over entire mesh (multi-elemental)
  */

  for (i=xs; i<xs+xm; i++) {
    for (j=0; j<appctx.param.N-1; j++) {
      x = (appctx.param.Le/2.0)*(appctx.SEMop.gll.nodes[j]+1.0)+appctx.param.Le*i;
      ind=i*(appctx.param.N-1)+j;
      wrk_ptr1[ind]=x;
      wrk_ptr2[ind]=.5*appctx.param.Le*appctx.SEMop.gll.weights[j];
      if (j==0) wrk_ptr2[ind]+=.5*appctx.param.Le*appctx.SEMop.gll.weights[j];
    }
  }
  ierr = DMDAVecRestoreArray(appctx.da,appctx.SEMop.grid,&wrk_ptr1);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(appctx.da,appctx.SEMop.mass,&wrk_ptr2);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create matrix data structure; set matrix evaluation routine.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMSetMatrixPreallocateOnly(appctx.da, PETSC_TRUE);CHKERRQ(ierr);
  ierr = DMCreateMatrix(appctx.da,&appctx.SEMop.stiff);CHKERRQ(ierr);
  ierr = DMCreateMatrix(appctx.da,&appctx.SEMop.advec);CHKERRQ(ierr);

  /*
   For linear problems with a time-dependent f(u,t) in the equation
   u_t = f(u,t), the user provides the discretized right-hand-side
   as a time-dependent matrix.
   */
  ierr = RHSLaplacian(appctx.ts,0.0,u,appctx.SEMop.stiff,appctx.SEMop.stiff,&appctx);CHKERRQ(ierr);
  ierr = RHSAdvection(appctx.ts,0.0,u,appctx.SEMop.advec,appctx.SEMop.advec,&appctx);CHKERRQ(ierr);
  ierr = MatAXPY(appctx.SEMop.stiff,-1.0,appctx.SEMop.advec,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatDuplicate(appctx.SEMop.stiff,MAT_COPY_VALUES,&appctx.SEMop.keptstiff);CHKERRQ(ierr);

  /* attach the null space to the matrix, this probably is not needed but does no harm */
  ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&nsp);CHKERRQ(ierr);
  ierr = MatSetNullSpace(appctx.SEMop.stiff,nsp);CHKERRQ(ierr);
  ierr = MatNullSpaceTest(nsp,appctx.SEMop.stiff,NULL);CHKERRQ(ierr);
  ierr = MatNullSpaceDestroy(&nsp);CHKERRQ(ierr);

  /* Create the TS solver that solves the ODE and its adjoint; set its options */
  ierr = TSCreate(PETSC_COMM_WORLD,&appctx.ts);CHKERRQ(ierr);
  ierr = TSSetSolutionFunction(appctx.ts,(PetscErrorCode (*)(TS,PetscReal,Vec, void *))ComputeReference,&appctx);CHKERRQ(ierr);
  ierr = TSSetProblemType(appctx.ts,TS_LINEAR);CHKERRQ(ierr);
  ierr = TSSetType(appctx.ts,TSRK);CHKERRQ(ierr);
  ierr = TSSetDM(appctx.ts,appctx.da);CHKERRQ(ierr);
  ierr = TSSetTime(appctx.ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(appctx.ts,appctx.initial_dt);CHKERRQ(ierr);
  ierr = TSSetMaxSteps(appctx.ts,appctx.param.steps);CHKERRQ(ierr);
  ierr = TSSetMaxTime(appctx.ts,appctx.param.Tend);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(appctx.ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetTolerances(appctx.ts,1e-7,NULL,1e-7,NULL);CHKERRQ(ierr);
  ierr = TSSetFromOptions(appctx.ts);CHKERRQ(ierr);
  /* Need to save initial timestep user may have set with -ts_dt so it can be reset for each new TSSolve() */
  ierr = TSGetTimeStep(appctx.ts,&appctx.initial_dt);CHKERRQ(ierr);
  ierr = TSSetRHSFunction(appctx.ts,NULL,TSComputeRHSFunctionLinear,&appctx);CHKERRQ(ierr);
  ierr = TSSetRHSJacobian(appctx.ts,appctx.SEMop.stiff,appctx.SEMop.stiff,TSComputeRHSJacobianConstant,&appctx);CHKERRQ(ierr);
  /*  ierr = TSSetRHSFunction(appctx.ts,NULL,RHSFunction,&appctx);CHKERRQ(ierr);
      ierr = TSSetRHSJacobian(appctx.ts,appctx.SEMop.stiff,appctx.SEMop.stiff,RHSJacobian,&appctx);CHKERRQ(ierr); */

  /* Set random initial conditions as initial guess, compute analytic reference solution and analytic (true) initial conditions */
  ierr = ComputeSolutionCoefficients(&appctx);CHKERRQ(ierr);
  ierr = InitialConditions(appctx.dat.ic,&appctx);CHKERRQ(ierr);
  ierr = ComputeReference(appctx.ts,appctx.param.Tend,appctx.dat.reference,&appctx);CHKERRQ(ierr);
  ierr = ComputeReference(appctx.ts,0.0,appctx.dat.true_solution,&appctx);CHKERRQ(ierr);

  /* Set up to save trajectory before TSSetFromOptions() so that TSTrajectory options can be captured */
  ierr = TSSetSaveTrajectory(appctx.ts);CHKERRQ(ierr);
  ierr = TSSetFromOptions(appctx.ts);CHKERRQ(ierr);

  /* Create TAO solver and set desired solution method  */
  ierr = TaoCreate(PETSC_COMM_WORLD,&tao);CHKERRQ(ierr);
  ierr = TaoSetMonitor(tao,MonitorError,&appctx,MonitorDestroy);CHKERRQ(ierr);
  ierr = TaoSetType(tao,TAOBQNLS);CHKERRQ(ierr);
  ierr = TaoSetInitialVector(tao,appctx.dat.ic);CHKERRQ(ierr);
  /* Set routine for function and gradient evaluation  */
  ierr = TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,(void *)&appctx);CHKERRQ(ierr);
  /* Check for any TAO command line options  */
  ierr = TaoSetTolerances(tao,1e-8,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = TaoSetFromOptions(tao);CHKERRQ(ierr);
  ierr = TaoSolve(tao);CHKERRQ(ierr);

  ierr = TaoDestroy(&tao);CHKERRQ(ierr);
  ierr = PetscFree(appctx.solutioncoefficients);CHKERRQ(ierr);
  ierr = MatDestroy(&appctx.SEMop.advec);CHKERRQ(ierr);
  ierr = MatDestroy(&appctx.SEMop.stiff);CHKERRQ(ierr);
  ierr = MatDestroy(&appctx.SEMop.keptstiff);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&appctx.dat.ic);CHKERRQ(ierr);
  ierr = VecDestroy(&appctx.dat.joe);CHKERRQ(ierr);
  ierr = VecDestroy(&appctx.dat.true_solution);CHKERRQ(ierr);
  ierr = VecDestroy(&appctx.dat.reference);CHKERRQ(ierr);
  ierr = VecDestroy(&appctx.SEMop.grid);CHKERRQ(ierr);
  ierr = VecDestroy(&appctx.SEMop.mass);CHKERRQ(ierr);
  ierr = VecDestroy(&appctx.dat.curr_sol);CHKERRQ(ierr);
  ierr = PetscFree2(appctx.SEMop.gll.nodes,appctx.SEMop.gll.weights);CHKERRQ(ierr);
  ierr = DMDestroy(&appctx.da);CHKERRQ(ierr);
  ierr = TSDestroy(&appctx.ts);CHKERRQ(ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
  */
    ierr = PetscFinalize();
    return ierr;
}

/*
    Computes the coefficients for the analytic solution to the PDE
*/
PetscErrorCode ComputeSolutionCoefficients(AppCtx *appctx)
{
  PetscErrorCode    ierr;
  PetscRandom       rand;
  PetscInt          i;

  PetscFunctionBegin;
  ierr = PetscMalloc1(appctx->ncoeff,&appctx->solutioncoefficients);CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rand);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rand,.9,1.0);CHKERRQ(ierr);
  for (i=0; i<appctx->ncoeff; i++) {
    ierr = PetscRandomGetValue(rand,&appctx->solutioncoefficients[i]);CHKERRQ(ierr);
  }
  ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*
   InitialConditions - Computes the (random) initial conditions for the Tao optimization solve (these are also initial conditions for the first TSSolve()

   Input Parameter:
   u - uninitialized solution vector (global)
   appctx - user-defined application context

   Output Parameter:
   u - vector with solution at initial time (global)
*/
PetscErrorCode InitialConditions(Vec u,AppCtx *appctx)
{
  PetscScalar       *s;
  const PetscScalar *xg;
  PetscErrorCode    ierr;
  PetscInt          i,j,lenglob;
  PetscReal         sum,val;
  PetscRandom       rand;

  PetscFunctionBegin;
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rand);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(rand,.9,1.0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(appctx->da,u,&s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(appctx->da,appctx->SEMop.grid,(void*)&xg);CHKERRQ(ierr);
  lenglob  = appctx->param.E*(appctx->param.N-1);
  for (i=0; i<lenglob; i++) {
    s[i]= 0;
    for (j=0; j<appctx->ncoeff; j++) {
      ierr =  PetscRandomGetValue(rand,&val);CHKERRQ(ierr);
      s[i] += val*PetscSinScalar(2*(j+1)*PETSC_PI*xg[i]);
    }
  }
  ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(appctx->da,u,&s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(appctx->da,appctx->SEMop.grid,(void*)&xg);CHKERRQ(ierr);
  /* make sure initial conditions do not contain the constant functions, since with periodic boundary conditions the constant functions introduce a null space */
  ierr = VecSum(u,&sum);CHKERRQ(ierr);
  ierr = VecShift(u,-sum/lenglob);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
   TrueSolution() computes the true solution for the Tao optimization solve which means they are the initial conditions for the objective function.

             InitialConditions() computes the initial conditions for the begining of the Tao iterations

   Input Parameter:
   u - uninitialized solution vector (global)
   appctx - user-defined application context

   Output Parameter:
   u - vector with solution at initial time (global)
*/
PetscErrorCode TrueSolution(Vec u,AppCtx *appctx)
{
  PetscScalar       *s;
  const PetscScalar *xg;
  PetscErrorCode    ierr;
  PetscInt          i,j,lenglob;
  PetscReal         sum;

  PetscFunctionBegin;
  ierr = DMDAVecGetArray(appctx->da,u,&s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(appctx->da,appctx->SEMop.grid,(void*)&xg);CHKERRQ(ierr);
  lenglob  = appctx->param.E*(appctx->param.N-1);
  for (i=0; i<lenglob; i++) {
    s[i]= 0;
    for (j=0; j<appctx->ncoeff; j++) {
      s[i] += appctx->solutioncoefficients[j]*PetscSinScalar(2*(j+1)*PETSC_PI*xg[i]);
    }
  }
  ierr = DMDAVecRestoreArray(appctx->da,u,&s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(appctx->da,appctx->SEMop.grid,(void*)&xg);CHKERRQ(ierr);
  /* make sure initial conditions do not contain the constant functions, since with periodic boundary conditions the constant functions introduce a null space */
  ierr = VecSum(u,&sum);CHKERRQ(ierr);
  ierr = VecShift(u,-sum/lenglob);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/* --------------------------------------------------------------------- */
/*
   Sets the desired profile for the final end time

   Input Parameters:
   t - final time
   obj - vector storing the desired profile
   appctx - user-defined application context

*/
PetscErrorCode ComputeReference(TS ts,PetscReal t,Vec obj,AppCtx *appctx)
{
  PetscScalar       *s,tc;
  const PetscScalar *xg;
  PetscErrorCode    ierr;
  PetscInt          i, j,lenglob;

  PetscFunctionBegin;
  ierr = DMDAVecGetArray(appctx->da,obj,&s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(appctx->da,appctx->SEMop.grid,(void*)&xg);CHKERRQ(ierr);
  lenglob  = appctx->param.E*(appctx->param.N-1);
  for (i=0; i<lenglob; i++) {
    s[i]= 0;
    for (j=0; j<appctx->ncoeff; j++) {
      tc    = -appctx->param.mu*(j+1)*(j+1)*4.0*PETSC_PI*PETSC_PI*t;
      s[i] += appctx->solutioncoefficients[j]*PetscSinScalar(2*(j+1)*PETSC_PI*(xg[i] + appctx->param.a*t))*PetscExpReal(tc);
    }
  }
  ierr = DMDAVecRestoreArray(appctx->da,obj,&s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(appctx->da,appctx->SEMop.grid,(void*)&xg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec globalin,Vec globalout,void *ctx)
{
  PetscErrorCode ierr;
  AppCtx          *appctx = (AppCtx*)ctx;

  PetscFunctionBegin;
  ierr = MatMult(appctx->SEMop.keptstiff,globalin,globalout);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RHSJacobian(TS ts,PetscReal t,Vec globalin,Mat A, Mat B,void *ctx)
{
  PetscErrorCode ierr;
  AppCtx         *appctx = (AppCtx*)ctx;

  PetscFunctionBegin;
  ierr = MatCopy(appctx->SEMop.keptstiff,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */

/*
   RHSLaplacian -   matrix for diffusion

   Input Parameters:
   ts - the TS context
   t - current time  (ignored)
   X - current solution (ignored)
   dummy - optional user-defined context, as set by TSetRHSJacobian()

   Output Parameters:
   AA - Jacobian matrix
   BB - optionally different matrix from which the preconditioner is built
   str - flag indicating matrix structure

   Scales by the inverse of the mass matrix (perhaps that should be pulled out)

*/
PetscErrorCode RHSLaplacian(TS ts,PetscReal t,Vec X,Mat A,Mat BB,void *ctx)
{
  PetscReal      **temp;
  PetscReal      vv;
  AppCtx         *appctx = (AppCtx*)ctx;     /* user-defined application context */
  PetscErrorCode ierr;
  PetscInt       i,xs,xn,l,j;
  PetscInt       *rowsDM;

  PetscFunctionBegin;
  /*
   Creates the element stiffness matrix for the given gll
   */
  ierr = PetscGaussLobattoLegendreElementLaplacianCreate(appctx->SEMop.gll.n,appctx->SEMop.gll.nodes,appctx->SEMop.gll.weights,&temp);CHKERRQ(ierr);

  /* scale by the size of the element */
  for (i=0; i<appctx->param.N; i++) {
    vv=-appctx->param.mu*2.0/appctx->param.Le;
    for (j=0; j<appctx->param.N; j++) temp[i][j]=temp[i][j]*vv;
  }

  ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(appctx->da,&xs,NULL,NULL,&xn,NULL,NULL);CHKERRQ(ierr);

  if (appctx->param.N-1 < 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Polynomial order must be at least 2");
  xs   = xs/(appctx->param.N-1);
  xn   = xn/(appctx->param.N-1);

  ierr = PetscMalloc1(appctx->param.N,&rowsDM);CHKERRQ(ierr);
  /*
   loop over local elements
   */
  for (j=xs; j<xs+xn; j++) {
    for (l=0; l<appctx->param.N; l++) {
      rowsDM[l] = 1+(j-xs)*(appctx->param.N-1)+l;
    }
    ierr = MatSetValuesLocal(A,appctx->param.N,rowsDM,appctx->param.N,rowsDM,&temp[0][0],ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = PetscFree(rowsDM);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecReciprocal(appctx->SEMop.mass);CHKERRQ(ierr);
  ierr = MatDiagonalScale(A,appctx->SEMop.mass,0);CHKERRQ(ierr);
  ierr = VecReciprocal(appctx->SEMop.mass);CHKERRQ(ierr);

  ierr = PetscGaussLobattoLegendreElementLaplacianDestroy(appctx->SEMop.gll.n,appctx->SEMop.gll.nodes,appctx->SEMop.gll.weights,&temp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    Almost identical to Laplacian

    Note that the element matrix is NOT scaled by the size of element like the Laplacian term.
 */
PetscErrorCode RHSAdvection(TS ts,PetscReal t,Vec X,Mat A,Mat BB,void *ctx)
{
  PetscReal      **temp;
  PetscReal      vv;
  AppCtx         *appctx = (AppCtx*)ctx;     /* user-defined application context */
  PetscErrorCode ierr;
  PetscInt       i,xs,xn,l,j;
  PetscInt       *rowsDM;

  PetscFunctionBegin;
  /*
   Creates the element stiffness matrix for the given gll
   */
  ierr = PetscGaussLobattoLegendreElementAdvectionCreate(appctx->SEMop.gll.n,appctx->SEMop.gll.nodes,appctx->SEMop.gll.weights,&temp);CHKERRQ(ierr);

  /* scale by the size of the element */
  for (i=0; i<appctx->param.N; i++) {
    vv = -appctx->param.a;
    for (j=0; j<appctx->param.N; j++) temp[i][j]=temp[i][j]*vv;
  }

  ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(appctx->da,&xs,NULL,NULL,&xn,NULL,NULL);CHKERRQ(ierr);

  if (appctx->param.N-1 < 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Polynomial order must be at least 2");
  xs   = xs/(appctx->param.N-1);
  xn   = xn/(appctx->param.N-1);

  ierr = PetscMalloc1(appctx->param.N,&rowsDM);CHKERRQ(ierr);
  /*
   loop over local elements
   */
  for (j=xs; j<xs+xn; j++) {
    for (l=0; l<appctx->param.N; l++) {
      rowsDM[l] = 1+(j-xs)*(appctx->param.N-1)+l;
    }
    ierr = MatSetValuesLocal(A,appctx->param.N,rowsDM,appctx->param.N,rowsDM,&temp[0][0],ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = PetscFree(rowsDM);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecReciprocal(appctx->SEMop.mass);CHKERRQ(ierr);
  ierr = MatDiagonalScale(A,appctx->SEMop.mass,0);CHKERRQ(ierr);
  ierr = VecReciprocal(appctx->SEMop.mass);CHKERRQ(ierr);

  ierr = PetscGaussLobattoLegendreElementAdvectionDestroy(appctx->SEMop.gll.n,appctx->SEMop.gll.nodes,appctx->SEMop.gll.weights,&temp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------ */
/*
   FormFunctionGradient - Evaluates the function and corresponding gradient.

   Input Parameters:
   tao - the Tao context
   ic   - the input vector
   ctx - optional user-defined context, as set when calling TaoSetObjectiveAndGradientRoutine()

   Output Parameters:
   f   - the newly evaluated function
   G   - the newly evaluated gradient

   Notes:

          The forward equation is
              M u_t = F(U)
          which is converted to
                u_t = M^{-1} F(u)
          in the user code since TS has no direct way of providing a mass matrix. The Jacobian of this is
                 M^{-1} J
          where J is the Jacobian of F. Now the adjoint equation is
                M v_t = J^T v
          but TSAdjoint does not solve this since it can only solve the transposed system for the
          Jacobian the user provided. Hence TSAdjoint solves
                 w_t = J^T M^{-1} w  (where w = M v)
          since there is no way to indicate the mass matrix as a separate entitity to TS. Thus one
          must be careful in initializing the "adjoint equation" and using the result. This is
          why
              G = -2 M(u(T) - u_d)
          below (instead of -2(u(T) - u_d)


*/
PetscErrorCode FormFunctionGradient(Tao tao,Vec ic,PetscReal *f,Vec G,void *ctx)
{
  AppCtx           *appctx = (AppCtx*)ctx;     /* user-defined application context */
  PetscErrorCode    ierr;
  Vec               temp;

  PetscFunctionBegin;
  ierr = TSSetTime(appctx->ts,0.0);CHKERRQ(ierr);
  ierr = TSSetStepNumber(appctx->ts,0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(appctx->ts,appctx->initial_dt);CHKERRQ(ierr);
  ierr = VecCopy(ic,appctx->dat.curr_sol);CHKERRQ(ierr);

  ierr = TSSolve(appctx->ts,appctx->dat.curr_sol);CHKERRQ(ierr);
  ierr = VecCopy(appctx->dat.curr_sol,appctx->dat.joe);CHKERRQ(ierr);

  /*     Compute the difference between the current ODE solution and target ODE solution */
  ierr = VecWAXPY(G,-1.0,appctx->dat.curr_sol,appctx->dat.reference);CHKERRQ(ierr);

  /*     Compute the objective/cost function   */
  ierr = VecDuplicate(G,&temp);CHKERRQ(ierr);
  ierr = VecPointwiseMult(temp,G,G);CHKERRQ(ierr);
  ierr = VecDot(temp,appctx->SEMop.mass,f);CHKERRQ(ierr);
  ierr = VecDestroy(&temp);CHKERRQ(ierr);

  /*     Compute initial conditions for the adjoint integration. See Notes above  */
  ierr = VecScale(G, -2.0);CHKERRQ(ierr);
  ierr = VecPointwiseMult(G,G,appctx->SEMop.mass);CHKERRQ(ierr);
  ierr = TSSetCostGradients(appctx->ts,1,&G,NULL);CHKERRQ(ierr);

  ierr = TSAdjointSolve(appctx->ts);CHKERRQ(ierr);
  /* ierr = VecPointwiseDivide(G,G,appctx->SEMop.mass);CHKERRQ(ierr);*/
  PetscFunctionReturn(0);
}

PetscErrorCode MonitorError(Tao tao,void *ctx)
{
  AppCtx         *appctx = (AppCtx*)ctx;
  Vec            temp,grad;
  PetscReal      nrm;
  PetscErrorCode ierr;
  PetscInt       its;
  PetscReal      fct,gnorm;

  PetscFunctionBegin;
  ierr  = VecDuplicate(appctx->dat.ic,&temp);CHKERRQ(ierr);
  ierr  = VecWAXPY(temp,-1.0,appctx->dat.ic,appctx->dat.true_solution);CHKERRQ(ierr);
  ierr  = VecPointwiseMult(temp,temp,temp);CHKERRQ(ierr);
  ierr  = VecDot(temp,appctx->SEMop.mass,&nrm);CHKERRQ(ierr);
  nrm   = PetscSqrtReal(nrm);
  ierr  = TaoGetGradientVector(tao,&grad);CHKERRQ(ierr);
  ierr  = VecPointwiseMult(temp,temp,temp);CHKERRQ(ierr);
  ierr  = VecDot(temp,appctx->SEMop.mass,&gnorm);CHKERRQ(ierr);
  gnorm = PetscSqrtReal(gnorm);
  ierr  = VecDestroy(&temp);CHKERRQ(ierr);
  ierr  = TaoGetIterationNumber(tao,&its);CHKERRQ(ierr);
  ierr  = TaoGetObjective(tao,&fct);CHKERRQ(ierr);
  if (!its) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%% Iteration Error Objective Gradient-norm\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"history = [\n");CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%3D %g %g %g\n",its,(double)nrm,(double)fct,(double)gnorm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MonitorDestroy(void **ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"];\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*TEST

   build:
     requires: !complex

   test:
     requires: !single
     args:  -ts_adapt_dt_max 3.e-3 -E 10 -N 8 -ncoeff 5 -tao_bqnls_mat_lmvm_scale_type none

   test:
     suffix: cn
     requires: !single
     args:  -ts_type cn -ts_dt .003 -pc_type lu -E 10 -N 8 -ncoeff 5 -tao_bqnls_mat_lmvm_scale_type none

   test:
     suffix: 2
     requires: !single
     args:  -ts_adapt_dt_max 3.e-3 -E 10 -N 8 -ncoeff 5  -a .1 -tao_bqnls_mat_lmvm_scale_type none


TEST*/
