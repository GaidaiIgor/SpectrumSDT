/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "power"

   Method: Power Iteration

   Algorithm:

       This solver implements the power iteration for finding dominant
       eigenpairs. It also includes the following well-known methods:
       - Inverse Iteration: when used in combination with shift-and-invert
         spectral transformation.
       - Rayleigh Quotient Iteration (RQI): also with shift-and-invert plus
         a variable shift.

       It can also be used for nonlinear inverse iteration on the problem
       A(x)*x=lambda*B(x)*x, where A and B are not constant but depend on x.

   References:

       [1] "Single Vector Iteration Methods in SLEPc", SLEPc Technical Report
           STR-2, available at https://slepc.upv.es.
*/

#include <slepc/private/epsimpl.h>                /*I "slepceps.h" I*/
#include <slepcblaslapack.h>
/* petsc headers */
#include <petscdm.h>
#include <petscsnes.h>

static PetscErrorCode EPSPowerFormFunction_Update(SNES,Vec,Vec,void*);
static PetscErrorCode SNESLineSearchPostheckFunction(SNESLineSearch linesearch,Vec x,Vec y,Vec w,PetscBool *changed_y,PetscBool *changed_w,void *ctx);
PetscErrorCode EPSSolve_Power(EPS);
PetscErrorCode EPSSolve_TS_Power(EPS);

typedef struct {
  EPSPowerShiftType shift_type;
  PetscBool         nonlinear;
  PetscBool         update;
  SNES              snes;
  PetscErrorCode    (*formFunctionB)(SNES,Vec,Vec,void*);
  void              *formFunctionBctx;
  PetscErrorCode    (*formFunctionA)(SNES,Vec,Vec,void*);
  void              *formFunctionActx;
  PetscErrorCode    (*formFunctionAB)(SNES,Vec,Vec,Vec,void*);
  PetscInt          idx;  /* index of the first nonzero entry in the iteration vector */
  PetscMPIInt       p;    /* process id of the owner of idx */
  PetscReal         norm0; /* norm of initial vector */
} EPS_POWER;

PetscErrorCode EPSSetUp_Power(EPS eps)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  PetscBool      flg,istrivial;
  STMatMode      mode;
  Mat            A,B,P;
  Vec            res;
  PetscContainer container;
  PetscErrorCode (*formFunctionA)(SNES,Vec,Vec,void*);
  PetscErrorCode (*formJacobianA)(SNES,Vec,Mat,Mat,void*);
  void           *ctx;
  SNESLineSearch linesearch;

  PetscFunctionBegin;
  if (eps->ncv!=PETSC_DEFAULT) {
    if (eps->ncv<eps->nev) SETERRQ(PetscObjectComm((PetscObject)eps),1,"The value of ncv must be at least nev");
  } else eps->ncv = eps->nev;
  if (eps->mpd!=PETSC_DEFAULT) { ierr = PetscInfo(eps,"Warning: parameter mpd ignored\n");CHKERRQ(ierr); }
  if (eps->max_it==PETSC_DEFAULT) eps->max_it = PetscMax(1000*eps->nev,100*eps->n);
  if (!eps->which) { ierr = EPSSetWhichEigenpairs_Default(eps);CHKERRQ(ierr); }
  if (eps->which!=EPS_LARGEST_MAGNITUDE && eps->which !=EPS_TARGET_MAGNITUDE) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Wrong value of eps->which");
  if (power->shift_type != EPS_POWER_SHIFT_CONSTANT) {
    if (power->nonlinear) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Variable shifts not allowed in nonlinear problems");
    ierr = PetscObjectTypeCompareAny((PetscObject)eps->st,&flg,STSINVERT,STCAYLEY,"");CHKERRQ(ierr);
    if (!flg) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Variable shifts only allowed in shift-and-invert or Cayley ST");
    ierr = STGetMatMode(eps->st,&mode);CHKERRQ(ierr);
    if (mode == ST_MATMODE_INPLACE) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"ST matrix mode inplace does not work with variable shifts");
  }
  if (eps->extraction) { ierr = PetscInfo(eps,"Warning: extraction type ignored\n");CHKERRQ(ierr); }
  if (eps->balance!=EPS_BALANCE_NONE) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Balancing not supported in this solver");
  if (eps->arbitrary) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Arbitrary selection of eigenpairs not supported in this solver");
  ierr = RGIsTrivial(eps->rg,&istrivial);CHKERRQ(ierr);
  if (!istrivial) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"This solver does not support region filtering");
  ierr = EPSAllocateSolution(eps,0);CHKERRQ(ierr);
  ierr = EPS_SetInnerProduct(eps);CHKERRQ(ierr);

  if (power->nonlinear) {
    if (eps->nev>1) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Nonlinear inverse iteration cannot compute more than one eigenvalue");
    ierr = EPSSetWorkVecs(eps,4);CHKERRQ(ierr);

    /* Set up SNES solver */
    /* If SNES was setup, we need to reset it so that it will be consistent with the current EPS */
    if (power->snes) { ierr = SNESReset(power->snes);CHKERRQ(ierr); }
    else { ierr = EPSPowerGetSNES(eps,&power->snes);CHKERRQ(ierr); }

    /* For nonlinear solver (Newton), we should scale the initial vector back.
       The initial vector will be scaled in EPSSetUp. */
    if (eps->IS) {
      ierr = VecNorm((eps->IS)[0],NORM_2,&power->norm0);CHKERRQ(ierr);
    }

    ierr = EPSGetOperators(eps,&A,&B);CHKERRQ(ierr);

    ierr = PetscObjectQueryFunction((PetscObject)A,"formFunction",&formFunctionA);CHKERRQ(ierr);
    if (!formFunctionA) SETERRQ(PetscObjectComm((PetscObject)eps),1,"For nonlinear inverse iteration you must compose a callback function 'formFunction' in matrix A");
    ierr = PetscObjectQuery((PetscObject)A,"formFunctionCtx",(PetscObject*)&container);CHKERRQ(ierr);
    if (container) {
      ierr = PetscContainerGetPointer(container,&ctx);CHKERRQ(ierr);
    } else ctx = NULL;
    ierr = MatCreateVecs(A,&res,NULL);CHKERRQ(ierr);
    power->formFunctionA = formFunctionA;
    power->formFunctionActx = ctx;
    if (power->update) {
      ierr = SNESSetFunction(power->snes,res,EPSPowerFormFunction_Update,ctx);CHKERRQ(ierr);
      ierr = PetscObjectQueryFunction((PetscObject)A,"formFunctionAB",&power->formFunctionAB);CHKERRQ(ierr);
    }
    else { ierr = SNESSetFunction(power->snes,res,formFunctionA,ctx);CHKERRQ(ierr); }
    ierr = VecDestroy(&res);CHKERRQ(ierr);

    ierr = PetscObjectQueryFunction((PetscObject)A,"formJacobian",&formJacobianA);CHKERRQ(ierr);
    if (!formJacobianA) SETERRQ(PetscObjectComm((PetscObject)eps),1,"For nonlinear inverse iteration you must compose a callback function 'formJacobian' in matrix A");
    ierr = PetscObjectQuery((PetscObject)A,"formJacobianCtx",(PetscObject*)&container);CHKERRQ(ierr);
    if (container) {
      ierr = PetscContainerGetPointer(container,&ctx);CHKERRQ(ierr);
    } else ctx = NULL;
    /* This allows users to compute a different preconditioning matrix than A.
     * It is useful when A and B are shell matrices.
     */
    ierr = STPrecondGetMatForPC(eps->st,&P);CHKERRQ(ierr);
    /* If users do not provide a matrix, we simply use A */
    ierr = SNESSetJacobian(power->snes,A,P? P:A,formJacobianA,ctx);CHKERRQ(ierr);
    ierr = SNESSetFromOptions(power->snes);CHKERRQ(ierr);
    ierr = SNESGetLineSearch(power->snes,&linesearch);CHKERRQ(ierr);
    if (power->update) {
      ierr = SNESLineSearchSetPostCheck(linesearch,SNESLineSearchPostheckFunction,ctx);CHKERRQ(ierr);
    }
    ierr = SNESSetUp(power->snes);CHKERRQ(ierr);
    if (B) {
      ierr = PetscObjectQueryFunction((PetscObject)B,"formFunction",&power->formFunctionB);CHKERRQ(ierr);
      ierr = PetscObjectQuery((PetscObject)B,"formFunctionCtx",(PetscObject*)&container);CHKERRQ(ierr);
      if (power->formFunctionB && container) {
        ierr = PetscContainerGetPointer(container,&power->formFunctionBctx);CHKERRQ(ierr);
      } else power->formFunctionBctx = NULL;
    }
  } else {
    if (eps->twosided) {
      ierr = EPSSetWorkVecs(eps,3);CHKERRQ(ierr);
    } else {
      ierr = EPSSetWorkVecs(eps,2);CHKERRQ(ierr);
    }
    ierr = DSSetType(eps->ds,DSNHEP);CHKERRQ(ierr);
    ierr = DSAllocate(eps->ds,eps->nev);CHKERRQ(ierr);
  }
  /* dispatch solve method */
  if (eps->twosided) {
    if (power->nonlinear) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Nonlinear inverse iteration does not have two-sided variant");
    if (power->shift_type == EPS_POWER_SHIFT_WILKINSON) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_SUP,"Two-sided variant does not support Wilkinson shifts");
    eps->ops->solve = EPSSolve_TS_Power;
  } else eps->ops->solve = EPSSolve_Power;
  PetscFunctionReturn(0);
}

/*
   Find the index of the first nonzero entry of x, and the process that owns it.
*/
static PetscErrorCode FirstNonzeroIdx(Vec x,PetscInt *idx,PetscMPIInt *p)
{
  PetscErrorCode    ierr;
  PetscInt          i,first,last,N;
  PetscLayout       map;
  const PetscScalar *xx;

  PetscFunctionBegin;
  ierr = VecGetSize(x,&N);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(x,&first,&last);CHKERRQ(ierr);
  ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
  for (i=first;i<last;i++) {
    if (PetscAbsScalar(xx[i-first])>10*PETSC_MACHINE_EPSILON) break;
  }
  ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
  if (i==last) i=N;
  ierr = MPI_Allreduce(&i,idx,1,MPIU_INT,MPI_MIN,PetscObjectComm((PetscObject)x));CHKERRQ(ierr);
  if (*idx==N) SETERRQ(PetscObjectComm((PetscObject)x),1,"Zero vector found");CHKERRQ(ierr);
  ierr = VecGetLayout(x,&map);CHKERRQ(ierr);
  ierr = PetscLayoutFindOwner(map,*idx,p);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Normalize a vector x with respect to a given norm as well as the
   sign of the first nonzero entry (assumed to be idx in process p).
*/
static PetscErrorCode Normalize(Vec x,PetscReal norm,PetscInt idx,PetscMPIInt p,PetscScalar *sign)
{
  PetscErrorCode    ierr;
  PetscScalar       alpha=1.0;
  PetscInt          first,last;
  PetscReal         absv;
  const PetscScalar *xx;

  PetscFunctionBegin;
  ierr = VecGetOwnershipRange(x,&first,&last);CHKERRQ(ierr);
  if (idx>=first && idx<last) {
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    absv = PetscAbsScalar(xx[idx-first]);
    if (absv>10*PETSC_MACHINE_EPSILON) alpha = xx[idx-first]/absv;
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
  }
  ierr = MPI_Bcast(&alpha,1,MPIU_SCALAR,p,PetscObjectComm((PetscObject)x));CHKERRQ(ierr);
  if (sign) *sign = alpha;
  alpha *= norm;
  ierr = VecScale(x,1.0/alpha);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerUpdateFunctionB(EPS eps,Vec x,Vec Bx)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  Mat            B;

  PetscFunctionBegin;
  ierr = STResetMatrixState(eps->st);CHKERRQ(ierr);
  ierr = EPSGetOperators(eps,NULL,&B);CHKERRQ(ierr);
  if (B) {
    if (power->formFunctionB) {
      ierr = (*power->formFunctionB)(power->snes,x,Bx,power->formFunctionBctx);CHKERRQ(ierr);
    } else {
      ierr = MatMult(B,x,Bx);CHKERRQ(ierr);
    }
  } else {
    ierr = VecCopy(x,Bx);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerUpdateFunctionA(EPS eps,Vec x,Vec Ax)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  Mat            A;

  PetscFunctionBegin;
  ierr = STResetMatrixState(eps->st);CHKERRQ(ierr);
  ierr = EPSGetOperators(eps,&A,NULL);CHKERRQ(ierr);
  if (A) {
    if (power->formFunctionA) {
      ierr = (*power->formFunctionA)(power->snes,x,Ax,power->formFunctionActx);CHKERRQ(ierr);
    } else {
      ierr = MatMult(A,x,Ax);CHKERRQ(ierr);
    }
  } else SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_NULL,"Matrix A is required for an eigenvalue problem");
  PetscFunctionReturn(0);
}

static PetscErrorCode SNESLineSearchPostheckFunction(SNESLineSearch linesearch,Vec x,Vec y,Vec w,PetscBool *changed_y,PetscBool *changed_w,void *ctx)
{
  PetscErrorCode ierr;
  SNES           snes;
  EPS            eps;
  Vec            oldx;

  PetscFunctionBegin;
  ierr = SNESLineSearchGetSNES(linesearch,&snes);CHKERRQ(ierr);
  ierr = PetscObjectQuery((PetscObject)snes,"eps",(PetscObject *)&eps);CHKERRQ(ierr);
  if (!eps) SETERRQ(PetscObjectComm((PetscObject)snes),PETSC_ERR_ARG_NULL,"No composed EPS");
  oldx = eps->work[3];
  ierr = VecCopy(x,oldx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerFormFunction_Update(SNES snes,Vec x,Vec y,void *ctx)
{
  PetscErrorCode ierr;
  EPS            eps;
  PetscReal      bx;
  Vec            Bx;
  EPS_POWER      *power;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)snes,"eps",(PetscObject *)&eps);CHKERRQ(ierr);
  if (!eps) SETERRQ(PetscObjectComm((PetscObject)snes),PETSC_ERR_ARG_NULL,"No composed EPS");
  power = (EPS_POWER*)eps->data;
  Bx = eps->work[2];
  if (power->formFunctionAB) {
    ierr = (*power->formFunctionAB)(snes,x,y,Bx,ctx);CHKERRQ(ierr);
  } else {
    ierr = EPSPowerUpdateFunctionA(eps,x,y);CHKERRQ(ierr);
    ierr = EPSPowerUpdateFunctionB(eps,x,Bx);CHKERRQ(ierr);
  }
  ierr = VecNorm(Bx,NORM_2,&bx);CHKERRQ(ierr);
  ierr = Normalize(Bx,bx,power->idx,power->p,NULL);CHKERRQ(ierr);
  ierr = VecAXPY(y,-1.0,Bx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Use SNES to compute y = A^{-1}*B*x for the nonlinear problem
*/
static PetscErrorCode EPSPowerApply_SNES(EPS eps,Vec x,Vec y)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  Vec            Bx;

  PetscFunctionBegin;
  ierr = VecCopy(x,y);CHKERRQ(ierr);
  if (power->update) {
    ierr = SNESSolve(power->snes,NULL,y);CHKERRQ(ierr);
  } else {
    Bx = eps->work[2];
    ierr = SNESSolve(power->snes,Bx,y);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
   Use nonlinear inverse power to compute an initial guess
*/
static PetscErrorCode EPSPowerComputeInitialGuess_Update(EPS eps)
{
  EPS            powereps;
  Mat            A,B,P;
  Vec            v1,v2;
  SNES           snes;
  DM             dm,newdm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = EPSCreate(PetscObjectComm((PetscObject)eps),&powereps);CHKERRQ(ierr);
  ierr = EPSGetOperators(eps,&A,&B);CHKERRQ(ierr);
  ierr = EPSSetType(powereps,EPSPOWER);CHKERRQ(ierr);
  ierr = EPSSetOperators(powereps,A,B);CHKERRQ(ierr);
  ierr = EPSSetTolerances(powereps,1e-6,4);CHKERRQ(ierr);
  ierr = EPSSetOptionsPrefix(powereps,((PetscObject)eps)->prefix);CHKERRQ(ierr);
  ierr = EPSAppendOptionsPrefix(powereps,"init_");CHKERRQ(ierr);
  ierr = EPSSetProblemType(powereps,EPS_GNHEP);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(powereps,EPS_TARGET_MAGNITUDE);CHKERRQ(ierr);
  ierr = EPSPowerSetNonlinear(powereps,PETSC_TRUE);CHKERRQ(ierr);
  ierr = STPrecondGetMatForPC(eps->st,&P);CHKERRQ(ierr);
  /* attach dm to initial solve */
  ierr = EPSPowerGetSNES(eps,&snes);CHKERRQ(ierr);
  ierr = SNESGetDM(snes,&dm);CHKERRQ(ierr);
  /* use  dmshell to temporarily store snes context */
  ierr = DMCreate(PetscObjectComm((PetscObject)eps),&newdm);CHKERRQ(ierr);
  ierr = DMSetType(newdm,DMSHELL);CHKERRQ(ierr);
  ierr = DMSetUp(newdm);CHKERRQ(ierr);
  ierr = DMCopyDMSNES(dm,newdm);CHKERRQ(ierr);
  ierr = EPSPowerGetSNES(powereps,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,dm);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(powereps);CHKERRQ(ierr);
  if (P) {
    ierr = STPrecondSetMatForPC(powereps->st,P);CHKERRQ(ierr);
  }
  ierr = EPSSolve(powereps);CHKERRQ(ierr);
  ierr = BVGetColumn(eps->V,0,&v2);CHKERRQ(ierr);
  ierr = BVGetColumn(powereps->V,0,&v1);CHKERRQ(ierr);
  ierr = VecCopy(v1,v2);CHKERRQ(ierr);
  ierr = BVRestoreColumn(powereps->V,0,&v1);CHKERRQ(ierr);
  ierr = BVRestoreColumn(eps->V,0,&v2);CHKERRQ(ierr);
  ierr = EPSDestroy(&powereps);CHKERRQ(ierr);
  /* restore context back to the old nonlinear solver */
  ierr = DMCopyDMSNES(newdm,dm);CHKERRQ(ierr);
  ierr = DMDestroy(&newdm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_Power(EPS eps)
{
  PetscErrorCode      ierr;
  EPS_POWER           *power = (EPS_POWER*)eps->data;
  PetscInt            k,ld;
  Vec                 v,y,e,Bx;
  Mat                 A;
  KSP                 ksp;
  PetscReal           relerr,norm,norm1,rt1,rt2,cs1;
  PetscScalar         theta,rho,delta,sigma,alpha2,beta1,sn1,*T,sign;
  PetscBool           breakdown;
  KSPConvergedReason  reason;
  SNESConvergedReason snesreason;

  PetscFunctionBegin;
  e = eps->work[0];
  y = eps->work[1];
  if (power->nonlinear) Bx = eps->work[2];
  else Bx = NULL;

  if (power->shift_type != EPS_POWER_SHIFT_CONSTANT) { ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr); }
  if (power->nonlinear) {
    ierr = PetscObjectCompose((PetscObject)power->snes,"eps",(PetscObject)eps);CHKERRQ(ierr);
    /* Compute an intial guess only when users do not provide one */
    if (power->update && !eps->nini) {
      ierr = EPSPowerComputeInitialGuess_Update(eps);CHKERRQ(ierr);
    }
  } else {
    ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  }
  if (!power->update) {
    ierr = EPSGetStartVector(eps,0,NULL);CHKERRQ(ierr);
  }
  if (power->nonlinear) {
    ierr = BVGetColumn(eps->V,0,&v);CHKERRQ(ierr);
    if (eps->nini) {
      /* We scale the initial vector back if the initial vector was provided by users */
      ierr = VecScale(v,power->norm0);CHKERRQ(ierr);
    }
    ierr = EPSPowerUpdateFunctionB(eps,v,Bx);CHKERRQ(ierr);
    ierr = VecNorm(Bx,NORM_2,&norm);CHKERRQ(ierr);
    ierr = FirstNonzeroIdx(Bx,&power->idx,&power->p);CHKERRQ(ierr);
    ierr = Normalize(Bx,norm,power->idx,power->p,NULL);CHKERRQ(ierr);
    ierr = BVRestoreColumn(eps->V,0,&v);CHKERRQ(ierr);
  }

  ierr = STGetShift(eps->st,&sigma);CHKERRQ(ierr);    /* original shift */
  rho = sigma;

  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;
    k = eps->nconv;

    /* y = OP v */
    ierr = BVGetColumn(eps->V,k,&v);CHKERRQ(ierr);
    if (power->nonlinear) {
      ierr = VecCopy(v,eps->work[3]);CHKERRQ(ierr);
      ierr = EPSPowerApply_SNES(eps,v,y);CHKERRQ(ierr);
      ierr = VecCopy(eps->work[3],v);CHKERRQ(ierr);
    } else {
      ierr = STApply(eps->st,v,y);CHKERRQ(ierr);
    }
    ierr = BVRestoreColumn(eps->V,k,&v);CHKERRQ(ierr);

    /* purge previously converged eigenvectors */
    if (power->nonlinear) {
      ierr = EPSPowerUpdateFunctionB(eps,y,Bx);CHKERRQ(ierr);
      ierr = VecNorm(Bx,NORM_2,&norm);CHKERRQ(ierr);
      ierr = Normalize(Bx,norm,power->idx,power->p,&sign);CHKERRQ(ierr);
    } else {
      ierr = DSGetArray(eps->ds,DS_MAT_A,&T);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(eps->V,0,k);CHKERRQ(ierr);
      ierr = BVOrthogonalizeVec(eps->V,y,T+k*ld,&norm,NULL);CHKERRQ(ierr);
    }

    /* theta = (v,y)_B */
    ierr = BVSetActiveColumns(eps->V,k,k+1);CHKERRQ(ierr);
    ierr = BVDotVec(eps->V,y,&theta);CHKERRQ(ierr);
    if (!power->nonlinear) {
      T[k+k*ld] = theta;
      ierr = DSRestoreArray(eps->ds,DS_MAT_A,&T);CHKERRQ(ierr);
    }

    /* Eigenvalue: 1/|Bx| */
    if (power->nonlinear) theta = 1.0/norm*sign;

    if (power->shift_type == EPS_POWER_SHIFT_CONSTANT) { /* direct & inverse iteration */

      /* approximate eigenvalue is the Rayleigh quotient */
      eps->eigr[eps->nconv] = theta;

      /* compute relative error as ||y-theta v||_2/|theta| */
      ierr = VecCopy(y,e);CHKERRQ(ierr);
      ierr = BVGetColumn(eps->V,k,&v);CHKERRQ(ierr);
      ierr = VecAXPY(e,power->nonlinear?-1.0:-theta,v);CHKERRQ(ierr);
      ierr = BVRestoreColumn(eps->V,k,&v);CHKERRQ(ierr);
      ierr = VecNorm(e,NORM_2,&relerr);CHKERRQ(ierr);
      if (power->nonlinear)
        relerr *= PetscAbsScalar(theta);
      else
        relerr /= PetscAbsScalar(theta);

    } else {  /* RQI */

      /* delta = ||y||_B */
      delta = norm;

      /* compute relative error */
      if (rho == 0.0) relerr = PETSC_MAX_REAL;
      else relerr = 1.0 / (norm*PetscAbsScalar(rho));

      /* approximate eigenvalue is the shift */
      eps->eigr[eps->nconv] = rho;

      /* compute new shift */
      if (relerr<eps->tol) {
        rho = sigma;  /* if converged, restore original shift */
        ierr = STSetShift(eps->st,rho);CHKERRQ(ierr);
      } else {
        rho = rho + PetscConj(theta)/(delta*delta);  /* Rayleigh quotient R(v) */
        if (power->shift_type == EPS_POWER_SHIFT_WILKINSON) {
          /* beta1 is the norm of the residual associated with R(v) */
          ierr = BVGetColumn(eps->V,k,&v);CHKERRQ(ierr);
          ierr = VecAXPY(v,-PetscConj(theta)/(delta*delta),y);CHKERRQ(ierr);
          ierr = BVRestoreColumn(eps->V,k,&v);CHKERRQ(ierr);
          ierr = BVScaleColumn(eps->V,k,1.0/delta);CHKERRQ(ierr);
          ierr = BVNormColumn(eps->V,k,NORM_2,&norm1);CHKERRQ(ierr);
          beta1 = norm1;

          /* alpha2 = (e'*A*e)/(beta1*beta1), where e is the residual */
          ierr = STGetMatrix(eps->st,0,&A);CHKERRQ(ierr);
          ierr = BVGetColumn(eps->V,k,&v);CHKERRQ(ierr);
          ierr = MatMult(A,v,e);CHKERRQ(ierr);
          ierr = VecDot(v,e,&alpha2);CHKERRQ(ierr);
          ierr = BVRestoreColumn(eps->V,k,&v);CHKERRQ(ierr);
          alpha2 = alpha2 / (beta1 * beta1);

          /* choose the eigenvalue of [rho beta1; beta1 alpha2] closest to rho */
          ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
          PetscStackCallBLAS("LAPACKlaev2",LAPACKlaev2_(&rho,&beta1,&alpha2,&rt1,&rt2,&cs1,&sn1));
          ierr = PetscFPTrapPop();CHKERRQ(ierr);
          if (PetscAbsScalar(rt1-rho) < PetscAbsScalar(rt2-rho)) rho = rt1;
          else rho = rt2;
        }
        /* update operator according to new shift */
        ierr = KSPSetErrorIfNotConverged(ksp,PETSC_FALSE);CHKERRQ(ierr);
        ierr = STSetShift(eps->st,rho);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
        if (reason) {
          ierr = PetscInfo(eps,"Factorization failed, repeat with a perturbed shift\n");CHKERRQ(ierr);
          rho *= 1+10*PETSC_MACHINE_EPSILON;
          ierr = STSetShift(eps->st,rho);CHKERRQ(ierr);
          ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
          if (reason) SETERRQ(PetscObjectComm((PetscObject)ksp),PETSC_ERR_CONV_FAILED,"Second factorization failed");
        }
        ierr = KSPSetErrorIfNotConverged(ksp,PETSC_TRUE);CHKERRQ(ierr);
      }
    }
    eps->errest[eps->nconv] = relerr;

    /* normalize */
    if (!power->nonlinear) { ierr = Normalize(y,norm,power->idx,power->p,NULL);CHKERRQ(ierr); }
    ierr = BVInsertVec(eps->V,k,y);CHKERRQ(ierr);

    if (power->update) {
      ierr = SNESGetConvergedReason(power->snes,&snesreason);CHKERRQ(ierr);
    }
    /* if relerr<tol, accept eigenpair */
    if (relerr<eps->tol || (power->update && snesreason > 0)) {
      eps->nconv = eps->nconv + 1;
      if (eps->nconv<eps->nev) {
        ierr = EPSGetStartVector(eps,eps->nconv,&breakdown);CHKERRQ(ierr);
        if (breakdown) {
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          ierr = PetscInfo(eps,"Unable to generate more start vectors\n");CHKERRQ(ierr);
          break;
        }
      }
    }
    ierr = EPSMonitor(eps,eps->its,eps->nconv,eps->eigr,eps->eigi,eps->errest,PetscMin(eps->nconv+1,eps->nev));CHKERRQ(ierr);
    ierr = (*eps->stopping)(eps,eps->its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);
  }

  if (power->nonlinear) {
    ierr = PetscObjectCompose((PetscObject)power->snes,"eps",NULL);CHKERRQ(ierr);
  } else {
    ierr = DSSetDimensions(eps->ds,eps->nconv,0,0,0);CHKERRQ(ierr);
    ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSolve_TS_Power(EPS eps)
{
  PetscErrorCode     ierr;
  EPS_POWER          *power = (EPS_POWER*)eps->data;
  PetscInt           k,ld;
  Vec                v,w,y,e,z;
  KSP                ksp;
  PetscReal          relerr=1.0,relerrl,delta;
  PetscScalar        theta,rho,alpha,sigma;
  PetscBool          breakdown,breakdownl;
  KSPConvergedReason reason;

  PetscFunctionBegin;
  e = eps->work[0];
  v = eps->work[1];
  w = eps->work[2];

  if (power->shift_type != EPS_POWER_SHIFT_CONSTANT) { ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr); }
  ierr = DSGetLeadingDimension(eps->ds,&ld);CHKERRQ(ierr);
  ierr = EPSGetStartVector(eps,0,NULL);CHKERRQ(ierr);
  ierr = EPSGetLeftStartVector(eps,0,NULL);CHKERRQ(ierr);
  ierr = BVBiorthonormalizeColumn(eps->V,eps->W,0,NULL);CHKERRQ(ierr);
  ierr = BVCopyVec(eps->V,0,v);CHKERRQ(ierr);
  ierr = BVCopyVec(eps->W,0,w);CHKERRQ(ierr);
  ierr = STGetShift(eps->st,&sigma);CHKERRQ(ierr);    /* original shift */
  rho = sigma;

  while (eps->reason == EPS_CONVERGED_ITERATING) {
    eps->its++;
    k = eps->nconv;

    /* y = OP v, z = OP' w */
    ierr = BVGetColumn(eps->V,k,&y);CHKERRQ(ierr);
    ierr = STApply(eps->st,v,y);CHKERRQ(ierr);
    ierr = BVRestoreColumn(eps->V,k,&y);CHKERRQ(ierr);
    ierr = BVGetColumn(eps->W,k,&z);CHKERRQ(ierr);
    ierr = STApplyHermitianTranspose(eps->st,w,z);CHKERRQ(ierr);
    ierr = BVRestoreColumn(eps->W,k,&z);CHKERRQ(ierr);

    /* purge previously converged eigenvectors */
    ierr = BVBiorthogonalizeColumn(eps->V,eps->W,k);CHKERRQ(ierr);

    /* theta = (w,y)_B */
    ierr = BVSetActiveColumns(eps->V,k,k+1);CHKERRQ(ierr);
    ierr = BVDotVec(eps->V,w,&theta);CHKERRQ(ierr);
    theta = PetscConj(theta);

    if (power->shift_type == EPS_POWER_SHIFT_CONSTANT) { /* direct & inverse iteration */

      /* approximate eigenvalue is the Rayleigh quotient */
      eps->eigr[eps->nconv] = theta;

      /* compute relative errors as ||y-theta v||_2/|theta| and ||z-conj(theta) w||_2/|theta|*/
      ierr = BVCopyVec(eps->V,k,e);CHKERRQ(ierr);
      ierr = VecAXPY(e,-theta,v);CHKERRQ(ierr);
      ierr = VecNorm(e,NORM_2,&relerr);CHKERRQ(ierr);
      ierr = BVCopyVec(eps->W,k,e);CHKERRQ(ierr);
      ierr = VecAXPY(e,-PetscConj(theta),w);CHKERRQ(ierr);
      ierr = VecNorm(e,NORM_2,&relerrl);CHKERRQ(ierr);
      relerr = PetscMax(relerr,relerrl)/PetscAbsScalar(theta);
    }

    /* normalize */
    ierr = BVSetActiveColumns(eps->V,k,k+1);CHKERRQ(ierr);
    ierr = BVGetColumn(eps->W,k,&z);CHKERRQ(ierr);
    ierr = BVDotVec(eps->V,z,&alpha);CHKERRQ(ierr);
    ierr = BVRestoreColumn(eps->W,k,&z);CHKERRQ(ierr);
    delta = PetscSqrtReal(PetscAbsScalar(alpha));
    if (delta==0.0) SETERRQ(PetscObjectComm((PetscObject)eps),1,"Breakdown in two-sided Power/RQI");
    ierr = BVScaleColumn(eps->V,k,1.0/PetscConj(alpha/delta));CHKERRQ(ierr);
    ierr = BVScaleColumn(eps->W,k,1.0/delta);CHKERRQ(ierr);
    ierr = BVCopyVec(eps->V,k,v);CHKERRQ(ierr);
    ierr = BVCopyVec(eps->W,k,w);CHKERRQ(ierr);

    if (power->shift_type == EPS_POWER_SHIFT_RAYLEIGH) { /* RQI */

      /* compute relative error */
      if (rho == 0.0) relerr = PETSC_MAX_REAL;
      else relerr = 1.0 / PetscAbsScalar(delta*rho);

      /* approximate eigenvalue is the shift */
      eps->eigr[eps->nconv] = rho;

      /* compute new shift */
      if (relerr<eps->tol) {
        rho = sigma;  /* if converged, restore original shift */
        ierr = STSetShift(eps->st,rho);CHKERRQ(ierr);
      } else {
        rho = rho + PetscConj(theta)/(delta*delta);  /* Rayleigh quotient R(v) */
        /* update operator according to new shift */
        ierr = KSPSetErrorIfNotConverged(ksp,PETSC_FALSE);CHKERRQ(ierr);
        ierr = STSetShift(eps->st,rho);CHKERRQ(ierr);
        ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
        if (reason) {
          ierr = PetscInfo(eps,"Factorization failed, repeat with a perturbed shift\n");CHKERRQ(ierr);
          rho *= 1+10*PETSC_MACHINE_EPSILON;
          ierr = STSetShift(eps->st,rho);CHKERRQ(ierr);
          ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
          if (reason) SETERRQ(PetscObjectComm((PetscObject)ksp),PETSC_ERR_CONV_FAILED,"Second factorization failed");
        }
        ierr = KSPSetErrorIfNotConverged(ksp,PETSC_TRUE);CHKERRQ(ierr);
      }
    }
    eps->errest[eps->nconv] = relerr;

    /* if relerr<tol, accept eigenpair */
    if (relerr<eps->tol) {
      eps->nconv = eps->nconv + 1;
      if (eps->nconv<eps->nev) {
        ierr = EPSGetStartVector(eps,eps->nconv,&breakdown);CHKERRQ(ierr);
        ierr = EPSGetLeftStartVector(eps,eps->nconv,&breakdownl);CHKERRQ(ierr);
        if (breakdown || breakdownl) {
          eps->reason = EPS_DIVERGED_BREAKDOWN;
          ierr = PetscInfo(eps,"Unable to generate more start vectors\n");CHKERRQ(ierr);
          break;
        }
        ierr = BVBiorthonormalizeColumn(eps->V,eps->W,eps->nconv,NULL);CHKERRQ(ierr);
      }
    }
    ierr = EPSMonitor(eps,eps->its,eps->nconv,eps->eigr,eps->eigi,eps->errest,PetscMin(eps->nconv+1,eps->nev));CHKERRQ(ierr);
    ierr = (*eps->stopping)(eps,eps->its,eps->max_it,eps->nconv,eps->nev,&eps->reason,eps->stoppingctx);CHKERRQ(ierr);
  }

  ierr = DSSetDimensions(eps->ds,eps->nconv,0,0,0);CHKERRQ(ierr);
  ierr = DSSetState(eps->ds,DS_STATE_RAW);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSStopping_Power(EPS eps,PetscInt its,PetscInt max_it,PetscInt nconv,PetscInt nev,EPSConvergedReason *reason,void *ctx)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  SNESConvergedReason snesreason;

  PetscFunctionBegin;
  if (power->update) {
    ierr = SNESGetConvergedReason(power->snes,&snesreason);CHKERRQ(ierr);
    if (snesreason < 0) {
      *reason = EPS_DIVERGED_BREAKDOWN;
      PetscFunctionReturn(0);
    }
  }
  ierr = EPSStoppingBasic(eps,its,max_it,nconv,nev,reason,ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSBackTransform_Power(EPS eps)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (power->nonlinear) eps->eigr[0] = 1.0/eps->eigr[0];
  else if (power->shift_type == EPS_POWER_SHIFT_CONSTANT) {
    ierr = EPSBackTransform_Default(eps);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetFromOptions_Power(PetscOptionItems *PetscOptionsObject,EPS eps)
{
  PetscErrorCode    ierr;
  EPS_POWER         *power = (EPS_POWER*)eps->data;
  PetscBool         flg,val;
  EPSPowerShiftType shift;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"EPS Power Options");CHKERRQ(ierr);

    ierr = PetscOptionsEnum("-eps_power_shift_type","Shift type","EPSPowerSetShiftType",EPSPowerShiftTypes,(PetscEnum)power->shift_type,(PetscEnum*)&shift,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSPowerSetShiftType(eps,shift);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-eps_power_nonlinear","Use nonlinear inverse iteration","EPSPowerSetNonlinear",power->nonlinear,&val,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSPowerSetNonlinear(eps,val);CHKERRQ(ierr); }

    ierr = PetscOptionsBool("-eps_power_update","Update residual monolithically","EPSPowerSetUpdate",power->update,&val,&flg);CHKERRQ(ierr);
    if (flg) { ierr = EPSPowerSetUpdate(eps,val);CHKERRQ(ierr); }

  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerSetShiftType_Power(EPS eps,EPSPowerShiftType shift)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  switch (shift) {
    case EPS_POWER_SHIFT_CONSTANT:
    case EPS_POWER_SHIFT_RAYLEIGH:
    case EPS_POWER_SHIFT_WILKINSON:
      if (power->shift_type != shift) {
        power->shift_type = shift;
        eps->state = EPS_STATE_INITIAL;
      }
      break;
    default:
      SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_OUTOFRANGE,"Invalid shift type");
  }
  PetscFunctionReturn(0);
}

/*@
   EPSPowerSetShiftType - Sets the type of shifts used during the power
   iteration. This can be used to emulate the Rayleigh Quotient Iteration
   (RQI) method.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  shift - the type of shift

   Options Database Key:
.  -eps_power_shift_type - Sets the shift type (either 'constant' or
                           'rayleigh' or 'wilkinson')

   Notes:
   By default, shifts are constant (EPS_POWER_SHIFT_CONSTANT) and the iteration
   is the simple power method (or inverse iteration if a shift-and-invert
   transformation is being used).

   A variable shift can be specified (EPS_POWER_SHIFT_RAYLEIGH or
   EPS_POWER_SHIFT_WILKINSON). In this case, the iteration behaves rather like
   a cubic converging method such as RQI.

   Level: advanced

.seealso: EPSPowerGetShiftType(), STSetShift(), EPSPowerShiftType
@*/
PetscErrorCode EPSPowerSetShiftType(EPS eps,EPSPowerShiftType shift)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveEnum(eps,shift,2);
  ierr = PetscTryMethod(eps,"EPSPowerSetShiftType_C",(EPS,EPSPowerShiftType),(eps,shift));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerGetShiftType_Power(EPS eps,EPSPowerShiftType *shift)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  *shift = power->shift_type;
  PetscFunctionReturn(0);
}

/*@
   EPSPowerGetShiftType - Gets the type of shifts used during the power
   iteration.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Input Parameter:
.  shift - the type of shift

   Level: advanced

.seealso: EPSPowerSetShiftType(), EPSPowerShiftType
@*/
PetscErrorCode EPSPowerGetShiftType(EPS eps,EPSPowerShiftType *shift)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(shift,2);
  ierr = PetscUseMethod(eps,"EPSPowerGetShiftType_C",(EPS,EPSPowerShiftType*),(eps,shift));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerSetNonlinear_Power(EPS eps,PetscBool nonlinear)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (power->nonlinear != nonlinear) {
    power->nonlinear = nonlinear;
    eps->useds = PetscNot(nonlinear);
    eps->state = EPS_STATE_INITIAL;
  }
  PetscFunctionReturn(0);
}

/*@
   EPSPowerSetNonlinear - Sets a flag to indicate that the problem is nonlinear.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  nonlinear - whether the problem is nonlinear or not

   Options Database Key:
.  -eps_power_nonlinear - Sets the nonlinear flag

   Notes:
   If this flag is set, the solver assumes that the problem is nonlinear,
   that is, the operators that define the eigenproblem are not constant
   matrices, but depend on the eigenvector: A(x)*x=lambda*B(x)*x. This is
   different from the case of nonlinearity with respect to the eigenvalue
   (use the NEP solver class for this kind of problems).

   The way in which nonlinear operators are specified is very similar to
   the case of PETSc's SNES solver. The difference is that the callback
   functions are provided via composed functions "formFunction" and
   "formJacobian" in each of the matrix objects passed as arguments of
   EPSSetOperators(). The application context required for these functions
   can be attached via a composed PetscContainer.

   Level: advanced

.seealso: EPSPowerGetNonlinear(), EPSSetOperators()
@*/
PetscErrorCode EPSPowerSetNonlinear(EPS eps,PetscBool nonlinear)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,nonlinear,2);
  ierr = PetscTryMethod(eps,"EPSPowerSetNonlinear_C",(EPS,PetscBool),(eps,nonlinear));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerGetNonlinear_Power(EPS eps,PetscBool *nonlinear)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  *nonlinear = power->nonlinear;
  PetscFunctionReturn(0);
}

/*@
   EPSPowerGetNonlinear - Returns a flag indicating if the problem is nonlinear.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Input Parameter:
.  nonlinear - the nonlinear flag

   Level: advanced

.seealso: EPSPowerSetUpdate(), EPSPowerSetNonlinear()
@*/
PetscErrorCode EPSPowerGetNonlinear(EPS eps,PetscBool *nonlinear)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(nonlinear,2);
  ierr = PetscUseMethod(eps,"EPSPowerGetNonlinear_C",(EPS,PetscBool*),(eps,nonlinear));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerSetUpdate_Power(EPS eps,PetscBool update)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (!power->nonlinear) SETERRQ(PetscObjectComm((PetscObject)eps),PETSC_ERR_ARG_INCOMP,"This option does not make sense for linear problems");
  power->update = update;
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSPowerSetUpdate - Sets a flag to indicate that the residual is updated monolithically
   for nonlinear problems. This potentially has a better convergence.

   Logically Collective on eps

   Input Parameters:
+  eps - the eigenproblem solver context
-  update - whether the residual is updated monolithically or not

   Options Database Key:
.  -eps_power_update - Sets the update flag

   Level: advanced

.seealso: EPSPowerGetUpdate(), EPSPowerGetNonlinear(), EPSSetOperators()
@*/
PetscErrorCode EPSPowerSetUpdate(EPS eps,PetscBool update)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidLogicalCollectiveBool(eps,update,2);
  ierr = PetscTryMethod(eps,"EPSPowerSetUpdate_C",(EPS,PetscBool),(eps,update));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerGetUpdate_Power(EPS eps,PetscBool *update)
{
  EPS_POWER *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  *update = power->update;
  PetscFunctionReturn(0);
}

/*@
   EPSPowerGetUpdate - Returns a flag indicating if the residual is updated monolithically
   for nonlinear problems.

   Not Collective

   Input Parameter:
.  eps - the eigenproblem solver context

   Input Parameter:
.  update - the update flag

   Level: advanced

.seealso: EPSPowerSetUpdate(), EPSPowerSetNonlinear()
@*/
PetscErrorCode EPSPowerGetUpdate(EPS eps,PetscBool *update)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidBoolPointer(update,2);
  ierr = PetscUseMethod(eps,"EPSPowerGetUpdate_C",(EPS,PetscBool*),(eps,update));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerSetSNES_Power(EPS eps,SNES snes)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  ierr = PetscObjectReference((PetscObject)snes);CHKERRQ(ierr);
  ierr = SNESDestroy(&power->snes);CHKERRQ(ierr);
  power->snes = snes;
  ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)power->snes);CHKERRQ(ierr);
  eps->state = EPS_STATE_INITIAL;
  PetscFunctionReturn(0);
}

/*@
   EPSPowerSetSNES - Associate a nonlinear solver object (SNES) to the
   eigenvalue solver (to be used in nonlinear inverse iteration).

   Collective on eps

   Input Parameters:
+  eps  - the eigenvalue solver
-  snes - the nonlinear solver object

   Level: advanced

.seealso: EPSPowerGetSNES()
@*/
PetscErrorCode EPSPowerSetSNES(EPS eps,SNES snes)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidHeaderSpecific(snes,SNES_CLASSID,2);
  PetscCheckSameComm(eps,1,snes,2);
  ierr = PetscTryMethod(eps,"EPSPowerSetSNES_C",(EPS,SNES),(eps,snes));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode EPSPowerGetSNES_Power(EPS eps,SNES *snes)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (!power->snes) {
    ierr = SNESCreate(PetscObjectComm((PetscObject)eps),&power->snes);CHKERRQ(ierr);
    ierr = PetscObjectIncrementTabLevel((PetscObject)power->snes,(PetscObject)eps,1);CHKERRQ(ierr);
    ierr = SNESSetOptionsPrefix(power->snes,((PetscObject)eps)->prefix);CHKERRQ(ierr);
    ierr = SNESAppendOptionsPrefix(power->snes,"eps_power_");CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)eps,(PetscObject)power->snes);CHKERRQ(ierr);
    ierr = PetscObjectSetOptions((PetscObject)power->snes,((PetscObject)eps)->options);CHKERRQ(ierr);
  }
  *snes = power->snes;
  PetscFunctionReturn(0);
}

/*@
   EPSPowerGetSNES - Retrieve the nonlinear solver object (SNES) associated
   with the eigenvalue solver.

   Not Collective

   Input Parameter:
.  eps - the eigenvalue solver

   Output Parameter:
.  snes - the nonlinear solver object

   Level: advanced

.seealso: EPSPowerSetSNES()
@*/
PetscErrorCode EPSPowerGetSNES(EPS eps,SNES *snes)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(snes,2);
  ierr = PetscUseMethod(eps,"EPSPowerGetSNES_C",(EPS,SNES*),(eps,snes));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSReset_Power(EPS eps)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (power->snes) { ierr = SNESReset(power->snes);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSDestroy_Power(EPS eps)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;

  PetscFunctionBegin;
  if (power->nonlinear) {
    ierr = SNESDestroy(&power->snes);CHKERRQ(ierr);
  }
  ierr = PetscFree(eps->data);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetShiftType_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetShiftType_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetNonlinear_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetNonlinear_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetUpdate_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetUpdate_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetSNES_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetSNES_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode EPSView_Power(EPS eps,PetscViewer viewer)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  PetscBool      isascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&isascii);CHKERRQ(ierr);
  if (isascii) {
    if (power->nonlinear) {
      ierr = PetscViewerASCIIPrintf(viewer,"  using nonlinear inverse iteration\n");CHKERRQ(ierr);
      if (power->update) {
        ierr = PetscViewerASCIIPrintf(viewer,"  updating the residual monolithically\n");CHKERRQ(ierr);
      }
      if (!power->snes) { ierr = EPSPowerGetSNES(eps,&power->snes);CHKERRQ(ierr); }
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = SNESView(power->snes,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"  %s shifts\n",EPSPowerShiftTypes[power->shift_type]);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSComputeVectors_Power(EPS eps)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  PetscReal      norm;
  PetscInt       i;

  PetscFunctionBegin;
  if (eps->twosided) {
    ierr = EPSComputeVectors_Twosided(eps);CHKERRQ(ierr);
    /* normalize (no need to take care of 2x2 blocks */
    for (i=0;i<eps->nconv;i++) {
      ierr = BVNormColumn(eps->V,i,NORM_2,&norm);CHKERRQ(ierr);
      ierr = BVScaleColumn(eps->V,i,1.0/norm);CHKERRQ(ierr);
      ierr = BVNormColumn(eps->W,i,NORM_2,&norm);CHKERRQ(ierr);
      ierr = BVScaleColumn(eps->W,i,1.0/norm);CHKERRQ(ierr);
    }
  } else if (!power->nonlinear) {
    ierr = EPSComputeVectors_Schur(eps);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EPSSetDefaultST_Power(EPS eps)
{
  PetscErrorCode ierr;
  EPS_POWER      *power = (EPS_POWER*)eps->data;
  KSP            ksp;
  PC             pc;

  PetscFunctionBegin;
  if (power->nonlinear) {
    eps->categ=EPS_CATEGORY_PRECOND;
    ierr = STGetKSP(eps->st,&ksp);CHKERRQ(ierr);
    /* Set ST as STPRECOND so it can carry one preconditioning matrix
     * It is useful when A and B are shell matrices
     */
    ierr = STSetType(eps->st,STPRECOND);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode EPSCreate_Power(EPS eps)
{
  EPS_POWER      *ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(eps,&ctx);CHKERRQ(ierr);
  eps->data = (void*)ctx;

  eps->useds = PETSC_TRUE;
  eps->hasts = PETSC_TRUE;
  eps->categ = EPS_CATEGORY_OTHER;

  eps->ops->setup          = EPSSetUp_Power;
  eps->ops->setfromoptions = EPSSetFromOptions_Power;
  eps->ops->reset          = EPSReset_Power;
  eps->ops->destroy        = EPSDestroy_Power;
  eps->ops->view           = EPSView_Power;
  eps->ops->backtransform  = EPSBackTransform_Power;
  eps->ops->computevectors = EPSComputeVectors_Power;
  eps->ops->setdefaultst   = EPSSetDefaultST_Power;
  eps->stopping            = EPSStopping_Power;

  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetShiftType_C",EPSPowerSetShiftType_Power);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetShiftType_C",EPSPowerGetShiftType_Power);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetNonlinear_C",EPSPowerSetNonlinear_Power);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetNonlinear_C",EPSPowerGetNonlinear_Power);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetUpdate_C",EPSPowerSetUpdate_Power);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetUpdate_C",EPSPowerGetUpdate_Power);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerSetSNES_C",EPSPowerSetSNES_Power);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)eps,"EPSPowerGetSNES_C",EPSPowerGetSNES_Power);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
