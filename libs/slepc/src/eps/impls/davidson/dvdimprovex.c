/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   SLEPc eigensolver: "davidson"

   Step: improve the eigenvectors X
*/

#include "davidson.h"
#include <slepcblaslapack.h>

/**** JD update step (I - Kfg'/(g'Kf)) K(A - sB) (I - Kfg'/(g'Kf)) t = (I - Kfg'/(g'Kf))r  ****/

typedef struct {
  PetscInt     size_X;
  KSP          ksp;                /* correction equation solver */
  Vec          friends;            /* reference vector for composite vectors */
  PetscScalar  theta[4],thetai[2]; /* the shifts used in the correction eq. */
  PetscInt     maxits;             /* maximum number of iterations */
  PetscInt     r_s,r_e;            /* the selected eigenpairs to improve */
  PetscInt     ksp_max_size;       /* the ksp maximum subvectors size */
  PetscReal    tol;                /* the maximum solution tolerance */
  PetscReal    lastTol;            /* last tol for dynamic stopping criterion */
  PetscReal    fix;                /* tolerance for using the approx. eigenvalue */
  PetscBool    dynamic;            /* if the dynamic stopping criterion is applied */
  dvdDashboard *d;                 /* the currect dvdDashboard reference */
  PC           old_pc;             /* old pc in ksp */
  BV           KZ;                 /* KZ vecs for the projector KZ*inv(X'*KZ)*X' */
  BV           U;                  /* new X vectors */
  PetscScalar  *XKZ;               /* X'*KZ */
  PetscScalar  *iXKZ;              /* inverse of XKZ */
  PetscInt     ldXKZ;              /* leading dimension of XKZ */
  PetscInt     size_iXKZ;          /* size of iXKZ */
  PetscInt     ldiXKZ;             /* leading dimension of iXKZ */
  PetscInt     size_cX;            /* last value of d->size_cX */
  PetscInt     old_size_X;         /* last number of improved vectors */
  PetscBLASInt *iXKZPivots;        /* array of pivots */
} dvdImprovex_jd;

/*
   Compute (I - KZ*iXKZ*X')*V where,
   V, the vectors to apply the projector,
   cV, the number of vectors in V,
*/
static PetscErrorCode dvd_improvex_apply_proj(dvdDashboard *d,Vec *V,PetscInt cV)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscInt       i,ldh,k,l;
  PetscScalar    *h;
  PetscBLASInt   cV_,n,info,ld;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
#endif

  PetscFunctionBegin;
  if (cV > 2) SETERRQ(PETSC_COMM_SELF,1,"Consistency broken");

  /* h <- X'*V */
  ierr = PetscMalloc1(data->size_iXKZ*cV,&h);CHKERRQ(ierr);
  ldh = data->size_iXKZ;
  ierr = BVGetActiveColumns(data->U,&l,&k);CHKERRQ(ierr);
  if (ldh!=k) SETERRQ(PETSC_COMM_SELF,1,"Consistency broken");
  ierr = BVSetActiveColumns(data->U,0,k);CHKERRQ(ierr);
  for (i=0;i<cV;i++) {
    ierr = BVDotVec(data->U,V[i],&h[ldh*i]);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    for (j=0; j<k; j++) h[ldh*i+j] = PetscConj(h[ldh*i+j]);
#endif
  }
  ierr = BVSetActiveColumns(data->U,l,k);CHKERRQ(ierr);

  /* h <- iXKZ\h */
  ierr = PetscBLASIntCast(cV,&cV_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(data->size_iXKZ,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(data->ldiXKZ,&ld);CHKERRQ(ierr);
  PetscValidScalarPointer(data->iXKZ,0);
  ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgetrs",LAPACKgetrs_("N",&n,&cV_,data->iXKZ,&ld,data->iXKZPivots,h,&n,&info));
  ierr = PetscFPTrapPop();CHKERRQ(ierr);
  SlepcCheckLapackInfo("getrs",info);

  /* V <- V - KZ*h */
  ierr = BVSetActiveColumns(data->KZ,0,k);CHKERRQ(ierr);
  for (i=0;i<cV;i++) {
    ierr = BVMultVec(data->KZ,-1.0,1.0,V[i],&h[ldh*i]);CHKERRQ(ierr);
  }
  ierr = BVSetActiveColumns(data->KZ,l,k);CHKERRQ(ierr);
  ierr = PetscFree(h);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   Compute (I - X*iXKZ*KZ')*V where,
   V, the vectors to apply the projector,
   cV, the number of vectors in V,
*/
static PetscErrorCode dvd_improvex_applytrans_proj(dvdDashboard *d,Vec *V,PetscInt cV)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscInt       i,ldh,k,l;
  PetscScalar    *h;
  PetscBLASInt   cV_, n, info, ld;
#if defined(PETSC_USE_COMPLEX)
  PetscInt       j;
#endif

  PetscFunctionBegin;
  if (cV > 2) SETERRQ(PETSC_COMM_SELF,1,"Consistency broken");

  /* h <- KZ'*V */
  ierr = PetscMalloc1(data->size_iXKZ*cV,&h);CHKERRQ(ierr);
  ldh = data->size_iXKZ;
  ierr = BVGetActiveColumns(data->U,&l,&k);CHKERRQ(ierr);
  if (ldh!=k) SETERRQ(PETSC_COMM_SELF,1,"Consistency broken");
  ierr = BVSetActiveColumns(data->KZ,0,k);CHKERRQ(ierr);
  for (i=0;i<cV;i++) {
    ierr = BVDotVec(data->KZ,V[i],&h[ldh*i]);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
    for (j=0;j<k;j++) h[ldh*i+j] = PetscConj(h[ldh*i+j]);
#endif
  }
  ierr = BVSetActiveColumns(data->KZ,l,k);CHKERRQ(ierr);

  /* h <- iXKZ\h */
  ierr = PetscBLASIntCast(cV,&cV_);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(data->size_iXKZ,&n);CHKERRQ(ierr);
  ierr = PetscBLASIntCast(data->ldiXKZ,&ld);CHKERRQ(ierr);
  PetscValidScalarPointer(data->iXKZ,0);
  ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgetrs",LAPACKgetrs_("C",&n,&cV_,data->iXKZ,&ld,data->iXKZPivots,h,&n,&info));
  ierr = PetscFPTrapPop();CHKERRQ(ierr);
  SlepcCheckLapackInfo("getrs",info);

  /* V <- V - U*h */
  ierr = BVSetActiveColumns(data->U,0,k);CHKERRQ(ierr);
  for (i=0;i<cV;i++) {
    ierr = BVMultVec(data->U,-1.0,1.0,V[i],&h[ldh*i]);CHKERRQ(ierr);
  }
  ierr = BVSetActiveColumns(data->U,l,k);CHKERRQ(ierr);
  ierr = PetscFree(h);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_improvex_jd_end(dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;

  PetscFunctionBegin;
  ierr = VecDestroy(&data->friends);CHKERRQ(ierr);

  /* Restore the pc of ksp */
  if (data->old_pc) {
    ierr = KSPSetPC(data->ksp, data->old_pc);CHKERRQ(ierr);
    ierr = PCDestroy(&data->old_pc);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_improvex_jd_d(dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;

  PetscFunctionBegin;
  /* Free local data and objects */
  ierr = PetscFree(data->XKZ);CHKERRQ(ierr);
  ierr = PetscFree(data->iXKZ);CHKERRQ(ierr);
  ierr = PetscFree(data->iXKZPivots);CHKERRQ(ierr);
  ierr = BVDestroy(&data->KZ);CHKERRQ(ierr);
  ierr = BVDestroy(&data->U);CHKERRQ(ierr);
  ierr = PetscFree(data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   y <- theta[1]A*x - theta[0]*B*x
   auxV, two auxiliary vectors
 */
PETSC_STATIC_INLINE PetscErrorCode dvd_aux_matmult(dvdImprovex_jd *data,const Vec *x,const Vec *y)
{
  PetscErrorCode ierr;
  PetscInt       n,i;
  const Vec      *Bx;
  Vec            *auxV;

  PetscFunctionBegin;
  n = data->r_e - data->r_s;
  for (i=0;i<n;i++) {
    ierr = MatMult(data->d->A,x[i],y[i]);CHKERRQ(ierr);
  }

  ierr = SlepcVecPoolGetVecs(data->d->auxV,2,&auxV);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (data->d->eigi[data->r_s+i] != 0.0) {
      if (data->d->B) {
        ierr = MatMult(data->d->B,x[i],auxV[0]);CHKERRQ(ierr);
        ierr = MatMult(data->d->B,x[i+1],auxV[1]);CHKERRQ(ierr);
        Bx = auxV;
      } else Bx = &x[i];

      /* y_i   <- [ t_2i+1*A*x_i   - t_2i*Bx_i + ti_i*Bx_i+1;
         y_i+1      t_2i+1*A*x_i+1 - ti_i*Bx_i - t_2i*Bx_i+1  ] */
      ierr = VecAXPBYPCZ(y[i],-data->theta[2*i],data->thetai[i],data->theta[2*i+1],Bx[0],Bx[1]);CHKERRQ(ierr);
      ierr = VecAXPBYPCZ(y[i+1],-data->thetai[i],-data->theta[2*i],data->theta[2*i+1],Bx[0],Bx[1]);CHKERRQ(ierr);
      i++;
    } else
#endif
    {
      if (data->d->B) {
        ierr = MatMult(data->d->B,x[i],auxV[0]);CHKERRQ(ierr);
        Bx = auxV;
      } else Bx = &x[i];
      ierr = VecAXPBY(y[i],-data->theta[i*2],data->theta[i*2+1],Bx[0]);CHKERRQ(ierr);
    }
  }
  ierr = SlepcVecPoolRestoreVecs(data->d->auxV,2,&auxV);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   y <- theta[1]'*A'*x - theta[0]'*B'*x
 */
PETSC_STATIC_INLINE PetscErrorCode dvd_aux_matmulttrans(dvdImprovex_jd *data,const Vec *x,const Vec *y)
{
  PetscErrorCode ierr;
  PetscInt       n,i;
  const Vec      *Bx;
  Vec            *auxV;

  PetscFunctionBegin;
  n = data->r_e - data->r_s;
  for (i=0;i<n;i++) {
    ierr = MatMultTranspose(data->d->A,x[i],y[i]);CHKERRQ(ierr);
  }

  ierr = SlepcVecPoolGetVecs(data->d->auxV,2,&auxV);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (data->d->eigi[data->r_s+i] != 0.0) {
      if (data->d->B) {
        ierr = MatMultTranspose(data->d->B,x[i],auxV[0]);CHKERRQ(ierr);
        ierr = MatMultTranspose(data->d->B,x[i+1],auxV[1]);CHKERRQ(ierr);
        Bx = auxV;
      } else Bx = &x[i];

      /* y_i   <- [ t_2i+1*A*x_i   - t_2i*Bx_i - ti_i*Bx_i+1;
         y_i+1      t_2i+1*A*x_i+1 + ti_i*Bx_i - t_2i*Bx_i+1  ] */
      ierr = VecAXPBYPCZ(y[i],-data->theta[2*i],-data->thetai[i],data->theta[2*i+1],Bx[0],Bx[1]);CHKERRQ(ierr);
      ierr = VecAXPBYPCZ(y[i+1],data->thetai[i],-data->theta[2*i],data->theta[2*i+1],Bx[0],Bx[1]);CHKERRQ(ierr);
      i++;
    } else
#endif
    {
      if (data->d->B) {
        ierr = MatMultTranspose(data->d->B,x[i],auxV[0]);CHKERRQ(ierr);
        Bx = auxV;
      } else Bx = &x[i];
      ierr = VecAXPBY(y[i],PetscConj(-data->theta[i*2]),PetscConj(data->theta[i*2+1]),Bx[0]);CHKERRQ(ierr);
    }
  }
  ierr = SlepcVecPoolRestoreVecs(data->d->auxV,2,&auxV);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PCApplyBA_dvd(PC pc,PCSide side,Vec in,Vec out,Vec w)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data;
  PetscInt       n,i;
  const Vec      *inx,*outx,*wx;
  Vec            *auxV;
  Mat            A;

  PetscFunctionBegin;
  ierr = PCGetOperators(pc,&A,NULL);CHKERRQ(ierr);
  ierr = MatShellGetContext(A,(void**)&data);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(in,NULL,&inx);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(out,NULL,&outx);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(w,NULL,&wx);CHKERRQ(ierr);
  n = data->r_e - data->r_s;
  ierr = SlepcVecPoolGetVecs(data->d->auxV,n,&auxV);CHKERRQ(ierr);
  switch (side) {
  case PC_LEFT:
    /* aux <- theta[1]A*in - theta[0]*B*in */
    ierr = dvd_aux_matmult(data,inx,auxV);CHKERRQ(ierr);

    /* out <- K * aux */
    for (i=0;i<n;i++) {
      ierr = data->d->improvex_precond(data->d,data->r_s+i,auxV[i],outx[i]);CHKERRQ(ierr);
    }
    break;
  case PC_RIGHT:
    /* aux <- K * in */
    for (i=0;i<n;i++) {
      ierr = data->d->improvex_precond(data->d,data->r_s+i,inx[i],auxV[i]);CHKERRQ(ierr);
    }

    /* out <- theta[1]A*auxV - theta[0]*B*auxV */
    ierr = dvd_aux_matmult(data,auxV,outx);CHKERRQ(ierr);
    break;
  case PC_SYMMETRIC:
    /* aux <- K^{1/2} * in */
    for (i=0;i<n;i++) {
      ierr = PCApplySymmetricRight(data->old_pc,inx[i],auxV[i]);CHKERRQ(ierr);
    }

    /* wx <- theta[1]A*auxV - theta[0]*B*auxV */
    ierr = dvd_aux_matmult(data,auxV,wx);CHKERRQ(ierr);

    /* aux <- K^{1/2} * in */
    for (i=0;i<n;i++) {
      ierr = PCApplySymmetricLeft(data->old_pc,wx[i],outx[i]);CHKERRQ(ierr);
    }
    break;
  default:
    SETERRQ(PETSC_COMM_SELF,1,"Unsupported KSP side");
  }
  /* out <- out - v*(u'*out) */
  ierr = dvd_improvex_apply_proj(data->d,(Vec*)outx,n);CHKERRQ(ierr);
  ierr = SlepcVecPoolRestoreVecs(data->d->auxV,n,&auxV);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PCApply_dvd(PC pc,Vec in,Vec out)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data;
  PetscInt       n,i;
  const Vec      *inx, *outx;
  Mat            A;

  PetscFunctionBegin;
  ierr = PCGetOperators(pc,&A,NULL);CHKERRQ(ierr);
  ierr = MatShellGetContext(A,(void**)&data);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(in,NULL,&inx);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(out,NULL,&outx);CHKERRQ(ierr);
  n = data->r_e - data->r_s;
  /* out <- K * in */
  for (i=0;i<n;i++) {
    ierr = data->d->improvex_precond(data->d,data->r_s+i,inx[i],outx[i]);CHKERRQ(ierr);
  }
  /* out <- out - v*(u'*out) */
  ierr = dvd_improvex_apply_proj(data->d,(Vec*)outx,n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PCApplyTranspose_dvd(PC pc,Vec in,Vec out)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data;
  PetscInt       n,i;
  const Vec      *inx, *outx;
  Vec            *auxV;
  Mat            A;

  PetscFunctionBegin;
  ierr = PCGetOperators(pc,&A,NULL);CHKERRQ(ierr);
  ierr = MatShellGetContext(A,(void**)&data);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(in,NULL,&inx);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(out,NULL,&outx);CHKERRQ(ierr);
  n = data->r_e - data->r_s;
  ierr = SlepcVecPoolGetVecs(data->d->auxV,n,&auxV);CHKERRQ(ierr);
  /* auxV <- in */
  for (i=0;i<n;i++) {
    ierr = VecCopy(inx[i],auxV[i]);CHKERRQ(ierr);
  }
  /* auxV <- auxV - u*(v'*auxV) */
  ierr = dvd_improvex_applytrans_proj(data->d,auxV,n);CHKERRQ(ierr);
  /* out <- K' * aux */
  for (i=0;i<n;i++) {
    ierr = PCApplyTranspose(data->old_pc,auxV[i],outx[i]);CHKERRQ(ierr);
  }
  ierr = SlepcVecPoolRestoreVecs(data->d->auxV,n,&auxV);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_dvd_jd(Mat A,Vec in,Vec out)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data;
  PetscInt       n;
  const Vec      *inx, *outx;
  PCSide         side;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,(void**)&data);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(in,NULL,&inx);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(out,NULL,&outx);CHKERRQ(ierr);
  n = data->r_e - data->r_s;
  /* out <- theta[1]A*in - theta[0]*B*in */
  ierr = dvd_aux_matmult(data,inx,outx);CHKERRQ(ierr);
  ierr = KSPGetPCSide(data->ksp,&side);CHKERRQ(ierr);
  if (side == PC_RIGHT) {
    /* out <- out - v*(u'*out) */
    ierr = dvd_improvex_apply_proj(data->d,(Vec*)outx,n);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMultTranspose_dvd_jd(Mat A,Vec in,Vec out)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data;
  PetscInt       n,i;
  const Vec      *inx,*outx,*r;
  Vec            *auxV;
  PCSide         side;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,(void**)&data);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(in,NULL,&inx);CHKERRQ(ierr);
  ierr = VecCompGetSubVecs(out,NULL,&outx);CHKERRQ(ierr);
  n = data->r_e - data->r_s;
  ierr = KSPGetPCSide(data->ksp,&side);CHKERRQ(ierr);
  if (side == PC_RIGHT) {
    /* auxV <- in */
    ierr = SlepcVecPoolGetVecs(data->d->auxV,n,&auxV);CHKERRQ(ierr);
    for (i=0;i<n;i++) {
      ierr = VecCopy(inx[i],auxV[i]);CHKERRQ(ierr);
    }
    /* auxV <- auxV - v*(u'*auxV) */
    ierr = dvd_improvex_applytrans_proj(data->d,auxV,n);CHKERRQ(ierr);
    r = auxV;
  } else r = inx;
  /* out <- theta[1]A*r - theta[0]*B*r */
  ierr = dvd_aux_matmulttrans(data,r,outx);CHKERRQ(ierr);
  if (side == PC_RIGHT) {
    ierr = SlepcVecPoolRestoreVecs(data->d->auxV,n,&auxV);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MatCreateVecs_dvd_jd(Mat A,Vec *right,Vec *left)
{
  PetscErrorCode ierr;
  Vec            *r,*l;
  dvdImprovex_jd *data;
  PetscInt       n,i;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,(void**)&data);CHKERRQ(ierr);
  n = data->ksp_max_size;
  if (right) {
    ierr = PetscMalloc1(n,&r);CHKERRQ(ierr);
  }
  if (left) {
    ierr = PetscMalloc1(n,&l);CHKERRQ(ierr);
  }
  for (i=0;i<n;i++) {
    ierr = MatCreateVecs(data->d->A,right?&r[i]:NULL,left?&l[i]:NULL);CHKERRQ(ierr);
  }
  if (right) {
    ierr = VecCreateCompWithVecs(r,n,data->friends,right);CHKERRQ(ierr);
    for (i=0;i<n;i++) {
      ierr = VecDestroy(&r[i]);CHKERRQ(ierr);
    }
  }
  if (left) {
    ierr = VecCreateCompWithVecs(l,n,data->friends,left);CHKERRQ(ierr);
    for (i=0;i<n;i++) {
      ierr = VecDestroy(&l[i]);CHKERRQ(ierr);
    }
  }

  if (right) {
    ierr = PetscFree(r);CHKERRQ(ierr);
  }
  if (left) {
    ierr = PetscFree(l);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_improvex_jd_start(dvdDashboard *d)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscInt       rA, cA, rlA, clA;
  Mat            A;
  PetscBool      t;
  PC             pc;
  Vec            v0[2];

  PetscFunctionBegin;
  data->size_cX = data->old_size_X = 0;
  data->lastTol = data->dynamic?0.5:0.0;

  /* Setup the ksp */
  if (data->ksp) {
    /* Create the reference vector */
    ierr = BVGetColumn(d->eps->V,0,&v0[0]);CHKERRQ(ierr);
    v0[1] = v0[0];
    ierr = VecCreateCompWithVecs(v0,data->ksp_max_size,NULL,&data->friends);CHKERRQ(ierr);
    ierr = BVRestoreColumn(d->eps->V,0,&v0[0]);CHKERRQ(ierr);
    ierr = PetscLogObjectParent((PetscObject)d->eps,(PetscObject)data->friends);CHKERRQ(ierr);

    /* Save the current pc and set a PCNONE */
    ierr = KSPGetPC(data->ksp, &data->old_pc);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)data->old_pc,PCNONE,&t);CHKERRQ(ierr);
    data->lastTol = 0.5;
    if (t) data->old_pc = 0;
    else {
      ierr = PetscObjectReference((PetscObject)data->old_pc);CHKERRQ(ierr);
      ierr = PCCreate(PetscObjectComm((PetscObject)d->eps),&pc);CHKERRQ(ierr);
      ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);
      ierr = PCSetOperators(pc,d->A,d->A);CHKERRQ(ierr);
      ierr = PCSetReusePreconditioner(pc,PETSC_TRUE);CHKERRQ(ierr);
      ierr = PCShellSetApply(pc,PCApply_dvd);CHKERRQ(ierr);
      ierr = PCShellSetApplyBA(pc,PCApplyBA_dvd);CHKERRQ(ierr);
      ierr = PCShellSetApplyTranspose(pc,PCApplyTranspose_dvd);CHKERRQ(ierr);
      ierr = KSPSetPC(data->ksp,pc);CHKERRQ(ierr);
      ierr = PCDestroy(&pc);CHKERRQ(ierr);
    }

    /* Create the (I-v*u')*K*(A-s*B) matrix */
    ierr = MatGetSize(d->A,&rA,&cA);CHKERRQ(ierr);
    ierr = MatGetLocalSize(d->A,&rlA,&clA);CHKERRQ(ierr);
    ierr = MatCreateShell(PetscObjectComm((PetscObject)d->A),rlA*data->ksp_max_size,clA*data->ksp_max_size,rA*data->ksp_max_size,cA*data->ksp_max_size,data,&A);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatMult_dvd_jd);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMultTranspose_dvd_jd);CHKERRQ(ierr);
    ierr = MatShellSetOperation(A,MATOP_CREATE_VECS,(void(*)(void))MatCreateVecs_dvd_jd);CHKERRQ(ierr);

    /* Try to avoid KSPReset */
    ierr = KSPGetOperatorsSet(data->ksp,&t,NULL);CHKERRQ(ierr);
    if (t) {
      Mat      M;
      PetscInt rM;
      ierr = KSPGetOperators(data->ksp,&M,NULL);CHKERRQ(ierr);
      ierr = MatGetSize(M,&rM,NULL);CHKERRQ(ierr);
      if (rM != rA*data->ksp_max_size) { ierr = KSPReset(data->ksp);CHKERRQ(ierr); }
    }
    ierr = KSPSetOperators(data->ksp,A,A);CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(data->ksp,PETSC_TRUE);CHKERRQ(ierr);
    ierr = KSPSetUp(data->ksp);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
  } else {
    data->old_pc = 0;
    data->friends = NULL;
  }
  ierr = BVSetActiveColumns(data->KZ,0,0);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(data->U,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Compute: u <- X, v <- K*(theta[0]*A+theta[1]*B)*X,
  kr <- K^{-1}*(A-eig*B)*X, being X <- V*pX[i_s..i_e-1], Y <- W*pY[i_s..i_e-1]
  where
  pX,pY, the right and left eigenvectors of the projected system
  ld, the leading dimension of pX and pY
*/
static PetscErrorCode dvd_improvex_jd_proj_cuv(dvdDashboard *d,PetscInt i_s,PetscInt i_e,Vec *kr,PetscScalar *theta,PetscScalar *thetai,PetscScalar *pX,PetscScalar *pY,PetscInt ld)
{
  PetscErrorCode ierr;
  PetscInt       n=i_e-i_s,size_KZ,V_new,rm,i,lv,kv,lKZ,kKZ;
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscScalar    *array;
  Mat            M;
  Vec            u[2],v[2];
  PetscBLASInt   s,ldXKZ,info;

  PetscFunctionBegin;
  /* Check consistency */
  ierr = BVGetActiveColumns(d->eps->V,&lv,&kv);CHKERRQ(ierr);
  V_new = lv - data->size_cX;
  if (V_new > data->old_size_X) SETERRQ(PETSC_COMM_SELF,1,"Consistency broken");
  data->old_size_X = n;
  data->size_cX = lv;

  /* KZ <- KZ(rm:rm+max_cX-1) */
  ierr = BVGetActiveColumns(data->KZ,&lKZ,&kKZ);CHKERRQ(ierr);
  rm = PetscMax(V_new+lKZ,0);
  if (rm > 0) {
    for (i=0;i<lKZ;i++) {
      ierr = BVCopyColumn(data->KZ,i+rm,i);CHKERRQ(ierr);
      ierr = BVCopyColumn(data->U,i+rm,i);CHKERRQ(ierr);
    }
  }

  /* XKZ <- XKZ(rm:rm+max_cX-1,rm:rm+max_cX-1) */
  if (rm > 0) {
    for (i=0;i<lKZ;i++) {
      ierr = PetscArraycpy(&data->XKZ[i*data->ldXKZ+i],&data->XKZ[(i+rm)*data->ldXKZ+i+rm],lKZ);CHKERRQ(ierr);
    }
  }
  lKZ = PetscMin(0,lKZ+V_new);
  ierr = BVSetActiveColumns(data->KZ,lKZ,lKZ+n);CHKERRQ(ierr);
  ierr = BVSetActiveColumns(data->U,lKZ,lKZ+n);CHKERRQ(ierr);

  /* Compute X, KZ and KR */
  ierr = BVGetColumn(data->U,lKZ,u);CHKERRQ(ierr);
  if (n>1) { ierr = BVGetColumn(data->U,lKZ+1,&u[1]);CHKERRQ(ierr); }
  ierr = BVGetColumn(data->KZ,lKZ,v);CHKERRQ(ierr);
  if (n>1) { ierr = BVGetColumn(data->KZ,lKZ+1,&v[1]);CHKERRQ(ierr); }
  ierr = d->improvex_jd_proj_uv(d,i_s,i_e,u,v,kr,theta,thetai,pX,pY,ld);CHKERRQ(ierr);
  ierr = BVRestoreColumn(data->U,lKZ,u);CHKERRQ(ierr);
  if (n>1) { ierr = BVRestoreColumn(data->U,lKZ+1,&u[1]);CHKERRQ(ierr); }
  ierr = BVRestoreColumn(data->KZ,lKZ,v);CHKERRQ(ierr);
  if (n>1) { ierr = BVRestoreColumn(data->KZ,lKZ+1,&v[1]);CHKERRQ(ierr); }

  /* XKZ <- U'*KZ */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,lKZ+n,lKZ+n,NULL,&M);CHKERRQ(ierr);
  ierr = BVMatProject(data->KZ,NULL,data->U,M);CHKERRQ(ierr);
  ierr = MatDenseGetArray(M,&array);CHKERRQ(ierr);
  for (i=lKZ;i<lKZ+n;i++) { /* upper part */
    ierr = PetscArraycpy(&data->XKZ[data->ldXKZ*i],&array[i*(lKZ+n)],lKZ);CHKERRQ(ierr);
  }
  for (i=0;i<lKZ+n;i++) { /* lower part */
    ierr = PetscArraycpy(&data->XKZ[data->ldXKZ*i+lKZ],&array[i*(lKZ+n)+lKZ],n);CHKERRQ(ierr);
  }
  ierr = MatDenseRestoreArray(M,&array);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);

  /* iXKZ <- inv(XKZ) */
  size_KZ = lKZ+n;
  ierr = PetscBLASIntCast(lKZ+n,&s);CHKERRQ(ierr);
  data->ldiXKZ = data->size_iXKZ = size_KZ;
  for (i=0;i<size_KZ;i++) {
    ierr = PetscArraycpy(&data->iXKZ[data->ldiXKZ*i],&data->XKZ[data->ldXKZ*i],size_KZ);CHKERRQ(ierr);
  }
  ierr = PetscBLASIntCast(data->ldiXKZ,&ldXKZ);CHKERRQ(ierr);
  ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);CHKERRQ(ierr);
  PetscStackCallBLAS("LAPACKgetrf",LAPACKgetrf_(&s,&s,data->iXKZ,&ldXKZ,data->iXKZPivots,&info));
  ierr = PetscFPTrapPop();CHKERRQ(ierr);
  SlepcCheckLapackInfo("getrf",info);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_improvex_jd_gen(dvdDashboard *d,PetscInt r_s,PetscInt r_e,PetscInt *size_D)
{
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscErrorCode ierr;
  PetscInt       i,j,n,maxits,maxits0,lits,s,ld,k,max_size_D,lV,kV;
  PetscScalar    *pX,*pY;
  PetscReal      tol,tol0;
  Vec            *kr,kr_comp,D_comp,D[2],kr0[2];
  PetscBool      odd_situation = PETSC_FALSE;

  PetscFunctionBegin;
  ierr = BVGetActiveColumns(d->eps->V,&lV,&kV);CHKERRQ(ierr);
  max_size_D = d->eps->ncv-kV;
  /* Quick exit */
  if ((max_size_D == 0) || r_e-r_s <= 0) {
   *size_D = 0;
    PetscFunctionReturn(0);
  }

  n = PetscMin(PetscMin(data->size_X, max_size_D), r_e-r_s);
  if (n == 0) SETERRQ(PETSC_COMM_SELF,1,"n == 0");
  if (data->size_X < r_e-r_s) SETERRQ(PETSC_COMM_SELF,1,"size_X < r_e-r_s");

  ierr = DSGetLeadingDimension(d->eps->ds,&ld);CHKERRQ(ierr);

  /* Restart lastTol if a new pair converged */
  if (data->dynamic && data->size_cX < lV)
    data->lastTol = 0.5;

  for (i=0,s=0;i<n;i+=s) {
    /* If the selected eigenvalue is complex, but the arithmetic is real... */
#if !defined(PETSC_USE_COMPLEX)
    if (d->eigi[r_s+i] != 0.0) {
      if (i+2 <= max_size_D) s=2;
      else break;
    } else
#endif
      s=1;

    data->r_s = r_s+i;
    data->r_e = r_s+i+s;
    ierr = SlepcVecPoolGetVecs(d->auxV,s,&kr);CHKERRQ(ierr);

    /* Compute theta, maximum iterations and tolerance */
    maxits = 0;
    tol = 1;
    for (j=0;j<s;j++) {
      ierr = d->improvex_jd_lit(d,r_s+i+j,&data->theta[2*j],&data->thetai[j],&maxits0,&tol0);CHKERRQ(ierr);
      maxits += maxits0;
      tol *= tol0;
    }
    maxits/= s;
    tol = data->dynamic?data->lastTol:PetscExpReal(PetscLogReal(tol)/s);

    /* Compute u, v and kr */
    k = r_s+i;
    ierr = DSVectors(d->eps->ds,DS_MAT_X,&k,NULL);CHKERRQ(ierr);
    k = r_s+i;
    ierr = DSVectors(d->eps->ds,DS_MAT_Y,&k,NULL);CHKERRQ(ierr);
    ierr = DSGetArray(d->eps->ds,DS_MAT_X,&pX);CHKERRQ(ierr);
    ierr = DSGetArray(d->eps->ds,DS_MAT_Y,&pY);CHKERRQ(ierr);
    ierr = dvd_improvex_jd_proj_cuv(d,r_s+i,r_s+i+s,kr,data->theta,data->thetai,pX,pY,ld);CHKERRQ(ierr);
    ierr = DSRestoreArray(d->eps->ds,DS_MAT_X,&pX);CHKERRQ(ierr);
    ierr = DSRestoreArray(d->eps->ds,DS_MAT_Y,&pY);CHKERRQ(ierr);

    /* Check if the first eigenpairs are converged */
    if (i == 0) {
      PetscInt oldnpreconv = d->npreconv;
      ierr = d->preTestConv(d,0,r_s+s,r_s+s,&d->npreconv);CHKERRQ(ierr);
      if (d->npreconv > oldnpreconv) break;
    }

    /* Test the odd situation of solving Ax=b with A=I */
#if !defined(PETSC_USE_COMPLEX)
    odd_situation = (data->ksp && data->theta[0] == 1. && data->theta[1] == 0. && data->thetai[0] == 0. && d->B == NULL)? PETSC_TRUE: PETSC_FALSE;
#else
    odd_situation = (data->ksp && data->theta[0] == 1. && data->theta[1] == 0. && d->B == NULL)? PETSC_TRUE: PETSC_FALSE;
#endif
    /* If JD */
    if (data->ksp && !odd_situation) {
      /* kr <- -kr */
      for (j=0;j<s;j++) {
        ierr = VecScale(kr[j],-1.0);CHKERRQ(ierr);
      }

      /* Compose kr and D */
      kr0[0] = kr[0];
      kr0[1] = (s==2 ? kr[1] : NULL);
      ierr = VecCreateCompWithVecs(kr0,data->ksp_max_size,data->friends,&kr_comp);CHKERRQ(ierr);
      ierr = BVGetColumn(d->eps->V,kV+i,&D[0]);CHKERRQ(ierr);
      if (s==2) { ierr = BVGetColumn(d->eps->V,kV+i+1,&D[1]);CHKERRQ(ierr); }
      else D[1] = NULL;
      ierr = VecCreateCompWithVecs(D,data->ksp_max_size,data->friends,&D_comp);CHKERRQ(ierr);
      ierr = VecCompSetSubVecs(data->friends,s,NULL);CHKERRQ(ierr);

      /* Solve the correction equation */
      ierr = KSPSetTolerances(data->ksp,tol,PETSC_DEFAULT,PETSC_DEFAULT,maxits);CHKERRQ(ierr);
      ierr = KSPSolve(data->ksp,kr_comp,D_comp);CHKERRQ(ierr);
      ierr = KSPGetIterationNumber(data->ksp,&lits);CHKERRQ(ierr);

      /* Destroy the composed ks and D */
      ierr = VecDestroy(&kr_comp);CHKERRQ(ierr);
      ierr = VecDestroy(&D_comp);CHKERRQ(ierr);
      ierr = BVRestoreColumn(d->eps->V,kV+i,&D[0]);CHKERRQ(ierr);
      if (s==2) { ierr = BVRestoreColumn(d->eps->V,kV+i+1,&D[1]);CHKERRQ(ierr); }

    /* If GD */
    } else {
      ierr = BVGetColumn(d->eps->V,kV+i,&D[0]);CHKERRQ(ierr);
      if (s==2) { ierr = BVGetColumn(d->eps->V,kV+i+1,&D[1]);CHKERRQ(ierr); }
      for (j=0;j<s;j++) {
        ierr = d->improvex_precond(d,r_s+i+j,kr[j],D[j]);CHKERRQ(ierr);
      }
      ierr = dvd_improvex_apply_proj(d,D,s);CHKERRQ(ierr);
      ierr = BVRestoreColumn(d->eps->V,kV+i,&D[0]);CHKERRQ(ierr);
      if (s==2) { ierr = BVRestoreColumn(d->eps->V,kV+i+1,&D[1]);CHKERRQ(ierr); }
    }
    /* Prevent that short vectors are discarded in the orthogonalization */
    if (i == 0 && d->eps->errest[d->nconv+r_s] > PETSC_MACHINE_EPSILON && d->eps->errest[d->nconv+r_s] < PETSC_MAX_REAL) {
      for (j=0;j<s;j++) {
        ierr = BVScaleColumn(d->eps->V,kV+i+j,1.0/d->eps->errest[d->nconv+r_s]);CHKERRQ(ierr);
      }
    }
    ierr = SlepcVecPoolRestoreVecs(d->auxV,s,&kr);CHKERRQ(ierr);
  }
  *size_D = i;
  if (data->dynamic) data->lastTol = PetscMax(data->lastTol/2.0,PETSC_MACHINE_EPSILON*10.0);
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_improvex_jd(dvdDashboard *d,dvdBlackboard *b,KSP ksp,PetscInt max_bs,PetscBool dynamic)
{
  PetscErrorCode ierr;
  dvdImprovex_jd *data;
  PetscBool      useGD;
  PC             pc;
  PetscInt       size_P;

  PetscFunctionBegin;
  /* Setting configuration constrains */
  ierr = PetscObjectTypeCompare((PetscObject)ksp,KSPPREONLY,&useGD);CHKERRQ(ierr);

  /* If the arithmetic is real and the problem is not Hermitian, then
     the block size is incremented in one */
#if !defined(PETSC_USE_COMPLEX)
  if (!DVD_IS(d->sEP,DVD_EP_HERMITIAN)) {
    max_bs++;
    b->max_size_P = PetscMax(b->max_size_P,2);
  } else
#endif
  {
    b->max_size_P = PetscMax(b->max_size_P,1);
  }
  b->max_size_X = PetscMax(b->max_size_X,max_bs);
  size_P = b->max_size_P;

  /* Setup the preconditioner */
  if (ksp) {
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = dvd_static_precond_PC(d,b,pc);CHKERRQ(ierr);
  } else {
    ierr = dvd_static_precond_PC(d,b,0);CHKERRQ(ierr);
  }

  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    ierr = PetscNewLog(d->eps,&data);CHKERRQ(ierr);
    data->dynamic = dynamic;
    ierr = PetscMalloc1(size_P*size_P,&data->XKZ);CHKERRQ(ierr);
    ierr = PetscMalloc1(size_P*size_P,&data->iXKZ);CHKERRQ(ierr);
    ierr = PetscMalloc1(size_P,&data->iXKZPivots);CHKERRQ(ierr);
    data->ldXKZ = size_P;
    data->size_X = b->max_size_X;
    d->improveX_data = data;
    data->ksp = useGD? NULL: ksp;
    data->d = d;
    d->improveX = dvd_improvex_jd_gen;
#if !defined(PETSC_USE_COMPLEX)
    if (!DVD_IS(d->sEP,DVD_EP_HERMITIAN)) data->ksp_max_size = 2;
    else
#endif
      data->ksp_max_size = 1;
    /* Create various vector basis */
    ierr = BVDuplicateResize(d->eps->V,size_P,&data->KZ);CHKERRQ(ierr);
    ierr = BVSetMatrix(data->KZ,NULL,PETSC_FALSE);CHKERRQ(ierr);
    ierr = BVDuplicate(data->KZ,&data->U);CHKERRQ(ierr);

    ierr = EPSDavidsonFLAdd(&d->startList,dvd_improvex_jd_start);CHKERRQ(ierr);
    ierr = EPSDavidsonFLAdd(&d->endList,dvd_improvex_jd_end);CHKERRQ(ierr);
    ierr = EPSDavidsonFLAdd(&d->destroyList,dvd_improvex_jd_d);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#if !defined(PETSC_USE_COMPLEX)
PETSC_STATIC_INLINE PetscErrorCode dvd_complex_rayleigh_quotient(Vec ur,Vec ui,Vec Axr,Vec Axi,Vec Bxr,Vec Bxi,PetscScalar *eigr,PetscScalar *eigi)
{
  PetscErrorCode ierr;
  PetscScalar    rAr,iAr,rAi,iAi,rBr,iBr,rBi,iBi,b0,b2,b4,b6,b7;

  PetscFunctionBegin;
  /* eigr = [(rAr+iAi)*(rBr+iBi) + (rAi-iAr)*(rBi-iBr)]/k
     eigi = [(rAi-iAr)*(rBr+iBi) - (rAr+iAi)*(rBi-iBr)]/k
     k    =  (rBr+iBi)*(rBr+iBi) + (rBi-iBr)*(rBi-iBr)    */
  ierr = VecDotBegin(Axr,ur,&rAr);CHKERRQ(ierr); /* r*A*r */
  ierr = VecDotBegin(Axr,ui,&iAr);CHKERRQ(ierr); /* i*A*r */
  ierr = VecDotBegin(Axi,ur,&rAi);CHKERRQ(ierr); /* r*A*i */
  ierr = VecDotBegin(Axi,ui,&iAi);CHKERRQ(ierr); /* i*A*i */
  ierr = VecDotBegin(Bxr,ur,&rBr);CHKERRQ(ierr); /* r*B*r */
  ierr = VecDotBegin(Bxr,ui,&iBr);CHKERRQ(ierr); /* i*B*r */
  ierr = VecDotBegin(Bxi,ur,&rBi);CHKERRQ(ierr); /* r*B*i */
  ierr = VecDotBegin(Bxi,ui,&iBi);CHKERRQ(ierr); /* i*B*i */
  ierr = VecDotEnd(Axr,ur,&rAr);CHKERRQ(ierr); /* r*A*r */
  ierr = VecDotEnd(Axr,ui,&iAr);CHKERRQ(ierr); /* i*A*r */
  ierr = VecDotEnd(Axi,ur,&rAi);CHKERRQ(ierr); /* r*A*i */
  ierr = VecDotEnd(Axi,ui,&iAi);CHKERRQ(ierr); /* i*A*i */
  ierr = VecDotEnd(Bxr,ur,&rBr);CHKERRQ(ierr); /* r*B*r */
  ierr = VecDotEnd(Bxr,ui,&iBr);CHKERRQ(ierr); /* i*B*r */
  ierr = VecDotEnd(Bxi,ur,&rBi);CHKERRQ(ierr); /* r*B*i */
  ierr = VecDotEnd(Bxi,ui,&iBi);CHKERRQ(ierr); /* i*B*i */
  b0 = rAr+iAi; /* rAr+iAi */
  b2 = rAi-iAr; /* rAi-iAr */
  b4 = rBr+iBi; /* rBr+iBi */
  b6 = rBi-iBr; /* rBi-iBr */
  b7 = b4*b4 + b6*b6; /* k */
  *eigr = (b0*b4 + b2*b6) / b7; /* eig_r */
  *eigi = (b2*b4 - b0*b6) / b7; /* eig_i */
  PetscFunctionReturn(0);
}
#endif

PETSC_STATIC_INLINE PetscErrorCode dvd_compute_n_rr(PetscInt i_s,PetscInt n,PetscScalar *eigr,PetscScalar *eigi,Vec *u,Vec *Ax,Vec *Bx)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    b0,b1;

  PetscFunctionBegin;
  for (i=0; i<n; i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (eigi[i_s+i] != 0.0) {
      PetscScalar eigr0=0.0,eigi0=0.0;
      ierr = dvd_complex_rayleigh_quotient(u[i],u[i+1],Ax[i],Ax[i+1],Bx[i],Bx[i+1],&eigr0,&eigi0);CHKERRQ(ierr);
      if (PetscAbsScalar(eigr[i_s+i]-eigr0)/PetscAbsScalar(eigr[i_s+i]) > 1e-10 || PetscAbsScalar(eigi[i_s+i]-eigi0)/PetscAbsScalar(eigi[i_s+i]) > 1e-10) {
        ierr = PetscInfo4(u[0],"The eigenvalue %g%+gi is far from its Rayleigh quotient value %g%+gi\n",(double)eigr[i_s+i],(double)eigi[i_s+i],(double)eigr0,(double)eigi0);CHKERRQ(ierr);
      }
      i++;
    } else
#endif
    {
      ierr = VecDotBegin(Ax[i],u[i],&b0);CHKERRQ(ierr);
      ierr = VecDotBegin(Bx[i],u[i],&b1);CHKERRQ(ierr);
      ierr = VecDotEnd(Ax[i],u[i],&b0);CHKERRQ(ierr);
      ierr = VecDotEnd(Bx[i],u[i],&b1);CHKERRQ(ierr);
      b0 = b0/b1;
      if (PetscAbsScalar(eigr[i_s+i]-b0)/PetscAbsScalar(eigr[i_s+i]) > 1e-10) {
        ierr = PetscInfo4(u[0],"The eigenvalue %g+%g is far from its Rayleigh quotient value %g+%g\n",(double)PetscRealPart(eigr[i_s+i]),(double)PetscImaginaryPart(eigr[i_s+i]),(double)PetscRealPart(b0),(double)PetscImaginaryPart(b0));CHKERRQ(ierr);
      }
    }
  }
  PetscFunctionReturn(0);
}

/*
  Compute: u <- X, v <- K*(theta[0]*A+theta[1]*B)*X,
  kr <- K^{-1}*(A-eig*B)*X, being X <- V*pX[i_s..i_e-1], Y <- W*pY[i_s..i_e-1]
  where
  pX,pY, the right and left eigenvectors of the projected system
  ld, the leading dimension of pX and pY
*/
PetscErrorCode dvd_improvex_jd_proj_uv_KZX(dvdDashboard *d,PetscInt i_s,PetscInt i_e,Vec *u,Vec *v,Vec *kr,PetscScalar *theta,PetscScalar *thetai,PetscScalar *pX,PetscScalar *pY,PetscInt ld)
{
  PetscErrorCode ierr;
  PetscInt       n = i_e-i_s,i;
  PetscScalar    *b;
  Vec            *Ax,*Bx,*r;
  Mat            M;
  BV             X;

  PetscFunctionBegin;
  ierr = BVDuplicateResize(d->eps->V,4,&X);CHKERRQ(ierr);
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,4,4,NULL,&M);CHKERRQ(ierr);
  /* u <- X(i) */
  ierr = dvd_improvex_compute_X(d,i_s,i_e,u,pX,ld);CHKERRQ(ierr);

  /* v <- theta[0]A*u + theta[1]*B*u */

  /* Bx <- B*X(i) */
  Bx = kr;
  if (d->BX) {
    for (i=i_s; i<i_e; ++i) {
      ierr = BVMultVec(d->BX,1.0,0.0,Bx[i-i_s],&pX[ld*i]);CHKERRQ(ierr);
    }
  } else {
    for (i=0;i<n;i++) {
      if (d->B) {
        ierr = MatMult(d->B, u[i], Bx[i]);CHKERRQ(ierr);
      } else {
        ierr = VecCopy(u[i], Bx[i]);CHKERRQ(ierr);
      }
    }
  }

  /* Ax <- A*X(i) */
  ierr = SlepcVecPoolGetVecs(d->auxV,i_e-i_s,&r);CHKERRQ(ierr);
  Ax = r;
  for (i=i_s; i<i_e; ++i) {
    ierr = BVMultVec(d->AX,1.0,0.0,Ax[i-i_s],&pX[ld*i]);CHKERRQ(ierr);
  }

  /* v <- Y(i) */
  for (i=i_s; i<i_e; ++i) {
    ierr = BVMultVec(d->W?d->W:d->eps->V,1.0,0.0,v[i-i_s],&pY[ld*i]);CHKERRQ(ierr);
  }

  /* Recompute the eigenvalue */
  ierr = dvd_compute_n_rr(i_s,n,d->eigr,d->eigi,v,Ax,Bx);CHKERRQ(ierr);

  for (i=0;i<n;i++) {
#if !defined(PETSC_USE_COMPLEX)
    if (d->eigi[i_s+i] != 0.0) {
      /* [r_i r_i+1 kr_i kr_i+1]*= [ theta_2i'    0            1        0
                                       0         theta_2i'     0        1
                                     theta_2i+1 -thetai_i   -eigr_i -eigi_i
                                     thetai_i    theta_2i+1  eigi_i -eigr_i ] */
      ierr = MatDenseGetArray(M,&b);CHKERRQ(ierr);
      b[0] = b[5] = PetscConj(theta[2*i]);
      b[2] = b[7] = -theta[2*i+1];
      b[6] = -(b[3] = thetai[i]);
      b[1] = b[4] = 0.0;
      b[8] = b[13] = 1.0/d->nX[i_s+i];
      b[10] = b[15] = -d->eigr[i_s+i]/d->nX[i_s+i];
      b[14] = -(b[11] = d->eigi[i_s+i]/d->nX[i_s+i]);
      b[9] = b[12] = 0.0;
      ierr = MatDenseRestoreArray(M,&b);CHKERRQ(ierr);
      ierr = BVInsertVec(X,0,Ax[i]);CHKERRQ(ierr);
      ierr = BVInsertVec(X,1,Ax[i+1]);CHKERRQ(ierr);
      ierr = BVInsertVec(X,2,Bx[i]);CHKERRQ(ierr);
      ierr = BVInsertVec(X,3,Bx[i+1]);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(X,0,4);CHKERRQ(ierr);
      ierr = BVMultInPlace(X,M,0,4);CHKERRQ(ierr);
      ierr = BVCopyVec(X,0,Ax[i]);CHKERRQ(ierr);
      ierr = BVCopyVec(X,1,Ax[i+1]);CHKERRQ(ierr);
      ierr = BVCopyVec(X,2,Bx[i]);CHKERRQ(ierr);
      ierr = BVCopyVec(X,3,Bx[i+1]);CHKERRQ(ierr);
      i++;
    } else
#endif
    {
      /* [Ax_i Bx_i]*= [ theta_2i'    1/nX_i
                        theta_2i+1  -eig_i/nX_i ] */
      ierr = MatDenseGetArray(M,&b);CHKERRQ(ierr);
      b[0] = PetscConj(theta[i*2]);
      b[1] = theta[i*2+1];
      b[4] = 1.0/d->nX[i_s+i];
      b[5] = -d->eigr[i_s+i]/d->nX[i_s+i];
      ierr = MatDenseRestoreArray(M,&b);CHKERRQ(ierr);
      ierr = BVInsertVec(X,0,Ax[i]);CHKERRQ(ierr);
      ierr = BVInsertVec(X,1,Bx[i]);CHKERRQ(ierr);
      ierr = BVSetActiveColumns(X,0,2);CHKERRQ(ierr);
      ierr = BVMultInPlace(X,M,0,2);CHKERRQ(ierr);
      ierr = BVCopyVec(X,0,Ax[i]);CHKERRQ(ierr);
      ierr = BVCopyVec(X,1,Bx[i]);CHKERRQ(ierr);
    }
  }
  for (i=0; i<n; i++) d->nX[i_s+i] = 1.0;

  /* v <- K^{-1} r = K^{-1}(theta_2i'*Ax + theta_2i+1*Bx) */
  for (i=0;i<n;i++) {
    ierr = d->improvex_precond(d,i_s+i,r[i],v[i]);CHKERRQ(ierr);
  }
  ierr = SlepcVecPoolRestoreVecs(d->auxV,i_e-i_s,&r);CHKERRQ(ierr);

  /* kr <- P*(Ax - eig_i*Bx) */
  ierr = d->calcpairs_proj_res(d,i_s,i_e,kr);CHKERRQ(ierr);
  ierr = BVDestroy(&X);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_improvex_jd_lit_const_0(dvdDashboard *d,PetscInt i,PetscScalar* theta,PetscScalar* thetai,PetscInt *maxits,PetscReal *tol)
{
  dvdImprovex_jd *data = (dvdImprovex_jd*)d->improveX_data;
  PetscReal      a;

  PetscFunctionBegin;
  a = SlepcAbsEigenvalue(d->eigr[i],d->eigi[i]);

  if (d->nR[i] < data->fix*a) {
    theta[0] = d->eigr[i];
    theta[1] = 1.0;
#if !defined(PETSC_USE_COMPLEX)
    *thetai = d->eigi[i];
#endif
  } else {
    theta[0] = d->target[0];
    theta[1] = d->target[1];
#if !defined(PETSC_USE_COMPLEX)
    *thetai = 0.0;
#endif
}

#if defined(PETSC_USE_COMPLEX)
  if (thetai) *thetai = 0.0;
#endif
  *maxits = data->maxits;
  *tol = data->tol;
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_improvex_jd_lit_const(dvdDashboard *d,dvdBlackboard *b,PetscInt maxits,PetscReal tol,PetscReal fix)
{
  dvdImprovex_jd  *data = (dvdImprovex_jd*)d->improveX_data;

  PetscFunctionBegin;
  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    data->maxits = maxits;
    data->tol = tol;
    data->fix = fix;
    d->improvex_jd_lit = dvd_improvex_jd_lit_const_0;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_improvex_jd_proj_uv(dvdDashboard *d,dvdBlackboard *b)
{
  PetscFunctionBegin;
  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    d->improvex_jd_proj_uv = dvd_improvex_jd_proj_uv_KZX;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_improvex_compute_X(dvdDashboard *d,PetscInt i_s,PetscInt i_e,Vec *u_,PetscScalar *pX,PetscInt ld)
{
  PetscErrorCode ierr;
  PetscInt       n = i_e - i_s,i;
  Vec            *u;

  PetscFunctionBegin;
  if (u_) u = u_;
  else if (d->correctXnorm) {
    ierr = SlepcVecPoolGetVecs(d->auxV,i_e-i_s,&u);CHKERRQ(ierr);
  }
  if (u_ || d->correctXnorm) {
    for (i=0; i<n; i++) {
      ierr = BVMultVec(d->eps->V,1.0,0.0,u[i],&pX[ld*(i+i_s)]);CHKERRQ(ierr);
    }
  }
  /* nX(i) <- ||X(i)|| */
  if (d->correctXnorm) {
    for (i=0; i<n; i++) {
      ierr = VecNormBegin(u[i],NORM_2,&d->nX[i_s+i]);CHKERRQ(ierr);
    }
    for (i=0; i<n; i++) {
      ierr = VecNormEnd(u[i],NORM_2,&d->nX[i_s+i]);CHKERRQ(ierr);
    }
#if !defined(PETSC_USE_COMPLEX)
    for (i=0;i<n;i++) {
      if (d->eigi[i_s+i] != 0.0) {
        d->nX[i_s+i] = d->nX[i_s+i+1] = PetscSqrtScalar(d->nX[i_s+i]*d->nX[i_s+i]+d->nX[i_s+i+1]*d->nX[i_s+i+1]);
        i++;
      }
    }
#endif
  } else {
    for (i=0;i<n;i++) d->nX[i_s+i] = 1.0;
  }
  if (d->correctXnorm && !u_) {
    ierr = SlepcVecPoolRestoreVecs(d->auxV,i_e-i_s,&u);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

