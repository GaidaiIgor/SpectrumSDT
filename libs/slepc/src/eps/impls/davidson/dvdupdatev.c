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

   Step: test for restarting, updateV, restartV
*/

#include "davidson.h"

typedef struct {
  PetscInt          min_size_V;        /* restart with this number of eigenvectors */
  PetscInt          plusk;             /* at restart, save plusk vectors from last iteration */
  PetscInt          mpd;               /* max size of the searching subspace */
  void              *old_updateV_data; /* old updateV data */
  PetscErrorCode    (*old_isRestarting)(dvdDashboard*,PetscBool*);  /* old isRestarting */
  Mat               oldU;              /* previous projected right igenvectors */
  Mat               oldV;              /* previous projected left eigenvectors */
  PetscInt          size_oldU;         /* size of oldU */
  PetscBool         allResiduals;      /* if computing all the residuals */
} dvdManagV_basic;

static PetscErrorCode dvd_updateV_start(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        i;

  PetscFunctionBegin;
  for (i=0;i<d->eps->ncv;i++) d->eigi[i] = 0.0;
  d->nR = d->real_nR;
  for (i=0;i<d->eps->ncv;i++) d->nR[i] = 1.0;
  d->nX = d->real_nX;
  for (i=0;i<d->eps->ncv;i++) d->errest[i] = 1.0;
  data->size_oldU = 0;
  d->nconv = 0;
  d->npreconv = 0;
  d->V_tra_s = d->V_tra_e = d->V_new_s = d->V_new_e = 0;
  d->size_D = 0;
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_isrestarting_fullV(dvdDashboard *d,PetscBool *r)
{
  PetscErrorCode  ierr;
  PetscInt        l,k;
  PetscBool       restart;
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;

  PetscFunctionBegin;
  ierr = BVGetActiveColumns(d->eps->V,&l,&k);CHKERRQ(ierr);
  restart = (k+2 > d->eps->ncv)? PETSC_TRUE: PETSC_FALSE;

  /* Check old isRestarting function */
  if (!restart && data->old_isRestarting) {
    ierr = data->old_isRestarting(d,&restart);CHKERRQ(ierr);
  }
  *r = restart;
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_managementV_basic_d(dvdDashboard *d)
{
  PetscErrorCode  ierr;
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;

  PetscFunctionBegin;
  /* Restore changes in dvdDashboard */
  d->updateV_data = data->old_updateV_data;

  /* Free local data */
  ierr = MatDestroy(&data->oldU);CHKERRQ(ierr);
  ierr = MatDestroy(&data->oldV);CHKERRQ(ierr);
  ierr = PetscFree(d->real_nR);CHKERRQ(ierr);
  ierr = PetscFree(d->real_nX);CHKERRQ(ierr);
  ierr = PetscFree(data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_updateV_conv_gen(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        npreconv,cMT,cMTX,lV,kV,nV;
  PetscErrorCode  ierr;
  Mat             Z;
  PetscBool       t;
#if !defined(PETSC_USE_COMPLEX)
  PetscInt        i;
#endif

  PetscFunctionBegin;
  npreconv = d->npreconv;
  /* Constrains the converged pairs to nev */
#if !defined(PETSC_USE_COMPLEX)
  /* Tries to maintain together conjugate eigenpairs */
  for (i=0; (i + (d->eigi[i]!=0.0?1:0) < npreconv) && (d->nconv + i < d->nev); i+= (d->eigi[i]!=0.0?2:1));
  npreconv = i;
#else
  npreconv = PetscMax(PetscMin(d->nev-d->nconv,npreconv),0);
#endif
  /* For GHEP without B-ortho, converge all of the requested pairs at once */
  ierr = PetscObjectTypeCompare((PetscObject)d->eps->ds,DSGHEP,&t);CHKERRQ(ierr);
  if (t && d->nconv+npreconv<d->nev) npreconv = 0;
  /* Quick exit */
  if (npreconv == 0) PetscFunctionReturn(0);

  ierr = BVGetActiveColumns(d->eps->V,&lV,&kV);CHKERRQ(ierr);
  nV  = kV - lV;
  cMT = nV - npreconv;
  /* Harmonics restarts wiht right eigenvectors, and other with the left ones.
     If the problem is standard or hermitian, left and right vectors are the same */
  if (!(d->W||DVD_IS(d->sEP,DVD_EP_STD)||DVD_IS(d->sEP,DVD_EP_HERMITIAN))) {
    /* ps.Q <- [ps.Q(0:npreconv-1) ps.Z(npreconv:size_H-1)] */
    ierr = DSGetMat(d->eps->ds,DS_MAT_Z,&Z);CHKERRQ(ierr);
    ierr = DSCopyMat(d->eps->ds,DS_MAT_Q,0,npreconv,Z,0,npreconv,nV,cMT,PETSC_FALSE);CHKERRQ(ierr);
    ierr = MatDestroy(&Z);CHKERRQ(ierr);
  }
  if (DVD_IS(d->sEP,DVD_EP_INDEFINITE)) {
    ierr = DSPseudoOrthogonalize(d->eps->ds,DS_MAT_Q,nV,d->nBds,&cMTX,d->nBds);CHKERRQ(ierr);
  } else {
    ierr = DSOrthogonalize(d->eps->ds,DS_MAT_Q,nV,&cMTX);CHKERRQ(ierr);
  }
  cMT = cMTX - npreconv;

  if (d->W) {
    ierr = DSOrthogonalize(d->eps->ds,DS_MAT_Z,nV,&cMTX);CHKERRQ(ierr);
    cMT = PetscMin(cMT,cMTX - npreconv);
  }

  /* Lock the converged pairs */
  d->eigr+= npreconv;
#if !defined(PETSC_USE_COMPLEX)
  if (d->eigi) d->eigi+= npreconv;
#endif
  d->nconv+= npreconv;
  d->errest+= npreconv;
  /* Notify the changes in V and update the other subspaces */
  d->V_tra_s = npreconv;          d->V_tra_e = nV;
  d->V_new_s = cMT;               d->V_new_e = d->V_new_s;
  /* Remove oldU */
  data->size_oldU = 0;

  d->npreconv-= npreconv;
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_updateV_restart_gen(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        lV,kV,nV,size_plusk,size_X,cMTX,cMTY,max_restart_size;
  Mat             Z;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  /* Select size_X desired pairs from V */
  /* The restarted basis should:
     - have at least one spot to add a new direction;
     - keep converged vectors, npreconv;
     - keep at least 1 oldU direction if possible.
  */
  ierr = BVGetActiveColumns(d->eps->V,&lV,&kV);CHKERRQ(ierr);
  nV = kV - lV;
  max_restart_size = PetscMax(0,PetscMin(d->eps->mpd - 1,d->eps->ncv - lV - 2));
  size_X = PetscMin(PetscMin(data->min_size_V+d->npreconv,max_restart_size - (max_restart_size - d->npreconv > 1 && data->plusk > 0 && data->size_oldU > 0 ? 1 : 0)), nV);

  /* Add plusk eigenvectors from the previous iteration */
  size_plusk = PetscMax(0,PetscMin(PetscMin(PetscMin(data->plusk,data->size_oldU),max_restart_size - size_X),nV - size_X));

  d->size_MT = nV;
  /* ps.Q <- orth([pX(0:size_X-1) [oldU(0:size_plusk-1); 0] ]) */
  /* Harmonics restarts wiht right eigenvectors, and other with the left ones.
     If the problem is standard or hermitian, left and right vectors are the same */
  if (!(d->W||DVD_IS(d->sEP,DVD_EP_STD)||DVD_IS(d->sEP,DVD_EP_HERMITIAN))) {
    ierr = DSGetMat(d->eps->ds,DS_MAT_Z,&Z);CHKERRQ(ierr);
    ierr = DSCopyMat(d->eps->ds,DS_MAT_Q,0,0,Z,0,0,nV,size_X,PETSC_FALSE);CHKERRQ(ierr);
    ierr = MatDestroy(&Z);CHKERRQ(ierr);
  }
  if (size_plusk > 0 && DVD_IS(d->sEP,DVD_EP_INDEFINITE)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unsupported plusk>0 in indefinite eigenvalue problems");
  if (size_plusk > 0) {
    ierr = DSCopyMat(d->eps->ds,DS_MAT_Q,0,size_X,data->oldU,0,0,nV,size_plusk,PETSC_FALSE);CHKERRQ(ierr);
  }
  if (DVD_IS(d->sEP,DVD_EP_INDEFINITE)) {
    ierr = DSPseudoOrthogonalize(d->eps->ds,DS_MAT_Q,size_X,d->nBds,&cMTX,d->nBds);CHKERRQ(ierr);
  } else {
    ierr = DSOrthogonalize(d->eps->ds,DS_MAT_Q,size_X+size_plusk,&cMTX);CHKERRQ(ierr);
  }

  if (d->W && size_plusk > 0) {
    /* ps.Z <- orth([ps.Z(0:size_X-1) [oldV(0:size_plusk-1); 0] ]) */
    ierr = DSCopyMat(d->eps->ds,DS_MAT_Z,0,size_X,data->oldV,0,0,nV,size_plusk,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DSOrthogonalize(d->eps->ds,DS_MAT_Z,size_X+size_plusk,&cMTY);CHKERRQ(ierr);
    cMTX = PetscMin(cMTX, cMTY);
  }
  if (cMTX > size_X+size_plusk) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Invalid number of columns to restart");

  /* Notify the changes in V and update the other subspaces */
  d->V_tra_s = 0;                     d->V_tra_e = cMTX;
  d->V_new_s = d->V_tra_e;            d->V_new_e = d->V_new_s;

  /* Remove oldU */
  data->size_oldU = 0;

  /* Remove npreconv */
  d->npreconv = 0;
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_updateV_testConv(dvdDashboard *d,PetscInt s,PetscInt pre,PetscInt e,PetscInt *nConv)
{
  PetscInt        i,j,b;
  PetscReal       norm;
  PetscErrorCode  ierr;
  PetscBool       conv, c;
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;

  PetscFunctionBegin;
  if (nConv) *nConv = s;
  for (i=s,conv=PETSC_TRUE;(conv || data->allResiduals) && (i < e);i+=b) {
#if !defined(PETSC_USE_COMPLEX)
    b = d->eigi[i]!=0.0?2:1;
#else
    b = 1;
#endif
    if (i+b-1 >= pre) {
      ierr = d->calcpairs_residual(d,i,i+b);CHKERRQ(ierr);
    }
    /* Test the Schur vector */
    for (j=0,c=PETSC_TRUE;j<b && c;j++) {
      norm = d->nR[i+j]/d->nX[i+j];
      c = d->testConv(d,d->eigr[i+j],d->eigi[i+j],norm,&d->errest[i+j]);
    }
    if (conv && c) { if (nConv) *nConv = i+b; }
    else conv = PETSC_FALSE;
  }
  pre = PetscMax(pre,i);

#if !defined(PETSC_USE_COMPLEX)
  /* Enforce converged conjugate complex eigenpairs */
  if (nConv) {
    for (j=0;j<*nConv;j++) if (d->eigi[j] != 0.0) j++;
    if (j>*nConv) (*nConv)--;
  }
#endif
  for (i=pre;i<e;i++) d->errest[i] = d->nR[i] = 1.0;
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_updateV_update_gen(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        size_D,s,lV,kV,nV;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  /* Select the desired pairs */
  ierr = BVGetActiveColumns(d->eps->V,&lV,&kV);CHKERRQ(ierr);
  nV = kV - lV;
  size_D = PetscMin(PetscMin(PetscMin(d->bs,nV),d->eps->ncv-nV),nV);
  if (size_D == 0) PetscFunctionReturn(0);

  /* Fill V with D */
  ierr = d->improveX(d,d->npreconv,d->npreconv+size_D,&size_D);CHKERRQ(ierr);

  /* If D is empty, exit */
  d->size_D = size_D;
  if (size_D == 0) PetscFunctionReturn(0);

  /* Get the residual of all pairs */
#if !defined(PETSC_USE_COMPLEX)
  s = (d->eigi[0]!=0.0)? 2: 1;
#else
  s = 1;
#endif
  ierr = BVGetActiveColumns(d->eps->V,&lV,&kV);CHKERRQ(ierr);
  nV = kV - lV;
  ierr = dvd_updateV_testConv(d,s,s,data->allResiduals?nV:size_D,NULL);CHKERRQ(ierr);

  /* Notify the changes in V */
  d->V_tra_s = 0;                 d->V_tra_e = 0;
  d->V_new_s = nV;                d->V_new_e = nV+size_D;

  /* Save the projected eigenvectors */
  if (data->plusk > 0) {
    ierr = MatZeroEntries(data->oldU);CHKERRQ(ierr);
    data->size_oldU = nV;
    ierr = DSCopyMat(d->eps->ds,DS_MAT_Q,0,0,data->oldU,0,0,nV,nV,PETSC_TRUE);CHKERRQ(ierr);
    if (d->W) {
      ierr = MatZeroEntries(data->oldV);CHKERRQ(ierr);
      ierr = DSCopyMat(d->eps->ds,DS_MAT_Z,0,0,data->oldV,0,0,nV,nV,PETSC_TRUE);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode dvd_updateV_extrapol(dvdDashboard *d)
{
  dvdManagV_basic *data = (dvdManagV_basic*)d->updateV_data;
  PetscInt        i;
  PetscBool       restart,t;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  /* TODO: restrict select pairs to each case */
  ierr = d->calcpairs_selectPairs(d, data->min_size_V+d->npreconv);CHKERRQ(ierr);

  /* If the subspaces doesn't need restart, add new vector */
  ierr = d->isRestarting(d,&restart);CHKERRQ(ierr);
  if (!restart) {
    d->size_D = 0;
    ierr = dvd_updateV_update_gen(d);CHKERRQ(ierr);

    /* If no vector were converged, exit */
    /* For GHEP without B-ortho, converge all of the requested pairs at once */
    ierr = PetscObjectTypeCompareAny((PetscObject)d->eps->ds,&t,DSGHEP,"");CHKERRQ(ierr);
    if (d->nconv+d->npreconv < d->nev && (t || d->npreconv == 0)) PetscFunctionReturn(0);
  }

  /* If some eigenpairs were converged, lock them  */
  if (d->npreconv > 0) {
    i = d->npreconv;
    ierr = dvd_updateV_conv_gen(d);CHKERRQ(ierr);

    /* If some eigenpair was locked, exit */
    if (i > d->npreconv) PetscFunctionReturn(0);
  }

  /* Else, a restarting is performed */
  ierr = dvd_updateV_restart_gen(d);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode dvd_managementV_basic(dvdDashboard *d,dvdBlackboard *b,PetscInt bs,PetscInt mpd,PetscInt min_size_V,PetscInt plusk,PetscBool harm,PetscBool allResiduals)
{
  PetscErrorCode  ierr;
  dvdManagV_basic *data;
#if !defined(PETSC_USE_COMPLEX)
  PetscBool       her_probl,std_probl;
#endif

  PetscFunctionBegin;
  /* Setting configuration constrains */
#if !defined(PETSC_USE_COMPLEX)
  /* if the last converged eigenvalue is complex its conjugate pair is also
     converged */
  her_probl = DVD_IS(d->sEP,DVD_EP_HERMITIAN)? PETSC_TRUE: PETSC_FALSE;
  std_probl = DVD_IS(d->sEP,DVD_EP_STD)? PETSC_TRUE: PETSC_FALSE;
  b->max_size_X = PetscMax(b->max_size_X,bs+((her_probl && std_probl)?0:1));
#else
  b->max_size_X = PetscMax(b->max_size_X,bs);
#endif

  b->max_size_V = PetscMax(b->max_size_V,mpd);
  min_size_V = PetscMin(min_size_V,mpd-bs);
  b->size_V = PetscMax(b->size_V,b->max_size_V+b->max_size_P+b->max_nev);
  b->max_size_oldX = plusk;

  /* Setup the step */
  if (b->state >= DVD_STATE_CONF) {
    ierr = PetscNewLog(d->eps,&data);CHKERRQ(ierr);
    data->mpd = b->max_size_V;
    data->min_size_V = min_size_V;
    d->bs = bs;
    data->plusk = plusk;
    data->allResiduals = allResiduals;

    d->eigr = d->eps->eigr;
    d->eigi = d->eps->eigi;
    d->errest = d->eps->errest;
    ierr = PetscMalloc1(d->eps->ncv,&d->real_nR);CHKERRQ(ierr);
    ierr = PetscMalloc1(d->eps->ncv,&d->real_nX);CHKERRQ(ierr);
    if (plusk > 0) { ierr = MatCreateSeqDense(PETSC_COMM_SELF,d->eps->ncv,d->eps->ncv,NULL,&data->oldU);CHKERRQ(ierr); }
    else data->oldU = NULL;
    if (harm && plusk>0) { ierr = MatCreateSeqDense(PETSC_COMM_SELF,d->eps->ncv,d->eps->ncv,NULL,&data->oldV);CHKERRQ(ierr); }
    else data->oldV = NULL;

    data->old_updateV_data = d->updateV_data;
    d->updateV_data = data;
    data->old_isRestarting = d->isRestarting;
    d->isRestarting = dvd_isrestarting_fullV;
    d->updateV = dvd_updateV_extrapol;
    d->preTestConv = dvd_updateV_testConv;
    ierr = EPSDavidsonFLAdd(&d->startList,dvd_updateV_start);CHKERRQ(ierr);
    ierr = EPSDavidsonFLAdd(&d->destroyList,dvd_managementV_basic_d);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

