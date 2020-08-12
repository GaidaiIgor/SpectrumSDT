/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <slepc/private/vecimplslepc.h>     /*I "slepcvec.h" I*/

/* Private MPI datatypes and operators */
static MPI_Datatype MPIU_NORM2=0, MPIU_NORM1_AND_2=0;
static PetscBool VecCompInitialized = PETSC_FALSE;
MPI_Op MPIU_NORM2_SUM=0;

/* Private functions */
PETSC_STATIC_INLINE void SumNorm2(PetscReal*,PetscReal*,PetscReal*,PetscReal*);
PETSC_STATIC_INLINE PetscReal GetNorm2(PetscReal,PetscReal);
PETSC_STATIC_INLINE void AddNorm2(PetscReal*,PetscReal*,PetscReal);
static PetscErrorCode VecCompSetSubVecs_Comp(Vec,PetscInt,Vec*);
static PetscErrorCode VecCompGetSubVecs_Comp(Vec,PetscInt*,const Vec**);

#include "veccomp0.h"

#define __WITH_MPI__
#include "veccomp0.h"

PETSC_STATIC_INLINE void SumNorm2(PetscReal *ssq0,PetscReal *scale0,PetscReal *ssq1,PetscReal *scale1)
{
  PetscReal q;
  if (*scale0 > *scale1) {
    q = *scale1/(*scale0);
    *ssq1 = *ssq0 + q*q*(*ssq1);
    *scale1 = *scale0;
  } else {
    q = *scale0/(*scale1);
    *ssq1 += q*q*(*ssq0);
  }
}

PETSC_STATIC_INLINE PetscReal GetNorm2(PetscReal ssq,PetscReal scale)
{
  return scale*PetscSqrtReal(ssq);
}

PETSC_STATIC_INLINE void AddNorm2(PetscReal *ssq,PetscReal *scale,PetscReal x)
{
  PetscReal absx,q;
  if (x != 0.0) {
    absx = PetscAbs(x);
    if (*scale < absx) {
      q = *scale/absx;
      *ssq = 1.0 + *ssq*q*q;
      *scale = absx;
    } else {
      q = absx/(*scale);
      *ssq += q*q;
    }
  }
}

SLEPC_EXTERN void MPIAPI SlepcSumNorm2_Local(void *in,void *out,PetscMPIInt *cnt,MPI_Datatype *datatype)
{
  PetscInt  i,count = *cnt;
  PetscReal *xin = (PetscReal*)in,*xout = (PetscReal*)out;

  PetscFunctionBegin;
  if (*datatype == MPIU_NORM2) {
    for (i=0;i<count;i++) {
      SumNorm2(&xin[i*2],&xin[i*2+1],&xout[i*2],&xout[i*2+1]);
    }
  } else if (*datatype == MPIU_NORM1_AND_2) {
    for (i=0;i<count;i++) {
      xout[i*3] += xin[i*3];
      SumNorm2(&xin[i*3+1],&xin[i*3+2],&xout[i*3+1],&xout[i*3+2]);
    }
  } else {
    (*PetscErrorPrintf)("Can only handle MPIU_NORM* data types");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  PetscFunctionReturnVoid();
}

static PetscErrorCode VecCompNormEnd(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Type_free(&MPIU_NORM2);CHKERRQ(ierr);
  ierr = MPI_Type_free(&MPIU_NORM1_AND_2);CHKERRQ(ierr);
  ierr = MPI_Op_free(&MPIU_NORM2_SUM);CHKERRQ(ierr);
  VecCompInitialized = PETSC_FALSE;
  PetscFunctionReturn(0);
}

static PetscErrorCode VecCompNormInit(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Type_contiguous(2,MPIU_REAL,&MPIU_NORM2);CHKERRQ(ierr);
  ierr = MPI_Type_commit(&MPIU_NORM2);CHKERRQ(ierr);
  ierr = MPI_Type_contiguous(3,MPIU_REAL,&MPIU_NORM1_AND_2);CHKERRQ(ierr);
  ierr = MPI_Type_commit(&MPIU_NORM1_AND_2);CHKERRQ(ierr);
  ierr = MPI_Op_create(SlepcSumNorm2_Local,PETSC_TRUE,&MPIU_NORM2_SUM);CHKERRQ(ierr);
  ierr = PetscRegisterFinalize(VecCompNormEnd);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VecDestroy_Comp(Vec v)
{
  Vec_Comp       *vs = (Vec_Comp*)v->data;
  PetscInt       i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
#if defined(PETSC_USE_LOG)
  PetscLogObjectState((PetscObject)v,"Length=%D",v->map->n);
#endif
  for (i=0;i<vs->nx;i++) {
    ierr = VecDestroy(&vs->x[i]);CHKERRQ(ierr);
  }
  if (--vs->n->friends <= 0) {
    ierr = PetscFree(vs->n);CHKERRQ(ierr);
  }
  ierr = PetscFree(vs->x);CHKERRQ(ierr);
  ierr = PetscFree(vs);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)v,"VecCompSetSubVecs_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)v,"VecCompGetSubVecs_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static struct _VecOps DvOps = {VecDuplicate_Comp, /* 1 */
            VecDuplicateVecs_Comp,
            VecDestroyVecs_Comp,
            VecDot_Comp_MPI,
            VecMDot_Comp_MPI,
            VecNorm_Comp_MPI,
            VecTDot_Comp_MPI,
            VecMTDot_Comp_MPI,
            VecScale_Comp,
            VecCopy_Comp, /* 10 */
            VecSet_Comp,
            VecSwap_Comp,
            VecAXPY_Comp,
            VecAXPBY_Comp,
            VecMAXPY_Comp,
            VecAYPX_Comp,
            VecWAXPY_Comp,
            VecAXPBYPCZ_Comp,
            VecPointwiseMult_Comp,
            VecPointwiseDivide_Comp,
            0, /* 20 */
            0,0,
            0 /*VecGetArray_Seq*/,
            VecGetSize_Comp,
            VecGetLocalSize_Comp,
            0/*VecRestoreArray_Seq*/,
            VecMax_Comp,
            VecMin_Comp,
            VecSetRandom_Comp,
            0, /* 30 */
            0,
            VecDestroy_Comp,
            VecView_Comp,
            0/*VecPlaceArray_Seq*/,
            0/*VecReplaceArray_Seq*/,
            VecDot_Comp_Seq,
            VecTDot_Comp_Seq,
            VecNorm_Comp_Seq,
            VecMDot_Comp_Seq,
            VecMTDot_Comp_Seq, /* 40 */
            0,
            VecReciprocal_Comp,
            VecConjugate_Comp,
            0,0,
            0/*VecResetArray_Seq*/,
            0,
            VecMaxPointwiseDivide_Comp,
            VecPointwiseMax_Comp,
            VecPointwiseMaxAbs_Comp,
            VecPointwiseMin_Comp,
            0,
            VecSqrtAbs_Comp,
            VecAbs_Comp,
            VecExp_Comp,
            VecLog_Comp,
            VecShift_Comp,
            0,
            0,
            0,
            VecDotNorm2_Comp_MPI
          };

PetscErrorCode VecDuplicateVecs_Comp(Vec w,PetscInt m,Vec *V[])
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(w,VEC_CLASSID,1);
  PetscValidPointer(V,3);
  if (m<=0) SETERRQ1(PetscObjectComm((PetscObject)w),PETSC_ERR_ARG_OUTOFRANGE,"m must be > 0: m = %D",m);
  ierr = PetscMalloc1(m,V);CHKERRQ(ierr);
  for (i=0;i<m;i++) { ierr = VecDuplicate(w,*V+i);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

PetscErrorCode VecDestroyVecs_Comp(PetscInt m,Vec v[])
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidPointer(v,1);
  if (m<=0) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"m must be > 0: m = %D",m);
  for (i=0;i<m;i++) { ierr = VecDestroy(&v[i]);CHKERRQ(ierr); }
  ierr = PetscFree(v);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode VecCreate_Comp_Private(Vec v,Vec *x,PetscInt nx,PetscBool x_to_me,Vec_Comp_N *n)
{
  Vec_Comp       *s;
  PetscErrorCode ierr;
  PetscInt       N=0,lN=0,i,k;

  PetscFunctionBegin;
  if (!VecCompInitialized) {
    VecCompInitialized = PETSC_TRUE;
    ierr = VecRegister(VECCOMP,VecCreate_Comp);CHKERRQ(ierr);
    ierr = VecCompNormInit();CHKERRQ(ierr);
  }

  /* Allocate a new Vec_Comp */
  if (v->data) { ierr = PetscFree(v->data);CHKERRQ(ierr); }
  ierr = PetscNewLog(v,&s);CHKERRQ(ierr);
  ierr = PetscMemcpy(v->ops,&DvOps,sizeof(DvOps));CHKERRQ(ierr);
  v->data        = (void*)s;
  v->petscnative = PETSC_FALSE;

  /* Allocate the array of Vec, if it is needed to be done */
  if (!x_to_me) {
    if (nx) { ierr = PetscMalloc1(nx,&s->x);CHKERRQ(ierr); }
    if (x) { ierr = PetscArraycpy(s->x,x,nx);CHKERRQ(ierr); }
  } else s->x = x;

  s->nx = nx;

  if (nx && x) {
    /* Allocate the shared structure, if it is not given */
    if (!n) {
      for (i=0;i<nx;i++) {
        ierr = VecGetSize(x[i],&k);CHKERRQ(ierr);
        N+= k;
        ierr = VecGetLocalSize(x[i],&k);CHKERRQ(ierr);
        lN+= k;
      }
      ierr = PetscNewLog(v,&n);CHKERRQ(ierr);
      s->n = n;
      n->n = nx;
      n->N = N;
      n->lN = lN;
      n->friends = 1;
    } else { /* If not, check in the vector in the shared structure */
      s->n = n;
      s->n->friends++;
    }

    /* Set the virtual sizes as the real sizes of the vector */
    ierr = VecSetSizes(v,s->n->lN,s->n->N);CHKERRQ(ierr);
  }

  ierr = PetscObjectChangeTypeName((PetscObject)v,VECCOMP);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)v,"VecCompSetSubVecs_C",VecCompSetSubVecs_Comp);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)v,"VecCompGetSubVecs_C",VecCompGetSubVecs_Comp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

SLEPC_EXTERN PetscErrorCode VecCreate_Comp(Vec V)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecCreate_Comp_Private(V,NULL,0,PETSC_FALSE,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   VecCreateComp - Creates a new vector containing several subvectors,
   each stored separately.

   Collective

   Input Parameters:
+  comm - communicator for the new Vec
.  Nx   - array of (initial) global sizes of child vectors
.  n    - number of child vectors
.  t    - type of the child vectors
-  Vparent - (optional) template vector

   Output Parameter:
.  V - new vector

   Notes:
   This is similar to PETSc's VecNest but customized for SLEPc's needs. In particular,
   the number of child vectors can be modified dynamically, with VecCompSetSubVecs().

   Level: developer

.seealso: VecCreateCompWithVecs(), VecCompSetSubVecs()
@*/
PetscErrorCode VecCreateComp(MPI_Comm comm,PetscInt *Nx,PetscInt n,VecType t,Vec Vparent,Vec *V)
{
  PetscErrorCode ierr;
  Vec            *x;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = VecCreate(comm,V);CHKERRQ(ierr);
  ierr = PetscMalloc1(n,&x);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)*V,n*sizeof(Vec));CHKERRQ(ierr);
  for (i=0;i<n;i++) {
    ierr = VecCreate(comm,&x[i]);CHKERRQ(ierr);
    ierr = VecSetSizes(x[i],PETSC_DECIDE,Nx[i]);CHKERRQ(ierr);
    ierr = VecSetType(x[i],t);CHKERRQ(ierr);
  }
  ierr = VecCreate_Comp_Private(*V,x,n,PETSC_TRUE,Vparent?((Vec_Comp*)Vparent->data)->n:NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   VecCreateCompWithVecs - Creates a new vector containing several subvectors,
   each stored separately, from an array of Vecs.

   Collective on x

   Input Parameters:
+  x - array of Vecs
.  n - number of child vectors
-  Vparent - (optional) template vector

   Output Parameter:
.  V - new vector

   Level: developer

.seealso: VecCreateComp()
@*/
PetscErrorCode VecCreateCompWithVecs(Vec *x,PetscInt n,Vec Vparent,Vec *V)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  PetscValidPointer(x,1);
  PetscValidHeaderSpecific(*x,VEC_CLASSID,1);
  PetscValidLogicalCollectiveInt(*x,n,2);
  ierr = VecCreate(PetscObjectComm((PetscObject)x[0]),V);CHKERRQ(ierr);
  for (i=0;i<n;i++) {
    ierr = PetscObjectReference((PetscObject)x[i]);CHKERRQ(ierr);
  }
  ierr = VecCreate_Comp_Private(*V,x,n,PETSC_FALSE,Vparent?((Vec_Comp*)Vparent->data)->n:NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VecDuplicate_Comp(Vec win,Vec *V)
{
  PetscErrorCode ierr;
  Vec            *x;
  PetscInt       i;
  Vec_Comp       *s = (Vec_Comp*)win->data;

  PetscFunctionBegin;
  SlepcValidVecComp(win,1);
  ierr = VecCreate(PetscObjectComm((PetscObject)win),V);CHKERRQ(ierr);
  ierr = PetscMalloc1(s->nx,&x);CHKERRQ(ierr);
  ierr = PetscLogObjectMemory((PetscObject)*V,s->nx*sizeof(Vec));CHKERRQ(ierr);
  for (i=0;i<s->nx;i++) {
    if (s->x[i]) {
      ierr = VecDuplicate(s->x[i],&x[i]);CHKERRQ(ierr);
    } else x[i] = NULL;
  }
  ierr = VecCreate_Comp_Private(*V,x,s->nx,PETSC_TRUE,s->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode VecCompGetSubVecs_Comp(Vec win,PetscInt *n,const Vec **x)
{
  Vec_Comp *s = (Vec_Comp*)win->data;

  PetscFunctionBegin;
  if (x) *x = s->x;
  if (n) *n = s->n->n;
  PetscFunctionReturn(0);
}

/*@C
   VecCompGetSubVecs - Returns the entire array of vectors defining a
   compound vector.

   Collective on win

   Input Parameter:
.  win - compound vector

   Output Parameters:
+  n - number of child vectors
-  x - array of child vectors

   Level: developer

.seealso: VecCreateComp()
@*/
PetscErrorCode VecCompGetSubVecs(Vec win,PetscInt *n,const Vec **x)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(win,VEC_CLASSID,1);
  ierr = PetscUseMethod(win,"VecCompGetSubVecs_C",(Vec,PetscInt*,const Vec**),(win,n,x));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode VecCompSetSubVecs_Comp(Vec win,PetscInt n,Vec *x)
{
  Vec_Comp       *s = (Vec_Comp*)win->data;
  PetscInt       i,N,nlocal;
  Vec_Comp_N     *nn;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!s) SETERRQ(PetscObjectComm((PetscObject)win),PETSC_ERR_ORDER,"Must call VecSetSizes first");
  if (!s->nx) {
    /* vector has been created via VecCreate+VecSetType+VecSetSizes, so allocate data structures */
    ierr = PetscMalloc1(n,&s->x);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)win,n*sizeof(Vec));CHKERRQ(ierr);
    ierr = VecGetSize(win,&N);CHKERRQ(ierr);
    if (N%n) SETERRQ2(PetscObjectComm((PetscObject)win),1,"Global dimension %D is not divisible by %D",N,n);
    ierr = VecGetLocalSize(win,&nlocal);CHKERRQ(ierr);
    if (nlocal%n) SETERRQ2(PetscObjectComm((PetscObject)win),1,"Local dimension %D is not divisible by %D",nlocal,n);
    s->nx = n;
    for (i=0;i<n;i++) {
      ierr = VecCreate(PetscObjectComm((PetscObject)win),&s->x[i]);CHKERRQ(ierr);
      ierr = VecSetSizes(s->x[i],nlocal/n,N/n);CHKERRQ(ierr);
      ierr = VecSetFromOptions(s->x[i]);CHKERRQ(ierr);
    }
    if (!s->n) {
      ierr = PetscNewLog(win,&nn);CHKERRQ(ierr);
      s->n = nn;
      nn->N = N;
      nn->lN = nlocal;
      nn->friends = 1;
    }
  } else if (n > s->nx) SETERRQ1(PetscObjectComm((PetscObject)win),PETSC_ERR_SUP,"Number of child vectors cannot be larger than %D",s->nx);
  if (x) {
    ierr = PetscArraycpy(s->x,x,n);CHKERRQ(ierr);
  }
  s->n->n = n;
  PetscFunctionReturn(0);
}

/*@C
   VecCompSetSubVecs - Resets the number of subvectors defining a compound vector,
   or replaces the subvectors.

   Collective on win

   Input Parameters:
+  win - compound vector
.  n - number of child vectors
-  x - array of child vectors

   Note:
   It is not possible to increase the number of subvectors with respect to the
   number set at its creation.

   Level: developer

.seealso: VecCreateComp(), VecCompGetSubVecs()
@*/
PetscErrorCode VecCompSetSubVecs(Vec win,PetscInt n,Vec *x)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(win,VEC_CLASSID,1);
  PetscValidLogicalCollectiveInt(win,n,2);
  ierr = PetscTryMethod(win,"VecCompSetSubVecs_C",(Vec,PetscInt,Vec*),(win,n,x));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VecAXPY_Comp(Vec v,PetscScalar alpha,Vec w)
{
  PetscErrorCode ierr;
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,3);
  for (i=0;i<vs->n->n;i++) {
    ierr = VecAXPY(vs->x[i],alpha,ws->x[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VecAYPX_Comp(Vec v,PetscScalar alpha,Vec w)
{
  PetscErrorCode ierr;
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,3);
  for (i=0;i<vs->n->n;i++) {
    ierr = VecAYPX(vs->x[i],alpha,ws->x[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VecAXPBY_Comp(Vec v,PetscScalar alpha,PetscScalar beta,Vec w)
{
  PetscErrorCode ierr;
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,4);
  for (i=0;i<vs->n->n;i++) {
    ierr = VecAXPBY(vs->x[i],alpha,beta,ws->x[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VecMAXPY_Comp(Vec v,PetscInt n,const PetscScalar *alpha,Vec *w)
{
  PetscErrorCode ierr;
  Vec_Comp       *vs = (Vec_Comp*)v->data;
  Vec            *wx;
  PetscInt       i,j;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  for (i=0;i<n;i++) SlepcValidVecComp(w[i],4);

  ierr = PetscMalloc1(n,&wx);CHKERRQ(ierr);

  for (j=0;j<vs->n->n;j++) {
    for (i=0;i<n;i++) wx[i] = ((Vec_Comp*)w[i]->data)->x[j];
    ierr = VecMAXPY(vs->x[j],n,alpha,wx);CHKERRQ(ierr);
  }

  ierr = PetscFree(wx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode VecWAXPY_Comp(Vec v,PetscScalar alpha,Vec w,Vec z)
{
  PetscErrorCode ierr;
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data,*zs = (Vec_Comp*)z->data;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,3);
  SlepcValidVecComp(z,4);
  for (i=0;i<vs->n->n;i++) {
    ierr = VecWAXPY(vs->x[i],alpha,ws->x[i],zs->x[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VecAXPBYPCZ_Comp(Vec v,PetscScalar alpha,PetscScalar beta,PetscScalar gamma,Vec w,Vec z)
{
  PetscErrorCode  ierr;
  Vec_Comp        *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data,*zs = (Vec_Comp*)z->data;
  PetscInt        i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,5);
  SlepcValidVecComp(z,6);
  for (i=0;i<vs->n->n;i++) {
    ierr = VecAXPBYPCZ(vs->x[i],alpha,beta,gamma,ws->x[i],zs->x[i]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VecGetSize_Comp(Vec v,PetscInt *size)
{
  Vec_Comp *vs = (Vec_Comp*)v->data;

  PetscFunctionBegin;
  PetscValidIntPointer(size,2);
  if (vs->n) {
    SlepcValidVecComp(v,1);
    *size = vs->n->N;
  } else *size = v->map->N;
  PetscFunctionReturn(0);
}

PetscErrorCode VecGetLocalSize_Comp(Vec v,PetscInt *size)
{
  Vec_Comp *vs = (Vec_Comp*)v->data;

  PetscFunctionBegin;
  PetscValidIntPointer(size,2);
  if (vs->n) {
    SlepcValidVecComp(v,1);
    *size = vs->n->lN;
  } else *size = v->map->n;
  PetscFunctionReturn(0);
}

PetscErrorCode VecMax_Comp(Vec v,PetscInt *idx,PetscReal *z)
{
  PetscErrorCode ierr;
  Vec_Comp       *vs = (Vec_Comp*)v->data;
  PetscInt       idxp,s=0,s0;
  PetscReal      zp,z0;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  if (!idx && !z) PetscFunctionReturn(0);

  if (vs->n->n > 0) {
    ierr = VecMax(vs->x[0],idx?&idxp:NULL,&zp);CHKERRQ(ierr);
  } else {
    zp = PETSC_MIN_REAL;
    if (idx) idxp = -1;
  }
  for (i=1;i<vs->n->n;i++) {
    ierr = VecGetSize(vs->x[i-1],&s0);CHKERRQ(ierr);
    s += s0;
    ierr = VecMax(vs->x[i],idx?&idxp:NULL,&z0);CHKERRQ(ierr);
    if (zp < z0) {
      if (idx) *idx = s+idxp;
      zp = z0;
    }
  }
  if (z) *z = zp;
  PetscFunctionReturn(0);
}

PetscErrorCode VecMin_Comp(Vec v,PetscInt *idx,PetscReal *z)
{
  PetscErrorCode ierr;
  Vec_Comp       *vs = (Vec_Comp*)v->data;
  PetscInt       idxp,s=0,s0;
  PetscReal      zp,z0;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  if (!idx && !z) PetscFunctionReturn(0);

  if (vs->n->n > 0) {
    ierr = VecMin(vs->x[0],idx?&idxp:NULL,&zp);CHKERRQ(ierr);
  } else {
    zp = PETSC_MAX_REAL;
    if (idx) idxp = -1;
  }
  for (i=1;i<vs->n->n;i++) {
    ierr = VecGetSize(vs->x[i-1],&s0);CHKERRQ(ierr);
    s += s0;
    ierr = VecMin(vs->x[i],idx?&idxp:NULL,&z0);CHKERRQ(ierr);
    if (zp > z0) {
      if (idx) *idx = s+idxp;
      zp = z0;
    }
  }
  if (z) *z = zp;
  PetscFunctionReturn(0);
}

PetscErrorCode VecMaxPointwiseDivide_Comp(Vec v,Vec w,PetscReal *m)
{
  PetscErrorCode ierr;
  Vec_Comp       *vs = (Vec_Comp*)v->data,*ws = (Vec_Comp*)w->data;
  PetscReal      work;
  PetscInt       i;

  PetscFunctionBegin;
  SlepcValidVecComp(v,1);
  SlepcValidVecComp(w,2);
  if (!m || vs->n->n == 0) PetscFunctionReturn(0);
  ierr = VecMaxPointwiseDivide(vs->x[0],ws->x[0],m);CHKERRQ(ierr);
  for (i=1;i<vs->n->n;i++) {
    ierr = VecMaxPointwiseDivide(vs->x[i],ws->x[i],&work);CHKERRQ(ierr);
    *m = PetscMax(*m,work);
  }
  PetscFunctionReturn(0);
}


#define __QUOTEME__(x) #x
#define __COMPOSE2__(A,B) A##B
#define __COMPOSE3__(A,B,C) A##B##C

#define __FUNC_TEMPLATE1__(NAME) \
PetscErrorCode __COMPOSE3__(Vec,NAME,_Comp)(Vec v) \
{ \
  PetscErrorCode  ierr; \
  Vec_Comp        *vs = (Vec_Comp*)v->data; \
  PetscInt        i; \
\
  PetscFunctionBegin; \
  SlepcValidVecComp(v,1); \
  for (i=0;i<vs->n->n;i++) { \
    ierr = __COMPOSE2__(Vec,NAME)(vs->x[i]);CHKERRQ(ierr); \
  } \
  PetscFunctionReturn(0);\
}

__FUNC_TEMPLATE1__(Conjugate)
__FUNC_TEMPLATE1__(Reciprocal)
__FUNC_TEMPLATE1__(SqrtAbs)
__FUNC_TEMPLATE1__(Abs)
__FUNC_TEMPLATE1__(Exp)
__FUNC_TEMPLATE1__(Log)

#define __FUNC_TEMPLATE2__(NAME,T0) \
PetscErrorCode __COMPOSE3__(Vec,NAME,_Comp)(Vec v,T0 __a) \
{ \
  PetscErrorCode  ierr; \
  Vec_Comp        *vs = (Vec_Comp*)v->data; \
  PetscInt        i; \
\
  PetscFunctionBegin; \
  SlepcValidVecComp(v,1); \
  for (i=0;i<vs->n->n;i++) { \
    ierr = __COMPOSE2__(Vec,NAME)(vs->x[i],__a);CHKERRQ(ierr); \
  } \
  PetscFunctionReturn(0);\
}

__FUNC_TEMPLATE2__(Set,PetscScalar)
__FUNC_TEMPLATE2__(View,PetscViewer)
__FUNC_TEMPLATE2__(Scale,PetscScalar)
__FUNC_TEMPLATE2__(SetRandom,PetscRandom)
__FUNC_TEMPLATE2__(Shift,PetscScalar)

#define __FUNC_TEMPLATE3__(NAME) \
PetscErrorCode __COMPOSE3__(Vec,NAME,_Comp)(Vec v,Vec w) \
{ \
  PetscErrorCode  ierr; \
  Vec_Comp        *vs = (Vec_Comp*)v->data,\
                  *ws = (Vec_Comp*)w->data; \
  PetscInt        i; \
\
  PetscFunctionBegin; \
  SlepcValidVecComp(v,1); \
  SlepcValidVecComp(w,2); \
  for (i=0;i<vs->n->n;i++) { \
    ierr = __COMPOSE2__(Vec,NAME)(vs->x[i],ws->x[i]);CHKERRQ(ierr); \
  } \
  PetscFunctionReturn(0);\
}

__FUNC_TEMPLATE3__(Copy)
__FUNC_TEMPLATE3__(Swap)

#define __FUNC_TEMPLATE4__(NAME) \
PetscErrorCode __COMPOSE3__(Vec,NAME,_Comp)(Vec v,Vec w,Vec z) \
{ \
  PetscErrorCode  ierr; \
  Vec_Comp        *vs = (Vec_Comp*)v->data, \
                  *ws = (Vec_Comp*)w->data, \
                  *zs = (Vec_Comp*)z->data; \
  PetscInt        i; \
\
  PetscFunctionBegin; \
  SlepcValidVecComp(v,1); \
  SlepcValidVecComp(w,2); \
  SlepcValidVecComp(z,3); \
  for (i=0;i<vs->n->n;i++) { \
    ierr = __COMPOSE2__(Vec,NAME)(vs->x[i],ws->x[i],zs->x[i]);CHKERRQ(ierr); \
  } \
  PetscFunctionReturn(0);\
}

__FUNC_TEMPLATE4__(PointwiseMax)
__FUNC_TEMPLATE4__(PointwiseMaxAbs)
__FUNC_TEMPLATE4__(PointwiseMin)
__FUNC_TEMPLATE4__(PointwiseMult)
__FUNC_TEMPLATE4__(PointwiseDivide)
