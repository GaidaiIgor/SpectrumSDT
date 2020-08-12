/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#if !defined(SLEPCBVIMPL_H)
#define SLEPCBVIMPL_H

#include <slepcbv.h>
#include <slepc/private/slepcimpl.h>

SLEPC_EXTERN PetscBool BVRegisterAllCalled;
SLEPC_EXTERN PetscErrorCode BVRegisterAll(void);

SLEPC_EXTERN PetscLogEvent BV_Create,BV_Copy,BV_Mult,BV_MultVec,BV_MultInPlace,BV_Dot,BV_DotVec,BV_Orthogonalize,BV_OrthogonalizeVec,BV_Scale,BV_Norm,BV_NormVec,BV_SetRandom,BV_MatMult,BV_MatMultVec,BV_MatProject;

typedef struct _BVOps *BVOps;

struct _BVOps {
  PetscErrorCode (*mult)(BV,PetscScalar,PetscScalar,BV,Mat);
  PetscErrorCode (*multvec)(BV,PetscScalar,PetscScalar,Vec,PetscScalar*);
  PetscErrorCode (*multinplace)(BV,Mat,PetscInt,PetscInt);
  PetscErrorCode (*multinplacetrans)(BV,Mat,PetscInt,PetscInt);
  PetscErrorCode (*dot)(BV,BV,Mat);
  PetscErrorCode (*dotvec)(BV,Vec,PetscScalar*);
  PetscErrorCode (*dotvec_local)(BV,Vec,PetscScalar*);
  PetscErrorCode (*dotvec_begin)(BV,Vec,PetscScalar*);
  PetscErrorCode (*dotvec_end)(BV,Vec,PetscScalar*);
  PetscErrorCode (*scale)(BV,PetscInt,PetscScalar);
  PetscErrorCode (*norm)(BV,PetscInt,NormType,PetscReal*);
  PetscErrorCode (*norm_local)(BV,PetscInt,NormType,PetscReal*);
  PetscErrorCode (*norm_begin)(BV,PetscInt,NormType,PetscReal*);
  PetscErrorCode (*norm_end)(BV,PetscInt,NormType,PetscReal*);
  PetscErrorCode (*matmult)(BV,Mat,BV);
  PetscErrorCode (*copy)(BV,BV);
  PetscErrorCode (*copycolumn)(BV,PetscInt,PetscInt);
  PetscErrorCode (*resize)(BV,PetscInt,PetscBool);
  PetscErrorCode (*getcolumn)(BV,PetscInt,Vec*);
  PetscErrorCode (*restorecolumn)(BV,PetscInt,Vec*);
  PetscErrorCode (*getarray)(BV,PetscScalar**);
  PetscErrorCode (*restorearray)(BV,PetscScalar**);
  PetscErrorCode (*getarrayread)(BV,const PetscScalar**);
  PetscErrorCode (*restorearrayread)(BV,const PetscScalar**);
  PetscErrorCode (*restoresplit)(BV,BV*,BV*);
  PetscErrorCode (*gramschmidt)(BV,PetscInt,Vec,PetscBool*,PetscScalar*,PetscScalar*,PetscReal*,PetscReal*);
  PetscErrorCode (*duplicate)(BV,BV);
  PetscErrorCode (*create)(BV);
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,BV);
  PetscErrorCode (*view)(BV,PetscViewer);
  PetscErrorCode (*destroy)(BV);
};

struct _p_BV {
  PETSCHEADER(struct _BVOps);
  /*------------------------- User parameters --------------------------*/
  Vec                t;            /* template vector */
  PetscInt           n,N;          /* dimensions of vectors (local, global) */
  PetscInt           m;            /* number of vectors */
  PetscInt           l;            /* number of leading columns */
  PetscInt           k;            /* number of active columns */
  PetscInt           nc;           /* number of constraints */
  BVOrthogType       orthog_type;  /* the method of vector orthogonalization */
  BVOrthogRefineType orthog_ref;   /* refinement method */
  PetscReal          orthog_eta;   /* refinement threshold */
  BVOrthogBlockType  orthog_block; /* the method of block orthogonalization */
  Mat                matrix;       /* inner product matrix */
  PetscBool          indef;        /* matrix is indefinite */
  BVMatMultType      vmm;          /* version of matmult operation */
  PetscBool          rrandom;      /* reproducible random vectors */

  /*---------------------- Cached data and workspace -------------------*/
  Vec                buffer;       /* buffer vector used in orthogonalization */
  Mat                Abuffer;      /* auxiliary seqdense matrix that wraps the buffer */
  Vec                Bx;           /* result of matrix times a vector x */
  PetscObjectId      xid;          /* object id of vector x */
  PetscObjectState   xstate;       /* state of vector x */
  Vec                cv[2];        /* column vectors obtained with BVGetColumn() */
  PetscInt           ci[2];        /* column indices of obtained vectors */
  PetscObjectState   st[2];        /* state of obtained vectors */
  PetscObjectId      id[2];        /* object id of obtained vectors */
  PetscScalar        *h,*c;        /* orthogonalization coefficients */
  Vec                omega;        /* signature matrix values for indefinite case */
  PetscBool          defersfo;     /* deferred call to setfromoptions */
  BV                 cached;       /* cached BV to store result of matrix times BV */
  PetscObjectState   bvstate;      /* state of BV when BVApplyMatrixBV() was called */
  BV                 L,R;          /* BV objects obtained with BVGetSplit() */
  PetscObjectState   lstate,rstate;/* state of L and R when BVGetSplit() was called */
  PetscInt           lsplit;       /* the value of l when BVGetSplit() was called */
  PetscInt           issplit;      /* >0 if this BV has been created by splitting (1=left, 2=right) */
  BV                 splitparent;  /* my parent if I am a split BV */
  PetscRandom        rand;         /* random number generator */
  Mat                Acreate;      /* matrix given at BVCreateFromMat() */
  Mat                Aget;         /* matrix returned for BVGetMat() */
  PetscBool          cuda;         /* true if GPU must be used in SVEC */
  PetscScalar        *work;
  PetscInt           lwork;
  void               *data;
};

/*
  BV_SafeSqrt - Computes the square root of a scalar value alpha, which is
  assumed to be z'*B*z. The result is
    if definite inner product:     res = sqrt(alpha)
    if indefinite inner product:   res = sgn(alpha)*sqrt(abs(alpha))
*/
PETSC_STATIC_INLINE PetscErrorCode BV_SafeSqrt(BV bv,PetscScalar alpha,PetscReal *res)
{
  PetscErrorCode ierr;
  PetscReal      absal,realp;

  PetscFunctionBegin;
  absal = PetscAbsScalar(alpha);
  realp = PetscRealPart(alpha);
  if (absal<PETSC_MACHINE_EPSILON) {
    ierr = PetscInfo(bv,"Zero norm, either the vector is zero or a semi-inner product is being used\n");CHKERRQ(ierr);
  }
#if defined(PETSC_USE_COMPLEX)
  if (PetscAbsReal(PetscImaginaryPart(alpha))>10*PETSC_MACHINE_EPSILON && PetscAbsReal(PetscImaginaryPart(alpha))/absal>100*PETSC_MACHINE_EPSILON) SETERRQ1(PetscObjectComm((PetscObject)bv),1,"The inner product is not well defined: nonzero imaginary part %g",PetscImaginaryPart(alpha));
#endif
  if (bv->indef) {
    *res = (realp<0.0)? -PetscSqrtReal(-realp): PetscSqrtReal(realp);
  } else {
    if (realp<-10*PETSC_MACHINE_EPSILON) SETERRQ(PetscObjectComm((PetscObject)bv),1,"The inner product is not well defined: indefinite matrix");
    *res = (realp<0.0)? 0.0: PetscSqrtReal(realp);
  }
  PetscFunctionReturn(0);
}

/*
  BV_IPMatMult - Multiply a vector x by the inner-product matrix, cache the
  result in Bx.
*/
PETSC_STATIC_INLINE PetscErrorCode BV_IPMatMult(BV bv,Vec x)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (((PetscObject)x)->id != bv->xid || ((PetscObject)x)->state != bv->xstate) {
    if (!bv->Bx) {
      ierr = MatCreateVecs(bv->matrix,&bv->Bx,NULL);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)bv->Bx);CHKERRQ(ierr);
    }
    ierr = MatMult(bv->matrix,x,bv->Bx);CHKERRQ(ierr);
    ierr = PetscObjectGetId((PetscObject)x,&bv->xid);CHKERRQ(ierr);
    ierr = PetscObjectStateGet((PetscObject)x,&bv->xstate);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  BV_IPMatMultBV - Multiply BV by the inner-product matrix, cache the
  result internally in bv->cached.
*/
PETSC_STATIC_INLINE PetscErrorCode BV_IPMatMultBV(BV bv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = BVGetCachedBV(bv,&bv->cached);CHKERRQ(ierr);
  if (((PetscObject)bv)->state != bv->bvstate || bv->l != bv->cached->l || bv->k != bv->cached->k) {
    ierr = BVSetActiveColumns(bv->cached,bv->l,bv->k);CHKERRQ(ierr);
    if (bv->matrix) {
      ierr = BVMatMult(bv,bv->matrix,bv->cached);CHKERRQ(ierr);
    } else {
      ierr = BVCopy(bv,bv->cached);CHKERRQ(ierr);
    }
    bv->bvstate = ((PetscObject)bv)->state;
  }
  PetscFunctionReturn(0);
}

/*
  BV_AllocateCoeffs - Allocate orthogonalization coefficients if not done already.
*/
PETSC_STATIC_INLINE PetscErrorCode BV_AllocateCoeffs(BV bv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!bv->h) {
    ierr = PetscMalloc2(bv->nc+bv->m,&bv->h,bv->nc+bv->m,&bv->c);CHKERRQ(ierr);
    ierr = PetscLogObjectMemory((PetscObject)bv,2*bv->m*sizeof(PetscScalar));CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  BV_AllocateSignature - Allocate signature coefficients if not done already.
*/
PETSC_STATIC_INLINE PetscErrorCode BV_AllocateSignature(BV bv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (bv->indef && !bv->omega) {
    if (bv->cuda) {
#if defined(PETSC_HAVE_CUDA)
      ierr = VecCreateSeqCUDA(PETSC_COMM_SELF,bv->nc+bv->m,&bv->omega);CHKERRQ(ierr);
#else
      SETERRQ(PetscObjectComm((PetscObject)bv),1,"Something wrong happened");
#endif
    } else {
      ierr = VecCreateSeq(PETSC_COMM_SELF,bv->nc+bv->m,&bv->omega);CHKERRQ(ierr);
    }
    ierr = PetscLogObjectParent((PetscObject)bv,(PetscObject)bv->omega);CHKERRQ(ierr);
    ierr = VecSet(bv->omega,1.0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  BVAvailableVec: First (0) or second (1) vector available for
  getcolumn operation (or -1 if both vectors already fetched).
*/
#define BVAvailableVec (((bv->ci[0]==-bv->nc-1)? 0: (bv->ci[1]==-bv->nc-1)? 1: -1))

/*
    Macros to test valid BV arguments
*/
#if !defined(PETSC_USE_DEBUG)

#define BVCheckSizes(h,arg) do {} while (0)
#define BVCheckOp(h,arg,op) do {} while (0)

#else

#define BVCheckSizes(h,arg) \
  do { \
    if (!(h)->m) SETERRQ1(PetscObjectComm((PetscObject)(h)),PETSC_ERR_ARG_WRONGSTATE,"BV sizes have not been defined: Parameter #%d",arg); \
  } while (0)

#define BVCheckOp(h,arg,op) \
  do { \
    if (!(h)->ops->op) SETERRQ1(PetscObjectComm((PetscObject)(h)),PETSC_ERR_SUP,"Operation not implemented in this BV type: Parameter #%d",arg); \
  } while (0)

#endif

SLEPC_INTERN PetscErrorCode BVView_Vecs(BV,PetscViewer);

SLEPC_INTERN PetscErrorCode BVAllocateWork_Private(BV,PetscInt);

SLEPC_INTERN PetscErrorCode BVMult_BLAS_Private(BV,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar,const PetscScalar*,const PetscScalar*,PetscScalar,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVMultVec_BLAS_Private(BV,PetscInt,PetscInt,PetscScalar,const PetscScalar*,const PetscScalar*,PetscScalar,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVMultInPlace_BLAS_Private(BV,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar*,const PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BVMultInPlace_Vecs_Private(BV,PetscInt,PetscInt,PetscInt,Vec*,const PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BVAXPY_BLAS_Private(BV,PetscInt,PetscInt,PetscScalar,const PetscScalar*,PetscScalar,PetscScalar*);
SLEPC_INTERN PetscErrorCode BVDot_BLAS_Private(BV,PetscInt,PetscInt,PetscInt,PetscInt,const PetscScalar*,const PetscScalar*,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BVDotVec_BLAS_Private(BV,PetscInt,PetscInt,const PetscScalar*,const PetscScalar*,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BVScale_BLAS_Private(BV,PetscInt,PetscScalar*,PetscScalar);
SLEPC_INTERN PetscErrorCode BVNorm_LAPACK_Private(BV,PetscInt,PetscInt,const PetscScalar*,NormType,PetscReal*,PetscBool);
SLEPC_INTERN PetscErrorCode BVMatCholInv_LAPACK_Private(BV,Mat,Mat);
SLEPC_INTERN PetscErrorCode BVMatTriInv_LAPACK_Private(BV,Mat,Mat);
SLEPC_INTERN PetscErrorCode BVMatSVQB_LAPACK_Private(BV,Mat,Mat);
SLEPC_INTERN PetscErrorCode BVOrthogonalize_LAPACK_TSQR(BV,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscInt);
SLEPC_INTERN PetscErrorCode BVOrthogonalize_LAPACK_TSQR_OnlyR(BV,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscInt);

/* reduction operation used in BVOrthogonalize */
SLEPC_EXTERN MPI_Op MPIU_TSQR;
SLEPC_EXTERN void MPIAPI SlepcGivensPacked(void*,void*,PetscMPIInt*,MPI_Datatype*);

#if defined(PETSC_HAVE_CUDA)

/* complex single */
#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#define cublasXgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n) cublasCgemm((a),(b),(c),(d),(e),(f),(cuComplex*)(g),(cuComplex*)(h),(i),(cuComplex*)(j),(k),(cuComplex*)(l),(cuComplex*)(m),(n))
#define cublasXgemv(a,b,c,d,e,f,g,h,i,j,k,l) cublasCgemv((a),(b),(c),(d),(cuComplex*)(e),(cuComplex*)(f),(g),(cuComplex*)(h),(i),(cuComplex*)(j),(cuComplex*)(k),(l))
#define cublasXscal(a,b,c,d,e) cublasCscal(a,b,(const cuComplex*)(c),(cuComplex*)(d),e)
#define cublasXnrm2(a,b,c,d,e) cublasScnrm2(a,b,(const cuComplex*)(c),d,e)
#define cublasXaxpy(a,b,c,d,e,f,g) cublasCaxpy((a),(b),(cuComplex*)(c),(cuComplex*)(d),(e),(cuComplex*)(f),(g))
#define cublasXdotc(a,b,c,d,e,f,g) cublasCdotc((a),(b),(const cuComplex *)(c),(d),(const cuComplex *)(e),(f),(cuComplex *)(g))
#else /* complex double */
#define cublasXgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n) cublasZgemm((a),(b),(c),(d),(e),(f),(cuDoubleComplex*)(g),(cuDoubleComplex*)(h),(i),(cuDoubleComplex*)(j),(k),(cuDoubleComplex*)(l),(cuDoubleComplex*)(m),(n))
#define cublasXgemv(a,b,c,d,e,f,g,h,i,j,k,l) cublasZgemv((a),(b),(c),(d),(cuDoubleComplex*)(e),(cuDoubleComplex*)(f),(g),(cuDoubleComplex*)(h),(i),(cuDoubleComplex*)(j),(cuDoubleComplex*)(k),(l))
#define cublasXscal(a,b,c,d,e) cublasZscal(a,b,(const cuDoubleComplex*)(c),(cuDoubleComplex*)(d),e)
#define cublasXnrm2(a,b,c,d,e) cublasDznrm2(a,b,(const cuDoubleComplex*)(c),d,e)
#define cublasXaxpy(a,b,c,d,e,f,g) cublasZaxpy((a),(b),(cuDoubleComplex*)(c),(cuDoubleComplex*)(d),(e),(cuDoubleComplex*)(f),(g))
#define cublasXdotc(a,b,c,d,e,f,g) cublasZdotc((a),(b),(const cuDoubleComplex *)(c),(d),(const cuDoubleComplex *)(e),(f),(cuDoubleComplex *)(g))
#endif
#else /* real single */
#if defined(PETSC_USE_REAL_SINGLE)
#define cublasXgemm cublasSgemm
#define cublasXgemv cublasSgemv
#define cublasXscal cublasSscal
#define cublasXnrm2 cublasSnrm2
#define cublasXaxpy cublasSaxpy
#define cublasXdotc cublasSdot
#else /* real double */
#define cublasXgemm cublasDgemm
#define cublasXgemv cublasDgemv
#define cublasXscal cublasDscal
#define cublasXnrm2 cublasDnrm2
#define cublasXaxpy cublasDaxpy
#define cublasXdotc cublasDdot
#endif
#endif

SLEPC_INTERN PetscErrorCode BV_CleanCoefficients_CUDA(BV,PetscInt,PetscScalar*);
SLEPC_INTERN PetscErrorCode BV_AddCoefficients_CUDA(BV,PetscInt,PetscScalar*,PetscScalar*);
SLEPC_INTERN PetscErrorCode BV_SetValue_CUDA(BV,PetscInt,PetscInt,PetscScalar*,PetscScalar);
SLEPC_INTERN PetscErrorCode BV_SquareSum_CUDA(BV,PetscInt,PetscScalar*,PetscReal*);
SLEPC_INTERN PetscErrorCode BV_ApplySignature_CUDA(BV,PetscInt,PetscScalar*,PetscBool);
SLEPC_INTERN PetscErrorCode BV_SquareRoot_CUDA(BV,PetscInt,PetscScalar*,PetscReal*);
SLEPC_INTERN PetscErrorCode BV_StoreCoefficients_CUDA(BV,PetscInt,PetscScalar*,PetscScalar*);

#endif /* PETSC_HAVE_CUDA */

/*
   BV_CleanCoefficients_Default - Sets to zero all entries of column j of the bv buffer
*/
PETSC_STATIC_INLINE PetscErrorCode BV_CleanCoefficients_Default(BV bv,PetscInt j,PetscScalar *h)
{
  PetscErrorCode ierr;
  PetscScalar    *hh=h,*a;
  PetscInt       i;

  PetscFunctionBegin;
  if (!h) {
    ierr = VecGetArray(bv->buffer,&a);CHKERRQ(ierr);
    hh = a + j*(bv->nc+bv->m);
  }
  for (i=0;i<bv->nc+j;i++) hh[i] = 0.0;
  if (!h) { ierr = VecRestoreArray(bv->buffer,&a);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   BV_AddCoefficients_Default - Add the contents of the scratch (0-th column) of the bv buffer
   into column j of the bv buffer
*/
PETSC_STATIC_INLINE PetscErrorCode BV_AddCoefficients_Default(BV bv,PetscInt j,PetscScalar *h,PetscScalar *c)
{
  PetscErrorCode ierr;
  PetscScalar    *hh=h,*cc=c;
  PetscInt       i;

  PetscFunctionBegin;
  if (!h) {
    ierr = VecGetArray(bv->buffer,&cc);CHKERRQ(ierr);
    hh = cc + j*(bv->nc+bv->m);
  }
  for (i=0;i<bv->nc+j;i++) hh[i] += cc[i];
  if (!h) { ierr = VecRestoreArray(bv->buffer,&cc);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   BV_SetValue_Default - Sets value in row j (counted after the constraints) of column k
   of the coefficients array
*/
PETSC_STATIC_INLINE PetscErrorCode BV_SetValue_Default(BV bv,PetscInt j,PetscInt k,PetscScalar *h,PetscScalar value)
{
  PetscErrorCode ierr;
  PetscScalar    *hh=h,*a;

  PetscFunctionBegin;
  if (!h) {
    ierr = VecGetArray(bv->buffer,&a);CHKERRQ(ierr);
    hh = a + k*(bv->nc+bv->m);
  }
  hh[bv->nc+j] = value;
  if (!h) { ierr = VecRestoreArray(bv->buffer,&a);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   BV_SquareSum_Default - Returns the value h'*h, where h represents the contents of the
   coefficients array (up to position j)
*/
PETSC_STATIC_INLINE PetscErrorCode BV_SquareSum_Default(BV bv,PetscInt j,PetscScalar *h,PetscReal *sum)
{
  PetscErrorCode ierr;
  PetscScalar    *hh=h;
  PetscInt       i;

  PetscFunctionBegin;
  *sum = 0.0;
  if (!h) { ierr = VecGetArray(bv->buffer,&hh);CHKERRQ(ierr); }
  for (i=0;i<bv->nc+j;i++) *sum += PetscRealPart(hh[i]*PetscConj(hh[i]));
  if (!h) { ierr = VecRestoreArray(bv->buffer,&hh);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   BV_ApplySignature_Default - Computes the pointwise product h*omega, where h represents
   the contents of the coefficients array (up to position j) and omega is the signature;
   if inverse=TRUE then the operation is h/omega
*/
PETSC_STATIC_INLINE PetscErrorCode BV_ApplySignature_Default(BV bv,PetscInt j,PetscScalar *h,PetscBool inverse)
{
  PetscErrorCode    ierr;
  PetscScalar       *hh=h;
  PetscInt          i;
  const PetscScalar *omega;

  PetscFunctionBegin;
  if (!(bv->nc+j)) PetscFunctionReturn(0);
  if (!h) { ierr = VecGetArray(bv->buffer,&hh);CHKERRQ(ierr); }
  ierr = VecGetArrayRead(bv->omega,&omega);CHKERRQ(ierr);
  if (inverse) for (i=0;i<bv->nc+j;i++) hh[i] /= PetscRealPart(omega[i]);
  else for (i=0;i<bv->nc+j;i++) hh[i] *= PetscRealPart(omega[i]);
  ierr = VecRestoreArrayRead(bv->omega,&omega);CHKERRQ(ierr);
  if (!h) { ierr = VecRestoreArray(bv->buffer,&hh);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   BV_SquareRoot_Default - Returns the square root of position j (counted after the constraints)
   of the coefficients array
*/
PETSC_STATIC_INLINE PetscErrorCode BV_SquareRoot_Default(BV bv,PetscInt j,PetscScalar *h,PetscReal *beta)
{
  PetscErrorCode ierr;
  PetscScalar    *hh=h;

  PetscFunctionBegin;
  if (!h) { ierr = VecGetArray(bv->buffer,&hh);CHKERRQ(ierr); }
  ierr = BV_SafeSqrt(bv,hh[bv->nc+j],beta);CHKERRQ(ierr);
  if (!h) { ierr = VecRestoreArray(bv->buffer,&hh);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
   BV_StoreCoefficients_Default - Copy the contents of the coefficients array to an array dest
   provided by the caller (only values from l to j are copied)
*/
PETSC_STATIC_INLINE PetscErrorCode BV_StoreCoefficients_Default(BV bv,PetscInt j,PetscScalar *h,PetscScalar *dest)
{
  PetscErrorCode ierr;
  PetscScalar    *hh=h,*a;
  PetscInt       i;

  PetscFunctionBegin;
  if (!h) {
    ierr = VecGetArray(bv->buffer,&a);CHKERRQ(ierr);
    hh = a + j*(bv->nc+bv->m);
  }
  for (i=bv->l;i<j;i++) dest[i-bv->l] = hh[bv->nc+i];
  if (!h) { ierr = VecRestoreArray(bv->buffer,&a);CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}

/*
  BV_GetEigenvector - retrieves k-th eigenvector from basis vectors V.
  The argument eigi is the imaginary part of the corresponding eigenvalue.
*/
PETSC_STATIC_INLINE PetscErrorCode BV_GetEigenvector(BV V,PetscInt k,PetscScalar eigi,Vec Vr,Vec Vi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
#if defined(PETSC_USE_COMPLEX)
  if (Vr) { ierr = BVCopyVec(V,k,Vr);CHKERRQ(ierr); }
  if (Vi) { ierr = VecSet(Vi,0.0);CHKERRQ(ierr); }
#else
  if (eigi > 0.0) { /* first value of conjugate pair */
    if (Vr) { ierr = BVCopyVec(V,k,Vr);CHKERRQ(ierr); }
    if (Vi) { ierr = BVCopyVec(V,k+1,Vi);CHKERRQ(ierr); }
  } else if (eigi < 0.0) { /* second value of conjugate pair */
    if (Vr) { ierr = BVCopyVec(V,k-1,Vr);CHKERRQ(ierr); }
    if (Vi) {
      ierr = BVCopyVec(V,k,Vi);CHKERRQ(ierr);
      ierr = VecScale(Vi,-1.0);CHKERRQ(ierr);
    }
  } else { /* real eigenvalue */
    if (Vr) { ierr = BVCopyVec(V,k,Vr);CHKERRQ(ierr); }
    if (Vi) { ierr = VecSet(Vi,0.0);CHKERRQ(ierr); }
  }
#endif
  PetscFunctionReturn(0);
}

/*
   BV_OrthogonalizeColumn_Safe - this is intended for cases where we know that
   the resulting vector is going to be numerically zero, so normalization or
   iterative refinement may cause problems in parallel (collective operation
   not being called by all processes)
*/
PETSC_STATIC_INLINE PetscErrorCode BV_OrthogonalizeColumn_Safe(BV bv,PetscInt j,PetscScalar *H,PetscReal *norm,PetscBool *lindep)
{
  PetscErrorCode     ierr;
  BVOrthogRefineType orthog_ref;

  PetscFunctionBegin;
  ierr = PetscInfo1(bv,"Orthogonalizing column %D without refinement\n",j);CHKERRQ(ierr);
  orthog_ref     = bv->orthog_ref;
  bv->orthog_ref = BV_ORTHOG_REFINE_NEVER;  /* avoid refinement */
  ierr = BVOrthogonalizeColumn(bv,j,H,NULL,NULL);CHKERRQ(ierr);
  bv->orthog_ref = orthog_ref;  /* restore refinement setting */
  if (norm)   *norm  = 0.0;
  if (lindep) *lindep = PETSC_TRUE;
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_CUDA)
#define BV_CleanCoefficients(a,b,c)   ((a)->cuda?BV_CleanCoefficients_CUDA:BV_CleanCoefficients_Default)((a),(b),(c))
#define BV_AddCoefficients(a,b,c,d)   ((a)->cuda?BV_AddCoefficients_CUDA:BV_AddCoefficients_Default)((a),(b),(c),(d))
#define BV_SetValue(a,b,c,d,e)        ((a)->cuda?BV_SetValue_CUDA:BV_SetValue_Default)((a),(b),(c),(d),(e))
#define BV_SquareSum(a,b,c,d)         ((a)->cuda?BV_SquareSum_CUDA:BV_SquareSum_Default)((a),(b),(c),(d))
#define BV_ApplySignature(a,b,c,d)    ((a)->cuda?BV_ApplySignature_CUDA:BV_ApplySignature_Default)((a),(b),(c),(d))
#define BV_SquareRoot(a,b,c,d)        ((a)->cuda?BV_SquareRoot_CUDA:BV_SquareRoot_Default)((a),(b),(c),(d))
#define BV_StoreCoefficients(a,b,c,d) ((a)->cuda?BV_StoreCoefficients_CUDA:BV_StoreCoefficients_Default)((a),(b),(c),(d))
#else
#define BV_CleanCoefficients(a,b,c)   BV_CleanCoefficients_Default((a),(b),(c))
#define BV_AddCoefficients(a,b,c,d)   BV_AddCoefficients_Default((a),(b),(c),(d))
#define BV_SetValue(a,b,c,d,e)        BV_SetValue_Default((a),(b),(c),(d),(e))
#define BV_SquareSum(a,b,c,d)         BV_SquareSum_Default((a),(b),(c),(d))
#define BV_ApplySignature(a,b,c,d)    BV_ApplySignature_Default((a),(b),(c),(d))
#define BV_SquareRoot(a,b,c,d)        BV_SquareRoot_Default((a),(b),(c),(d))
#define BV_StoreCoefficients(a,b,c,d) BV_StoreCoefficients_Default((a),(b),(c),(d))
#endif /* PETSC_HAVE_CUDA */

#endif
