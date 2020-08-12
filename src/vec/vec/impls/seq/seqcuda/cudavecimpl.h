#if !defined(__CUDAVECIMPL)
#define __CUDAVECIMPL

#include <petscvec.h>
#include <petsccublas.h>
#include <petsc/private/vecimpl.h>

typedef struct {
  PetscScalar  *GPUarray;           /* this always holds the GPU data */
  PetscScalar  *GPUarray_allocated; /* if the array was allocated by PETSc this is its pointer */
  cudaStream_t stream;              /* A stream for doing asynchronous data transfers */
} Vec_CUDA;

PETSC_INTERN PetscErrorCode VecDotNorm2_SeqCUDA(Vec,Vec,PetscScalar*, PetscScalar*);
PETSC_INTERN PetscErrorCode VecPointwiseDivide_SeqCUDA(Vec,Vec,Vec);
PETSC_INTERN PetscErrorCode VecWAXPY_SeqCUDA(Vec,PetscScalar,Vec,Vec);
PETSC_INTERN PetscErrorCode VecMDot_SeqCUDA(Vec,PetscInt,const Vec[],PetscScalar*);
PETSC_EXTERN PetscErrorCode VecSet_SeqCUDA(Vec,PetscScalar);
PETSC_INTERN PetscErrorCode VecMAXPY_SeqCUDA(Vec,PetscInt,const PetscScalar*,Vec*);
PETSC_INTERN PetscErrorCode VecAXPBYPCZ_SeqCUDA(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
PETSC_INTERN PetscErrorCode VecPointwiseMult_SeqCUDA(Vec,Vec,Vec);
PETSC_INTERN PetscErrorCode VecPlaceArray_SeqCUDA(Vec,const PetscScalar*);
PETSC_INTERN PetscErrorCode VecResetArray_SeqCUDA(Vec);
PETSC_INTERN PetscErrorCode VecReplaceArray_SeqCUDA(Vec,const PetscScalar*);
PETSC_INTERN PetscErrorCode VecDot_SeqCUDA(Vec,Vec,PetscScalar*);
PETSC_INTERN PetscErrorCode VecTDot_SeqCUDA(Vec,Vec,PetscScalar*);
PETSC_INTERN PetscErrorCode VecScale_SeqCUDA(Vec,PetscScalar);
PETSC_EXTERN PetscErrorCode VecCopy_SeqCUDA(Vec,Vec);
PETSC_INTERN PetscErrorCode VecSwap_SeqCUDA(Vec,Vec);
PETSC_EXTERN PetscErrorCode VecAXPY_SeqCUDA(Vec,PetscScalar,Vec);
PETSC_INTERN PetscErrorCode VecAXPBY_SeqCUDA(Vec,PetscScalar,PetscScalar,Vec);
PETSC_INTERN PetscErrorCode VecDuplicate_SeqCUDA(Vec,Vec*);
PETSC_INTERN PetscErrorCode VecConjugate_SeqCUDA(Vec xin);
PETSC_INTERN PetscErrorCode VecNorm_SeqCUDA(Vec,NormType,PetscReal*);
PETSC_INTERN PetscErrorCode VecCUDACopyToGPU(Vec);
PETSC_INTERN PetscErrorCode VecCUDAAllocateCheck(Vec);
PETSC_EXTERN PetscErrorCode VecCreate_SeqCUDA(Vec);
PETSC_INTERN PetscErrorCode VecCreate_SeqCUDA_Private(Vec,const PetscScalar*);
PETSC_INTERN PetscErrorCode VecCreate_MPICUDA(Vec);
PETSC_INTERN PetscErrorCode VecCreate_MPICUDA_Private(Vec,PetscBool,PetscInt,const PetscScalar*);
PETSC_INTERN PetscErrorCode VecCreate_CUDA(Vec);
PETSC_INTERN PetscErrorCode VecDestroy_SeqCUDA(Vec);
PETSC_INTERN PetscErrorCode VecDestroy_MPICUDA(Vec);
PETSC_INTERN PetscErrorCode VecAYPX_SeqCUDA(Vec,PetscScalar,Vec);
PETSC_INTERN PetscErrorCode VecSetRandom_SeqCUDA(Vec,PetscRandom);
PETSC_INTERN PetscErrorCode VecGetLocalVector_SeqCUDA(Vec,Vec);
PETSC_INTERN PetscErrorCode VecRestoreLocalVector_SeqCUDA(Vec,Vec);
PETSC_INTERN PetscErrorCode VecGetArrayWrite_SeqCUDA(Vec,PetscScalar**);
PETSC_INTERN PetscErrorCode VecCopy_SeqCUDA_Private(Vec xin,Vec yin);
PETSC_INTERN PetscErrorCode VecSetRandom_SeqCUDA_Private(Vec xin,PetscRandom r);
PETSC_INTERN PetscErrorCode VecDestroy_SeqCUDA_Private(Vec v);
PETSC_INTERN PetscErrorCode VecResetArray_SeqCUDA_Private(Vec vin);
PETSC_INTERN PetscErrorCode VecCUDACopyToGPU_Public(Vec);
PETSC_INTERN PetscErrorCode VecCUDAAllocateCheck_Public(Vec);
PETSC_INTERN PetscErrorCode VecCUDACopyToGPUSome(Vec,PetscCUDAIndices,ScatterMode);
PETSC_INTERN PetscErrorCode VecCUDACopyFromGPUSome(Vec,PetscCUDAIndices,ScatterMode);

PETSC_INTERN PetscErrorCode VecScatterCUDAIndicesCreate_PtoP(PetscInt, PetscInt*,PetscInt, PetscInt*,PetscCUDAIndices*);
PETSC_INTERN PetscErrorCode VecScatterCUDAIndicesCreate_StoS(PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt*,PetscInt*,PetscCUDAIndices*);
PETSC_INTERN PetscErrorCode VecScatterCUDAIndicesDestroy(PetscCUDAIndices*);
PETSC_INTERN PetscErrorCode VecScatterCUDA_StoS(Vec,Vec,PetscCUDAIndices,InsertMode,ScatterMode);

typedef enum {VEC_SCATTER_CUDA_STOS, VEC_SCATTER_CUDA_PTOP} VecCUDAScatterType;
typedef enum {VEC_SCATTER_CUDA_GENERAL, VEC_SCATTER_CUDA_STRIDED} VecCUDASequentialScatterMode;

struct  _p_VecScatterCUDAIndices_PtoP {
  PetscInt ns;
  PetscInt sendLowestIndex;
  PetscInt nr;
  PetscInt recvLowestIndex;
};

struct  _p_VecScatterCUDAIndices_StoS {
  /* from indices data */
  PetscInt *fslots;
  PetscInt fromFirst;
  PetscInt fromStep;
  VecCUDASequentialScatterMode fromMode;

  /* to indices data */
  PetscInt *tslots;
  PetscInt toFirst;
  PetscInt toStep;
  VecCUDASequentialScatterMode toMode;

  PetscInt n;
  PetscInt MAX_BLOCKS;
  PetscInt MAX_CORESIDENT_THREADS;
  cudaStream_t stream;
};

struct  _p_PetscCUDAIndices {
  void * scatter;
  VecCUDAScatterType scatterType;
};

/* complex single */
#if defined(PETSC_USE_COMPLEX)
#if defined(PETSC_USE_REAL_SINGLE)
#define cublasXaxpy(a,b,c,d,e,f,g) cublasCaxpy((a),(b),(cuComplex*)(c),(cuComplex*)(d),(e),(cuComplex*)(f),(g))
#define cublasXscal(a,b,c,d,e)     cublasCscal((a),(b),(cuComplex*)(c),(cuComplex*)(d),(e))
#define cublasXdotu(a,b,c,d,e,f,g) cublasCdotu((a),(b),(cuComplex*)(c),(d),(cuComplex*)(e),(f),(cuComplex*)(g))
#define cublasXdot(a,b,c,d,e,f,g)  cublasCdotc((a),(b),(cuComplex*)(c),(d),(cuComplex*)(e),(f),(cuComplex*)(g))
#define cublasXswap(a,b,c,d,e,f)   cublasCswap((a),(b),(cuComplex*)(c),(d),(cuComplex*)(e),(f))
#define cublasXnrm2(a,b,c,d,e)     cublasScnrm2((a),(b),(cuComplex*)(c),(d),(e))
#define cublasIXamax(a,b,c,d,e)    cublasIcamax((a),(b),(cuComplex*)(c),(d),(e))
#define cublasXasum(a,b,c,d,e)     cublasScasum((a),(b),(cuComplex*)(c),(d),(e))
#define cublasXgemv(a,b,c,d,e,f,g,h,i,j,k,l) cublasCgemv((a),(b),(c),(d),(cuComplex*)(e),(cuComplex*)(f),(g),(cuComplex*)(h),(i),(cuComplex*)(j),(cuComplex*)(k),(l))
#define cublasXgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n) cublasCgemm((a),(b),(c),(d),(e),(f),(cuComplex*)(g),(cuComplex*)(h),(i),(cuComplex*)(j),(k),(cuComplex*)(l),(cuComplex*)(m),(n))
#else /* complex double */
#define cublasXaxpy(a,b,c,d,e,f,g) cublasZaxpy((a),(b),(cuDoubleComplex*)(c),(cuDoubleComplex*)(d),(e),(cuDoubleComplex*)(f),(g))
#define cublasXscal(a,b,c,d,e)     cublasZscal((a),(b),(cuDoubleComplex*)(c),(cuDoubleComplex*)(d),(e))
#define cublasXdotu(a,b,c,d,e,f,g) cublasZdotu((a),(b),(cuDoubleComplex*)(c),(d),(cuDoubleComplex*)(e),(f),(cuDoubleComplex*)(g))
#define cublasXdot(a,b,c,d,e,f,g)  cublasZdotc((a),(b),(cuDoubleComplex*)(c),(d),(cuDoubleComplex*)(e),(f),(cuDoubleComplex*)(g))
#define cublasXswap(a,b,c,d,e,f)   cublasZswap((a),(b),(cuDoubleComplex*)(c),(d),(cuDoubleComplex*)(e),(f))
#define cublasXnrm2(a,b,c,d,e)     cublasDznrm2((a),(b),(cuDoubleComplex*)(c),(d),(e))
#define cublasIXamax(a,b,c,d,e)    cublasIzamax((a),(b),(cuDoubleComplex*)(c),(d),(e))
#define cublasXasum(a,b,c,d,e)     cublasDzasum((a),(b),(cuDoubleComplex*)(c),(d),(e))
#define cublasXgemv(a,b,c,d,e,f,g,h,i,j,k,l) cublasZgemv((a),(b),(c),(d),(cuDoubleComplex*)(e),(cuDoubleComplex*)(f),(g),(cuDoubleComplex*)(h),(i),(cuDoubleComplex*)(j),(cuDoubleComplex*)(k),(l))
#define cublasXgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n) cublasZgemm((a),(b),(c),(d),(e),(f),(cuDoubleComplex*)(g),(cuDoubleComplex*)(h),(i),(cuDoubleComplex*)(j),(k),(cuDoubleComplex*)(l),(cuDoubleComplex*)(m),(n))
#endif
#else /* real single */
#if defined(PETSC_USE_REAL_SINGLE)
#define cublasXaxpy  cublasSaxpy
#define cublasXscal  cublasSscal
#define cublasXdotu  cublasSdot
#define cublasXdot   cublasSdot
#define cublasXswap  cublasSswap
#define cublasXnrm2  cublasSnrm2
#define cublasIXamax cublasIsamax
#define cublasXasum  cublasSasum
#define cublasXgemv  cublasSgemv
#define cublasXgemm  cublasSgemm
#else /* real double */
#define cublasXaxpy  cublasDaxpy
#define cublasXscal  cublasDscal
#define cublasXdotu  cublasDdot
#define cublasXdot   cublasDdot
#define cublasXswap  cublasDswap
#define cublasXnrm2  cublasDnrm2
#define cublasIXamax cublasIdamax
#define cublasXasum  cublasDasum
#define cublasXgemv  cublasDgemv
#define cublasXgemm  cublasDgemm
#endif
#endif

#endif
