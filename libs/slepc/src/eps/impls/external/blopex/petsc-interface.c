/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1) LGPL Version 2.1 or above.See www.gnu.org. */
/* @@@ Copyright 2010 BLOPEX team https://github.com/lobpcg/blopex     */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#include <petscvec.h>
#include <petscblaslapack.h>
#include <interpreter.h>
#include <temp_multivector.h>
#include <fortran_matrix.h>

static PetscRandom LOBPCG_RandomContext = NULL;

#if !defined(PETSC_USE_COMPLEX)
BlopexInt PETSC_dpotrf_interface (char *uplo,BlopexInt *n,double *a,BlopexInt * lda,BlopexInt *info)
{
  PetscBLASInt n_,lda_,info_;

  /* type conversion */
  n_ = *n;
  lda_ = *lda;
  info_ = *info;

  LAPACKpotrf_(uplo,&n_,(PetscScalar*)a,&lda_,&info_);

  *info = info_;
  return 0;
}

BlopexInt PETSC_dsygv_interface (BlopexInt *itype,char *jobz,char *uplo,BlopexInt *n,double *a,BlopexInt *lda,double *b,BlopexInt *ldb,double *w,double *work,BlopexInt *lwork,BlopexInt *info)
{
  PetscBLASInt itype_,n_,lda_,ldb_,lwork_,info_;

  itype_ = *itype;
  n_ = *n;
  lda_ = *lda;
  ldb_ = *ldb;
  lwork_ = *lwork;
  info_ = *info;

  LAPACKsygv_(&itype_,jobz,uplo,&n_,(PetscScalar*)a,&lda_,(PetscScalar*)b,&ldb_,(PetscScalar*)w,(PetscScalar*)work,&lwork_,&info_);

  *info = info_;
  return 0;
}
#else
BlopexInt PETSC_zpotrf_interface (char *uplo,BlopexInt *n,komplex *a,BlopexInt* lda,BlopexInt *info)
{
  PetscBLASInt n_,lda_,info_;

  /* type conversion */
  n_ = *n;
  lda_ = (PetscBLASInt)*lda;

  LAPACKpotrf_(uplo,&n_,(PetscScalar*)a,&lda_,&info_);

  *info = info_;
  return 0;
}

BlopexInt PETSC_zsygv_interface (BlopexInt *itype,char *jobz,char *uplo,BlopexInt *n,komplex *a,BlopexInt *lda,komplex *b,BlopexInt *ldb,double *w,komplex *work,BlopexInt *lwork,double *rwork,BlopexInt *info)
{
  PetscBLASInt itype_,n_,lda_,ldb_,lwork_,info_;

  itype_ = *itype;
  n_ = *n;
  lda_ = *lda;
  ldb_ = *ldb;
  lwork_ = *lwork;
  info_ = *info;

  LAPACKsygv_(&itype_,jobz,uplo,&n_,(PetscScalar*)a,&lda_,(PetscScalar*)b,&ldb_,(PetscReal*)w,(PetscScalar*)work,&lwork_,(PetscReal*)rwork,&info_);

  *info = info_;
  return 0;
}
#endif

void *PETSC_MimicVector(void *vvector)
{
  PetscErrorCode  ierr;
  Vec temp;

  ierr = VecDuplicate((Vec)vvector,&temp);CHKERRABORT(PETSC_COMM_SELF,ierr);
  return (void*)temp;
}

BlopexInt PETSC_DestroyVector(void *vvector)
{
  PetscErrorCode ierr;
  Vec v = (Vec)vvector;

  ierr = VecDestroy(&v);CHKERRQ(ierr);
  return 0;
}

BlopexInt PETSC_InnerProd(void *x,void *y,void *result)
{
  PetscErrorCode ierr;

  ierr = VecDot((Vec)x,(Vec)y,(PetscScalar*)result);CHKERRQ(ierr);
  return 0;
}

BlopexInt PETSC_CopyVector(void *x,void *y)
{
  PetscErrorCode  ierr;

  ierr = VecCopy((Vec)x,(Vec)y);CHKERRQ(ierr);
  return 0;
}

BlopexInt PETSC_ClearVector(void *x)
{
  PetscErrorCode  ierr;

  ierr = VecSet((Vec)x,0.0);CHKERRQ(ierr);
  return 0;
}

BlopexInt PETSC_SetRandomValues(void* v,BlopexInt seed)
{
  PetscErrorCode ierr;

  /* note: without previous call to LOBPCG_InitRandomContext LOBPCG_RandomContext will be null,
    and VecSetRandom will use internal petsc random context */

  ierr = VecSetRandom((Vec)v,LOBPCG_RandomContext);CHKERRQ(ierr);
  return 0;
}

BlopexInt PETSC_ScaleVector(double alpha,void *x)
{
  PetscErrorCode ierr;

  ierr = VecScale((Vec)x,alpha);CHKERRQ(ierr);
  return 0;
}

BlopexInt PETSC_Axpy(void *alpha,void *x,void *y)
{
  PetscErrorCode ierr;

  ierr = VecAXPY((Vec)y,*(PetscScalar*)alpha,(Vec)x);CHKERRQ(ierr);
  return 0;
}

BlopexInt PETSC_VectorSize(void *x)
{
  PetscInt N;
  VecGetSize((Vec)x,&N);
  return N;
}

int LOBPCG_InitRandomContext(MPI_Comm comm,PetscRandom rand)
{
  PetscErrorCode ierr;
  /* PetscScalar rnd_bound = 1.0; */

  if (rand) {
    ierr = PetscObjectReference((PetscObject)rand);CHKERRQ(ierr);
    ierr = PetscRandomDestroy(&LOBPCG_RandomContext);CHKERRQ(ierr);
    LOBPCG_RandomContext = rand;
  } else {
    ierr = PetscRandomCreate(comm,&LOBPCG_RandomContext);CHKERRQ(ierr);
  }
  return 0;
}

int LOBPCG_SetFromOptionsRandomContext(void)
{
  PetscErrorCode ierr;
  ierr = PetscRandomSetFromOptions(LOBPCG_RandomContext);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
  ierr = PetscRandomSetInterval(LOBPCG_RandomContext,(PetscScalar)PetscCMPLX(-1.0,-1.0),(PetscScalar)PetscCMPLX(1.0,1.0));CHKERRQ(ierr);
#else
  ierr = PetscRandomSetInterval(LOBPCG_RandomContext,(PetscScalar)-1.0,(PetscScalar)1.0);CHKERRQ(ierr);
#endif
  return 0;
}

int LOBPCG_DestroyRandomContext(void)
{
  PetscErrorCode ierr;

  ierr = PetscRandomDestroy(&LOBPCG_RandomContext);CHKERRQ(ierr);
  return 0;
}

int PETSCSetupInterpreter(mv_InterfaceInterpreter *i)
{
  i->CreateVector = PETSC_MimicVector;
  i->DestroyVector = PETSC_DestroyVector;
  i->InnerProd = PETSC_InnerProd;
  i->CopyVector = PETSC_CopyVector;
  i->ClearVector = PETSC_ClearVector;
  i->SetRandomValues = PETSC_SetRandomValues;
  i->ScaleVector = PETSC_ScaleVector;
  i->Axpy = PETSC_Axpy;
  i->VectorSize = PETSC_VectorSize;

  /* Multivector part */

  i->CreateMultiVector = mv_TempMultiVectorCreateFromSampleVector;
  i->CopyCreateMultiVector = mv_TempMultiVectorCreateCopy;
  i->DestroyMultiVector = mv_TempMultiVectorDestroy;

  i->Width = mv_TempMultiVectorWidth;
  i->Height = mv_TempMultiVectorHeight;
  i->SetMask = mv_TempMultiVectorSetMask;
  i->CopyMultiVector = mv_TempMultiVectorCopy;
  i->ClearMultiVector = mv_TempMultiVectorClear;
  i->SetRandomVectors = mv_TempMultiVectorSetRandom;
  i->Eval = mv_TempMultiVectorEval;

#if defined(PETSC_USE_COMPLEX)
  i->MultiInnerProd = mv_TempMultiVectorByMultiVector_complex;
  i->MultiInnerProdDiag = mv_TempMultiVectorByMultiVectorDiag_complex;
  i->MultiVecMat = mv_TempMultiVectorByMatrix_complex;
  i->MultiVecMatDiag = mv_TempMultiVectorByDiagonal_complex;
  i->MultiAxpy = mv_TempMultiVectorAxpy_complex;
  i->MultiXapy = mv_TempMultiVectorXapy_complex;
#else
  i->MultiInnerProd = mv_TempMultiVectorByMultiVector;
  i->MultiInnerProdDiag = mv_TempMultiVectorByMultiVectorDiag;
  i->MultiVecMat = mv_TempMultiVectorByMatrix;
  i->MultiVecMatDiag = mv_TempMultiVectorByDiagonal;
  i->MultiAxpy = mv_TempMultiVectorAxpy;
  i->MultiXapy = mv_TempMultiVectorXapy;
#endif

  return 0;
}
