/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
   Various types of linearization for quadratic eigenvalue problem
*/

#include <slepc/private/pepimpl.h>
#include "linear.h"

/*
    Given the quadratic problem (l^2*M + l*C + K)*x = 0 the linearization is
    A*z = l*B*z for z = [  x  ] and A,B defined as follows:
                        [ l*x ]

            N:
                     A = [-bK      aI    ]     B = [ aI+bC   bM    ]
                         [-aK     -aC+bI ]         [ bI      aM    ]


            S:
                     A = [ bK      aK    ]     B = [ aK-bC  -bM    ]
                         [ aK      aC-bM ]         [-bM     -aM    ]


            H:
                     A = [ aK     -bK    ]     B = [ bM      aK+bC ]
                         [ aC+bM   aK    ]         [-aM      bM    ]
 */

/* - - - N - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

PetscErrorCode MatCreateExplicit_Linear_NA(MPI_Comm comm,PEP_LINEAR *ctx,Mat *A)
{
  PetscErrorCode ierr;
  PetscInt       M,N,m,n,i,Istart,Iend;
  Mat            Id,T=NULL;
  PetscReal      a=ctx->alpha,b=ctx->beta;
  PetscScalar    scalt=1.0;

  PetscFunctionBegin;
  ierr = MatGetSize(ctx->M,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(ctx->M,&m,&n);CHKERRQ(ierr);
  ierr = MatCreate(PetscObjectComm((PetscObject)ctx->M),&Id);CHKERRQ(ierr);
  ierr = MatSetSizes(Id,m,n,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(Id);CHKERRQ(ierr);
  ierr = MatSetUp(Id);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(Id,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(Id,i,i,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(Id,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Id,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (a!=0.0 && b!=0.0) {
    ierr = MatDuplicate(ctx->C,MAT_COPY_VALUES,&T);CHKERRQ(ierr);
    ierr = MatScale(T,-a*ctx->dsfactor*ctx->sfactor);CHKERRQ(ierr);
    ierr = MatShift(T,b);CHKERRQ(ierr);
  } else {
    if (a==0.0) { T = Id; scalt = b; }
    else { T = ctx->C; scalt = -a*ctx->dsfactor*ctx->sfactor; }
  }
  ierr = MatCreateTile(-b*ctx->dsfactor,ctx->K,a,Id,-ctx->dsfactor*a,ctx->K,scalt,T,A);CHKERRQ(ierr);
  ierr = MatDestroy(&Id);CHKERRQ(ierr);
  if (a!=0.0 && b!=0.0) {
    ierr = MatDestroy(&T);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MatCreateExplicit_Linear_NB(MPI_Comm comm,PEP_LINEAR *ctx,Mat *B)
{
  PetscErrorCode ierr;
  PetscInt       M,N,m,n,i,Istart,Iend;
  Mat            Id,T=NULL;
  PetscReal      a=ctx->alpha,b=ctx->beta;
  PetscScalar    scalt=1.0;

  PetscFunctionBegin;
  ierr = MatGetSize(ctx->M,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(ctx->M,&m,&n);CHKERRQ(ierr);
  ierr = MatCreate(PetscObjectComm((PetscObject)ctx->M),&Id);CHKERRQ(ierr);
  ierr = MatSetSizes(Id,m,n,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(Id);CHKERRQ(ierr);
  ierr = MatSetUp(Id);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(Id,&Istart,&Iend);CHKERRQ(ierr);
  for (i=Istart;i<Iend;i++) {
    ierr = MatSetValue(Id,i,i,1.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(Id,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Id,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (a!=0.0 && b!=0.0) {
    ierr = MatDuplicate(ctx->C,MAT_COPY_VALUES,&T);CHKERRQ(ierr);
    ierr = MatScale(T,b*ctx->dsfactor*ctx->sfactor);CHKERRQ(ierr);
    ierr = MatShift(T,a);CHKERRQ(ierr);
  } else {
    if (b==0.0) { T = Id; scalt = a; }
    else { T = ctx->C; scalt = b*ctx->dsfactor*ctx->sfactor; }
  }
  ierr = MatCreateTile(scalt,T,b*ctx->dsfactor*ctx->sfactor*ctx->sfactor,ctx->M,b,Id,a*ctx->sfactor*ctx->sfactor*ctx->dsfactor,ctx->M,B);CHKERRQ(ierr);
  ierr = MatDestroy(&Id);CHKERRQ(ierr);
  if (a!=0.0 && b!=0.0) {
    ierr = MatDestroy(&T);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* - - - S - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

PetscErrorCode MatCreateExplicit_Linear_SA(MPI_Comm comm,PEP_LINEAR *ctx,Mat *A)
{
  PetscErrorCode ierr;
  Mat            T=NULL;
  PetscScalar    scalt=1.0;
  PetscReal      a=ctx->alpha,b=ctx->beta;

  PetscFunctionBegin;
  if (a!=0.0 && b!=0.0) {
    ierr = MatDuplicate(ctx->C,MAT_COPY_VALUES,&T);CHKERRQ(ierr);
    ierr = MatScale(T,a*ctx->dsfactor*ctx->sfactor);CHKERRQ(ierr);
    ierr = MatAXPY(T,-b*ctx->dsfactor*ctx->sfactor*ctx->sfactor,ctx->M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  } else {
    if (a==0.0) { T = ctx->M; scalt = -b*ctx->dsfactor*ctx->sfactor*ctx->sfactor; }
    else { T = ctx->C; scalt = a*ctx->dsfactor*ctx->sfactor; }
  }
  ierr = MatCreateTile(b*ctx->dsfactor,ctx->K,a*ctx->dsfactor,ctx->K,a*ctx->dsfactor,ctx->K,scalt,T,A);CHKERRQ(ierr);
  if (a!=0.0 && b!=0.0) {
    ierr = MatDestroy(&T);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MatCreateExplicit_Linear_SB(MPI_Comm comm,PEP_LINEAR *ctx,Mat *B)
{
  PetscErrorCode ierr;
  Mat            T=NULL;
  PetscScalar    scalt=1.0;
  PetscReal      a=ctx->alpha,b=ctx->beta;

  PetscFunctionBegin;
  if (a!=0.0 && b!=0.0) {
    ierr = MatDuplicate(ctx->C,MAT_COPY_VALUES,&T);CHKERRQ(ierr);
    ierr = MatScale(T,-b*ctx->dsfactor*ctx->sfactor);CHKERRQ(ierr);
    ierr = MatAXPY(T,a*ctx->dsfactor,ctx->K,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  } else {
    if (b==0.0) { T = ctx->K; scalt = a*ctx->dsfactor; }
    else { T = ctx->C; scalt = -b*ctx->dsfactor*ctx->sfactor; }
  }
  ierr = MatCreateTile(scalt,T,-b*ctx->dsfactor*ctx->sfactor*ctx->sfactor,ctx->M,-b*ctx->dsfactor*ctx->sfactor*ctx->sfactor,ctx->M,-a*ctx->dsfactor*ctx->sfactor*ctx->sfactor,ctx->M,B);CHKERRQ(ierr);
  if (a!=0.0 && b!=0.0) {
    ierr = MatDestroy(&T);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* - - - H - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

PetscErrorCode MatCreateExplicit_Linear_HA(MPI_Comm comm,PEP_LINEAR *ctx,Mat *A)
{
  PetscErrorCode ierr;
  Mat            T=NULL;
  PetscScalar    scalt=1.0;
  PetscReal      a=ctx->alpha,b=ctx->beta;

  PetscFunctionBegin;
  if (a!=0.0 && b!=0.0) {
    ierr = MatDuplicate(ctx->C,MAT_COPY_VALUES,&T);CHKERRQ(ierr);
    ierr = MatScale(T,a*ctx->dsfactor*ctx->sfactor);CHKERRQ(ierr);
    ierr = MatAXPY(T,b*ctx->dsfactor*ctx->sfactor*ctx->sfactor,ctx->M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  } else {
    if (a==0.0) { T = ctx->M; scalt = b*ctx->dsfactor*ctx->sfactor*ctx->sfactor; }
    else { T = ctx->C; scalt = a*ctx->dsfactor*ctx->sfactor; }
  }
  ierr = MatCreateTile(a*ctx->dsfactor,ctx->K,-b*ctx->dsfactor,ctx->K,scalt,T,a*ctx->dsfactor,ctx->K,A);CHKERRQ(ierr);
  if (a!=0.0 && b!=0.0) {
    ierr = MatDestroy(&T);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MatCreateExplicit_Linear_HB(MPI_Comm comm,PEP_LINEAR *ctx,Mat *B)
{
  PetscErrorCode ierr;
  Mat            T=NULL;
  PetscScalar    scalt=1.0;
  PetscReal      a=ctx->alpha,b=ctx->beta;

  PetscFunctionBegin;
  if (a!=0.0 && b!=0.0) {
    ierr = MatDuplicate(ctx->C,MAT_COPY_VALUES,&T);CHKERRQ(ierr);
    ierr = MatScale(T,b*ctx->dsfactor*ctx->sfactor);CHKERRQ(ierr);
    ierr = MatAXPY(T,a*ctx->dsfactor,ctx->K,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  } else {
    if (b==0.0) { T = ctx->K; scalt = a*ctx->dsfactor; }
    else { T = ctx->C; scalt = b*ctx->dsfactor*ctx->sfactor; }
  }
  ierr = MatCreateTile(b*ctx->dsfactor*ctx->sfactor*ctx->sfactor,ctx->M,scalt,T,-a*ctx->dsfactor*ctx->sfactor*ctx->sfactor,ctx->M,b*ctx->dsfactor*ctx->sfactor*ctx->sfactor,ctx->M,B);CHKERRQ(ierr);
  if (a!=0.0 && b!=0.0) {
    ierr = MatDestroy(&T);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
