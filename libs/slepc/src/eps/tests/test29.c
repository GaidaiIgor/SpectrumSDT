/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2020, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Illustrates the computation of left eigenvectors for generalized eigenproblems.\n\n"
  "The command line options are:\n"
  "  -f1 <filename> -f2 <filename>, PETSc binary files containing A and B\n\n";

#include <slepceps.h>

/*
   User-defined routines
*/
PetscErrorCode ComputeResidualNorm(Mat,Mat,PetscBool,PetscScalar,PetscScalar,Vec,Vec,Vec*,PetscReal*);

int main(int argc,char **argv)
{
  Mat            A,B;
  EPS            eps;
  EPSType        type;
  PetscInt       i,nconv;
  PetscBool      twosided,flg;
  PetscReal      nrmr,nrml=0.0,re,im,lev;
  PetscScalar    *kr,*ki;
  Vec            t,*xr,*xi,*yr,*yi,*z;
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, Ax=kBx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nGeneralized eigenproblem stored in file.\n\n");CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-f1",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate a file name for matrix A with the -f1 option");

#if defined(PETSC_USE_COMPLEX)
  ierr = PetscPrintf(PETSC_COMM_WORLD," Reading COMPLEX matrices from binary files...\n");CHKERRQ(ierr);
#else
  ierr = PetscPrintf(PETSC_COMM_WORLD," Reading REAL matrices from binary files...\n");CHKERRQ(ierr);
#endif
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatLoad(A,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = PetscOptionsGetString(NULL,NULL,"-f2",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
    ierr = MatSetFromOptions(B);CHKERRQ(ierr);
    ierr = MatLoad(B,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD," Matrix B was not provided, setting B=I\n\n");CHKERRQ(ierr);
    B = NULL;
  }
  ierr = MatCreateVecs(A,NULL,&t);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);

  /* use a two-sided algorithm to compute left eigenvectors as well */
  ierr = EPSSetTwoSided(eps,PETSC_TRUE);CHKERRQ(ierr);

  /* allow user to change settings at run time */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  ierr = EPSGetTwoSided(eps,&twosided);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Get number of converged approximate eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);
  ierr = PetscMalloc2(nconv,&kr,nconv,&ki);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(t,3,&z);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(t,nconv,&xr);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(t,nconv,&xi);CHKERRQ(ierr);
  if (twosided) {
    ierr = VecDuplicateVecs(t,nconv,&yr);CHKERRQ(ierr);
    ierr = VecDuplicateVecs(t,nconv,&yi);CHKERRQ(ierr);
  }

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k            ||Ax-kBx||         ||y'A-y'Bk||\n"
         "   ---------------- ------------------ ------------------\n");CHKERRQ(ierr);

    for (i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,i,&kr[i],&ki[i],xr[i],xi[i]);CHKERRQ(ierr);
      if (twosided) {
        ierr = EPSGetLeftEigenvector(eps,i,yr[i],yi[i]);CHKERRQ(ierr);
      }
      /*
         Compute the residual norms associated to each eigenpair
      */
      ierr = ComputeResidualNorm(A,B,PETSC_FALSE,kr[i],ki[i],xr[i],xi[i],z,&nrmr);CHKERRQ(ierr);
      if (twosided) {
        ierr = ComputeResidualNorm(A,B,PETSC_TRUE,kr[i],ki[i],yr[i],yi[i],z,&nrml);CHKERRQ(ierr);
      }

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr[i]);
      im = PetscImaginaryPart(kr[i]);
#else
      re = kr[i];
      im = ki[i];
#endif
      if (im!=0.0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," %8f%+8fi %12g       %12g\n",(double)re,(double)im,(double)nrmr,(double)nrml);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g       %12g\n",(double)re,(double)nrmr,(double)nrml);CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    /*
       Check bi-orthogonality of eigenvectors
    */
    if (twosided) {
      ierr = VecCheckOrthogonality(xr,nconv,yr,nconv,B,NULL,&lev);CHKERRQ(ierr);
      if (lev<100*PETSC_MACHINE_EPSILON) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of bi-orthogonality of eigenvectors < 100*eps\n\n");CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Level of bi-orthogonality of eigenvectors: %g\n\n",(double)lev);CHKERRQ(ierr);
      }
    }
  }

  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = VecDestroy(&t);CHKERRQ(ierr);
  ierr = PetscFree2(kr,ki);CHKERRQ(ierr);
  ierr = VecDestroyVecs(3,&z);CHKERRQ(ierr);
  ierr = VecDestroyVecs(nconv,&xr);CHKERRQ(ierr);
  ierr = VecDestroyVecs(nconv,&xi);CHKERRQ(ierr);
  if (twosided) {
    ierr = VecDestroyVecs(nconv,&yr);CHKERRQ(ierr);
    ierr = VecDestroyVecs(nconv,&yi);CHKERRQ(ierr);
  }
  ierr = SlepcFinalize();
  return ierr;
}

/*
   ComputeResidualNorm - Computes the norm of the residual vector
   associated with an eigenpair.

   Input Parameters:
     trans - whether A' must be used instead of A
     kr,ki - eigenvalue
     xr,xi - eigenvector
     z     - three work vectors (the second one not referenced in complex scalars)
*/
PetscErrorCode ComputeResidualNorm(Mat A,Mat B,PetscBool trans,PetscScalar kr,PetscScalar ki,Vec xr,Vec xi,Vec *z,PetscReal *norm)
{
  PetscErrorCode ierr;
  Vec            u,w=NULL;
  PetscScalar    alpha;
#if !defined(PETSC_USE_COMPLEX)
  Vec            v;
  PetscReal      ni,nr;
#endif
  PetscErrorCode (*matmult)(Mat,Vec,Vec) = trans? MatMultHermitianTranspose: MatMult;

  PetscFunctionBegin;
  u = z[0];
  if (B) w = z[2];

#if !defined(PETSC_USE_COMPLEX)
  v = z[1];
  if (ki == 0 || PetscAbsScalar(ki) < PetscAbsScalar(kr*PETSC_MACHINE_EPSILON)) {
#endif
    ierr = (*matmult)(A,xr,u);CHKERRQ(ierr);                          /* u=A*x */
    if (PetscAbsScalar(kr) > PETSC_MACHINE_EPSILON) {
      if (B) { ierr = (*matmult)(B,xr,w);CHKERRQ(ierr); }             /* w=B*x */
      else w = xr;
      alpha = trans? -PetscConj(kr): -kr;
      ierr = VecAXPY(u,alpha,w);CHKERRQ(ierr);                        /* u=A*x-k*B*x */
    }
    ierr = VecNorm(u,NORM_2,norm);CHKERRQ(ierr);
#if !defined(PETSC_USE_COMPLEX)
  } else {
    ierr = (*matmult)(A,xr,u);CHKERRQ(ierr);                          /* u=A*xr */
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      if (B) { ierr = (*matmult)(B,xr,v);CHKERRQ(ierr); }             /* v=B*xr */
      else { ierr = VecCopy(xr,v);CHKERRQ(ierr); }
      ierr = VecAXPY(u,-kr,v);CHKERRQ(ierr);                          /* u=A*xr-kr*B*xr */
      if (B) { ierr = (*matmult)(B,xi,w);CHKERRQ(ierr); }             /* w=B*xi */
      else w = xi;
      ierr = VecAXPY(u,trans?-ki:ki,w);CHKERRQ(ierr);                 /* u=A*xr-kr*B*xr+ki*B*xi */
    }
    ierr = VecNorm(u,NORM_2,&nr);CHKERRQ(ierr);
    ierr = (*matmult)(A,xi,u);CHKERRQ(ierr);                          /* u=A*xi */
    if (SlepcAbsEigenvalue(kr,ki) > PETSC_MACHINE_EPSILON) {
      ierr = VecAXPY(u,-kr,w);CHKERRQ(ierr);                          /* u=A*xi-kr*B*xi */
      ierr = VecAXPY(u,trans?ki:-ki,v);CHKERRQ(ierr);                 /* u=A*xi-kr*B*xi-ki*B*xr */
    }
    ierr = VecNorm(u,NORM_2,&ni);CHKERRQ(ierr);
    *norm = SlepcAbsEigenvalue(nr,ni);
  }
#endif
  PetscFunctionReturn(0);
}

/*TEST

   testset:
      args: -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -eps_nev 4 -st_type sinvert -eps_target -190000
      filter: grep -v "method" | sed -e "s/[+-]0.0*i//" | sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: double !complex !define(PETSC_USE_64BIT_INDICES)
      test:
         suffix: 1
      test:
         suffix: 1_rqi
         TODO: gives differences after adding EPSGetLeftStartVector()
         args: -eps_type power -eps_power_shift_type rayleigh
         output_file: output/test29_1.out
      test:
         suffix: 1_rqi_singular
         args: -eps_type power -eps_power_shift_type rayleigh -eps_nev 1 -eps_target -195500

   test:
      suffix: 2
      args: -f1 ${DATAFILESPATH}/matrices/complex/mhd1280a.petsc -f2 ${DATAFILESPATH}/matrices/complex/mhd1280b.petsc -eps_nev 6 -eps_tol 1e-11
      filter: sed -e "s/[0-9]\.[0-9]*e[+-]\([0-9]*\)/removed/g"
      requires: complex datafilespath !define(PETSC_USE_64BIT_INDICES)
      timeoutfactor: 2

TEST*/
